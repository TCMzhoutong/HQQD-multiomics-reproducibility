from pathlib import Path
from io import StringIO
import time

import pandas as pd
from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.selenium_manager import SeleniumManager


BASE = Path(__file__).resolve().parent
INPUT = BASE / "04.GutMicrobe_metabolite_targets" / "02_metabolite_SMILES.csv"
DOWNLOAD_DIR = BASE / "metabolite_swisstarget_downloads"
ORIG_OUT = BASE / "04.GutMicrobe_metabolite_targets" / "03_metabolite_target_orig.csv"
CLEAN_OUT = BASE / "04.GutMicrobe_metabolite_targets" / "03_metabolite_target.csv"
UNIQUE_OUT = BASE / "05.disease_targets" / "03_metabolite_targets.csv"


def make_driver():
    options = webdriver.ChromeOptions()
    options.page_load_strategy = "eager"
    options.add_argument("--headless=new")
    options.add_experimental_option("excludeSwitches", ["enable-logging"])
    prefs = {
        "download.default_directory": str(DOWNLOAD_DIR),
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": True,
    }
    options.add_experimental_option("prefs", prefs)
    paths = SeleniumManager().binary_paths(
        ["--browser", "chrome", "--skip-driver-in-path", "--avoid-stats"]
    )
    service = webdriver.ChromeService(executable_path=paths["driver_path"])
    driver = webdriver.Chrome(options=options, service=service)
    try:
        driver.command_executor._client_config.timeout = 600
    except Exception:
        pass
    return driver


def normalize_cid(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def predict_one(driver, row, index):
    smiles = str(row["SMILES"]).strip()
    cid = normalize_cid(row.get("PubChem_CID"))
    file_id = cid or f"compound_{index + 1}"
    name = row.get("CM_Name") or row.get("Name") or file_id
    out_file = DOWNLOAD_DIR / f"{file_id}_SwissTargetPrediction.csv"

    if out_file.exists() and out_file.stat().st_size > 0:
        try:
            existing_rows = len(pd.read_csv(out_file))
        except Exception:
            existing_rows = 0
        if existing_rows >= 50:
            print(f"[skip] {index + 1}: {name} -> existing {out_file.name} ({existing_rows} rows)")
            return out_file
        print(f"[redo] {index + 1}: {name} -> existing {out_file.name} has only {existing_rows} rows")

    print(f"[run] {index + 1}: {name} ({file_id})")
    driver.get("https://www.swisstargetprediction.ch/index.php")
    try:
        alert = WebDriverWait(driver, 5).until(EC.alert_is_present())
        print(f"      accepted site alert: {alert.text}")
        alert.accept()
    except TimeoutException:
        pass
    wait = WebDriverWait(driver, 90)
    wait.until(EC.presence_of_element_located((By.ID, "smilesBox")))

    state = driver.execute_script(
        """
        const smilesBox = document.getElementById('smilesBox');
        smilesBox.value = arguments[0];
        if (typeof checkForm === 'function') {
            checkForm();
        } else {
            document.getElementById('submitButton').disabled = false;
        }
        return {
          value: smilesBox.value,
          disabled: document.getElementById('submitButton').disabled
        };
        """,
        smiles,
    )
    print(f"      filled={bool(state['value'])} submit_disabled={state['disabled']}")
    if state["disabled"]:
        raise RuntimeError(f"Submit button stayed disabled for {name}")

    driver.execute_script("document.getElementById('submitButton').click();")
    print("      submitted")

    try:
        WebDriverWait(driver, 420).until(
            lambda d: "result.php" in d.current_url
            and len(d.find_elements(By.ID, "resultTable")) > 0
        )
    except TimeoutException:
        print(f"      timeout url={driver.current_url}")
        raise

    print(f"      result url={driver.current_url}")
    try:
        select = Select(driver.find_element(By.NAME, "resultTable_length"))
        select.select_by_visible_text("All")
        WebDriverWait(driver, 30).until(
            lambda d: "Showing 1 to" in d.find_element(By.ID, "resultTable_info").text
            and "of" in d.find_element(By.ID, "resultTable_info").text
        )
        print(f"      {driver.find_element(By.ID, 'resultTable_info').text}")
    except Exception as exc:
        print(f"      could not switch table length to All: {exc}")
    tables = pd.read_html(StringIO(driver.page_source), attrs={"id": "resultTable"})
    if not tables:
        raise RuntimeError(f"No resultTable parsed for {name}")

    df = tables[0]
    df.to_csv(out_file, index=False, encoding="utf-8-sig")
    print(f"      saved {out_file.name}: {len(df)} rows")
    return out_file


def rebuild_metabolite_target_files(input_df):
    files = []
    for index, row in input_df.iterrows():
        cid = normalize_cid(row.get("PubChem_CID"))
        file_id = cid or f"compound_{index + 1}"
        file_path = DOWNLOAD_DIR / f"{file_id}_SwissTargetPrediction.csv"
        if file_path.exists() and file_path.stat().st_size > 0:
            files.append((index, row, file_path))
        else:
            print(f"[warn] missing target file for {row.get('CM_Name') or row.get('Name')}: {file_path.name}")

    rows = []
    for _, row, file_path in files:
        target_df = pd.read_csv(file_path)
        target_df["Probability*"] = pd.to_numeric(target_df["Probability*"], errors="coerce")
        target_df = target_df[target_df["Probability*"] > 0].copy()
        for _, target in target_df.iterrows():
            rows.append(
                {
                    "CM_Name": row.get("CM_Name") or row.get("Name"),
                    "Pubchem_CID": normalize_cid(row.get("PubChem_CID")),
                    "SMILES": row.get("SMILES"),
                    "Target": target.get("Target"),
                    "Common name": target.get("Common name"),
                    "Uniprot ID": target.get("Uniprot ID"),
                    "Target Class": target.get("Target Class"),
                    "Probability*": target.get("Probability*"),
                }
            )

    orig_df = pd.DataFrame(rows)
    orig_df.to_csv(ORIG_OUT, index=False, encoding="utf-8-sig")

    if orig_df.empty:
        clean_df = orig_df.copy()
    else:
        expanded_rows = []
        for _, row in orig_df.iterrows():
            common_name = "" if pd.isna(row["Common name"]) else str(row["Common name"])
            uniprot_id = "" if pd.isna(row["Uniprot ID"]) else str(row["Uniprot ID"])
            common_names = [x.strip() for x in common_name.split() if x.strip()]
            uniprot_ids = [x.strip() for x in uniprot_id.split() if x.strip()]
            pairs = list(zip(common_names, uniprot_ids)) if len(common_names) == len(uniprot_ids) and common_names else [(common_name, uniprot_id)]
            for common, uniprot in pairs:
                item = row.copy()
                item["Common name"] = common
                item["Uniprot ID"] = uniprot
                expanded_rows.append(item)
        clean_df = pd.DataFrame(expanded_rows)

    clean_df.to_csv(CLEAN_OUT, index=False, encoding="utf-8-sig")

    if "Common name" in clean_df.columns:
        unique_targets = sorted(clean_df["Common name"].dropna().astype(str).str.strip().replace("", pd.NA).dropna().unique())
    else:
        unique_targets = []
    pd.DataFrame({"targets": unique_targets}).to_csv(UNIQUE_OUT, index=False, encoding="utf-8-sig")

    print(f"[done] {ORIG_OUT.name}: {len(orig_df)} rows")
    print(f"[done] {CLEAN_OUT.name}: {len(clean_df)} rows")
    print(f"[done] {UNIQUE_OUT.name}: {len(unique_targets)} unique targets")


def main():
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    input_df = pd.read_csv(INPUT, encoding="utf-8-sig")
    input_df = input_df[input_df["SMILES"].notna() & (input_df["SMILES"].astype(str).str.strip() != "")].copy()
    print(f"[input] {len(input_df)} metabolites from {INPUT}")

    driver = make_driver()
    try:
        for index, row in input_df.iterrows():
            predict_one(driver, row, index)
            time.sleep(2)
    finally:
        driver.quit()

    rebuild_metabolite_target_files(input_df)


if __name__ == "__main__":
    main()
