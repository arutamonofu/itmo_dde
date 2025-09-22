import requests
from pathlib import Path
import time
import pandas as pd

CHEMBL_TARGET_ID = "CHEMBL4822"
PUBCHEM_ACCESSION = "P56817"

def download_chembl_data(target_chembl_id: str, output_path: str | Path):
    """
    Скачивает данные по биологической активности для заданной мишени из базы данных ChEMBL.
    Собирает данные по IC50 для малых молекул с использованием пагинации.

    Args:
        target_chembl_id (str): Идентификатор мишени в ChEMBL (например, "CHEMBL230").
        output_path (str | Path): Путь для сохранения итогового CSV файла.
    """    
    # Запрос
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    PARAMS = {
        "target_chembl_id": target_chembl_id,
        "standard_type": "IC50",
        "limit": 1000
    }

    all_activities = []
    page = 0

    # Цикл загрузки пачками по limit строк
    while True:
        PARAMS["offset"] = PARAMS["limit"] * page
        
        try:
            response = requests.get(BASE_URL, params=PARAMS, timeout=30)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429:
                time.sleep(10)
                continue
            else:
                print(f"Ошибка сети при запросе: {e}")
                break

        data = response.json()
        activities = data.get("activities")

        # Если пустая страница, значит всё загружено
        if not activities:
            break
        
        all_activities.extend(activities)
        page += 1
        time.sleep(10) # Небольшая задержка между запросами

    # Ничего не скачалось
    if not all_activities:
        return

    # Сохранение
    df = pd.DataFrame(all_activities)
    df.to_csv(output_path, index=False)


def download_pubchem_data(target_accession: str, output_path: str | Path):
    """
    Скачивает данные по BioAssay для заданной мишени из PubChem,
    получает для них SMILES и объединяет в один файл.

    Args:
        target_accession (str): Идентификатор белка (например, "P04058").
        output_path (str | Path): Путь для сохранения итогового CSV файла.
    """
    try:
        # 1. Загрузка всех BioAssay по мишени
        bioassay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{target_accession}/concise/json"
        response = requests.get(bioassay_url, timeout=60)
        response.raise_for_status()
        concise_data = response.json()
        
        df_bioassay = pd.DataFrame(
            data=[row["Cell"] for row in concise_data["Table"]["Row"]],
            columns=concise_data["Table"]["Columns"]["Column"]
        )

        # 2. Загрузка SMILES для всех CID в загруженных BioAssay
        cids_string = ",".join(df_bioassay["CID"].astype(str).unique())
        smiles_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/SMILES/json"
        post_data = {"cid": cids_string}
        
        response = requests.post(smiles_url, data=post_data, timeout=120)
        response.raise_for_status()
        smiles_data = response.json()
        
        df_smiles = pd.DataFrame(smiles_data["PropertyTable"]["Properties"])

        # 3. Объединение данных
        df_bioassay["CID"] = pd.to_numeric(df_bioassay["CID"], errors="coerce")
        df_bioassay.dropna(subset=["CID"], inplace=True)
        df_bioassay["CID"] = df_bioassay["CID"].astype(int)

        df_merged = pd.merge(df_bioassay, df_smiles, how='left', on='CID')
        df_merged.to_csv(output_path, index=False)

    except requests.exceptions.RequestException as e:
        print(f"X Ошибка при загрузке данных из PubChem: {e}")
    except KeyError as e:
        print(f"X Ошибка: неверный формат ответа от API. Не найден ключ: {e}")

# Определение путей для сохранения
chembl_raw_path = f"chembl_{CHEMBL_TARGET_ID}_raw.csv"
pubchem_raw_path = f"pubchem_{PUBCHEM_ACCESSION}_raw.csv"

# Скачивание
print("Загрузка с ChEMBL...")
download_chembl_data(
    target_chembl_id=CHEMBL_TARGET_ID,
    output_path=chembl_raw_path
)
print(f"✓ Сохранено в: {chembl_raw_path}")

print("Загрузка с PubChem...")
download_pubchem_data(
    target_accession=PUBCHEM_ACCESSION,
    output_path=pubchem_raw_path
)
print(f"✓ Сохранено в: {pubchem_raw_path}")