# /etl/extract.py

import yaml
import requests
import time
import pandas as pd
from pathlib import Path
import argparse

def download_chembl_data(target_chembl_id: str) -> pd.DataFrame | None:
    """
    Скачивает сырые данные по биологической активности для заданной мишени из базы данных ChEMBL.
    Собирает данные по IC50 для малых молекул с использованием пагинации.

    Args:
        target_chembl_id (str): Идентификатор мишени в ChEMBL (например, "CHEMBL230").

    Returns:
        pd.DataFrame | None: DataFrame с сырыми данными или None в случае ошибки.
    """

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    PARAMS = {
        "target_chembl_id": target_chembl_id,
        "standard_type": "IC50",
        "limit": 1000,
    }

    all_activities = []
    page = 0
    
    while True:
        PARAMS["offset"] = PARAMS["limit"] * page
        try:
            response = requests.get(BASE_URL, params=PARAMS, timeout=30)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429:
                print("Слишком много запросов. Ожидание 10 секунд...")
                time.sleep(10)
                continue
            else:
                print(f"X Ошибка сети при запросе: {e}")
                return None

        data = response.json()
        activities = data.get("activities")

        if not activities:
            break

        all_activities.extend(activities)
        page += 1
        time.sleep(1)

    if not all_activities:
        return None

    return pd.DataFrame(all_activities)


def download_pubchem_data(target_pubchem_id: str) -> pd.DataFrame | None:
    """
    Скачивает сырые данные по BioAssay для заданной мишени из PubChem,
    получает для них SMILES и объединяет в один файл.

    Args:
        target_pubchem_id (str): Идентификатор белка (например, "P04058").

    Returns:
        pd.DataFrame | None: DataFrame с сырыми данными или None в случае ошибки.
    """
    
    try:
        bioassay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{target_pubchem_id}/concise/json"
        response = requests.get(bioassay_url, timeout=60)
        response.raise_for_status()
        concise_data = response.json()

        df_bioassay = pd.DataFrame(
            data=[row["Cell"] for row in concise_data["Table"]["Row"]],
            columns=concise_data["Table"]["Columns"]["Column"],
        )

        cids_string = ",".join(df_bioassay["CID"].astype(str).unique())
        smiles_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/SMILES/json"
        post_data = {"cid": cids_string}

        response = requests.post(smiles_url, data=post_data, timeout=120)
        response.raise_for_status()
        smiles_data = response.json()
        df_smiles = pd.DataFrame(smiles_data["PropertyTable"]["Properties"])

        df_bioassay["CID"] = pd.to_numeric(df_bioassay["CID"], errors="coerce")
        df_bioassay.dropna(subset=["CID"], inplace=True)
        df_bioassay["CID"] = df_bioassay["CID"].astype(int)
        df_merged = pd.merge(df_bioassay, df_smiles, how="left", on="CID")
        
        return df_merged

    except requests.exceptions.RequestException as e:
        print(f"X Ошибка при загрузке данных из PubChem: {e}")
        return None
    except KeyError as e:
        print(f"X Ошибка: неверный формат ответа от API. Не найден ключ: {e}")
        return None


def main():
    print("Загрузка конфигурации...")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_path", help="Путь к файлу конфигурации .yml")
    args = parser.parse_args()
    config_path = Path(args.config_path)

    try:
        project_root = Path(__file__).resolve().parent.parent

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        
        target_chembl_id = config["data"]["target_chembl"]
        target_pubchem_id = config["data"]["target_pubchem"]
        
        data_dir = project_root / config["data"]["data_dir"]
        data_file_raw_chembl = config["data"]["data_file_raw_chembl"]
        data_file_raw_pubchem = config["data"]["data_file_raw_pubchem"]
        
    except FileNotFoundError:
        print(f"X Ошибка: файл 'config.yml' не найден по пути {config_path}")
        return
    except KeyError as e:
        print(f"X Ошибка: в 'config.yml' отсутствует необходимый ключ: {e}")
        return
        
    print("Загрузка конфигурации завершена\n")
    data_dir.mkdir(parents=True, exist_ok=True)

    print("Скачивание данных с ChEMBL...")
    df_chembl_raw = download_chembl_data(target_chembl_id)
    if df_chembl_raw is not None:
        output_path = data_dir / data_file_raw_chembl
        df_chembl_raw.to_csv(output_path, index=False)
        print(f"Сырые данные ChEMBL сохранены в: {output_path}\n")
    else:
        print("Не удалось загрузить данные из ChEMBL. Пропускаем сохранение.\n")

    print("Скачивание данных с PubChem...")
    df_pubchem_raw = download_pubchem_data(target_pubchem_id)
    if df_pubchem_raw is not None:
        output_path = data_dir / data_file_raw_pubchem
        df_pubchem_raw.to_csv(output_path, index=False)
        print(f"Сырые данные PubChem сохранены в: {output_path}\n")
    else:
        print("Не удалось загрузить данные из PubChem. Пропускаем сохранение.\n")


if __name__ == "__main__":
    main()