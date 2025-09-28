import yaml
import numpy as np
import requests
from pathlib import Path
import time
import pandas as pd
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem import Descriptors


def download_chembl_data(target_chembl_id: str) -> pd.DataFrame:
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
        "limit": 1000,
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
            if (
                isinstance(e, requests.exceptions.HTTPError)
                and e.response.status_code == 429
            ):
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
        time.sleep(10)  # Небольшая задержка между запросами

    # Ничего не скачалось
    if not all_activities:
        return

    df = pd.DataFrame(all_activities)

    return df


def download_pubchem_data(target_pubchem_id: str) -> pd.DataFrame:
    """
    Скачивает данные по BioAssay для заданной мишени из PubChem,
    получает для них SMILES и объединяет в один файл.

    Args:
        target_accession (str): Идентификатор белка (например, "P04058").
        output_path (str | Path): Путь для сохранения итогового CSV файла.
    """
    try:
        # 1. Загрузка всех BioAssay по мишени
        bioassay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{target_pubchem_id}/concise/json"
        response = requests.get(bioassay_url, timeout=60)
        response.raise_for_status()
        concise_data = response.json()

        df_bioassay = pd.DataFrame(
            data=[row["Cell"] for row in concise_data["Table"]["Row"]],
            columns=concise_data["Table"]["Columns"]["Column"],
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
        df_merged = pd.merge(df_bioassay, df_smiles, how="left", on="CID")

        # # Сохранение
        # output_path = output_dir + f"pubchem_{target_pubchem_id}_raw.csv"
        # df_merged.to_csv(output_path, index=False)

        return df_merged

    except requests.exceptions.RequestException as e:
        print(f"X Ошибка при загрузке данных из PubChem: {e}")
    except KeyError as e:
        print(f"X Ошибка: неверный формат ответа от API. Не найден ключ: {e}")


def is_valid_smiles(smiles: str) -> bool:
    """Проверяет, является ли строка SMILES валидной."""
    if not isinstance(smiles, str) or pd.isna(smiles):
        return False
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except (
        ValueError,
        Chem.rdchem.AtomValenceException,
        Chem.rdchem.KekulizeException,
    ):
        return False


def preprocess_chembl_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Очищает и преобразует сырые данные из ChEMBL.
    - Фильтрует по типу активности IC50.
    - Вычисляет pIC50.
    - Валидирует SMILES.
    - Удаляет дубликаты.
    """

    # 1. Фильтрация IC50
    df = df[df["standard_type"] == "IC50"].copy()
    df = df.dropna(subset=["standard_value", "canonical_smiles"])

    # 2. Валидация SMILES
    df["is_valid"] = df["canonical_smiles"].apply(is_valid_smiles)
    df = df[df["is_valid"]].drop(columns=["is_valid"])

    # 3. Преобразование в pIC50
    sv_numeric_nm = pd.to_numeric(df["standard_value"], errors="coerce")
    sv_molar = sv_numeric_nm * 1e-9
    safe_sv_molar = sv_molar.where(sv_molar > 0, 1e-12)
    pIC50_calculated = -np.log10(safe_sv_molar)
    df["pIC50"] = df["pchembl_value"].fillna(pIC50_calculated)
    df["pIC50"] = pd.to_numeric(df["pIC50"], errors="coerce")

    # 4. Удаление дубликатов по SMILES, оставляем запись с наибольшей активностью
    df_clean = df[["canonical_smiles", "pIC50"]].copy()
    df_clean = df_clean.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles"
    )

    return df_clean


def preprocess_pubchem_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Очищает и преобразует сырые данные из PubChem.
    - Фильтрует по типу активности IC50 и статусу 'Active'.
    - Вычисляет pIC50 из 'Activity Value [uM]'.
    - Валидирует SMILES.
    - Удаляет дубликаты.
    """

    # 1. Фильтрация IC50 и активности
    df = df[df["Activity Name"] == "IC50"].copy()
    df = df[df["Activity Outcome"] == "Active"].copy()

    # 2. Валидация SMILES
    df["is_valid"] = df["SMILES"].apply(is_valid_smiles)
    df = df[df["is_valid"]].drop(columns=["is_valid"])

    # 3. Преобразование в pIC50
    df["Activity Value [uM]"] = pd.to_numeric(
        df["Activity Value [uM]"], errors="coerce"
    )
    df = df.dropna(subset=["Activity Value [uM]", "SMILES"])
    activity_uM = pd.to_numeric(df["Activity Value [uM]"], errors="coerce")
    safe_activity_uM = activity_uM.where(activity_uM > 0, 1e-6)
    df["pIC50"] = 6 - np.log10(safe_activity_uM)

    # 4. Удаление дубликатов
    df = df.rename(columns={"SMILES": "canonical_smiles"})
    df_clean = df[["canonical_smiles", "pIC50"]].copy()
    df_clean = df_clean.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles"
    )

    return df_clean


def combine_preprocessed_datasets(
    df_chembl: pd.DataFrame,
    df_pubchem: pd.DataFrame,
) -> pd.DataFrame:
    """
    Объединяет очищенные датасеты из ChEMBL и PubChem в один.
    """
    # Объединение
    df_combined = pd.concat([df_chembl, df_pubchem], ignore_index=True)

    # Удаление дубликатов
    df_combined = df_combined.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles"
    )

    return df_combined


def calculate_rdkit_descriptors(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """
    Вычисляет ~200 дескрипторов RDKit для каждой SMILES-строки в DataFrame.

    Args:
        df (pd.DataFrame): DataFrame, содержащий колонку со SMILES.
        smiles_col (str): Название колонки со SMILES.

    Returns:
        pd.DataFrame: Исходный DataFrame, объединенный с новыми дескрипторами.
    """

    mols = [
        Chem.MolFromSmiles(smi)
        for smi in tqdm(df[smiles_col], desc="Конвертация SMILES в Mol")
    ]

    desc_func_list = {
        "MolWt": Descriptors.MolWt,
        "MolLogP": Descriptors.MolLogP,
        "NumHDonors": Descriptors.NumHDonors,
        "NumHAcceptors": Descriptors.NumHAcceptors,
        "NumRotatableBonds": Descriptors.NumRotatableBonds,
        "RingCount": Descriptors.RingCount,
        "NumAromaticRings": Descriptors.NumAromaticRings,
        "TPSA": Descriptors.TPSA,
    }

    desc_list = []
    for mol in tqdm(mols, desc="Генерация дескрипторов RDKit"):
        if mol is not None:
            descriptors = {name: func(mol) for name, func in desc_func_list.items()}
            desc_list.append(descriptors)
        else:
            desc_list.append({name: None for name in desc_func_list.keys()})
    df_desc = pd.DataFrame(desc_list, index=df.index)

    return pd.concat([df, df_desc], axis=1)


def convert_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    optimal_dtypes = {
        "MolWt": "float32",
        "MolLogP": "float32",
        "TPSA": "float32",
        "NumHDonors": "uint8",
        "NumHAcceptors": "uint8",
        "NumRotatableBonds": "uint8",
        "RingCount": "uint8",
        "NumAromaticRings": "uint8",
    }
    df = df.astype(optimal_dtypes)

    return df


def main():
    # Загрузка конфигурации

    print("Загрузка конфигурации...")
    with open("config.yml", "r") as f:
        config = yaml.safe_load(f)
    target_chembl_id = config["data"]["target_chembl"]
    target_pubchem_id = config["data"]["target_pubchem"]
    data_dir = Path(config["data"]["data_dir"])
    print("Загрузка конфигурации завершена\n")

    # Скачивание

    print("Скачивание данных с ChEMBL...")
    df_chembl_raw = download_chembl_data(target_chembl_id)
    print("Скачивание данных с ChEMBL завершено\n")

    print("Скачивание данных с PubChem...")
    df_pubchem_raw = download_pubchem_data(target_pubchem_id)
    print("Скачивание данных с PubChem завершено\n")

    # Предобработка

    print("Предобработка данных с ChEMBL...")
    df_chembl_preprocessed = preprocess_chembl_data(df_chembl_raw)
    print("Предобработка данных с ChEMBL завершена\n")

    print("Предобработка данных с PubChem...")
    df_pubchem_preprocessed = preprocess_pubchem_data(df_pubchem_raw)
    print("Предобработка данных с PubChem завершена\n")

    print("Объединение датасетов ChEMBL и PubChem...")
    df_combined_preprocessed = combine_preprocessed_datasets(
        df_chembl_preprocessed, df_pubchem_preprocessed
    )
    print("Объединение датасетов ChEMBL и PubChem завершено\n")

    # Расчёт молекулярных дескрипторов

    print("Расчёт молекулярных дескрипторов...")
    df_with_features = calculate_rdkit_descriptors(
        df_combined_preprocessed, "canonical_smiles"
    )
    print("Расчёт молекулярных дескрипторов завершён\n")

    # Приведение типов данных к оптимальным

    df_converted = convert_dtypes(df_with_features)

    # Сохранение

    print("Сохранение файла...")
    data_dir.mkdir(parents=True, exist_ok=True)
    df_converted.to_parquet(data_dir / "preprocessed.parquet", index=False)
    print("Сохранение файла завершено")


if __name__ == "__main__":
    main()