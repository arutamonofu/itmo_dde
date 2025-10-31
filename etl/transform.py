# /etl/transform.py

import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem import Descriptors
import argparse


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
    df.dropna(subset=["pIC50"], inplace=True)
    
    # 4. Удаление дубликатов по SMILES, оставляем запись с наибольшей активностью
    df_clean = df[["canonical_smiles", "pIC50"]].copy()
    df_clean = df_clean.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles", keep="first"
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

    df = df.rename(columns={"SMILES": "canonical_smiles"}) # Переименуем сразу для унификации

    # 1. Фильтрация IC50 и активности
    df = df[df["Activity Name"] == "IC50"].copy()
    df = df[df["Activity Outcome"] == "Active"].copy()

    # 2. Валидация SMILES
    df = df.dropna(subset=["canonical_smiles"])
    df["is_valid"] = df["canonical_smiles"].apply(is_valid_smiles)
    df = df[df["is_valid"]].drop(columns=["is_valid"])

    # 3. Преобразование в pIC50
    df["Activity Value [uM]"] = pd.to_numeric(df["Activity Value [uM]"], errors="coerce")
    df = df.dropna(subset=["Activity Value [uM]", "canonical_smiles"])
    activity_uM = pd.to_numeric(df["Activity Value [uM]"], errors="coerce")
    safe_activity_uM = activity_uM.where(activity_uM > 0, 1e-6)
    df["pIC50"] = 6 - np.log10(safe_activity_uM)
    
    # 4. Удаление дубликатов
    df_clean = df[["canonical_smiles", "pIC50"]].copy()
    df_clean = df_clean.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles", keep="first"
    )

    return df_clean

def combine_preprocessed_datasets(
    df_chembl: pd.DataFrame,
    df_pubchem: pd.DataFrame
) -> pd.DataFrame:
    """Объединяет очищенные датасеты из ChEMBL и PubChem в один."""
    
    # Объединение
    df_combined = pd.concat([df_chembl, df_pubchem], ignore_index=True)

    # Удаление дубликатов
    df_combined = df_combined.sort_values("pIC50", ascending=False).drop_duplicates(
        "canonical_smiles"
    )

    return df_combined

def calculate_rdkit_descriptors(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """Вычисляет набор дескрипторов RDKit для каждой SMILES-строки в DataFrame."""
    
    mols = [Chem.MolFromSmiles(smi) for smi in df[smiles_col]]
    desc_func_list = {
        "MolWt": Descriptors.MolWt, "MolLogP": Descriptors.MolLogP,
        "NumHDonors": Descriptors.NumHDonors, "NumHAcceptors": Descriptors.NumHAcceptors,
        "NumRotatableBonds": Descriptors.NumRotatableBonds, "RingCount": Descriptors.RingCount,
        "NumAromaticRings": Descriptors.NumAromaticRings, "TPSA": Descriptors.TPSA,
    }
    desc_list = []
    for mol in tqdm(mols, desc="Генерация дескрипторов RDKit"):
        if mol:
            desc_list.append({name: func(mol) for name, func in desc_func_list.items()})
        else:
            desc_list.append({name: None for name in desc_func_list.keys()})
    df_desc = pd.DataFrame(desc_list, index=df.index)
    return pd.concat([df, df_desc], axis=1).dropna()

def convert_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Приводит типы данных в DataFrame к более оптимальным для экономии памяти."""
    
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
    print("Загрузка конфигурации...")

    parser = argparse.ArgumentParser()
    parser.add_argument("config_path", help="Путь к файлу конфигурации .yml")
    args = parser.parse_args()
    config_path = Path(args.config_path)

    try:
        # 1. Определяем корень проекта
        project_root = Path(__file__).resolve().parent.parent

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        
        # 3. Формируем абсолютный путь к директории с данными
        data_dir = project_root / config["data"]["data_dir"]
        data_file_raw_chembl = config["data"]["data_file_raw_chembl"]
        data_file_raw_pubchem = config["data"]["data_file_raw_pubchem"]
        data_file_preprocessed = config["data"]["data_file_preprocessed"]
        
    except FileNotFoundError:
        print(f"X Ошибка: файл 'config.yml' не найден по пути {config_path}")
        return
    except KeyError as e:
        print(f"X Ошибка: в 'config.yml' отсутствует необходимый ключ: {e}")
        return
    print("Загрузка конфигурации завершена\n")

    try:
        path_chembl_raw = data_dir / data_file_raw_chembl
        path_pubchem_raw = data_dir / data_file_raw_pubchem
        df_chembl_raw = pd.read_csv(path_chembl_raw)
        df_pubchem_raw = pd.read_csv(path_pubchem_raw)
    except FileNotFoundError as e:
        print(f"X Ошибка: не найден сырой файл данных: {e}. "
              "Убедитесь, что скрипт extract.py был успешно запущен.")
        return

    print("Предобработка данных с ChEMBL...")
    df_chembl_preprocessed = preprocess_chembl_data(df_chembl_raw)
    print("Предобработка данных с ChEMBL завершена\n")

    print("Предобработка данных с PubChem ---")
    df_pubchem_preprocessed = preprocess_pubchem_data(df_pubchem_raw)
    print("Предобработка данных с PubChem завершена\n")

    print("Объединение датасетов ChEMBL и PubChem...")
    df_combined = combine_preprocessed_datasets(
        df_chembl_preprocessed, df_pubchem_preprocessed
    )
    print("Объединение датасетов ChEMBL и PubChem завершено\n")

    print("Расчёт молекулярных дескрипторов...")
    df_with_features = calculate_rdkit_descriptors(
        df_combined, "canonical_smiles"
    )
    print("Расчёт молекулярных дескрипторов завершён\n")

    df_final = convert_dtypes(df_with_features)
    
    print("Сохранение предобработанного датасета...")
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / data_file_preprocessed
    df_final.to_parquet(output_path, index=False)
    print("Сохранение предобработанного датасета завершено")


if __name__ == "__main__":
    main()