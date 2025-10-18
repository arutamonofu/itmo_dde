# /scripts/write_to_db.py

import os
import sys
import pandas as pd
import yaml
from pathlib import Path
from dotenv import load_dotenv
from sqlalchemy import create_engine, inspect

def get_engine():
    db_url = os.getenv('DB_URL')
    db_port = os.getenv('DB_PORT')
    db_user = os.getenv('DB_USER')
    db_password = os.getenv('DB_PASSWORD')
    db_root_base = os.getenv('DB_ROOT_BASE')

    url = f"postgresql+psycopg2://{db_user}:{db_password}@{db_url}:{db_port}/{db_root_base}"
    return create_engine(url)

def upload_data(engine, table_name, data_path, sample_size):
    data = pd.read_parquet(data_path)[:sample_size]
    data = data.rename({col:col.lower() for col in data.columns}, axis=1)
    data.to_sql(table_name, engine, schema='public', if_exists='replace', index=False, method='multi')

def inspect_table(engine, table_name):
    inspector = inspect(engine)
    if inspector.has_table(table_name):
        print(f"Таблица '{table_name}' найдена. Структура:")
        columns = inspector.get_columns(table_name)
        for col in columns:
            print(f"- {col['name']}: {col['type']}")
    else:
        print(f"Таблица '{table_name}' не найдена!")

def main():
    print("Загрузка конфигурации...")
    project_root = Path(__file__).resolve().parent.parent

    dotenv_path = project_root / '.env'
    load_dotenv(dotenv_path=dotenv_path)

    table_name = os.getenv('TABLE_NAME')

    config_path = project_root / 'config.yml'
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    data_dir = config.get("data", {}).get("data_dir")
    data_file = config.get("data", {}).get("data_file")
    sample_size = config.get("db", {}).get("sample_size")
    
    data_path = project_root / data_dir / data_file
    print("Загрузка конфигурации завершена\n")

    print("Подключение к базе данных...")
    engine = get_engine()
    print("Подключение к базе данных завершено\n")

    print("Загрузка данных в базу...")
    upload_data(engine, table_name, data_path, sample_size)
    print("Загрузка данных в базу завершена\n")

    print("Инспектирование загруженной таблицы...")
    inspect_table(engine, table_name)


if __name__ == '__main__':
    main()