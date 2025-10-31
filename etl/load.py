# /etl/load.py

import os
import pandas as pd
import yaml
from pathlib import Path
from dotenv import load_dotenv
from sqlalchemy import (create_engine, inspect, text, Table, Column, Integer, 
                        Float, String, MetaData)
import argparse

def get_db_engine():
    db_url = os.getenv('DB_URL')
    db_port = os.getenv('DB_PORT')
    db_user = os.getenv('DB_USER')
    db_password = os.getenv('DB_PASSWORD')
    db_name = os.getenv('DB_NAME')

    connection_url = f"postgresql+psycopg2://{db_user}:{db_password}@{db_url}:{db_port}/{db_name}"
    return create_engine(connection_url)

def upload_data(engine, table_name: str, data_path: Path, sample_size: int):
    data = pd.read_parquet(data_path)
    data = data.head(sample_size)
    data.columns = [col.lower() for col in data.columns]

    metadata = MetaData()
    molecules_table = Table(
        table_name,
        metadata,
        Column('id', Integer, primary_key=True),  # Автоинкрементный ID
        Column('canonical_smiles', String, nullable=False),
        Column('pic50', Float),
        Column('molwt', Float),
        Column('mollogp', Float),
        Column('numhdonors', Integer),
        Column('numhacceptors', Integer),
        Column('numrotatablebonds', Integer),
        Column('ringcount', Integer),
        Column('numaromaticrings', Integer),
        Column('tpsa', Float),
        schema='public'
    )

    with engine.connect() as connection:
        connection.execute(text(f'DROP TABLE IF EXISTS public."{table_name}" CASCADE'))
        connection.commit()

    metadata.create_all(engine)


    data.to_sql(
        table_name,
        engine,
        schema='public',
        if_exists='append',
        index=False,
        method='multi'
    )

def inspect_table(engine, table_name: str):
    inspector = inspect(engine)
    if inspector.has_table(table_name, schema='public'):
        print(f"Таблица '{table_name}' найдена. Структура:")
        columns = inspector.get_columns(table_name, schema='public')
        for col in columns:
            print(f"- {col['name']}: {col['type']}")
        try:
            with engine.connect() as connection:
                query = text(f'SELECT COUNT(*) FROM public."{table_name}"')              
                result = connection.execute(query)
                row_count = result.scalar_one() 
            print(f"\nСтрок в таблице: {row_count}")
        except Exception as e:
            print(f"\nX Ошибка при подсчете строк в таблице: {e}")
    else:
        print(f"Таблица '{table_name}' не найдена!")

def main():
    print("Загрузка конфигурации...")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_path", help="Путь к файлу конфигурации .yml")
    args = parser.parse_args()
    config_path = Path(args.config_path)
    project_root = Path(__file__).resolve().parent.parent

    dotenv_path = project_root / '.env'
    load_dotenv(dotenv_path=dotenv_path)
    
    table_name = os.getenv('DB_TABLE')

    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        data_dir = config["data"]["data_dir"]
        data_file_preprocessed = config["data"]["data_file_preprocessed"]
        sample_size = config["db"]["sample_size"]

        data_path = project_root / data_dir / data_file_preprocessed

    except FileNotFoundError:
        print(f"X Ошибка: Файл 'config.yml' не найден по пути {config_path}.")
        return
    except (KeyError, ValueError) as e:
        print(f"X Ошибка конфигурации: {e}")
        return
    
    print("Конфигурация успешно загружена.\n")

    try:
        print("Подключение к базе данных...")
        engine = get_db_engine()
        print("Подключение успешно установлено\n")

        print("Загрузка данных в базу...")
        upload_data(engine, table_name, data_path, sample_size)
        print("Загрузка данных в базу завершена\n")

        print("Инспектирование загруженной таблицы...")
        inspect_table(engine, table_name)

    except Exception as e:
        print(f"\nX Произошла ошибка во время работы с базой данных: {e}")


if __name__ == '__main__':
    main()