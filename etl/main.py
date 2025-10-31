import subprocess
import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Запуск ETL-пайплайна для обработки данных по молекулам")
    parser.add_argument(
        '--config',
        type=str,
        default='config.yml',
        help='Путь к файлу конфигурации .yml (по умолчанию: config.yml)'
    )
    args = parser.parse_args()
    config_path = args.config

    print(f"Используется файл конфигурации: {config_path}\n")

    print("Запуск скрипта извлечения...")
    extract_script_path = os.path.join("etl", "extract.py")
    try:
        subprocess.run([sys.executable, extract_script_path, config_path], check=True)
        print("Скрипт извлечения выполнен")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка выполнения extract.py: {e}")
        sys.exit(1)

    print("\n\nЗапуск скрипта преобразования...")
    transform_script_path = os.path.join("etl", "transform.py")
    try:
        subprocess.run([sys.executable, transform_script_path, config_path], check=True)
        print("Скрипт преобразования выполнен")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка выполнения transform.py: {e}")
        sys.exit(1)

    print("\n\nЗапуск скрипта загрузки данных...")
    load_script_path = os.path.join("etl", "load.py")
    try:
        subprocess.run([sys.executable, load_script_path, config_path], check=True)
        print("Скрипт загрузки выполнен")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка выполнения load.py: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()