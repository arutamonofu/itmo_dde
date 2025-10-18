import subprocess
import sys
import os

def main():
    # Подготовка данных
    print("Запуск скрипта подготовки данных...")
    data_loader_path = os.path.join("scripts", "data_loader.py")
    subprocess.run([sys.executable, data_loader_path], check=True)
    print("Скрипт подготовка данных выполнен")
    print("Запуск скрипта загрузки данных в базу")
    db_writer_path = os.path.join("scripts", "write_to_db.py")
    subprocess.run([sys.executable, db_writer_path], check=True)
    print("Скрипт загрузки данных в базу выполнен")

if __name__ == "__main__":
    main()