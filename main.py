import subprocess
import sys
import os

def main():
    # Подготовка данных
    print("Запуск скрипта подготовки данных...")
    data_loader_path = os.path.join("scripts", "data_loader.py")
    subprocess.run([sys.executable, data_loader_path], check=True)
    print("Скрипт подготовка данных выполнен")

if __name__ == "__main__":
    main()