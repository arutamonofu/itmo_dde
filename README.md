# Учебный проект в рамках дисциплины «Инжиниринг управления данными» ИТМО

Данные — молекулы, для которых будет необходимо сгенерировать и обработать численные признаки (молекулярные дескрипторы).

Источники данных — химические базы данных (выгрузка по API):

- ChEMBL: https://www.ebi.ac.uk/chembl/ ([документация API](https://www.ebi.ac.uk/chembl/api/data/docs)):
  - Пример запроса (10765 строк): https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL4822&standard_type=IC50&limit=1000
- PubChem: https://pubchem.ncbi.nlm.nih.gov/ ([документация API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial))
