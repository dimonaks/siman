rm dist/*
python3.8 setup.py sdist 
python3.8 -m twine upload dist/*
