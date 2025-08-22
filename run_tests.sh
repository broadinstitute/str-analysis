python3 -m pip install --upgrade coverage
#python3 -m unittest discover -p "*tests.py"
coverage run --omit="setup.py,*/site-packages/*,__init__.py" -m unittest discover -p "*tests.py"
coverage html
