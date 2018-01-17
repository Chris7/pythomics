testenv:
	pip install -r requirements.txt
	pip install -e .

test:
	nosetests --with-coverage --cover-erase --cover-package=pythomics -s tests

