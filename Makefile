testenv:
	pip install -r requirements/testing.txt
	pip install -e .

test:
	py.test --cov=pythomics -s tests

