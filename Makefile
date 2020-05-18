release/major release/minor release/patch release/rc:
	bumpversion $(@F)
	git push
	git push --tags

testenv:
	pip install -r requirements/testing.txt
	pip install -e .

test:
	py.test --cov=pythomics --cov-report=xml -s tests

