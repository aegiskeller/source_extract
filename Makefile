install:
	pip3 install --upgrade pip &&\
		pip3 install -r requirements.txt

install-astrometry:
	./install_dependencies.sh

install-indexfiles:
	./install_indexfiles.sh

test:
	python -m pytest -vv tests/test_main.py

format:
	black .

lint:
	pylint --disable=R,C *.py

all: install lint format test