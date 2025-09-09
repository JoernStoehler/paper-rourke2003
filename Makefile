PY := python3
PIP := pip3

install:
	$(PIP) install -r requirements.txt

build:
	$(PY) scripts/build_cffi.py

test: build
	PYTHONPATH=src pytest -q

jupyter: build
	jupyter lab --ip=0.0.0.0 --no-browser

clean:
	rm -rf build/
	rm -f src/optproj/_cext.*
	find . -name '__pycache__' -type d -exec rm -rf {} +
	find . -name '*.pyc' -delete
