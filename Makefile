VENV_DIR := .venv
PY := $(VENV_DIR)/bin/python3
PIP := $(VENV_DIR)/bin/pip3

.PHONY: all install build test jupyter clean

all: install build

$(VENV_DIR):
	python3 -m venv $(VENV_DIR)

install: $(VENV_DIR)
	$(PIP) install -r requirements.txt

build: $(VENV_DIR)
	$(PY) scripts/build_cffi.py

test: build
	PYTHONPATH=src $(PY) -m pytest -q

jupyter: build
	$(PY) -m jupyter lab --ip=0.0.0.0 --no-browser

clean:
	rm -rf build/
	rm -f src/optproj/_cext.*
	find . -name '__pycache__' -type d -exec rm -rf {} +
	find . -name '*.pyc' -delete