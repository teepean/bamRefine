.PHONY: build

install:
	bash install_dependencies.sh
build:
	python3 setup.py build_ext --inplace

