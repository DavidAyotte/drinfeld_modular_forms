# This Makefile is for convenience as a reminder and shortcut for the most used commands

# Check if SageMath is installed.
# ifeq ($(shell which sage >/dev/null 2>&1; echo $$?), 1)
# $(error The 'sage' command was not found. Make sure you have SageMath installed (cf. https://www.sagemath.org/).)
# endif

# Package name
PACKAGE = drinfeld_modular_forms

# Package source directory
PACKAGE_SRC = src/$(PACKAGE)

# change to your sage command if needed
SAGE = sage

all: install

install:
	$(SAGE) -pip install --upgrade --no-index -v .

uninstall:
	$(SAGE) -pip uninstall $(PACKAGE)

develop:
	$(SAGE) -pip install --upgrade -e .

test:  # run the doctest
	$(SAGE) -t $(PACKAGE_SRC)

coverage:  # check the documentation coverage
	$(SAGE) -coverage $(PACKAGE_SRC)/*

doc:
	cd docs && $(SAGE) -sh -c "make html"

build-dist:  # build the distribution for PyPI package
	$(sage) -sh -c "python3 -m build"

clean: clean-doc

clean-doc:
	cd docs && rm -R html && rm -R doctrees

.PHONY: all install develop coverage clean clean-doc doc test
