# This Makefile is for convenience as a reminder and shortcut for the most used commands

# Check if SageMath is installed.
ifeq ($(shell which sage >/dev/null 2>&1; echo $$?), 1)
$(error The 'sage' command was not found. Make sure you have SageMath installed (cf. https://www.sagemath.org/).)
endif

# Package folder
PACKAGE = drinfeld_modular_forms

# change to your sage command if needed
SAGE = sage

all: install

install:
	$(SAGE) -pip install --upgrade --no-index -v .

uninstall:
	$(SAGE) -pip uninstall $(PACKAGE)

develop:
	$(SAGE) -pip install --upgrade -e .

coverage:
	$(SAGE) -coverage $(PACKAGE)/*

doc:
	cd docs && $(SAGE) -sh -c "make html"

doc-pdf:
	cd docs && $(SAGE) -sh -c "make latexpdf"

clean: clean-doc

clean-doc:
	cd docs && $(SAGE) -sh -c "make clean"

.PHONY: all install develop coverage clean clean-doc doc doc-pdf
