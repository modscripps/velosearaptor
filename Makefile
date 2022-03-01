.PHONY: test docs servedocs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

test: ## run tests quickly with the default Python
	pytest -W ignore::DeprecationWarning

docs: ## generate documentation using pdoc
	rm -rf docs
	pdoc -d numpy -o docs gadcp
	$(BROWSER) docs/index.html

servedocs: ## compile the docs & watch for changes
	pdoc -d numpy gadcp
	# $(BROWSER) http://localhost:8080
