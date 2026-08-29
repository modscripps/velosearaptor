.PHONY: check format format-check test docs ghdocs servedocs help
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

BROWSER := uv run python -c "$$BROWSER_PYSCRIPT"

help:
	@uv run python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

check: ## check style
	uv run ruff check src/velosearaptor/ tests/ notebooks/

format: ## format code using ruff
	uv run ruff format src/velosearaptor/ tests/ notebooks/

format-check: ## check code style using ruff format --diff
	uv run ruff format --diff src/velosearaptor/
	uv run ruff format --diff tests/
	uv run ruff format --diff notebooks/

LOGO := https://github.com/modscripps/velosearaptor/raw/main/logo/velosearaptor.png
PDOC_OPTS := -d numpy -t .pdoc-theme-gv --math --logo $(LOGO)

test: ## run tests quickly with the default Python
	uv run pytest

docs: ## generate documentation using pdoc
	rm -rf docs
	uv run pdoc $(PDOC_OPTS) -o docs src/velosearaptor/
	$(BROWSER) docs/index.html

ghdocs: ## generate documentation for GitHub Pages
	rm -rf docs
	PDOC_ALLOW_EXEC=1 pdoc $(PDOC_OPTS) -o docs src/velosearaptor/

servedocs: ## compile the docs & watch for changes
	uv run pdoc $(PDOC_OPTS) src/velosearaptor/
