# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext

.PHONY: default
default: cclib update attributes coverage html copy
copy:
	cp -rvf _build/html/* ../.

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"

clean:
	-rm -rf $(BUILDDIR)/*
	-rm -rf cclib*
	-rm -rf attributes*.rst
	-rm -rf converage*.rst

PRODUCTION_VERSION = v1.5.3
CCLIB_PROD = $(BUILDDIR)/cclib_prod
CCLIB_DEV = $(BUILDDIR)/cclib_dev
SRC_PROD = $(CCLIB_PROD)/cclib
SRC_DEV = $(CCLIB_DEV)/cclib

# Have the production and master in two separate directories so that
# dependencies are triggered correctly. If we had one, it would not
# work correctly, since git updates timestamps on each checkout, which
# would trigger some targets unnecessarily.
.PHONY: cclib
cclib: $(CCLIB_PROD) $(CCLIB_DEV)
$(CCLIB_PROD):
	@echo "Checking out prod..."
	git clone git@github.com:cclib/cclib.git $@
	cd $(CCLIB_PROD); git checkout $(PRODUCTION_VERSION)
$(CCLIB_DEV):
	@echo "Checking out dev..."
	git clone -b master git@github.com:cclib/cclib.git $@

.PHONY: update update_prod update_dev
update: update_prod update_dev
update_prod:
	@echo "Updating prod..."
	cd $(CCLIB_PROD); git checkout $(PRODUCTION_VERSION); git fetch --all --tags
	cd $(CCLIB_PROD); find . -name "*.pyc" -o -name "__pycache__" -delete
update_dev:
	@echo "Updating dev..."
	cd $(CCLIB_DEV); git pull origin; git fetch --all --tags
	cd $(CCLIB_DEV); find . -name "*.pyc" -o -name "__pycache__" -delete

# Since we have two separate directories, 'checking out' a branch actually
# means making a symlink to the appropriate directory, because scripts should be
# impervious to this Makefile and will import from always the same name (cclib).
.PHONY: checkout_prod checkout_dev
checkout_prod:
	-rm -rf cclib
	ln -s $(SRC_PROD) cclib
checkout_dev:
	-rm -rf cclib
	ln -s $(SRC_DEV) cclib

# We need three layers of targets in order to always get the checkout
# to execute, but the actual file-generating target to run only if the
# appropriate files have been updated (resolved via the right symlink).
.PHONY: attributes
attributes: attributes.rst attributes_dev.rst
attributes.rst: attributes.py $(SRC_PROD)/parser/data.py
	@echo "Generating prod attributes..."
	@$(MAKE) --no-print-directory checkout_prod
	python attributes.py > $@
attributes_dev.rst: attributes.py $(SRC_DEV)/parser/data.py
	@echo "Generating dev attributes..."
	@$(MAKE) --no-print-directory checkout_dev
	python attributes.py > $@

# Same as above, we need three layers of targets, since this script
# especially takes some time as it performs all the tests.
.PHONY: coverage
coverage: coverage.rst coverage_dev.rst
coverage.rst: coverage.py $(wildcard $(SRC_PROD)/parser/*parser.py) $(wildcard $(CCLIB_PROD)/test/*.py)
	@echo "Generating prod coverage..."
	@$(MAKE) --no-print-directory checkout_prod
	python coverage.py $(CCLIB_PROD) > $@
coverage_dev.rst: coverage.py $(wildcard $(SRC_DEV)/parser/*parser.py) $(wildcard $(CCLIB_DEV)/test/*.py)
	@echo "Generating dev coverage..."
	@$(MAKE) --no-print-directory checkout_dev
	python coverage.py $(CCLIB_DEV) > $@

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/cclib.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/cclib.qhc"

devhelp:
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/cclib"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/cclib"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."
