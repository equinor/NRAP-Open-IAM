# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = openiam
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	# Compile user's guide and copy to website source
	cd ../user; \
	make latexpdf; \
	make html; \
	cp -f build/latex/OpenIAM_user.pdf ../website/source/doc; \
	# Compile developer's guide and copy to website source
	cd ../developer; \
	make latexpdf; \
	make html; \
	cp -f build/latex/OpenIAM_dev.pdf ../website/source/doc; \
	# Compile QAQC docs and copy to website source
	cd ../qaqc; \
	make latexpdf; \
	make html; \
	cp -f build/latex/OpenIAM_QAQC.pdf ../website/source/doc; \
	# Compile website
	cd ../website;
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	# Copy users, developers, and qaqc folders to website folder
	cp -rf ../user/build/html build/html/user_html; \
	cp -rf ../developer/build/html build/html/dev_html; \
	cp -rf ../qaqc/build/html build/html/qaqc_html; \

clean:
	rm -rf build/*
	rm -rf ../user/build/*
	rm -rf ../developer/build/*
	rm -rf ../qaqc/build/*
