PREFIX = $(PWD)/venv

install: $(PREFIX)

$(PREFIX):
	mamba env create -p $(PREFIX) -f environment.yaml

tests:
	snakemake -d test -c 1 --use-conda -s test/Snakefile $(EXTRA)
