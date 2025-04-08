make test:
	make test_tool tool=visualization/jaspar2logo profile=test
	make test_tool tool=CRISPR/ICE profile=test
	make test_tool tool=dorado/duplex profile=test
	make test_tool tool=dorado/simplex profile=test
	make test_tool tool=dorado/demux profile=test
	make test_tool tool=strainline/strainline profile=test
	make test_tool tool=strainline/clipqs profile=test
	make test_tool tool=phylo/FastTree profile=test
	make test_tool tool=msa/muscle profile=test
	make test_tool tool=phylo/phytreeviz profile=test
	make test_tool tool=phylo/reroot profile=test
	make test_tool tool=pod5/convert_fast5 profile=test
	make test_tool tool=pod5/split_by_channel profile=test
	make test_tool tool=picard/addorreplacereadgroups profile=test
	make test_tool tool=seqkit/primercheck profile=test
	make test_tool tool=hiv/intactness profile=test
	make test_tool tool=huggingface/hiv-bert profile=test
	
test_tool:
	snakemake --profile profiles/$(profile) -d $(tool)/test/ -s $(tool)/test/Snakefile

clean:
	find . -type d -name ".snakemake" -exec rm -rf {} +
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name "scratch" -exec rm -rf {} +
	find . -type d -name "*.log" -exec rm -rf {} +

very_clean:
	make clean
	find . -type d -name "venv" -exec rm -rf {} +

venv:
	mamba env create -f environment.yaml --prefix ./venv

strainline/venv:
	cd strainline
	make venv
