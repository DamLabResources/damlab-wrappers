make test:
	make test_tool tool=visualization/jaspar2logo profile=test
	make test_tool tool=CRISPR/ICE profile=test
	make test_tool tool=dorado/duplex profile=test
	make test_tool tool=dorado/demux profile=test

test_tool:
	snakemake --profile profiles/$(profile) -d $(tool)/test/ -s $(tool)/test/Snakefile

clean:
	find . -type d -name ".snakemake" -exec rm -rf {} +
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
