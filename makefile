make test:
	make test_tool tool=visualization/jaspar2logo profile=test
	make test_tool tool=CRISPR/ICE profile=test

test_tool:
	snakemake --profile profiles/$(profile) -d $(tool)/test/ -s $(tool)/test/Snakefile
