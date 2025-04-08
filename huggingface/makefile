PREFIX = $(PWD)/venv

# Default target
all: help

# Help target
help:
	@echo "Available targets:"
	@echo "  install-picotte    - Install Picotte-specific environment with CUDA 11.6 support"
	@echo "  install-modern     - Install modern environment with latest PyTorch"
	@echo "  clean              - Remove installed environments"

# Picotte environment (CUDA 11.6)
install-picotte: $(PREFIX)/picotte

$(PREFIX)/picotte:
	mkdir -p $(PREFIX)
	mamba env create -p $(PREFIX)/picotte -f hiv-bert/environment.picotte.yaml

# Modern environment
install-modern: $(PREFIX)/modern

$(PREFIX)/modern:
	mkdir -p $(PREFIX)
	mamba env create -p $(PREFIX)/modern -f hiv-bert/environment.yaml

# Clean up
clean:
	rm -rf $(PREFIX)

.PHONY: all help install-picotte install-modern clean 