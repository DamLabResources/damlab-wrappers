PREFIX = $(PWD)/venv
ENV_FILE = environment.picotte.yaml

# Default target
all: $(PREFIX)

# Create the environment on a GPU node
$(PREFIX):
	mkdir -p $(PREFIX)
	srun --partition=gpu --gres=gpu:1 --time=01:00:00 --mem=16G \
		mamba env create -p $(PREFIX) -f $(ENV_FILE)

# Clean up
clean:
	rm -rf $(PREFIX)

# Rebuild the environment
rebuild: clean all

.PHONY: all clean rebuild 