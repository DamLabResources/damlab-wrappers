PREFIX = $(PWD)/venv
ENV_FILE = env.yaml

# Default target
all: $(PREFIX)

# Create the environment on a GPU node
$(PREFIX):
	mkdir -p $(PREFIX)
	srun --partition=gpu --gres=gpu:1 --time=01:00:00 --mem=16G \
		mamba env create -p $(PREFIX) -f $(ENV_FILE)

# Run tests
test: $(PREFIX)
	srun --partition=gpu --gres=gpu:1 --time=01:00:00 --mem=16G \
		python -m pytest tests.py

# Clean up
clean:
	rm -rf $(PREFIX)
	rm -rf .build

# Rebuild the environment
rebuild: clean all

.PHONY: all clean rebuild test 