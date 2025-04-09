# Snakemake Wrappers Guide

This directory contains example wrappers that demonstrate how to create Snakemake wrappers for both shell commands and Python APIs. These examples serve as templates and educational resources for creating new wrappers in the damlab-wrappers repository.

## What are Snakemake Wrappers?

Snakemake wrappers are reusable components that encapsulate the execution of bioinformatics tools within Snakemake workflows. They provide a standardized interface for running tools, handling parameters, and managing dependencies.

### Benefits of Using Wrappers

- **Reproducibility**: Wrappers ensure tools are run with consistent parameters and versions
- **Reusability**: Once created, wrappers can be used across multiple workflows
- **Maintainability**: Centralizing tool execution logic makes workflows easier to maintain
- **Documentation**: Wrappers include built-in documentation about inputs, outputs, and parameters

## Wrapper Types

This repository demonstrates two main types of wrappers:

1. **Shell Command Wrappers** (`example/shell/`): For wrapping command-line tools
2. **Python API Wrappers** (`example/api/`): For wrapping Python libraries and APIs

## Wrapper Structure

All wrappers follow a consistent directory structure:

```
wrapper_name/
├── environment.yaml       # Conda environment definition
├── README.md             # Documentation
├── wrapper.py            # Main wrapper script
├── CHANGELOG.md          # Version history and changes
└── test/                 # Test directory
    ├── env.yaml          # Test environment
    ├── Snakefile         # Test workflow
    ├── tests.py          # Test script
    └── test_data.*       # Test data files
```

### Key Components

#### environment.yaml

Defines the conda environment required to run the wrapper:

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - tool_name>=version
  - snakemake-wrapper-utils>=0.3
```

#### README.md

Documents the wrapper's functionality, inputs, outputs, and parameters:

```markdown
# Tool Name Wrapper

Description of the tool and its purpose.

## Input
* Input file formats and requirements

## Output
* Output file formats

## Parameters
* `param1` (optional, default: value)
    Description of parameter
* `param2` (required)
    Description of parameter

## Example
```python
rule example:
    input:
        "input.file"
    output:
        "output.file"
    params:
        param1="value",
        param2="value"
    wrapper:
        "file:path/to/wrapper"
```

## Author
* Author Name

## Software Requirements
* [Tool Name](https://tool-url) (tested with vX.Y.Z)
```

#### wrapper.py

The main wrapper script that executes the tool:

- For shell wrappers: Uses `snakemake.shell` to execute commands
- For API wrappers: Uses Python libraries directly to process data

#### CHANGELOG.md

Documents version history and changes:

```markdown
# Changelog

## [1.2.3] - 2023-04-09
### Fixed
- Bug in parameter handling

## [1.2.0] - 2023-03-15
### Added
- New optional parameter `param2`

## [1.1.0] - 2023-02-01
### Added
- Support for new input format
- Improved algorithm performance

## [1.0.0] - 2023-01-01
### Added
- Initial release
```

#### test/

Contains tests to verify the wrapper works correctly:

- `env.yaml`: Test environment with additional testing dependencies
- `Snakefile`: Test workflow that runs the wrapper with different parameters
- `tests.py`: Test script that verifies outputs
- `test_data.*`: Sample data files for testing

## Versioning Wrappers

Wrappers should follow semantic versioning (MAJOR.MINOR.PATCH) to ensure reproducibility and compatibility:

### Version Number Format: MAJOR.MINOR.PATCH

- **MAJOR**: Breaking changes that affect input/output formats or metrics
- **MINOR**: New features or parameters that maintain backward compatibility
- **PATCH**: Bug fixes and minor improvements

### Implementation

1. **Version in wrapper.py**:
   ```python
   __version__ = "1.2.3"
   ```

2. **Version in README.md**:
   ```markdown
   Version: 1.2.3
   ```

3. **CHANGELOG.md**:
   Document all changes between versions

4. **Git Tags**:
   Tag each release with the version number

### Versioning Guidelines

- **MAJOR Version (1.0.0 → 2.0.0)**: Increment when input/output formats change
- **MINOR Version (1.0.0 → 1.1.0)**: Increment when new optional parameters are added
- **PATCH Version (1.0.0 → 1.0.1)**: Increment when bug fixes are made

### Version Specification in Workflows

When using wrappers in workflows, specify the version explicitly:

```python
rule example:
    input:
        "input.fasta"
    output:
        "output.fasta"
    params:
        param1="value"
    wrapper:
        "file:path/to/wrapper@1.2.3"  # Specify version
```

## Creating a New Wrapper

### Step 1: Choose the Wrapper Type

- Use a **shell wrapper** for command-line tools with a well-defined CLI
- Use an **API wrapper** for Python libraries or tools without a CLI

### Step 2: Create the Directory Structure

```bash
mkdir -p wrapper_name/test
touch wrapper_name/environment.yaml
touch wrapper_name/README.md
touch wrapper_name/wrapper.py
touch wrapper_name/CHANGELOG.md
touch wrapper_name/test/env.yaml
touch wrapper_name/test/Snakefile
touch wrapper_name/test/tests.py
```

### Step 3: Define Dependencies

Create the `environment.yaml` file with the required dependencies.

### Step 4: Document the Wrapper

Create the `README.md` file with comprehensive documentation.

### Step 5: Implement the Wrapper

Write the `wrapper.py` script to execute the tool.

### Step 6: Create Tests

Implement tests to verify the wrapper works correctly.

## Using Wrappers in Workflows

Wrappers can be used in Snakemake workflows by specifying the path to the wrapper:

```python
rule example:
    input:
        "input.file"
    output:
        "output.file"
    params:
        param1="value",
        param2="value"
    wrapper:
        "file:path/to/wrapper"
```

## Best Practices

1. **Version Dependencies**: Always specify version numbers for dependencies
2. **Document Parameters**: Clearly document all parameters with descriptions and defaults
3. **Provide Examples**: Include example usage in the README
4. **Write Tests**: Create comprehensive tests for all functionality
5. **Handle Errors**: Implement proper error handling and logging
6. **Follow Conventions**: Adhere to the established wrapper structure and conventions
7. **Version Wrappers**: Use semantic versioning for wrappers

## Resources

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [Snakemake Wrapper Repository](https://github.com/snakemake/snakemake-wrappers)
- [Conda Documentation](https://docs.conda.io/)
- [Bioconda Documentation](https://bioconda.github.io/)
- [Semantic Versioning](https://semver.org/)

## Example Wrappers

This directory contains two example wrappers:

1. **Shell Command Wrapper** (`shell/`): Wraps the `seqkit translate` command
2. **Python API Wrapper** (`api/`): Wraps the Biopython `Seq.translate()` method

Both wrappers accomplish the same task (translating DNA sequences to protein sequences) but use different approaches, highlighting the trade-offs between shell and API wrappers. 