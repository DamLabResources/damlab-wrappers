import click # type: ignore
import logging
import sys

log = logging.getLogger('multiqc.modules.damlab')

# Define the CLI option
debug_option = click.option(
    '--my-testing-debug',
    is_flag=True,
    help='Print debug information for DAMlab modules and exit'
)

def check_debug(config):
    """Check if debug flag is set and print debug info if needed."""
    if config.kwargs.get('my_testing_debug'):
        print("\nDAMlab MultiQC Plugin Debug Info:")
        print("================================")
        print(f"Python version: {sys.version}")
        print(f"Python path: {sys.path}")
        print("\nMultiQC Configuration:")
        print("---------------------")
        print(f"Search paths: {config.search_patterns}")
        print(f"Module order: {config.module_order}")
        print(f"Available modules: {config.available_modules}")
        print("\nExiting due to --my-testing-debug flag...")
        sys.exit(0) 