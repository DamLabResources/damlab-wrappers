from multiqc.utils import config # type: ignore


def add_config():
    """Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Add to the search patterns used by modules
    if "strainline" not in config.sp:
        config.update_dict(config.sp, {"strainline": {"fn": "*.strainline.yaml"}})
    
    # Add to the module order
    config.module_order = ["strainline"] + config.module_order
