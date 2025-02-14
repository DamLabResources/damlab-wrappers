from multiqc import config # type: ignore

import logging
log = logging.getLogger('multiqc.modules.strainline')
log.critical('Importing hook')


def add_config():
    """Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Add to the search patterns used by modules
    if "strainline" not in config.sp:
        config.update_dict(config.sp, {"strainline": {"fn": "*.strainline.*",
                                                      "contents": "# Strainline MultiQC Log",
                                                      "num_lines": 10000}})
        log.critical('updated in damlabqc hook {}'.format(config.sp['strainline']))
        log.critical('Other SP example {}'.format(config.sp['afterqc']))
