"""MultiQC plugin for DAMlab wrappers."""

try:
    from .modules.dorado import DoradoModule
    from .modules.strainline import StrainlineModule
except ImportError as e:
    # Handle import errors more gracefully
    import logging
    log = logging.getLogger('multiqc.modules.damlab')
    log.debug(f"Could not load module: {str(e)}") 