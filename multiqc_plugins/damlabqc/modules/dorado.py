from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule # type: ignore
from multiqc.plots import linegraph, bargraph # type: ignore
import logging
import os
import yaml # type: ignore
from collections import OrderedDict
import sys
from ..cli import check_debug

# Rest of file remains the same... 