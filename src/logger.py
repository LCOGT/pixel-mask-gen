# needed to make sphinx work, but super hacky. Please fix
try:
    import lcogt_logging
except ModuleNotFoundError:
    pass

import logging
import os
import sys

# Taken from: https://github.com/LCOGT/lcogt_logging/blob/master/example.py

def setup_custom_logger():
    logger = logging.getLogger(__name__)
    helpful_info = {'pid': os.getpid()}
    formatter = lcogt_logging.LCOGTFormatter(extra_tags=helpful_info)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.NOTSET)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger

