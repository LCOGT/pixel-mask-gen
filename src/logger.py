import lcogt_logging
import logging
import os
import sys

def setup_custom_logger(name='pixel-mask-gen'):
    # Taken from: https://github.com/LCOGT/lcogt_logging/blob/master/example.py

    logger = logging.getLogger(__name__)
    helpful_info = {'pid': os.getpid()}
    formatter = lcogt_logging.LCOGTFormatter(extra_tags=helpful_info)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.NOTSET)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    logger.setLevel(logging.DEBUG)

    return logger

global logger_obj
logger_obj = setup_custom_logger()
#logger_obj.setLevel(logging.DEBUG)
