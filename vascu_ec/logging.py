import logging
import sys

LOGGER_NAME = 'vascu_ec'


def get_logger():
    return logging.getLogger(LOGGER_NAME)  # root logger


def configure_logger(loglevel=None, stream_handler=None, formatter_string=None):
    logger = logging.getLogger(LOGGER_NAME)

    # create formatter
    if not formatter_string:
        formatter = get_default_formatter()
    else:
        formatter = logging.Formatter(formatter_string)

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(loglevel)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # additional stream handler?
    if stream_handler:
        ch = logging.StreamHandler(stream_handler)
        ch.setLevel(loglevel.name)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger


def get_default_formatter():
    return logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')