from __future__ import annotations

import logging


def configure_logging(level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("ppt_patient_parser")
    if logger.handlers:
        logger.setLevel(level.upper())
        return logger

    logger.setLevel(level.upper())
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger
