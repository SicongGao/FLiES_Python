import logging
import logging.config

INPUT_MODE = 1  # 1: input from files; 0: input manually
OUTPUT_PATH = "../result/"
LOG_PATH = "../config/logger.config"


def log_init():
    logging.config.fileConfig(LOG_PATH)
    logging.info("Finished log initialize.")

