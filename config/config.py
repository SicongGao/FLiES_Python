import logging
import logging.config

INPUT_MODE = 1  # 1: input from files; 0: input manually
OUTPUT_PATH = "../result/"
LOG_PATH = "../config/logger.config"
CPU_MANUAL = 1  # 1: read value from CPUS; 0: rely on the machine cpu amount
CPUS = 32


def log_init():
    logging.config.fileConfig(LOG_PATH)
    logging.info("Finished log initialize.")

