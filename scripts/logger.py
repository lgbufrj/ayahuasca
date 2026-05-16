import logging
import sys

RESET  = "\033[0m"
BOLD   = "\033[1m"
GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
CYAN   = "\033[36m"
GREY   = "\033[90m"

LEVEL_COLOURS = {
    logging.DEBUG:    GREY,
    logging.INFO:     CYAN,
    logging.WARNING:  YELLOW,
    logging.ERROR:    RED,
    logging.CRITICAL: RED + BOLD,
}


class ColourFormatter(logging.Formatter):
    """Coloured, human-friendly log formatter."""

    FMT = "{asctime}  {levelname:<8}  {message}"
    DATEFMT = "%H:%M:%S"

    def format(self, record: logging.LogRecord) -> str:  # noqa: A003
        colour = LEVEL_COLOURS.get(record.levelno, RESET)
        record.levelname = f"{colour}{record.levelname}{RESET}"
        record.msg        = f"{colour}{record.msg}{RESET}"
        return super().format(record)


def build_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(
        ColourFormatter(
            fmt=ColourFormatter.FMT,
            datefmt=ColourFormatter.DATEFMT,
            style="{",
        )
    )
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
    return logger