import logging
import sys

logging.basicConfig(
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    level=logging.INFO,
    handlers=[
        logging.FileHandler("TestLog.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('TestLogger')


# logging with showtraceback if working in Notebook
def exception_handler_IP(self, etype, evalue, tb, tb_offset=None):
    logger.error("An Error occured while executing CIDER!", exc_info=True)
    self.showtraceback((etype, evalue, tb), tb_offset=tb_offset)


# logging with excepthook if working in command line
def handle_exception(exc_type, exc_value, exc_traceback):
    logger.exception(
        "An Error occured while executing CIDER",
        exc_info=(exc_type, exc_value, exc_traceback)
    )
    return


try:
    __IPYTHON__
    from IPython import get_ipython
    ip = get_ipython()
    ip.set_custom_exc((BaseException,), exception_handler_IP)
except NameError:
    sys.excepthook = handle_exception

