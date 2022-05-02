import logging, sys

# def setupLogging(debug=False, output=None):
#     """
#         Setup for logging
#         
#         Inputs:
#         - debug (boolean for debug modes yes or no)
#         - output (full path of location to store output log file)
#           (default = currentworkingdir/Typing.log
#     """
#     if output==None:
#         logLevel = logging.DEBUG if debug else logging.INFO
#         logFormat = "%(asctime)s [%(levelname)s] %(message)s"
#         logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
#         logging.info("Running %s" % " ".join(sys.argv))
#     else:
#         logLevel = logging.DEBUG if debug else logging.INFO
#         logFormat = "%(asctime)s [%(levelname)s] %(message)s"
#         logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat,
#                             filename=output,
#                             filemode='w')
#         console = logging.StreamHandler()
#         console.setLevel(logLevel)
#         formatter = logging.Formatter(logFormat)
#         console.setFormatter(formatter)
#         # add the handler to the root logger
#         logging.getLogger('').addHandler(console)
#         logging.info("Running %s" % " ".join(sys.argv))

# def setupLogging(debug=False, output=None):
#     logging.basicConfig(level=logging.DEBUG,
#                         filename="example2.log",
#                         format='%(asctime)s %(message)s',
#                         handlers=[logging.StreamHandler()])
    
def setupLogging(debug=False):    
    # create logger
    logger = logging.getLogger('test')
    
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    
    # create a file handler
    handler = logging.FileHandler('test.log')
    
    if debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
        handler.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # add formatter to ch
    ch.setFormatter(formatter)
    
    # create a logging format
    handler.setFormatter(formatter)
    
    # add ch to logger
    logger.addHandler(ch)
    
    # add the handlers to the logger
    logger.addHandler(handler)
    
    return logger
