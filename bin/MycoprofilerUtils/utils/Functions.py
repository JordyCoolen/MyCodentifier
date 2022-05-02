#!/usr/bin/env python

"""
    Set of functions regular used
"""
import subprocess, os, logging, sys, time
from string import Template
from utils.CommandRunner import exe

def check_path(program='nano'):
    """
        Checks if program is in PATH
    """
    try: 
        subprocess.call([program],
                        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except OSError:
        logging.warning("\n\tModule %s is not in PATH" % (program))
        sys.exit() 

def multi_run_wrapper(args):
    """
        Needed to unpack arguments of pairs 
    """
    return add(*args)

def add(x, y):
    """
        Needed to unpack arguments of pairs 
    """
    return x + y

def create_folder(path, name):
    """
        Creates folder with name including
        warning if already exists
    """
    folder = os.path.join(path, name)
    
    if not os.path.exists(folder):
        # create folder
        try:
            os.makedirs(folder)
        except OSError:
            logging.error('Cannot make %s' %(folder))
            sys.exit()
            pass

    return folder

def is_file_done_writing(inputfile, delay = 10):
    """
        Function will see if inputfile size is increasing, if
        inputfile size did not increase during delay time the
        the function will return True, if inputfile does not 
        exist function returns None.
        
        Inputs:
            - inputfile (full path to inputfile)
            - delay (time in seconds between new size measurement)
            
        Output:
            - retcode   True: inputfile did not increase in size (so finished transfering)
                        None: inputfile does not exist
    """
    if os.path.exists(inputfile):
        logging.debug('File exists: %s' %(inputfile))
        statsize = os.stat(inputfile).st_size
        oldsize = os.stat(inputfile).st_size - 1
        while statsize > oldsize:
            oldsize = statsize
            time.sleep(delay)
            statsize = os.stat(inputfile).st_size
        retcode = True
    else:
        retcode = None 
    return retcode

def list_filetype(directory, fileextension):
    """
        Creates list of all files with fileextension in folder
        
        Input
            - directory
            - fileextension (extension without . Example: fasta)
        Output
            - list of files with fileextension
    """
    try:
        listdir = os.listdir(directory)
    except OSError:
        logging.error('Path does not exist')
        return []
        
    newlist = []
    logging.debug('Files in directory:\n%s' % (listdir))
    for f in listdir:
        if f.endswith('.%s' % (fileextension)):
            newlist.append(f)
    
    # no files with fileextension found handling
    if newlist == []:
        logging.warning('No %s files found' % (fileextension))
        return []
            
    return newlist

def get_pairs(inputlist):
    """
        Get pairs from inputlist
        
        Input:
            - [1,2,3,4,5,6]
        Output:
            - [[(1,2),(3,4),(5,6)]]
        
    """
    inputlist = sorted(inputlist)
    newlist = []
    i = iter(inputlist)
    newlist.append(zip(i, i))
    return newlist

def remove_files(inputlist):
    """
        Removes all files from disk that are in the list if possible
        
        Input:
            - inputlist of files (files given as full path)
            
        Output:
            - all files are removed from disk
    """
    
    # removes files in list
    for f in inputlist:
        try:
            os.remove(f)
        except OSError:
            logging.warning('%s not found' % (f))
    
def check_retcode(code):
    """
        Takes return code as input
        and exits if its not 0
    """
    
    if code != 0:
        logging.error('Command failed')
        sys.exit()

class Config():
    """
        Reads the config.xml file and obtains
        paths.
        
        Input read:
            - fullPath of Config.xml
            
        Outputs:
            - dbname (pathoscope reference db name)
            - dbpath (path to pathoscope reference db)
            - ncbipath (path to ncbi reference files)
            - deps (list of all dependencies)
    """
        
    def read(self, fullPath):
        from xml.dom import minidom
        import ast
        logging.debug(fullPath)
        
        # read xml
        xmldoc = minidom.parse(fullPath)
        itemlist = xmldoc.getElementsByTagName('item')
        
        # get xml attributes and fill object
        self.dbname = itemlist[0].attributes['dbname'].value
        self.dbpath = itemlist[1].attributes['dbpath'].value
        self.ncbipath = itemlist[2].attributes['ncbipath'].value
        self.abricatedb = itemlist[3].attributes['abricatedb'].value
        self.deps = ast.literal_eval(itemlist[4].attributes['deps'].value)
        self.wDir = itemlist[5].attributes['wDir'].value

def reverseComp(string):
    """
        Takes a DNA string (lower and/or uppercase) and reverse complement it
        
        Input:
            - DNA string
            
        Output:
            - reverse complemented DNA string
    """
    
    # define the complementation
    complement = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    
    # perform the reverse complementation on string
    result = string.translate(complement)[::-1]
    
    return result

def create_JSON(jobId, location, finishedSuccesfully, errorMessage='', samples=[]):
    """
        JSON object that needs to be returned to the
        BUZZZ.
        
        Inputs:
        - jobId (id of the job needed for BUZZ to follow flow)
        - errorMessage (error message if present default='')
        - location (full path to output folder)
        - samples (can be ignored default=[])
        - finishedSuccessfully (True or False if step was successfully finished)
        
        Outputs:
        - jsonString (string containing parameters in JSON format)
    """
    from json import JSONEncoder
    
    jsonString = JSONEncoder().encode({
        "jobId":jobId,
        "errorMessage":errorMessage,
        "location":location,
        "samples":samples,
        "finishedSuccesfully":finishedSuccesfully})
    
    return jsonString

def sent_to_BUZZ(jsonString, destination, port):
    """
        STOMP protocol to sent messages to BUZZZ
        
        Inputs:
        - jsonString (generated using create_JSON)
        - destination (which box it should be sent to)
        
        Outputs:
        - return 0 if finished
    """
    
    import stomp
    
    c = stomp.Connection([(destination, port)])
    c.start()
    c.connect(wait=True)
    c.send("JobFinishedQueueInput", jsonString)
    c.disconnect()
    
    return 0

def absoluteFilePaths(directory, extension):
    '''
        Yields full paths of files in directory
        with containing extension
    '''
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            if f.endswith(extension):
                yield os.path.abspath(os.path.join(dirpath, f))

def make_dict_with_list_as_value_same_length(dictio, fill):
    '''
        Function will make all list values in a dictionary
        of the same length by filling the list with fill.
        
        Inputs
        - dictio (dictionary with format 
        {key: [data,data,data], key2: [data,data]}
        - fill (element to fill lists, should be a string)
        
        Output
        - dictio (dictionary filled)
        {key: [data,data,data], key2: [data,data,fill]}
        
    '''
    
    maxi = max([len(v) for _,v in dictio.items()])
    for _,v in dictio.items():
        if len(v) < maxi:
            diff = maxi - len(v)
            add = diff*fill
            v.extend(add)
            
    return dictio

def extract_gz(inputfile):
    """
        Extracts .gz file
        
        Inputs:
            - inputfile (full path of .gz file)
            
        Outputs:
            - unpacked file path
    """
    
    # create output file name
    output = inputfile.replace(".gz", "")
    
    # template
    extractfastq = Template("gunzip -c ${input} > ${output}")
        
    # Creating command
    cmd = extractfastq.substitute({"input":inputfile,
                                    "output":output})
    
    c, o, e = exe(cmd)
    logging.debug(c)
    logging.debug(o)
    logging.debug(e)
    check_retcode(c)
    
    return output

def create_gz(inputfile):
    """
        Create .gz file
        
        Inputs:
            - inputfile (full path of file to create archive)
            
        Outputs:
            - file archive with .gz extension
    """
    
    # create output file name
    output = inputfile + '.gz'
    
    # template
    creategz = Template("gzip -f ${input} > ${output}")
        
    # Creating command
    cmd = creategz.substitute({"input":inputfile,
                                    "output":output})
    
    c, o, e = exe(cmd)
    logging.debug(c)
    logging.debug(o)
    logging.debug(e)
    check_retcode(c)
    
    return output
    
    
    
    
    
    
    
    
    
    
    
    
    
    
