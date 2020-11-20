import sys
import os
from subprocess import *


class Basic(object):
    def __init__(self):
        pass

    @staticmethod
    def run(cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir, stdout=PIPE)
        p.wait()
        return(p)
        
    @staticmethod    
    def read_arguments(arg_file=None):
        """
        this method used to parse the config file 
        arg_file: config file
        return: dict
        """
        arguments = {}
        handle = open(arg_file, "r")
        for line in handle:
            #print(line)
            line = line.rstrip("\r\n")
            if (not line):
                continue
            if (line[0] == "#"):
                continue
            line = line.split("=")
            arg = line[0].strip()
            value = line[1].strip()
            #print(os.path.isfile(value))
            if not os.path.isfile(value):
                print("excute file %s dons't exsits!" % value)
                exit(1)
            arguments[arg] = value
        handle.close()
        return(arguments)
        