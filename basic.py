import sys
from subprocess import *


class Basic(object):
    def __init__(self):
        pass

    @staticmethod
    def run(cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir, stdout=PIPE)
        p.wait()
        return p