#!/usr/bin/env python
import sys
import getopt
import logging

from cclib.parser import guesstype

def moreusage():
    """More detailed usage information"""
    print """Usage:  python ccget.py <attribute> [<attribute>] <compchemlogfile>
where <attribute> is one of the available attributes parsed by cclib
from <compchemlogfile>.
For a list of attributes parsed from the file, type:
     ccget --list     [or -l]"""

def usage():
    """Display usage information"""
    print """Usage:  python ccget.py <attribute> [<attribute>] <compchemlogfile>
Try     ccget --help    for more information"""


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hl", ["help","list"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    showattr = False
    for o, a in opts:
        if o in ("-h", "--help"):
            moreusage()
            sys.exit()
        if o in ("-l", "--list"):
            showattr = True
    if (not showattr and len(args)<2) or (showattr and len(args)!=1): # Need at least one attribute and the filename
        usage()
        sys.exit()

    log = guesstype(args[-1])
    log.logger.setLevel(logging.ERROR)
    log.parse()
    if showattr:
        print "cclib can parse the following attributes from %s:" % args[-1]
        for x in ['aonames','aooveralps','atomcoords','atomnos',
                  'etenergies','etoscs','etrotats','etsecs','etsyms',
                  'fonames','fooverlaps','geotargets','geovalues',
                  'homos','mocoeffs','moenergies','mosyms','natom','nbasis',
                  'nmo','scfenergies','scftargets','scfvalues',
                  'vibfreqs','vibirs','vibramans','vibsyms']:
            if hasattr(log,x):
                print "  %s" % x
    else:
        for arg in args[:-1]:
            print getattr(log,arg)
    return log 

if __name__ == "__main__":
    t = main()
