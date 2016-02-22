#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os
import csv

simplenodelist = []
otudict = {}
edgelookupdict = {}
otulookupdict = {}
root = ''
publication = ''
def main():
    eulerfile = sys.argv[1]
    convfile = sys.argv[2]
    file_name, extension = os.path.splitext(convfile)
    outputfile = '%s_nl.txt' % file_name
    eulertxt = ''
    convdict = {}
    if os.path.exists(eulerfile):
        f = open(eulerfile, 'r')
        eulertxt = f.read()
    else:
        print "can't find the euler file"
        return;
    if os.path.exists(convfile):
        with open(convfile, 'rb') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                if len(row) < 2 or row[1] == '' or row[1].isspace():
                    print "conversion file does not have new names"
                    return
                convdict[row[0]] = row[1].replace(' ','_')
    else:
        print "can't find the conversion file"
        return;

    convkeys = convdict.keys()
    # sort by length of strings, so we're replacing longest strings to shortest strings
    convkeys.sort(cmp=lambda x,y: cmp(len(y), len(x)))
    for key in convkeys:
        eulertxt = eulertxt.replace(key, convdict[key])
    outf = open(outputfile, 'w')
    outf.write(eulertxt)
    outf.close()
    
    
if __name__ == '__main__':
    main()
