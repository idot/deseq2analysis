#!/usr/bin/env python 

import sys
import os.path
import argparse

# --file a --file b ...  all input files are merged  
# --column=X as indicator of which column to use for merge, default = 1, can be anywhere
# --round if output should be rounded to integer
# --outfile output file path 

def parseArgs():
    parser = argparse.ArgumentParser(description='Merge one column from multiple files. Column 0 from 1. file will be in column 0.')
    parser.add_argument('--file', action='append', help='input file to merge')
    parser.add_argument('--column', type=int , default=1, help='column to use for merge, 0 based')
    parser.add_argument('--round', action='store_true', help='if values should be rounded')
    parser.add_argument('--header', action='store_true', help='if a header row should be created')
    parser.add_argument('--outfile', help='output file path')
    args = parser.parse_args()
    return args


# reads a column (0 indexed) from file 
# 
#
def readColumn(ifile, column, doRound, doHeader):
    result = []
    with open(ifile) as f:
      lineNr = 0
      for line in f.readlines():
         try:
            items = line.strip().split("\t")
            #if lineNr == 0 and (len(items) < column + 1 or (len(items) > 2 and items[column] == '')): #first row add filename if missing in header was 'autodetect from htseq format' now explicit doHeader
            if doHeader:
               fn = os.path.basename(ifile)
               result.append(fn)
            else:
               c = items[column].strip()
               if doRound:
                  c = str(int(round(float(c))))
               result.append(c)
            lineNr += 1
         except:
            print("error parsing file lineNr: "+str(lineNr)+" column: "+str(column)+" "+line) #+" "+sys.exc_info()[0]
            raise
    return result

def checkRowcountsEqual(arrs):
   rowcounts = map(lambda a: len(a), arrs)
   rowset = set(rowcounts)
   if not len(rowset) == 1:
      sys.stderr.write("numbers of rows differ in files: "+(",".join([str(c) for c in rowcounts]))+"\n")
      sys.exit(1)

def mergeColumns(files, column: int, doRound: bool, doHeader: bool):
   names = readColumn(files[0], 0, False)
   counts = map(lambda l: readColumn(l, column, doRound, doHeader), files)
   counts.insert(0, names)
   checkRowcountsEqual(counts)
   resarr = zip(*counts)
   return resarr

def writeOut(resarr, outfile):
    with open(outfile, 'w') as f:
       for larr in resarr:
          f.write("\t".join(larr)+"\n")


def doIt():
    args = parseArgs()
    merged = mergeColumns(args.file, args.column, args.round, args.header)
    writeOut(merged, args.outfile) 

doIt()


