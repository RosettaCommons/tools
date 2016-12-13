#!/usr/bin/env python

## @file   ss_pred_converter.py
#  @brief  Converts various secondary structure profile formats to a PsiPred-SS2 file
#  @author Dominik Gront (dgront@chem.uw.edu.pl)

import sys,re
import getopt

class ResiduePrediction:
  id = -1 
  aa = 'X'
  ss = 'C'
  h = 0.0
  e = 0.0
  c = 1.0
 
  def __init__(self,id,aa,ss,h,e,c):
    self.id = id
    self.aa = aa
    self.ss = ss
    self.h = h
    self.e = e
    self.c = c

def usage():
    usage = """Converts Juffo and SAM predictions to PsiPred-SS2 file format
    -h --help                     Prints this message
    -s --sam filename             Converts SAM    -> PsiPred SS2
    -j --juffo filename           Converts JUFFO  -> PsiPred SS2
    -t --talos filename           Converts TALOS+ -> PsiPred SS2
    -p --porter filename          Converts PORTER -> PsiPred SS2
    """
    print usage

def isFloat(s):
   try: return float(s) or True
   except (ValueError, TypeError), e: return False

def readTalos(fileName) :
  out = []
  for line in open(fileName) :
    entries = re.split(" +",line.strip())
    if len(entries) != 9 : continue
    if isFloat(entries[5]) and isFloat(entries[6]) and isFloat(entries[4]) :
      e = float(entries[5])
      h = float(entries[4])
      c = float(entries[6])
      ss = entries[8]
      if ss == 'S' : ss = 'E'
      if ss == 'U' : ss = 'L'
      out.append( ResiduePrediction(int(entries[0]),entries[1],ss,h,e,c) )

  return out

def readPorter(fileName) :
  out = []
  f = open(fileName)
  skipIt = f.readline().strip()
  seq = f.readline().strip()
  sec = f.readline().strip()
  h = re.split("\s+",f.readline().strip())
  e = re.split("\s+",f.readline().strip())  
  c = re.split("\s+",f.readline().strip())  

  for i in range(len(seq)) :
      out.append( ResiduePrediction(i+1,seq[i],sec[i],float(h[i]),float(e[i]),float(c[i])) )

  return out

def readJuffo(fileName) :
  out = []
  for line in open(fileName) :
    entries = re.split(" +",line.strip())
    if len(entries) != 6 : continue
    if isFloat(entries[5]) and isFloat(entries[3]) and isFloat(entries[4]) :
      e = float(entries[5])
      h = float(entries[4])
      c = float(entries[3])
      ss = entries[2]
      if ss == 'S' : ss = 'E'
      if ss == 'U' : ss = 'L'
      out.append( ResiduePrediction(int(entries[0]),entries[1],ss,h,e,c) )

  return out

def readSam(fileName) :
  out = []
  plusOne = 0
  for line in open(fileName) :
    if line[0] == '#' : continue
    entries = re.split("\s+",line.strip())
    if len(entries) != 5 : 
      continue
    if isFloat(entries[2]) and isFloat(entries[3]) and isFloat(entries[4]) :
      e = float(entries[2])
      h = float(entries[3])
      c = float(entries[4])
      if int(entries[0]) == 0 : plusOne = 1
      if e >= h and e >= c :
        out.append( ResiduePrediction(int(entries[0])+plusOne,entries[1],'E',h,e,c) )
      if h > e and h >= c :
        out.append( ResiduePrediction(int(entries[0])+plusOne,entries[1],'H',h,e,c) )
      if c > e and c > h :
        out.append( ResiduePrediction(int(entries[0])+plusOne,entries[1],'L',h,e,c) )

  return out

def printSS2(ss_entries) :
  print "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n"
  for e in ss_entries:
    print "%4d %c %c   %5.3f  %5.3f  %5.3f" % (e.id,e.aa,e.ss,e.c,e.h,e.e)

def main():
    if len(sys.argv) == 1 : usage()
    # parse command line options
    try:
      opts, args = getopt.getopt(sys.argv[1:], "hs:j:t:p:", ["help","sam=","juffo=","talos=","porter="])
    except getopt.error, msg:
      print msg
      print "for help use --help"
      sys.exit(2)
    # process options
    for o, a in opts:
      if o in ("-h", "--help") :
        usage()
        sys.exit(0)
      if o in ("-s", "--sam") :
        printSS2( readSam(a) )
        continue
      if o in ("-j", "--juffo") :
        printSS2( readJuffo(a) )
        continue
      if o in ("-t", "--talos") :
        printSS2( readTalos(a) )
        continue
      if o in ("-p", "--porter") :
        printSS2( readPorter(a) )
        continue
      assert False, "unhandled option: "+o

if __name__ == "__main__":
    main()
