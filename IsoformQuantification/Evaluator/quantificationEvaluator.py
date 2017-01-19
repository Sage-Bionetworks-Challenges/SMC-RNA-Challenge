#! /usr/bin/env python

import sys
import os
import math
import getopt
import warnings
import numpy
from collections import Counter 
from collections import defaultdict
from scipy import stats

def usage():
    print """
    quantificationEvaluator -t <truth-tsv>  -i <input-tsv>
    
    Requested Parameters:
        -t/--truth-tsv    [string:    path to truth tsv]
        -i/--input-tsv         [string:    path to input tsv]

    Version:                    1.0.0
          """


#parameters
inFile = ''               #input tsv
inTruthFile = ''          #input truth tsv

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"ht:i:",["help",
                                                 "truth-tsv=",
                                                 "input-tsv="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-i", "--input-tsv"):
            global inFile
            inFile= arg
        elif opt in ("-t","--truth-tsv"):
            global inTruthFile
            inTruthFile = arg



input_values_dic = defaultdict(lambda: 0.0)

def getInputDic():
    infile = "%s" % inFile
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        elif not line.startswith("ENST"):
            continue
        else:
            tmp=line.split("\t")
            name=tmp[0]
            value=float(tmp[1][0:len(tmp[1])-1])
            input_values_dic[name]=value

truth_values = []
input_values = []

def getBothValues():
    infile = "%s" % inTruthFile
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        elif not line.startswith("ENST"):
            continue
        else:
            tmp=line.split("\t")
            name=tmp[0]
            value=float(tmp[1][0:len(tmp[1])-1])
            truth_values.append(value)
            if input_values_dic.get(name) is not None:
                input_values.append(input_values_dic[name])
            else:
                input_values.append(0)

def calculateCor():
    cor,p_value=stats.spearmanr(truth_values,input_values)
    pearson,pearson_pvalue=stats.pearsonr(truth_values,input_values)
    log_pearson,log_pearson_pvalue=stats.pearsonr(numpy.log(numpy.add(truth_values,0.01)),numpy.log(numpy.add(input_values,0.01)))

    final = "spearman\tpearson\tlog_pearson\n%s\t%s\t%s" % (cor,pearson,log_pearson) 

    print(final)
    return(final)

def main(argv):
    getParameters(argv[1:])
    if inFile=='' or inTruthFile=='':
        usage()
        return 1
    getInputDic()
    getBothValues()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final = calculateCor()
        with open("result.out",'w') as results:
            results.write(final)
            results.close()
    return(0)
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
