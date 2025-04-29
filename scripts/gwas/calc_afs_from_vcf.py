#!/usr/bin/env python
import sys
import os
import subprocess
#import re
######################################################
#print sys.argv[1]
#f = open(sys.stdin)
for line in sys.stdin:
#while True:
#    line = f.readline()
    if not line: break
    geno_split = line.split('\t')
    REF = 0
    ALT = 0
    for i in geno_split:
        if i == "0/0" or i == "0|0":
            REF += 2
        elif i == "1/1" or i == "1|1":
            ALT += 2
        elif i == "0/1" or i == "0|1" or i == "1|0":
            ALT += 1
    print(REF,"\t",ALT)
#f.close
