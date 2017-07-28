#!/usr/bin/env python3
import sys, linecache
infile = sys.argv[1]
levels_num = len(linecache.getline(infile, 3).strip().split("\t")[0].split(";"))

def biom2txt2spf():
    with open(infile) as hd:
        next(hd)
        for line in hd:
            cols = line.strip().split("\t")
            sample_num = len(cols)
            if line.startswith("#"):
                header = cols[1:]
                new_header = ["Class" + str(i+1) for i in list(range(levels_num))]
                print("\t".join(new_header) + "\t" + "\t".join(header))
            else:
                levels = cols[0].split(";")
                new_header = ["Class" + str(i) for i in list(range(levels_num))]
                print("\t".join(levels) + "\t" + "\t".join(cols[1:]))
biom2txt2spf()
