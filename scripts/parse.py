#!/usr/bin/python3
import os
import re
import argparse

parser = argparse.ArgumentParser();
parser.add_argument('file');
args = parser.parse_args()

f = open(args.file)
for line in f.readlines():
    pattern = re.compile(" *[0-9]+:(.*)\((.*)\).*\[(.*)\]")
    if pattern.match(line):
        addr=pattern.findall(line)[0]
        os.system("addr2line -Cfie " + str(addr[0]) + " " + str(addr[1]))
        print()


