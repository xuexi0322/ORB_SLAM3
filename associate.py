#!/usr/bin/python
# Software License Agreement (BSD License)
#
# Copyright (c) 2013, Juergen Sturm, TUM
# All rights reserved.

import sys
import numpy

def read_file_list(filename):
    file = open(filename)
    data = file.read()
    lines = data.replace(","," ").replace("\t"," ").split("\n") 
    list = [[v.strip() for v in line.split(" ") if v.strip()!=""] for line in lines if len(line)>0 and line[0]!="#"] 
    list = [(float(l[0]),l[1:]) for l in list if len(l)>1]
    return dict(list)

def associate(first_list, second_list, offset, max_difference):
    first_keys = list(first_list.keys())
    second_keys = list(second_list.keys())
    potential_matches = [(abs(a - (b + offset)), a, b) 
                         for a in first_keys 
                         for b in second_keys 
                         if abs(a - (b + offset)) < max_difference]
    potential_matches.sort()
    matches = []
    for diff, a, b in potential_matches:
        if a in first_keys and b in second_keys:
            first_keys.remove(a)
            second_keys.remove(b)
            matches.append((a, b))
    
    matches.sort()
    return matches

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: associate.py first_file second_file [offset] [max_difference]")
        sys.exit(1)
        
    first_list = read_file_list(sys.argv[1])
    second_list = read_file_list(sys.argv[2])

    offset = 0.0
    if len(sys.argv) > 3:
        offset = float(sys.argv[3])
    max_difference = 0.02
    if len(sys.argv) > 4:
        max_difference = float(sys.argv[4])
    
    matches = associate(first_list, second_list, offset, max_difference)    
    
    for a,b in matches:
        print("%f %s %f %s"%(a," ".join(first_list[a]),b-offset," ".join(second_list[b])))
