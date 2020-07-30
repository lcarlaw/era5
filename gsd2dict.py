"""
Simple script to convert GSD airport metadata from 
https://rucsoundings.noaa.gov/metar.short to a python-readable dictionary.

This script then update airports.py with the proper dictionary. 
"""

import pprint
metadata = {}

with open('airports.txt') as temp: data = temp.readlines()
knt = 0
for line in data:

    try:
        line = line.strip().split(' ')
        line = filter(None, line)
        print line
    
        # Grab the header information. 
        if knt == 0:
            header = line

        else:
            num_entries = len(line)
            n = 6
            desc = ''
            while n < num_entries:
                desc += line[n]
                n += 1
            metadata[line[0]] = line[2], line[3], line[4], desc
        knt += 1
    except:
        pass
f = open('./airports.py', 'w')
f.write('metadata = ')
pprint.pprint(metadata, stream=f)
