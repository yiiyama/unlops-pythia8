#!/usr/bin/env python

"""Merge multiple LHE files into several larger files."""

import os
import sys
import subprocess
import shutil

# The path of a file that lists the source LHE files.
source_list = sys.argv[1]
# Start line number within the list.
start = int(sys.argv[2])
# End line number within the list.
end = int(sys.argv[3])
# Number of source files to merge into each output.
unit = int(sys.argv[4])

output_path = '/eos/cms/store/user/yiiyama/unlops_lhe_merged/%s/%s'

channel = os.path.basename(source_list).replace('.list', '')

lhes = []

with open(source_list) as source:
    iline = 0
    for line in source:
        if iline >= start:
            lhes.append(line.strip())
        
        iline += 1
            
        if iline == end:
            break

ifile = 0
while start + ifile != end:
    outname = os.environ['TMPDIR'] + ('/%s_%d.lhe' % (channel, start + ifile))
    
    output = open(outname, 'w')

    try:
        path = lhes[ifile]
    except IndexError:
        break

    print path
    with open(path) as lhe:
        for line in lhe:
            if '</LesHouchesEvents>' in line:
                break
    
            output.write(line)

    ifile += 1

    while ifile % unit != 0:
        try:
            path = lhes[ifile]
        except IndexError:
            break

        print path
        with open(path) as lhe:
            events_block = False
            for line in lhe:
                if '</LesHouchesEvents>' in line:
                    break
    
                if events_block:
                    output.write(line)
                    continue
                
                if '</init>' in line:
                    events_block = True

        ifile += 1
    
    output.write('</LesHouchesEvents>\n')
    output.close()

    shutil.copyfile(outname, output_path % (channel, os.path.basename(outname)))
    os.unlink(outname)
