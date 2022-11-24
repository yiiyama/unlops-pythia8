#!/usr/bin/env python

"""Extract a specified number of events from an LHE file."""

import sys

input_name = sys.argv[1]
event_start = int(sys.argv[2])
num_events = int(sys.argv[3])
output_name = sys.argv[4]

input_lhe = open(input_name)
output_lhe = open(output_name, 'w')

for line in input_lhe:
    output_lhe.write(line)

    if '</init>' in line:
        break

ievt = 0
for line in input_lhe:
    if ievt >= event_start:
        output_lhe.write(line)

    if '</event>' in line:
        ievt += 1

    if ievt == event_start + num_events:
        break

output_lhe.write('</LesHouchesEvents>\n')

output_lhe.close()
input_lhe.close()
