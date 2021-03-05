#!/usr/bin/env python3
# This script reduces the states within a BED file returned from ChromHMM and assigns new colors for visualization with UCSC

import re
import sys

input = open(sys.argv[1], 'r')
output = open(sys.argv[2], 'w')
strongRE_input = sys.argv[3]
weakRE_input =  sys.argv[4]
poisedRE_input = sys.argv[5]
noRE_input = sys.argv[6]

strongRE_list = strongRE_input.split(",") # state1
weakRE_list = weakRE_input.split(",") # state2
poisedRE_list = poisedRE_input.split(",") # state3
noRE_list = noRE_input.split(",") # state4


for line in input:
	if re.match('track', line):
		output.write(line)
	else:
		tab = re.split('\t', line)
		state = tab[3]
		color = tab[8]
		if state in strongRE_list:
			color = '030,080,255'
			state = 1
		elif state in weakRE_input:
			color = '050,200,255'
			state = 2
		elif state in poisedRE_input:
			color = '189,188,188'
			state = 3
		elif state in noRE_input:
			color = '255,255,255'
			state = 4
		else:
			break
		output.write(tab[0] + '\t' + tab[1] + '\t' + tab[2] + '\t' + str(state) + '\t' + tab[4] + '\t' + tab[5] + '\t' + tab[6] + '\t' + tab[7] + '\t' + color + '\n')

input.close()
output.close()
