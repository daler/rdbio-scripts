#!/usr/bin/python
import os
import sys
try:
    seq = sys.argv[1]
except IndexError:
    print 'Usage: %s 1.seq (creates a new PNG, 1.png)'
    sys.exit(1)
output = os.path.splitext(seq)[0]+'.png'
cmds = ['weblogo',
        '-F','png_print',
        '-f',seq,
        '--size','large',
        '-o',output,
        '--sequence-type','dna',
        '--color red "A" A',
        '--color green "T" T',
        '--color blue "C" C',
        '--color orange "G" G']
os.system(' '.join(cmds))
