"""
WIG scaler . . . script to scale a pair of WIG files.
"""
import optparse 
import sys

op = optparse.OptionParser()
op.add_option('-i',dest='input',help='Input file to scale to')
options,args = op.parse_args()

# sum the values in each file
sys.stderr.write('summing values in %s\n' % options.input)
count = 0
for line in open(options.input):
    if line.startswith('track'):
        continue
    if line.startswith('fixedStep'):
        continue
    count += float(line.strip())

# scale it down!
count /= 1e6

for line in open(options.input):
    if line.startswith('track'):
        sys.stdout.write(line)
        continue
    if line.startswith('fixedStep'):
        sys.stdout.write(line)
        continue

    value = float(line.strip())
    value /= count
    sys.stdout.write('%s\n'%value)

