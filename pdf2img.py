#!python2
# -*- coding: utf-8 -*-

# ***************************
# Antoine Juckler
#
# pdf2img.py
# Last modified: 2015/12/24
# ***************************

import getopt
import sys
import os

from subprocess import Popen

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:gb', ['input=', 'gif', 'stop', 'bunch', 'delay'])
except getopt.GetoptError as err:
	print(err)
	sys.exit(2)

gif = False
bunch = False
delay = 60

for o, a in opts:
	if o in ('-i', '--input'):
		infile = str(a)
		if infile.endswith('.pdf'):
			infile = infile[:-len('.pdf')]
	elif o in ('-g', '--gif'):
		gif = True
	elif o in ('-b', '--bunch'):
		bunch = True
	elif o in ('--delay'):
		delay = eval(a)
	else:
		print 'Unknown option: ' + format(o)
		sys.exit(2)

# Run
if bunch:
	maxval = 20
else:
	maxval = 1

names = []
for i in range(0, maxval):
	if maxval != 1:
		filename = infile + str(i + 1)
	else:
		filename = infile
	if os.path.isfile(filename + '.pdf'):
		names.append(filename + '.pdf')
		args = ['convert_magick', filename + '.pdf', filename + '.png']
		Popen(args).communicate()[0]
		print filename + ' converted'

if gif:
	namestr = ''
	for i in range(0, len(names) - 1):
		namestr += names[i] + ' '
	namestr += names[-1]

	args = ['convert_magick', '-delay', str(delay), '-density', '150', '-dispose', 'background']
	args.extend(names)
	args.append(infile + '.gif')
	Popen(args).communicate()[0]
	print 'gif generated'
