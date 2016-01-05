#!python2
# -*- coding: utf-8 -*-

# ***************************
# Antoine Juckler
#
# pdf2img.py
# Last modified: 2016/01/05
# ***************************

import getopt
import sys
import os
import re

from subprocess import Popen

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:gb', ['input=', 'gif', 'stop', 'bunch', 'delay=', 'density=', 'fixed'])
except getopt.GetoptError as err:
	print(err)
	sys.exit(2)

gif = False
bunch = False
delay = 60
density = 150
infile = ""
fixed = False

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
	elif o in ('--density'):
		density = eval(a)
	elif o in ('--fixed'):
		fixed = True
	else:
		print 'Unknown option: ' + format(o)
		sys.exit(2)

if infile == "":
	print 'No input file specified'
	sys.exit(2)

# Run
if bunch:
	maxval = 20
else:
	maxval = 1

if infile.endswith('_BW'):
	infile = infile[:-len('_BW')]
	BW = '_BW'
else:
	BW = ''


if fixed and re.search("\_[0-9]+$", infile) != None:
	m = re.search("\_[0-9]+$", infile)
	fixed_dist = m.group(0)
	infile = infile[:-len(fixed_dist)]
else:
	fixed_dist = ''


names = []
for i in range(0, maxval):
	if maxval != 1:
		filename = infile + str(i) + fixed_dist + BW
	else:
		filename = infile + fixed_dist + BW
	if os.path.isfile(filename + '.pdf'):
		names.append('img/' + filename + '.png')
		args = ['convert_magick', '-density', str(density), filename + '.pdf', names[-1]]
		Popen(args).communicate()[0]
		print filename + ' converted'

if gif:
	args = ['convert_magick', '-delay', str(delay), '-density', str(density), '-dispose', 'background']
	args.extend(names)
	args.append('img/' + infile + fixed_dist + BW + '.gif')
	Popen(args).communicate()[0]
	print 'gif generated'
