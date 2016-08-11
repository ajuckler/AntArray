#!python2
# -*- coding: utf-8 -*-

# ***************************
# Antoine Juckler
#
# pdf2img.py
# Last modified: 2016/01/27
# ***************************

import getopt
import sys
import os
import re
import math

from subprocess import Popen

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:gbt:h',
                                ['input=', 'gif', 'stop', 'bunch', 'delay=',
                                    'density=', 'fixed', 'theta=', 'repeat',
                                    'start=', 'max=', 'help'])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

gif = False
bunch = False
delay = 60
density = 150
infile = ""
fixed = False
theta = 0
ang = False
repeat = False
fixed_dist = ''
start = 0
maxval = 20

for o, a in opts:
    if o in ('-h', '--help'):
        help = open('py_help.txt')
        for l in help.readlines():
            print l
        help.close()
        sys.exit()
    elif o in ('-i', '--input'):
        infile = str(a)
        print infile
        print a
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
    elif o in('-t', '--theta'):
        if a:
            theta = eval(a.replace('pi', 'math.pi'))
        ang = True
    elif o in ('--repeat'):
        repeat = True
    elif o in ('--start'):
        start = eval(a)
    elif o in ('--max'):
        maxval = eval(a)
    else:
        print 'Unknown option: ' + format(o)
        sys.exit(2)

if infile == "":
    print 'No input file specified'
    sys.exit(2)

# if bunch and theta == 0 and 'theta' in infile:
#   print 'Bunch parameter specified, but no theta value given'
#   sys.exit(2)

# Run
if not os.path.exists("img"):  # Check if dir exists and create it
    os.makedirs("img")

if not bunch:
    maxval = 1

if infile.endswith('_BW'):  # Check if BW pattern
    infile = infile[:-len('_BW')]
    BW = '_BW'
else:
    BW = ''

names = []
currang = 0
startang = ''
infile1 = ''
infile2 = ''

if theta != 0 and bunch:  # Create list of files, for theta patterns
    startang = re.search("theta\_([0-9\.]+)\_", infile).group(1)
    infile1 = infile[:infile.find("theta_") + len("theta_")]
    infile2 = infile[infile.find("theta_") + len("theta_" + startang):]

    startang = eval(startang)

    for i in range(0, maxval):
        currang = startang + theta * i
        if currang > math.pi / 2:
            break
        elif currang > 1:
            currang_s = str("{0:.2f}".format(currang))
        else:
            currang_s = str("{0:.3f}".format(currang))
        if "." in currang_s:
            while currang_s[-1] == '0':
                currang_s = currang_s[:-1]
            if currang_s[-1] == '.':
                currang_s = currang_s[:-1]

        filename = infile1 + currang_s + infile2 + BW
        if os.path.isfile(filename + '.pdf'):
            names.append(filename)
else:  # Create list of files, for non-theta patterns
    if fixed and re.search("\_[0-9]+$", infile) not None:
        m = re.search("\_[0-9]+$", infile)
        fixed_dist = m.group(0)
        infile = infile[:-len(fixed_dist)]

    for i in range(0, maxval):
        if maxval != 1:
            filename = infile + str(i + start) + fixed_dist + BW
        else:
            filename = infile + fixed_dist + BW

        if os.path.isfile(filename + '.pdf'):
            names.append(filename)

for i in range(0, len(names)):  # Convert to PNG
    args = ['convert_magick', '-density',
            str(density), names[i] + '.pdf', 'img/' + names[i] + '.png']
    Popen(args).communicate()[0]
    print names[i] + ' converted'

if gif:  # Convert to GIF
    if repeat and theta != 0:  # If theta, add more files to make a cyclic gif
        namesadd = []
        invdir = True
        effang = currang
        maxang = currang
        pos = len(names) - 1
        while effang < 2 * math.pi + startang - theta:
            if invdir:
                pos -= 1
            else:
                pos += 1

            if pos < 0:
                invdir = not invdir
                pos += 2
            elif pos > len(names) - 1:
                invdir = not invdir
                pos -= 2

            effang += theta

            namesadd.append(names[pos])
        names.extend(namesadd)

    for i in range(0, len(names)):
        names[i] = 'img/' + names[i] + '.png'

    args = ['convert_magick', '-delay', str(delay),
            '-density', str(density), '-dispose', 'background']
    args.extend(names)
    args.append('img/' + infile + fixed_dist + BW + '.gif')
    Popen(args).communicate()[0]
    print 'gif generated'
