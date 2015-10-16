import os
import glob
import sys
from DataFormats.FWLite import Events, Handle
import ROOT

from optparse import OptionParser

from top_xs_TreeMaker_muon import *

#### This bit allows us to run the analyzer using command-line options
parser = OptionParser()

parser.add_option('-i', '--dirs', metavar='F', type='string', action='store',
					default='/uscms_data/d3/dfehling/NTUPLES/June/dataA/',
					dest='dirs',
					help='Input Directories (glob format)')

parser.add_option('-N', '--Nevents', metavar='N', type='int', action='store',
					default = -1,
					dest='Nevents',
					help='numEvents')

parser.add_option('--sec', metavar='N', type='int', action='store',
					default = 1,
					dest='sec',
					help='Section number')

parser.add_option('--totalSec', metavar='N', type='int', action='store',
					default = 10,
					dest='totalSec',
					help='Total number of sections')

# parser.add_option('--seed', metavar='N', type='int', action='store',
# 					default=12345,
# 					dest='seed',
# 					help='Random Seed')

parser.add_option('-o', '--outfile', metavar='N', type='string', action='store',
					default='output',
					dest='outfile',
					help='output file')

# parser.add_option('--mistagFile', metavar='N', type='string', action='store',
# 					default='',
# 					dest='mistagFile',
# 					help='mistag file')

parser.add_option('--modMassFile', metavar='N', type='string', action='store',
					default='',
					dest='modMassFile',
					help='mod mass file')

parser.add_option('--triggerFile', metavar='N', type='string', action='store',
					default='',
					dest='triggerFile',
					help='trigger file')

(options, args) = parser.parse_args()


files = sorted(glob.glob( options.dirs + "*.root" ))
#print files


totalSection = options.totalSec
section = options.sec

numFiles = len(files) / (totalSection)

# To make sure the last block is not much bigger than the rest
while(len(files)-numFiles*(totalSection-1) >= totalSection
	and len(files)-numFiles*(totalSection-1) > numFiles):
	numFiles = numFiles+1
if numFiles == totalSection-1:
	numFiles = numFiles+1
	totalSection = totalSection-1
if section == totalSection:
	secFiles = files[(section-1)*numFiles:]
else :
	secFiles = files[(section-1)*numFiles:(section)*numFiles]

files = secFiles
# print files

events = Events (files)
ntotal = events.size()

analyzer = tree_maker(options.outfile, options.modMassFile, options.triggerFile)

count = 0
print "Start looping"
for event in events:
	count = count + 1
	if count % 10000 == 0 or count == 1:
			percentDone = float(count) / float(ntotal) * 100.0
			print 'Processing Job {0:2.0f} {1:10.0f}/{2:10.0f} : {3:5.2f} %'.format(section, count, ntotal, percentDone )
			
	error = analyzer.analyze(event)
	analyzer.reset()

	if count > options.Nevents and options.Nevents > 0: 
		break   

del analyzer
