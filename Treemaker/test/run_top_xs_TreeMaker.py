import os
import glob
import sys
from DataFormats.FWLite import Events, Handle
import ROOT
import math

from optparse import OptionParser

#### This bit allows us to run the analyzer using command-line options
parser = OptionParser()

parser.add_option('-i', '--dirs', metavar='F', type='string', action='store',
                    default='/uscms_data/d3/dfehling/NTUPLES/dataA/',
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

parser.add_option('-o', '--outfile', metavar='N', type='string', action='store',
                    default='output',
                    dest='outfile',
                    help='output file')

parser.add_option('--triggerFile', metavar='N', type='string', action='store',
                    default='',
                    dest='triggerFile',
                    help='trigger file')

parser.add_option('--pileupFile', metavar='N', type='string', action='store',
                    default='',
                    dest='pileupFile',
                    help='pileup file')

parser.add_option('--config', metavar='N', type='string', action='store',
                    default='',
                    dest='config',
                    help='leave blank to run hadronic analyzer. Use "muon" to run over muon NTUPLES')

parser.add_option('--isMC', metavar='N', type='int', action='store',
                    default = 0,
                    dest='isMC',
                    help='are we running over MC(1) or not(0)')

parser.add_option('--unfoldWeight', metavar='N', type='float', action='store',
                    default=1,
                    dest='unfoldWeight',
                    help='mc weight to be used when making the response matrix')

parser.add_option('--invMassCut', metavar='N', type='int', action='store',
                    default=-1,
                    dest='invMassCut',
                    help='invariant Mass cut. Really just -1 or 700 for full powheg to be able to stitch')

parser.add_option('--doUnfold', metavar='N', type='int', action='store',
                    default=0,
                    dest='doUnfold',
                    help='Whether or not to unfold. Needs to be off for JER/JES currently')

parser.add_option('--includeTrigger', metavar='N', type='int', action='store',
                    default=0,
                    dest='includeTrigger',
                    help='Whether to include Trigger reweighting. (1 for nominal weighting, >1 for Trigger up, <0 for Trigger down, 0 for off')

parser.add_option('--includePileup', metavar='N', type='int', action='store',
                    default=0,
                    dest='includePileup',
                    help='Whether to include Pileup reweighting. 0 for no weighting. Otherwise make sure to include the correct file for systematics.')

parser.add_option('--includePDF', metavar='N', type='int', action='store',
                    default=0,
                    dest='includePDF',
                    help='Whether to include PDF reweighting. (0 for nominal weighting, >0 for PDF up, <0 for PDF down')

parser.add_option('--useSyst', metavar='N', type='string', action='store',
                    default="",
                    dest='useSyst',
                    help='Which systematic (JERup,JERdn,JESup,JERdn) to run over. Should be blank for data')

parser.add_option('--isMCatNLO', metavar='N', type='int', action='store',
                    default=0,
                    dest='isMCatNLO',
                    help='Whether we are running on MC@NLO and need to use the generator weights')

parser.add_option('--useCondor', metavar='N', type='int', action='store',
                    default=0,
                    dest='useCondor',
                    help='Whether we are running on the grid. Affects path to files')



(options, args) = parser.parse_args()

if options.config=='':
    from top_xs_TreeMaker import *
if options.config=='muon':
    from top_xs_TreeMaker_muon import *

files = sorted(glob.glob( options.dirs + "*.root" ))

totalSection = options.totalSec
section = options.sec
'''
To make sure the last job is not overly long, we round up the integer division so each job is 1 file longer. 
This may result in the last requested job not being filled, so we exit in that case
'''
numFiles = int(math.ceil( float(len(files)) / (totalSection) ) )
limit = min((section*numFiles), len(files))

if (section-1)*numFiles>=len(files):
    exit("Last section %i would be empty, exiting" % section)

files = files[(section-1)*numFiles:limit]
#print files

events = Events (files)
ntotal = events.size()
   
analyzer = tree_maker(options.outfile, options.triggerFile, options.includeTrigger, options.pileupFile, options.includePileup, options.includePDF, 
                        options.isMC, options.unfoldWeight, options.invMassCut, options.doUnfold, options.useSyst, options.isMCatNLO, options.useCondor)

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
print 'Finished   Job {0:2.0f} {1:10.0f}/{2:10.0f} : {3:5.0f} %'.format(section, count, ntotal, float(count)/float(ntotal)*100. )
