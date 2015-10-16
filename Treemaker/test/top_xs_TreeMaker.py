import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
import math
from math import *
import sys

from Analysis.Tools.JetTools import *

class tree_maker:
	def __init__(self, outputname, modMassFileStr, triggerFileStr):
		# load all the event info:
		# self.out_info = 0
		self.name = outputname
		# self.seed = seed
		# self.makeMistag = makeMistag
		self.modMassFileStr = modMassFileStr
		self.triggerFileStr = triggerFileStr

		#General Quantities
		#N Primary Vertices
		self.npvHandle = Handle( "unsigned int" )
		self.npvLabel  = ( "jhuGen", "npv" )

		#MET Pt and Phi - why not eta?
		self.metPtHandle = Handle ( "double" )
		self.metPtLabel  = ( "jhuGen", "metpt" )
		self.metPhiHandle = Handle ( "double" )
		self.metPhiLabel  = ( "jhuGen", "metphi" )
		
		#Pruned Jet Collection
		#We need this for CSV values and subjet CSV values
		self.prunedHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
		self.prunedLabel  = ( "jhuCa8pp", "PrunedCA8CORR" )

		#Btagging CSV values
		self.CSVHandle = Handle( "std::vector<double>" )
		self.CSVLabel  = ( "jhuCa8pp", "PrunedCA8csv" )

		#Subjet btagging CSV values
		self.subjet1CSVHandle = Handle( "std::vector<double>" )
		self.subjet1CSVLabel  = ( "jhuCa8tt", "TopTaggedCA8sub0csv" )
		self.subjet2CSVHandle = Handle( "std::vector<double>" )
		self.subjet2CSVLabel  = ( "jhuCa8tt", "TopTaggedCA8sub1csv" )
		self.subjet3CSVHandle = Handle( "std::vector<double>" )
		self.subjet3CSVLabel  = ( "jhuCa8tt", "TopTaggedCA8sub2csv" )
		self.subjet4CSVHandle = Handle( "std::vector<double>" )
		self.subjet4CSVLabel  = ( "jhuCa8tt", "TopTaggedCA8sub3csv" )

		#TopTagged Jet collection
		#Eventually this is the one we will use most likely
		self.topTaggedHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
		self.topTaggedLabel  = ( "jhuCa8tt", "TopTaggedCA8CORR")
		self.topTaggedJECHandle = Handle ( "std::vector<double>" )
		self.topTaggedJECLabel = ( "jhuCa8tt", "TopTaggedCA8JEC" )

		#Number of Subjets
		self.NsubjetsHandle = Handle( "std::vector<unsigned int>" )
		self.NsubjetsLabel  = ( "jhuCa8tt", "TopTaggedCA8nsub" )

		#Min-pairwise Mass
		self.minMassHandle = Handle ( "std::vector<double>")
		self.minMassLabel  = ( "jhuCa8tt", "TopTaggedCA8topTagMinMass")

		#Top Mass
		self.topTagTopMassHandle = Handle ( "std::vector<double>")
		self.topTagTopMassLabel  = ( "jhuCa8tt", "TopTaggedCA8topTagTopMass")

		#Unpruned Jet collection
		#We need this to get N-subjettiness
		self.unprunedHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > ")
		self.unprunedLabel  = ( "jhuCa8", "UnprunedCA8CORR" )

		#N-Subjettiness(tau)
		self.t1Handle = Handle( "std::vector<double>" )
		self.t1Label  = ("jhuCa8", "UnprunedCA8tau1")
		self.t2Handle = Handle( "std::vector<double>" )
		self.t2Label  = ("jhuCa8", "UnprunedCA8tau2")
		self.t3Handle = Handle( "std::vector<double>" )
		self.t3Label  = ("jhuCa8", "UnprunedCA8tau3")
		self.t4Handle = Handle( "std::vector<double>" )
		self.t4Label  = ("jhuCa8", "UnprunedCA8tau4")

		#Trigger
		self.triggerHandle = Handle( "edm::TriggerResults" )
		self.triggerLabel  = ( "TriggerResults","","HLT" )

		self.__book__()

	def __book__(self):

		if (self.modMassFileStr == ''):
			self.doModMass = False
		else:
			self.modMassFile = ROOT.TFile(self.modMassFileStr + ".root")
			self.modMassFile.cd()
			self.modMass = ROOT.TH1F()
			self.modMass = self.modMassFile.Get("jetMassOneTag_MassWindow").Clone()
			self.modMass.SetName('modMass')
			ROOT.SetOwnership( self.modMass, False )
			self.doModMass = True
		
		if (self.triggerFileStr == ''):
			self.doTrigger = False
		else:
			self.triggerFile = ROOT.TFile(self.triggerFileStr + ".root")
			self.triggerFile.cd()
			self.trigger = self.triggerFile.Get("TRIGGER_EFF").Clone()
			self.trigger.SetName('trigger')
			ROOT.SetOwnership( self.trigger, False )
			self.doTrigger = True

		print "Booking Histograms and Trees..."
		self.f = ROOT.TFile( self.name + ".root", "recreate" )
		self.f.cd()
		self.treeVars = ROOT.TTree('treeVars', 'treeVars')

		self.run = array('i', [-1])
		self.event = array('l', [-1])
		self.lumi = array('i', [-1])
		
		self.index = array('i', [-1])
		self.trigWt = array('f', [-1.0])

		self.npv = array('i', [-1])
		self.MET = array('f', [-1.0])
		self.jet1pt = array('f', [-1.0])
		self.jet2pt = array('f', [-1.0])
		self.jet1eta = array('f', [-10.0])
		self.jet2eta = array('f', [-10.0])
		self.jet1phi = array('f', [-10.0])
		self.jet2phi = array('f', [-10.0])
		self.jet1mass = array('f', [-1.0])
		self.jet2mass = array('f', [-1.0])
		self.jet1topMass = array('f', [-1.0])
		self.jet2topMass = array('f', [-1.0])
		# self.jet1topMassRaw = array('f', [-1.0])
		# self.jet2topMassRaw = array('f', [-1.0])
		self.jet1csv = array('f', [-1.0])
		self.jet2csv = array('f', [-1.0])
		self.jet1maxSubjetCSV = array('f', [-1.0])
		self.jet2maxSubjetCSV = array('f', [-1.0])
		self.jet1tau1 = array('f', [-1.0])
		self.jet1tau2 = array('f', [-1.0])
		self.jet1tau3 = array('f', [-1.0])
		self.jet1tau4 = array('f', [-1.0])
		self.jet2tau1 = array('f', [-1.0])
		self.jet2tau2 = array('f', [-1.0])
		self.jet2tau3 = array('f', [-1.0])
		self.jet2tau4 = array('f', [-1.0])
		self.deltaY = array('f', [-10.0])
		self.deltaPhi = array('f', [-10.0])

		self.jet1tau32 = array('f', [-1.0])
		self.jet2tau32 = array('f', [-1.0])
		self.jet1tau31 = array('f', [-1.0])
		self.jet2tau31 = array('f', [-1.0])
		self.jet1tau21 = array('f', [-1.0])
		self.jet2tau21 = array('f', [-1.0])
		self.jet1nSubj = array('i', [-1])
		self.jet2nSubj = array('i', [-1])
		self.jet1minMass = array('f', [-1.0])
		self.jet2minMass = array('f', [-1.0])

		# self.cutflow = array('f', [-1.0])

		# self.x = array('f',[-1.0])
		self.jet1bTagged = array('i', [-1])
		self.jet2bTagged = array('i', [-1])
		self.jet1topTagged = array('i', [-1])
		self.jet2topTagged = array('i', [-1])
		# self.taggedPt = array('f', [-1.0])
		# self.taggedMass = array('f', [-1.0])

		# self.mistagWt = array('f', [-1.0])
		# self.mistagWt1B = array('f', [-1.0])
		# self.mistagWt2B = array('f', [-1.0])

		self.htSum = array('f', [-1.0])
		self.triggerEff = array('f', [1.0])

		self.jet1JEC = array('f', [-1.0])
		self.jet2JEC = array('f', [-1.0])


		self.treeVars.Branch('run', self.run, 'run/I')
		self.treeVars.Branch('event', self.event, 'event/L')
		self.treeVars.Branch('lumi', self.lumi, 'lumi/I')

		self.treeVars.Branch('npv', self.npv, 'npv/I')
		self.treeVars.Branch('index', self.index, 'index/I')
		self.treeVars.Branch('trigWt', self.trigWt, 'trigWt/F')
		self.treeVars.Branch('MET', self.MET, 'MET/F')

		self.treeVars.Branch('jet1pt', self.jet1pt, 'jet1pt/F')
		self.treeVars.Branch('jet2pt', self.jet2pt, 'jet2pt/F')
		self.treeVars.Branch('jet1eta', self.jet1eta, 'jet1eta/F')
		self.treeVars.Branch('jet2eta', self.jet2eta, 'jet2eta/F')
		self.treeVars.Branch('jet1phi', self.jet1phi, 'jet1phi/F')
		self.treeVars.Branch('jet2phi', self.jet2phi, 'jet2phi/F')
		self.treeVars.Branch('jet1mass', self.jet1mass, 'jet1mass/F')
		self.treeVars.Branch('jet2mass', self.jet2mass, 'jet2mass/F')
		self.treeVars.Branch('jet1topMass', self.jet1topMass, 'jet1topMass/F')
		self.treeVars.Branch('jet2topMass', self.jet2topMass, 'jet2topMass/F')
		# self.treeVars.Branch('jet1topMassRaw', self.jet1topMassRaw, 'jet1topMassRaw/F')
		# self.treeVars.Branch('jet2topMassRaw', self.jet2topMassRaw, 'jet2topMassRaw/F')
		self.treeVars.Branch('jet1csv', self.jet1csv, 'jet1csv/F')
		self.treeVars.Branch('jet2csv', self.jet2csv, 'jet2csv/F')
		self.treeVars.Branch('jet1maxSubjetCSV', self.jet1maxSubjetCSV, 'jet1maxSubjetCSV/F')
		self.treeVars.Branch('jet2maxSubjetCSV', self.jet2maxSubjetCSV, 'jet2maxSubjetCSV/F')
		self.treeVars.Branch('jet1tau1', self.jet1tau1, 'jet1tau1/F')
		self.treeVars.Branch('jet1tau2', self.jet1tau2, 'jet1tau2/F')
		self.treeVars.Branch('jet1tau3', self.jet1tau3, 'jet1tau3/F')
		self.treeVars.Branch('jet1tau4', self.jet1tau4, 'jet1tau4/F')
		self.treeVars.Branch('jet2tau1', self.jet2tau1, 'jet2tau1/F')
		self.treeVars.Branch('jet2tau2', self.jet2tau2, 'jet2tau2/F')
		self.treeVars.Branch('jet2tau3', self.jet2tau3, 'jet2tau3/F')
		self.treeVars.Branch('jet2tau4', self.jet2tau4, 'jet2tau4/F')
		self.treeVars.Branch('jet1tau32', self.jet1tau32, 'jet1tau32/F')
		self.treeVars.Branch('jet2tau32', self.jet2tau32, 'jet2tau32/F')
		self.treeVars.Branch('jet1tau31', self.jet1tau31, 'jet1tau31/F')
		self.treeVars.Branch('jet2tau31', self.jet2tau31, 'jet2tau31/F')
		self.treeVars.Branch('jet1tau21', self.jet1tau21, 'jet1tau21/F')
		self.treeVars.Branch('jet2tau21', self.jet2tau21, 'jet2tau21/F')
		self.treeVars.Branch('jet1nSubj', self.jet1nSubj, 'jet1nSubj/I')
		self.treeVars.Branch('jet2nSubj', self.jet2nSubj, 'jet2nSubj/I')
		self.treeVars.Branch('jet1minMass', self.jet1minMass, 'jet1minMass/F')
		self.treeVars.Branch('jet2minMass', self.jet2minMass, 'jet2minMass/F')

		self.treeVars.Branch('deltaY', self.deltaY, 'deltaY/F')
		self.treeVars.Branch('deltaPhi', self.deltaPhi, 'deltaPhi/F')
		# self.treeVars.Branch('cutflow', self.cutflow, 'cutflow/F')

		# self.treeVars.Branch('x', self.x, 'x/F')
		self.treeVars.Branch('jet1bTagged', self.jet1bTagged, 'jet1bTagged/I')
		self.treeVars.Branch('jet2bTagged', self.jet2bTagged, 'jet2bTagged/I')
		self.treeVars.Branch('jet1topTagged', self.jet1topTagged, 'jet1topTagged/I')
		self.treeVars.Branch('jet2topTagged', self.jet2topTagged, 'jet2topTagged/I')
		# self.treeVars.Branch('taggedPt', self.taggedPt, 'taggedPt/F')
		# self.treeVars.Branch('taggedMass', self.taggedMass, 'taggedMass/F')

		# self.treeVars.Branch('mistagWt', self.mistagWt, 'mistagWt/F')
		# self.treeVars.Branch('mistagWt1B', self.mistagWt1B, 'mistagWt1B/F')
		# self.treeVars.Branch('mistagWt2B', self.mistagWt2B, 'mistagWt2B/F')

		self.treeVars.Branch('htSum', self.htSum, 'htSum/F')
		self.treeVars.Branch('triggerEff', self.triggerEff, 'triggerEff/F')

		self.treeVars.Branch('jet1JEC', self.jet1JEC, 'jet1JEC/F')
		self.treeVars.Branch('jet2JEC', self.jet2JEC, 'jet2JEC/F')

		self.invarmass = array('f', [-1.0])
		self.jetangle = array('f', [-10.0])
		self.treeVars.Branch('invariant_mass', self.invarmass, 'invarmass/F')
		self.treeVars.Branch('angle_between_jets', self.jetangle, 'jetangle/F')


		# self.topTagPt         = ROOT.TH1D("topTagPt",         "Top Tag Pt",       400,  0,  2000 )
		# self.topProbePt       = ROOT.TH1D("topProbePt",       "Top Probe Pt",     400,  0,  2000 )
		# self.testTagPt        = ROOT.TH1D("testTagPt",        "Top Tag Pt",       400,  0,  2000 )
		# self.testProbePt      = ROOT.TH1D("testProbePt",      "Top Probe Pt",     400,  0,  2000 )
		# self.lowmMinTagWP1Pt  = ROOT.TH1D("lowmMinTagPtWP1",  "Top Tag Pt WP1",   400,  0,  2000 )
		# self.lowmMinTagWP2Pt  = ROOT.TH1D("lowmMinTagPtWP2",  "Top Tag Pt WP2",   400,  0,  2000 )
		# self.lowmMinTagWP3Pt  = ROOT.TH1D("lowmMinTagPtWP3",  "Top Tag Pt WP3",   400,  0,  2000 )
		# self.lowmMinTagWP4Pt  = ROOT.TH1D("lowmMinTagPtWP4",  "Top Tag Pt WP4",   400,  0,  2000 )
		# self.lowmMinTagWP5Pt  = ROOT.TH1D("lowmMinTagPtWP5",  "Top Tag Pt WP5",   400,  0,  2000 )
		# self.lowmMinTagTau31WP1Pt = ROOT.TH1D("lowmMinTagTau31WP1", "Top Tag Pt Tau31 WP1", 400, 0, 2000 )
		# self.lowmMinTagTau31WP2Pt = ROOT.TH1D("lowmMinTagTau31WP2", "Top Tag Pt Tau31 WP2", 400, 0, 2000 )
		# self.lowmMinTagPt     = ROOT.TH1D("lowmMinTagPt",     "Top Tag Pt",       400,  0,  2000 )
		# self.lowmMinProbePt   = ROOT.TH1D("lowmMinProbePt",   "Top Probe Pt",     400,  0,  2000 )
		# self.lowmMinBTagPt    = ROOT.TH1D("lowmMinBTagPt",    "Top Tag Pt",       400,  0,  2000 )
		# self.lowmMinBProbePt  = ROOT.TH1D("lowmMinBProbePt",  "Top Probe Pt",     400,  0,  2000 )
		# self.lowmMin2BTagPt   = ROOT.TH1D("lowmMin2BTagPt",   "Top Tag Pt",       400,  0,  2000 )
		# self.lowmMin2BProbePt = ROOT.TH1D("lowmMin2BProbePt", "Top Probe Pt",     400,  0,  2000 )
		# self.antiBTagPt       = ROOT.TH1D("antiBTagPt",       "Top Tag Pt",       400,  0,  2000 ) 
		# self.antiBProbePt     = ROOT.TH1D("antiBProbePt",     "Top Probe Pt",     400,  0,  2000 )

		# self.lowmMinTagTau31Tau    = ROOT.TH1D("lowmMinTagTau31Tau",    "Top Tag Tau31",     20, 0, 1 )
		# self.lowmMinTagTau31WP1Tau = ROOT.TH1D("lowmMinTagTau31WP1Tau", "Top Tag Tau31 WP1", 20, 0, 1 )
		# self.lowmMinTagTau31WP2Tau = ROOT.TH1D("lowmMinTagTau31WP2Tau", "Top Tag Tau31 WP2", 20, 0, 1 )

		self.cutflow = ROOT. TH1D("cutflow", "cutfow", 10, 0, 10 )
		self.cutflow.Sumw2()

		# self.topTagPt.Sumw2()
		# self.topProbePt.Sumw2()
		# self.testTagPt.Sumw2()
		# self.testProbePt.Sumw2()
		# self.lowmMinTagPt.Sumw2()
		# self.lowmMinTagWP1Pt.Sumw2()
		# self.lowmMinTagWP2Pt.Sumw2()
		# self.lowmMinTagWP3Pt.Sumw2()
		# self.lowmMinTagWP4Pt.Sumw2()
		# self.lowmMinTagWP5Pt.Sumw2()
		# self.lowmMinTagTau31WP1Pt.Sumw2()
		# self.lowmMinTagTau31WP2Pt.Sumw2()
		# self.lowmMinProbePt.Sumw2()
		# self.lowmMinBTagPt.Sumw2()
		# self.lowmMinBProbePt.Sumw2()
		# self.lowmMin2BTagPt.Sumw2()
		# self.lowmMin2BProbePt.Sumw2()
		# self.antiBTagPt.Sumw2()
		# self.antiBProbePt.Sumw2()
		# self.lowmMinTagTau31Tau.Sumw2()
		# self.lowmMinTagTau31WP1Tau.Sumw2()
		# self.lowmMinTagTau31WP2Tau.Sumw2()

	def analyze(self, event):
		# // To get the trigger names
		# edm::Handle< edm::TriggerResults > h_trigresults
		# edm::InputTag triggerResultsSrc_("TriggerResults", "", "HLT")
		# iEvent.getByLabel( triggerResultsSrc_, h_trigresults )
		# const edm::TriggerResults* thing = h_trigresults.product()
		# edm::TriggerNames const & trig_names = iEvent.triggerNames(*thing)
		# std::vector<std::string> const & trig_strings = trig_names.triggerNames()

		event.getByLabel (self.triggerLabel, self.triggerHandle)

		self.run[0] = event.object().id().run()
		self.event[0] = event.object().id().event()
		self.lumi[0] = event.object().id().luminosityBlock()

		event.getByLabel (self.npvLabel, self.npvHandle)
		event.getByLabel (self.metPtLabel, self.metPtHandle)
		event.getByLabel (self.metPhiLabel, self.metPhiHandle)

		event.getByLabel (self.prunedLabel, self.prunedHandle)
		event.getByLabel (self.unprunedLabel, self.unprunedHandle)
		event.getByLabel (self.topTaggedLabel, self.topTaggedHandle)

		event.getByLabel (self.CSVLabel, self.CSVHandle)

		event.getByLabel (self.subjet1CSVLabel, self.subjet1CSVHandle)
		event.getByLabel (self.subjet2CSVLabel, self.subjet2CSVHandle)
		event.getByLabel (self.subjet3CSVLabel, self.subjet3CSVHandle)
		event.getByLabel (self.subjet4CSVLabel, self.subjet4CSVHandle)

		event.getByLabel (self.t1Label, self.t1Handle)
		event.getByLabel (self.t2Label, self.t2Handle)
		event.getByLabel (self.t3Label, self.t3Handle)
		event.getByLabel (self.t4Label, self.t4Handle)

		event.getByLabel (self.NsubjetsLabel, self.NsubjetsHandle)

		event.getByLabel (self.minMassLabel, self.minMassHandle)
		event.getByLabel (self.topTagTopMassLabel, self.topTagTopMassHandle)

		event.getByLabel (self.topTaggedJECLabel, self.topTaggedJECHandle)	

		npv = self.npvHandle.product()[0]
		metPt = self.metPtHandle.product()[0]
		metPhi = self.metPhiHandle.product()[0]

		CSVVals = self.CSVHandle.product()

		subjet1CSV = self.subjet1CSVHandle.product()
		subjet2CSV = self.subjet2CSVHandle.product()
		subjet3CSV = self.subjet3CSVHandle.product()
		subjet4CSV = self.subjet4CSVHandle.product()

		Tau1  =  self.t1Handle.product()
		Tau2  =  self.t2Handle.product()
		Tau3  =  self.t3Handle.product()
		Tau4  =  self.t4Handle.product()

		nSubjets = self.NsubjetsHandle.product()

		minMass = self.minMassHandle.product()
		topTagTopMass = self.topTagTopMassHandle.product()

		topTaggedJEC = self.topTaggedJECHandle.product()

		HTsum = 0

		#Total events seen
		self.index[0]=0;
		self.cutflow.Fill(1);

		if (self.prunedHandle.isValid() and self.unprunedHandle.isValid() and self.topTaggedHandle.isValid()) :
			unpj = self.unprunedHandle.product()
			pj = self.prunedHandle.product()
			ttpj = self.topTaggedHandle.product()

			#Events with valid handle
			self.cutflow.Fill(2);

			## approximately 3% of the time we have more pj than unpj. Why? Shouldn't this be impossible?
			## The pruning can separate 2 subjets into standalone jets so this makes sense
			# if len(unpj)<len(pj):
			# 	self.numbtag=self.numbtag+1

			#Two jets, including top tagged:
			if len(pj) < 2 or len(unpj) < 2 or len(ttpj) < 2:
				return
		
			#Events with 2 or more valid jets
			self.cutflow.Fill(3);

			#Make sure there's at least a top condidiate by checking the pT
			nTopCand = 0
			for i in range(0,len(ttpj) ) :
				if ttpj[i].pt() > 400. :
					nTopCand = nTopCand + 1
					HTsum += ttpj[i].pt()
			if nTopCand < 2:
				return
			
			#Events with 2 top candidates (pt > 400)
			self.cutflow.Fill(4);

			self.htSum[0] = HTsum
			self.MET[0] = metPt
			self.npv[0] = npv

			#Reorder to make sure the highest pT is first. This is after any JEC. Is this correct?
			#pt_sorted_jets = ttpj
			pt_sorted_jets = ReorderByPt(ttpj)

			#We need these for the substructure
			ca1 = ROOT.TLorentzVector()
			ca2 = ROOT.TLorentzVector()

			#We create the CA jet to match and fill the jet parameters			
			ca1.SetPtEtaPhiM(pt_sorted_jets[0].Pt(), pt_sorted_jets[0].Eta(), pt_sorted_jets[0].Phi(), pt_sorted_jets[0].M())
			ca2.SetPtEtaPhiM(pt_sorted_jets[1].Pt(), pt_sorted_jets[1].Eta(), pt_sorted_jets[1].Phi(), pt_sorted_jets[1].M())

			#Match unpruned jets with pruned - so we have both subjet btagging and nsubjettiness
			#This returns the index of the jet in the first collection that matches within dr = 0.4 to the jet of the second argument
			jet1matchIndex = MatchCol(unpj, ca1)
			jet2matchIndex = MatchCol(unpj, ca2)
			#This is temporary so we get the correct pj index for CSV and subjet CSV
			jet1matchIndex_pj = MatchCol(pj, ca1)
			jet2matchIndex_pj = MatchCol(pj, ca2)
			#Again temporary to match toptagged for Nsubjets, minMass, and topMass
			jet2matchIndex_tt = MatchCol(ttpj, ca2)
			jet1matchIndex_tt = MatchCol(ttpj, ca1)

			#Make sure matching is correct for jet1
			if jet1matchIndex==-1 or jet1matchIndex_pj==-1 or jet1matchIndex_tt ==-1:
				return
			
			#Make sure matching is correct for jet2
			if jet2matchIndex==-1 or jet2matchIndex_pj==-1 or jet2matchIndex_tt ==-1:
				return
			
			#Events with valid ca-other jet matches
			self.cutflow.Fill(5)

			#Make sure the subjets are filled
			if len(nSubjets)==0:
				return

			#Events with filled subjets
			self.cutflow.Fill(6)

			#Nsubjettiness
			if Tau2[jet1matchIndex] == 0 or Tau2[jet2matchIndex] == 0:
				return
			if Tau1[jet1matchIndex] == 0 or Tau1[jet2matchIndex] == 0:
				return

			self.jet1tau1[0] = Tau1[jet1matchIndex]
			self.jet1tau2[0] = Tau2[jet1matchIndex]
			self.jet1tau3[0] = Tau3[jet1matchIndex]
			self.jet1tau4[0] = Tau4[jet1matchIndex]
			self.jet2tau1[0] = Tau1[jet2matchIndex]
			self.jet2tau2[0] = Tau2[jet2matchIndex]
			self.jet2tau3[0] = Tau3[jet2matchIndex]
			self.jet2tau4[0] = Tau4[jet2matchIndex]

			#Events with nonzero tau2 and tau1
			self.cutflow.Fill(7)
			
			self.jet1pt[0] = ca1.Pt()
			self.jet2pt[0] = ca2.Pt()
			self.jet1eta[0] = ca1.Eta()
			self.jet2eta[0] = ca2.Eta()
			self.jet1phi[0] = ca1.Phi()
			self.jet2phi[0] = ca2.Phi()
			self.jet1mass[0] = ca1.M()
			self.jet2mass[0] = ca2.M()

			#Top Tagging Parameters
			self.jet1nSubj[0] = nSubjets[jet1matchIndex_tt]
			self.jet2nSubj[0] = nSubjets[jet2matchIndex_tt]
			self.jet1minMass[0] = minMass[jet1matchIndex_tt]
			self.jet2minMass[0] = minMass[jet2matchIndex_tt]
			#### This should rightfully be jet1mass and jet2mass ####
			self.jet1topMass[0] = topTagTopMass[jet1matchIndex_tt]
			self.jet2topMass[0] = topTagTopMass[jet2matchIndex_tt]
			# self.jet1topMassRaw[0] = topTagTopMass[jet1matchIndex_tt]
			# self.jet2topMassRaw[0] = topTagTopMass[jet2matchIndex_tt]

			self.jet1JEC[0] = topTaggedJEC[jet1matchIndex_tt]
			self.jet2JEC[0] = topTaggedJEC[jet2matchIndex_tt]
			
			#Invariant Mass			
			self.invarmass[0] = (ca1+ca2).M()

			#Angular Parameters
			self.jetangle[0] = ca1.DeltaR(ca2)
			self.deltaY[0] = (ca1.Rapidity() - ca2.Rapidity())
			self.deltaPhi[0] = ca1.DeltaPhi(ca2)
			
			#Nsubjettiness
			self.jet1tau32[0] = Tau3[jet1matchIndex] / Tau2[jet1matchIndex]
			self.jet2tau32[0] = Tau3[jet2matchIndex] / Tau2[jet2matchIndex]

			self.jet1tau31[0] = Tau3[jet1matchIndex] / Tau1[jet1matchIndex]
			self.jet2tau31[0] = Tau3[jet2matchIndex] / Tau1[jet2matchIndex]

			self.jet1tau21[0] = Tau2[jet1matchIndex] / Tau1[jet1matchIndex]
			self.jet2tau21[0] = Tau2[jet2matchIndex] / Tau1[jet2matchIndex]

			#Fill the btagging information
			self.jet1csv[0] = CSVVals[jet1matchIndex_pj]
			self.jet2csv[0] = CSVVals[jet2matchIndex_pj]
			
			jet1subjetCSVs = []
			jet1subjetCSVs.append(subjet1CSV[jet1matchIndex_tt])
			jet1subjetCSVs.append(subjet2CSV[jet1matchIndex_tt])
			jet1subjetCSVs.append(subjet3CSV[jet1matchIndex_tt])
			jet1subjetCSVs.append(subjet4CSV[jet1matchIndex_tt])

			jet2subjetCSVs = []
			jet2subjetCSVs.append(subjet1CSV[jet2matchIndex_tt])
			jet2subjetCSVs.append(subjet2CSV[jet2matchIndex_tt])
			jet2subjetCSVs.append(subjet3CSV[jet2matchIndex_tt])
			jet2subjetCSVs.append(subjet4CSV[jet2matchIndex_tt])

			self.jet1maxSubjetCSV[0] = max(jet1subjetCSVs)
			self.jet2maxSubjetCSV[0] = max(jet2subjetCSVs)

			self.jet1bTagged[0] = self.jet1maxSubjetCSV[0] > 0.679
			self.jet2bTagged[0] = self.jet2maxSubjetCSV[0] > 0.679

			self.jet1topTagged[0] = self.jet1topMass[0] > 140.0 and self.jet1topMass[0] < 250.0 and self.jet1minMass[0] > 50.0 and self.jet1nSubj[0] > 2 and self.jet1pt[0] > 400.
			self.jet2topTagged[0] = self.jet2topMass[0] > 140.0 and self.jet2topMass[0] < 250.0 and self.jet2minMass[0] > 50.0 and self.jet2nSubj[0] > 2 and self.jet2pt[0] > 400.

			if self.jet1topTagged[0] and self.jet2topTagged[0]:	
				self.index[0] = 1
				self.cutflow.Fill(7)

			self.treeVars.Fill()

			# if self.seed < 0:
			# 	temp_x = ROOT.TRandom3(0)
			# else:
			# 	#This is not truly correct, but it is better than what we had before
			# 	temp_x = ROOT.TRandom3(self.run[0]+self.event[0]+self.lumi[0]+self.seed)
			# x = temp_x.Uniform(1.0)

			# self.x[0] = x

			# Signal
				# self.cutflow[0] = 4.0
				# if self.x[0] < 0.5:
				# 	self.taggedPt[0] = self.jet1pt[0]
				# 	self.taggedMass[0] = self.jet1topMass[0];
				# else:
				# 	self.taggedPt[0] = self.jet2pt[0]
				# 	self.taggedMass[0] = self.jet2topMass[0];

				#Signal Region

			# if self.doTrigger:
			# 	triggerBin = self.trigger.FindBin(HTsum)
			# 	if self.trigger.IsBinOverflow(triggerBin):
			# 		triggerBin = self.trigger.GetNbinsX()
			# 	elif triggerBin == 0:
			# 		triggerBin = 1
				
			# 	self.triggerEff = self.trigger.GetBinContent( triggerBin )

			# if x < 0.5:
			# 	if self.jet1topTagged[0]:
					
			# 		if self.doModMass:
			# 			self.jet2topMass[0] = self.modMass.GetRandom()

			# 		self.index[0] = 2

			# 		self.taggedPt[0] = self.jet1pt[0]
			# 		self.taggedMass[0] = self.jet1topMass[0]

			# 		#Jet 1 top tagged in background estimate
			# 		self.cutflow.Fill(8)
			# 		self.cutflow.Fill(10)
			# 		self.treeVars.Fill()

			# # if x >= 0.5:
			# else:
			# 	if self.jet2topTagged[0]:

			# 		if self.doModMass:
			# 			self.jet1topMass[0] = self.modMass.GetRandom()

			# 		self.index[0] = 2

			# 		self.taggedPt[0] = self.jet2pt[0]
			# 		self.taggedMass[0] = self.jet2topMass[0]

			# 		#Jet 2 top tagged in background estimate
			# 		self.cutflow.Fill(9)
			# 		self.cutflow.Fill(10)
			# 		self.treeVars.Fill()

			# self.makeMistag()
			return

	def makeMistag(self):


		bTag1 = self.jet1bTagged[0]
		bTag2 = self.jet2bTagged[0]

		bTagLoose1 = self.jet1maxSubjetCSV[0] > 0.244
		bTagLoose2 = self.jet2maxSubjetCSV[0] > 0.244
		
		#Cuts - courtesy of JP
		#Why don't we also require nSubj>2?
		failMinMass1    = self.jet1topMass[0] > 140 and self.jet1topMass[0] < 250 and self.jet1minMass[0] < 50  
		failMinMass2    = self.jet2topMass[0] > 140 and self.jet2topMass[0] < 250 and self.jet2minMass[0] < 50
		failMinMassLow1 = self.jet1topMass[0] > 140 and self.jet1topMass[0] < 250 and self.jet1minMass[0] < 30
		failMinMassLow2 = self.jet2topMass[0] > 140 and self.jet2topMass[0] < 250 and self.jet2minMass[0] < 30

		fail2Nsub1 = self.jet1topMass[0] > 140 and self.jet1topMass[0] < 250 and self.jet1nSubj[0] < 2
		fail2Nsub2 = self.jet2topMass[0] > 140 and self.jet2topMass[0] < 250 and self.jet2nSubj[0] < 2
		fail3Nsub1 = self.jet1topMass[0] > 140 and self.jet1topMass[0] < 250 and self.jet1nSubj[0] < 3
		fail3Nsub2 = self.jet2topMass[0] > 140 and self.jet2topMass[0] < 250 and self.jet2nSubj[0] < 3

		topTag1WP1 = self.jet1topTagged[0] and self.jet1tau32[0] < 0.7  and bTagLoose1
		topTag2WP1 = self.jet2topTagged[0] and self.jet2tau32[0] < 0.7  and bTagLoose2
		topTag1WP2 = self.jet1topTagged[0] and self.jet1tau32[0] < 0.6  and bTagLoose1
		topTag2WP2 = self.jet2topTagged[0] and self.jet2tau32[0] < 0.6  and bTagLoose2
		topTag1WP3 = self.jet1topTagged[0] and self.jet1tau32[0] < 0.55 and bTag1
		topTag2WP3 = self.jet2topTagged[0] and self.jet2tau32[0] < 0.55 and bTag2
		topTag1WP4 = self.jet1topTagged[0] and self.jet1tau32[0] < 0.5  and bTag1 and self.jet1minMass[0] > 65
		topTag2WP4 = self.jet2topTagged[0] and self.jet2tau32[0] < 0.5  and bTag2 and self.jet2minMass[0] > 65
		topTag1WP5 = self.jet1topTagged[0] and self.jet1tau32[0] < 0.4  and bTag1 and self.jet1minMass[0] > 55
		topTag2WP5 = self.jet2topTagged[0] and self.jet2tau32[0] < 0.4  and bTag2 and self.jet2minMass[0] > 55

		topTag1Tau31WP1 = self.jet1topTagged[0] and self.jet1tau31[0] < 0.55 and bTag1
		topTag2Tau31WP1 = self.jet2topTagged[0] and self.jet2tau31[0] < 0.55 and bTag2
		topTag1Tau31WP2 = self.jet1topTagged[0] and self.jet1tau31[0] < 0.40 and bTag1
		topTag2Tau31WP2 = self.jet2topTagged[0] and self.jet2tau31[0] < 0.40 and bTag2

		weight = 1.0

		if self.x[0] < 0.5 :
			if not self.jet1topTagged[0]:
				self.topProbePt.Fill( self.jet2pt[0], weight )    
				if self.jet2topTagged[0]:
					self.topTagPt.Fill( self.jet2pt[0], weight )
			if failMinMass1:
				self.testProbePt.Fill( self.jet2pt[0], weight )
				if self.jet2topTagged[0]:
					self.testTagPt.Fill( self.jet2pt[0], weight )
			if failMinMassLow1:
				self.lowmMinProbePt.Fill( self.jet2pt[0], weight )
				if topTag2WP1:
					self.lowmMinTagWP1Pt.Fill( self.jet2pt[0], weight)
				if topTag2WP2:
					self.lowmMinTagWP2Pt.Fill( self.jet2pt[0], weight)
				if topTag2WP3:
					self.lowmMinTagWP3Pt.Fill( self.jet2pt[0], weight)
				if topTag2WP4:
					self.lowmMinTagWP4Pt.Fill( self.jet2pt[0], weight)
				if topTag2WP5:
					self.lowmMinTagWP5Pt.Fill( self.jet2pt[0], weight)
				if self.jet2topTagged[0]:
					self.lowmMinTagPt.Fill( self.jet2pt[0], weight )
					self.lowmMinTagTau31Tau.Fill( self.jet2tau31[0], weight)
				#I have a feeling that this is being done wrong
				# if not self.jet1bTagged[0]:
				# 	self.lowmMinBProbePt.Fill( self.jet2pt[0], weight )
				# 	if self.jet2bTagged[0] and self.jet2topTagged[0]:
				# 		self.lowmMinBTagPt.Fill( self.jet2pt[0], weight ) 
				if not self.jet1bTagged[0]:
					self.lowmMinBProbePt.Fill( self.jet2pt[0], weight )
					if self.jet2bTagged[0] and self.jet2topTagged[0]:
						self.lowmMinBTagPt.Fill( self.jet2pt[0], weight ) 
				if self.jet2bTagged[0] and self.jet1bTagged[0]:
					self.lowmMin2BProbePt.Fill( self.jet2pt[0], weight )
					if self.jet2topTagged[0]:
						self.lowmMin2BTagPt.Fill( self.jet2pt[0], weight )
				if topTag2Tau31WP1:
					self.lowmMinTagTau31WP1Pt.Fill( self.jet2pt[0], weight )
					self.lowmMinTagTau31WP1Tau.Fill( self.jet2tau31[0], weight )
				if topTag2Tau31WP2:
					self.lowmMinTagTau31WP2Pt.Fill( self.jet2pt[0], weight )
					self.lowmMinTagTau31WP2Tau.Fill( self.jet2tau31[0], weight )
			if self.jet1topTagged[0] and self.jet2topTagged[0] and not self.jet1bTagged[0]:
				self.antiBProbePt.Fill( self.jet2pt[0], weight )
				if self.jet2bTagged[0]:
					self.antiBTagPt.Fill( self.jet2pt[0], weight)

		if self.x[0] >= 0.5 :
			if not self.jet2topTagged[0]:
				self.topProbePt.Fill( self.jet1pt[0], weight )
				if self.jet1topTagged[0]:
					self.topTagPt.Fill( self.jet1pt[0], weight )
			if failMinMass2:
				self.testProbePt.Fill( self.jet1pt[0], weight )
				if self.jet1topTagged[0]:
					self.testTagPt.Fill( self.jet1pt[0], weight )
			if failMinMassLow2:
				self.lowmMinProbePt.Fill( self.jet1pt[0], weight )
				if topTag1WP1:
					self.lowmMinTagWP1Pt.Fill( self.jet1pt[0], weight)
				if topTag1WP2:
					self.lowmMinTagWP2Pt.Fill( self.jet1pt[0], weight)
				if topTag1WP3:
					self.lowmMinTagWP3Pt.Fill( self.jet1pt[0], weight)
				if topTag1WP4:
					self.lowmMinTagWP4Pt.Fill( self.jet1pt[0], weight)
				if topTag1WP5:
					self.lowmMinTagWP5Pt.Fill( self.jet1pt[0], weight)
				if self.jet1topTagged[0]:
					self.lowmMinTagPt.Fill( self.jet1pt[0], weight )
					self.lowmMinTagTau31Tau.Fill( self.jet1tau31[0], weight)
				#I have a feeling that this is being done wrong
				# if not self.jet2bTagged[0]:
				# 	self.lowmMinBProbePt.Fill( self.jet1pt[0], weight )
				# 	if self.jet1bTagged[0] and self.jet1topTagged[0]:
				# 		self.lowmMinBTagPt.Fill( self.jet1pt[0], weight ) 
				if not self.jet2bTagged[0]:
					self.lowmMinBProbePt.Fill( self.jet1pt[0], weight )
					if self.jet1bTagged[0] and self.jet1topTagged[0]:
						self.lowmMinBTagPt.Fill( self.jet1pt[0], weight ) 

				if self.jet2bTagged[0] and self.jet1bTagged[0]:
					self.lowmMin2BProbePt.Fill( self.jet1pt[0], weight )
					if self.jet1topTagged[0]:
						self.lowmMin2BTagPt.Fill( self.jet1pt[0], weight )
				if topTag1Tau31WP1:
					self.lowmMinTagTau31WP1Pt.Fill( self.jet1pt[0], weight )
					self.lowmMinTagTau31WP1Tau.Fill( self.jet1tau31[0], weight )
				if topTag1Tau31WP2:
					self.lowmMinTagTau31WP2Pt.Fill( self.jet1pt[0], weight )
					self.lowmMinTagTau31WP1Tau.Fill( self.jet1tau31[0], weight )
			if self.jet1topTagged[0] and self.jet2topTagged[0] and not self.jet2bTagged[0]:
				self.antiBProbePt.Fill( self.jet1pt[0], weight )
				if self.jet1bTagged[0]:
					self.antiBTagPt.Fill( self.jet1pt[0], weight)      


	def reset(self):

		self.run[0] = -1
		self.event[0] = -1
		self.lumi[0] = -1

		self.npv[0] = -1
		self.index[0] = -1
		self.trigWt[0] = -1.0
		self.MET[0] = -1.0

		self.jet1pt[0] = -1.0
		self.jet2pt[0] = -1.0
		self.jet1eta[0] = -10.0
		self.jet2eta[0] = -10.0
		self.jet1phi[0] = -10.0
		self.jet2phi[0] = -10.0
		self.jet1mass[0] = -1.0
		self.jet2mass[0] = -1.0
		self.jet1topMass[0] = -1.0
		self.jet2topMass[0] = -1.0
		# self.jet1topMassRaw[0] = -1.0
		# self.jet2topMassRaw[0] = -1.0
		self.jet1csv[0] = -1.0
		self.jet2csv[0] = -1.0
		self.jet1maxSubjetCSV[0] = -1.0
		self.jet2maxSubjetCSV[0] = -1.0
		self.jet1tau1[0] = -1.0
		self.jet1tau2[0] = -1.0
		self.jet1tau3[0] = -1.0
		self.jet1tau4[0] = -1.0
		self.jet2tau1[0] = -1.0
		self.jet2tau2[0] = -1.0
		self.jet2tau3[0] = -1.0
		self.jet2tau4[0] = -1.0
		self.jet1tau32[0] = -1.0
		self.jet2tau32[0] = -1.0
		self.jet1tau31[0] = -1.0
		self.jet2tau31[0] = -1.0
		self.jet1tau21[0] = -1.0
		self.jet2tau21[0] = -1.0
		self.jet1nSubj[0] = -1
		self.jet2nSubj[0] = -1
		self.jet1minMass[0] = -1.0
		self.jet2minMass[0] = -1.0
		self.deltaY[0] = -10.0
		self.deltaPhi[0] = -10.0
		# self.cutflow[0] = -1.0

		# self.x[0] = -1.0
		self.jet1bTagged[0] = -1
		self.jet2bTagged[0] = -1
		self.jet1topTagged[0] = -1
		self.jet2topTagged[0] = -1
		# self.taggedPt[0] = -1.0
		# self.taggedMass[0] = -1.0

		# self.mistagWt[0] = -1.0
		# self.mistagWt1B[0] = -1.0
		# self.mistagWt2B[0] = -1.0

		self.htSum[0] = -1.0
		self.triggerEff[0] = -1.0

		self.jet1JEC[0] = -1.0
		self.jet2JEC[0] = -1.0

		self.invarmass[0] = -1.0
		self.jetangle[0] = -10.0

	def __del__(self):	
		# print self.out_info
		self.f.cd()
		self.f.Write()
		self.f.Close()
		if self.doModMass:
			self.modMassFile.Close()
		if self.doTrigger:
			self.triggerFile.Close()