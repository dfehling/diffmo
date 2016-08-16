import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
import math
from math import *
import sys
# from FWCore.Common import TriggerNames

from Analysis.Tools.JetTools import *

class tree_maker:
    def __init__(self, outputname, triggerFileStr, isMC, unfoldWeight, invMassCut, doUnfold):
        # load all the event info:
        # self.out_info = 0
        self.name = outputname
        self.triggerFileStr = triggerFileStr
        self.isMC = isMC
        self.unfoldWeight = unfoldWeight
        self.invMassCut = invMassCut
        self.doUnfold = doUnfold

        #General Quantities
        #N Primary Vertices
        self.npvHandle = Handle( "unsigned int" )
        # We need the true distribution to properly PU reweight
        self.npvLabel  = ( "jhuGen", "npv" )
        # However this is always 0 for some reason
        if self.isMC:
            self.npvLabel  = ( "jhuGen", "npvTrue" )

        #MET Pt and Phi - why not eta?
        self.metPtHandle = Handle ( "double" )
        self.metPtLabel  = ( "jhuGen", "metpt" )
        self.metPhiHandle = Handle ( "double" )
        self.metPhiLabel  = ( "jhuGen", "metphi" )
        
        #Pruned Jet Collection
        #We need this for CSV values and subjet CSV values
        self.prunedHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
        self.prunedLabel  = ( "jhuCa8pp", "PrunedCA8CORR" )

        #Pruned Uncorrected Jet Collection
        #Do we need to compare the uncorrected jets to truth?
        self.prunedUCHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
        self.prunedUCLabel  = ( "jhuCa8pp", "PrunedCA8" )

        #Btagging CSV values
        self.CSVHandle = Handle( "std::vector<double>" )
        self.CSVLabel  = ( "jhuCa8pp", "PrunedCA8csv" )

        #Subjet btagging CSV values
        self.subjet1CSVHandle = Handle( "std::vector<double>" )
        self.subjet1CSVLabel  = ( "jhuCa8pp", "PrunedCA8sub0csv" )
        self.subjet2CSVHandle = Handle( "std::vector<double>" )
        self.subjet2CSVLabel  = ( "jhuCa8pp", "PrunedCA8sub1csv" )
        self.subjet3CSVHandle = Handle( "std::vector<double>" )
        self.subjet3CSVLabel  = ( "jhuCa8pp", "PrunedCA8sub2csv" )
        self.subjet4CSVHandle = Handle( "std::vector<double>" )
        self.subjet4CSVLabel  = ( "jhuCa8pp", "PrunedCA8sub3csv" )

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

        #GenParticles for MCunfolding
        self.genParticlesHandle = Handle( "vector<reco::GenParticle>" )
        self.genParticlesLabel = ( "prunedGenParticles" )

        self.__book__()

    def __book__(self):
    
        # if (self.triggerFileStr == ''):
        #     self.doTrigger = False
        # else:
        #     self.triggerFile = ROOT.TFile(self.triggerFileStr + ".root")
        #     self.triggerFile.cd()
        #     self.trigger = self.triggerFile.Get("TRIGGER_EFF").Clone()
        #     self.trigger.SetName('trigger')
        #     ROOT.SetOwnership( self.trigger, False )
        #     self.doTrigger = True

        print "Booking Histograms and Trees..."
        self.f = ROOT.TFile( self.name + ".root", "recreate" )
        self.f.cd()
        self.treeVars = ROOT.TTree('treeVars', 'treeVars')
        if self.isMC:
            self.treeVars.SetWeight(self.unfoldWeight)

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

        self.jet1topTagged = array('i', [-1])
        self.jet2topTagged = array('i', [-1])
        self.jet1bTagged = array('i', [-1])
        self.jet2bTagged = array('i', [-1])

        self.htSum = array('f', [-1.0])
        self.triggerEff = array('f', [1.0])

        self.pass400pt = array('f', [-1.0])
        self.pass450pt = array('f', [-1.0])
        self.fail400pt = array('f', [-1.0])
        self.fail450pt = array('f', [-1.0])

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

        self.treeVars.Branch('deltaY', self.deltaY, 'deltaY/F')
        self.treeVars.Branch('deltaPhi', self.deltaPhi, 'deltaPhi/F')

        self.treeVars.Branch('jet1topTagged', self.jet1topTagged, 'jet1topTagged/I')
        self.treeVars.Branch('jet2topTagged', self.jet2topTagged, 'jet2topTagged/I')
        self.treeVars.Branch('jet1bTagged', self.jet1bTagged, 'jet1bTagged/I')
        self.treeVars.Branch('jet2bTagged', self.jet2bTagged, 'jet2bTagged/I')

        self.treeVars.Branch('htSum', self.htSum, 'htSum/F')
        self.treeVars.Branch('triggerEff', self.triggerEff, 'triggerEff/F')

        self.Mtt = array('f', [-1.0])
        self.jetangle = array('f', [-10.0])
        self.treeVars.Branch('Mtt', self.Mtt, 'Mtt/F')
        self.treeVars.Branch('angle_between_jets', self.jetangle, 'jetangle/F')

        self.genPartonJet1pt = array('f', [-1.0])
        self.genPartonJet2pt = array('f', [-1.0])
        self.genPartonJet1eta = array('f', [-10.0])
        self.genPartonJet2eta = array('f', [-10.0])
        self.genPartonJet1phi = array('f', [-10.0])
        self.genPartonJet2phi = array('f', [-10.0])
        self.genPartonJet1mass = array('f', [-1.0])
        self.genPartonJet2mass = array('f', [-1.0])
        self.genPartonJetMtt = array('f', [-1.0])
        self.genPartonJet1id = array('i', [0])

        self.treeVars.Branch('genPartonJet1pt', self.genPartonJet1pt, 'genPartonJet1pt/F')
        self.treeVars.Branch('genPartonJet2pt', self.genPartonJet2pt, 'genPartonJet2pt/F')
        self.treeVars.Branch('genPartonJet1eta', self.genPartonJet1eta, 'genPartonJet1eta/F')
        self.treeVars.Branch('genPartonJet2eta', self.genPartonJet2eta, 'genPartonJet2eta/F')
        self.treeVars.Branch('genPartonJet1phi', self.genPartonJet1phi, 'genPartonJet1phi/F')
        self.treeVars.Branch('genPartonJet2phi', self.genPartonJet2phi, 'genPartonJet2phi/F')
        self.treeVars.Branch('genPartonJet1mass', self.genPartonJet1mass, 'genPartonJet1mass/F')
        self.treeVars.Branch('genPartonJet2mass', self.genPartonJet2mass, 'genPartonJet2mass/F')
        self.treeVars.Branch('genPartonJetMtt',   self.genPartonJetMtt, 'genPartonJetMtt/F')
        self.treeVars.Branch('genPartonJet1id',   self.genPartonJet1id, 'genPartonJet1id/I')

        # self.isGenLeptonic = array('i', [-1])
        self.isGenHadronic = array('i', [-1])
        self.passKinCuts = array('i', [-1])
        self.passFullSel = array('i', [-1])
        # self.treeVars.Branch('isGenLeptonic', self.isGenLeptonic, 'isGenLeptonic/I')
        self.treeVars.Branch('isGenHadronic', self.isGenHadronic, 'isGenHadronic/I')
        self.treeVars.Branch('passKinCuts', self.passKinCuts, 'passKinCuts/I')
        self.treeVars.Branch('passFullSel', self.passFullSel, 'passFullSel/I')

        self.genPartonJet1matchPt = array('f', [-1.0])
        self.genPartonJet2matchPt = array('f', [-1.0])
        self.genPartonJet1match = array('i', [-10])
        self.genPartonJet2match = array('i', [-10])
        self.treeVars.Branch('genPartonJet1matchPt', self.genPartonJet1matchPt, 'genPartonJet1matchPt/F')
        self.treeVars.Branch('genPartonJet2matchPt', self.genPartonJet2matchPt, 'genPartonJet2matchPt/F')
        self.treeVars.Branch('genPartonJet1match', self.genPartonJet1match, 'genPartonJet1match/I')
        self.treeVars.Branch('genPartonJet2match', self.genPartonJet2match, 'genPartonJet2match/I')
        

        self.ucJet1pt = array('f', [-1.0])
        self.ucJet2pt = array('f', [-1.0])
        self.ucJet1eta = array('f', [-10.0])
        self.ucJet2eta = array('f', [-10.0])
        self.ucJet1phi = array('f', [-10.0])
        self.ucJet2phi = array('f', [-10.0])
        self.ucJet1mass = array('f', [-1.0])
        self.ucJet2mass = array('f', [-1.0])

        self.treeVars.Branch('ucJet1pt', self.ucJet1pt, 'ucJet1pt/F')
        self.treeVars.Branch('ucJet2pt', self.ucJet2pt, 'ucJet2pt/F')
        self.treeVars.Branch('ucJet1eta', self.ucJet1eta, 'ucJet1eta/F')
        self.treeVars.Branch('ucJet2eta', self.ucJet2eta, 'ucJet2eta/F')
        self.treeVars.Branch('ucJet1phi', self.ucJet1phi, 'ucJet1phi/F')
        self.treeVars.Branch('ucJet2phi', self.ucJet2phi, 'ucJet2phi/F')
        self.treeVars.Branch('ucJet1mass', self.ucJet1mass, 'ucJet1mass/F')
        self.treeVars.Branch('ucJet2mass', self.ucJet2mass, 'ucJet2mass/F')

        self.treeVars.Branch('pass400pt', self.pass400pt, 'pass400pt/F')
        self.treeVars.Branch('pass450pt', self.pass450pt, 'pass450pt/F')
        self.treeVars.Branch('fail400pt', self.fail400pt, 'fail400pt/F')
        self.treeVars.Branch('fail450pt', self.fail450pt, 'fail450pt/F')


        if self.isMC == True and self.doUnfold == True and self.unfoldWeight != 1.0:
            ROOT.gSystem.Load("RooUnfold-1.1.1/libRooUnfold")
            # dummy histogram used only to specify dimensions for reponse matrix
            # ptbins = array('d',[0.0,200.0,400.0,500.0,600.0,700.0,800.0,1200.0,2000.0])
            ptbins = array('d',[400.0,500.0,600.0,700.0,800.0,1600.0])
            h_bins = ROOT.TH1F("bins", ";;", len(ptbins)-1, ptbins)
            # self.doUnfold = True
            self.response = ROOT.RooUnfoldResponse(h_bins, h_bins)
            self.response.SetName('response_pt')
        else:
            self.doUnfold = False

        self.cutflow = ROOT. TH1D("cutflow", "cutfow", 10, 0, 10 )
        self.cutflow.Sumw2()

        # self.allEvents = ROOT.TTree('allEvents', 'allEvents')
        # self.allJet1pt = array('f', [-1.0])
        # self.allEvents.Branch('allJet1pt', self.allJet1pt, 'allJet1pt/F')
        # self.allEvents.Branch('allJet2pt', self.allJet2pt, 'allJet2pt/F')
        # self.allEvents.Branch('allJet1eta', self.allJet1eta, 'allJet1eta/F')
        # self.allEvents.Branch('allJet2eta', self.allJet2eta, 'allJet2eta/F')
        # self.allEvents.Branch('allGenJet1pt', self.allGenJet1pt, 'allGenJet1pt/F')
        # self.allEvents.Branch('allGenJet2pt', self.allGenJet2pt, 'allGenJet2pt/F')
        # self.allEvents.Branch('allGenJet1eta', self.allGenJet1eta, 'allGenJet1eta/F')
        # self.allEvents.Branch('allGenJet2eta', self.allGenJet2eta, 'allGenJet2eta/F')

        # We don't want to fill anything if the event is outside the mtt range we're looking for
        self.noFill = 0

    def analyze(self, event):
        # // To get the trigger names
        # edm::Handle< edm::TriggerResults > h_trigresults
        # edm::InputTag triggerResultsSrc_("TriggerResults", "", "HLT")
        # iEvent.getByLabel( triggerResultsSrc_, h_trigresults )
        # const edm::TriggerResults* thing = h_trigresults.product()
        # edm::TriggerNames const & trig_names = iEvent.triggerNames(*thing)
        # std::vector<std::string> const & trig_strings = trig_names.triggerNames()

        # if self.doTrigger == True:
        if self.isMC == False:
            event.getByLabel (self.triggerLabel, self.triggerHandle)
            trigNames = event.object().triggerNames(self.triggerHandle.product())

            path400 = "HLT_HT400"
            path450 = "HLT_HT450"

            index400 = trigNames.triggerIndex(path400)
            index450 = trigNames.triggerIndex(path450)

            for version in ["_v1","_v2","_v3","_v4","_v5","_v6","_v7"]:
                newpath400 = path400+version
                newindex400 = trigNames.triggerIndex(newpath400)
                if newindex400==trigNames.size():
                    continue
                else:
                    newpath450 = path450+version
                    break

            index400 = trigNames.triggerIndex(newpath400)
            index450 = trigNames.triggerIndex(newpath450)

            # print index400,index450
            # print trigNames.triggerName(index400), trigNames.triggerName(index450)

            pass400 = self.triggerHandle.product().accept(index400)
            pass450 = self.triggerHandle.product().accept(index450)

        self.run[0] = event.object().id().run()
        self.event[0] = event.object().id().event()
        self.lumi[0] = event.object().id().luminosityBlock()

        event.getByLabel (self.npvLabel, self.npvHandle)
        event.getByLabel (self.metPtLabel, self.metPtHandle)
        event.getByLabel (self.metPhiLabel, self.metPhiHandle)

        event.getByLabel (self.prunedLabel, self.prunedHandle)
        event.getByLabel (self.prunedUCLabel, self.prunedUCHandle)
        event.getByLabel (self.unprunedLabel, self.unprunedHandle)

        event.getByLabel (self.CSVLabel, self.CSVHandle)

        event.getByLabel (self.subjet1CSVLabel, self.subjet1CSVHandle)
        event.getByLabel (self.subjet2CSVLabel, self.subjet2CSVHandle)
        event.getByLabel (self.subjet3CSVLabel, self.subjet3CSVHandle)
        event.getByLabel (self.subjet4CSVLabel, self.subjet4CSVHandle)

        event.getByLabel (self.t1Label, self.t1Handle)
        event.getByLabel (self.t2Label, self.t2Handle)
        event.getByLabel (self.t3Label, self.t3Handle)
        event.getByLabel (self.t4Label, self.t4Handle)

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

        HTsum = 0

        #Total events seen
        self.index[0]=0;
        self.cutflow.Fill(1);

        unpj = self.unprunedHandle.product()
        pj = self.prunedHandle.product()
        ucPJ = self.prunedUCHandle.product()


        #Reorder to make sure the highest pT is first. This is after any JEC. Is this correct?
        #pt_sorted_jets = pj
        pt_sorted_jets = ReorderByPt(pj)

        #The very first thing we want to do is, if we are running over MC, to get the genParticle info
        #and save this. If we are running over the Powheg full sample, we need to only save events
        #with the generated mTT<700

        #This will only possibly be set true if isMC is true so we don't need both after the next bit
        # passParton = False
        # doUnfold = False

        # if self.isMC == True and self.unfoldWeight > 0:
        if self.doUnfold == True:
            # doUnfold = True
            event.getByLabel( self.genParticlesLabel, self.genParticlesHandle )
            genParticles  = self.genParticlesHandle.product()
            
            # if self.genParticlesHandle.isValid() == False :
            #     print "DEBUG: genParticlesPtHandle is not valid! continuing..." 
            #     return
            
            
            genT = ROOT.TLorentzVector()
            genTbar = ROOT.TLorentzVector()
            genJet1 = ROOT.TLorentzVector()
            genJet2 = ROOT.TLorentzVector()
            
            # hadTDecay = 0     # 1 = hadronic, 0 = leptonic
            # hadTbarDecay = 0    # 1 = hadronic, 0 = leptonic
            # self.isGenLeptonic[0] = 0
            # self.noFill = 0
            
            # loop over gen particles
            # we want to loop over all the gen particles and save some info
            # we want to find out if we have both a top and an antitop
            # if the decay chain includes a lepton, the decay is semileptonic - we don't want this
            #TODO add in cut for bottom quark to check btagging
            # Right now, we aren't worried about fakes
            self.isGenHadronic[0] = 1
            for igen in xrange( len(genParticles) ) :

                if  genParticles[igen].status() != 3 :
                    continue
                if  abs(genParticles[igen].pdgId()) == 11 or abs(genParticles[igen].pdgId()) == 13 or abs(genParticles[igen].pdgId()) == 15 :
                    self.isGenHadronic[0] = 0
                    continue
                    # self.noFill = 1
                    # return

                if  abs(genParticles[igen].pdgId()) != 6 :
                    continue
                
                # if  abs(genParticles[igen].pdgId()) > 16 :
                #     continue

                if genParticles[igen].pdgId() == 6 :
                    gen = ROOT.TLorentzVector()
                    gen.SetPtEtaPhiM( genParticles[igen].pt(), genParticles[igen].eta(), genParticles[igen].phi(), genParticles[igen].mass() )
                    genT = gen
                    # self.isGenHadronic
                    # hadTDecay = 1
                elif genParticles[igen].pdgId() == -6 :
                    gen = ROOT.TLorentzVector()
                    gen.SetPtEtaPhiM( genParticles[igen].pt(), genParticles[igen].eta(), genParticles[igen].phi(), genParticles[igen].mass() )
                    genTbar = gen

            #We only want the hadronic events
            # if isGenLeptonic:
            #     return


                    # hadTbarDecay = 1
                # If there is an antilepton (e+, mu+, tau+) then the T is leptonic
                # elif ( genParticles[igen].pdgId() == -11 or genParticles[igen].pdgId() == -13 or genParticles[igen].pdgId() == -15) :
                #     hadTDecay = 0
                # # If there is an lepton (e-, mu-, tau-) then the Tbar is leptonic
                # elif ( genParticles[igen].pdgId() == 11 or genParticles[igen].pdgId() == 13 or genParticles[igen].pdgId() == 15) :                
                #     hadTbarDecay = 0

            #This might crash, but I want to try to save mtt for all possible events, even if there's a lepton in the event

            # If we have both a t and a tbar, then we've passed the parton level selection
            # If these events have a mTT > 700 for the full sample, we return...fill NOTHING
            # if hadTDecay == True and hadTbarDecay == True:
                #Events generated with ttbar
                # self.cutflow.Fill(11)
                # passParton = True
                # genPartonJetMtt = (genT+genTbar).M()
                #Allows us to get the mTT<700

            #Sort the t and tbar by pT
            # if genT.Pt() > genTbar.Pt():
            genJet1 = genT
            genJet2 = genTbar
            #     self.genPartonJet1id[0] = 6
            # else:
            #     genJet1 = genTbar
            #     genJet2 = genT
            #     self.genPartonJet1id[0] = -6

            if genJet2.Pt() > 0 and genJet1.Pt() > 0:
                genPartonJetMtt = (genT+genTbar).M()
                self.genPartonJetMtt[0] = genPartonJetMtt
            
            if self.invMassCut > 0 and genPartonJetMtt >= self.invMassCut:
                self.noFill = 1
                return
            else:
                self.noFill = 0

            #Events passing massCut (full powheg only difference)
            # if self.isGenLeptonic[0] == 0:
            #     self.isGenHadronic[0] = 1
            # elif self.isGenLeptonic[0] == 1:
            #     self.isGenHadronic[0]
            self.cutflow.Fill(12)
            self.index[0] = 2
            #Fill our generator level cuts. we need to remember to fill the tree before returning in the Miss case
            
            if genJet1.Pt() > 0:
                self.genPartonJet1pt[0] = genJet1.Pt()
                self.genPartonJet1eta[0] = genJet1.Eta()
                self.genPartonJet1phi[0] = genJet1.Phi()
                self.genPartonJet1mass[0] = genJet1.M()

            if genJet2.Pt() > 0:
                self.genPartonJet2pt[0] = genJet2.Pt()
                self.genPartonJet2eta[0] = genJet2.Eta()
                self.genPartonJet2phi[0] = genJet2.Phi()
                self.genPartonJet2mass[0] = genJet2.M()
            # self.genPartonJetMtt[0] = genPartonJetMtt
        # else:
        #     self.isGenHadronic[0] = 0

        #Fill the pT of the leading jet before ANY reco cuts so we can get the efficiency to get the diff xsection
        # self.allEvents.Fill(pt_sorted_jets.Pt());
        # if len(pj) > 0:
        #     self.allJet1pt[0] = pt_sorted_jets[0].Pt()
        #     self.allEvents.Fill()
        #     #Events with valid handle
        #     self.cutflow.Fill(2);

        self.passKinCuts[0] = 0

        #Two jets:
        if len(pj) < 2 or len(unpj) < 2 or len(ucPJ) < 2:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return
    
        #Events with 2 or more valid jets
        self.cutflow.Fill(3);

        #Make sure there's at least a top condidiate by checking the pT
        nTopCand = 0
        for i in range(0,len(pj) ) :
            if pj[i].pt() > 400. :
                nTopCand = nTopCand + 1
                HTsum += pj[i].pt()
        if nTopCand < 2:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return
        #make sure we are within correct eta
        elif abs(pt_sorted_jets[0].Eta()) > 2.4 or abs(pt_sorted_jets[1].Eta()) > 2.4:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return

        #Events with 2 top candidates (pt > 400) and eta < 2.4
        self.cutflow.Fill(4)

        self.htSum[0] = HTsum
        self.MET[0] = metPt
        self.npv[0] = npv

        #We need these for the substructure
        ca1 = ROOT.TLorentzVector()
        ca2 = ROOT.TLorentzVector()

        #We create the CA jet to match and fill the jet parameters          
        ca1.SetPtEtaPhiM(pt_sorted_jets[0].Pt(), pt_sorted_jets[0].Eta(), pt_sorted_jets[0].Phi(), pt_sorted_jets[0].M())
        ca2.SetPtEtaPhiM(pt_sorted_jets[1].Pt(), pt_sorted_jets[1].Eta(), pt_sorted_jets[1].Phi(), pt_sorted_jets[1].M())

        #Fill the trigger info
        if self.isMC == False:
            if pass400:
                self.pass400pt[0] = ca1.Pt()
            else:
                self.fail400pt[0] = ca1.Pt()
            if pass450:
                self.pass450pt[0] = ca1.Pt()
            else:
                self.fail450pt[0] = ca1.Pt()

        #Match unpruned jets with pruned - so we have both subjet btagging and nsubjettiness
        #This returns the index of the jet in the first collection that matches within dr = 0.4 to the jet of the second argument
        jet1matchIndex = MatchCol(unpj, ca1)
        jet2matchIndex = MatchCol(unpj, ca2)
        #This is temporary so we get the correct pj index for CSV and subjet CSV
        jet1matchIndex_pj = MatchCol(pj, ca1)
        jet2matchIndex_pj = MatchCol(pj, ca2)
        #This is temporary so we get the uncorrected pj index
        jet1matchIndex_ucPJ = MatchCol(ucPJ, ca1)
        jet2matchIndex_ucPJ = MatchCol(ucPJ, ca2)

        #Let's see what matches with the truth info
        if self.doUnfold == True:
            maxDR = 0.4
            j = -1
            if genJet1.Pt() > 0:
                for i in range(len(pt_sorted_jets)):
                    C = ROOT.TLorentzVector()
                    C.SetPtEtaPhiM( pt_sorted_jets[i].Pt(), pt_sorted_jets[i].Eta(), pt_sorted_jets[i].Phi(), pt_sorted_jets[i].M() )
                    dr = abs(genJet1.DeltaR(C))
                    if dr < maxDR :
                        j = i
                        break
            # if dr > maxDR:
            self.genPartonJet1match[0] = j
            if self.genPartonJet1match[0]>-1:
                self.genPartonJet1matchPt[0] = genParticles[self.genPartonJet1match[0]].pt()

            maxDR = 0.4
            j = -1
            if genJet2.Pt() > 0:
                for i in range(len(pt_sorted_jets)):
                    C = ROOT.TLorentzVector()
                    C.SetPtEtaPhiM( pt_sorted_jets[i].Pt(), pt_sorted_jets[i].Eta(), pt_sorted_jets[i].Phi(), pt_sorted_jets[i].M() )
                    dr = abs(genJet2.DeltaR(C))
                    if dr < maxDR :
                        j = i
                        break
            self.genPartonJet2match[0] = j
            if self.genPartonJet2match[0]>-1:
                self.genPartonJet2matchPt[0] = genParticles[self.genPartonJet2match[0]].pt()
            
        #At this point, we have 2 jets which have passed the basic kinematic cuts of pt>400 and eta<2.4
        self.passKinCuts[0] = 1
        self.passFullSel[0] = 0
        self.jet1pt[0] = ca1.Pt()
        self.jet2pt[0] = ca2.Pt()
        self.jet1eta[0] = ca1.Eta()
        self.jet2eta[0] = ca2.Eta()
        self.jet1phi[0] = ca1.Phi()
        self.jet2phi[0] = ca2.Phi()
        self.jet1mass[0] = ca1.M()
        self.jet2mass[0] = ca2.M()

        #Make sure matching is correct for jet1
        if jet1matchIndex==-1 or jet1matchIndex_pj==-1 or jet1matchIndex_ucPJ==-1:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return
        
        #Make sure matching is correct for jet2
        if jet2matchIndex==-1 or jet2matchIndex_pj==-1 or jet2matchIndex_ucPJ==-1:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return
        
        #Events with valid ca-other jet matches
        self.cutflow.Fill(5)

        # #Events with filled subjets
        # self.cutflow.Fill(6)

        #Nsubjettiness
        if Tau2[jet1matchIndex] == 0 or Tau2[jet2matchIndex] == 0:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
            return
        if Tau1[jet1matchIndex] == 0 or Tau1[jet2matchIndex] == 0:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
                # self.treeVars.Fill()
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
        
        #Let's fill the uncorrected jet info
        self.ucJet1pt[0] = ucPJ[jet1matchIndex_ucPJ].Pt()
        self.ucJet2pt[0] = ucPJ[jet2matchIndex_ucPJ].Pt()
        self.ucJet1eta[0] = ucPJ[jet1matchIndex_ucPJ].Eta()
        self.ucJet2eta[0] = ucPJ[jet2matchIndex_ucPJ].Eta()
        self.ucJet1phi[0] = ucPJ[jet1matchIndex_ucPJ].Phi()
        self.ucJet2phi[0] = ucPJ[jet2matchIndex_ucPJ].Phi()
        self.ucJet1mass[0] = ucPJ[jet1matchIndex_ucPJ].M()
        self.ucJet2mass[0] = ucPJ[jet2matchIndex_ucPJ].M()
        
        #Invariant Mass         
        self.Mtt[0] = (ca1+ca2).M()

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
        jet1subjetCSVs.append(subjet1CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet2CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet3CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet4CSV[jet1matchIndex_pj])

        jet2subjetCSVs = []
        jet2subjetCSVs.append(subjet1CSV[jet2matchIndex_pj])
        jet2subjetCSVs.append(subjet2CSV[jet2matchIndex_pj])
        jet2subjetCSVs.append(subjet3CSV[jet2matchIndex_pj])
        jet2subjetCSVs.append(subjet4CSV[jet2matchIndex_pj])

        self.jet1maxSubjetCSV[0] = max(jet1subjetCSVs)
        self.jet2maxSubjetCSV[0] = max(jet2subjetCSVs)

        self.jet1bTagged[0] = self.jet1maxSubjetCSV[0] > 0.679
        self.jet2bTagged[0] = self.jet2maxSubjetCSV[0] > 0.679

        # self.jet1topTagged[0] = self.jet1mass[0] > 140.0 and self.jet1mass[0] < 250.0 and self.jet1pt[0] > 400.
        # self.jet2topTagged[0] = self.jet2mass[0] > 140.0 and self.jet2mass[0] < 250.0 and self.jet2pt[0] > 400.

        #This is what we are actually using
        self.jet1topTagged[0] = self.jet1mass[0] > 100.0 and self.jet1pt[0] > 400. and self.jet1tau21[0] > 0.1 and self.jet1tau32[0] > 0. and self.jet1tau32[0] < 0.55
        self.jet2topTagged[0] = self.jet2mass[0] > 140.0 and self.jet2mass[0] < 250.0 and self.jet2pt[0] > 400.

        if self.jet1topTagged[0] and self.jet2topTagged[0]: 
            self.index[0] = 1
            self.passFullSel[0] = 1
            self.cutflow.Fill(7)
            if self.doUnfold:
                if self.isGenHadronic[0] == 1:
                    self.response.Fill(self.jet1pt[0], self.genPartonJet1pt[0], self.unfoldWeight)
                else:
                    self.response.Fake(self.genPartonJet1pt[0], self.unfoldWeight)
        elif self.doUnfold and self.isGenHadronic[0] == 1:
            self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)

        # self.treeVars.Fill()
        return

    def reset(self):
        # self.allEvents.Fill()
        if not self.noFill:
            self.treeVars.Fill()

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
        self.deltaY[0] = -10.0
        self.deltaPhi[0] = -10.0

        self.jet1bTagged[0] = -1
        self.jet2bTagged[0] = -1
        self.jet1topTagged[0] = -1
        self.jet2topTagged[0] = -1

        self.htSum[0] = -1.0
        self.triggerEff[0] = -1.0

        self.Mtt[0] = -1.0
        self.jetangle[0] = -10.0

        self.genPartonJet1pt[0] = -1.0
        self.genPartonJet2pt[0] = -1.0
        self.genPartonJet1eta[0] = -10.0
        self.genPartonJet2eta[0] = -10.0
        self.genPartonJet1phi[0] = -10.0
        self.genPartonJet2phi[0] = -10.0
        self.genPartonJet1mass[0] = -1.0
        self.genPartonJet2mass[0] = -1.0
        self.genPartonJetMtt[0] = -1.0
        self.genPartonJet1id[0] = 0

        self.genPartonJet1matchPt[0] = -1.0
        self.genPartonJet2matchPt[0] = -1.0
        self.genPartonJet1match[0] = -10
        self.genPartonJet2match[0] = -10
        
        self.ucJet1pt[0] = -1.0
        self.ucJet2pt[0] = -1.0
        self.ucJet1eta[0] = -10.0
        self.ucJet2eta[0] = -10.0
        self.ucJet1phi[0] = -10.0
        self.ucJet2phi[0] = -10.0
        self.ucJet1mass[0] = -1.0
        self.ucJet2mass[0] = -1.0

        self.isGenHadronic[0] = -1
        self.passKinCuts[0] = -1
        self.passFullSel[0] = -1

        self.pass400pt[0] = -1
        self.pass450pt[0] = -1
        self.fail400pt[0] = -1
        self.fail450pt[0] = -1

        # self.allJet1pt[0] = -1.0

    def __del__(self):  
        # print self.out_info
        self.f.cd()
        if self.doUnfold == True:
            self.response.Write()
        self.f.Write()
        self.f.Close()
        # if self.doTrigger:
        #     self.triggerFile.Close()