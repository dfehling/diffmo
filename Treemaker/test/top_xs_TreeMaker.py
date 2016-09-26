import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
import math
from math import *
import sys
# from FWCore.Common import TriggerNames

from Analysis.Tools.JetTools import *

class tree_maker:
    def __init__(self, outputname, triggerFileStr, useTrigger, pileupFileStr, usePileup, usePDF, isMC, unfoldWeight, invMassCut, doUnfold):
        # load all the event info:
        # self.out_info = 0
        self.name = outputname
        self.triggerFileStr = triggerFileStr
        self.useTrigger = useTrigger
        self.pileupFileStr = pileupFileStr
        self.usePileup = usePileup
        self.usePDF = usePDF
        self.isMC = isMC
        self.unfoldWeight = unfoldWeight
        self.invMassCut = invMassCut
        self.doUnfold = doUnfold

        self.btagSF = 1.08513350715
        self.nsubSF = 0.814651566377

        if not self.isMC:
            self.useTrigger = False
            self.usePileup = False
            self.usePDF = False

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

        #PDF
        self.pdfWeightHandle = Handle( "std::vector<double>" )
        self.pdfWeightLabel = ( "pdfWeights", "CT10" )

        self.__book__()

    def __book__(self):
    
        if self.isMC:
            if (self.triggerFileStr != ''):
                self.triggerFile = ROOT.TFile(self.triggerFileStr + ".root")
                self.triggerFile.cd()
                self.trigger = self.triggerFile.Get("trigger").Clone()
                self.trigger.SetName('trigger')
                ROOT.SetOwnership( self.trigger, False )
                self.useTrigger = True
                
            if (self.pileupFileStr != ''):
                self.pileupFile = ROOT.TFile(self.pileupFileStr + ".root")
                self.pileupFile.cd()
                self.pileup = self.pileupFile.Get("pileup").Clone()
                self.pileup.SetName('pileup')
                ROOT.SetOwnership( self.pileup, False )
                self.usePileup = True

        print "Booking Histograms and Trees..."
        self.f = ROOT.TFile( self.name + ".root", "recreate" )
        self.f.cd()
        self.treeVars = ROOT.TTree('treeVars', 'treeVars')
        if self.isMC and self.doUnfold:
            self.treeVars.SetWeight(self.unfoldWeight)

        self.run = array('i', [-1])
        self.event = array('l', [-1])
        self.lumi = array('i', [-1])
        
        self.index = array('i', [-1])

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
        self.nJets = array('i',[-1])
                
        self.jet1topTagged = array('i', [-1])
        self.jet2topTagged = array('i', [-1])
        self.jet1bTagged = array('i', [-1])
        self.jet2bTagged = array('i', [-1])

        self.htSum = array('f', [-1.0])
        self.triggerWeight = array('f', [1.0])
        self.pileupWeight = array('f', [1.0])
        self.pdfWeight = array('f', [1.0])
        self.pdfWeightNom = array('f', [1.0])
        self.pdfWeightUp = array('f', [1.0])
        self.pdfWeightDown = array('f', [1.0])
        self.unfoldWeightUsed = array('f', [1.0])

        self.pass400pt = array('f', [-1.0])
        self.pass750pt = array('f', [-1.0])
        self.fail400pt = array('f', [-1.0])
        self.fail750pt = array('f', [-1.0])

        self.treeVars.Branch('run', self.run, 'run/I')
        self.treeVars.Branch('event', self.event, 'event/L')
        self.treeVars.Branch('lumi', self.lumi, 'lumi/I')

        self.treeVars.Branch('npv', self.npv, 'npv/I')
        self.treeVars.Branch('index', self.index, 'index/I')
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
        self.treeVars.Branch('nJets', self.nJets, 'nJets/I')
        self.treeVars.Branch('deltaY', self.deltaY, 'deltaY/F')
        self.treeVars.Branch('deltaPhi', self.deltaPhi, 'deltaPhi/F')

        self.treeVars.Branch('jet1topTagged', self.jet1topTagged, 'jet1topTagged/I')
        self.treeVars.Branch('jet2topTagged', self.jet2topTagged, 'jet2topTagged/I')
        self.treeVars.Branch('jet1bTagged', self.jet1bTagged, 'jet1bTagged/I')
        self.treeVars.Branch('jet2bTagged', self.jet2bTagged, 'jet2bTagged/I')

        self.treeVars.Branch('htSum', self.htSum, 'htSum/F')
        self.treeVars.Branch('triggerWeight', self.triggerWeight, 'triggerWeight/F')
        self.treeVars.Branch('pileupWeight', self.pileupWeight, 'pileupWeight/F')
        self.treeVars.Branch('pdfWeight', self.pdfWeight, 'pdfWeight/F')
        self.treeVars.Branch('pdfWeightNom', self.pdfWeightNom, 'pdfWeightNom/F')
        self.treeVars.Branch('pdfWeightUp', self.pdfWeightUp, 'pdfWeightUp/F')
        self.treeVars.Branch('pdfWeightDown', self.pdfWeightDown, 'pdfWeightDown/F')
        self.treeVars.Branch('unfoldWeightUsed', self.unfoldWeightUsed, 'unfoldWeightUsed/F')

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
        self.treeVars.Branch('pass750pt', self.pass750pt, 'pass750pt/F')
        self.treeVars.Branch('fail400pt', self.fail400pt, 'fail400pt/F')
        self.treeVars.Branch('fail750pt', self.fail750pt, 'fail750pt/F')


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

        # We don't want to fill anything if the event is outside the mtt range we're looking for
        self.noFill = 0

    def analyze(self, event):

        self.run[0] = event.object().id().run()
        self.event[0] = event.object().id().event()
        self.lumi[0] = event.object().id().luminosityBlock()

        event.getByLabel (self.npvLabel, self.npvHandle)
        event.getByLabel (self.metPtLabel, self.metPtHandle)
        event.getByLabel (self.metPhiLabel, self.metPhiHandle)
        npv = self.npvHandle.product()[0]
        metPt = self.metPtHandle.product()[0]
        metPhi = self.metPhiHandle.product()[0]

       
        #Save information about trigger paths so we can calculate the trigger sf later
        event.getByLabel (self.triggerLabel, self.triggerHandle)
        trigNames = event.object().triggerNames(self.triggerHandle.product())

        path400 = "HLT_HT400"
        path750 = "HLT_HT750"
        index400 = trigNames.triggerIndex(path400)
        index750 = trigNames.triggerIndex(path750)

        for version in ["_v1","_v2","_v3","_v4","_v5","_v6","_v7"]:
            newpath400 = path400+version
            newindex400 = trigNames.triggerIndex(newpath400)
            if newindex400==trigNames.size():
                continue
            else:
                newpath750 = path750+version
                break

        index400 = trigNames.triggerIndex(newpath400)
        index750 = trigNames.triggerIndex(newpath750)
        pass400 = self.triggerHandle.product().accept(index400)
        pass750 = self.triggerHandle.product().accept(index750)


        #Pileup SF. If data, do nothing. Otherwise you should be passed the appropriate histogram for nominal, up, or down.
        #Simply get the weight from it
        if self.isMC == False:
            self.pileupWeight[0] = 1.0
        elif self.usePileup:
            self.pileupWeight[0] = self.pileup.GetBinContent(self.pileup.FindBin(npv))
        else:
            self.pileupWeight[0] = 1.0

        #Include pileup reweighting for the response matrix
        weight = 0
        weight = self.unfoldWeight * self.pileupWeight[0]

        #PDF
        #If present and in MC, calculate PDF to scale up and down. Otherwise return weight = 1
        if self.usePDF == False:
            self.pdfWeight[0] = 1.0
        else:
            event.getByLabel (self.pdfWeightLabel, self.pdfWeightHandle)
            pdfWeights = self.pdfWeightHandle.product()
            self.pdfWeightNom[0] = pdfWeights[0]
            # print "Event weight for central PDF:",pdfWeights[0]

            tempPdfWeight1 = 0
            tempPdfWeight2 = 0
            tempPdfWeightUp = 0
            tempPdfWeightDown = 0
            for pdfIndex in xrange(1,len(pdfWeights),2):
                tempPdfWeight1 = pdfWeights[pdfIndex]  -1.0
                tempPdfWeight2 = pdfWeights[pdfIndex+1]-1.0

                # if tempPdfWeight1>1 or tempPdfWeight2>1:
                #     print tempPdfWeight1,tempPdfWeight2

                if(tempPdfWeight1>tempPdfWeight2):
                    if(tempPdfWeight1<0.):
                        tempPdfWeight1=0.
                    if(tempPdfWeight2>0.):
                        tempPdfWeight2=0.
                    tempPdfWeightUp += tempPdfWeight1*tempPdfWeight1
                    tempPdfWeightDown += tempPdfWeight2*tempPdfWeight2
                else:
                    if(tempPdfWeight2<0.):
                        tempPdfWeight2=0.
                    if(tempPdfWeight1>0.):
                        tempPdfWeight1=0.
                    tempPdfWeightUp += tempPdfWeight2*tempPdfWeight2
                    tempPdfWeightDown += tempPdfWeight1*tempPdfWeight1

            # if sqrt(tempPdfWeightDown) > .9:
            #     print sqrt(tempPdfWeightUp),sqrt(tempPdfWeightDown)

            self.pdfWeightUp[0] = 1.0 + sqrt(tempPdfWeightUp)
            self.pdfWeightDown[0] = 1.0 - sqrt(tempPdfWeightDown)
            if self.usePDF>1:
                self.pdfWeight[0] = self.pdfWeightUp[0]
            elif self.usePDF<0:
                self.pdfWeight[0] = self.pdfWeightDown[0]
            else:
                self.pdfWeight[0] = 1.0

            if self.pdfWeight[0]>10.0:
                return

        #Include pdf reweighting for the response matrix
        weight = weight * self.pdfWeight[0]

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
        self.nJets[0] = len(pt_sorted_jets)

        #The very first thing we want to do is, if we are running over MC, to get the genParticle info
        #and save this. If we are running over the Powheg full sample, we need to only save events
        #with the generated mTT<700

        #This will only possibly be set true if isMC is true so we don't need both after the next bit
        # passParton = False
        # doUnfold = False

        # if self.isMC == True and self.unfoldWeight > 0:
        if self.doUnfold == True:
            event.getByLabel( self.genParticlesLabel, self.genParticlesHandle )
            genParticles  = self.genParticlesHandle.product()
            
            genT = ROOT.TLorentzVector()
            genTbar = ROOT.TLorentzVector()
            genJet1 = ROOT.TLorentzVector()
            genJet2 = ROOT.TLorentzVector()
            
            # loop over gen particles
            # we want to loop over all the gen particles and save some info
            # we want to find out if we have both a top and an antitop
            # if the decay chain includes a lepton, the decay is semileptonic - we don't want this
            #TODO add in cut for bottom quark to check btagging
            # Right now, we aren't worried about fakes
            self.isGenHadronic[0] = 1
            for igen in xrange( len(genParticles) ) :

                #Make sure we have a stable particle
                if  genParticles[igen].status() != 3 :
                    continue
                #If we have a lepton, note that we're not hadronic
                if  abs(genParticles[igen].pdgId()) == 11 or abs(genParticles[igen].pdgId()) == 13 or abs(genParticles[igen].pdgId()) == 15 :
                    self.isGenHadronic[0] = 0
                    continue
                #Make sure we have a top quark/antiquark
                if  abs(genParticles[igen].pdgId()) != 6 :
                    continue

                #Take the top quark as the first parton/jet
                if genParticles[igen].pdgId() == 6 :
                    gen = ROOT.TLorentzVector()
                    gen.SetPtEtaPhiM( genParticles[igen].pt(), genParticles[igen].eta(), genParticles[igen].phi(), genParticles[igen].mass() )
                    genT = gen
                #Take the top antiquark as the second parton/jet
                elif genParticles[igen].pdgId() == -6 :
                    gen = ROOT.TLorentzVector()
                    gen.SetPtEtaPhiM( genParticles[igen].pt(), genParticles[igen].eta(), genParticles[igen].phi(), genParticles[igen].mass() )
                    genTbar = gen

            genJet1 = genT
            genJet2 = genTbar
            
            #If we have 2 valid genJets, save the Mtt spectrum which we need for splicing ttbar MC
            if genJet2.Pt() > 0 and genJet1.Pt() > 0:
                genPartonJetMtt = (genT+genTbar).M()
                self.genPartonJetMtt[0] = genPartonJetMtt
            
            #If we only want the <700 Mtt events, don't save anything else
            if self.invMassCut > 0 and genPartonJetMtt >= self.invMassCut:
                self.noFill = 1
                return
            else:
                self.noFill = 0

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

        #At this point we have not passed any reco cuts
        self.passKinCuts[0] = 0

        #############################################################################3
        #A note on response matrix weights. At this time, we aren't sure we have a jet, so we can't get the trigger efficiency.
        #This might make us off by a little.
        #############################################################################3
        self.unfoldWeightUsed[0]=weight

        #Two jets:
        if len(pj) < 2 or len(unpj) < 2 or len(ucPJ) < 2:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], weight)
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
                self.response.Miss(self.genPartonJet1pt[0], weight)
            return
        
        #Make sure we are within correct eta
        elif abs(pt_sorted_jets[0].Eta()) > 2.4 or abs(pt_sorted_jets[1].Eta()) > 2.4:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], weight)
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
        if pass400:
            self.pass400pt[0] = ca1.Pt()
        else:
            self.fail400pt[0] = ca1.Pt()
        if pass750:
            self.pass750pt[0] = ca1.Pt()
        else:
            self.fail750pt[0] = ca1.Pt()

        #Trigger SF. If data, do nothing. Otherwise you should be passed the appropriated histogram for nominal, up, or down.
        #Simply get the weight from it
        if self.isMC == False:
            self.triggerWeight[0] = 1.0
        elif self.useTrigger:
            bin = 0
            bin = self.trigger.FindBin(ca1.Pt())
            self.triggerWeight[0] = self.trigger.GetBinContent(bin)
            # print bin,ca1.Pt(),self.trigger.GetBinContent(bin),self.triggerWeight[0]
            if self.useTrigger > 1:
                # self.triggerWeight[0] = self.trigger.GetBinContent(bin) + self.trigger.GetBinError(bin)
                self.triggerWeight[0] = self.triggerWeight[0] + self.trigger.GetBinError(bin)
            elif self.useTrigger < 0:
                # self.triggerWeight[0] = self.trigger.GetBinContent(bin) - self.trigger.GetBinError(bin)
                self.triggerWeight[0] = self.triggerWeight[0] - self.trigger.GetBinError(bin)
            if self.triggerWeight[0] < 0.0:
                self.triggerWeight[0] = 0.0
        else:
            self.triggerWeight[0] = 1.0

        #Include trigger reweighting for the response matrix
        weight = weight * self.triggerWeight[0]
        self.unfoldWeightUsed[0]=weight

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
            return
        
        #Make sure matching is correct for jet2
        if jet2matchIndex==-1 or jet2matchIndex_pj==-1 or jet2matchIndex_ucPJ==-1:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
            return
        
        #Events with valid ca-other jet matches
        self.cutflow.Fill(5)

        # #Events with filled subjets
        # self.cutflow.Fill(6)

        #Nsubjettiness
        if Tau2[jet1matchIndex] == 0 or Tau2[jet2matchIndex] == 0:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
            return
        if Tau1[jet1matchIndex] == 0 or Tau1[jet2matchIndex] == 0:
            if self.doUnfold and self.isGenHadronic[0] == 1:
                self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
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

        #nSubjettiness and btag weights for response matrix
        numBtags = self.jet1bTagged[0] + self.jet2bTagged[0]
        # btagWt = math.pow(self.btagSF,numBtags)*math.pow(1-self.btagSF,2-numBtags)
        btagWt = math.pow(self.btagSF,numBtags)
        isNsubTag = 0
        if self.jet1tau32[0] < 0.55:
            isNsubTag = 1
        nsubWt = math.pow(self.nsubSF,isNsubTag)
        #Final weight to be put into response matrix. Shouldn't matter if btagged or pass/fail tau32
        weight = weight * self.unfoldWeight * btagWt * nsubWt
        self.unfoldWeightUsed[0]=weight

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
                    # self.response.Fill(self.jet1pt[0], self.genPartonJet1pt[0], self.unfoldWeight)
                    self.response.Fill(self.jet1pt[0], self.genPartonJet1pt[0], weight)
                else:
                    # self.response.Fake(self.genPartonJet1pt[0], self.unfoldWeight)
                    self.response.Fake(self.genPartonJet1pt[0], weight)
        elif self.doUnfold and self.isGenHadronic[0] == 1:
            # self.response.Miss(self.genPartonJet1pt[0], self.unfoldWeight)
            self.response.Miss(self.genPartonJet1pt[0], weight)

        return

    def reset(self):
        if not self.noFill:
            self.treeVars.Fill()

        self.run[0] = -1
        self.event[0] = -1
        self.lumi[0] = -1

        self.npv[0] = -1
        self.index[0] = -1
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
        self.nJets[0] = -1
        self.deltaY[0] = -10.0
        self.deltaPhi[0] = -10.0

        self.jet1bTagged[0] = -1
        self.jet2bTagged[0] = -1
        self.jet1topTagged[0] = -1
        self.jet2topTagged[0] = -1

        self.htSum[0] = -1.0
        self.triggerWeight[0] = 1.0
        self.pdfWeight[0] = 1.0
        self.pileupWeight[0] = 1.0
        self.pdfWeightNom[0] = 1.0
        self.pdfWeightUp[0] = 1.0
        self.pdfWeightDown[0] = 1.0
        self.unfoldWeightUsed[0] = 1.0

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

        self.pass400pt[0] = -1.0
        self.pass750pt[0] = -1.0
        self.fail400pt[0] = -1.0
        self.fail750pt[0] = -1.0

    def __del__(self):  
        # print self.out_info
        self.f.cd()
        if self.doUnfold == True:
            self.response.Write()
        self.f.Write()
        self.f.Close()
        if self.useTrigger:
            self.triggerFile.Close()
        if self.usePileup:
            self.pileupFile.Close()