import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
import math
from math import *
import sys

from Analysis.Tools.JetTools import *

class tree_maker:
    # def __init__(self, outputname, triggerFileStr, useTrigger, pileupFileStr, usePileup, usePDF, isMC, unfoldWeight, invMassCut, doUnfold):
    def __init__(self, outputname, triggerFileStr, useTrigger, pileupFileStr, usePileup, usePDF, isMC, unfoldWeight, invMassCut, doUnfold, syst, isMCatNLO):

        # load all the event info:
        # self.out_info = 0
        self.name = outputname
        self.triggerFileStr = triggerFileStr
        self.useTrigger = useTrigger
        self.pileupFileStr = pileupFileStr
        self.usePileup = usePileup
        # self.pdfFileStr = pdfFileStr
        # self.usePDF = usePDF
        self.isMC = isMC
        # self.unfoldWeight = unfoldWeight
        self.invMassCut = invMassCut
        # self.doUnfold = doUnfold
        self.isMCatNLO = isMCatNLO

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

        #Muons
        self.muonHandle = Handle( "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
        self.muonLabel = ( "jhuMuonPFlow", "muon")

        self.muonChargeHandle = Handle( "vector<int>" )
        self.muonChargeLabel = ( "jhuMuonPFlow", "muoncharge" )

        self.muonIsoHandle = Handle( "vector<double>" )
        self.muonIsoLabel = ( "jhuMuonPFlow", "muoniso" )

        self.muonIsLooseHandle = Handle( "vector<unsigned int>" )
        self.muonIsLooseLabel = ( "jhuMuonPFlow", "muonisloose" )

        self.muonIsTightHandle = Handle( "vector<unsigned int>" )
        self.muonIsTightLabel = ( "jhuMuonPFlow", "muonistight" )

        #MC@NLO weights
        self.geninfoHandle = Handle( "GenEventInfoProduct" )
        self.geninfoLabel = ( "generator" )

        self.__book__()

    def __book__(self):

        if (self.pileupFileStr == ''):
            self.usePileup = False
        else:
            self.pileupFile = ROOT.TFile(self.pileupFileStr + ".root")
            self.pileupFile.cd()
            self.pileup = self.pileupFile.Get("pileup").Clone()
            self.pileup.SetName('pileup')
            ROOT.SetOwnership( self.pileup, False )
            self.usePileup = True

        if (self.triggerFileStr != ''):
            self.triggerFile = ROOT.TFile(self.triggerFileStr + ".root")
            self.triggerFile.cd()
            self.trigger = self.triggerFile.Get("trigger").Clone()
            self.trigger.SetName('trigger')
            ROOT.SetOwnership( self.trigger, False )
            self.useTrigger = True

        print "Booking Histograms and Trees..."
        self.f = ROOT.TFile( self.name + ".root", "recreate" )
        self.f.cd()
        self.treeVars = ROOT.TTree('treeVars', 'treeVars')

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
        self.jet1tau1 = array('f', [-1.0])
        self.jet1tau2 = array('f', [-1.0])
        self.jet1tau3 = array('f', [-1.0])
        self.jet1tau4 = array('f', [-1.0])
        self.nJets = array('i', [-1])
        self.deltaY = array('f', [-10.0])
        self.deltaPhi = array('f', [-10.0])

        self.jet1tau32 = array('f', [-1.0])
        self.jet1tau31 = array('f', [-1.0])
        self.jet1tau21 = array('f', [-1.0])
        self.jet1nSubj = array('i', [-1])
       
        self.jet1topTagged = array('i', [-1])
        self.jet1bTagged = array('i', [-1])
        self.jet2bTagged = array('i', [-1])
        self.jet1bTaggedSub = array('i', [-1])
        self.jet2bTaggedSub = array('i', [-1])

        self.htSum = array('f', [-1.0])
        self.triggerWeight = array('f', [1.0])
        self.pileupWeight = array('f', [1.0])
        self.mcatnloWeight = array('f', [0.0])
        self.mcatnloWeightNorm = array('f', [0.0])

        self.passTriggerPt = array('f', [-1.0])
        self.failTriggerPt = array('f', [-1.0])
        
        self.muonPt = array('f', [-1.0])
        self.muonEta = array('f', [-10.0])
        self.muonPhi = array('f', [-10.0])
        self.muonMass = array('f', [-1.0])
        self.muonCharge = array('f', [-1.0])
        self.muonIso = array('f', [-1.0])
        self.muonIsLoose = array('i', [-1])
        self.muonIsTight = array('i', [-1])
        self.muonDr = array('f', [-1.0])
        self.muonPtRel = array('f', [-1.0])
        self.passMuon2D = array('i', [-1])
        self.nMuons = array('i', [-1])

        self.wPt = array('f', [-1.0])
        self.wEta = array('f', [-10.0])
        self.wPhi = array('f', [-10.0])
        self.wMass = array('f', [-1.0])
        self.wDr = array('f', [-10.0])
        self.lepTopPt = array('f', [-1.0])
        self.lepTopEta = array('f', [-1.0])
        self.lepTopPhi = array('f', [-1.0])
        self.lepTopMass = array('f', [-1.0])
        self.lepTopTagged = array('i', [-1])

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
        self.treeVars.Branch('jet1tau1', self.jet1tau1, 'jet1tau1/F')
        self.treeVars.Branch('jet1tau2', self.jet1tau2, 'jet1tau2/F')
        self.treeVars.Branch('jet1tau3', self.jet1tau3, 'jet1tau3/F')
        self.treeVars.Branch('jet1tau4', self.jet1tau4, 'jet1tau4/F')
        self.treeVars.Branch('jet1tau32', self.jet1tau32, 'jet1tau32/F')
        self.treeVars.Branch('jet1tau31', self.jet1tau31, 'jet1tau31/F')
        self.treeVars.Branch('jet1tau21', self.jet1tau21, 'jet1tau21/F')
        self.treeVars.Branch('jet1nSubj', self.jet1nSubj, 'jet1nSubj/I')
        self.treeVars.Branch('nJets', self.nJets, 'nJets/I')

        self.treeVars.Branch('deltaY', self.deltaY, 'deltaY/F')
        self.treeVars.Branch('deltaPhi', self.deltaPhi, 'deltaPhi/F')

        self.treeVars.Branch('jet1topTagged', self.jet1topTagged, 'jet1topTagged/I')
        self.treeVars.Branch('jet1bTagged', self.jet1bTagged, 'jet1bTagged/I')
        self.treeVars.Branch('jet2bTagged', self.jet2bTagged, 'jet2bTagged/I')
        self.treeVars.Branch('jet1bTaggedSub', self.jet1bTaggedSub, 'jet1bTaggedSub/I')

        self.treeVars.Branch('htSum', self.htSum, 'htSum/F')
        self.treeVars.Branch('triggerWeight', self.triggerWeight, 'triggerWeight/F')
        self.treeVars.Branch('pileupWeight', self.pileupWeight, 'pileupWeight/F')
        self.treeVars.Branch('mcatnloWeight', self.mcatnloWeight, 'mcatnloWeight/F')
        self.treeVars.Branch('mcatnloWeightNorm', self.mcatnloWeightNorm, 'mcatnloWeightNorm/F')

        self.treeVars.Branch('passTriggerPt', self.passTriggerPt, 'passTriggerPt/F')
        self.treeVars.Branch('failTriggerPt', self.failTriggerPt, 'failTriggerPt/F')

        self.treeVars.Branch('muonPt', self.muonPt, 'muonPt/F')
        self.treeVars.Branch('muonEta', self.muonEta, 'muonEta/F')
        self.treeVars.Branch('muonPhi', self.muonPhi, 'muonPhi/F')
        self.treeVars.Branch('muonMass', self.muonMass, 'muonMass/F')
        self.treeVars.Branch('muonCharge', self.muonCharge, 'muonCharge/I')
        self.treeVars.Branch('muonIso', self.muonIso, 'muonIso/F')
        self.treeVars.Branch('muonIsLoose', self.muonIsLoose, 'muonIsLoose/I')
        self.treeVars.Branch('muonIsTight', self.muonIsTight, 'muonIsTight/I')
        self.treeVars.Branch('muonDr', self.muonDr, 'muonDr/F')
        self.treeVars.Branch('muonPtRel', self.muonPtRel, 'muonPtRel/F')
        self.treeVars.Branch('passMuon2D', self.passMuon2D, 'passMuon2D/I')
        self.treeVars.Branch('nMuons', self.nMuons, 'nMuons/I')

        self.treeVars.Branch('wPt', self.wPt, 'wPt/F')
        self.treeVars.Branch('wEta', self.wEta, 'wEta/F')
        self.treeVars.Branch('wPhi', self.wPhi, 'wPhi/F')
        self.treeVars.Branch('wMass', self.wMass, 'wMass/F')
        self.treeVars.Branch('wDr', self.wDr, 'wDr/F')
        self.treeVars.Branch('lepTopPt', self.lepTopPt, 'lepTopPt/F')
        self.treeVars.Branch('lepTopEta', self.lepTopEta, 'lepTopEta/F')
        self.treeVars.Branch('lepTopPhi', self.lepTopPhi, 'lepTopPhi/F')
        self.treeVars.Branch('lepTopMass', self.lepTopMass, 'lepTopMass/F')
        self.treeVars.Branch('lepTopTagged', self.lepTopTagged, 'lepTopTagged/I')

        self.invarmass = array('f', [-1.0])
        self.jetangle = array('f', [-10.0])
        self.treeVars.Branch('invariant_mass', self.invarmass, 'invarmass/F')
        self.treeVars.Branch('angle_between_jets', self.jetangle, 'jetangle/F')

        self.cutflow = ROOT. TH1D("cutflow", "cutfow", 10, 0, 10 )
        self.cutflow.Sumw2()

    def analyze(self, event):
        #Fill trigger info
        event.getByLabel (self.triggerLabel, self.triggerHandle)
        trigNames = event.object().triggerNames(self.triggerHandle.product())

        trigPath = "HLT_Mu40_eta2p1"

        trigIndex = trigNames.triggerIndex(trigPath)

        for version in ["_v1","_v2","_v3","_v4","_v5","_v6","_v7","_v8","_v9","_v10","_v11"]:
            newPath = trigPath+version
            newIndex = trigNames.triggerIndex(newPath)
            if newIndex==trigNames.size():
                continue
            else:
                break

        trigIndex = trigNames.triggerIndex(newPath)

        passTrig = self.triggerHandle.product().accept(trigIndex)

        self.run[0] = event.object().id().run()
        self.event[0] = event.object().id().event()
        self.lumi[0] = event.object().id().luminosityBlock()

        event.getByLabel (self.npvLabel, self.npvHandle)
        event.getByLabel (self.metPtLabel, self.metPtHandle)
        event.getByLabel (self.metPhiLabel, self.metPhiHandle)
        npv = self.npvHandle.product()[0]
        metPt = self.metPtHandle.product()[0]
        metPhi = self.metPhiHandle.product()[0]

        event.getByLabel (self.prunedLabel, self.prunedHandle)
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

        event.getByLabel (self.muonLabel, self.muonHandle)
        event.getByLabel (self.muonChargeLabel, self.muonChargeHandle)
        event.getByLabel (self.muonIsoLabel, self.muonIsoHandle)
        event.getByLabel (self.muonIsLooseLabel, self.muonIsLooseHandle)
        event.getByLabel (self.muonIsTightLabel, self.muonIsTightHandle)


        CSVVals = self.CSVHandle.product()

        subjet1CSV = self.subjet1CSVHandle.product()
        subjet2CSV = self.subjet2CSVHandle.product()
        subjet3CSV = self.subjet3CSVHandle.product()
        subjet4CSV = self.subjet4CSVHandle.product()

        Tau1  =  self.t1Handle.product()
        Tau2  =  self.t2Handle.product()
        Tau3  =  self.t3Handle.product()
        Tau4  =  self.t4Handle.product()

        muons = self.muonHandle.product()
        muonCharge = self.muonChargeHandle.product()
        muonIso = self.muonIsoHandle.product()
        muonIsLoose = self.muonIsLooseHandle.product()
        muonIsTight = self.muonIsTightHandle.product()

        #Pileup Weight
        if self.isMC == False:
            self.pileupWeight[0] = 1.0
        elif self.usePileup:
            self.pileupWeight[0] = 0
            self.pileupWeight[0] = self.pileup.GetBinContent(self.pileup.FindBin(npv))
        else:
            self.pileupWeight[0] = 1.0

        #Include MC@NLO weights
        if self.isMCatNLO:
            event.getByLabel(self.geninfoLabel, self.geninfoHandle)
            geninfo = self.geninfoHandle.product()
            geninfo_weight = geninfo.weight()
            self.mcatnloWeight[0] = geninfo_weight

            if self.mcatnloWeight[0] < 0:
                self.mcatnloWeightNorm[0] = -1
            else:
                self.mcatnloWeightNorm[0] = 1


        if len(muons) == 0:
            return

        #Total events seen
        self.index[0]=0;
        self.cutflow.Fill(1);

        ##Making sure the first muon is the largest
        max_pt = muons[0].Pt()
        n_muons = 0
        for i in range(0,len(muons) ):
            temp_pt = muons[i].Pt()
            if (temp_pt > 45):
                n_muons = n_muons + 1
            if( temp_pt > max_pt):
                # print muons[i].Pt()
                print "Muon %i has greater pT (%f) than the first muon pT (%f)" %(i,temp_pt,max_pt)
                
        if n_muons == 0:
            return

        self.muonPt[0] = muons[0].Pt()
        self.muonEta[0] = muons[0].Eta()
        self.muonPhi[0] = muons[0].Phi()
        self.muonMass[0] = .1056583715


        self.muonCharge[0] = muonCharge[0]
        self.muonIso[0] = muonIso[0]
        self.muonIsLoose[0] = muonIsLoose[0]
        self.muonIsTight[0] = muonIsTight[0]

        self.nMuons[0] = n_muons

        if n_muons != 1:
            return

        # #Events with ONE valid muon
        self.cutflow.Fill(2)

        HTsum = 0

        unpj = self.unprunedHandle.product()
        pj = self.prunedHandle.product()

        self.nJets[0] = len(pj)

        # Two jets, including top tagged:
        if len(pj) < 2 or len(unpj) < 2:
            return
        # if self.nJets[0] < 2:
            # return
    
        #Events with 2 or more valid jets
        self.cutflow.Fill(3);

        #Make sure there's at least a top condidiate by checking the pT
        nTopCand = 0
        for i in range(0,len(pj) ) :
            if pj[i].pt() > 400. :
                nTopCand = nTopCand + 1
                HTsum += pj[i].pt()
        if nTopCand < 1:
            return
        
        #Events with 2 top candidates (pt > 400)
        self.cutflow.Fill(4);

        self.htSum[0] = HTsum
        self.MET[0] = metPt
        self.npv[0] = npv

        #Reorder to make sure the highest pT is first. This is after any JEC. Is this correct?
        #pt_sorted_jets = pj
        pt_sorted_jets = ReorderByPt(pj)

        #We need these for the substructure
        ca1 = ROOT.TLorentzVector()
        ca2 = ROOT.TLorentzVector()

        #We create the CA jet to match and fill the jet parameters          
        ca1.SetPtEtaPhiM(pt_sorted_jets[0].Pt(), pt_sorted_jets[0].Eta(), pt_sorted_jets[0].Phi(), pt_sorted_jets[0].M())
        ca2.SetPtEtaPhiM(pt_sorted_jets[1].Pt(), pt_sorted_jets[1].Eta(), pt_sorted_jets[1].Phi(), pt_sorted_jets[1].M())

        #Fill the trigger info
        if passTrig:
            self.passTriggerPt[0] = ca1.Pt()
        else:
            self.failTriggerPt[0] = ca1.Pt()
        
        #Fill trigger weight if present
        if self.isMC == False:
            self.triggerWeight[0] = 1.0
        elif self.useTrigger:
            bin = 0
            bin = self.trigger.FindBin(ca1.Pt())
            self.triggerWeight[0] = self.trigger.GetBinContent(bin)
        else:
            self.triggerWeight[0] = getMuonSF(self.muonEta[0])
            
        #Match unpruned jets with pruned - so we have both subjet btagging and nsubjettiness
        #This returns the index of the jet in the first collection that matches within dr = 0.4 to the jet of the second argument
        jet1matchIndex = MatchCol(unpj, ca1)
        jet2matchIndex = MatchCol(unpj, ca2)
        #This is temporary so we get the correct pj index for CSV and subjet CSV
        jet1matchIndex_pj = MatchCol(pj, ca1)
        jet2matchIndex_pj = MatchCol(pj, ca2)

        #Make sure matching is correct for jet1
        if jet1matchIndex==-1 or jet1matchIndex_pj==-1:
            return
        
        #Make sure matching is correct for jet2
        if jet2matchIndex==-1 or jet2matchIndex_pj==-1:
            return
        
        #Events with valid ca-other jet matches
        self.cutflow.Fill(5)

        #Nsubjettiness
        # if Tau2[jet1matchIndex] == 0 or Tau2[jet2matchIndex] == 0:
        #     return
        # if Tau1[jet1matchIndex] == 0 or Tau1[jet2matchIndex] == 0:
        #     return
        if Tau2[jet1matchIndex] == 0 or Tau1[jet1matchIndex] == 0:
            return
        

        self.jet1tau1[0] = Tau1[jet1matchIndex]
        self.jet1tau2[0] = Tau2[jet1matchIndex]
        self.jet1tau3[0] = Tau3[jet1matchIndex]
        self.jet1tau4[0] = Tau4[jet1matchIndex]

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

        #Fill our W candidate variables 
        newmet = ROOT.TLorentzVector()
        Wcand = ROOT.TLorentzVector()
        newmet.SetPtEtaPhiM(metPt,0,metPhi,0)
        phivec = [math.cos(metPhi), math.sin(metPhi)]
        P_muon = math.sqrt((muons[0].Px()*muons[0].Px())+(muons[0].Py()*muons[0].Py())+(muons[0].Pz()*muons[0].Pz()))
        P_phi = (muons[0].Px()*phivec[0])+(muons[0].Py()*phivec[1])
        b = (80.4*80.4) + (P_muon*P_muon) - (muons[0].E()*muons[0].E()) + (2*metPt*P_phi)
        Pz_met = 0
        # Pz_met = muons[0].Pz()*b/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
        # newmet.SetPz(Pz_met)
        # newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))

        #New 9/6 stolen from Marc after fixing phi bug
        arg = (muons[0].E()*muons[0].E()) * ((4*metPt*metPt*((muons[0].Pz()*muons[0].Pz())-(muons[0].E()*muons[0].E())))+(b*b))
        if arg <= 0:
                Pz_met = muons[0].Pz()*b/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
                newmet.SetPz(Pz_met)
                newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))
                # return [newmet+muons[0], newmet+muons[0]]
        else:
                Pz_met = ((muons[0].Pz()*b)+math.sqrt(arg))/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
                # Pz_met_m = ((muons[0].Pz()*b)-math.sqrt(arg))/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
                newmet.SetPz(Pz_met)
                newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))
                # newmet_m.SetPz(Pz_met_m)
                # newmet_m.SetE(math.sqrt(newmet_m.Px()*newmet_m.Px()+newmet_m.Py()*newmet_m.Py()+newmet_m.Pz()*newmet_m.Pz()))
                # return [newmet_p+muons[0], newmet_m+muons[0]]


        muon = ROOT.TLorentzVector()
        # jet2 = ROOT.TLorentzVector()
        muon.SetPtEtaPhiM(muons[0].Pt(),muons[0].Eta(),muons[0].Phi(),muons[0].M())
        # jet2.SetPtEtaPhiM(self.jet2pt[0],self.jet2eta[0],self.jet2phi[0],self.jet2mass[0])
        Wcand = newmet + muon
        topCand = Wcand + ca2

        # print 2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz()))
        # print 2*(muons[0].Pt()*muons[0].Pt())
        # print sqrt((muons[0].E()*muons[0].E()) - P_muon*P_muon)
        # print muons[0].M()

        # test_newmet = ROOT.TLorentzVector()
        # test_Wcand = ROOT.TLorentzVector()
        # test_newmet.SetPtEtaPhiM(metPt,0,metPhi,0)
        # test_Wcand = test_newmet + muon

        # print Wcand.E(),Wcand.Pt(),Wcand.Pz(),Wcand.M()
        # print test_Wcand.E(),test_Wcand.Pt(),test_Wcand.Pz(),test_Wcand.M()
        


        self.wPt[0] = Wcand.Pt()
        self.wEta[0] = Wcand.Eta()
        self.wPhi[0] = Wcand.Phi()
        self.wMass[0] = Wcand.M()
        #there isn't a MET eta, so let's take the dr of 2nd jet (hopefully btagged) wrt our Wcand
        # self.wDr[0] = reco.DeltaR(self.wEta[0], self.wPhi[0], self.jet2eta[0], self.jet2phi[0])
        self.wDr[0] = ca2.DeltaR(muon)

        self.lepTopPt[0] = topCand.Pt()
        self.lepTopEta[0] = topCand.Eta()
        self.lepTopPhi[0] = topCand.Phi()
        self.lepTopMass[0] = topCand.M()
        
        #2D cut
        muJetIndex=-1
        muJetIndex = ClosestJet(pj,muon)
        muonJet = ROOT.TLorentzVector()
        muonJet.SetPtEtaPhiM(pj[muJetIndex].Pt(),pj[muJetIndex].Eta(),pj[muJetIndex].Phi(),pj[muJetIndex].M())
        self.muonDr[0] = muon.DeltaR(muonJet)
        self.muonPtRel[0] = muon.Perp(muonJet.Vect())
        if self.muonDr[0] > 0.5 or self.muonPtRel[0] > 25.:
            self.passMuon2D[0] = 1
        else:
            self.passMuon2D[0] = 0
      
        #Invariant Mass         
        self.invarmass[0] = (ca1+ca2).M()

        #Angular Parameters
        self.jetangle[0] = ca1.DeltaR(ca2)
        self.deltaY[0] = (ca1.Rapidity() - ca2.Rapidity())
        self.deltaPhi[0] = ca1.DeltaPhi(ca2)
        
        #Nsubjettiness
        self.jet1tau32[0] = Tau3[jet1matchIndex] / Tau2[jet1matchIndex]
        self.jet1tau31[0] = Tau3[jet1matchIndex] / Tau1[jet1matchIndex]
        self.jet1tau21[0] = Tau2[jet1matchIndex] / Tau1[jet1matchIndex]
 
        #Fill the btagging information
        self.jet1csv[0] = CSVVals[jet1matchIndex_pj]
        self.jet2csv[0] = CSVVals[jet2matchIndex_pj]
        
        jet1subjetCSVs = []
        jet1subjetCSVs.append(subjet1CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet2CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet3CSV[jet1matchIndex_pj])
        jet1subjetCSVs.append(subjet4CSV[jet1matchIndex_pj])

        self.jet1maxSubjetCSV[0] = max(jet1subjetCSVs)

        self.jet1bTagged[0] = self.jet1csv[0] > 0.679
        self.jet2bTagged[0] = self.jet2csv[0] > 0.679

        self.jet1bTaggedSub[0] = self.jet1maxSubjetCSV[0] > 0.679
 
        self.jet1topTagged[0] = self.jet1mass[0] > 140.0 and self.jet1mass[0] < 250.0 and self.jet1pt[0] > 400.

        self.lepTopTagged[0] = self.lepTopMass[0] > 140.0 and self.lepTopMass[0] < 250.0 and self.jet2bTagged[0]

        if self.jet1topTagged[0] and self.lepTopTagged[0]: 
            self.index[0] = 1
            self.cutflow.Fill(8)

        self.treeVars.Fill()
        return

    def reset(self):

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
        self.jet1tau1[0] = -1.0
        self.jet1tau2[0] = -1.0
        self.jet1tau3[0] = -1.0
        self.jet1tau4[0] = -1.0
        self.jet1tau32[0] = -1.0
        self.jet1tau31[0] = -1.0
        self.jet1tau21[0] = -1.0
        self.jet1nSubj[0] = -1
        self.deltaY[0] = -10.0
        self.deltaPhi[0] = -10.0
        self.nJets[0] = -1
 
        self.jet1bTagged[0] = -1
        self.jet2bTagged[0] = -1
        self.jet1bTaggedSub[0] = -1
        self.jet2bTaggedSub[0] = -1
        self.jet1topTagged[0] = -1
 
        self.htSum[0] = -1.0
        self.triggerWeight[0] = 1.0
        self.pileupWeight[0] = 1.0
        self.mcatnloWeight[0] = 0.0
        self.mcatnloWeightNorm[0] = 0.0
        self.passTriggerPt[0] = -1.0
        self.failTriggerPt[0] = -1.0

        self.muonPt[0] = -1.0
        self.muonEta[0] = -10.0
        self.muonPhi[0] = -10.0
        self.muonMass[0] = -1.0
        self.muonCharge[0] = -1.0
        self.muonIso[0] = -1.0
        self.muonIsLoose[0] = -1
        self.muonIsTight[0] = -1
        self.muonDr[0] = -1.0
        self.muonPtRel[0] = -1.0
        self.passMuon2D[0] = -1
        self.nMuons[0] = -1

        self.wPt[0] = -1.0
        self.wEta[0] = -10.0
        self.wPhi[0] = -10.0
        self.wMass[0] = -1.0
        self.wDr[0] = -10.0

        self.lepTopPt[0] = -1.0
        self.lepTopEta[0] = -10.0
        self.lepTopPhi[0] = -10.0
        self.lepTopMass[0] = -1.0

        self.invarmass[0] = -1.0
        self.jetangle[0] = -10.0

    def __del__(self):  
        # print self.out_info
        self.f.cd()
        self.f.Write()
        self.f.Close()
        if self.usePileup:
            self.pileupFile.Close()

# Muon trigger * ID SF
def getMuonSF(muEta) :

    ## from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs#22Jan2013_ReReco_of_2012_data_re
    ## ID: https://indico.cern.ch/getFile.py/access?contribId=1&resId=2&materialId=slides&confId=257630
    ## trigger: https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=257000

    #Sal's numbers
    # muSF = 1.0
    # if abs(muEta) < 0.9 :
    #     muSF = 0.9827 * 0.9925
    # elif abs(muEta) < 1.2 :
    #     muSF = 0.9622 * 0.9928
    # else :
    #     muSF = 0.9906 * 0.9960

    #From here: https://indico.cern.ch/event/257000/contributions/1586340/attachments/451105/625487/SingleMuTriggerEff_10June_updated.pdf
    muSF = 1.0
    if abs(muEta) < 0.9 :
        muSF = 0.98299 * 0.9925
    elif abs(muEta) < 1.2 :
        muSF = 0.9604 * 0.9928
    else :
        muSF = 0.9955 * 0.9960


    return float(muSF)