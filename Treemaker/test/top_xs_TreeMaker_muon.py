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
        self.out_info = 0
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
        self.unprunedHandle = Handle("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > " )
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
        # self.jet2topMass = array('f', [-1.0])
        self.jet1csv = array('f', [-1.0])
        self.jet2csv = array('f', [-1.0])
        self.jet1maxSubjetCSV = array('f', [-1.0])
        # self.jet2maxSubjetCSV = array('f', [-1.0])
        self.jet1tau1 = array('f', [-1.0])
        self.jet1tau2 = array('f', [-1.0])
        self.jet1tau3 = array('f', [-1.0])
        self.jet1tau4 = array('f', [-1.0])
        # self.jet2tau1 = array('f', [-1.0])
        # self.jet2tau2 = array('f', [-1.0])
        # self.jet2tau3 = array('f', [-1.0])
        # self.jet2tau4 = array('f', [-1.0])
        self.deltaY = array('f', [-10.0])
        self.deltaPhi = array('f', [-10.0])

        self.jet1tau32 = array('f', [-1.0])
        # self.jet2tau32 = array('f', [-1.0])
        self.jet1tau31 = array('f', [-1.0])
        # self.jet2tau31 = array('f', [-1.0])
        self.jet1tau21 = array('f', [-1.0])
        # self.jet2tau21 = array('f', [-1.0])
        self.jet1nSubj = array('i', [-1])
        # self.jet2nSubj = array('i', [-1])
        self.jet1minMass = array('f', [-1.0])
        # self.jet2minMass = array('f', [-1.0])

        self.jet1bTagged = array('i', [-1])
        self.jet2bTagged = array('i', [-1])
        self.jet1bTaggedSub = array('i', [-1])
        # self.jet2bTaggedSub = array('i', [-1])
        self.jet1topTagged = array('i', [-1])
        # self.jet2topTagged = array('i', [-1])
        self.taggedPt = array('f', [-1.0])
        self.taggedMass = array('f', [-1.0])

        self.htSum = array('f', [-1.0])
        self.triggerEff = array('f', [1.0])

        self.jet1JEC = array('f', [-1.0])
        # self.jet2JEC = array('f', [-1.0])

        self.muonPt = array('f', [-1.0])
        self.muonEta = array('f', [-10.0])
        self.muonPhi = array('f', [-10.0])
        self.muonMass = array('f', [-1.0])
        self.muonCharge = array('f', [-1.0])
        self.muonIso = array('f', [-1.0])
        self.muonIsLoose = array('i', [-1])
        self.muonIsTight = array('i', [-1])

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
        # self.treeVars.Branch('jet2topMass', self.jet2topMass, 'jet2topMass/F')
        self.treeVars.Branch('jet1csv', self.jet1csv, 'jet1csv/F')
        self.treeVars.Branch('jet2csv', self.jet2csv, 'jet2csv/F')
        self.treeVars.Branch('jet1maxSubjetCSV', self.jet1maxSubjetCSV, 'jet1maxSubjetCSV/F')
        # self.treeVars.Branch('jet2maxSubjetCSV', self.jet2maxSubjetCSV, 'jet2maxSubjetCSV/F')
        self.treeVars.Branch('jet1tau1', self.jet1tau1, 'jet1tau1/F')
        self.treeVars.Branch('jet1tau2', self.jet1tau2, 'jet1tau2/F')
        self.treeVars.Branch('jet1tau3', self.jet1tau3, 'jet1tau3/F')
        self.treeVars.Branch('jet1tau4', self.jet1tau4, 'jet1tau4/F')
        # self.treeVars.Branch('jet2tau1', self.jet2tau1, 'jet2tau1/F')
        # self.treeVars.Branch('jet2tau2', self.jet2tau2, 'jet2tau2/F')
        # self.treeVars.Branch('jet2tau3', self.jet2tau3, 'jet2tau3/F')
        # self.treeVars.Branch('jet2tau4', self.jet2tau4, 'jet2tau4/F')
        self.treeVars.Branch('jet1tau32', self.jet1tau32, 'jet1tau32/F')
        # self.treeVars.Branch('jet2tau32', self.jet2tau32, 'jet2tau32/F')
        self.treeVars.Branch('jet1tau31', self.jet1tau31, 'jet1tau31/F')
        # self.treeVars.Branch('jet2tau31', self.jet2tau31, 'jet2tau31/F')
        self.treeVars.Branch('jet1tau21', self.jet1tau21, 'jet1tau21/F')
        # self.treeVars.Branch('jet2tau21', self.jet2tau21, 'jet2tau21/F')
        self.treeVars.Branch('jet1nSubj', self.jet1nSubj, 'jet1nSubj/I')
        # self.treeVars.Branch('jet2nSubj', self.jet2nSubj, 'jet2nSubj/I')
        self.treeVars.Branch('jet1minMass', self.jet1minMass, 'jet1minMass/F')
        # self.treeVars.Branch('jet2minMass', self.jet2minMass, 'jet2minMass/F')

        self.treeVars.Branch('deltaY', self.deltaY, 'deltaY/F')
        self.treeVars.Branch('deltaPhi', self.deltaPhi, 'deltaPhi/F')
        
        self.treeVars.Branch('jet1bTagged', self.jet1bTagged, 'jet1bTagged/I')
        self.treeVars.Branch('jet2bTagged', self.jet2bTagged, 'jet2bTagged/I')
        self.treeVars.Branch('jet1bTaggedSub', self.jet1bTaggedSub, 'jet1bTaggedSub/I')
        # self.treeVars.Branch('jet2bTaggedSub', self.jet2bTaggedSub, 'jet2bTaggedSub/I')
        self.treeVars.Branch('jet1topTagged', self.jet1topTagged, 'jet1topTagged/I')
        # self.treeVars.Branch('jet2topTagged', self.jet2topTagged, 'jet2topTagged/I')
        self.treeVars.Branch('taggedPt', self.taggedPt, 'taggedPt/F')
        self.treeVars.Branch('taggedMass', self.taggedMass, 'taggedMass/F')

        self.treeVars.Branch('htSum', self.htSum, 'htSum/F')
        self.treeVars.Branch('triggerEff', self.triggerEff, 'triggerEff/F')

        self.treeVars.Branch('jet1JEC', self.jet1JEC, 'jet1JEC/F')
        # self.treeVars.Branch('jet2JEC', self.jet2JEC, 'jet2JEC/F')

        self.treeVars.Branch('muonPt', self.muonPt, 'muonPt/F')
        self.treeVars.Branch('muonEta', self.muonEta, 'muonEta/F')
        self.treeVars.Branch('muonPhi', self.muonPhi, 'muonPhi/F')
        self.treeVars.Branch('muonMass', self.muonMass, 'muonMass/F')
        self.treeVars.Branch('muonCharge', self.muonCharge, 'muonCharge/I')
        self.treeVars.Branch('muonIso', self.muonIso, 'muonIso/F')
        self.treeVars.Branch('muonIsLoose', self.muonIsLoose, 'muonIsLoose/I')
        self.treeVars.Branch('muonIsTight', self.muonIsTight, 'muonIsTight/I')

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

        event.getByLabel (self.muonLabel, self.muonHandle)
        event.getByLabel (self.muonChargeLabel, self.muonChargeHandle)
        event.getByLabel (self.muonIsoLabel, self.muonIsoHandle)
        event.getByLabel (self.muonIsLooseLabel, self.muonIsLooseHandle)
        event.getByLabel (self.muonIsTightLabel, self.muonIsTightHandle)

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

        # if not (self.muonHandle.isValid()): # and self.muonChargeHandle.isValid()):
        #     return
        muons = self.muonHandle.product()
        muonCharge = self.muonChargeHandle.product()
        muonIso = self.muonIsoHandle.product()
        muonIsLoose = self.muonIsLooseHandle.product()
        muonIsTight = self.muonIsTightHandle.product()

        if len(muons) == 0:
            return

        #Total events seen
        self.index[0]=0
        self.cutflow.Fill(1)

        unpj = self.unprunedHandle.product()
        pj = self.prunedHandle.product()
        ttpj = self.topTaggedHandle.product()


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

        if n_muons != 1:
            return

        # #Events with ONE valid muon
        self.cutflow.Fill(2)

        HTsum = 0

        # newmet = ROOT.TLorentzVector()
        # Wcand = ROOT.TLorentzVector()
        # # newmet_m = ROOT.TLorentzVector()
        # # newmet_p = ROOT.TLorentzVector()
        # newmet.SetPtEtaPhiM(metPt,0,metPhi,0)
        # # newmet_m.SetPtEtaPhiM(metPt,0,metPhi,0)
        # # newmet_p.SetPtEtaPhiM(metPt,0,metPhi,0)
        # phivec = [math.cos(metPt), math.sin(metPhi)]
        # P_muon = math.sqrt((muons[0].Px()*muons[0].Px())+(muons[0].Py()*muons[0].Py())+(muons[0].Pz()*muons[0].Pz()))
        # P_phi = (muons[0].Px()*phivec[0])+(muons[0].Py()*phivec[1])
        # b = (80.4*80.4) + (P_muon*P_muon) - (muons[0].E()*muons[0].E()) + (2*metPt*P_phi)
        # # arg = (muons[0].E()*muons[0].E()) * ((4*met.Pt()*met.Pt()*((muons[0].Pz()*muons[0].Pz())-(muons[0].E()*muons[0].E())))+(b*b))
        # # if arg <= 0:
        # Pz_met = muons[0].Pz()*b/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
        # newmet.SetPz(Pz_met)
        # newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))

        # Wcand = newmet+muons[0]

        # return [newmet+lep, newmet+lep]


        # else:
        #         Pz_met_p = ((lep.Pz()*b)+math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
        #         Pz_met_m = ((lep.Pz()*b)-math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
        #         newmet_p.SetPz(Pz_met_p)
        #         newmet_p.SetE(math.sqrt(newmet_p.Px()*newmet_p.Px()+newmet_p.Py()*newmet_p.Py()+newmet_p.Pz()*newmet_p.Pz()))
        #         newmet_m.SetPz(Pz_met_m)
        #         newmet_m.SetE(math.sqrt(newmet_m.Px()*newmet_m.Px()+newmet_m.Py()*newmet_m.Py()+newmet_m.Pz()*newmet_m.Pz()))
        #         return [newmet_p+lep, newmet_m+lep]



        #Two normal jets and 1 top tagged:
        if len(pj) < 2 or len(unpj) < 2 or len(ttpj) < 1:
            # print len(pj),len(unpj),len(ttpj)
            return
    
        #Events with 2 or more valid jets
        self.cutflow.Fill(3);

        #Make sure there's at least a top condidiate by checking the pT
        #Reorder to make sure the highest pT is first. This is after any JEC. Is this correct?
        pt_sorted_ttJets = ReorderByPt(ttpj)
        pt_sorted_pjJets = ReorderByPt(pj)

        nTopCand = 0
        
        for i in range(0,len(pt_sorted_ttJets) ) :
            if pt_sorted_ttJets[i].Pt() > 400. :
                nTopCand = nTopCand + 1
                HTsum += pt_sorted_ttJets[i].Pt()
        if nTopCand < 1:
            return

        #Events with at least 1 hadronic top candidate (pt > 400)
        self.cutflow.Fill(4);

        self.htSum[0] = HTsum
        self.MET[0] = metPt
        self.npv[0] = npv


        #We need these for the substructure
        ca1 = ROOT.TLorentzVector()
        ca2 = ROOT.TLorentzVector()

        #We create the CA jet to match and fill the jet parameters          
        ca1.SetPtEtaPhiM(pt_sorted_ttJets[0].Pt(), pt_sorted_ttJets[0].Eta(), pt_sorted_ttJets[0].Phi(), pt_sorted_ttJets[0].M())
        ca2.SetPtEtaPhiM(pt_sorted_pjJets[1].Pt(), pt_sorted_pjJets[1].Eta(), pt_sorted_pjJets[1].Phi(), pt_sorted_pjJets[1].M())

        #Match unpruned jets with pruned - so we have both subjet btagging and nsubjettiness
        #This returns the index of the jet in the first collection that matches within dr = 0.4 to the jet of the second argument
        jet1matchIndex = MatchCol(unpj, ca1)
        jet2matchIndex = MatchCol(unpj, ca2)
        #This is temporary so we get the correct pj index for CSV and subjet CSV
        jet1matchIndex_pj = MatchCol(pj, ca1)
        jet2matchIndex_pj = MatchCol(pj, ca2)
        #Again temporary to match toptagged for Nsubjets, minMass, and topMass
        jet1matchIndex_tt = MatchCol(ttpj, ca1)
        # jet2matchIndex_tt = MatchCol(ttpj, ca2)

        #Make sure matching is correct for jet1
        if jet1matchIndex==-1 or jet1matchIndex_pj==-1 or jet1matchIndex_tt ==-1:
            return
        
        #Make sure matching is correct for jet2
        if jet2matchIndex==-1 or jet2matchIndex_pj==-1: # or jet2matchIndex_tt ==-1:
            return

        #Events with valid ca-other jet matches
        self.cutflow.Fill(5)

        #Make sure the toptagged jet is not somehow the 2nd pruned jet
        if jet1matchIndex_pj == jet2matchIndex_pj:
            return

        #Events without duplicate pruned jet
        self.cutflow.Fill(9)

        #Make sure the subjets are filled
        if len(nSubjets)==0:
            return

        #Events with filled subjets
        self.cutflow.Fill(6)

        #Nsubjettiness
        if Tau2[jet1matchIndex] == 0: # or Tau2[jet2matchIndex] == 0:
            return
        if Tau1[jet1matchIndex] == 0: # or Tau1[jet2matchIndex] == 0:
            return

        self.jet1tau1[0] = Tau1[jet1matchIndex]
        self.jet1tau2[0] = Tau2[jet1matchIndex]
        self.jet1tau3[0] = Tau3[jet1matchIndex]
        self.jet1tau4[0] = Tau4[jet1matchIndex]
        # self.jet2tau1[0] = Tau1[jet2matchIndex]
        # self.jet2tau2[0] = Tau2[jet2matchIndex]
        # self.jet2tau3[0] = Tau3[jet2matchIndex]
        # self.jet2tau4[0] = Tau4[jet2matchIndex]

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
       
        self.treeVars.Branch('lepTopPt', self.lepTopPt, 'lepTopPt/F')
        self.treeVars.Branch('lepTopMass', self.lepTopMass, 'lepTopMass/F')
        self.treeVars.Branch('lepTopTagged', self.lepTopTagged, 'lepTopTagged/I')


        #Fill our W candidate variables 
        newmet = ROOT.TLorentzVector()
        Wcand = ROOT.TLorentzVector()
        newmet.SetPtEtaPhiM(metPt,0,metPhi,0)
        phivec = [math.cos(metPt), math.sin(metPhi)]
        P_muon = math.sqrt((muons[0].Px()*muons[0].Px())+(muons[0].Py()*muons[0].Py())+(muons[0].Pz()*muons[0].Pz()))
        P_phi = (muons[0].Px()*phivec[0])+(muons[0].Py()*phivec[1])
        b = (80.4*80.4) + (P_muon*P_muon) - (muons[0].E()*muons[0].E()) + (2*metPt*P_phi)
        Pz_met = muons[0].Pz()*b/(2*((muons[0].E()*muons[0].E()) -(muons[0].Pz()*muons[0].Pz())))
        newmet.SetPz(Pz_met)
        newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))

        muon = ROOT.TLorentzVector()
        jet2 = ROOT.TLorentzVector()
        muon.SetPtEtaPhiM(muons[0].Pt(),muons[0].Eta(),muons[0].Phi(),muons[0].M())
        jet2.SetPtEtaPhiM(self.jet2pt[0],self.jet2eta[0],self.jet2phi[0],self.jet2mass[0])
        Wcand = newmet+muon
        topCand = Wcand + ca2

        self.wPt[0] = Wcand.Pt()
        self.wEta[0] = Wcand.Eta()
        self.wPhi[0] = Wcand.Phi()
        self.wMass[0] = Wcand.M()
        #there isn't a MET eta, so let's take the dr of 2nd jet (hopefully btagged) wrt our Wcand
        # self.wDr[0] = reco.DeltaR(self.wEta[0], self.wPhi[0], self.jet2eta[0], self.jet2phi[0])
        self.wDr[0] = jet2.DeltaR(muon)

        self.lepTopPt[0] = topCand.Pt()
        self.lepTopEta[0] = topCand.Eta()
        self.lepTopPhi[0] = topCand.Phi()
        self.lepTopMass[0] = topCand.M()

        #Top Tagging Parameters
        self.jet1nSubj[0] = nSubjets[jet1matchIndex_tt]
        # self.jet2nSubj[0] = nSubjets[jet2matchIndex_tt]
        self.jet1minMass[0] = minMass[jet1matchIndex_tt]
        # self.jet2minMass[0] = minMass[jet2matchIndex_tt]
        #### This should rightfully be jet1mass and jet2mass ####
        self.jet1topMass[0] = topTagTopMass[jet1matchIndex_tt]
        # self.jet2topMass[0] = topTagTopMass[jet2matchIndex_tt]
        
        self.jet1JEC[0] = topTaggedJEC[jet1matchIndex_tt]
        # self.jet2JEC[0] = topTaggedJEC[jet2matchIndex_tt]
        
        #Invariant Mass         
        self.invarmass[0] = (ca1+ca2).M()

        #Angular Parameters
        self.jetangle[0] = ca1.DeltaR(ca2)
        self.deltaY[0] = (ca1.Rapidity() - ca2.Rapidity())
        self.deltaPhi[0] = ca1.DeltaPhi(ca2)
        
        #Nsubjettiness
        self.jet1tau32[0] = Tau3[jet1matchIndex] / Tau2[jet1matchIndex]
        # self.jet2tau32[0] = Tau3[jet2matchIndex] / Tau2[jet2matchIndex]

        self.jet1tau31[0] = Tau3[jet1matchIndex] / Tau1[jet1matchIndex]
        # self.jet2tau31[0] = Tau3[jet2matchIndex] / Tau1[jet2matchIndex]

        self.jet1tau21[0] = Tau2[jet1matchIndex] / Tau1[jet1matchIndex]
        # self.jet2tau21[0] = Tau2[jet2matchIndex] / Tau1[jet2matchIndex]

        #Fill the btagging information
        self.jet1csv[0] = CSVVals[jet1matchIndex_pj]
        self.jet2csv[0] = CSVVals[jet2matchIndex_pj]
        
        jet1subjetCSVs = []
        jet1subjetCSVs.append(subjet1CSV[jet1matchIndex_tt])
        jet1subjetCSVs.append(subjet2CSV[jet1matchIndex_tt])
        jet1subjetCSVs.append(subjet3CSV[jet1matchIndex_tt])
        jet1subjetCSVs.append(subjet4CSV[jet1matchIndex_tt])

        # jet2subjetCSVs = []
        # jet2subjetCSVs.append(subjet1CSV[jet2matchIndex_tt])
        # jet2subjetCSVs.append(subjet2CSV[jet2matchIndex_tt])
        # jet2subjetCSVs.append(subjet3CSV[jet2matchIndex_tt])
        # jet2subjetCSVs.append(subjet4CSV[jet2matchIndex_tt])

        self.jet1maxSubjetCSV[0] = max(jet1subjetCSVs)
        # self.jet2maxSubjetCSV[0] = max(jet2subjetCSVs)

        self.jet1bTagged[0] = self.jet1csv[0] > 0.679
        self.jet2bTagged[0] = self.jet2csv[0] > 0.679

        self.jet1bTaggedSub[0] = self.jet1maxSubjetCSV[0] > 0.679
        # self.jet2bTaggedSub[0] = self.jet2maxSubjetCSV[0] > 0.679

        self.jet1topTagged[0] = self.jet1topMass[0] > 140.0 and self.jet1topMass[0] < 250.0 and self.jet1minMass[0] > 50.0 and self.jet1nSubj[0] > 2 and self.jet1pt[0] > 400.
        # self.jet2topTagged[0] = self.jet2topMass[0] > 140.0 and self.jet2topMass[0] < 250.0 and self.jet2minMass[0] > 50.0 and self.jet2nSubj[0] > 2 and self.jet2pt[0] > 400.

        self.lepTopTagged[0] = self.lepTopMass[0] > 140.0 and self.lepTopMass[0] < 250.0 and self.jet2bTagged[0]

        if self.jet1topTagged[0] and self.lepTopTagged[0]: 
            self.index[0] = 1
            self.cutflow.Fill(8)

        self.treeVars.Fill()


        # if self.doTrigger:
        #   triggerBin = self.trigger.FindBin(HTsum)
        #   if self.trigger.IsBinOverflow(triggerBin):
        #       triggerBin = self.trigger.GetNbinsX()
        #   elif triggerBin == 0:
        #       triggerBin = 1
            
        #   self.triggerEff = self.trigger.GetBinContent( triggerBin )

        return

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
        # self.jet2topMass[0] = -1.0
        self.jet1csv[0] = -1.0
        self.jet2csv[0] = -1.0
        self.jet1maxSubjetCSV[0] = -1.0
        # self.jet2maxSubjetCSV[0] = -1.0
        self.jet1tau1[0] = -1.0
        self.jet1tau2[0] = -1.0
        self.jet1tau3[0] = -1.0
        self.jet1tau4[0] = -1.0
        # self.jet2tau1[0] = -1.0
        # self.jet2tau2[0] = -1.0
        # self.jet2tau3[0] = -1.0
        # self.jet2tau4[0] = -1.0
        self.jet1tau32[0] = -1.0
        # self.jet2tau32[0] = -1.0
        self.jet1tau31[0] = -1.0
        # self.jet2tau31[0] = -1.0
        self.jet1tau21[0] = -1.0
        # self.jet2tau21[0] = -1.0
        self.jet1nSubj[0] = -1
        # self.jet2nSubj[0] = -1
        self.jet1minMass[0] = -1.0
        # self.jet2minMass[0] = -1.0
        self.deltaY[0] = -10.0
        self.deltaPhi[0] = -10.0
        
        self.jet1bTagged[0] = -1
        self.jet2bTagged[0] = -1
        self.jet1bTaggedSub[0] = -1
        # self.jet2bTaggedSub[0] = -1
        self.jet1topTagged[0] = -1
        # self.jet2topTagged[0] = -1
        self.taggedPt[0] = -1.0
        self.taggedMass[0] = -1.0

        self.htSum[0] = -1.0
        self.triggerEff[0] = -1.0

        self.jet1JEC[0] = -1.0
        # self.jet2JEC[0] = -1.0

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
        print self.out_info
        self.f.cd()
        self.f.Write()
        self.f.Close()
        if self.doModMass:
            self.modMassFile.Close()
        if self.doTrigger:
            self.triggerFile.Close()