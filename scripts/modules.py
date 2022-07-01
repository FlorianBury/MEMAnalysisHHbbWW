import os
import copy
import ROOT
from itertools import combinations

from memhep.computer import MEMComputer 
from memhep.types import Event

class DoubleLeptonChannel(MEMComputer):
    def __init__(self,args):
        super(DoubleLeptonChannel,self).__init__(args)
        self._loadHistograms()

    def _loadHistograms(self):
        flavourPath = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data','MEM_data','Btag_FR2.root')
        flavourFile = ROOT.TFile(flavourPath)
        
        ak4_b = flavourFile.Get("N_ak4_truth_bjets")
        ak4_c = flavourFile.Get("N_ak4_truth_cjets")
        ak4_l = flavourFile.Get("N_ak4_truth_lightjets")

        ak4_den = ak4_b.Clone()
        ak4_den.Add(ak4_c)
        ak4_den.Add(ak4_l)

        ak4_ratio = ak4_b.Clone()
        ak4_ratio.Divide(ak4_den)

        self.ak4_weight = copy.deepcopy(
        btagFile_b = F.Get('bjets_truth')
        h_c = F.Get('cjets_truth')
        h_l = F.Get('lightjets_truth')
        h_j = h_b.Clone()
        h_j.Add(h_c)
        h_j.Add(h_l)

        self.weight = copy.deepcopy(h_b)
        self.weight.Divide(h_l)
        F.Close()

    def _getJetBProb(self,eta,pt,btag):
        assert abs(eta) <= 2.5
        assert pt >= 0.
        assert btag >= 0. and btag <= 1.
        return self.weight.GetBinContent(self.weight.FindBin(abs(eta),pt,btag))

    def defineBranches(self):
        return {'considered_jets': int}

    def addBranchValues(self,event):
        if event.boosted1b_tag == 1:
            return {'considered_jets': 0}
        else:
            jets_E = [event.j1_E,event.j2_E,event.j3_E,event.j4_E]
            return {'considered_jets': sum([1 for E in jets_E if E>0])}

    def _addLeptons(self,event,memEvent):
        if event.l1_charge == -1:
            memEvent.addParticle('lepton1',event.l1_Px,event.l1_Py,event.l1_Pz,event.l1_E)
            memEvent.addParticle('lepton2',event.l2_Px,event.l2_Py,event.l2_Pz,event.l2_E)
        else:
            memEvent.addParticle('lepton1',event.l2_Px,event.l2_Py,event.l2_Pz,event.l2_E)
            memEvent.addParticle('lepton2',event.l1_Px,event.l1_Py,event.l1_Pz,event.l1_E)
       

    def processEvent(self,event):
        if False:
#        if event.boosted1b_tag == 1:
            memEvent = Event()
            self._addLeptons(event,memEvent)
            memEvent.addParticle('bjet1',event.fatj_sub1_Px,fatj_sub1_Py,fatj_sub1_Pz,fatj_sub1_E)
            memEvent.addParticle('bjet2',event.fatj_sub2_Px,fatj_sub2_Py,fatj_sub2_Pz,fatj_sub2_E)
            memEvent.addParticle('isr',
                   event.l1_Px+event.l2_Px+event.fatj_sub1_Px+event.fatj_sub2_Px,
                   event.l1_Py+event.l2_Py+event.fatj_sub1_Py+event.fatj_sub2_Py,
                   event.l1_Pz+event.l2_Pz+event.fatj_sub1_Pz+event.fatj_sub2_Pz,
                   event.l1_E+event.l2_E+event.fatj_sub1_E+event.fatj_sub2_E)

            return memEvent
            
        else:
            memEvents = {}
            idx_jets = [1,2,3,4]
            for pair in combinations(idx_jets,2):
                jet1 = f'j{pair[0]}' 
                jet2 = f'j{pair[1]}' 
                if getattr(event,f'{jet1}_E') == -9999. or getattr(event,f'{jet2}_E') == -9999.:
                    continue
                memEvent = Event()
                self._addLeptons(event,memEvent)
                memEvent.addParticle('bjet1',getattr(event,f'{jet1}_Px'),getattr(event,f'{jet1}_Py'),getattr(event,f'{jet1}_Pz'),getattr(event,f'{jet1}_E'))
                memEvent.addParticle('bjet2',getattr(event,f'{jet2}_Px'),getattr(event,f'{jet2}_Py'),getattr(event,f'{jet2}_Pz'),getattr(event,f'{jet2}_E'))
            
                memEvent.addParticle('isr',
                       event.l1_Px+event.l2_Px+getattr(event,f'{jet1}_Px')+getattr(event,f'{jet2}_Px'),
                       event.l1_Py+event.l2_Py+getattr(event,f'{jet1}_Py')+getattr(event,f'{jet2}_Py'),
                       event.l1_Pz+event.l2_Pz+getattr(event,f'{jet1}_Pz')+getattr(event,f'{jet2}_Pz'),
                       event.l1_E+event.l2_E+getattr(event,f'{jet1}_E')+getattr(event,f'{jet2}_E'))
                memEvent.setMET(event.met_Px,event.met_Py,event.met_Pz,event.met_E)

                memEvents[memEvent] = self._getJetBProb(getattr(event,f'{jet1}_eta'),
                                                        getattr(event,f'{jet1}_pt'),
                                                        getattr(event,f'{jet1}_btag')) * \
                                      self._getJetBProb(getattr(event,f'{jet2}_eta'),
                                                        getattr(event,f'{jet2}_pt'),
                                                        getattr(event,f'{jet2}_btag'))
                                        
            return memEvents


class ttbar_fullyLeptonic(DoubleLeptonChannel):
    @property
    def configFile(self):
        return '/home/ucl/cp3/fbury/HHMEM/MEMWeight/confs/ttbar_fullyLeptonic.lua'


class dy_to_llbb(DoubleLeptonChannel):
    @property
    def configFile(self):
        return '/home/ucl/cp3/fbury/HHMEM/MEMWeight/confs/dy_to_llbb.lua'


