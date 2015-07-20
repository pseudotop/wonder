from __future__ import division
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,sys,gc
from deltaR import *
from copy import deepcopy
from array import array
import math

metsLabel, mets = "catMETs", Handle("std::vector<cat::MET>")
jetsLabel, jets = "catJets", Handle("std::vector<cat::Jet>")
muonsLabel, muons = "catMuons", Handle("std::vector<cat::Muon>")
genParticlesLabel, genParticles = "prunedGenParticles", Handle("std::vector<reco::GenParticle>") 
electronsLabel, electrons = "catElectrons", Handle("std::vector<cat::Electron>")
vertexsLabel, vertexs = "offlineSlimmedPrimaryVertices", Handle("std::vector<reco::Vertex>")
cut_lep_pt = 20.
cut_jet_pt = 30.
cut_jet_eta = 4.7

def muonSelect(m,toggle=False):
    if toggle:
        if not m.isMidiumMuon():
            return False
    if not toggle:
        if not m.isTightMuon():
            return False
    if m.pt() <= cut_lep_pt:
        return False
    if abs(m.eta()) >= 2.4:
        return False
    if m.relIso(0.4) >= 0.12:
        return False
    return True

def jetSelect(j):
    if not j.LooseId():
        return False
    if j.pileupJetId() < 0.9:
        return False
    if j.pt() < cut_jet_pt:
        return False
    if abs(j.eta()) > cut_jet_eta:
        return False
    return True
  
def jetcat_preSelect(selcjets, met):
    if len(selcjets) > 1:
        leadingjet = selcjets[0]
        subjet = selcjets[1]
        missingpt = met.product()[0].pt()
        if (leadingjet.pt()>40 and subjet.pt()>30 and missingpt<40):
            return 3
        #else :
        #    if (leadingjet.pt()>30) :
        #        return 2
        #    else :
        #        return 1
          
    if len(selcjets) == 1:
        #leadingjet = selcjets[0]
        #if (leadingjet.pt()>40) :
        return 2
        #else :
        #    return 1
    if len(selcjets) == 0:
        return 1
        
def jetcat_Select(selcjets, met, higgs):
    sel = jetcat_preSelect(selcjets, met)
    if sel==3:
        return jetcat_2(selcjets,higgs)
    if sel==2:
        return jetcat_1(higgs)
    if sel==1:
        return jetcat_0(higgs)


def jetcat_0(higgs):        
    if higgs.Pt()>10 :
        return 1
    if higgs.Pt()<=10 :
        return 2

def jetcat_1(higgs):
    if higgs.Pt()>10 :
        return 3
    if higgs.Pt()<=10 :
        return 4

def jetcat_2(selcjets,higgs):
    M_jets = ROOT.TLorentzVector(selcjets[0].px(), selcjets[0].py(), selcjets[0].pz(), selcjets[0].energy()) + ROOT.TLorentzVector(selcjets[1].px(), selcjets[1].py(), selcjets[1].pz(), selcjets[1].energy())
    delta_eta = selcjets[0].eta()-selcjets[1].eta()
    VBF_Tight = (M_jets.M() > 650 and abs(delta_eta) > 3.5)
    ggF_Tight = (M_jets.M() > 250 and higgs.Pt() > 50)
    if (VBF_Tight or ggF_Tight):
        if not ggF_Tight: #not ggF_Tight but only VBF_Tight
            return 5
        if not VBF_Tight: #also contrast of above
            return 6
        if (VBF_Tight and ggF_Tight):
            return 7
    else:
        return 8

def jetcat_GC(mu1, mu2):
    eta_mu = [abs(mu1.Eta()),abs(mu2.Eta())]
    GC = 0
    for i in range(2):
        if (eta_mu[i] < 0.8):
            GC += 1
        if (eta_mu[i] > 0.8 and eta_mu[i] < 1.6):
            GC += 10
        if (eta_mu[i] > 1.6 and eta_mu[i] < 2.4):
            GC += 100
    return GC



if not os.path.isdir("./h2mu_test"):
    os.mkdir("./h2mu_test")
os.chdir("./h2mu_test")
path = "/pnfs/user/jlee/cattuple/cat74v1/DYJetsToLL_M-50_13TeV-madgraph-pythia8/cat74v1_Phys14DR-PU20bx25_PHYS14_25_V1-v1/150616_191810/0000/"
f_list = os.listdir(path)
target = []
for x in f_list:
  if 'root' in x:
    target.append(x)
# toggle default is False.
# if toggle=True, selection generates 'isMidiumMuon'. else if toggle=False, 'isTightMuon'.
toggle = True

f_name = "test_muonselc"
out_rt = ROOT.TFile(f_name+".root", "RECREATE")
out_txt = open(f_name+".txt", "w")

br_v1= ["mu1_pt","mu1_eta","mu1_phi","mu1_energy","mu2_pt","mu2_eta","mu2_phi","mu2_energy",\
    "jet1_pt","jet1_eta","jet1_phi","jet1_energy","jet2_pt","jet2_eta","jet2_phi","jet2_energy",\
    "higgs_M","0jet", "1jet", "2jet", "0jet_GC", "1jet_GC",\
    "mid_mu1_pt","mid_mu1_eta","mid_mu1_phi","mid_mu1_energy","mid_mu2_pt","mid_mu2_eta","mid_mu2_phi","mid_mu2_energy",\
    "mid_jet1_pt","mid_jet1_eta","mid_jet1_phi","mid_jet1_energy","mid_jet2_pt","mid_jet2_eta","mid_jet2_phi","mid_jet2_energy",\
    "mid_higgs_M","mid_0jet", "mid_1jet", "mid_2jet", "mid_0jet_GC", "mid_1jet_GC"]
br_v2=["gen_mu_pt","gen_mu_eta","gen_mu_phi","gen_mu_energy","reco_mu_pt","reco_mu_eta","reco_mu_phi","reco_mu_energy",\
    "resolution",\
    "mid_reco_mu_pt","mid_reco_mu_eta","mid_reco_mu_phi","mid_reco_mu_energy",\
    "mid_resolution"]

br_l1 = []
br_l2 = []
tr1 = ROOT.TTree("h2mu","h2mu")
tr2 = ROOT.TTree("h2mu-genNreco","h2mu-genNreco")
    
for i,x in enumerate(br_v1):
    br_l1.append(array("d", [0.0]))
    tr1.Branch(x, br_l1[i], x+"/D")
for i,x in enumerate(br_v2):
    br_l2.append(array("d", [0.0]))
    tr2.Branch(x, br_l2[i], x+"/D")


pass_ps = [0,0,0,0,0,0]
mid_arr1 = 22
mid_arr2 = 9

for i,x in enumerate(target):
    if (i > 0):
        break
    path_tot = path + x
    events = Events(path_tot)
    percent = ((i+1)/len(target)*100)
    print "%d/%d -- %.1f%%" % (i+1,len(target),percent)

    nmuons = 0
    njets = 0
    npass = 0
    totaliev = 0
    gc.collect()
    for iev,event in enumerate(events):
        selectedmuons=[]
        mid_selectedmuons=[]
        selectedgenmuons=[]
        selectedjets=[]
        br_r1 = []
        for br in range(len(br_l1)):
          br_r1.append(-100)
       
        event.getByLabel(vertexsLabel, vertexs)
        event.getByLabel(genParticlesLabel,genParticles)
        event.getByLabel(muonsLabel,muons)
        event.getByLabel(jetsLabel,jets)
        event.getByLabel(metsLabel, mets)
        vert = vertexs.product()[0]
        genParticle_list = genParticles.product()
        muon_list = muons.product()
        jet_list = jets.product()
       
        for g,m in enumerate(muon_list):
            #midium
            if muonSelect(m,toggle):
                mid_selectedmuons.append(m)
            #tight
            if muonSelect(m):
                selectedmuons.append(m)
            nmuons +=1
      
        for ig, g in enumerate(genParticle_list):
            br_r2 = []
            for br2 in range(4):
                br_r2.append(-100)
            if ( abs(g.pdgId())!= 13 or g.pt() <= cut_lep_pt): continue
            isfromHiggs = False
            for i in range(g.numberOfMothers()):
                #pdgID : Z-boson=23, Higgs-boson=25
                if ( g.mother(i).pdgId()== 23 ):
                    isfromHiggs = True
            if not isfromHiggs: continue
            br_r2[0:3] = [g.pt(), g.eta(), g.phi(), g.energy()]
            for y in range(4):
                br_l2[y][0] = br_r2[y]
            #h_geneta.Fill(g.eta())
            #all
            #if (abs(g.eta()) < 2.4):
            #    h_genpt[0].Fill(g.pt())
            #endcap
            #if (abs(g.eta()) < 2.4 and abs(g.eta()) > 1.5):
            #    h_genpt[1].Fill(g.pt())
            #me0
            #if (abs(g.eta()) < 3.0 and abs(g.eta()) > 2.5):
                #h_genpt[2].Fill(g.pt())     
            for m in selectedmuons:
                dr = deltaR(g.eta(), g.phi(), m.eta(), m.phi())
                if dr < 0.1:
                    for y in range(4,9):
                        br_r2.append(-100)
                    #h_recoeta.Fill(g.eta())
                    resolution = (m.pt()-g.pt())/g.pt()
                    br_r2[4:8] = [g.pt(), g.eta(), g.phi(), g.energy(),resolution]
                    for y in range(4,9):
                        br_l2[y][0] = br_r2[y]
                    #all
                    #if (abs(g.eta()) < 2.4):
                     # h_recopt[0].Fill(g.pt())
                     # h_resolution[0].Fill(resolution)
                    #endcap
                    #if (abs(g.eta()) < 2.4 and abs(g.eta()) > 1.5):
                     # h_recopt[1].Fill(g.pt())
                     # h_resolution[1].Fill(resolution)
                    break
            if (len(mid_selectedmuons)!=0):
                for m in mid_selectedmuons:
                    dr = deltaR(g.eta(), g.phi(), m.eta(), m.phi())
                    if dr < 0.1:
                        for y in range(4+mid_arr2,9+mid_arr2):
                            br_r2.append(-100)
                        #h_recoeta.Fill(g.eta())
                        resolution = (m.pt()-g.pt())/g.pt()
                        br_r2[4+mid_arr2:8+mid_arr2] = [g.pt(), g.eta(), g.phi(), g.energy(),resolution]
                        for y in range(4+mid_arr2,9+mid_arr2):
                            br_l2[y][0] = br_r2[y]
                        break
            tr2.Fill()
        for j in jet_list:
            if jetSelect(j):
                selectedjets.append(j)
       
        if len(selectedmuons) < 2:
            continue
                
        if selectedmuons[0].charge()*selectedmuons[1].charge() == 1:
            continue
        mu1 = ROOT.TLorentzVector()
        mu2 = ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.105)
        mu2.SetPtEtaPhiM(selectedmuons[1].pt(), selectedmuons[1].eta(), selectedmuons[1].phi(), 0.105)
        br_r1[0:3] = [mu1.Pt(), mu1.Eta(), mu1.Phi(), mu1.Energy()]
        br_r1[4:7] = [mu2.Pt(), mu2.Eta(), mu2.Phi(), mu2.Energy()]
        
        
        #genmu1 = selectedgenmuons[0]
        #genmu2 = selectedgenmuons[1]
        #br_r[8:11] = [genmu1.pt(), genmu1.eta(), genmu1.phi(), genmu1.energy()]
        #br_r[12:15] = [genmu2.pt(), genmu2.eta(), genmu2.phi(), genmu2.energy()]
     
        if len(selectedjets) > 0:
            br_r1[8:11] = [selectedjets[0].pt(), selectedjets[0].eta(), selectedjets[0].phi(), selectedjets[0].energy()]
        if len(selectedjets) > 1:
            br_r1[12:15] = [selectedjets[1].pt(), selectedjets[1].eta(), selectedjets[1].phi(), selectedjets[1].energy()]

        higgs = mu1+mu2
        br_r1[16] = higgs.M()
      
        jetcat_Id = jetcat_Select(selectedjets, mets, higgs) 
        jetcat_range = 0 
        if jetcat_Id < 5:
          jetcat_range = jetcat_GC(mu1, mu2)
        if (jetcat_Id==1 or jetcat_Id==2):
            br_r1[17] = jetcat_Id
            br_r1[20] = jetcat_range
        if (jetcat_Id==3 or jetcat_Id==4):
            br_r1[18] = jetcat_Id
            br_r1[21] = jetcat_range
        if (jetcat_Id>4):
            br_r1[19] = jetcat_Id
        #print "jetcat_Id : ", jetcat_Id, "jetcat_range : ", jetcat_range   
        if (len(mid_selectedmuons)!=0):
            if len(mid_selectedmuons) < 2:
                continue
            if mid_selectedmuons[0].charge()*mid_selectedmuons[1].charge() == 1:
                continue
            mid_mu1 = ROOT.TLorentzVector()
            mid_mu2 = ROOT.TLorentzVector()
            mid_mu1.SetPtEtaPhiM(mid_selectedmuons[0].pt(), mid_selectedmuons[0].eta(), mid_selectedmuons[0].phi(), 0.105)
            mid_mu2.SetPtEtaPhiM(mid_selectedmuons[1].pt(), mid_selectedmuons[1].eta(), mid_selectedmuons[1].phi(), 0.105)
            br_r1[0+mid_arr1:3+mid_arr1] = [mid_mu1.Pt(), mid_mu1.Eta(), mid_mu1.Phi(), mid_mu1.Energy()]
            br_r1[4+mid_arr1:7+mid_arr1] = [mid_mu2.Pt(), mid_mu2.Eta(), mid_mu2.Phi(), mid_mu2.Energy()]
            if len(selectedjets) > 0:
                br_r1[8+mid_arr1:11+mid_arr1] = [selectedjets[0].pt(), selectedjets[0].eta(), selectedjets[0].phi(), selectedjets[0].energy()]
            if len(selectedjets) > 1:
                br_r1[12+mid_arr1:15+mid_arr1] = [selectedjets[1].pt(), selectedjets[1].eta(), selectedjets[1].phi(), selectedjets[1].energy()]
         
            higgs = mid_mu1+mid_mu2
            br_r1[16+mid_arr1] = higgs.M()
          
            jetcat_Id = jetcat_Select(selectedjets, mets, higgs) 
            jetcat_range = 0 
            if jetcat_Id < 5:
              jetcat_range = jetcat_GC(mid_mu1, mid_mu2)
            if (jetcat_Id==1 or jetcat_Id==2):
                br_r1[17+mid_arr1] = jetcat_Id
                br_r1[20+mid_arr1] = jetcat_range
            if (jetcat_Id==3 or jetcat_Id==4):
                br_r1[18+mid_arr1] = jetcat_Id
                br_r1[21+mid_arr1] = jetcat_range
            if (jetcat_Id>4):
                br_r1[19+mid_arr1] = jetcat_Id
        for y in range(len(br_l1)):
            br_l1[y][0] = br_r1[y]
        tr1.Fill()
        npass +=1
      
tr1.Print()
tr2.Print()

out_txt.write("=" * 50)
out_txt.write("\n muons cut flow \n")
out_txt.write("=" * 50)
out_txt.write("\n total, isTightMuon, pt_cut : %d, eta_cut : 2.4, relIso(0.4) >= 0.12\n" % (cut_lep_pt))
out_rt.Write()
out_rt.Close()
out_txt.close()
