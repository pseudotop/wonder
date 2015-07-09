from __future__ import division
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,sys,gc
from deltaR import *
import array

if not os.path.isdir("./h2mu"):
    os.mkdir("./h2mu")
os.chdir("./h2mu")
path = "/pnfs/user/jlee/cattuple/cat74v1/DYJetsToLL_M-50_13TeV-madgraph-pythia8/cat74v1_Phys14DR-PU20bx25_PHYS14_25_V1-v1/150616_191810/0000/"
f_list = os.listdir(path)
target = []
for x in f_list:
  if 'root' in x:
    target.append(x)
f_name = "test"   
out_rt = ROOT.TFile(f_name+".root", "RECREATE")
out_txt = open(f_name+".txt", "w")
############################################################
##################### HISTO SETUP ##########################
############################################################
bins_N = 19
# 0 to 80, 10 bins; 80 to 150, 7 bins; 150 to 200, 2bins;
Edges = [0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200]
Edges_bin = array.array('d', Edges)
eta_range = 3
eta_bin = 10*eta_range

h_resolution = []
h_resolution.append(ROOT.TH1F("resolution_hist_all", "resolution hist; (p_{t}^{gen}-p_{t}^{reco})/p_{t}^{gen}; number of entries", 100, -0.2, 0.2))
h_resolution.append(ROOT.TH1F("resolution_hist_endcap", "resolution hist; (p_{t}^{gen}-p_{t}^{reco})/p_{t}^{gen}; number of entries", 100, -0.2, 0.2))

h_geneta = ROOT.TH1F("eta_of_genMuons", "eta of genMuons; #eta; number of entries", eta_bin, -eta_range, eta_range)
h_recoeta = ROOT.TH1F("eta_of_recoMuons", "eta of recoMuons; #eta; number of entries", eta_bin, -eta_range, eta_range)
h_H2eta = ROOT.TH1F("eta_of_HtoMuons", "eta of H2Muons; #eta; number of entries", 30, -3, 3)

h_genpt = []
h_genpt.append(ROOT.TH1D("pT_of_genMuons_all", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_genpt.append(ROOT.TH1D("pT_of_genMuons_endcap", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
#h_genpt.append(ROOT.TH1D("pT_of_genMuons_ME0", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))

h_recopt = []
h_recopt.append(ROOT.TH1D("pT_of_recoMuons_all", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_recopt.append(ROOT.TH1D("pT_of_recoMuons_endcap", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
#h_recopt.append(ROOT.TH1D("pT_of_recoMuons_ME0", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_H2pt = ROOT.TH1F("pT_of_H2Muons", "pT of H2Muons; pT(GeV/c); number of entries", 30, 0, 200)

h_eff_pt = []
h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_all", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))
h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_endcap", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))
#h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_ME0", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))

h_eff_eta = ROOT.TH1F("eta_efficiency_of_recoMuons_per_genMuons", "eta-efficiency of recoMuons/genMuons; #eta; ratio", eta_bin, -eta_range, eta_range)
#h_eff2_eta = ROOT.TH1F("eta_efficiency_of_H2Muons_per_genMuons", "eta-efficiency of HiggsToMuons/genMuons; #eta; ratio", 30, -3, 3)
#h_eff2_pt = ROOT.TH1F("pT_efficiency_of_H2Muons_per_genMuons", "pT-efficiency of HiggsToMuons/genMuons; pT(GeV/c); ratio", 30, 0, 200)
###########################################################

pass_ps = [0,0,0,0,0,0]
for i,x in enumerate(target):
  path_tot = path + x
  events = Events(path_tot)
  percent = ((i+1)/len(target)*100)
  print "%d/%d -- %.1f%%" % (i+1,len(target),percent)

  cut_lep_pt = 20.
  cut_lep_eta = 2.1
  cut_lep_iso = 10.12
  cut_jet_pt = 30.
  cut_jet_eta = 4.7

  metsLabel, mets = "catMETs", Handle("std::vector<cat::MET>")
  jetsLabel, jets = "catJets", Handle("std::vector<cat::Jet>")
  muonsLabel, muons = "catMuons", Handle("std::vector<cat::Muon>")
  genParticlesLabel, genParticles = "prunedGenParticles", Handle("std::vector<reco::GenParticle>") 
  electronsLabel, electrons = "catElectrons", Handle("std::vector<cat::Electron>")
  vertexsLabel, vertexs = "offlineSlimmedPrimaryVertices", Handle("std::vector<reco::Vertex>")
  nmuons = 0
  njets = 0
  npass = 0
  chems=[0, 0, 0, 0, 0]
  chmms=[0, 0, 0, 0, 0]
  chees=[0, 0, 0, 0, 0]
  totaliev = 0
  gc.collect()
  for iev,event in enumerate(events):
      if iev%10000 == 0:
          print round(iev/50000.*100), "%"
      selectedelecs=[]
      selectedmuons=[]
      selectedjets=[]

      event.getByLabel(vertexsLabel, vertexs)
      vert = vertexs.product()[0]
      event.getByLabel(genParticlesLabel,genParticles)
      event.getByLabel(electronsLabel,electrons)
      event.getByLabel(muonsLabel,muons)
      genParticle_list = genParticles.product()
      electron_list = electrons.product()
      muon_list = muons.product()

      for g,m in enumerate(electron_list):
          if not m.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"):
              continue
          if not m.passConversionVeto():
              continue
          if not m.isPF():
              continue
          #if m.gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS') > 0:
          #    continue
          if m.pt() <= cut_lep_pt:
              continue
          if abs(m.scEta()) <= 1.4442:
              if m.relIso(0.3) >= 0.1649:
                  continue
          if abs(m.scEta()) >= 1.566:
              if m.relIso(0.3) >= 0.2075:
                  continue
          if abs(m.scEta()) > 1.4442 and abs(m.scEta()) < 1.566:
              continue
          if abs(m.eta()) >= 2.5:
              continue
          selectedelecs.append(m)

      for g,m in enumerate(muon_list):
          pass_ps[0] += 1
          if not m.isTightMuon():
              continue
          pass_ps[1] += 1
          #print m.pt()
          if m.pt() <= cut_lep_pt:
              continue
          pass_ps[2] += 1
          if abs(m.eta()) >= 2.4:
              continue
          pass_ps[3] += 1
          if m.relIso(0.4) >= 0.12:
              continue
          pass_ps[4] += 1
          selectedmuons.append(m)
          nmuons +=1

      for g in genParticle_list:
          if ( abs(g.pdgId())!= 13 or g.pt() <= cut_lep_pt): continue
          isfromHiggs = False
          for i in range(g.numberOfMothers()):
              if ( g.mother(i).pdgId()== 23 ):
                isfromHiggs = True
          if not isfromHiggs: continue
          h_geneta.Fill(g.eta())
          #all
          if (abs(g.eta()) < 2.4):
              h_genpt[0].Fill(g.pt())
          #endcap
          if (abs(g.eta()) < 2.4 and abs(g.eta()) > 1.5):
              h_genpt[1].Fill(g.pt())
          #me0
          #if (abs(g.eta()) < 3.0 and abs(g.eta()) > 2.5):
              #h_genpt[2].Fill(g.pt())     
          for m in selectedmuons:
               dr = deltaR(g.eta(), g.phi(), m.eta(), m.phi())
               if dr < 0.1:
                   h_recoeta.Fill(g.eta())
                   resolution = (m.pt()-g.pt())/g.pt()
                   #all
                   if (abs(g.eta()) < 2.4):
                     h_recopt[0].Fill(g.pt())
                     h_resolution[0].Fill(resolution)
                   #endcap
                   if (abs(g.eta()) < 2.4 and abs(g.eta()) > 1.5):
                     h_recopt[1].Fill(g.pt())
                     h_resolution[1].Fill(resolution)
                   break
  #step1
      lep1 = ROOT.TLorentzVector()
      lep2 = ROOT.TLorentzVector()
      ch=[]
      charge = 1
      if len(selectedmuons) + len(selectedelecs) == 2:
          if len(selectedmuons) == 1:
              lep1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.1056)
              lep2.SetPtEtaPhiM(selectedelecs[0].pt(), selectedelecs[0].eta(), selectedelecs[0].phi(), 0.00051)
              charge = selectedmuons[0].charge()*selectedelecs[0].charge()
              ch = "chem"
          elif len(selectedmuons) == 2:
              lep1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.1056)
              lep2.SetPtEtaPhiM(selectedmuons[1].pt(), selectedmuons[1].eta(), selectedmuons[1].phi(), 0.1056)
              charge = selectedmuons[0].charge()*selectedmuons[1].charge()
              ch = "chmm"
          else:
              lep1.SetPtEtaPhiM(selectedelecs[0].pt(), selectedelecs[0].eta(), selectedelecs[0].phi(), 0.00051)
              lep2.SetPtEtaPhiM(selectedelecs[1].pt(), selectedelecs[1].eta(), selectedelecs[1].phi(), 0.00051)
              charge = selectedelecs[0].charge()*selectedelecs[1].charge()
              ch = "chee"
          
      if (lep1+lep2).M() <= 20:
          continue
      if charge > 0:
          continue

      if ch == "chem": chems[0] +=1
      if ch == "chmm": chmms[0] +=1
      if ch == "chee": chees[0] +=1

  #step2
      if ch == []:
          continue
      if ch != "chem":
          if (lep1+lep2).M() > 76 and (lep1+lep2).M() < 106:
              continue

      if ch == "chem": chems[1] +=1
      if ch == "chmm": chmms[1] +=1
      if ch == "chee": chees[1] +=1

  #step3    
      event.getByLabel(jetsLabel,jets)
      for g,j in enumerate(jets.product()):
          if not j.LooseId():
              continue
          #if j.pileupJetId() < 0.9:
          #    continue
          if j.pt() <= cut_jet_pt:
              continue
          if abs(j.eta()) >= 2.4:
              continue 
          jet = ROOT.TLorentzVector(j.px(), j.py(), j.pz(), j.energy())
          dr1 = jet.DeltaR(lep1)
          dr2 = jet.DeltaR(lep2)
          if dr1 <= 0.4:
              continue
          if dr2 <= 0.4:
              continue
          selectedjets.append(j)
          njets += 1

      if len(selectedjets) < 2:
          continue

      if ch == "chem": chems[2] +=1
      if ch == "chmm": chmms[2] +=1
      if ch == "chee": chees[2] +=1

  #step4
      event.getByLabel(metsLabel,mets)
      met = mets.product()[0]

      if ch == "chem": chems[3] +=1

      if ch == "chmm" or ch == "chee":
          if met.pt() <= 40.:
              continue
      if ch == "chmm": chmms[3] +=1
      if ch == "chee": chees[3] +=1

  #step5
      btag = 0
      for g,j in enumerate(selectedjets):
          jets_CSVInclV2 = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")
          if jets_CSVInclV2 <= 0.814:
              continue
          btag += 1

      if btag == 0:
          continue
   
      if ch == "chem": chems[4] +=1
      if ch == "chmm": chmms[4] +=1
      if ch == "chee": chees[4] +=1

      """
      print "met pt", met.pt(), "met eta", met.eta(), "met phi", met.phi()

      if len(selectedmuons) < 2:
          continue

      for m in selectedmuons:
          print "muon pt", m.pt(), "muon eta", m.eta(), "muon charge", m.charge()
              
      if selectedmuons[0].charge()*selectedmuons[1].charge() == 1:
          continue
      mu1 = selectedmuons[0]
      mu2 = selectedmuons[1]
      
      higgs = ROOT.TLorentzVector(mu1.px(), mu1.py(), mu1.pz(), mu1.energy()) + ROOT.TLorentzVector(mu2.px(), mu2.py(), mu2.pz(), mu2.energy())

      print higgs.M()
      
      npass +=1
      """
      
  print
  print "<num of passed>"
  print "    ", "s1 ", "s2 ", "s3 ", "s4 ", "s5 "
  print "chee" , 
  for i in range(0,5):
      print chees[i],
  print
  print "chem" , 
  for i in range(0,5):
      print chems[i],
  print
  print "chmm" , 
  for i in range(0,5):
      print chmms[i],
  print
  print
print pass_ps
for i in range(2):
  h_resolution[i].Draw()

  h_eff_pt[i].SetStats(0)
  h_recopt[i].Sumw2()
  h_genpt[i].Sumw2()
  h_eff_pt[i].Divide(h_recopt[i], h_genpt[i], 1.0, 1.0, "B")
  h_eff_pt[i].Draw()

h_eff_eta.SetStats(0)
h_recoeta.Sumw2()
h_geneta.Sumw2()
h_eff_eta.Divide(h_recoeta, h_geneta, 1.0, 1.0, "B")
h_eff_eta.Draw()


out_txt.write("=" * 50)
out_txt.write("\n muons cut flow \n")
out_txt.write("=" * 50)
out_txt.write("\n total, isTightMuon, pt_cut : %d, eta_cut : 2.4, relIso(0.4) >= 0.12\n" % (cut_lep_pt))
out_txt.write("%d, %d, %d, %d, %d\n" % (pass_ps[0],pass_ps[1],pass_ps[2],pass_ps[3],pass_ps[4]))
out_rt.Write()
out_rt.Close()
out_txt.close()
