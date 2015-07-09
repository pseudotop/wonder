import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,sys,gc
from deltaR import *
import array


    
def cloneHisto(ori_hist,new_hist_name,title,xmin,xmax):
    hist = ori_hist.Clone(new_hist_name)
    #e.g. title = "histogram title;x axis;y axis"
    hist.SetTitle(title)
    nx = xmax-xmin
    hist.SetBins(nx,xmin,xmax)
    return hist    

def createHisto(hist_name,title,xmin,xmax):
    bin = xmax-xmin
    histc = ROOT.TH1F(hist_name,title,bin,xmin,xmax)
    return histc

Phase_target = ["crab_h2mu_ggh_M125GeV_14TeV_2019WithGem_PU50", "crab_h2mu_ggh_M125GeV_14TeV_2019WithGem_PU140","crab_h2mu_ggh_M125GeV_14TeV_cfi_2023SHCal_PU140"]# "crab_h2mu_ggh_M125GeV_14TeV_2023WithGem_PU50", "crab_h2mu_ggh_M125GeV_14TeV_2023WithGem_PU140"]
Dataset_8TeV = "/pnfs/user/jlee/datagem/GluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/"
Dataset_catTuple = ["crab_h2mu_ggh_M125GeV_14TeV_cfi_2019WithGem1000fb_PU140_catTuple.root","crab_h2mu_ggh_M125GeV_14TeV_cfi_2019WithGem_PU50_catTuple.root","crab_h2mu_ggh_M125GeV_14TeV_cfi_2023SHCal_PU140_catTuple.root","crab_h2mu_VBF_M125GeV_14TeV_cfi_2019WithGem1000fb_PU140_catTuple.root","crab_h2mu_VBF_M125GeV_14TeV_cfi_2019WithGem_PU50_catTuple.root","crab_h2mu_VBF_M125GeV_14TeV_cfi_2023SHCal_PU140_catTuple.root", "crab_h2mu_Wh_M125GeV_14TeV_cfi_2019WithGem1000fb_PU140_catTuple.root","crab_h2mu_Wh_M125GeV_14TeV_cfi_2019WithGem_PU50_catTuple.root","crab_h2mu_Wh_M125GeV_14TeV_cfi_2023SHCal_PU140_catTuple.root", "crab_h2mu_Zh_M125GeV_14TeV_cfi_2019WithGem1000fb_PU140_catTuple.root","crab_h2mu_Zh_M125GeV_14TeV_cfi_2019WithGem_PU50_catTuple.root","crab_h2mu_Zh_M125GeV_14TeV_cfi_2023SHCal_PU140_catTuple.root", "crab_h2mu_tth_M125GeV_14TeV_cfi_2019WithGem1000fb_PU140_catTuple.root","crab_h2mu_tth_M125GeV_14TeV_cfi_2019WithGem_PU50_catTuple.root","crab_h2mu_tth_M125GeV_14TeV_cfi_2023SHCal_PU140_catTuple.root","GluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_catTuple.root","crab_DYToMuMu_M_110To140_TuneZ2star_14TeV_pythia6_cfi_2023SHCal_PU140_catTuple.root",""]
print("=" * 50)
print("higgs production.")
print("1. ggh")
print("2. VBF")
print("3. Wh")
print("4. Zh")
print("5. tth")
print("6. TuneZ2star")
sel = input("which higgs production do you want to use? : ")
if sel:
  print("=" * 50)
  print("1. %s" % Phase_target[0])
  print("2. %s" % Phase_target[1])
  print("3. %s" % Phase_target[2])
  #print("4. %s" % Phase_target[3])
  print("4. %s" % Dataset_8TeV)
  print("5. %s" % Dataset_catTuple[sel*3-3])
  print("6. %s" % Dataset_catTuple[sel*3-2])
  print("7. %s" % Dataset_catTuple[sel*3-1])
  print("=" * 50)
  print("ME0 exist in 2023 year.")
 
num = input("which program do you want to use? : ")

filename = "new"
#Events("DIR/h2mu_ggh----.root"), DIR list : Phase1_PU_140_1000fb_1  Phase1_PU_50  Phase2_PU_140  
#events = Events("/pnfs/user/h2mu_TP/Phase2_PU_140/h2mu_ggh_M125GeV_14TeV_Phase2_PU_140_RECO_10.root")

exist_me0 = False


path = ""
if num < (len(Phase_target)+1):
    if (num == 1 or num == 2):
      path = "/pnfs/user/jlee/datagem/%s/results/" % Phase_target[num-1]
    if (num == 3):
      path = "/pnfs/user/jlee/datagem/tp/%s/results/" % Phase_target[num-1]
      exist_me0 = True
    f_list = os.listdir(path)
    target = []
    if (num==1 or num==2):
      for x in f_list:
          if 'root' in x:
              tmp2 = x.split("_")[1]
              if tmp2 == "ggh":
                  target.append(x)
    if (num==3):
      for x in f_list:
          if 'root' in x:
              target.append(x)
    print target
    out_root_tot = ROOT.TFile(Phase_target[num-1]+"_"+filename+".root", "RECREATE")
    out_num_event_tot = open(Phase_target[num-1]+"_"+filename+"_num_of_events.txt","w")

if (num == 4):
    path = Dataset_8TeV
    f_list = os.listdir(path)
    target = []
    for x in f_list:
        if 'root' in x:
            target.append(x)
    print target
    out_root_tot = ROOT.TFile("Dataset_8TeV_noip"+"_"+filename+".root", "RECREATE")
    out_num_event_tot = open("Dataset_8TeV_noip"+"_"+filename+"_num_of_events.txt","w")

if (num > 4):
    path = "/pnfs/user/jlee/datagem/h2muCatuple/"
    target = []
    target.append(Dataset_catTuple[sel*3-3+num-5])
    print target
    out_root_tot = ROOT.TFile(Dataset_catTuple[sel*3-3+num-5]+"_"+filename+".root", "RECREATE")
    out_num_event_tot = open(Dataset_catTuple[sel*3-3+num-5]+"_"+filename+"_num_of_events.txt","w")

if not os.path.isdir("./h2mu_ggh"):
    os.mkdir("./h2mu_ggh")
os.chdir("./h2mu_ggh")


c_tot = ROOT.TCanvas("c_tot", "c_tot", 600, 400)

# histogram
higgs_bin = 400
hist_tot = ROOT.TH1F("hist_tot", "Higgs Mass Distribution", higgs_bin, 0., 200.)
hist_tot.GetXaxis().SetTitle("M_{H}(GeV/c^{2})")
hist_tot.GetYaxis().SetTitle("# of Higgs")
hist_tot.SetLineColor(2) #red
hist_tot.SetLineWidth(2)

h_barrel = ROOT.TH1F("h_barrel", "Higgs Mass Distribution", higgs_bin, 0., 200.)
h_barrel.GetXaxis().SetTitle("M_{H}(GeV/c^{2})")
h_barrel.GetYaxis().SetTitle("# of Higgs")
h_barrel.SetLineColor(2) #red
h_barrel.SetLineWidth(2)

h_endcap = ROOT.TH1F("h_endcap", "Higgs Mass Distribution", higgs_bin, 0., 200.)
h_endcap.GetXaxis().SetTitle("M_{H}(GeV/c^{2})")
h_endcap.GetYaxis().SetTitle("# of Higgs")
h_endcap.SetLineColor(2) #red
h_endcap.SetLineWidth(2)

h_cutflow = ROOT.TH1F("h_cutflow", "# of muon-candidates per each cut-processes", 11, 0, 11)
pr_num = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th", "11th"]
h_cutflow.GetXaxis().SetTitle("cut process")
h_cutflow.GetYaxis().SetTitle("remainder after cut")
h_cutflow.GetYaxis().SetTitleOffset(1.4)
h_cutflow.SetFillColor(38)
sbl1 = 0
while sbl1 < 11:
    sbl1 += 1
    h_cutflow.GetXaxis().SetBinLabel(sbl1, pr_num[sbl1-1])

hist2 = ROOT.TH1F("hist2", "# of muons per each event Distribution", 9, 0, 9)
mu_num = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
hist2.GetXaxis().SetTitle("# of muons per an event")
hist2.GetYaxis().SetTitle("# of events happening")
hist2.SetFillColor(38)
sbl = 0
while sbl < 9:
    sbl += 1
    hist2.GetXaxis().SetBinLabel(sbl, mu_num[sbl-1])

hist_tot2 = ROOT.TH1F("hist_tot2", "Higgs Mass Distribution +-5GeV Cut", higgs_bin, 0., 200.)
hist_tot2.GetXaxis().SetTitle("p_{t}(GeV/c)")
hist_tot2.GetYaxis().SetTitle("# of Higgs")
hist_tot2.SetLineColor(9) #dark purple?
hist_tot2.SetLineWidth(2)

h_resolution = []
h_resolution.append(ROOT.TH1F("resolution_hist_all", "resolution hist; (p_{t}^{gen}-p_{t}^{reco})/p_{t}^{gen}; number of entries", 100, -0.2, 0.2))
h_resolution.append(ROOT.TH1F("resolution_hist_endcap", "resolution hist; (p_{t}^{gen}-p_{t}^{reco})/p_{t}^{gen}; number of entries", 100, -0.2, 0.2))
h_resolution.append(ROOT.TH1F("resolution_hist_ME0", "resolution hist; (p_{t}^{gen}-p_{t}^{reco})/p_{t}^{gen}; number of entries", 100, -0.2, 0.2))

bins_N = 19
# 0 to 80, 10 bins; 80 to 150, 7 bins; 150 to 200, 2bins;
Edges = [0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200]
Edges_bin = array.array('d', Edges)
eta_range = 4
eta_bin = 10*eta_range


h_geneta = ROOT.TH1F("eta_of_genMuons", "eta of genMuons; #eta; number of entries", eta_bin, -eta_range, eta_range)
h_recoeta = ROOT.TH1F("eta_of_recoMuons", "eta of recoMuons; #eta; number of entries", eta_bin, -eta_range, eta_range)
h_H2eta = ROOT.TH1F("eta_of_HtoMuons", "eta of H2Muons; #eta; number of entries", 30, -3, 3)
h_genpt = []
h_genpt.append(ROOT.TH1D("pT_of_genMuons_all", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_genpt.append(ROOT.TH1D("pT_of_genMuons_endcap", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_genpt.append(ROOT.TH1D("pT_of_genMuons_ME0", "pT of genMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_recopt = []
h_recopt.append(ROOT.TH1D("pT_of_recoMuons_all", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_recopt.append(ROOT.TH1D("pT_of_recoMuons_endcap", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_recopt.append(ROOT.TH1D("pT_of_recoMuons_ME0", "pT of recoMuons; pT(GeV/c); number of entries", bins_N, Edges_bin))
h_H2pt = ROOT.TH1F("pT_of_H2Muons", "pT of H2Muons; pT(GeV/c); number of entries", 30, 0, 200)

h_eff_pt = []
h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_all", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))
h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_endcap", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))
h_eff_pt.append(ROOT.TH1D("pT_efficiency_of_recoMuons_per_genMuons_ME0", "pT-efficiency of recoMuons/genMuons; pT(GeV/c); Efficiency", bins_N, Edges_bin))
h_eff_eta = ROOT.TH1F("eta_efficiency_of_recoMuons_per_genMuons", "eta-efficiency of recoMuons/genMuons; #eta; ratio", eta_bin, -eta_range, eta_range)
#h_eff2_eta = ROOT.TH1F("eta_efficiency_of_H2Muons_per_genMuons", "eta-efficiency of HiggsToMuons/genMuons; #eta; ratio", 30, -3, 3)
#h_eff2_pt = ROOT.TH1F("pT_efficiency_of_H2Muons_per_genMuons", "pT-efficiency of HiggsToMuons/genMuons; pT(GeV/c); ratio", 30, 0, 200)

hist_l = []
hist_l.append(cloneHisto(hist_tot,"h_0jet_tight","'0' Jet Category-Tight;pT(#mu#mu) > 10GeV;# of di-muons",0,200)) 
hist_l.append(cloneHisto(hist_tot,"h_0jet_loose","'0' Jet Category-Loose;pT(#mu#mu) #leq 10GeV;# of di-muons",0,11)) 
hist_l.append(cloneHisto(hist_tot,"h_1jet_tight","'1' Jet Category-Tight;pT(#mu#mu) > 10GeV;# of di-muons",0,200))
hist_l.append(cloneHisto(hist_tot,"h_1jet_loose","'1' Jet Category-Tight;pT(#mu#mu) #leq 10GeV;# of di-muons",0,11))
hist_l.append(createHisto("h_2jet_VBF","'2' Jet Category-VBF Tight;M(jj) > 650GeV;# of di-muons",650,2000))
hist_l.append(createHisto("h_2jet_ggF","'2' Jet Category-ggF Tight;M(jj) > 250GeV;# of di-muons",250,650))
hist_l.append(createHisto("h_2jet_loose","'2' Jet Category-Loose; M(jj) GeV;# of di-muons",0,250))

h_higgs = []
h_higgs.append(cloneHisto(hist_tot,"higgs_0jet_tight","Higgs distribution for '0' Jet Category-Tight;M_{H}(GeV/c^{2});# of Higgs ",0,200)) 
h_higgs.append(cloneHisto(hist_tot,"higgs_0jet_loose","Higgs distribution for '0' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200)) 
h_higgs.append(cloneHisto(hist_tot,"higgs_1jet_tight","Higgs distribution for '1' Jet Category-Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h_higgs.append(cloneHisto(hist_tot,"higgs_1jet_loose","Higgs distribution for '1' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200))
h_higgs.append(cloneHisto(hist_tot,"higgs_2jet_VBF","Higgs distribution for '2' Jet Category-VBF Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h_higgs.append(cloneHisto(hist_tot,"higgs_2jet_ggF","Higgs distribution for '2' Jet Category-ggF Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h_higgs.append(cloneHisto(hist_tot,"higgs_2jet_loose","Higgs distribution for '2' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200))

h2_higgs = []
h2_higgs.append(cloneHisto(hist_tot,"higgs2_0jet_tight","Higgs distribution for '0' Jet Category-Tight;M_{H}(GeV/c^{2});# of Higgs ",0,200)) 
h2_higgs.append(cloneHisto(hist_tot,"higgs2_0jet_loose","Higgs distribution for '0' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200)) 
h2_higgs.append(cloneHisto(hist_tot,"higgs2_1jet_tight","Higgs distribution for '1' Jet Category-Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h2_higgs.append(cloneHisto(hist_tot,"higgs2_1jet_loose","Higgs distribution for '1' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200))
h2_higgs.append(cloneHisto(hist_tot,"higgs2_2jet_VBF","Higgs distribution for '2' Jet Category-VBF Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h2_higgs.append(cloneHisto(hist_tot,"higgs2_2jet_ggF","Higgs distribution for '2' Jet Category-ggF Tight;M_{H}(GeV/c^{2});# of Higgs",0,200))
h2_higgs.append(cloneHisto(hist_tot,"higgs2_2jet_loose","Higgs distribution for '2' Jet Category-Loose;M_{H}(GeV/c^{2});# of Higgs",0,200))

geo_Cat = []
geo_Cat_list = []
cat_count = 0
cat_name = ["higgs_0jet_tight", "higgs_0jet_loose", "higgs_1jet_tight", "higgs_1jet_loose"]
cat_title = ["'0' Jet Category-Tight", "'0' Jet Category-Loose", "'1' Jet Category-Tight", "'1' Jet Category-Loose"]
while cat_count < 4:
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_BB", cat_title[cat_count]+" : BB"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_BO", cat_title[cat_count]+" : BO"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_BE", cat_title[cat_count]+" : BE"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_OO", cat_title[cat_count]+" : OO"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_OE", cat_title[cat_count]+" : OE"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat_list.append(cloneHisto(hist_tot, cat_name[cat_count]+"_EE", cat_title[cat_count]+" : EE"+";M_{H}(GeV/c^{2});# of Higgs ",0,200))
  geo_Cat.append(geo_Cat_list)
  cat_count += 1
cut_muon_pt = 15.
cut_muon_eta = 2.4
cut_muon_iso = 10.12
cut_jet_pt = 30.
cut_jet_eta = 4.7

process = [0,0,0,0,0,0,0,0,0,0,0,0,0]
higgs_proc = [0,0,0,0,0,0,0]
higgs_5cut = [0,0,0,0,0,0,0]
cut = [0,0,0,0,0,0,0,0,0,0]
sum_tot = [0,0,0,0,0,0]
sum_iev = 0
sum_a = 0
geo_Count = []
count = 0
while count < 4:
  geo_Count.append([0,0,0,0,0,0])
  count += 1
for giev,x in enumerate(target):
    path_tot = path + x
    #print("=" * 50) 
    #print path_tot
    #print("=" * 50) 
    events = Events(path_tot)
    print "%d/%d" % (giev+1,len(target))
    #div_dir = path_tot.split('/')
    #out_f = div_dir[-1] 
    #out_root = ROOT.TFile(out_f,"RECREATE")
    #out_num_event = open(out_f[:-5]+"_num_of_events.txt","w")

    # histogram
    #hist = ROOT.TH1F("hist", "Higgs Mass Distribution", 200, 0., 200.)
    #hist.GetXaxis().SetTitle("GeV")
    #hist.GetYaxis().SetTitle("# of Higgs")

    metLabel, met = "catMETs", Handle("std::vector<cat::MET>")
    #jetsLabel, jets = "ak5PFJets", Handle("std::vector<reco::PFJet>")
    jetsLabel, jets = "catJets", Handle("std::vector<cat::Jet>")
    muonsLabel, muons = "catMuons", Handle("std::vector<cat::Muon>")
    #muonsLabel, muons = "muons", Handle("std::vector<reco::Muon>")
    #vertexsLabel, vertexs = "offlinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    #if (sel == 6 and num == 5):
    #  vertexsLabel, vertexs = "goodOfflinePrimaryVertices", Handle("std::vector<cat::Vertex>")
    #else:
    vertexsLabel, vertexs = "goodOfflinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    genParticlesLabel, genParticles = "genParticles", Handle("std::vector<reco::GenParticle>")
    if not (sel == 6 and num == 5):
      me0Label, me0 = "me0SegmentMatching", Handle("std::vector<reco::ME0Muon>")
    
    nmuons = 0
    njets = 0
    npass = 0
    n2muons = 0
    nopp_charge = 0
    nhiggs_cut = 0
    gc.collect()
    for iev,event in enumerate(events):
        #print iev
        #not included ME0
        selectedmuons=[]
        #included ME0
        selectedmuons_ME0=[]
        selectedjets=[]

        event.getByLabel(vertexsLabel, vertexs)
        event.getByLabel(genParticlesLabel, genParticles)
        event.getByLabel(muonsLabel,muons)
        event.getByLabel(jetsLabel,jets)
        event.getByLabel(metLabel,met)
        #vert = vertexs.product()[0]
        #print met.product()[0].pt()
        genParticle_list = genParticles.product()
        muon_list = muons.product()
        #if not vert:
        #    print "no vertax"

        a = [0,0,0,0,0]
        for j,m in enumerate(muon_list):
            #if m.isPFMuon() and ( m.isGlobalMuon() or m.isTrackerMuon()):
            sum_a += 1
            cut[0] += 1
            if not m.isTightMuon():
                continue
            cut[1] += 1
            if m.pt() < cut_muon_pt:
                continue
            cut[2] += 1
            if abs(m.eta()) > cut_muon_eta:
                continue
            cut[3] += 1
            a[4] += 1
            selectedmuons.append(m)
            selectedmuons_ME0.append(m)
            nmuons +=1
           
        #ME0 without 2019 
        if (exist_me0):
          event.getByLabel(me0Label,me0)
          for j,m in enumerate(me0.product()):
              if m.pt() < 5.:
                  continue
              #selectedmuons.append(m)
              selectedmuons_ME0.append(m)
        #genparticle compare with recomuons 
        for g in genParticle_list :
          if ( abs(g.pdgId())!=13 or g.pt() < cut_muon_pt ): continue 
          isfromHiggs = False
          for i in range(g.numberOfMothers()):
            if ( g.mother(i).pdgId()==25):
              isfromHiggs = True
              h_H2eta.Fill(g.eta())
              h_H2pt.Fill(g.pt())
          if not isfromHiggs: continue
          h_geneta.Fill(g.eta())
          #all
          if (abs(g.eta()) < 2.4):
            h_genpt[0].Fill(g.pt())
          #endcap
          if (abs(g.eta()) < 2.4 and abs(g.eta()) > 1.5):
            h_genpt[1].Fill(g.pt())
          #me0
          if (abs(g.eta()) < 3.0 and abs(g.eta()) > 2.5):
            h_genpt[2].Fill(g.pt())
          if (exist_me0):
            event.getByLabel(me0Label,me0)
            for k,n in enumerate(me0.product()):
              if n.pt() < 5.:
                continue
              resolution = (n.pt()-g.pt())/g.pt()
              dr = deltaR(g.eta(), g.phi(), n.eta(), n.phi())
              if dr < 0.1:
                if (abs(g.eta()) < 3.0 and abs(g.eta()) > 2.5):
                  h_recopt[2].Fill(g.pt())
                  h_resolution[2].Fill(resolution)
          
          for m in selectedmuons_ME0 :
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
        # the number of muons per event#  
        ev_mu = a[4]
        hist2.Fill(mu_num[ev_mu],1)
        ################################
        #print a
        #print cut
        for g,j in enumerate(jets.product()):
            cut[4] += 1
            if not j.LooseId():
                continue
            cut[5] += 1
            if j.pileupJetId() < 0.9:
                continue
            cut[6] += 1
            if j.pt() < cut_jet_pt:
                continue
            cut[7] += 1
            if abs(j.eta()) > cut_jet_eta:
                continue
            cut[8] += 1
            selectedjets.append(j)
            njets += 1
        if len(selectedmuons) < 2:
            continue
        n2muons += 1
        #for m in selectedmuons:
        #    print "muon pt", m.pt(), "muon eta", m.eta(), "muon charge", m.charge()
            
        if selectedmuons[0].charge()*selectedmuons[1].charge() == 1:
            continue
        nopp_charge += 1
        mu1 = ROOT.TLorentzVector()
        mu2 = ROOT.TLorentzVector()
        mu1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.105)
        mu2.SetPtEtaPhiM(selectedmuons[1].pt(), selectedmuons[1].eta(), selectedmuons[1].phi(), 0.105)
        higgs = mu1+mu2

        #print higgs.M()
        h_M = higgs.M()
        
        if h_M > 120 and h_M < 130:
            higgs_cut = higgs.M()
            hist_tot2.Fill(higgs_cut)
            if higgs_cut:
                nhiggs_cut += 1
        #### Pre-Selection Cut ####
        jetcategory = []
        jetcate_nomiss = []
        if len(selectedjets) > 1:
          process[0] += 1  
          if selectedjets[0].pt() > 40 :
            leadingjet = selectedjets[0] #already sorted by pt.
            jetcategory.append(leadingjet)
            jetcate_nomiss.append(leadingjet)
            process[1] += 1
          if selectedjets[1].pt() > 30 :
            subjet = selectedjets[1]
            jetcategory.append(subjet)
            jetcate_nomiss.append(subjet)
            process[2] += 1
          missingpt = met.product()[0].pt()
          if missingpt < 40 :
            jetcategory.append(missingpt)
            process[3] += 1
        if len(selectedjets)==1:
          process[4] += 1
          if selectedjets[0].pt() > 40 :
            leadingjet = selectedjets[0] #already sorted by pt.
            jetcategory.append(leadingjet)
            jetcate_nomiss.append(leadingjet)
          missingpt = met.product()[0].pt()
          if missingpt < 40:
            jetcategory.append(missingpt)
                        
        ##########################
        eta_mu = [abs(mu1.Eta()),abs(mu2.Eta())]
        B = []
        O = []
        E = []
        GC = [False, False, False, False, False, False]
        for i in range(2):
          B.append(eta_mu[i] < 0.8)
          O.append(eta_mu[i] > 0.8 and eta_mu[i] < 1.6)
          E.append(eta_mu[i] > 1.6 and eta_mu[i] < 2.4)
        if (B[0] and B[1]) :
          GC[0] = True
        if (B[0] and O[1]) or (B[1] and O[0]) :
          GC[1] = True  
        if (B[0] and E[1]) or (B[1] and E[0]) :  
          GC[2] = True
        if (O[0] and O[1]) :  
          GC[3] = True
        if (O[0] and E[1]) or (O[1] and E[0]) :  
          GC[4] = True
        if (E[0] and E[1]) :  
          GC[5] = True
          
       
        #### 0,1 Jet Category ####
        if len(jetcate_nomiss)==0:
          if (higgs.Pt()>10):
            process[5] += 1
            if h_M:
              h_higgs[0].Fill(h_M)
              higgs_proc[0] += 1
            if (h_M > 120 and h_M < 130):
              h2_higgs[0].Fill(h_M)
              higgs_5cut[0] += 1
            for i in range(6):
              if GC[i] :
                geo_Cat[0][i].Fill(h_M)
                geo_Count[0][i] += 1
            hist_l[0].Fill(higgs.Pt())  
          if (higgs.Pt()<=10):
            process[6] += 1
            if h_M:
              h_higgs[1].Fill(h_M)
              higgs_proc[1] += 1
            if (h_M > 120 and h_M < 130):
              h2_higgs[1].Fill(h_M)
              higgs_5cut[1] += 1
            for i in range(6):
              if GC[i] :
                geo_Cat[1][i].Fill(h_M)
                geo_Count[1][i] += 1
            hist_l[1].Fill(higgs.Pt())
        if len(jetcate_nomiss)==1:
          if (higgs.Pt()>10):
            process[7] += 1
            if h_M:
              h_higgs[2].Fill(h_M)
              higgs_proc[2] += 1
            if (h_M > 120 and h_M < 130):
              h2_higgs[2].Fill(h_M)
              higgs_5cut[2] += 1
            for i in range(6):
              if GC[i] :
                geo_Cat[2][i].Fill(h_M)
                geo_Count[2][i] += 1
            hist_l[2].Fill(higgs.Pt())
          if (higgs.Pt()<=10):
            process[8] += 1
            if h_M:
              h_higgs[3].Fill(h_M)
              higgs_proc[3] += 1
            if (h_M > 120 and h_M < 130):
              h2_higgs[3].Fill(h_M)
              higgs_5cut[3] += 1
            for i in range(6):
              if GC[i] :
                geo_Cat[3][i].Fill(h_M)
                geo_Count[3][i] += 1
            hist_l[3].Fill(higgs.Pt())
        ##########################
        
        ####  2 Jet Category  ####
        if len(jetcategory) > 2: #missing.pt() is/are included.
          process[9] += 1
          M_jets = ROOT.TLorentzVector(jetcategory[0].px(), jetcategory[0].py(), jetcategory[0].pz(), jetcategory[0].energy()) + ROOT.TLorentzVector(jetcategory[1].px(), jetcategory[1].py(), jetcategory[1].pz(), jetcategory[1].energy())
          delta_eta = jetcategory[0].eta()-jetcategory[1].eta()
          VBF_Tight = (M_jets.M() > 650 and abs(delta_eta) > 3.5)
          ggF_Tight = (M_jets.M() > 250 and higgs.Pt() > 50)
          if VBF_Tight or ggF_Tight:
          #VBF Tight
            if VBF_Tight:
              if h_M:
                h_higgs[4].Fill(h_M)
                higgs_proc[4] += 1
              if (h_M > 120 and h_M < 130):
                h2_higgs[4].Fill(h_M)
                higgs_5cut[4] += 1
              hist_l[4].Fill(M_jets.M())
              process[10] += 1
          #ggF Tight
            if ggF_Tight:
              if h_M:
                h_higgs[5].Fill(h_M)
                higgs_proc[5] += 1
              if (h_M > 120 and h_M < 130):
                h2_higgs[5].Fill(h_M)
                higgs_5cut[5] += 1
              hist_l[5].Fill(M_jets.M())
              process[11] += 1
          #Loose
          else:
            if h_M:
              h_higgs[6].Fill(h_M)
              higgs_proc[6] += 1
            if (h_M > 120 and h_M < 130):
              h2_higgs[6].Fill(h_M)
              higgs_5cut[6] += 1
            hist_l[6].Fill(M_jets.M())
            process[12] += 1
        ##########################

          
        ###########################         
        #hist.Fill(h_M)
        if abs(mu1.Eta()) < 0.8 and abs(mu2.Eta()) < 0.8 :
            h_barrel.Fill(h_M)
        if abs(mu1.Eta()) > 1.5 and abs(mu2.Eta()) > 1.5 :
            h_endcap.Fill(h_M)        
        hist_tot.Fill(h_M)

        npass +=1
        #print "nmuons",nmuons
        #print "njets",njets
        #print "npass",npass
    print process
    sum_iev += iev+1
    sum_tot[0] += nmuons
    sum_tot[1] += njets
    sum_tot[2] += npass
    sum_tot[3] += nhiggs_cut    
    sum_tot[4] += n2muons
    sum_tot[5] += nopp_charge
    sum_cut = []
    sum_cut += cut
    #hist.Draw()
    #out_root.Write()
    #out_root.Close()    

print sum_cut
print sum_tot
hist_tot.Draw()
hist_tot2.Draw()
#hist2.LabelsDeflate()
#hist2.LabelsOption("u","X")
hist2.Draw()
h_cutflow.Draw()

hist_l[0].Draw()
hist_l[1].Draw()
hist_l[2].Draw()
hist_l[3].Draw()
hist_l[4].Draw()
hist_l[5].Draw()
hist_l[6].Draw()

for i in range(3):
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
##############################################################
"""
h_eff2_pt.SetStats(0)
h_H2pt.Sumw2()
h_genpt.Sumw2()
h_eff2_pt.Divide(h_H2pt, h_genpt, 1.0, 1.0, "B")
h_eff2_pt.Draw()

h_eff2_eta.SetStats(0)
h_H2eta.Sumw2()
h_geneta.Sumw2()
h_eff2_eta.Divide(h_H2eta, h_geneta, 1.0, 1.0, "B")
h_eff2_eta.Draw()
"""
##############################################################
out_root_tot.Write()
out_root_tot.Close()

if not (num == 4):
    if num < 4 :
      Name = Phase_target[num-1]+"_"+filename
    if num > 4 :
      Name = Dataset_catTuple[sel*3-3+num-5]+"_"+filename
else:
    Name = "8TeV_data_set"+"_"+filename


out_num_event_tot.write("=" * 50)
out_num_event_tot.write("[[%s_no_ip]]" % Name)
out_num_event_tot.write("=" * 50)
#out_num_event_tot.write("\n # of events : %d\n total # of muons : %d\n total # of jets : %d\n total # of pass : %d\n total # of higgs_cut : %d\n total # of (muons > 1) : %d\n total # of opposite charge : %d\n" % (sum_iev,sum_tot[0],sum_tot[1],sum_tot[2],sum_tot[3],sum_tot[4],sum_tot[5]))
out_num_event_tot.write("\n # of events : %d\n total # of (muons > 1) : %d\n total # of opposite charge : %d\n total # of higgs_cut : %d\n " % (sum_iev,sum_tot[4],sum_tot[5],sum_tot[3]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n")
out_num_event_tot.write("catMuons label.\n # of muon candidates : %d\n muons # passed Tightmuon() : %d\n muons # passed pt_cut : %d\n muons # passed eta_cut : %d\n" % (cut[0], cut[1], cut[2], cut[3]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n")
out_num_event_tot.write("catJets label.\n # of jet candidates : %d\n jets # passed LooseId() : %d\n jets # passed pileupJetId() : %d\n jets # passed pt_cut : %d\n jets # passed eta_cut : %d\n" % (cut[4], cut[5], cut[6], cut[7], cut[8]))
#out_num_event_tot.write("in isTightMuon func. total tight-muon processes\n # of muon candidates : %d\n # of cut muon_pt : %d\n # of cut muon_eta : %d\n # of ParticleFlow : %d\n # of GlobalMuon : %d\n" % (sum_a,cut[0],cut[9],cut[1],cut[2]))
#out_num_event_tot.write(" # of global fit kai^2 : %d\n # of at least 1 valid muon hit associated with the global muon : %d\n # of at least chambers in different stations with matching segments : %d\n # of that track is null : %d\n" % (cut[3],cut[4],cut[5],cut[6]))
#out_num_event_tot.write(" # of pixel hits >=1 : %d\n # of tracker layers > 5 : %d\n" % (cut[7],cut[8]))
out_num_event_tot.write("=" * 50)
for i in range(13):
  out_num_event_tot.write("\n Jet-process%d : %d\n" % (i, process[i]))
out_num_event_tot.write("=" * 50)
for i in range(7):
  out_num_event_tot.write("\n higgs # : %d\n" % (higgs_proc[i]))
  out_num_event_tot.write(" higgs-cut # : %d\n" % (higgs_5cut[i]))
out_num_event_tot.write("=" * 50)
GC_Name = ["BB","BO","BE","OO","OE","EE"]
out_num_event_tot.write("\n '0'Jet Category - Tight\n")
for i in range(6):
  out_num_event_tot.write("%s count : %d\n" % (GC_Name[i], geo_Count[0][i]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n '0'Jet Category - Loose\n")
for i in range(6):
  out_num_event_tot.write("%s count : %d\n" % (GC_Name[i], geo_Count[1][i]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n '1'Jet Category - Tight\n")
for i in range(6):
  out_num_event_tot.write("%s count : %d\n" % (GC_Name[i], geo_Count[2][i]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n '1'Jet Category - Loose\n")
for i in range(6):
  out_num_event_tot.write("%s count : %d\n" % (GC_Name[i], geo_Count[3][i]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n")
out_num_event_tot.close()
