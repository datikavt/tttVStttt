import uproot
import pandas as pd
import vector
import numpy as np
import awkward as ak
vector.register_awkward()
import os  # file handling
import ROOT
import argparse  # I/O
import energyflow 
import mplhep as hep #CAT
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
hep.style.use("CMS") 
import cmsstyle as CMS # For ROOT plotting

parser = argparse.ArgumentParser(description='selecting events from a ROOT file.')
parser.add_argument('-f','--file',type=str, help='specify the filepath that you want to start pairing from')
parser.add_argument('-p','--process',type=str,help='specify the process ex. tttW, tttt, tttJ',default='tttW')
args = parser.parse_args()

fileData=uproot.open(args.file)
events=fileData["Events"]
Same_Mother_Top=np.array(events["Same_Mother_top"])
is_7=ak.fill_none((Same_Mother_Top==7),False)
Same_Mother_Top_is_7=Same_Mother_Top[is_7]
print(len(Same_Mother_Top_is_7)) 


top_mass=172.76 #GeV
electron_mass=0.000511 #GeV
muon_mass=0.10566 #GeV
top_masses=np.ones(len(np.array(events["gen_top_1_eta"])))*top_mass


#top kinematics
top1_etas=np.array(events["gen_top_1_eta"])
top1_phis=np.array(events["gen_top_1_phi"])
top1_pts=np.array(events["gen_top_1_pt"])
top1_ys=energyflow.ys_from_pts_etas_ms(top1_pts,top1_etas, top_masses)

top2_etas=np.array(events["gen_top_2_eta"])
top2_phis=np.array(events["gen_top_2_phi"])
top2_pts=np.array(events["gen_top_2_pt"])
top2_ys=energyflow.ys_from_pts_etas_ms(top2_pts,top2_etas, top_masses)

top3_etas=np.array(events["gen_top_3_eta"])
top3_phis=np.array(events["gen_top_3_phi"])
top3_pts=np.array(events["gen_top_3_pt"])
top3_ys=energyflow.ys_from_pts_etas_ms(top3_pts,top3_etas, top_masses)

if args.process=='tttt':
    top4_etas=np.array(events["gen_top_4_eta"])
    top4_phis=np.array(events["gen_top_4_phi"])
    top4_pts=np.array(events["gen_top_4_pt"])
    top4_ys=energyflow.ys_from_pts_etas_ms(top4_pts,top4_etas, top_masses)

#lepton kinematics + pdgIds
lepton1_pts=np.where(np.array(events["gen_lep_1_pt"]),np.array(events["gen_lep_1_pt"]),0)
lepton1_etas=np.where(np.array(events["gen_lep_1_eta"]),np.array(events["gen_lep_1_eta"]),0) #Where we have lepton1_etass keep it otherwise set it to 0
lepton1_phis=np.where(np.array(events["gen_lep_1_phi"]),np.array(events["gen_lep_1_phi"]),0) 
lepton1_Ids=np.where(np.array(events["gen_lep_1_pdgId"]),np.array(events["gen_lep_1_pdgId"]),0)
lepton1_pts=lepton1_pts.astype('float') 
lepton1_etas=lepton1_etas.astype('float')
lepton1_phis=lepton1_phis.astype('float')
lepton1_Ids=lepton1_Ids.astype('float')
lepton1_ms = np.where(lepton1_Ids == 11, electron_mass, muon_mass) #Notice that if the particle is missing we will give it muon mass we can use np.select if we want 0s
lepton1_ys=energyflow.ys_from_pts_etas_ms(lepton1_pts, lepton1_etas, lepton1_ms)

lepton2_pts=np.where(np.array(events["gen_lep_2_pt"]),np.array(events["gen_lep_2_pt"]),0)
lepton2_etas=np.where(np.array(events["gen_lep_2_eta"]),np.array(events["gen_lep_2_eta"]),0)
lepton2_phis=np.where(np.array(events["gen_lep_2_phi"]),np.array(events["gen_lep_2_phi"]),0) 
lepton2_Ids=np.where(np.array(events["gen_lep_2_pdgId"]),np.array(events["gen_lep_2_pdgId"]),0)
lepton2_pts = lepton2_pts.astype('float')
lepton2_etas= lepton2_etas.astype('float')
lepton2_phis=lepton2_phis.astype('float')
lepton2_Ids=lepton2_Ids.astype('float')
lepton2_ms = np.where(lepton2_Ids == 11, electron_mass, muon_mass)  
lepton2_ys=energyflow.ys_from_pts_etas_ms(lepton2_pts, lepton2_etas, lepton2_ms)


lepton3_pts=np.where(np.array(events["gen_lep_3_pt"]),np.array(events["gen_lep_3_pt"]),0)
lepton3_etas=np.where(np.array(events["gen_lep_3_eta"]),np.array(events["gen_lep_3_eta"]),0)
lepton3_phis=np.where(np.array(events["gen_lep_3_phi"]),np.array(events["gen_lep_3_phi"]),0) 
lepton3_Ids=np.where(np.array(events["gen_lep_3_pdgId"]),np.array(events["gen_lep_3_pdgId"]),0)
lepton3_pts = lepton3_pts.astype('float')
lepton3_etas= lepton3_etas.astype('float')
lepton3_phis=lepton3_phis.astype('float')
lepton3_Ids=lepton3_Ids.astype('float')
lepton3_ms = np.where(lepton3_Ids == 11, electron_mass, muon_mass)  
lepton3_ys=energyflow.ys_from_pts_etas_ms(lepton3_pts, lepton3_etas, lepton3_ms)


if args.process=='tttt':
    lepton4_pts=np.where(np.array(events["gen_lep_4_pt"]),np.array(events["gen_lep_4_pt"]),0)
    lepton4_etas=np.where(np.array(events["gen_lep_4_eta"]),np.array(events["gen_lep_4_eta"]),0)
    lepton4_phis=np.where(np.array(events["gen_lep_4_phi"]),np.array(events["gen_lep_4_phi"]),0) 
    lepton4_Ids=np.where(np.array(events["gen_lep_4_pdgId"]),np.array(events["gen_lep_4_pdgId"]),0)
    lepton4_pts = lepton4_pts.astype('float')
    lepton4_etas= lepton4_etas.astype('float')
    lepton4_phis=lepton4_phis.astype('float')
    lepton4_Ids=lepton4_Ids.astype('float')
    lepton4_ms = np.where(lepton4_Ids == 11, electron_mass, muon_mass)  
    lepton4_ys=energyflow.ys_from_pts_etas_ms(lepton4_pts, lepton4_etas, lepton4_ms)


top1_ptyphims=np.array([top1_pts,top1_ys,top1_phis,top_masses])
top2_ptyphims=np.array([top2_pts,top2_ys,top2_phis,top_masses])
top3_ptyphims=np.array([top3_pts,top3_ys,top3_phis,top_masses])    
lepton1_ptyphipms=np.array([lepton1_pts,lepton1_ys,lepton1_phis,lepton1_ms])
lepton2_ptyphipms=np.array([lepton2_pts,lepton2_ys,lepton2_phis,lepton2_ms])
lepton3_ptyphipms=np.array([lepton3_pts,lepton3_ys,lepton3_phis,lepton3_ms])
if args.process=='tttt':
    top4_ptyphims=np.array([top4_pts,top4_ys,top4_phis,top_masses])
    lepton4_ptyphipms=np.array([lepton4_pts,lepton4_ys,lepton4_phis,lepton4_ms])

# Convert hadronic coordinates to cartesian
top1_p4s=energyflow.p4s_from_ptyphims(np.transpose(top1_ptyphims)) #Notice that we need to transpose! #it now contains
                                                                                                         #[[E, px, py, pz]
#                                                                                                          [E, px, py, pz]
#                                                                                                           .
#                                                                                                           .
#                                                                                                           .
#                                                                                                          [E, px, py, pz]]
top2_p4s=energyflow.p4s_from_ptyphims(np.transpose(top2_ptyphims))
top3_p4s=energyflow.p4s_from_ptyphims(np.transpose(top3_ptyphims))
lepton1_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton1_ptyphipms))
lepton2_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton2_ptyphipms))
lepton3_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton3_ptyphipms))
# construct vectors in p4s
top1 = vector.obj(px=top1_p4s[:,1], py=top1_p4s[:,2], pz=top1_p4s[:,3], mass=top_mass)
top2 = vector.obj(px=top2_p4s[:,1], py=top2_p4s[:,2], pz=top2_p4s[:,3], mass=top_mass)
top3 = vector.obj(px=top3_p4s[:,1], py=top3_p4s[:,2], pz=top3_p4s[:,3], mass=top_mass)

lep1= vector.obj(px=lepton1_p4s[:,1], py=lepton1_p4s[:,2], pz=lepton1_p4s[:,3], mass=lepton1_ms)
lep2= vector.obj(px=lepton2_p4s[:,1], py=lepton2_p4s[:,2], pz=lepton2_p4s[:,3], mass=lepton2_ms)
lep3= vector.obj(px=lepton3_p4s[:,1], py=lepton3_p4s[:,2], pz=lepton3_p4s[:,3], mass=lepton3_ms)
#For boosting with vector package check: https://github.com/scikit-hep/vector/blob/main/src/vector/_methods.py
#Check test.py if you want to see more details on the numbers
#Boosts the lepton to top1 rest frame
lep1_in_top_frame=lep1.boostCM_of_p4(top1) #Remember these are all in p4s
lep2_in_top_frame=lep2.boostCM_of_p4(top2) 
lep3_in_top_frame=lep3.boostCM_of_p4(top3)
if args.process=='tttt':
    top4_p4s=energyflow.p4s_from_ptyphims(np.transpose(top4_ptyphims))
    lepton4_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton4_ptyphipms))
    top4 = vector.obj(px=top4_p4s[:,1], py=top4_p4s[:,2], pz=top4_p4s[:,3], mass=top_mass)
    lep4 = vector.obj(px=lepton4_p4s[:,1], py=lepton4_p4s[:,2], pz=lepton4_p4s[:,3], mass=lepton4_ms)
    lep4_in_top_frame=lep4.boostCM_of_p4(top4)


#Calculates the cos_phi for the indicated pair of leptons and for the indicated_event
def Calculate_cos_phi(event_number):
    if Same_Mother_Top[event_number]==1:
        #Dot product of the vectors
        lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
        lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
        cos_phi_e=(lep1_in_top_frame_3D @ lep2_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep2_in_top_frame_3D.mag)
        return cos_phi_e
    if Same_Mother_Top[event_number]==2:
        lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
        lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
        cos_phi_e=(lep1_in_top_frame_3D @ lep3_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep3_in_top_frame_3D.mag)
        return cos_phi_e
    if Same_Mother_Top[event_number]==3:
        lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
        lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
        cos_phi_e=(lep2_in_top_frame_3D @ lep3_in_top_frame_3D)/(lep2_in_top_frame_3D.mag * lep3_in_top_frame_3D.mag)
        return cos_phi_e
    if Same_Mother_Top[event_number]==4 or Same_Mother_Top[event_number]==0: 
        cos_phi_e=0 # Put 0s for where the Mother is not clear
        return cos_phi_e
    if Same_Mother_Top[event_number]==5:
        lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
        lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
        cos_phi_e=(lep1_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
        return cos_phi_e
    if Same_Mother_Top[event_number]==6:
        lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
        lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
        cos_phi_e=(lep2_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep2_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
        return cos_phi_e
    if Same_Mother_Top[event_number]==7:
        lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
        lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
        cos_phi_e=(lep3_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep3_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
        return cos_phi_e
    
    
bins=np.linspace(-1,1,9)
cos_phi_with_0s=[] #Will store cos_phis of every event and 0s for events without tt pair 
cos_phi_hist_r=ROOT.TH1F("cos_phi_histogram","cos_phi_histogram",8,-1,1) #define the cos_phi_hist
for i in range(len(Same_Mother_Top)):
    cos_phi_event=Calculate_cos_phi(i)
    cos_phi_with_0s.append(cos_phi_event)
    if Same_Mother_Top[i]==4 or Same_Mother_Top[i]==0: 
        continue
    else:
        cos_phi_hist_r.Fill(cos_phi_event) 
    
Good_Mother=ak.fill_none(((Same_Mother_Top!=0 ) & (Same_Mother_Top!=4)),False)
cos_phi_with_0s=np.array(cos_phi_with_0s)
cos_phi_without_0s=cos_phi_with_0s[Good_Mother]
print(len(cos_phi_with_0s))
print(len(cos_phi_without_0s))
midbins=(bins[:-1] + bins[1:])/2
bin_width=bins[1:]-bins[:-1]
cos_phi_hist_with_0s=np.histogram(np.array(cos_phi_with_0s),bins)
cos_phi_hist_without_0s=np.histogram(np.array(cos_phi_without_0s),bins)
tttW_style = {'histtype': 'step', 'density': False, 'lw': 1, 'zorder': 2,'label': 'tttW'}
#From CAT
fig, ax = plt.subplots()
hep.cms.label("", ax=ax, loc=0)

ax.hist(cos_phi_with_0s,8)
# ax.hist(cos_phi_with_0s,midbins,*tttW_style)
# Calculate D with event's that introduce 0 s 
D_with_0s=-3*np.average(np.asarray(cos_phi_with_0s))
D_without_0s=-3*np.average(np.asarray(cos_phi_without_0s))
print("D_with_0s is " + str(D_with_0s))
print("D without 0s is " + str(D_without_0s))
fig.savefig('Cos_phi_with_0_'+str(args.process) +'.pdf', bbox_inches='tight')


#ROOT plotting
CMS.SetExtraText("Simulation Preliminary")
CMS.SetLumi("")
canv = CMS.cmsCanvas('', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)
hdf = CMS.GetcmsCanvasHist(canv)
hdf.GetYaxis().SetMaxDigits(2)
# Shift multiplier position
ROOT.TGaxis.SetExponentOffset(-0.10, 0.01, "Y")
cos_phi_hist_r.Draw("HIST")
canv.SaveAs(os.path.join("plots","ROOT_cos_phi:"+str(args.process)+";n_toys"+".pdf"))
