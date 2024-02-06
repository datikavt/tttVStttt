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
import math
top_mass=172.76 #GeV
electron_mass=0.000511 #GeV
muon_mass=0.10566 #GeV



def Get_cosphi(file,process):
    fileData=uproot.open(file)
    events=fileData["Events"]
    Possible_Pair_1=np.array(events["Possible_Pair_1"])
    Possible_Pair_2=np.array(events["Possible_Pair_2"])


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

    if process=='tttt':
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


    if  process=='tttt':
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
    if process=='tttt':
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
    if process=='tttt':
        top4_p4s=energyflow.p4s_from_ptyphims(np.transpose(top4_ptyphims))
        lepton4_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton4_ptyphipms))
        top4 = vector.obj(px=top4_p4s[:,1], py=top4_p4s[:,2], pz=top4_p4s[:,3], mass=top_mass)
        lep4 = vector.obj(px=lepton4_p4s[:,1], py=lepton4_p4s[:,2], pz=lepton4_p4s[:,3], mass=lepton4_ms)
        lep4_in_top_frame=lep4.boostCM_of_p4(top4)


    #Calculates the cos_phi for the indicated pair of leptons and for the indicated_event
    def Calculate_cos_phi(event_number,Pair):
        if Pair[event_number]==1:
            #Dot product of the vectors
            lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
            lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
            cos_phi_e=(lep1_in_top_frame_3D @ lep2_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep2_in_top_frame_3D.mag)
            return cos_phi_e
        if Pair[event_number]==2:
            lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
            lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
            cos_phi_e=(lep1_in_top_frame_3D @ lep3_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep3_in_top_frame_3D.mag)
            return cos_phi_e
        if Pair[event_number]==3:
            lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
            lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
            cos_phi_e=(lep2_in_top_frame_3D @ lep3_in_top_frame_3D)/(lep2_in_top_frame_3D.mag * lep3_in_top_frame_3D.mag)
            return cos_phi_e
        if Pair[event_number]==4 or Pair[event_number]==0:  # Put 0s for where the Mother is not clear
            return "continue"
        if Pair[event_number]==5:
            lep1_in_top_frame_3D=vector.obj(px=lep1_in_top_frame.px[event_number],py=lep1_in_top_frame.py[event_number],pz=lep1_in_top_frame.pz[event_number])
            lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
            cos_phi_e=(lep1_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep1_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
            return cos_phi_e
        if Pair[event_number]==6:
            lep2_in_top_frame_3D=vector.obj(px=lep2_in_top_frame.px[event_number],py=lep2_in_top_frame.py[event_number],pz=lep2_in_top_frame.pz[event_number])
            lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
            cos_phi_e=(lep2_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep2_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
            return cos_phi_e
        if Pair[event_number]==7:
            lep3_in_top_frame_3D=vector.obj(px=lep3_in_top_frame.px[event_number],py=lep3_in_top_frame.py[event_number],pz=lep3_in_top_frame.pz[event_number])
            lep4_in_top_frame_3D=vector.obj(px=lep4_in_top_frame.px[event_number],py=lep4_in_top_frame.py[event_number],pz=lep4_in_top_frame.pz[event_number])
            cos_phi_e=(lep3_in_top_frame_3D @ lep4_in_top_frame_3D)/(lep3_in_top_frame_3D.mag * lep4_in_top_frame_3D.mag)
            return cos_phi_e

    cos_phi=[] #Will store cos_phis of every possible ttbar pair

    for i in range(len(Possible_Pair_1)):
        cos_phi_pair=Calculate_cos_phi(i,Possible_Pair_1)
        if cos_phi_pair != "continue":
            cos_phi.append(cos_phi_pair)
        else:
            continue    
        cos_phi_pair=Calculate_cos_phi(i,Possible_Pair_2)
        cos_phi.append(cos_phi_pair)
    
    return cos_phi
    


def Get_cosphi_lab(file,process):
    
    fileData=uproot.open(file)
    events=fileData["Events"]
    Possible_Pair_1=np.array(events["Possible_Pair_1"])
    Possible_Pair_2=np.array(events["Possible_Pair_2"])


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

    if process=='tttt':
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


    if  process=='tttt':
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
    if process=='tttt':
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

    if process=='tttt':
        top4_p4s=energyflow.p4s_from_ptyphims(np.transpose(top4_ptyphims))
        lepton4_p4s=energyflow.p4s_from_ptyphims(np.transpose(lepton4_ptyphipms))
        top4 = vector.obj(px=top4_p4s[:,1], py=top4_p4s[:,2], pz=top4_p4s[:,3], mass=top_mass)
        lep4 = vector.obj(px=lepton4_p4s[:,1], py=lepton4_p4s[:,2], pz=lepton4_p4s[:,3], mass=lepton4_ms)

    def Calculate_cos_phi_lab(event_number,Pair):
        if Pair[event_number]==1:
            #Dot product of the vectors
            lep1_3D=vector.obj(px=lep1.px[event_number],py=lep1.py[event_number],pz=lep1.pz[event_number])
            lep2_3D=vector.obj(px=lep2.px[event_number],py=lep2.py[event_number],pz=lep2.pz[event_number])
            cos_phi_e=(lep1_3D @ lep2_3D)/(lep1_3D.mag * lep2_3D.mag)
            return cos_phi_e
        if Pair[event_number]==2:
            lep1_3D=vector.obj(px=lep1.px[event_number],py=lep1.py[event_number],pz=lep1.pz[event_number])
            lep3_3D=vector.obj(px=lep3.px[event_number],py=lep3.py[event_number],pz=lep3.pz[event_number])
            cos_phi_e=(lep1_3D @ lep3_3D)/(lep1_3D.mag * lep3_3D.mag)
            return cos_phi_e
        if Pair[event_number]==3:
            lep2_3D=vector.obj(px=lep2.px[event_number],py=lep2.py[event_number],pz=lep2.pz[event_number])
            lep3_3D=vector.obj(px=lep3.px[event_number],py=lep3.py[event_number],pz=lep3.pz[event_number])
            cos_phi_e=(lep2_3D @ lep3_3D)/(lep2_3D.mag * lep3_3D.mag)
            return cos_phi_e
        if Pair[event_number]==4 or Pair[event_number]==0: 
            cos_phi_e=0 # Put 0s for where the Mother is not clear
            return cos_phi_e
        if Pair[event_number]==5:
            lep1_3D=vector.obj(px=lep1.px[event_number],py=lep1.py[event_number],pz=lep1.pz[event_number])
            lep4_3D=vector.obj(px=lep4.px[event_number],py=lep4.py[event_number],pz=lep4.pz[event_number])
            cos_phi_e=(lep1_3D @ lep4_3D)/(lep1_3D.mag * lep4_3D.mag)
            return cos_phi_e
        if Pair[event_number]==6:
            lep2_3D=vector.obj(px=lep2.px[event_number],py=lep2.py[event_number],pz=lep2.pz[event_number])
            lep4_3D=vector.obj(px=lep4.px[event_number],py=lep4.py[event_number],pz=lep4.pz[event_number])
            cos_phi_e=(lep2_3D @ lep4_3D)/(lep2_3D.mag * lep4_3D.mag)
            return cos_phi_e
        if Pair[event_number]==7:
            lep3_3D=vector.obj(px=lep3.px[event_number],py=lep3.py[event_number],pz=lep3.pz[event_number])
            lep4_3D=vector.obj(px=lep4.px[event_number],py=lep4.py[event_number],pz=lep4.pz[event_number])
            cos_phi_e=(lep3_3D @ lep4_3D)/(lep3_3D.mag * lep4_3D.mag)
            return cos_phi_e
    
    cos_phi_lab=[] #Will store cos_phis_labs of every possible ttbar pair

    for i in range(len(Possible_Pair_1)):
        cos_phi_pair=Calculate_cos_phi_lab(i,Possible_Pair_1)
        cos_phi_lab.append(cos_phi_pair)
        cos_phi_pair=Calculate_cos_phi_lab(i,Possible_Pair_2)
        cos_phi_lab.append(cos_phi_pair)

    return cos_phi_lab


def normalize_angle(angle):
    # Normalize angle to be between 0 and 
    normalized_angle = angle % (2*math.pi)
    # Make sure the angle is between 0 and 1 radian
    if normalized_angle > math.pi:
        normalized_angle = 2*math.pi - normalized_angle
    return normalized_angle


def Get_delta_phi(file,process):
    
    fileData=uproot.open(file)
    events=fileData["Events"]
    
    Possible_Pair_1=np.array(events["Possible_Pair_1"])
    Possible_Pair_2=np.array(events["Possible_Pair_2"])
    
    lepton1_phis=np.where(np.array(events["gen_lep_1_phi"]),np.array(events["gen_lep_1_phi"]),0) 
    lepton2_phis=np.where(np.array(events["gen_lep_2_phi"]),np.array(events["gen_lep_2_phi"]),0)
    lepton3_phis=np.where(np.array(events["gen_lep_3_phi"]),np.array(events["gen_lep_3_phi"]),0)
    lepton1_phis=lepton1_phis.astype('float')
    lepton2_phis=lepton2_phis.astype('float') 
    lepton3_phis=lepton3_phis.astype('float')
    if process=='tttt':
        lepton4_phis=np.where(np.array(events["gen_lep_4_phi"]),np.array(events["gen_lep_4_phi"]),0) 
        lepton4_phis=lepton3_phis.astype('float')
    
    def Get_delta_phi_pair(event_number,Pair):
        if Pair[event_number]==1:
            phi1=(lepton1_phis[event_number]) % (2*math.pi)
            phi2=(lepton2_phis[event_number]) % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi1-phi2))
            return delta_phi_pair
        if Pair[event_number]==2:
            phi1=(lepton1_phis[event_number])  % (2*math.pi)
            phi3=(lepton3_phis[event_number])  % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi1-phi3))
            return delta_phi_pair
        if Pair[event_number]==3:
            phi2=(lepton2_phis[event_number]) % (2*math.pi)
            phi3=(lepton3_phis[event_number]) % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi2-phi3))
            return delta_phi_pair
        if Pair[event_number]==4 or Pair[event_number]==0: 
            return False
        if Pair[event_number]==5:
            phi1=(lepton1_phis[event_number])  % (2*math.pi)
            phi4=(lepton4_phis[event_number])  % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi1-phi4))
            return delta_phi_pair
        if Pair[event_number]==6:
            phi2=(lepton2_phis[event_number])  % (2*math.pi)
            phi4=(lepton4_phis[event_number])  % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi2-phi4))
            return delta_phi_pair
        if Pair[event_number]==7:
            phi3=(lepton3_phis[event_number])  % (2*math.pi)
            phi4=(lepton4_phis[event_number])  % (2*math.pi)
            delta_phi_pair=normalize_angle(abs(phi3-phi4))
            return delta_phi_pair
    
    delta_phi=[]
    for i in range(len(Possible_Pair_1)):
        delta_phi_pair=Get_delta_phi_pair(i,Possible_Pair_1)
        if delta_phi_pair!=False:
            delta_phi.append(delta_phi_pair)
        delta_phi_pair=Get_delta_phi_pair(i,Possible_Pair_2)
        if delta_phi_pair!=False: #not required but just to be sure
            delta_phi.append(delta_phi_pair)
        
    return delta_phi


def Calculate_stat_err(cosphi):
    ntoys=1000
    fluctuated_cosphi_w=np.random.poisson(1,size=(len(cosphi),ntoys)) 
    av_fluctuated_D=0
    fluctuated_D_proxies=[]
    for i in range(ntoys):
        fluctuated_cosphi=cosphi*fluctuated_cosphi_w[:,i] 
        fluctuated_D_proxy=-3*np.average(np.asarray(fluctuated_cosphi))
        fluctuated_D_proxies.append(fluctuated_D_proxy)
        av_fluctuated_D+=fluctuated_D_proxy
    
    av_fluctuated_D=av_fluctuated_D/ntoys
    var=0
    for i in range(ntoys):
        var+=(fluctuated_D_proxies[i]-av_fluctuated_D)*(fluctuated_D_proxies-av_fluctuated_D)
    
    std_err=np.sqrt(var/ntoys)
    
    return std_err
        