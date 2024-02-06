import numpy as np
import vector
import awkward as ak
import uproot
electron_mass=0.000511 #GeV
muon_mass=0.10566 #GeV

from Get_kinematics import Get_kinematics
from datasets import get_datasets

def lep_mass_pair_selection(event_number,lep_1,lep_2):
    if ((lep_1['pt'][event_number]!=0) and (lep_2['pt'][event_number])!=0):
            if abs(lep_1['Id'][event_number])==11:
                lep_1_mass=electron_mass
            else:
                lep_1_mass=muon_mass
            if abs(lep_2['Id'][event_number])==11:
                lep_2_mass=electron_mass
            else:
                lep_2_mass=muon_mass
            lep_1_vector=vector.obj(pt=lep_1['pt'][event_number],phi=lep_1['phi'][event_number],eta=lep_1['eta'][event_number],mass=lep_1_mass)
            lep_2_vector=vector.obj(pt=lep_2['pt'][event_number],phi=lep_2['phi'][event_number],eta=lep_2['eta'][event_number],mass=lep_2_mass)
            lep_Pair=lep_1_vector+lep_2_vector
            if ((lep_Pair.mass>12)):
                return True
            else:
                return False
    else:
        return False
    
def Match_A_Pair(event_number,Possible_Pair):
    if Possible_Pair[event_number]==1:
        return ['lep1','lep2']
    elif Possible_Pair[event_number]==2:
        return ['lep1','lep3']
    elif Possible_Pair[event_number]==3: 
        return ['lep2','lep3']
    elif Possible_Pair[event_number]==4: 
        return []
    elif Possible_Pair[event_number]==5:
        return ['lep1','lep4']
    elif Possible_Pair[event_number]==6:
        return ['lep2','lep4']
    elif Possible_Pair[event_number]==7:
        return ['lep3','lep4']
    else: return []
    
def Mass_selector(filepath,process):
    kinematics=Get_kinematics(filepath,process)[0]
    Pairs=Get_kinematics(filepath,process)[1]
    Possible_Pair_1=Pairs['Possible_pairs']['Possible_Pair_1']
    Possible_Pair_2=Pairs['Possible_pairs']['Possible_Pair_2']
    is_Selected=[]
    for i in range(len(Possible_Pair_1)):
        Matched_pair_1=Match_A_Pair(i,Possible_Pair_1)
        Matched_pair_2=Match_A_Pair(i,Possible_Pair_2)
        if (len(Matched_pair_1)==0 or len(Matched_pair_2)==0):
            is_Selected.append(False)
            continue
        if (lep_mass_pair_selection(i,kinematics[Matched_pair_1[0]],kinematics[Matched_pair_1[1]]) & lep_mass_pair_selection(i,kinematics[Matched_pair_2[0]],kinematics[Matched_pair_2[1]])):
             is_Selected.append(True)
        else:
             is_Selected.append(False)
    
    return is_Selected
             
             
        
def Store_Selected_Kinematics(filepath,process,fileLastFolder,fileName):
    kinematics=Get_kinematics(filepath,process)[0]
    Pairs=Get_kinematics(filepath,process)[1]
    top_kinematics = Get_kinematics(filepath,process)[2]
    is_Selected=Mass_selector(filepath,process)
    gen_lep_1_pt=kinematics['lep1']['pt'][is_Selected]
    gen_lep_1_phi=kinematics['lep1']['phi'][is_Selected]
    gen_lep_1_eta=kinematics['lep1']['eta'][is_Selected]
    gen_lep_1_pdgId=kinematics['lep1']['Id'][is_Selected]
    
    gen_lep_2_pt=kinematics['lep2']['pt'][is_Selected]
    gen_lep_2_phi=kinematics['lep2']['phi'][is_Selected]
    gen_lep_2_eta=kinematics['lep2']['eta'][is_Selected]
    gen_lep_2_pdgId=kinematics['lep2']['Id'][is_Selected]
    
    gen_lep_3_pt=kinematics['lep3']['pt'][is_Selected]
    gen_lep_3_phi=kinematics['lep3']['phi'][is_Selected]
    gen_lep_3_eta=kinematics['lep3']['eta'][is_Selected]
    gen_lep_3_pdgId=kinematics['lep3']['Id'][is_Selected]
    
    if process =='tttt':
        gen_lep_4_pt=kinematics['lep4']['pt'][is_Selected]
        gen_lep_4_phi=kinematics['lep4']['phi'][is_Selected]
        gen_lep_4_eta=kinematics['lep4']['eta'][is_Selected]
        gen_lep_4_pdgId=kinematics['lep4']['Id'][is_Selected]
        
    #Top Kinematics
    gen_top_1_pt=top_kinematics['top1']['pt'][is_Selected]
    gen_top_1_phi=top_kinematics['top1']['phi'][is_Selected]
    gen_top_1_eta=top_kinematics['top1']['eta'][is_Selected]
    
    gen_top_2_pt=top_kinematics['top2']['pt'][is_Selected]
    gen_top_2_phi=top_kinematics['top2']['phi'][is_Selected]
    gen_top_2_eta=top_kinematics['top2']['eta'][is_Selected]
    
    gen_top_3_pt=top_kinematics['top3']['pt'][is_Selected]
    gen_top_3_phi=top_kinematics['top3']['phi'][is_Selected]
    gen_top_3_eta=top_kinematics['top3']['eta'][is_Selected]
    
    if process=='tttt':
        gen_top_4_pt=top_kinematics['top4']['pt'][is_Selected]
        gen_top_4_phi=top_kinematics['top4']['phi'][is_Selected]
        gen_top_4_eta=top_kinematics['top4']['eta'][is_Selected]
        
    #Pairs
    Possible_Pair_1=Pairs['Possible_pairs']['Possible_Pair_1'][is_Selected]
    Possible_Pair_2=Pairs['Possible_pairs']['Possible_Pair_2'][is_Selected]
    
    if process=='tttt':
        with uproot.recreate("selected_Events/"+str(process)+"/Mass_Pair_requirement/"+str(fileLastFolder)+"/"+str(fileName)) as output_file:
            output_file["Events"] = {
                            "gen_lep_1_pt": gen_lep_1_pt, "gen_lep_1_phi": gen_lep_1_phi,"gen_lep_1_eta": gen_lep_1_eta, "gen_lep_1_pdgId": gen_lep_1_pdgId,
                            "gen_lep_2_pt": gen_lep_2_pt, "gen_lep_2_phi": gen_lep_2_phi,"gen_lep_2_eta": gen_lep_2_eta, "gen_lep_2_pdgId": gen_lep_2_pdgId,
                            "gen_lep_3_pt": gen_lep_3_pt, "gen_lep_3_phi": gen_lep_3_phi,"gen_lep_3_eta": gen_lep_3_eta, "gen_lep_3_pdgId": gen_lep_3_pdgId,
                            "gen_lep_4_pt": gen_lep_4_pt, "gen_lep_4_phi": gen_lep_4_phi,"gen_lep_4_eta": gen_lep_4_eta, "gen_lep_4_pdgId": gen_lep_4_pdgId,
                            "gen_top_1_pt": gen_top_1_pt, "gen_top_1_phi": gen_top_1_phi,"gen_top_1_eta": gen_top_1_eta, 
                            "gen_top_2_pt": gen_top_2_pt, "gen_top_2_phi": gen_top_2_phi,"gen_top_2_eta": gen_top_2_eta, 
                            "gen_top_3_pt": gen_top_3_pt, "gen_top_3_phi": gen_top_3_phi,"gen_top_3_eta": gen_top_3_eta,
                            "gen_top_4_pt": gen_top_4_pt, "gen_top_4_phi": gen_top_4_phi,"gen_top_4_eta": gen_top_4_eta,  
                            "Possible_Pair_1":Possible_Pair_1,"Possible_Pair_2":Possible_Pair_2,
                        }
    else:   
        with uproot.recreate("selected_Events/"+str(process)+"/Mass_Pair_requirement/"+str(fileLastFolder)+"/"+str(fileName)) as output_file:
            output_file["Events"] = {
                            "gen_lep_1_pt": gen_lep_1_pt, "gen_lep_1_phi": gen_lep_1_phi,"gen_lep_1_eta": gen_lep_1_eta, "gen_lep_1_pdgId": gen_lep_1_pdgId,
                            "gen_lep_2_pt": gen_lep_2_pt, "gen_lep_2_phi": gen_lep_2_phi,"gen_lep_2_eta": gen_lep_2_eta, "gen_lep_2_pdgId": gen_lep_2_pdgId,
                            "gen_lep_3_pt": gen_lep_3_pt, "gen_lep_3_phi": gen_lep_3_phi,"gen_lep_3_eta": gen_lep_3_eta, "gen_lep_3_pdgId": gen_lep_3_pdgId,
                            "gen_top_1_pt": gen_top_1_pt, "gen_top_1_phi": gen_top_1_phi,"gen_top_1_eta": gen_top_1_eta, 
                            "gen_top_2_pt": gen_top_2_pt, "gen_top_2_phi": gen_top_2_phi,"gen_top_2_eta": gen_top_2_eta, 
                            "gen_top_3_pt": gen_top_3_pt, "gen_top_3_phi": gen_top_3_phi,"gen_top_3_eta": gen_top_3_eta, 
                            "Possible_Pair_1":Possible_Pair_1,"Possible_Pair_2":Possible_Pair_2,
                        }
    output_file.close()
        
datasets=get_datasets()

for process,prcs in datasets.items():
        for i in range(prcs['nrepos']):
            repo=(prcs['repos']).split(',')[i]
            Nfiles=len((prcs[repo]).split(','))
            for f in range(Nfiles):
                file=(prcs[repo].split(','))[f] 
                filepath='selected_Events/'+str(process)+'/'+str(repo)+'/'+str(file)
                Store_Selected_Kinematics(filepath,process,repo,file)
