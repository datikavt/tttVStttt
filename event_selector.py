import uproot
import pandas as pd
import vector
import numpy as np
import awkward as ak
vector.register_awkward()
import os  # file handling
import argparse  # I/O
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


parser = argparse.ArgumentParser(description='selecting events from a ROOT file.')
parser.add_argument('-f','--file',type=str, help='specify the filepath that you want to select events from')
parser.add_argument('-p','--process',type=str,help='specify the process ex. tttW, tttt, tttJ')
args = parser.parse_args()

print('INFO: selecting events from the file:' + str(args.file))

fileData=uproot.open(args.file)

events = NanoEventsFactory.from_root(
    args.file,
    schemaclass=NanoAODSchema,
).events()


# ------------------------------------------------------------
# use `GenParticle` behavior
gen_part = ak.with_name(events["GenPart"], "GenParticle") 

# derived gen-particle properties
gen_part["genPartIdxMotherMasked"] = ak.mask(gen_part.genPartIdxMother, gen_part.genPartIdxMother >= 0) #Containing only elements where the gen_part has a mother
gen_part["absPdgId"] = abs(gen_part.pdgId)
gen_part["index"] = ak.local_index(gen_part, axis=1) #Index the gen_part i.e. to know the index of it in an event.
                                                    #  for example [0.5,200,21],[0.4,300,210] after indexing [0,1,2],[0,1,2]

def parent(gen_part):
        """
        Return an array with the same structure as `gen_part` that contains the parent particle
        at the same positions. Entries are masked if not parent particle exists.
        """
        return gen_part[gen_part.genPartIdxMotherMasked]

def is_descended_from(gen_part, idx):
        """
        Return an array with the same structure as `gen_part` that indicated whether the
        respective generator particle is part of the decay chain of the particle at position `idx`.
        """
        idx = ak.fill_none(idx, -1)
        result = ak.zeros_like(gen_part.index, dtype=bool)

        # apply `parents` repeatedly and check index
        gen_parent = gen_part
        while ak.any(~ak.is_none(gen_parent.index, axis=1)):
            result = result | ak.fill_none(gen_parent.index == idx, False)
            gen_parent = gen_part[gen_parent.genPartIdxMotherMasked]

        return result

# extract bool flags
is_hard_proc = gen_part.hasFlags("isHardProcess")
from_hard_proc = gen_part.hasFlags("fromHardProcess")
is_last_copy = gen_part.hasFlags("isLastCopy") # the last copy of the particle in the chain with the same pdg id
                                             # (and therefore is more likely, but not guaranteed,
                                             # to carry the final physical momentum)

# find top quarks produced during the hard scattering
gen_part_hp = gen_part[is_hard_proc]
gen_hp_top = gen_part_hp[gen_part_hp.absPdgId == 6]
gen_hp_top = ak.pad_none(gen_hp_top, 4)
gen_hp_top_1 = gen_hp_top[:, 0]
gen_hp_top_2 = gen_hp_top[:, 1]
gen_hp_top_3=gen_hp_top[:,2]
if args.process=='tttt':
    gen_hp_top_4=gen_hp_top[:,3]

# find last copy of top quarks before decay. 
gen_top = gen_part[(gen_part.absPdgId == 6) & is_last_copy]
gen_top = ak.pad_none(gen_top, 4)
gen_top_1 = gen_top[:, 0]
gen_top_2 = gen_top[:, 1]
gen_top_3 = gen_top[:,2]
if args.process=='tttt':
    gen_top_4 = gen_top[:,3]

# identify top quark decay products
is_top_1_decay = is_descended_from(gen_part, gen_top_1.index)
is_top_2_decay = is_descended_from(gen_part, gen_top_2.index)
is_top_3_decay= is_descended_from(gen_part, gen_top_3.index)
gen_part["is_top_decay"] = (is_top_1_decay | is_top_2_decay | is_top_3_decay)
if args.process=='tttt':
    is_top_4_decay= is_descended_from(gen_part, gen_top_4.index)
    gen_part["is_top_decay"] = (is_top_1_decay | is_top_2_decay | is_top_3_decay | is_top_4_decay)


# identify W bosons coming from top decay!
is_w = (gen_part.absPdgId == 24) & from_hard_proc & is_last_copy
gen_w_1 = ak.firsts(gen_part[is_w & is_top_1_decay])
gen_w_2 = ak.firsts(gen_part[is_w & is_top_2_decay])
gen_w_3= ak.firsts(gen_part[is_w & is_top_3_decay])
if args.process=='tttt':
    gen_w_4= ak.firsts(gen_part[is_w & is_top_4_decay])

# identify W boson decay products
is_w_1_decay = is_descended_from(gen_part, gen_w_1.index)
is_w_2_decay = is_descended_from(gen_part, gen_w_2.index)
is_w_3_decay = is_descended_from(gen_part, gen_w_3.index)
gen_part["is_w_decay"] = (is_w_1_decay | is_w_2_decay | is_w_3_decay)
if args.process=='tttt':
    is_w_4_decay = is_descended_from(gen_part, gen_w_4.index)
    gen_part["is_w_decay"] = (is_w_1_decay | is_w_2_decay | is_w_3_decay | is_w_4_decay )

# charged leptons (e, mu only)
is_lepton = is_last_copy & from_hard_proc & ((gen_part.absPdgId == 11) | (gen_part.absPdgId == 13))
has_enough_pt=ak.fill_none((gen_part.pt>5),False)
lepton_condition1=ak.fill_none((is_lepton & is_w_1_decay & has_enough_pt),False)
lepton_condition2=ak.fill_none((is_lepton & is_w_2_decay & has_enough_pt),False)
lepton_condition3=ak.fill_none((is_lepton & is_w_3_decay & has_enough_pt),False)

gen_lep_1 = gen_part[lepton_condition1] #We require the pt of the lepton to be at least 5 GeV
gen_lep_2 = gen_part[lepton_condition2]
gen_lep_3 = gen_part[lepton_condition3]

if args.process=='tttt':
    lepton_condition4=ak.fill_none((is_lepton & is_w_4_decay & has_enough_pt),False)
    gen_lep_4 = gen_part[lepton_condition4]

#Select events with at least 3 generator level leptons that come from  top-W! ------------------------------------------------------------------------------
n_gen_lep_1, n_gen_lep_2, n_gen_lep_3 = ak.num(gen_lep_1), ak.num(gen_lep_2), ak.num(gen_lep_3) #Count the number of elements in gen_leptons 
is_3L_channel = ak.fill_none(  #Fill none's with False
        (
            (n_gen_lep_1+n_gen_lep_2+n_gen_lep_3) == 3
            
        ),
        False,
    )


if args.process=='tttt':
    n_gen_lep_4 = ak.num(gen_lep_4)
    is_3L_channel = ak.fill_none(  #Fill none's with False
        (
            (n_gen_lep_1+n_gen_lep_2+n_gen_lep_3+n_gen_lep_4 ) == 3
        ),
        False,
    )


L3_Channel_events=events[is_3L_channel]

# L3_gen_part = ak.with_name(L3_Channel_events["GenPart"], "GenParticle")
gen_part=gen_part[is_3L_channel] #save the selected GenPart 



gen_lep_1_selected = gen_lep_1[is_3L_channel]
gen_lep_1_pt= (gen_lep_1_selected.pt)
gen_lep_1_phi=(gen_lep_1_selected.phi)
gen_lep_1_eta=gen_lep_1_selected.eta
gen_lep_1_pdgId=(gen_lep_1_selected.absPdgId)

gen_lep_2_selected = gen_lep_2[is_3L_channel]
gen_lep_2_pt= (gen_lep_2_selected.pt)
gen_lep_2_phi=(gen_lep_2_selected.pt)
gen_lep_2_eta=(gen_lep_2_selected.eta)
gen_lep_2_pdgId=(gen_lep_2_selected.absPdgId)

gen_lep_3_selected = gen_lep_3[is_3L_channel]
gen_lep_3_pt= (gen_lep_3_selected.pt)
gen_lep_3_phi=(gen_lep_3_selected.phi)
gen_lep_3_eta=(gen_lep_3_selected.eta)
gen_lep_3_pdgId=(gen_lep_3_selected.absPdgId)

if args.process=='tttt':
    gen_lep_4_selected = gen_lep_4[is_3L_channel]
    gen_lep_4_pt= (gen_lep_4_selected.pt)
    gen_lep_4_phi=(gen_lep_4_selected.phi)
    gen_lep_4_eta=(gen_lep_4_selected.eta)
    gen_lep_4_pdgId=(gen_lep_4_selected.absPdgId)

gen_top_1_selected = gen_top_1[is_3L_channel]
gen_top_1_pt=np.array(gen_top_1_selected.pt) #we need to use sth otherwise ROOT can't save it as it's an ak-like array
gen_top_1_phi=np.array(gen_top_1_selected.phi)
gen_top_1_eta=np.array(gen_top_1_selected.eta)

gen_top_2_selected = gen_top_2[is_3L_channel]
gen_top_2_pt=np.array(gen_top_2_selected.pt)
gen_top_2_phi=np.array(gen_top_2_selected.phi)
gen_top_2_eta=np.array(gen_top_2_selected.eta)

gen_top_3_selected = gen_top_3[is_3L_channel]
gen_top_3_pt=np.array(gen_top_3_selected.pt)
gen_top_3_phi=np.array(gen_top_3_selected.phi)
gen_top_3_eta=np.array(gen_top_3_selected.eta)

if args.process=='tttt':
    gen_top_4_selected = gen_top_4[is_3L_channel]
    gen_top_4_pt=np.array(gen_top_4_selected.pt)
    gen_top_4_phi=np.array(gen_top_4_selected.phi)
    gen_top_4_eta=np.array(gen_top_4_selected.eta)

# Pairing the correct leptons to the ttbar
gen_top_1_Mothers=parent(gen_top_1_selected)
gen_top_2_Mothers=parent(gen_top_2_selected)
gen_top_3_Mothers=parent(gen_top_3_selected)

if args.process=='tttt':
    gen_top_4_Mothers=parent(gen_top_4_selected)

#Remember gen_part is already 3L selected

def is_same_mother(gen_top_Mother1, gen_top_Mother2 ):
        if gen_top_Mother1.index == gen_top_Mother2.index:
            return True
        else:
            return False



#Need to construct an array that will write 
Same_Mother_Top=[]
#I'll go for a loop but surely there's a better way
if args.process=='tttt':
     for i in range(len(L3_Channel_events)):
        #Make sure there's no error and they don't all have same mother
        gen_top_1_Mother=gen_top_1_Mothers[i]
        gen_top_2_Mother=gen_top_2_Mothers[i]
        gen_top_3_Mother=gen_top_3_Mothers[i]
        gen_top_4_Mother=gen_top_4_Mothers[i]
        if (is_same_mother(gen_top_1_Mother,gen_top_2_Mother) & is_same_mother(gen_top_1_Mother,gen_top_3_Mother)) | (is_same_mother(gen_top_2_Mother,gen_top_3_Mother) & is_same_mother(gen_top_1_Mother,gen_top_2_Mother)) | (is_same_mother(gen_top_2_Mother,gen_top_3_Mother) & is_same_mother(gen_top_1_Mother,gen_top_2_Mother)):
            print("error same mother for two pair combination in event " + str(i))
            Same_Mother_Top.append(0)
            continue
        if is_same_mother(gen_top_1_Mother,gen_top_2_Mother):
            Same_Mother_Top.append(1)
            continue
        if is_same_mother(gen_top_1_Mother,gen_top_3_Mother):
            Same_Mother_Top.append(2)
            continue
        if is_same_mother(gen_top_2_Mother,gen_top_3_Mother):
            Same_Mother_Top.append(3)
            continue
        if is_same_mother(gen_top_1_Mother,gen_top_4_Mother):
            Same_Mother_Top.append(5)
            continue
        if is_same_mother(gen_top_2_Mother,gen_top_4_Mother):
            Same_Mother_Top.append(6)
            continue
        if is_same_mother(gen_top_3_Mother,gen_top_4_Mother):
            Same_Mother_Top.append(7)
            continue    
        Same_Mother_Top.append(4) #append 4 if there was not found a same mother tops.
else:     
    for i in range(len(L3_Channel_events)):
         #Make sure there's no error and they don't all have same mother
        gen_top_1_Mother=gen_top_1_Mothers[i]
        gen_top_2_Mother=gen_top_2_Mothers[i]
        gen_top_3_Mother=gen_top_3_Mothers[i]
        if (is_same_mother(gen_top_1_Mother,gen_top_2_Mother) & is_same_mother(gen_top_1_Mother,gen_top_3_Mother)) | (is_same_mother(gen_top_2_Mother,gen_top_3_Mother) & is_same_mother(gen_top_1_Mother,gen_top_2_Mother)) | (is_same_mother(gen_top_2_Mother,gen_top_3_Mother) & is_same_mother(gen_top_1_Mother,gen_top_2_Mother)):
            print("error same mother for two pair combination in event " + str(i))
            Same_Mother_Top.append(0)
            continue
        if is_same_mother(gen_top_1_Mother,gen_top_2_Mother):
            Same_Mother_Top.append(1)
            continue
        if is_same_mother(gen_top_1_Mother,gen_top_3_Mother):
            Same_Mother_Top.append(2)
            continue
        if is_same_mother(gen_top_2_Mother,gen_top_3_Mother):
            Same_Mother_Top.append(3)
            continue
        Same_Mother_Top.append(4) #append 4 if there was not found a same mother tops.
    
print("len of L3_events is " +str(len(L3_Channel_events)) )
print("Same mother top length is " + str(len(Same_Mother_Top)))

if args.process=='tttt':
    with uproot.recreate(str(args.process)+"_selected_Events.root") as output_file:
        output_file["Events"] = {
                            "gen_lep_1_pt": gen_lep_1_pt, "gen_lep_1_phi": gen_lep_1_phi,"gen_lep_1_eta": gen_lep_1_eta, "gen_lep_1_pdgId": gen_lep_1_pdgId,
                            "gen_lep_2_pt": gen_lep_2_pt, "gen_lep_2_phi": gen_lep_2_phi,"gen_lep_2_eta": gen_lep_2_eta, "gen_lep_2_pdgId": gen_lep_2_pdgId,
                            "gen_lep_3_pt": gen_lep_3_pt, "gen_lep_3_phi": gen_lep_3_phi,"gen_lep_3_eta": gen_lep_3_eta, "gen_lep_3_pdgId": gen_lep_3_pdgId,
                            "gen_lep_4_pt": gen_lep_4_pt, "gen_lep_4_phi": gen_lep_4_phi,"gen_lep_4_eta": gen_lep_4_eta, "gen_lep_4_pdgId": gen_lep_4_pdgId,
                            "gen_top_1_pt": gen_top_1_pt, "gen_top_1_phi": gen_top_1_phi,"gen_top_1_eta": gen_top_1_eta, 
                            "gen_top_2_pt": gen_top_2_pt, "gen_top_2_phi": gen_top_2_phi,"gen_top_2_eta": gen_top_2_eta, 
                            "gen_top_3_pt": gen_top_3_pt, "gen_top_3_phi": gen_top_3_phi,"gen_top_3_eta": gen_top_3_eta,
                            "gen_top_4_pt": gen_top_4_pt, "gen_top_4_phi": gen_top_4_phi,"gen_top_4_eta": gen_top_4_eta,  
                            "Same_Mother_top":Same_Mother_Top,
                        }

else:   
    with uproot.recreate(str(args.process)+"_selected_Events.root")  as output_file:
        output_file["Events"] = {
                            "gen_lep_1_pt": gen_lep_1_pt, "gen_lep_1_phi": gen_lep_1_phi,"gen_lep_1_eta": gen_lep_1_eta, "gen_lep_1_pdgId": gen_lep_1_pdgId,
                            "gen_lep_2_pt": gen_lep_2_pt, "gen_lep_2_phi": gen_lep_2_phi,"gen_lep_2_eta": gen_lep_2_eta, "gen_lep_2_pdgId": gen_lep_2_pdgId,
                            "gen_lep_3_pt": gen_lep_3_pt, "gen_lep_3_phi": gen_lep_3_phi,"gen_lep_3_eta": gen_lep_3_eta, "gen_lep_3_pdgId": gen_lep_3_pdgId,
                            "gen_top_1_pt": gen_top_1_pt, "gen_top_1_phi": gen_top_1_phi,"gen_top_1_eta": gen_top_1_eta, 
                            "gen_top_2_pt": gen_top_2_pt, "gen_top_2_phi": gen_top_2_phi,"gen_top_2_eta": gen_top_2_eta, 
                            "gen_top_3_pt": gen_top_3_pt, "gen_top_3_phi": gen_top_3_phi,"gen_top_3_eta": gen_top_3_eta, 
                            "Same_Mother_top":Same_Mother_Top,
                        }


    
fileData.close()
output_file.close()
    