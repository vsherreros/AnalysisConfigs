import uproot
import awkward as ak
import numpy as np
import vector
from pocket_coffea.lib.deltaR_matching import object_matching
from coffea.nanoevents.methods.vector import LorentzVector

file = uproot.open(
	"/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/07798B6D-27E2-5940-9C12-E4E568EBD19E.root"
#	"root://cmsxrootd-kit-disk.gridka.de:1094///store/mc/RunIISummer20UL18NanoAODv9/TTbb_4f_TTTo2L2Nu_TuneCP5-Powheg-Openloops-Pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/130000/63DDFFBA-1381-C742-B38E-36F5DFB5AACB.root"
)
	
events = file["Events"]

ev=10000
GP = events.arrays(filter_name="GenPart_*", entry_stop=ev)
GJ = events.arrays(filter_name="GenJet_*", entry_stop=ev)
#El = events.arrays(filter_name="Electron", entry_stop=ev)
#Mu = events.arrays(filter_name="Muon_*", entry_stop=ev)

#print("Mu", Mu.tolist())


def find_mother(progenitorId):
	mother_idx = GP.GenPart_genPartIdxMother
	matched_particle = ak.full_like(GP.GenPart_pdgId, False) #Boolean array that will change for True if the particle comes from given progenitor
	control_array = ak.full_like(GP.GenPart_pdgId, False)  # To check if all mothers were analyzed

	while ak.any(control_array == 0):
		#Control if all mother's particle were studied
		first_mother = mother_idx == ak.full_like(GP.GenPart_pdgId, -1)
		control_array = control_array | first_mother
		
		#Check if mother is top
		new_mother = GP.GenPart_pdgId[mother_idx] == ak.full_like(GP.GenPart_pdgId, progenitorId)
		matched_particle = matched_particle  |  new_mother
	
		#Define index for next (previous mother)
		mother_idx = GP.GenPart_genPartIdxMother[mother_idx]
			
	#Boolean array with particles from top	
	return matched_particle 


"""
def find_mother(progenitorId, mother_idx):
    # Caso base: si el índice de la madre es -1, hemos llegado al final del árbol de decaimiento
    if ak.any(mother_idx == -1):
        return ak.Array(False, dtype=bool)
    
    # Verifica si la madre es el progenitor deseado
    is_progenitor = GP.GenPart_pdgId[mother_idx] == progenitorId
    
    # Si encontramos el progenitor, devolvemos True
    if ak.any(is_progenitor):
        return is_progenitor
    
    # Si no, seguimos buscando recursivamente en la madre de esta madre
    return find_mother(progenitorId, GP.GenPart_genPartIdxMother[mother_idx])
"""




### BUILD MASK

top_b, top_anti_b, higgs_b, higgs_anti_b, leptons, anti_leptons, neutrinos, anti_neutrinos = [], [], [], [], [], [], [], []
particle_mass = {"b": 4.18, "electron": 0.0005110, "nu_e": 0.0, "muon": 0.105658, "nu_mu": 0.0}
mother_idx = GP.GenPart_genPartIdxMother

# Built the mask for particle identification
from_hard_process = (GP.GenPart_statusFlags & 256) == 256  
is_last_copy = (GP.GenPart_statusFlags & 4096) == 4096 # 4096: first copy / 8192: last copy
#is_first_copy = (GP.GenPart_statusFlags & 4096) == 4096

is_b = GP.GenPart_pdgId == 5
is_anti_b = GP.GenPart_pdgId == -5

is_lepton = (GP.GenPart_pdgId == 11) | (GP.GenPart_pdgId == 13)
is_anti_lepton = (GP.GenPart_pdgId == -11) | (GP.GenPart_pdgId == -13)
is_neutrino = (GP.GenPart_pdgId == 12) | (GP.GenPart_pdgId == 14)
is_anti_neutrino = (GP.GenPart_pdgId == -12) | (GP.GenPart_pdgId == -14)
		
is_b_jet = abs(GJ.GenJet_partonFlavour) == 5 		
	

	
#print(GJ.tolist())
	
#print(is_b_jet.tolist())		

#print(GJ[is_b_jet].tolist())		

		
# comes from top or anti-top
is_from_top = find_mother(6)
is_from_anti_top = find_mother(-6)
not_from_top = (is_from_top | is_from_anti_top) == 0

#is_from_Wp = find_mother(24)
#is_from_Wm = find_mother(-24)
#is_from_higgs = find_mother(25)


#### GENERATE BRANCHES FOR IDENTIFIED PARTICLES

# b quark from top
b_from_top_mask = (from_hard_process & is_last_copy & is_b & is_from_top) != 0

GP["b_from_t"] = ak.zip(
    {
        "pt": GP.GenPart_pt[b_from_top_mask],
        "eta": GP.GenPart_eta[b_from_top_mask],
        "phi": GP.GenPart_phi[b_from_top_mask],
        "mass": ak.full_like(GP.GenPart_mass[b_from_top_mask], particle_mass["b"])
    },
    with_name="PtEtaPhiMCandidate"
)
	
# anti-b quark from anti-top
anti_b_from_anti_top_mask = (from_hard_process & is_last_copy & is_anti_b & is_from_anti_top) != 0

GP["anti_b_from_anti_t"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_b_from_anti_top_mask],
        "eta": GP.GenPart_eta[anti_b_from_anti_top_mask],
        "phi": GP.GenPart_phi[anti_b_from_anti_top_mask],
        "mass": ak.full_like(GP.GenPart_mass[anti_b_from_anti_top_mask], particle_mass["b"])
    },
    with_name="PtEtaPhiMCandidate"
)




# b quark from Higgs
b_from_higgs_mask = (from_hard_process & is_last_copy & is_b & not_from_top) != 0

GP["b_from_H"] = ak.zip(
    {
        "pt": GP.GenPart_pt[b_from_higgs_mask],
        "eta": GP.GenPart_eta[b_from_higgs_mask],
        "phi": GP.GenPart_phi[b_from_higgs_mask],
        "mass": ak.full_like(GP.GenPart_mass[b_from_higgs_mask], particle_mass["b"])
    },
    with_name="PtEtaPhiMCandidate"
)




# anti-b quark from anti-Higgs
anti_b_from_higgs_mask = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0

GP["anti_b_from_H"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_b_from_higgs_mask],
        "eta": GP.GenPart_eta[anti_b_from_higgs_mask],
        "phi": GP.GenPart_phi[anti_b_from_higgs_mask],
        "mass": ak.full_like(GP.GenPart_mass[anti_b_from_higgs_mask], particle_mass["b"])
    },
    with_name="PtEtaPhiMCandidate"
)


# lepton
lepton_mask = (from_hard_process & is_last_copy & is_lepton & is_from_anti_top) != 0

#GP["leptons"] = ak.zip(
leptons = ak.zip(
    {
        "pt": GP.GenPart_pt[lepton_mask],
        "eta": GP.GenPart_eta[lepton_mask],
        "phi": GP.GenPart_phi[lepton_mask],
        "mass": GP.GenPart_mass[lepton_mask]
    },
    with_name="PtEtaPhiMCandidate"
)

# anti-lepton
anti_lepton_mask = (from_hard_process & is_last_copy & is_anti_lepton & is_from_top) != 0

GP["anti_leptons"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_lepton_mask],
        "eta": GP.GenPart_eta[anti_lepton_mask],
        "phi": GP.GenPart_phi[anti_lepton_mask],
        "mass": GP.GenPart_mass[anti_lepton_mask]
    },
    with_name="PtEtaPhiMCandidate"
)

# neutrino
neutrino_mask = (from_hard_process & is_last_copy & is_neutrino & is_from_top) != 0

GP["neutrinos"] = ak.zip(
    {
        "pt": GP.GenPart_pt[neutrino_mask],
        "eta": GP.GenPart_eta[neutrino_mask],
        "phi": GP.GenPart_phi[neutrino_mask],
        "mass": GP.GenPart_mass[neutrino_mask]
    },
    with_name="PtEtaPhiMCandidate"
)

# anti-neutrino
anti_neutrino_mask = (from_hard_process & is_last_copy & is_anti_neutrino & is_from_anti_top) != 0

GP["anti_neutrinos"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_neutrino_mask],
        "eta": GP.GenPart_eta[anti_neutrino_mask],
        "phi": GP.GenPart_phi[anti_neutrino_mask],
        "mass": GP.GenPart_mass[anti_neutrino_mask]
    },
    with_name="PtEtaPhiMCandidate"
)



b_jets = ak.zip(
	{
		"pt": GJ[is_b_jet]["GenJet_pt"],
        "eta": GJ[is_b_jet]["GenJet_eta"],
        "phi": GJ[is_b_jet]["GenJet_phi"],
        "mass": GJ[is_b_jet]["GenJet_mass"]
	},
    with_name="PtEtaPhiMCandidate"
)

#for line in b_jets:
#	print(line.tolist())	

#print(GP.b_from_H.tolist())

#print("\nb from H last", GP.b_from_H_last.pt.tolist(), "\nb from H first", GP.b_from_H_first.pt.tolist())#, "\nleptons", GP.leptons.pt.tolist())
	
####### LORENTZ VECTORS




# For each particle type creates PtEtaPhiMCandidate object

#Dileptonic condition
all_particles_present_mask = (
    (ak.num(GP["b_from_t"]) > 0)
    & (ak.num(GP["anti_b_from_anti_t"]) > 0)
    & (ak.num(GP["b_from_H"]) > 0)
    & (ak.num(GP["anti_b_from_H"]) > 0)
    & (ak.num(leptons) > 0)
    & (ak.num(GP["anti_leptons"]) > 0)
    & (ak.num(GP["neutrinos"]) > 0)
    & (ak.num(GP["anti_neutrinos"]) > 0)
#    & (ak.num(b_jets) > 3)
)





#Improve with a loop with the keys
b_from_t_AllParticles = GP["b_from_t"][all_particles_present_mask]
anti_b_from_anti_t_AllParticles = GP["anti_b_from_anti_t"][all_particles_present_mask]
b_from_H_AllParticles = GP["b_from_H"][all_particles_present_mask]
anti_b_from_H_AllParticles = GP["anti_b_from_H"][all_particles_present_mask]

leptons_AllParticles = leptons[all_particles_present_mask]
anti_leptons_AllParticles = GP["anti_leptons"][all_particles_present_mask]
neutrinos_AllParticles = GP["neutrinos"][all_particles_present_mask]
anti_neutrinos_AllParticles = GP["anti_neutrinos"][all_particles_present_mask]

b_jets_AllParticles = b_jets[all_particles_present_mask]
b_jets_LV = vector.zip( {"pt": b_jets_AllParticles['pt'], "eta": b_jets_AllParticles['eta'], "phi": b_jets_AllParticles['phi'], "mass": b_jets_AllParticles['mass']})

#b_jets_LV = vector.backends.awkward.MomentumAwkward4D({"pt": pt, "eta": eta, "phi": phi, "mass": mass})

# When selecting events with the 8 particles, we reduce the size to that number of events
		
		
pt_index_sort = ak.argsort(b_from_t_AllParticles["pt"], ascending=False) 
b_from_t_clean =  b_from_t_AllParticles[pt_index_sort][:,0]#ak.pad_none(ak.with_name(b_from_t_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1) #[:,0] #ak.pad_none(ak.sort(GP.b_from_t.pt.tolist(), axis=1, ascending=False), 1, clip=True) # select events with all particles, sort by pt and take the higher pt

pt_index_sort = ak.argsort(anti_b_from_anti_t_AllParticles["pt"], ascending=False)
anti_b_from_anti_t_clean = anti_b_from_anti_t_AllParticles[pt_index_sort][:,0]

pt_index_sort = ak.argsort(b_from_H_AllParticles["pt"], ascending=False)
b_from_H_clean = b_from_H_AllParticles[pt_index_sort][:,0]

#print("b_jets:",b_jets_LV.pt)
#print("b_H:", (b_from_H_LV).tolist())
pt_index_sort = ak.argsort(anti_b_from_H_AllParticles["pt"], ascending=False)
anti_b_from_H_clean = anti_b_from_H_AllParticles[pt_index_sort][:,0]
#anti_b_from_H_clean_A = ak.pad_none(ak.with_name(anti_b_from_H_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1, clip=True) 


pt_index_sort = ak.argsort(leptons_AllParticles["pt"], ascending=False)
leptons_clean = leptons_AllParticles[pt_index_sort][:,0]

pt_index_sort = ak.argsort(anti_leptons_AllParticles["pt"], ascending=False)
anti_leptons_clean = anti_leptons_AllParticles[pt_index_sort][:,0]

pt_index_sort = ak.argsort(neutrinos_AllParticles["pt"], ascending=False)
neutrinos_clean = neutrinos_AllParticles[pt_index_sort][:,0]

pt_index_sort = ak.argsort(anti_neutrinos_AllParticles["pt"], ascending=False)
anti_neutrinos_clean = anti_neutrinos_AllParticles[pt_index_sort][:,0]


# b_from_t_LV = ak.to_numpy(b_from_t_clean).view(vector.MomentumNumpy4D)
# anti_b_from_anti_t_LV = ak.to_numpy(anti_b_from_anti_t_clean).view(vector.MomentumNumpy4D)
# b_from_H_LV = ak.to_numpy(b_from_H_clean).view(vector.MomentumNumpy4D)
# anti_b_from_H_LV = ak.to_numpy(anti_b_from_H_clean).view(vector.MomentumNumpy4D)



b_from_t_LV = vector.zip( {"pt": b_from_t_clean['pt'], "eta": b_from_t_clean['eta'], "phi": b_from_t_clean['phi'], "mass": b_from_t_clean['mass']})
anti_b_from_anti_t_LV = vector.zip( {"pt": anti_b_from_anti_t_clean['pt'], "eta": anti_b_from_anti_t_clean['eta'], "phi": anti_b_from_anti_t_clean['phi'], "mass": anti_b_from_anti_t_clean['mass']})
b_from_H_LV = vector.zip( {"pt": b_from_H_clean['pt'], "eta": b_from_H_clean['eta'], "phi": b_from_H_clean['phi'], "mass": b_from_H_clean['mass']})
anti_b_from_H_LV = vector.zip( {"pt": anti_b_from_H_clean['pt'], "eta": anti_b_from_H_clean['eta'], "phi": anti_b_from_H_clean['phi'], "mass": anti_b_from_H_clean['mass']})

# leptons_LV = vector.zip( {"pt": leptons_clean['pt'], "eta": leptons_clean['eta'], "phi": leptons_clean['phi'], "mass": leptons_clean['mass']})
# anti_leptons_LV = vector.zip( {"pt": anti_leptons_clean['pt'], "eta": anti_leptons_clean['eta'], "phi": anti_leptons_clean['phi'], "mass": anti_leptons_clean['mass']})
# neutrinos_LV = vector.zip( {"pt": neutrinos_clean['pt'], "eta": neutrinos_clean['eta'], "phi": neutrinos_clean['phi'], "mass": neutrinos_clean['mass']})
# anti_neutrinos_LV = vector.zip( {"pt": anti_neutrinos_clean['pt'], "eta": anti_neutrinos_clean['eta'], "phi": anti_neutrinos_clean['phi'], "mass": anti_neutrinos_clean['mass']})


"""
dR_b_t_mask =  b_jets_LV.deltaR(b_from_t_LV)<= 0.5 #& ak.min(b_jets_LV.deltaR(b_from_t_LV))
dR_anti_b_t_mask =  b_jets_LV.deltaR(anti_b_from_anti_t_LV)<= 0.5
dR_b_H_mask =  b_jets_LV.deltaR(b_from_H_LV)<= 0.5
dR_anti_b_H_mask =  b_jets_LV.deltaR(anti_b_from_H_LV)<= 0.5
"""



#############################  JET MATCHING  #############################################
min_dR = 0.5

dR_b_t = b_jets_LV.deltaR(b_from_t_LV) #dR for b_t
min_bt_mask = ( dR_b_t == ak.min(dR_b_t, axis=1) ) & (dR_b_t < min_dR) #mask to find the nice b-jets
b_top_jet = b_jets_LV[min_bt_mask] #b-jet assignment

#remove the assigned jets
b_jets_LV = b_jets_LV[~min_bt_mask]

dR_anti_b_t = b_jets_LV.deltaR(anti_b_from_anti_t_LV)
min_anti_bt_mask = ( dR_anti_b_t == ak.min(dR_anti_b_t, axis=1) ) & (dR_anti_b_t < min_dR)
anti_b_top_jet = b_jets_LV[min_anti_bt_mask]

#remove the assigned jets
b_jets_LV = b_jets_LV[~min_anti_bt_mask] 

dR_b_H = b_jets_LV.deltaR(b_from_H_LV)
min_bH_mask = ( dR_b_H == ak.min(dR_b_H, axis=1) ) & (dR_b_H  < min_dR)
b_H_jet = b_jets_LV[min_bH_mask]

#remove the assigned jets
b_jets_LV = b_jets_LV[~min_bH_mask] 

dR_anti_b_H = b_jets_LV.deltaR(anti_b_from_H_LV)
min_anti_bH_mask = ( dR_anti_b_H == ak.min(dR_anti_b_H, axis=1) ) & (dR_anti_b_H  < min_dR)
anti_b_H_jet = b_jets_LV[min_anti_bH_mask]

##########################################################################################

#Mask to select events with 4 kind of b-jets
all_b_jets_mask = (
    (ak.num(b_top_jet) > 0)
    & (ak.num(anti_b_top_jet) > 0)
    & (ak.num(b_H_jet) > 0)
    & (ak.num(anti_b_H_jet) > 0)
	)

#print(len(ak.num(b_top_jet)), len(ak.num(anti_b_top_jet)), len(ak.num(b_H_jet)), len(ak.num(anti_b_H_jet)))

b_top_jet = b_top_jet[all_b_jets_mask]
anti_b_top_jet = anti_b_top_jet[all_b_jets_mask]
b_H_jet = b_H_jet[all_b_jets_mask]
anti_b_H_jet = anti_b_H_jet[all_b_jets_mask]

leptons_clean = leptons_clean[all_b_jets_mask]
anti_leptons_clean = anti_leptons_clean[all_b_jets_mask]
neutrinos_clean = neutrinos_clean[all_b_jets_mask]
anti_neutrinos_clean = anti_neutrinos_clean[all_b_jets_mask]


#for event in leptons_LV:
#	print("Ev len", (event))
#print(ak.num(leptons_LV, axis=1))
#print("############", ak.any( ak.num(leptons_LV) > 1))

#print(ak.any(all_b_jets_mask) == 0)
#for event in leptons_LV:
#	print(ak.count(event))
#print(ak.any(leptons_LV) == None)
#print(len(anti_leptons_LV), len(all_b_jets_mask), len(neutrinos_LV))
#print(len(b_top_jet))

#leptons_LV = leptons_LV[all_b_jets_mask]
#anti_leptons_LV = anti_leptons_LV[all_b_jets_mask]
#neutrinos_LV = neutrinos_LV[all_b_jets_mask]
#anti_neutrinos_LV = anti_neutrinos_LV[all_b_jets_mask]

leptons_LV = ak.to_numpy(leptons_clean).view(vector.MomentumNumpy4D)
anti_leptons_LV = ak.to_numpy(anti_leptons_clean).view(vector.MomentumNumpy4D)
neutrinos_LV = ak.to_numpy(neutrinos_clean).view(vector.MomentumNumpy4D)
anti_neutrinos_LV = ak.to_numpy(anti_neutrinos_clean).view(vector.MomentumNumpy4D)

"""
any_length_gt_1 = any(len(subarray) < 1 for subarray in b_H)
print(any_length_gt_1)
indices = ak.where([len(subarray) < 1 for subarray in b_H])[0]

print("Índices donde la longitud es mayor que 1:", len(indices))

for index in indices:
    print(f"Contenido en el índice {index}: {(b_H[index]).tolist()}")


#print( ak.any(len(b_H) > 1, axis=0) )


#if ak.any(len(b_H) > 1, axis=1):
	#print("Oh no")
#print(b_top.tolist(), "\n",anti_b_top.tolist(),"\n",b_H.tolist(), "\n",anti_b_H.tolist)
"""
####CREATE LV



#COMBINE LV
Wp_LV = anti_leptons_LV + neutrinos_LV
Wn_LV = leptons_LV + anti_neutrinos_LV
top_LV = Wp_LV + b_top_jet
anti_top_LV = Wn_LV + anti_b_top_jet
#higgs_LV = b_from_H_LV + anti_b_from_H_LV
H = b_H_jet + anti_b_H_jet
#b = b_from_H_LV
#print(H.mass)


#print(h.deltaR(b))

print("###########################")
# cos (theta*)
b_ZMF = b_H_jet.boostCM_of_p4(H)
cos_theta_b_h = H.to_beta3().unit().dot(b_ZMF.to_beta3().unit())

### LVector boosting ### 
tt_LV = top_LV + anti_top_LV
top_ZMF = top_LV.boostCM_of_p4(tt_LV)
anti_top_ZMF = anti_top_LV.boostCM_of_p4(tt_LV)
leptons_ZMF = leptons_LV.boostCM_of(tt_LV)
anti_leptons_ZMF = anti_leptons_LV.boostCM_of(tt_LV)

lepton_ZMF = leptons_ZMF.boostCM_of_p4(anti_top_ZMF)
anti_lepton_ZMF = anti_leptons_ZMF.boostCM_of_p4(top_ZMF)

L = lepton_ZMF.to_beta3().unit()
L_bar = anti_lepton_ZMF.to_beta3().unit()

### Orthonormal Basis ###
pvect_p = vector.VectorObject3D(x=0, y=0, z=1)
kvect_p = top_ZMF.to_beta3().unit()
y_p = pvect_p.dot(kvect_p)
r_p = np.sqrt(1 - y_p*y_p)
rvect_p = np.sign(y_p) * 1 / r_p * (pvect_p - y_p * kvect_p)
nvect_p = np.sign(y_p) * 1 / r_p * (pvect_p.cross(kvect_p))


### Spin Observables ###
c_hel = L.dot(L_bar)
LLbar_DeltaPhi = abs(L.deltaphi(L_bar))
LLbar_DeltaEta = abs(L.deltaeta(L_bar))
b1k_plus_b2k = L.dot(kvect_p) +  L_bar.dot(kvect_p)
b1r_plus_b2r = L.dot(rvect_p) +  L_bar.dot(rvect_p)
b1n_plus_b2n = L.dot(nvect_p) +  L_bar.dot(nvect_p)
b1k_minus_b2k = L.dot(kvect_p) -  L_bar.dot(kvect_p)
b1r_minus_b2r = L.dot(rvect_p) -  L_bar.dot(rvect_p)
b1n_minus_b2n = L.dot(nvect_p) -  L_bar.dot(nvect_p) 
Ckk = L.dot(kvect_p) * L_bar.dot(kvect_p)
Crr = L.dot(rvect_p) * L_bar.dot(rvect_p)
Cnn = L.dot(nvect_p) * L_bar.dot(nvect_p)
Crk_plus_Ckr = L.dot(rvect_p) * L_bar.dot(kvect_p) + L.dot(kvect_p) * L_bar.dot(rvect_p)
Crn_plus_Cnr = L.dot(rvect_p) * L_bar.dot(nvect_p) + L.dot(nvect_p) * L_bar.dot(rvect_p)
Ckn_plus_Cnk = L.dot(kvect_p) * L_bar.dot(nvect_p) + L.dot(nvect_p) * L_bar.dot(kvect_p)
Crk_minus_Ckr = L.dot(rvect_p) * L_bar.dot(kvect_p) - L.dot(kvect_p) * L_bar.dot(rvect_p)
Crn_minus_Cnr = L.dot(rvect_p) * L_bar.dot(nvect_p) - L.dot(nvect_p) * L_bar.dot(rvect_p)
Ckn_minus_Cnk = L.dot(kvect_p) * L_bar.dot(nvect_p) - L.dot(nvect_p) * L_bar.dot(kvect_p)

####  min DeltaPhi in LAB  #####

#calculate Phi between two particles
TTbar_phi = abs( top_LV.deltaphi(anti_top_LV) )
THiggs_phi = abs( top_LV.deltaphi(H) )
TbarHiggs_phi = abs( anti_top_LV.unit().deltaphi(H.unit()) )

# def min_deltaPhi()
Phi_combined = np.concatenate((TTbar_phi , THiggs_phi , TbarHiggs_phi ), axis=1)
min_Phi = np.min(Phi_combined, axis=1)


####   min DeltaPhi in ttH_ZMF   #####

ttH_ZMF = top_LV + anti_top_LV + H

top_ttH_ZMF = top_LV.boostCM_of_p4(tt_LV)
anti_top_ttH_ZMF = anti_top_LV.boostCM_of_p4(tt_LV)
H_ttH_ZMF = H.boostCM_of_p4(tt_LV)

#print(top_ttH_ZMF.phi.tolist())
#print(anti_top_ttH_ZMF.phi.tolist())
#print(H_ttH_ZMF.phi.tolist())

TTbar_phi_ttH_ZMF = abs( top_ttH_ZMF.deltaphi(anti_top_ttH_ZMF) )
THiggs_phi_ttH_ZMF = abs( top_ttH_ZMF.deltaphi(H_ttH_ZMF) )
TbarHiggs_phi_ttH_ZMF = abs( anti_top_ttH_ZMF.deltaphi(H_ttH_ZMF) )

TTbar_eta_ttH_ZMF = abs( top_ttH_ZMF.deltaeta(anti_top_ttH_ZMF) )
THiggs_eta_ttH_ZMF = abs( top_ttH_ZMF.deltaeta(H_ttH_ZMF) )
TbarHiggs_eta_ttH_ZMF = abs( anti_top_ttH_ZMF.deltaeta(H_ttH_ZMF) )


print(THiggs_phi_ttH_ZMF)

Phi_combined_ttH_ZMF = np.concatenate((TTbar_phi_ttH_ZMF , THiggs_phi_ttH_ZMF , TbarHiggs_phi_ttH_ZMF ), axis=1)

#print(Phi_combined_ttH_ZMF.tolist())
min_Phi_ttH_ZMF = np.min(Phi_combined_ttH_ZMF, axis=1)

total_Phi =  TTbar_phi_ttH_ZMF + THiggs_phi_ttH_ZMF + TbarHiggs_phi_ttH_ZMF
data_list = ak.to_list(total_Phi)

filter_phi = [x for x in THiggs_phi_ttH_ZMF if x is not None]
#print(filter_phi)


import hist
import matplotlib.pyplot as plt

hist.Hist(hist.axis.Regular(20, 0, 4, label=r"$c_{hel}$ []")).fill(
    ak.flatten(filter_phi)
).plot()

plt.show()





"""


Higgs = ak.zip(
    {
        "pt": higgs_LV.pt,
        "eta": higgs_LV.eta,
        "phi": higgs_LV.phi,
        "mass": higgs_LV.mass,

	},
	with_name="PtEtaPhiMRatioCandidate"
)








import hist
import matplotlib.pyplot as plt

hist.Hist(hist.axis.Regular(50, 0, 2, label="mass [GeV]")).fill(
    Higgs.ratio
).plot()

plt.show()




import hist
import matplotlib.pyplot as plt

hist.Hist(hist.axis.Regular(50, 0, 140, label="mass [GeV]")).fill(
    Higgs.mass
).plot()

plt.show()






GP["Higgs"] = ak.with_field(
    {
        "pt": higgs_LV.pt,
        "eta": higgs_LV.eta,
        "phi": higgs_LV.phi,
        "mass": higgs_LV.mass
	},
	with_name="PtEtaPhiMCandidate"
), len(GP.GenPart_pt), axis=0)


print("\n W+:")
for line in Wp_LV:
    print(line)
    
print("\n W-:")
for line in Wn_LV:
    print(line)    
    
print("\n top:")
for line in top_LV:
    print(line)
    
print("\n anti-top:")
for line in anti_top_LV:
    print(line)    
 
    
print("\n Higgs:")
for line in higgs_LV["mass"]:
    print(line) 

print("Higgs events:",len(higgs_LV["mass"]),   "\nMean mass:", np.mean(np.array(higgs_LV["mass"])))

	


for line in Mu:
	print(line)



def get_dilepton(electrons, muons, transverse=False):

	fields = {
		"pt": 0.,
		"eta": 0.,
		"phi": 0.,
		"mass": 0.,
		"charge": 0.,
	}
	
	leptons = ak.pad_none(ak.with_name(ak.concatenate([ muons[:, 0:2], electrons[:, 0:2]], axis=1), "PtEtaPhiMCandidate"), 2)
	nlep =  ak.num(leptons[~ak.is_none(leptons, axis=1)])
	ll = leptons[:,0] + leptons[:,1]

	for var in fields.keys():
		fields[var] = ak.where(
			(nlep == 2),
			getattr(ll, var),
			fields[var]
		)
		
	fields["deltaR"] = ak.where(
		(nlep == 2), leptons[:,0].delta_r(leptons[:,1]), -1)

	if transverse:
		fields["eta"] = ak.zeros_like(fields["pt"])
	dileptons = ak.zip(fields, with_name="PtEtaPhiMCandidate")

	return dileptons
	
def lepton_selection(events, lepton_flavour, params):

	leptons = events[lepton_flavour]
	cuts = params.object_preselection[lepton_flavour]
	# Requirements on pT and eta
	passes_eta = abs(leptons.eta) < cuts["eta"]
	passes_pt = leptons.pt > cuts["pt"]

	if lepton_flavour == "Electron":
		# Requirements on SuperCluster eta, isolation and id
		etaSC = abs(leptons.deltaEtaSC + leptons.eta)
		passes_SC = np.invert((etaSC >= 1.4442) & (etaSC <= 1.5660))
		passes_iso = leptons.pfRelIso03_all < cuts["iso"]
		passes_id = leptons[cuts['id']] == True

		good_leptons = passes_eta & passes_pt & passes_SC & passes_iso & passes_id

	elif lepton_flavour == "Muon":
		# Requirements on isolation and id
		passes_iso = leptons.pfRelIso04_all < cuts["iso"]
		passes_id = leptons[cuts['id']] == True

		good_leptons = passes_eta & passes_pt & passes_iso & passes_id

	return leptons[good_leptons]	


params = {
  Muon:
    pt: 15
    eta: 2.4
    iso: 0.25  #PFIsoLoose
    id: tightId

  Electron:
    pt: 15
    eta: 2.4
    iso: 0.06
    id: mvaFall17V2Iso_WP80


Mu["MuonGood"] = lepton_selection(Mu, "Muon", params)
El["ElectronGood"] = lepton_selection(El, "Electron", params)  	

ll = get_dilepton(El.ElectronGood, Mu.MuonGood)
		
"""      
		





        
       

