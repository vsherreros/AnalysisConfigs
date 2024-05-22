import uproot
import awkward as ak
import numpy as np
import vector


file = uproot.open(
	"/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/50000/07798B6D-27E2-5940-9C12-E4E568EBD19E.root"
)
	
events = file["Events"]

ev=10000
GP = events.arrays(filter_name="GenPart_*", entry_stop=ev)
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
mother_idx = GP.GenPart_genPartIdxMother

# Built the mask for particle identification
from_hard_process = (GP.GenPart_statusFlags & 256) == 256  
is_last_copy = (GP.GenPart_statusFlags & 8192) == 8192 # 4096: first copy / 8192: last copy
is_first_copy = (GP.GenPart_statusFlags & 4096) == 4096

is_b = GP.GenPart_pdgId == 5
is_anti_b = GP.GenPart_pdgId == -5

is_lepton = (GP.GenPart_pdgId == 11) | (GP.GenPart_pdgId == 13)
is_anti_lepton = (GP.GenPart_pdgId == -11) | (GP.GenPart_pdgId == -13)
is_neutrino = (GP.GenPart_pdgId == 12) | (GP.GenPart_pdgId == 14)
is_anti_neutrino = (GP.GenPart_pdgId == -12) | (GP.GenPart_pdgId == -14)
		
#is_from_higgs = find_mother(25)
		
# comes from top or anti-top
is_from_top = find_mother(6)
is_from_anti_top = find_mother(-6)

not_from_top = (is_from_top | is_from_anti_top) == 0
#is_from_Wp = find_mother(24)
#is_from_Wm = find_mother(-24)


#### GENERATE BRANCHES FOR IDENTIFIED PARTICLES

# b quark from top
b_from_top_mask = (from_hard_process & is_last_copy & is_b & is_from_top) != 0

GP["b_from_t"] = ak.zip(
    {
        "pt": GP.GenPart_pt[b_from_top_mask],
        "eta": GP.GenPart_eta[b_from_top_mask],
        "phi": GP.GenPart_phi[b_from_top_mask],
        "mass": GP.GenPart_mass[b_from_top_mask]
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
        "mass": GP.GenPart_mass[anti_b_from_anti_top_mask]
    },
    with_name="PtEtaPhiMCandidate"
)

# b quark from Higgs
b_from_higgs_mask_last = (from_hard_process & is_last_copy & is_b & not_from_top) != 0

GP["b_from_H_last"] = ak.zip(
    {
        "pt": GP.GenPart_pt[b_from_higgs_mask_last],
        "eta": GP.GenPart_eta[b_from_higgs_mask_last],
        "phi": GP.GenPart_phi[b_from_higgs_mask_last],
        "mass": GP.GenPart_mass[b_from_higgs_mask_last]
    },
    with_name="PtEtaPhiMCandidate"
)


b_from_higgs_mask_first = (from_hard_process & is_first_copy & is_b & not_from_top) != 0

GP["b_from_H_first"] = ak.zip(
    {
        "pt": GP.GenPart_pt[b_from_higgs_mask_first],
        "eta": GP.GenPart_eta[b_from_higgs_mask_first],
        "phi": GP.GenPart_phi[b_from_higgs_mask_first],
        "mass": GP.GenPart_mass[b_from_higgs_mask_first]
    },
    with_name="PtEtaPhiMCandidate"
)

#print("Mother idx:", GP.GenPart_genPartIdxMother[0] )
#print("Event:", GP.GenPart_pt[0][27:29] .tolist())
#print("Mask:", ak.argmax(b_from_higgs_mask[0]).tolist() )
#print(GP.b_from_H.pt.tolist())

# anti-b quark from anti-Higgs
anti_b_from_higgs_mask_last = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0

GP["anti_b_from_H_last"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_b_from_higgs_mask_last],
        "eta": GP.GenPart_eta[anti_b_from_higgs_mask_last],
        "phi": GP.GenPart_phi[anti_b_from_higgs_mask_last],
        "mass": GP.GenPart_mass[anti_b_from_higgs_mask_last]
    },
    with_name="PtEtaPhiMCandidate"
)

anti_b_from_higgs_mask_first = (from_hard_process & is_first_copy & is_anti_b & not_from_top) != 0

GP["anti_b_from_H_first"] = ak.zip(
    {
        "pt": GP.GenPart_pt[anti_b_from_higgs_mask_first],
        "eta": GP.GenPart_eta[anti_b_from_higgs_mask_first],
        "phi": GP.GenPart_phi[anti_b_from_higgs_mask_first],
        "mass": GP.GenPart_mass[anti_b_from_higgs_mask_first]
    },
    with_name="PtEtaPhiMCandidate"
)

# lepton
lepton_mask = (from_hard_process & is_last_copy & is_lepton & is_from_anti_top) != 0

GP["leptons"] = ak.zip(
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
	
	
#print(GP.b_from_H.tolist())	
	
#print("\nb from H last", GP.b_from_H_last.pt.tolist(), "\nb from H first", GP.b_from_H_first.pt.tolist())#, "\nleptons", GP.leptons.pt.tolist())
	
####### LORENTZ VECTORS


# For each particle type creates PtEtaPhiMCandidate object

#Dileptonic condition
all_particles_present_mask = (
    (ak.num(GP["b_from_t"]) > 0)
    & (ak.num(GP["anti_b_from_anti_t"]) > 0)
    & (ak.num(GP["b_from_H_last"]) > 0)
    & (ak.num(GP["anti_b_from_H_last"]) > 0)
    & (ak.num(GP["b_from_H_first"]) > 0)
    & (ak.num(GP["anti_b_from_H_first"]) > 0)
    & (ak.num(GP["leptons"]) > 0)
    & (ak.num(GP["anti_leptons"]) > 0)
    & (ak.num(GP["neutrinos"]) > 0)
    & (ak.num(GP["anti_neutrinos"]) > 0)
)


all_particles_present_first_mask = (
    (ak.num(GP["b_from_t"]) > 0)
    & (ak.num(GP["anti_b_from_anti_t"]) > 0)
    & (ak.num(GP["b_from_H_first"]) > 0)
    & (ak.num(GP["anti_b_from_H_first"]) > 0)
    & (ak.num(GP["leptons"]) > 0)
    & (ak.num(GP["anti_leptons"]) > 0)
    & (ak.num(GP["neutrinos"]) > 0)
    & (ak.num(GP["anti_neutrinos"]) > 0)
)

#print("first", all_particles_present_first_mask, "\nlast", all_particles_present_mask)

#ratio = GP.b_from_H_last.pt / GP.b_from_H_first.pt 
#print(ratio)
#for line in GP.b_from_t:
#	if len(line) > 1:
#		print(line.tolist(), len(line))

#print(all_particles_present_mask)

#Improve with a loop with the keys
b_from_t_AllParticles = GP["b_from_t"][all_particles_present_mask]
anti_b_from_anti_t_AllParticles = GP["anti_b_from_anti_t"][all_particles_present_mask]

b_from_H_last_AllParticles = GP["b_from_H_last"][all_particles_present_mask]
anti_b_from_H_last_AllParticles = GP["anti_b_from_H_last"][all_particles_present_mask]
b_from_H_first_AllParticles = GP["b_from_H_first"][all_particles_present_mask]
anti_b_from_H_first_AllParticles = GP["anti_b_from_H_first"][all_particles_present_mask]

leptons_AllParticles = GP["leptons"][all_particles_present_mask]
anti_leptons_AllParticles = GP["anti_leptons"][all_particles_present_mask]
neutrinos_AllParticles = GP["neutrinos"][all_particles_present_mask]
anti_neutrinos_AllParticles = GP["anti_neutrinos"][all_particles_present_mask]

# When selecting events with the 8 particles, we reduce the size to that number of events
"""
print("all mask")
for line in b_from_t_AllParticles:
	if len(line) > 1:
		print(line.tolist(), len(line))
		

for line in anti_b_from_anti_t_AllParticles:
	if len(line) > 1:
		print(line.tolist(), len(line))		


for line in b_from_H_AllParticles:
	if len(line) > 1:
		print("b_H:",line.tolist(), len(line))
"""


"""
for line in leptons_AllParticles:
	if len(line) > 1:
		print(line.tolist(), len(line))
"""		
		
pt_index_sort = ak.argsort(b_from_t_AllParticles["pt"], ascending=False) 
b_from_t_clean = ak.pad_none(ak.with_name(b_from_t_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1) #[:,0] #ak.pad_none(ak.sort(GP.b_from_t.pt.tolist(), axis=1, ascending=False), 1, clip=True) # select events with all particles, sort by pt and take the higher pt
b_from_t_LV = ak.to_numpy(b_from_t_clean).view(vector.MomentumNumpy4D)
#print(b_from_t_clean[0].tolist())
#print(b_from_t_LV[0])

pt_index_sort = ak.argsort(anti_b_from_anti_t_AllParticles["pt"], ascending=False)
anti_b_from_anti_t_clean = anti_b_from_anti_t_AllParticles[pt_index_sort][:,0]
anti_b_from_anti_t_LV = ak.to_numpy(anti_b_from_anti_t_clean).view(vector.MomentumNumpy4D)
#print(anti_b_from_anti_t_clean[0].tolist)
#print(anti_b_from_anti_t_clean.type)

pt_index_sort = ak.argsort(b_from_H_last_AllParticles["pt"], ascending=False)
b_from_H_clean_last = b_from_H_last_AllParticles[pt_index_sort][:,0]
#b_from_H_clean_A = ak.pad_none(ak.with_name(b_from_H_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1, clip=True) 
b_from_H_last_LV = ak.to_numpy(b_from_H_clean_last).view(vector.MomentumNumpy4D)

pt_index_sort = ak.argsort(anti_b_from_H_last_AllParticles["pt"], ascending=False)
anti_b_from_H_clean_last = anti_b_from_H_last_AllParticles[pt_index_sort][:,0]
#anti_b_from_H_clean_A = ak.pad_none(ak.with_name(anti_b_from_H_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1, clip=True) 
anti_b_from_H_last_LV = ak.to_numpy(anti_b_from_H_clean_last).view(vector.MomentumNumpy4D)


#for line in b_from_H_first_AllParticles["pt"]:
	#print(line)
	
#print("Sorted")	

pt_index_sort = ak.argsort(b_from_H_first_AllParticles["pt"], ascending=False)
#for line in b_from_H_first_AllParticles[pt_index_sort]:
#	print(line)
#print("zero:", b_from_H_first_AllParticles[pt_index_sort][0])
b_from_H_clean_first = b_from_H_first_AllParticles[pt_index_sort][:,0]
#b_from_H_clean_A = ak.pad_none(ak.with_name(b_from_H_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1, clip=True) 
b_from_H_first_LV = ak.to_numpy(b_from_H_clean_first).view(vector.MomentumNumpy4D)


pt_index_sort = ak.argsort(anti_b_from_H_first_AllParticles["pt"], ascending=False)
anti_b_from_H_clean_first = anti_b_from_H_first_AllParticles[pt_index_sort][:,0]
#anti_b_from_H_clean_A = ak.pad_none(ak.with_name(anti_b_from_H_AllParticles[pt_index_sort], "PtEtaPhiMCandidate"), 1, clip=True) 
anti_b_from_H_first_LV = ak.to_numpy(anti_b_from_H_clean_first).view(vector.MomentumNumpy4D)



#sum_R = b_from_H_clean_R + anti_b_from_H_clean_R
#sum_A = b_from_H_clean_A[:,0] + anti_b_from_H_clean_A[:,0]

# pt_index_sort = ak.argsort(leptons_AllParticles["pt"], ascending=False)
# leptons_clean = leptons_AllParticles[pt_index_sort][:,0]
# leptons_LV = ak.to_numpy(leptons_clean).view(vector.MomentumNumpy4D)

# pt_index_sort = ak.argsort(anti_leptons_AllParticles["pt"], ascending=False)
# anti_leptons_clean = anti_leptons_AllParticles[pt_index_sort][:,0]
# anti_leptons_LV = ak.to_numpy(anti_leptons_clean).view(vector.MomentumNumpy4D)

# pt_index_sort = ak.argsort(neutrinos_AllParticles["pt"], ascending=False)
# neutrinos_clean = neutrinos_AllParticles[pt_index_sort][:,0]
# neutrinos_LV = ak.to_numpy(neutrinos_clean).view(vector.MomentumNumpy4D)

# pt_index_sort = ak.argsort(anti_neutrinos_AllParticles["pt"], ascending=False)
# anti_neutrinos_clean = anti_neutrinos_AllParticles[pt_index_sort][:,0]
# anti_neutrinos_LV = ak.to_numpy(anti_neutrinos_clean).view(vector.MomentumNumpy4D)




#COMBINE LV
#Wp_LV = anti_leptons_LV + neutrinos_LV
#Wn_LV = leptons_LV + anti_neutrinos_LV
#top_LV = Wp_LV + b_from_t_LV
#anti_top_LV = Wn_LV + anti_b_from_anti_t_LV
higgs_last_LV = b_from_H_last_LV + anti_b_from_H_last_LV
higgs_first_LV = b_from_H_first_LV + anti_b_from_H_first_LV





higgs_first_LV.pt[0]=0
higgs_last_LV.pt[0]=0.5




#print("\nH-last", higgs_last_LV.pt.tolist(), "\nH-first", higgs_first_LV.pt.tolist())
import numpy as np

#print((higgs_first_LV.pt == 0) & (higgs_last_LV.pt == 0) )

higgs_last_pt_fixed = np.where( ( (higgs_first_LV.pt == 0) & (higgs_last_LV.pt == 0) ), 1, higgs_last_LV.pt)
higgs_last_pt_fixed = np.where( ( (higgs_first_LV.pt == 0) & (higgs_last_LV.pt >= 0) ), np.inf, higgs_last_LV.pt)
#print(higgs_last_pt_fixed )

higgs_first_pt_fixed = np.where(higgs_first_LV.pt == 0, 1, higgs_first_LV.pt)
#print(higgs_first_pt_fixed )
ratio = higgs_last_pt_fixed / higgs_first_pt_fixed



#rint(ratio) 





#print("Higgs_pt:", higgs_LV.pt)


Higgs = ak.zip(
    {
        "pt_last": higgs_last_LV.pt,
        "eta_last": higgs_last_LV.eta,
        "phi_last": higgs_last_LV.phi,
        "mass_last": higgs_last_LV.mass,
        "pt_first": higgs_first_LV.pt,
		"eta_first": higgs_first_LV.eta,
		"phi_first": higgs_first_LV.phi,
		"mass_first": higgs_first_LV.mass,
        "ratio": ratio
	},
	with_name="PtEtaPhiMRatioCandidate"
)

#print(ratio.tolist())


import hist
import matplotlib.pyplot as plt

hist.Hist(hist.axis.Regular(50, 0, 2, label="R_{last/first}")).fill(
    Higgs.ratio
).plot()

plt.show()


"""
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
		





        
       

