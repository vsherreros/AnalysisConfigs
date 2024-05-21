import sys
import awkward as ak
import numba
import numpy as np

import vector

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.workflows.tthbb_base_processor import BaseProcessorABC

				

class MCProcessor(BaseProcessorABC):
	
	def __init__(self, cfg: Configurator):
		super().__init__(cfg)
		
	def apply_object_preselection(self, variation):
		pass

	def count_objects(self, variation):
		pass
     
	def do_gen_selection(self):
        # Initialize empty lists for particle types
		#top_b, top_anti_b, higgs_b, higgs_anti_b, leptons, anti_leptons, neutrinos, anti_neutrinos = [], [], [], [], [], [], [], []

		
		def find_mother(self, progenitorId):
			mother_idx = self.events.GenPart.genPartIdxMother
			matched_particle = ak.full_like(self.events.GenPart.pdgId, False) #Boolean array that will change for True if the particle comes from given progenitor
			control_array = ak.full_like(self.events.GenPart.pdgId, False)  # To check if all mothers were analyzed
			
			i=0
			while ak.any(control_array == 0): #While... all mothers were not analyzed
				i+=1
				
				#Control if all mother's particle were studied
				first_mother = mother_idx == ak.full_like(self.events.GenPart.pdgId, -1)
				control_array = control_array | first_mother
				
				#Check if mother is top
				new_mother = self.events.GenPart.pdgId[mother_idx] == ak.full_like(self.events.GenPart.pdgId, progenitorId)
				matched_particle = matched_particle  |  new_mother
				
				#Define index for next (previous mother)
				mother_idx = self.events.GenPart.genPartIdxMother[mother_idx]
				
				if i==100:
					print("Warning: Posible error finding mother of particle")
					break		
			return matched_particle 
		
		
		
		mother_idx = self.events.GenPart.genPartIdxMother

        # Built the mask for particle identification
		from_hard_process = (self.events.GenPart.statusFlags & 256) == 256  
		is_last_copy = (self.events.GenPart.statusFlags & 8192) == 8192  # 4096: first copy / 8192: last copy
		is_first_copy = (self.events.GenPart.statusFlags & 4096) == 4096  

		is_b = self.events.GenPart.pdgId == 5
		is_anti_b = self.events.GenPart.pdgId == -5

		is_lepton = (self.events.GenPart.pdgId == 11) | (self.events.GenPart.pdgId == 13)
		is_anti_lepton = (self.events.GenPart.pdgId == -11) | (self.events.GenPart.pdgId == -13)
		is_neutrino = (self.events.GenPart.pdgId == 12) | (self.events.GenPart.pdgId == 14)
		is_anti_neutrino = (self.events.GenPart.pdgId == -12) | (self.events.GenPart.pdgId == -14)
		
		
		
		# comes from top or anti-top
		is_from_top = find_mother(self, 6)
		is_from_anti_top = find_mother(self, -6)
		
		not_from_top = (is_from_top | is_from_anti_top) == 0
		#is_from_higgs = find_mother(self, 25)
        
		#is_from_Wp = find_mother(24)
		#is_from_Wm = find_mother(-24)


#### GENERATE BRANCHES FOR IDENTIFIED PARTICLES

		# b quark from top
		b_from_top_mask = (from_hard_process & is_last_copy & is_b & is_from_top) != 0
		
		


		self.events["b_from_t"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[b_from_top_mask],
				"eta": self.events.GenPart.eta[b_from_top_mask],
				"phi": self.events.GenPart.phi[b_from_top_mask],
				"mass": self.events.GenPart.mass[b_from_top_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)
		
	
	
		# anti-b quark from anti-top
		anti_b_from_anti_top_mask = (from_hard_process & is_last_copy & is_anti_b & is_from_anti_top) != 0

		self.events["anti_b_from_anti_t"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[anti_b_from_anti_top_mask],
				"eta": self.events.GenPart.eta[anti_b_from_anti_top_mask],
				"phi": self.events.GenPart.phi[anti_b_from_anti_top_mask],
				"mass": self.events.GenPart.mass[anti_b_from_anti_top_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)


		# b quark from Higgs
		# b_from_higgs_mask = (from_hard_process & is_last_copy & is_b & not_from_top) != 0

		# self.events["b_from_H"] = ak.zip(
			# {
				# "pt": self.events.GenPart.pt[b_from_higgs_mask],
				# "eta": self.events.GenPart.eta[b_from_higgs_mask],
				# "phi": self.events.GenPart.phi[b_from_higgs_mask],
				# "mass": self.events.GenPart.mass[b_from_higgs_mask]
			# },
			# with_name="PtEtaPhiMCandidate"
		# )


		# # anti-b quark from anti-Higgs
		# anti_b_from_higgs_mask = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0

		# self.events["anti_b_from_H"] = ak.zip(
			# {
				# "pt": self.events.GenPart.pt[anti_b_from_higgs_mask],
				# "eta": self.events.GenPart.eta[anti_b_from_higgs_mask],
				# "phi": self.events.GenPart.phi[anti_b_from_higgs_mask],
				# "mass": self.events.GenPart.mass[anti_b_from_higgs_mask]
			# },
			# with_name="PtEtaPhiMCandidate"
		# )

		b_from_higgs_mask_last = (from_hard_process & is_last_copy & is_b & not_from_top) != 0

		self.events["b_from_H_last"] = ak.zip(
		    {
		        "pt": self.events.GenPart.pt[b_from_higgs_mask_last],
		        "eta": self.events.GenPart.eta[b_from_higgs_mask_last],
		        "phi": self.events.GenPart.phi[b_from_higgs_mask_last],
		        "mass": self.events.GenPart.mass[b_from_higgs_mask_last]
		    },
		    with_name="PtEtaPhiMCandidate"
		)
		
		
		b_from_higgs_mask_first = (from_hard_process & is_first_copy & is_b & not_from_top) != 0
		
		self.events["b_from_H_first"] = ak.zip(
		    {
		        "pt": self.events.GenPart.pt[b_from_higgs_mask_first],
		        "eta": self.events.GenPart.eta[b_from_higgs_mask_first],
		        "phi": self.events.GenPart.phi[b_from_higgs_mask_first],
		        "mass": self.events.GenPart.mass[b_from_higgs_mask_first]
		    },
		    with_name="PtEtaPhiMCandidate"
		)
		
		
		# anti-b quark from anti-Higgs
		anti_b_from_higgs_mask_last = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0
		
		self.events["anti_b_from_H_last"] = ak.zip(
		    {
		        "pt": self.events.GenPart.pt[anti_b_from_higgs_mask_last],
		        "eta": self.events.GenPart.eta[anti_b_from_higgs_mask_last],
		        "phi": self.events.GenPart.phi[anti_b_from_higgs_mask_last],
		        "mass": self.events.GenPart.mass[anti_b_from_higgs_mask_last]
		    },
		    with_name="PtEtaPhiMCandidate"
		)
		
		anti_b_from_higgs_mask_first = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0
		
		self.events["anti_b_from_H_first"] = ak.zip(
		    {
		        "pt": self.events.GenPart.pt[anti_b_from_higgs_mask_first],
		        "eta": self.events.GenPart.eta[anti_b_from_higgs_mask_first],
		        "phi": self.events.GenPart.phi[anti_b_from_higgs_mask_first],
		        "mass": self.events.GenPart.mass[anti_b_from_higgs_mask_first]
		    },
		    with_name="PtEtaPhiMCandidate"
		)


		# lepton
		lepton_mask = (from_hard_process & is_last_copy & is_lepton & is_from_anti_top) != 0

		self.events["leptons"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[lepton_mask],
				"eta": self.events.GenPart.eta[lepton_mask],
				"phi": self.events.GenPart.phi[lepton_mask],
				"mass": self.events.GenPart.mass[lepton_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)


		# anti-lepton
		anti_lepton_mask = (from_hard_process & is_last_copy & is_anti_lepton & is_from_top) != 0

		self.events["anti_leptons"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[anti_lepton_mask],
				"eta": self.events.GenPart.eta[anti_lepton_mask],
				"phi": self.events.GenPart.phi[anti_lepton_mask],
				"mass": self.events.GenPart.mass[anti_lepton_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)


		# neutrino
		neutrino_mask = (from_hard_process & is_last_copy & is_neutrino & is_from_top) != 0

		self.events["neutrinos"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[neutrino_mask],
				"eta": self.events.GenPart.eta[neutrino_mask],
				"phi": self.events.GenPart.phi[neutrino_mask],
				"mass": self.events.GenPart.mass[neutrino_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)


		# anti-neutrino
		anti_neutrino_mask = (from_hard_process & is_last_copy & is_anti_neutrino & is_from_anti_top) != 0

		self.events["anti_neutrinos"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[anti_neutrino_mask],
				"eta": self.events.GenPart.eta[anti_neutrino_mask],
				"phi": self.events.GenPart.phi[anti_neutrino_mask],
				"mass": self.events.GenPart.mass[anti_neutrino_mask]
			},
			with_name="PtEtaPhiMCandidate"
		)	
		
		####### LORENTZ VECTORS

		#Selecting events with desired particles
		all_particles_present_mask = (
			(ak.num(self.events["b_from_t"]) > 0)
			& (ak.num(self.events["anti_b_from_anti_t"]) > 0)
		    & (ak.num(self.events["b_from_H_last"]) > 0)
		    & (ak.num(self.events["anti_b_from_H_last"]) > 0)
		    & (ak.num(self.events["b_from_H_first"]) > 0)
		    & (ak.num(self.events["anti_b_from_H_first"]) > 0)
		    & (ak.num(self.events["leptons"]) > 0)
			& (ak.num(self.events["anti_leptons"]) > 0)
			& (ak.num(self.events["neutrinos"]) > 0)
		    & (ak.num(self.events["anti_neutrinos"]) > 0)
		)
		


		#TODO: Improve with a loop with the keys
		b_from_top_AllParticles = self.events["b_from_t"][all_particles_present_mask]
		anti_b_from_anti_top_AllParticles = self.events["anti_b_from_anti_t"][all_particles_present_mask]
		#b_from_H_AllParticles = self.events["b_from_H"][all_particles_present_mask]
		#anti_b_from_H_AllParticles = self.events["anti_b_from_H"][all_particles_present_mask]
		
		b_from_H_last_AllParticles = self.events["b_from_H_last"][all_particles_present_mask]
		anti_b_from_H_last_AllParticles = self.events["anti_b_from_H_last"][all_particles_present_mask]
		b_from_H_first_AllParticles = self.events["b_from_H_first"][all_particles_present_mask]
		anti_b_from_H_first_AllParticles = self.events["anti_b_from_H_first"][all_particles_present_mask]
		
		leptons_AllParticles = self.events["leptons"][all_particles_present_mask]
		anti_leptons_AllParticles = self.events["anti_leptons"][all_particles_present_mask]
		neutrinos_AllParticles = self.events["neutrinos"][all_particles_present_mask]
		anti_neutrinos_AllParticles = self.events["anti_neutrinos"][all_particles_present_mask]
		
		### Generate LV
		pt_index_sort = ak.argsort(b_from_top_AllParticles["pt"], ascending=False)
		b_from_t_clean = b_from_top_AllParticles[pt_index_sort][:,0] # select events with all particles, sort by pt and take the higher pt
		b_from_t_LV = ak.to_numpy(b_from_t_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_b_from_anti_top_AllParticles["pt"], ascending=False)
		anti_b_from_anti_t_clean = anti_b_from_anti_top_AllParticles[pt_index_sort][:,0]
		anti_b_from_anti_t_LV = ak.to_numpy(anti_b_from_anti_t_clean).view(vector.MomentumNumpy4D)

		# pt_index_sort = ak.argsort(b_from_H_AllParticles["pt"], ascending=False)
		# b_from_H_clean = b_from_H_AllParticles[pt_index_sort][:,0]
		# b_from_H_LV = ak.to_numpy(b_from_H_clean).view(vector.MomentumNumpy4D)

		# pt_index_sort = ak.argsort(anti_b_from_H_AllParticles["pt"], ascending=False)
		# anti_b_from_H_clean = anti_b_from_H_AllParticles[pt_index_sort][:,0]
		# anti_b_from_H_LV = ak.to_numpy(anti_b_from_H_clean).view(vector.MomentumNumpy4D)
		
		
		pt_index_sort = ak.argsort(b_from_H_last_AllParticles["pt"], ascending=False)
		b_from_H_clean_last = b_from_H_last_AllParticles[pt_index_sort][:,0]
		b_from_H_last_LV = ak.to_numpy(b_from_H_clean_last).view(vector.MomentumNumpy4D)
		
		pt_index_sort = ak.argsort(anti_b_from_H_last_AllParticles["pt"], ascending=False)
		anti_b_from_H_clean_last = anti_b_from_H_last_AllParticles[pt_index_sort][:,0]
		anti_b_from_H_last_LV = ak.to_numpy(anti_b_from_H_clean_last).view(vector.MomentumNumpy4D)
		
		pt_index_sort = ak.argsort(b_from_H_first_AllParticles["pt"], ascending=False)
		b_from_H_clean_first = b_from_H_first_AllParticles[pt_index_sort][:,0]
		b_from_H_first_LV = ak.to_numpy(b_from_H_clean_first).view(vector.MomentumNumpy4D)
		
		pt_index_sort = ak.argsort(anti_b_from_H_first_AllParticles["pt"], ascending=False)
		anti_b_from_H_clean_first = anti_b_from_H_first_AllParticles[pt_index_sort][:,0]
		anti_b_from_H_first_LV = ak.to_numpy(anti_b_from_H_clean_first).view(vector.MomentumNumpy4D)
		

		pt_index_sort = ak.argsort(leptons_AllParticles["pt"], ascending=False)
		leptons_clean = leptons_AllParticles[pt_index_sort][:,0]
		leptons_LV = ak.to_numpy(leptons_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_leptons_AllParticles["pt"], ascending=False)
		anti_leptons_clean = anti_leptons_AllParticles[pt_index_sort][:,0]
		anti_leptons_LV = ak.to_numpy(anti_leptons_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(neutrinos_AllParticles["pt"], ascending=False)
		neutrinos_clean = neutrinos_AllParticles[pt_index_sort][:,0]
		neutrinos_LV = ak.to_numpy(neutrinos_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_neutrinos_AllParticles["pt"], ascending=False)
		anti_neutrinos_clean = anti_neutrinos_AllParticles[pt_index_sort][:,0]
		anti_neutrinos_LV = ak.to_numpy(anti_neutrinos_clean).view(vector.MomentumNumpy4D)

		

		#COMBINE LV to create top and Higgs

		Wp_LV = anti_leptons_LV + neutrinos_LV
		Wn_LV = leptons_LV + anti_neutrinos_LV
		top_LV = Wp_LV + b_from_t_LV
		anti_top_LV = Wn_LV + anti_b_from_anti_t_LV
		higgs_last_LV = b_from_H_last_LV + anti_b_from_H_last_LV
		higgs_first_LV = b_from_H_first_LV + anti_b_from_H_first_LV
		
		higgs_last_pt_fixed = np.where( ( (higgs_first_LV.pt == 0) & (higgs_last_LV.pt == 0) ), 1, higgs_last_LV.pt)
		higgs_last_pt_fixed = np.where( ( (higgs_first_LV.pt == 0) & (higgs_last_LV.pt >= 0) ), np.inf, higgs_last_LV.pt)
		
		higgs_first_pt_fixed = np.where(higgs_first_LV.pt == 0, 1, higgs_first_LV.pt)

		ratio = higgs_last_pt_fixed / higgs_first_pt_fixed
		
		#ratio = higgs_last_LV.pt / higgs_first_LV.pt
		
		
		#Higgs branch
		self.events["Higgs"] = ak.pad_none(ak.zip( 
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
		) , len(self.events), axis=0)
		
	
			
		#Boolean array with particles from asked mother	

		
	def count_higgs(self):
		self.events["nHiggs"] = ak.num(self.events.Higgs, axis=0)
         # use count since we have None
		

		
	def process_extra_after_presel(self, variation) -> ak.Array:
		self.do_gen_selection()
		self.count_higgs()	
