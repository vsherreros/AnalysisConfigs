#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  workflow_mc_copy.py
#  
#  Copyright 2024 Victor Serrano <victorsh@naf-cms16.desy.de>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import sys
import awkward as ak
import numpy as np
import numba

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
		
		
		particle_mass = {"b": 4.18, "electron": 0.0005110, "nu_e": 0.0, "muon": 0.105658, "nu_mu": 0.0}
		mother_idx = self.events.GenPart.genPartIdxMother

		##################  PARTICLE IDENTIFICATION  #####################################

        # Built the mask for particle identification
		from_hard_process = (self.events.GenPart.statusFlags & 256) == 256  
		is_last_copy = (self.events.GenPart.statusFlags & 8192) == 8192  # 4096: first copy / 8192: last copy in FSR
		#is_first_copy = (self.events.GenPart.statusFlags & 4096) == 4096  

		is_b = self.events.GenPart.pdgId == 5
		is_anti_b = self.events.GenPart.pdgId == -5

		is_lepton = (self.events.GenPart.pdgId == 11) | (self.events.GenPart.pdgId == 13)
		is_anti_lepton = (self.events.GenPart.pdgId == -11) | (self.events.GenPart.pdgId == -13)
		is_neutrino = (self.events.GenPart.pdgId == 12) | (self.events.GenPart.pdgId == 14)
		is_anti_neutrino = (self.events.GenPart.pdgId == -12) | (self.events.GenPart.pdgId == -14)
		
		is_b_jet = abs(self.events.GenJet.partonFlavour) == 5 	
		
		# comes from top or anti-top
		is_from_top = find_mother(self, 6)
		is_from_anti_top = find_mother(self, -6)
		not_from_top = (is_from_top | is_from_anti_top) == 0
	


		#### GENERATE BRANCHES FOR IDENTIFIED PARTICLES

		# b quark from top
		b_from_top_mask = (from_hard_process & is_last_copy & is_b & is_from_top) != 0
		
		self.events["b_from_t"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[b_from_top_mask],
				"eta": self.events.GenPart.eta[b_from_top_mask],
				"phi": self.events.GenPart.phi[b_from_top_mask],
				"mass": ak.full_like(self.events.GenPart.mass[b_from_top_mask], particle_mass["b"])
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
				"mass": ak.full_like(self.events.GenPart.mass[anti_b_from_anti_top_mask], particle_mass["b"])
			},
			with_name="PtEtaPhiMCandidate"
		)


		#b quark from Higgs
		b_from_higgs_mask = (from_hard_process & is_last_copy & is_b & not_from_top) != 0

		self.events["b_from_H"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[b_from_higgs_mask],
				"eta": self.events.GenPart.eta[b_from_higgs_mask],
				"phi": self.events.GenPart.phi[b_from_higgs_mask],
				"mass": ak.full_like(self.events.GenPart.mass[b_from_higgs_mask], particle_mass["b"])
			},
			with_name="PtEtaPhiMCandidate"
		)


		# anti-b quark from Higgs
		anti_b_from_higgs_mask = (from_hard_process & is_last_copy & is_anti_b & not_from_top) != 0

		self.events["anti_b_from_H"] = ak.zip(
			{
				"pt": self.events.GenPart.pt[anti_b_from_higgs_mask],
				"eta": self.events.GenPart.eta[anti_b_from_higgs_mask],
				"phi": self.events.GenPart.phi[anti_b_from_higgs_mask],
				"mass": ak.full_like(self.events.GenPart.mass[anti_b_from_higgs_mask], particle_mass["b"])
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
		
		#b-jets
		b_jets = ak.zip(
			{
				"pt": self.events.GenJet.pt[is_b_jet],
				"eta": self.events.GenJet.eta[is_b_jet],
				"phi": self.events.GenJet.phi[is_b_jet],
				"mass": self.events.GenJet.mass[is_b_jet]
			},
			with_name="PtEtaPhiMCandidate"
		)
		
		####### LORENTZ VECTORS

		#Selecting events with desired particles
		all_particles_present_mask = (
			(ak.num(self.events["b_from_t"]) > 0)
			& (ak.num(self.events["anti_b_from_anti_t"]) > 0)
		    & (ak.num(self.events["b_from_H"]) > 0)
		    & (ak.num(self.events["anti_b_from_H"]) > 0)
		    & (ak.num(self.events["leptons"]) > 0)
			& (ak.num(self.events["anti_leptons"]) > 0)
			& (ak.num(self.events["neutrinos"]) > 0)
		    & (ak.num(self.events["anti_neutrinos"]) > 0)
		    #    & (ak.num(b_jets) > 3)
		)
		
		b_jets_AllParticles = b_jets[all_particles_present_mask]
		b_jets_LV = vector.zip( {"pt": b_jets_AllParticles['pt'], "eta": b_jets_AllParticles['eta'], "phi": b_jets_AllParticles['phi'], "mass": b_jets_AllParticles['mass']})


		#TODO: Improve with a loop with the keys ??
		b_from_top_AllParticles = self.events["b_from_t"][all_particles_present_mask]
		anti_b_from_anti_top_AllParticles = self.events["anti_b_from_anti_t"][all_particles_present_mask]
		b_from_H_AllParticles = self.events["b_from_H"][all_particles_present_mask]
		anti_b_from_H_AllParticles = self.events["anti_b_from_H"][all_particles_present_mask]
		
		leptons_AllParticles = self.events["leptons"][all_particles_present_mask]
		anti_leptons_AllParticles = self.events["anti_leptons"][all_particles_present_mask]
		neutrinos_AllParticles = self.events["neutrinos"][all_particles_present_mask]
		anti_neutrinos_AllParticles = self.events["anti_neutrinos"][all_particles_present_mask]
		
		### Generate LV for b quarks and selecting the one with highest pt
		pt_index_sort = ak.argsort(b_from_top_AllParticles["pt"], ascending=False)
		b_from_t_clean = b_from_top_AllParticles[pt_index_sort][:,0] # select events with all particles, sort by pt and take the higher pt
		b_from_t_LV = ak.to_numpy(b_from_t_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_b_from_anti_top_AllParticles["pt"], ascending=False)
		anti_b_from_anti_t_clean = anti_b_from_anti_top_AllParticles[pt_index_sort][:,0]
		anti_b_from_anti_t_LV = ak.to_numpy(anti_b_from_anti_t_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(b_from_H_AllParticles["pt"], ascending=False)
		b_from_H_clean = b_from_H_AllParticles[pt_index_sort][:,0]
		b_from_H_LV = ak.to_numpy(b_from_H_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_b_from_H_AllParticles["pt"], ascending=False)
		anti_b_from_H_clean = anti_b_from_H_AllParticles[pt_index_sort][:,0]
		anti_b_from_H_LV = ak.to_numpy(anti_b_from_H_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(leptons_AllParticles["pt"], ascending=False)
		leptons_clean = leptons_AllParticles[pt_index_sort][:,0]
		#leptons_LV = ak.to_numpy(leptons_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_leptons_AllParticles["pt"], ascending=False)
		anti_leptons_clean = anti_leptons_AllParticles[pt_index_sort][:,0]
		#anti_leptons_LV = ak.to_numpy(anti_leptons_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(neutrinos_AllParticles["pt"], ascending=False)
		neutrinos_clean = neutrinos_AllParticles[pt_index_sort][:,0]
		#neutrinos_LV = ak.to_numpy(neutrinos_clean).view(vector.MomentumNumpy4D)

		pt_index_sort = ak.argsort(anti_neutrinos_AllParticles["pt"], ascending=False)
		anti_neutrinos_clean = anti_neutrinos_AllParticles[pt_index_sort][:,0]
		#anti_neutrinos_LV = ak.to_numpy(anti_neutrinos_clean).view(vector.MomentumNumpy4D)


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
		
		b_top_jet = b_top_jet[all_b_jets_mask]
		anti_b_top_jet = anti_b_top_jet[all_b_jets_mask]
		b_H_jet = b_H_jet[all_b_jets_mask]
		anti_b_H_jet = anti_b_H_jet[all_b_jets_mask]
		
		leptons_clean = leptons_clean[all_b_jets_mask]
		anti_leptons_clean = anti_leptons_clean[all_b_jets_mask]
		neutrinos_clean = neutrinos_clean[all_b_jets_mask]
		anti_neutrinos_clean = anti_neutrinos_clean[all_b_jets_mask]
		
		#Lorentz vector for leptons
		leptons_LV = ak.to_numpy(leptons_clean).view(vector.MomentumNumpy4D)
		anti_leptons_LV = ak.to_numpy(anti_leptons_clean).view(vector.MomentumNumpy4D)
		neutrinos_LV = ak.to_numpy(neutrinos_clean).view(vector.MomentumNumpy4D)
		anti_neutrinos_LV = ak.to_numpy(anti_neutrinos_clean).view(vector.MomentumNumpy4D)
		
		#COMBINE LV to create top and Higgs
		Wp_LV = anti_leptons_LV + neutrinos_LV
		Wn_LV = leptons_LV + anti_neutrinos_LV
		top_LV = Wp_LV + b_top_jet
		anti_top_LV = Wn_LV + anti_b_top_jet
		H_jets = b_H_jet + anti_b_H_jet
		
		
		#################  SPIN VARIABLES  ############################
		
		#cos(theta) - H/Z/g system
		b_ZMF = b_H_jet.boostCM_of(H_jets)
		cos_theta_b_h = H_jets.to_beta3().unit().dot(b_ZMF.to_beta3().unit())
		
		
		#Boosting vectors
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
		
		
		### Spin Observables for tt_bar ###
		
		LLbar_DeltaPhi = abs(leptons_LV.deltaphi(anti_leptons_LV)) ### In LAB system not in tt ZMF
		LLbar_DeltaEta = abs(leptons_LV.deltaeta(anti_leptons_LV)) ### In LAB system not in tt ZMF
		c_hel = L.dot(L_bar)
		c_han = ( L.dot(kvect_p) * L_bar.dot(kvect_p) ) - ( L.dot(rvect_p) * L_bar.dot(rvect_p) ) - ( L.dot(nvect_p) *  L_bar.dot(nvect_p) )
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
	
	
		#min DeltaPhi in LAB

		#calculate Phi between two particles
		TTbar_phi = abs( top_LV.deltaphi(anti_top_LV) )
		THiggs_phi = abs( top_LV.deltaphi(H_jets) )
		TbarHiggs_phi = abs( anti_top_LV.deltaphi(H_jets) )
		
		TTbar_dR = abs( top_LV.deltaR(anti_top_LV) ) ##
		THiggs_dR = abs( top_LV.deltaR(H_jets) ) ##
		TbarHiggs_dR = abs( anti_top_LV.deltaR(H_jets) ) ##
		
		Phi_combined = np.concatenate((TTbar_phi , THiggs_phi , TbarHiggs_phi ), axis=1)
		min_phi = np.min(Phi_combined, axis=1)
		
		dR_combined = np.concatenate((TTbar_dR, THiggs_dR, TbarHiggs_dR), axis=1) 
		min_dR = np.min(dR_combined, axis=1) ##
		
		####   min DeltaPhi in tt_ZMF   #####
		# Boosting to ttH frame

		#top_ZMF 
		#anti_top_ZMF  
		# H_ZMF = H_jets.boostCM_of_p4(tt_LV)
		
		# TTbar_phi_tt_ZMF = abs( top_ZMF.deltaphi(anti_top_ZMF) )
		# THiggs_phi_tt_ZMF = abs( top_ZMF.deltaphi(H_ZMF) )
		# TbarHiggs_phi_tt_ZMF = abs( anti_top_ZMF.deltaphi(H_ZMF) )
		
		# TTbar_eta_tt_ZMF = abs( top_ZMF.deltaeta(anti_top_ZMF) )
		# THiggs_eta_tt_ZMF = abs( top_ZMF.deltaeta(H_ZMF) )
		# TbarHiggs_eta_tt_ZMF = abs( anti_top_ZMF.deltaeta(H_ZMF) )
		
		# Phi_combined_tt_ZMF = np.concatenate((TTbar_phi_tt_ZMF , THiggs_phi_tt_ZMF , TbarHiggs_phi_tt_ZMF ), axis=1)
		# min_phi_tt_ZMF = np.min(Phi_combined_tt_ZMF, axis=1)
		
		# Eta_combined_tt_ZMF = np.concatenate((TTbar_eta_tt_ZMF, THiggs_eta_tt_ZMF, TbarHiggs_eta_tt_ZMF), axis=1)
		# min_eta_tt_ZMF = np.min(Eta_combined_tt_ZMF, axis=1)

		
		
		####   min DeltaPhi in ttH_ZMF   #####
		
		# Boosting to ttH frame
		ttH_ZMF = top_LV + anti_top_LV + H_jets
		
		top_ttH_ZMF = top_LV.boostCM_of_p4(ttH_ZMF)
		anti_top_ttH_ZMF = anti_top_LV.boostCM_of_p4(ttH_ZMF)
		H_ttH_ZMF = H_jets.boostCM_of_p4(ttH_ZMF)
		
		TTbar_phi_ttH_ZMF = abs( top_ttH_ZMF.deltaphi(anti_top_ttH_ZMF) )
		THiggs_phi_ttH_ZMF = abs( top_ttH_ZMF.deltaphi(H_ttH_ZMF) )
		TbarHiggs_phi_ttH_ZMF = abs( anti_top_ttH_ZMF.deltaphi(H_ttH_ZMF) )
		
		TTbar_dR_ttH_ZMF = abs( top_ttH_ZMF.deltaR(anti_top_ttH_ZMF) )
		THiggs_dR_ttH_ZMF = abs( top_ttH_ZMF.deltaR(H_ttH_ZMF) )
		TbarHiggs_dR_ttH_ZMF = abs( anti_top_ttH_ZMF.deltaR(H_ttH_ZMF) )
		
		Phi_combined_ttH_ZMF = np.concatenate((TTbar_phi_ttH_ZMF , THiggs_phi_ttH_ZMF , TbarHiggs_phi_ttH_ZMF ), axis=1)
		min_phi_ttH_ZMF = np.min(Phi_combined_ttH_ZMF, axis=1)
		
		dR_combined_ttH_ZMF = np.concatenate((TTbar_dR_ttH_ZMF, THiggs_dR_ttH_ZMF, TbarHiggs_dR_ttH_ZMF), axis=1)
		min_dR_ttH_ZMF = np.min(dR_combined_ttH_ZMF, axis=1)

		
		
		
		self.events["Spin"] = ak.pad_none(ak.zip( 
			{
				# "cos_theta": cos_theta_b_h,
				# "c_hel": c_hel,
				"c_han": c_han,
				# "LLbar_DeltaPhi": LLbar_DeltaPhi ,
				# "LLbar_DeltaEta": LLbar_DeltaEta,
				# "b1k_plus_b2k": b1k_plus_b2k,
		        # "b1r_plus_b2r": b1r_plus_b2r,
		        # "b1n_plus_b2n": b1n_plus_b2n,
		        # "b1k_minus_b2k": b1k_minus_b2k,
		        # "b1r_minus_b2r": b1r_minus_b2r,
		        # "b1n_minus_b2n": b1n_minus_b2n,
		        # "Ckk": Ckk,
		        # "Crr": Crr,
		        # "Cnn": Cnn,
		        # "Crk_plus_Ckr": Crk_plus_Ckr,
		        # "Crn_plus_Cnr": Crn_plus_Cnr,
		        # "Ckn_plus_Cnk": Ckn_plus_Cnk,
		        # "Crk_minus_Ckr": Crk_minus_Ckr,
		        # "Crn_minus_Cnr": Crn_minus_Cnr,
		        # "Ckn_minus_Cnk": Ckn_minus_Cnk,

			},
			with_name="SpinCandidate"
		) , len(self.events), axis=0)
		
		
		self.events["Angular"] = ak.pad_none(ak.zip( 
			{
				"min_phi": min_phi,
		        "min_phi_ttH_ZMF": min_phi_ttH_ZMF,
		        "TTbar_phi": TTbar_phi,
		        "TTbar_phi_ttH_ZMF": TTbar_phi_ttH_ZMF,
		        "THiggs_phi": THiggs_phi,
		        "THiggs_phi_ttH_ZMF": THiggs_phi_ttH_ZMF,
		        "TbarHiggs_phi": TbarHiggs_phi,
		        "TbarHiggs_phi_ttH_ZMF": TbarHiggs_phi_ttH_ZMF,

				"min_dR" : min_dR, 
				"min_dR_ttH_ZMF" : min_dR_ttH_ZMF,
				"TTbar_dR" : TTbar_dR,
				"TTbar_dR_ttH_ZMF" : TTbar_dR_ttH_ZMF,
				"THiggs_dR" : THiggs_dR,
				"THiggs_dR_ttH_ZMF" : THiggs_dR_ttH_ZMF,
				"TbarHiggs_dR" : TbarHiggs_dR,
				"TbarHiggs_dR_ttH_ZMF" : TbarHiggs_dR_ttH_ZMF,
				
			},
			with_name="AngularCandidate"
		) , len(self.events), axis=0)
		
		
		#Higgs branch
#		self.events["Higgs"] = ak.pad_none(ak.zip( 
#			{
#				"pt": higgs_LV.pt,
#				"eta": higgs_LV.eta,
#				"phi": higgs_LV.phi,
#				"mass": higgs_LV.mass,
#			},
#			with_name="PtEtaPhiMCandidate"
#		) , len(self.events), axis=0)
		
		
		#Higgs branch
		self.events["Higgs"] = ak.pad_none(ak.zip( 
			{
				"pt": H_jets.pt,
				"eta": H_jets.eta,
				"phi": H_jets.phi,
				"mass": H_jets.mass,
			},
			with_name="PtEtaPhiMCandidate"
		) , len(self.events), axis=0)
		

		
	def count_higgs(self):
		self.events["nHiggs"] = ak.num(self.events.Higgs, axis=0)
         # use count since we have None
		

		
	def process_extra_after_presel(self, variation) -> ak.Array:
		self.do_gen_selection()
		self.count_higgs()	
