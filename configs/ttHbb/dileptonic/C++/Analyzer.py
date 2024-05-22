#!/usr/bin/env python3
 
import ROOT
import numpy as np
import copy
from argparse import ArgumentParser

parser = ArgumentParser(description="Select events with Gen-Lev informatino to reconstruct obervables")
parser.add_argument("--proc", help="Signal or Background. Default is Background")
parser.add_argument("--sig_med", help="For signal, specify type of mediator")
parser.add_argument("--sig_mass", help="For signal, specify mass")
args = parser.parse_args()


nGenPart = np.empty((1), dtype="int32")
GenPart_pdgId = np.empty((1000), dtype="int32")
GenPart_genPartIdxMother = np.empty((1000), dtype="int32")
GenPart_status = np.empty((1000), dtype="int32")
GenPart_statusFlags = np.empty((1000), dtype="int32")
GenPart_mass = np.empty((1000), dtype="float32")
GenPart_eta = np.empty((1000), dtype="float32")
GenPart_phi = np.empty((1000), dtype="float32")
GenPart_pt = np.empty((1000), dtype="float32")
GenMET_pt = np.empty((1), dtype="float32")
GenMET_phi = np.empty((1), dtype="float32")
TTbarDMJets_Scalar_Mphi_10 = np.empty((1), dtype="bool")
TTbarDMJets_Pseudoscalar_Mphi_10 = np.empty((1), dtype="bool")
TTbarDMJets_Scalar_Mphi_500 = np.empty((1), dtype="bool")
TTbarDMJets_Pseudoscalar_Mphi_500 = np.empty((1), dtype="bool")


CHEL = np.empty((1), dtype="float32")
LLbar_DeltaPhi = np.empty((1), dtype="float32")
b1k_plus_b2k = np.empty((1), dtype="float32")
b1r_plus_b2r = np.empty((1), dtype="float32")
b1n_plus_b2n = np.empty((1), dtype="float32")
b1k_minus_b2k = np.empty((1), dtype="float32")
b1r_minus_b2r = np.empty((1), dtype="float32")
b1n_minus_b2n = np.empty((1), dtype="float32")
Ckk = np.empty((1), dtype="float32")
Crr = np.empty((1), dtype="float32")
Cnn = np.empty((1), dtype="float32")
Crk_plus_Ckr = np.empty((1), dtype="float32")
Crn_plus_Cnr = np.empty((1), dtype="float32")
Ckn_plus_Cnk = np.empty((1), dtype="float32")
Crk_minus_Ckr = np.empty((1), dtype="float32")
Crn_minus_Cnr = np.empty((1), dtype="float32")
Ckn_minus_Cnk = np.empty((1), dtype="float32")

pdgId_dic = {6 : 't', -6 : 't_bar', 5 : 'b', -5 : 'b_bar', 24 : 'w_plus', -24 : 'w_minus', 11 : 'l', -11 : 'l_bar', 13 : 'l', -13 : 'l_bar', 12 : 'v', -12 : 'v_bar', 14 : 'v', -14 : 'v_bar'}

filein = 'None'
fileout = 'None'
massbool = False

print("\n")

if args.proc is None or args.proc == 'Background':
	print("Processing 'Background' ... \n")
	filein = 'TTbar_Dilepton'
	fileout = 'Background_TTree_Observables'
else:
	if args.sig_med is None or args.sig_med=='Scalar':
		print("Processing 'Signal' 'Scalar'... \n")
		filein = 'TTbarDMJets_Dilepton_Scalar'
		if args.sig_mass is None:
			print("Processing M=10 ... \n")
			fileout = 'Signal_Scalar_M10_TTree_Observables'
			massbool = TTbarDMJets_Scalar_Mphi_10
		else:
			print(f'Processing M={args.sig_mass} ...\n')
			fileout = 'Signal_Scalar_M'+args.sig_mass+'_TTree_Observables'
			if args.sig_mass == '10':
				massbool = TTbarDMJets_Scalar_Mphi_10
			elif args.sig_mass == '500':
				massbool = TTbarDMJets_Scalar_Mphi_500
	else:
		print("Processing 'Signal' 'Pseudoscalar'... \n")
		filein = 'TTbarDMJets_Dilepton_Pseudoscalar'
		if args.sig_mass is None:
			print("Processing M=10 ... \n")
			fileout = 'Signal_Pseudoscalar_M10_TTree_Observables'
			massbool = TTbarDMJets_Pseudoscalar_Mphi_10
		else:
			print(f'Processing M={args.sig_mass} ... \n')
			fileout = 'Signal_Pseudoscalar_M'+args.sig_mass+'_TTree_Observables'
			if args.sig_mass == '10':
				massbool = TTbarDMJets_Pseudoscalar_Mphi_10
			elif args.sig_mass == '500':
				massbool = TTbarDMJets_Pseudoscalar_Mphi_500		
				

file_in = ROOT.TFile(filein+'.root', 'read')
file_out = ROOT.TFile(fileout+'.root', 'recreate')

tree = file_in.Get('Events')
tree.SetBranchAddress('nGenPart', nGenPart)
tree.SetBranchAddress('GenPart_pdgId', GenPart_pdgId)
tree.SetBranchAddress('GenPart_genPartIdxMother', GenPart_genPartIdxMother)
tree.SetBranchAddress('GenPart_status', GenPart_status)
tree.SetBranchAddress('GenPart_statusFlags', GenPart_statusFlags)
tree.SetBranchAddress('GenPart_mass', GenPart_mass)
tree.SetBranchAddress('GenPart_eta', GenPart_eta)
tree.SetBranchAddress('GenPart_phi', GenPart_phi)
tree.SetBranchAddress('GenPart_pt', GenPart_pt)
tree.SetBranchAddress('MET_pt', GenMET_pt)
tree.SetBranchAddress('MET_phi', GenMET_phi)
if args.proc == 'Signal' and args.sig_med == 'Scalar':
	tree.SetBranchAddress('GenModel__TTbarDMJets_Dilepton_scalar_LO_Mchi_1_Mphi_50_TuneCP5_13TeV_madgraph_mcatnlo_pythia8', TTbarDMJets_Scalar_Mphi_10)
	tree.SetBranchAddress('GenModel__TTbarDMJets_Dilepton_scalar_LO_Mchi_1_Mphi_500_TuneCP5_13TeV_madgraph_mcatnlo_pythia8', TTbarDMJets_Scalar_Mphi_500)
if args.proc == 'Signal' and args.sig_med == 'Pseudoscalar':
	tree.SetBranchAddress('GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_50_TuneCP5_13TeV_madgraph_mcatnlo_pythia8', TTbarDMJets_Pseudoscalar_Mphi_10)
	tree.SetBranchAddress('GenModel__TTbarDMJets_Dilepton_pseudoscalar_LO_Mchi_1_Mphi_500_TuneCP5_13TeV_madgraph_mcatnlo_pythia8', TTbarDMJets_Pseudoscalar_Mphi_500)

Tree_Obs = ROOT.TTree('Tree_Obs','Tree_Obs')
Tree_Obs.Branch('CHEL',CHEL,'CHEL')
Tree_Obs.Branch('LLbar_DeltaPhi',LLbar_DeltaPhi,'LLbar_DeltaPhi')
Tree_Obs.Branch('b1k_plus_b2k',b1k_plus_b2k,'b1k_plus_b2k')
Tree_Obs.Branch('b1r_plus_b2r',b1r_plus_b2r,'b1r_plus_b2r')
Tree_Obs.Branch('b1n_plus_b2n',b1n_plus_b2n,'b1n_plus_b2n')
Tree_Obs.Branch('b1k_minus_b2k',b1k_minus_b2k,'b1k_minus_b2k')
Tree_Obs.Branch('b1r_minus_b2r',b1r_minus_b2r,'b1r_minus_b2r')
Tree_Obs.Branch('b1n_minus_b2n',b1n_minus_b2n,'b1n_minus_b2n')
Tree_Obs.Branch('Ckk',Ckk,'Ckk')
Tree_Obs.Branch('Crr',Crr,'Crr')
Tree_Obs.Branch('Cnn',Cnn,'Cnn')
Tree_Obs.Branch('Crk_plus_Ckr',Crk_plus_Ckr,'Crk_plus_Ckr')
Tree_Obs.Branch('Crn_plus_Cnr',Crn_plus_Cnr,'Crn_plus_Cnr')
Tree_Obs.Branch('Ckn_plus_Cnk',Ckn_plus_Cnk,'Ckn_plus_Cnk')
Tree_Obs.Branch('Crk_minus_Ckr',Crk_minus_Ckr,'Crk_minus_Ckr')
Tree_Obs.Branch('Crn_minus_Cnr',Crn_minus_Cnr,'Crn_minus_Cnr')
Tree_Obs.Branch('Ckn_minus_Cnk',Ckn_minus_Cnk,'Ckn_minus_Cnk')

print(f'Number of events in the total sample: {tree.GetEntries()} \n')

counter = 0
for iEvent in range(100):
	tree.GetEntry(iEvent)
	
	if args.proc=='Signal':
		#selecting mass point
		if not massbool:
			continue
	
	iPart_dic = {}
	
	#find particles
	for iPart in range(nGenPart[0]):
		
		if(len(iPart_dic.keys()) >= 10):
			break
		else:
			if not all(x in iPart_dic.keys() for x in ['t', 't_bar']):
				if(abs(GenPart_pdgId[iPart])==6 and GenPart_statusFlags[iPart]==10497):
					iPart_dic[pdgId_dic[GenPart_pdgId[iPart]]] = iPart
			if not all(x in iPart_dic.keys() for x in ['w_plus', 'w_minus']):
				if(abs(GenPart_pdgId[iPart])==24 and abs(GenPart_pdgId[GenPart_genPartIdxMother[iPart]])==6):
					iPart_dic[pdgId_dic[GenPart_pdgId[iPart]]] = iPart
			if not all(x in iPart_dic.keys() for x in ['b', 'b_bar']):
				if(abs(GenPart_pdgId[iPart])==5 and abs(GenPart_pdgId[GenPart_genPartIdxMother[iPart]])==6):
					iPart_dic[pdgId_dic[GenPart_pdgId[iPart]]] = iPart
			if not all(x in iPart_dic.keys() for x in ['l', 'l_bar']):
				if((abs(GenPart_pdgId[iPart])==11 or abs(GenPart_pdgId[iPart])==13) and abs(GenPart_pdgId[GenPart_genPartIdxMother[iPart]])==24):
					iPart_dic[pdgId_dic[GenPart_pdgId[iPart]]] = iPart
			if not all(x in iPart_dic.keys() for x in ['v', 'v_bar']):
				if((abs(GenPart_pdgId[iPart])==12 or abs(GenPart_pdgId[iPart])==14) and abs(GenPart_pdgId[GenPart_genPartIdxMother[iPart]])==24):
					iPart_dic[pdgId_dic[GenPart_pdgId[iPart]]] = iPart
			
	if(len(iPart_dic.keys()) is not 10):
		continue
	
	#construct TLorenzVector
	T = ROOT.TLorentzVector()
	T.SetPtEtaPhiM(GenPart_pt[iPart_dic['t']],GenPart_eta[iPart_dic['t']],GenPart_phi[iPart_dic['t']],GenPart_mass[iPart_dic['t']])
	T_bar = ROOT.TLorentzVector()
	T_bar.SetPtEtaPhiM(GenPart_pt[iPart_dic['t_bar']],GenPart_eta[iPart_dic['t_bar']],GenPart_phi[iPart_dic['t_bar']],GenPart_mass[iPart_dic['t_bar']])
	W_plus = ROOT.TLorentzVector()
	W_plus.SetPtEtaPhiM(GenPart_pt[iPart_dic['w_plus']],GenPart_eta[iPart_dic['w_plus']],GenPart_phi[iPart_dic['w_plus']],GenPart_mass[iPart_dic['w_plus']])
	W_minus = ROOT.TLorentzVector()
	W_minus.SetPtEtaPhiM(GenPart_pt[iPart_dic['w_minus']],GenPart_eta[iPart_dic['w_minus']],GenPart_phi[iPart_dic['w_minus']],GenPart_mass[iPart_dic['w_minus']])
	B = ROOT.TLorentzVector()
	B.SetPtEtaPhiM(GenPart_pt[iPart_dic['b']],GenPart_eta[iPart_dic['b']],GenPart_phi[iPart_dic['b']],GenPart_mass[iPart_dic['b']])
	B_bar = ROOT.TLorentzVector()
	B_bar.SetPtEtaPhiM(GenPart_pt[iPart_dic['b_bar']],GenPart_eta[iPart_dic['b_bar']],GenPart_phi[iPart_dic['b_bar']],GenPart_mass[iPart_dic['b_bar']])
	L = ROOT.TLorentzVector()
	L.SetPtEtaPhiM(GenPart_pt[iPart_dic['l']],GenPart_eta[iPart_dic['l']],GenPart_phi[iPart_dic['l']],GenPart_mass[iPart_dic['l']])
	L_bar = ROOT.TLorentzVector()
	L_bar.SetPtEtaPhiM(GenPart_pt[iPart_dic['l_bar']],GenPart_eta[iPart_dic['l_bar']],GenPart_phi[iPart_dic['l_bar']],GenPart_mass[iPart_dic['l_bar']])
	Nu = ROOT.TLorentzVector()
	Nu.SetPtEtaPhiM(GenPart_pt[iPart_dic['v']],GenPart_eta[iPart_dic['v']],GenPart_phi[iPart_dic['v']],GenPart_mass[iPart_dic['v']])
	Nu_bar = ROOT.TLorentzVector()
	Nu_bar.SetPtEtaPhiM(GenPart_pt[iPart_dic['v_bar']],GenPart_eta[iPart_dic['v_bar']],GenPart_phi[iPart_dic['v_bar']],GenPart_mass[iPart_dic['v_bar']])
	Gen_MET = ROOT.TLorentzVector()
	Gen_MET.SetPtEtaPhiM(GenMET_pt, 0., GenMET_phi,0.)
	
	
	#obersvables
	# -------------------------------------------------- #
	RecoTop = T #B + L_bar + Nu
	RecoAntiTop = T_bar #B_bar + L + Nu_bar
	# -------------------------------------------------- #
	TTbar_boostv3 = -(RecoTop+RecoAntiTop).BoostVector()
	RecoTop_ZMF = copy.deepcopy(RecoTop)
	RecoAntiTop_ZMF = copy.deepcopy(RecoAntiTop)
	L_ZMF_AntiTop = copy.deepcopy(L)
	L_bar_ZMF_Top = copy.deepcopy(L_bar)
	
	L_ZMF_AntiTop.Boost(TTbar_boostv3)
	L_bar_ZMF_Top.Boost(TTbar_boostv3)
	RecoTop_ZMF.Boost(TTbar_boostv3)
	RecoAntiTop_ZMF.Boost(TTbar_boostv3)
	
	L_ZMF_AntiTop.Boost(-RecoAntiTop_ZMF.BoostVector())
	L_bar_ZMF_Top.Boost(-RecoTop_ZMF.BoostVector())
	
	L1 = L_ZMF_AntiTop.Vect().Unit()
	L2 = L_bar_ZMF_Top.Vect().Unit() 
	# ------------------------------------------------- #
	pvect_p = ROOT.TVector3(0.,0.,1.)
	kvect_p = RecoTop_ZMF.Vect().Unit()
	y_p = pvect_p.Unit().Dot(kvect_p)
	r_p = ROOT.TMath.Sqrt(1- y_p*y_p)
	rvect_p = np.sign(y_p) * 1/r_p * (pvect_p - y_p * kvect_p)
	nvect_p = np.sign(y_p) * 1/r_p * (pvect_p.Cross(kvect_p))
	
	
	CHEL[0] = L1.Dot(L2)
	LLbar_DeltaPhi[0] = abs(L.DeltaPhi(L_bar))
	b1k_plus_b2k[0] = L1.Dot(kvect_p) +  L2.Dot(kvect_p)
	b1r_plus_b2r[0] = L1.Dot(rvect_p) +  L2.Dot(rvect_p)
	b1n_plus_b2n[0] = L1.Dot(nvect_p) +  L2.Dot(nvect_p)
	b1k_minus_b2k[0] = L1.Dot(kvect_p) -  L2.Dot(kvect_p)
	b1r_minus_b2r[0] = L1.Dot(rvect_p) -  L2.Dot(rvect_p)
	b1n_minus_b2n[0] = L1.Dot(nvect_p) -  L2.Dot(nvect_p) 
	Ckk[0] = L1.Dot(kvect_p) * L2.Dot(kvect_p)
	Crr[0] = L1.Dot(rvect_p) * L2.Dot(rvect_p)
	Cnn[0] = L1.Dot(nvect_p) * L2.Dot(nvect_p)
	Crk_plus_Ckr[0] = L1.Dot(rvect_p) * L2.Dot(kvect_p) + L1.Dot(kvect_p) * L2.Dot(rvect_p)
	Crn_plus_Cnr[0] = L1.Dot(rvect_p) * L2.Dot(nvect_p) + L1.Dot(nvect_p) * L2.Dot(rvect_p)
	Ckn_plus_Cnk[0] = L1.Dot(kvect_p) * L2.Dot(nvect_p) + L1.Dot(nvect_p) * L2.Dot(kvect_p)
	Crk_minus_Ckr[0] = L1.Dot(rvect_p) * L2.Dot(kvect_p) - L1.Dot(kvect_p) * L2.Dot(rvect_p)
	Crn_minus_Cnr[0] = L1.Dot(rvect_p) * L2.Dot(nvect_p) - L1.Dot(nvect_p) * L2.Dot(rvect_p)
	Ckn_minus_Cnk[0] = L1.Dot(kvect_p) * L2.Dot(nvect_p) - L1.Dot(nvect_p) * L2.Dot(kvect_p)

	
	Tree_Obs.Fill()
	
	counter +=1

file_out.cd()
Tree_Obs.Write()
file_out.Close()
		
print(f'Number of events selected: {counter}')
print('Done')

