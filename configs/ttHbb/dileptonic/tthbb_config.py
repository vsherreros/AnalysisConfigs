from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_nObj_eq, get_nBtagEq, get_nBtagMin, get_HLTsel
from pocket_coffea.parameters.histograms import *
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.columns_manager import ColOut

import cloudpickle
import workflow_spin
from workflow_spin import MCProcessor
cloudpickle.register_pickle_by_value(workflow_spin)


# Local imports of functions
from preselection_cuts import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

# adding object preselection
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection_dileptonic.yaml",
                                                  # f"{localdir}/params/btagsf_calibration.yaml",
                                                  #f"{localdir}/params/triggers.yaml",
                                                  # f"{localdir}/params/overrides.yaml",
                                                  update=True)


cfg = Configurator(
	#Parameters
    parameters = parameters,

    #Datasets
    datasets = {
        "jsons": [f"{localdir}/datasets/ttHTobb_ttTo2L2Nu.json",
				  f"{localdir}/datasets/ttbb_4f_TTTo2L2Nu.json",
				  f"{localdir}/datasets/ttZToQQ.json"
                 ],
                 
        "filter" : {
            "samples": ["ttHTobb_ttTo2L2Nu",
						"ttbb_4f_TTTo2L2Nu",
						"ttZToQQ"],
            "samples_exclude" : [],
            "year": ['2018']
        },
        "subsamples":{
        }
    },

    #Workflow
    #workflow = ttHbbBaseProcessor,
	workflow = MCProcessor,

    #Cuts and Categories
    skim = [passthrough],

    preselections = [passthrough],

    categories = {
        "baseline": [passthrough]
    },


    weights = {
        "common": {
             "inclusive": [#"genWeight","lumi","XS",
                          # "pileup",
                          # "sf_ele_reco", "sf_ele_id",
                          # "sf_mu_id","sf_mu_iso",
                          # "sf_btag", "sf_jet_puId", 
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },

    variations = {
        "weights": {
            "common": {
                "inclusive": [  #"pileup",
                                #"sf_ele_reco", "sf_ele_id",
                                #"sf_mu_id", "sf_mu_iso", "sf_jet_puId",
                                #"sf_btag"                               
                              ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
    },

    
   variables = {
#        **ele_hists(coll="ElectronGood", pos=0),
#        **muon_hists(coll="MuonGood", pos=0),
#        **count_hist(name="nElectronGood", coll="ElectronGood",bins=4, start=0, stop=4),
#        **count_hist(name="nMuonGood", coll="MuonGood",bins=4, start=0, stop=4),
#        **count_hist(name="nJets", coll="JetGood",bins=10, start=4, stop=14),
#        **count_hist(name="nBJets", coll="BJetGood",bins=12, start=2, stop=14),
#        **jet_hists(coll="JetGood", pos=0),
#        **jet_hists(coll="JetGood", pos=1),
        "m_bb" : HistConf( [Axis( coll="Higgs", field="mass", bins=30, start=0, stop=200, label=r"$m_{b\bar{b}}$ [GeV]")]),
		"pt_bb" : HistConf( [Axis( coll="Higgs", field="pt", bins=30, start=0, stop=400, label=r"$p_{T,b\bar{b}}$ [GeV]")]),
#		"eta_bb" : HistConf( [Axis( coll="Higgs", field="eta", bins=30, start=-5, stop=5, label=r"$\eta_{b\bar{b}}$ [GeV]")]),
#		"phi_bb" : HistConf( [Axis( coll="Higgs", field="phi", bins=30, start=-3.1416, stop=3.1416, label=r"$\phi_{b\bar{b}}$ [GeV]")]),
#		"ratio" : HistConf( [Axis( coll="Higgs", field="ratio", bins=30, start=0, stop=2, label=r"$\frac{p_{T,last}}{p_{T,first}}$")])
		# "cos_theta" : HistConf( [Axis( coll="Spin", field="cos_theta", bins=20, start=-1, stop=1, label=r"$\cos \theta^*$")]),
		# "c_hel" : HistConf( [Axis( coll="Spin", field="c_hel", bins=20, start=-1, stop=1, label=r"$c_{hel}$")]),
		"c_han" : HistConf( [Axis( coll="Spin", field="c_han", bins=20, start=-1, stop=1, label=r"$c_{han}$")]),
		# "LLbar_DeltaPhi": HistConf([Axis(coll="Spin", field="LLbar_DeltaPhi", bins=20, start=0, stop=3.14, label=r"$|\phi_{\ell}-\phi_{\bar{\ell}}}|^{Lab}$")]),
		# "LLbar_DeltaEta": HistConf([Axis(coll="Spin", field="LLbar_DeltaEta", bins=20, start=0, stop=10, label=r"$|\eta_{\ell}-\eta_{\bar{\ell}}}|^{Lab}$")]),
		# "b1k_plus_b2k": HistConf([Axis(coll="Spin", field="b1k_plus_b2k", bins=20, start=-2, stop=2, label=r"$b_{1}^k + b_{2}^k$")]),
		# "b1r_plus_b2r": HistConf([Axis(coll="Spin", field="b1r_plus_b2r", bins=20, start=-2, stop=2, label=r"$b_{1}^r + b_{2}^r$")]),
		# "b1n_plus_b2n": HistConf([Axis(coll="Spin", field="b1n_plus_b2n", bins=20, start=-2, stop=2, label=r"$b_{1}^n + b_{2}^n$")]),
		# "b1k_minus_b2k": HistConf([Axis(coll="Spin", field="b1k_minus_b2k", bins=20, start=-2, stop=2, label=r"$b_{1}^k - b_{2}^k$")]),
		# "b1r_minus_b2r": HistConf([Axis(coll="Spin", field="b1r_minus_b2r", bins=20, start=-2, stop=2, label=r"$b_{1}^r - b_{2}^r$")]),
		# "b1n_minus_b2n": HistConf([Axis(coll="Spin", field="b1n_minus_b2n", bins=20, start=-2, stop=2, label=r"$b_{1}^n - b_{2}^n$")]),
		# "Ckk": HistConf([Axis(coll="Spin", field="Ckk", bins=20, start=-1, stop=1, label=r"$C_{kk}$")]),
		# "Crr": HistConf([Axis(coll="Spin", field="Crr", bins=20, start=-1, stop=1, label=r"$C_{rr}$")]),
		# "Cnn": HistConf([Axis(coll="Spin", field="Cnn", bins=20, start=-1, stop=1, label=r"$C_{nn}$")]),
		# "Crk_plus_Ckr": HistConf([Axis(coll="Spin", field="Crk_plus_Ckr", bins=20, start=-1, stop=1, label=r"$C_{rk} + C_{kr}$")]),
		# "Crn_plus_Cnr": HistConf([Axis(coll="Spin", field="Crn_plus_Cnr", bins=20, start=-1, stop=1, label=r"$C_{rn} + C_{nr}$")]),
		# "Ckn_plus_Cnk": HistConf([Axis(coll="Spin", field="Ckn_plus_Cnk", bins=20, start=-1, stop=1, label=r"$C_{kn} + C_{nk}$")]),
		# "Crk_minus_Ckr": HistConf([Axis(coll="Spin", field="Crk_minus_Ckr", bins=20, start=-1, stop=1, label=r"$C_{rk} - C_{kr}$")]),
		# "Crn_minus_Cnr": HistConf([Axis(coll="Spin", field="Crn_minus_Cnr", bins=20, start=-1, stop=1, label=r"$C_{rn} - C_{nr}$")]),
		# "Ckn_minus_Cnk": HistConf([Axis(coll="Spin", field="Ckn_minus_Cnk", bins=20, start=-1, stop=1, label=r"$C_{kn} - C_{nk}$")]),
		
		# "min_phi": HistConf([Axis(coll="Angular", field="min_phi", bins=20, start=0, stop=3.14, label=r"$min [ |\Delta\phi({t\bar{t}H})|]^{Lab}$")]),					
        # "min_phi_ttH_ZMF": HistConf([Axis(coll="Angular", field="min_phi_ttH_ZMF", bins=20, start=0, stop=3.14, label=r"$min [ |\Delta\phi({t\bar{t}H})|]^{t\bar{t}H}$")]),
        # "min_phi_tt_ZMF": HistConf([Axis(coll="Angular", field="min_phi_tt_ZMF", bins=20, start=0, stop=3.14, label=r"$min [ |\Delta\phi({t\bar{t}H})|]^{t\bar{t}}$")]),
        # "TTbar_phi": HistConf([Axis(coll="Angular", field="TTbar_phi", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({t\bar{t}})|^{Lab}$")]),					
        # "TTbar_phi_ttH_ZMF": HistConf([Axis(coll="Angular", field="TTbar_phi_ttH_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({t\bar{t}})|^{t\bar{t}H}$")]),
        # "TTbar_phi_tt_ZMF": HistConf([Axis(coll="Angular", field="TTbar_phi_tt_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({t\bar{t}})|^{t\bar{t}}$")]),
        # "THiggs_phi": HistConf([Axis(coll="Angular", field="THiggs_phi", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({tH})|^{Lab}$")]),					
        # "THiggs_phi_ttH_ZMF": HistConf([Axis(coll="Angular", field="THiggs_phi_ttH_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({tH})|^{t\bar{t}H}$")]),
        # "THiggs_phi_tt_ZMF": HistConf([Axis(coll="Angular", field="THiggs_phi_tt_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({tH})|^{t\bar{t}}$")]),
        # "TbarHiggs_phi": HistConf([Axis(coll="Angular", field="TbarHiggs_phi", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({\bar{t}H})|^{Lab}$")]),					
        # "TbarHiggs_phi_ttH_ZMF": HistConf([Axis(coll="Angular", field="TbarHiggs_phi_ttH_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({\bar{t}H})|^{t\bar{t}H}$")]),
        # "TbarHiggs_phi_tt_ZMF": HistConf([Axis(coll="Angular", field="TbarHiggs_phi_tt_ZMF", bins=20, start=0, stop=3.14, label=r"$|\Delta\phi({\bar{t}H})|^{t\bar{t}}$")]),
        
        #"min_eta": HistConf([Axis(coll="Angular", field="min_eta", bins=20, start=0, stop=5, label=r"$min [ |\Delta\eta({t\bar{t}H})|]^{Lab}$")]),					
		#"min_eta_ttH_ZMF": HistConf([Axis(coll="Angular", field="min_eta_ttH_ZMF", bins=20, start=0, stop=5, label=r"$min [ |\Delta\eta({t\bar{t}H})|]^{t\bar{t}H}$")]),
		#"min_eta_tt_ZMF": HistConf([Axis(coll="Angular", field="min_eta_tt_ZMF", bins=20, start=0, stop=5, label=r"$min [ |\Delta\eta({t\bar{t}H})|]^{t\bar{t}}$")]),
		#"TTbar_eta": HistConf([Axis(coll="Angular", field="TTbar_eta", bins=20, start=0, stop=5, label=r"$|\Delta\eta({t\bar{t}})|^{Lab}$")]),					
		#"TTbar_eta_ttH_ZMF": HistConf([Axis(coll="Angular", field="TTbar_eta_ttH_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({t\bar{t}})|^{t\bar{t}H}$")]),
		#"TTbar_eta_tt_ZMF": HistConf([Axis(coll="Angular", field="TTbar_eta_tt_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({t\bar{t}})|^{t\bar{t}}$")]),
		#"THiggs_eta": HistConf([Axis(coll="Angular", field="THiggs_eta", bins=20, start=0, stop=5, label=r"$|\Delta\eta({tH})|^{Lab}$")]),					
		#"THiggs_eta_ttH_ZMF": HistConf([Axis(coll="Angular", field="THiggs_eta_ttH_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({tH})|^{t\bar{t}H}$")]),
		#"THiggs_eta_tt_ZMF": HistConf([Axis(coll="Angular", field="THiggs_eta_tt_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({tH})|^{t\bar{t}}$")]),
		#"TbarHiggs_eta": HistConf([Axis(coll="Angular", field="TbarHiggs_eta", bins=20, start=0, stop=5, label=r"$|\Delta\eta({\bar{t}H})|^{Lab}$")]),					
		#"TbarHiggs_eta_ttH_ZMF": HistConf([Axis(coll="Angular", field="TbarHiggs_eta_ttH_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({\bar{t}H})|^{t\bar{t}H}$")]),
		#"TbarHiggs_eta_tt_ZMF": HistConf([Axis(coll="Angular", field="TbarHiggs_eta_tt_ZMF", bins=20, start=0, stop=5, label=r"$|\Delta\eta({\bar{t}H})|^{t\bar{t}}$")]),
		
		"min_dR": HistConf([Axis(coll="Angular", field="min_dR", bins=20, start=0, stop=5, label=r"$min [ \Delta R({t\bar{t}H})]^{Lab}$")]),					
		"min_dR_ttH_ZMF": HistConf([Axis(coll="Angular", field="min_dR_ttH_ZMF", bins=20, start=0, stop=5, label=r"$min [ \Delta R({t\bar{t}H})]^{t\bar{t}H}$")]),
		"TTbar_dR": HistConf([Axis(coll="Angular", field="TTbar_dR", bins=20, start=0, stop=6, label=r"$\Delta R({t\bar{t}})^{Lab}$")]),					
		"TTbar_dR_ttH_ZMF": HistConf([Axis(coll="Angular", field="TTbar_dR_ttH_ZMF", bins=20, start=0, stop=6, label=r"$\Delta R({t\bar{t}})^{t\bar{t}H}$")]),
		"THiggs_dR": HistConf([Axis(coll="Angular", field="THiggs_dR", bins=20, start=0, stop=5, label=r"$\Delta R({tH})^{Lab}$")]),					
		"THiggs_dR_ttH_ZMF": HistConf([Axis(coll="Angular", field="THiggs_dR_ttH_ZMF", bins=20, start=0, stop=6, label=r"$\Delta R({tH})^{t\bar{t}H}$")]),
		"TbarHiggs_dR": HistConf([Axis(coll="Angular", field="TbarHiggs_dR", bins=20, start=0, stop=6, label=r"$\Delta R({\bar{t}H})^{Lab}$")]),					
		"TbarHiggs_dR_ttH_ZMF": HistConf([Axis(coll="Angular", field="TbarHiggs_dR_ttH_ZMF", bins=20, start=0, stop=6, label=r"$\Delta R({\bar{t}H})^{t\bar{t}H}$")]),

        
        #**jet_hists(coll="JetGood", pos=2),
    }
)


run_options = {
        "executor"       : "parsl-condor",
        "env"            : "conda",
        "workers"        : 1,
        "scaleout"       : 10,
        "queue"          : "standard",
        "walltime"       : "00:30:00",
        "mem_per_worker" : "4GB", # GB
        "disk_per_worker" : "2GB", # GB
        "exclusive"      : False,
        "chunk"          : 200000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False,
        
    }
   



