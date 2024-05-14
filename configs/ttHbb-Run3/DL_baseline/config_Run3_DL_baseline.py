from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_functions import get_nObj_eq, get_nObj_min, get_nObj_less, get_HLTsel, get_nBtagMin, get_nElectron, get_nMuon
from pocket_coffea.parameters.histograms import *
from pocket_coffea.parameters.cuts import passthrough

import workflow
from workflow import ttHbb_Run3Test

import custom_cut_functions
import custom_cuts
from custom_cut_functions import *
from custom_cuts import *
# ~ from params.binning import bins
# ~ from params.axis_settings import axis_settings
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")


lepton_categories_dict = {
    "ee" : [custom_get_dilepton(selec = "ee")],
    "em" : [custom_get_dilepton(selec = "em")],
    "mm" : [custom_get_dilepton(selec = "mm")],
}

year = "2022_postEE"
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection_DL_Run3.yaml",
                                                  f"{localdir}/params/triggers_DL.yaml",
                                                  update=True)
parameters["run_period"] = "Run3"

cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons":[
			f"{localdir}/datasets/Run3_MC_ttBkg.json",
			f"{localdir}/datasets/Run3_DATA_DoubleEle.json",
			f"{localdir}/datasets/Run3_DATA_DoubleMuon.json",
			f"{localdir}/datasets/Run3_DATA_MuonEG.json",
			f"{localdir}/datasets/Run3_MC_SIG.json",
			f"{localdir}/datasets/Run3_MC_otherBKG.json",
			
		],
        "filter" : {
            "samples": [
				"DATA_EGamma",
				"DATA_DoubleMuon",
				"DATA_MuonEG",
				"TTTo2L2Nu",
				"TTToLNu2Q",
				"DYJetsToLL",
				"WtoLNu+jets",
				"TWto2L2Nu",
				"WWto2L2Nu",
				"WZto3LNu",
				"ZZto2L2Nu",
				"TTLL",
				"TTNuNu",
				"TTLNu",
				"TTZ-ZtoQQ",
				"TTH_Hto2B",
            ],
            "samples_exclude" : [],
            "year": [year]
            # ~ "year": ['2022_preEE','2022_postEE','2023_preBPix','2023_postBPix']
        },
        "subsamples": {
            'DATA_MuonEG'  : {'rmOverlap' : [
                get_HLTsel(primaryDatasets=["MuonEG"])] #i.e. HLT_eleXX_muXX
            },
            'DATA_DoubleMuon'  : {'rmOverlap' : [
                get_HLTsel(primaryDatasets=["DoubleMuon", "SingleMuon"]), #i.e. HLT_muXX_muXX OR HLT_muXX
                get_HLTsel(primaryDatasets=["MuonEG"], invert=True)] #i.e. NOT selected by HLT_eleXX_muXX
            },
            'DATA_EGamma'  : {'rmOverlap' : [
                get_HLTsel(primaryDatasets=["DoubleEle", "SingleEle"]), #i.e. HLT_eleXX_eleXX OR HLT_eleXX
                get_HLTsel(primaryDatasets=["MuonEG", "DoubleMuon", "SingleMuon"], invert=True)] #i.e. NOT selected by one above
            },
            'TTTo2L2Nu'  : {
                'ttB' : [get_ttB_id(ttBid = "ttB")],
			    'ttC' : [get_ttB_id(ttBid = "ttC")],
                'ttLF' : [get_ttB_id(ttBid = "ttLF")]
            },
            'TTToLNu2Q'  : {
                'ttB' : [get_ttB_id(ttBid = "ttB")],
			    'ttC' : [get_ttB_id(ttBid = "ttC")],
                'ttLF' : [get_ttB_id(ttBid = "ttLF")]
            }
        }
    },

    workflow = ttHbb_Run3Test,
    # ~ workflow_options = {"parton_jet_min_dR": 0.3},
    
    skim = [
        get_nObj_min(2, 15., "Jet"),
        get_HLTsel(primaryDatasets=["SingleEle", "SingleMuon", "MuonEG", "DoubleEle", "DoubleMuon"])
        ],
    preselections = [dileptonic_presel],
    categories = {
        "baseline": [passthrough],
        **lepton_categories_dict,
    },

    weights = {
        "common": {
			"inclusive": [
                "genWeight",
                "lumi",
                "XS",
                "pileup",
                "sf_ele_reco",
                "sf_ele_id",
                # "sf_ele_trigger",
                "sf_mu_id",
                "sf_mu_iso", 
                # "sf_mu_trigger",
                # "sf_btag", 
                # "sf_btag_calib",
                # "sf_jet_puId"
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
                "inclusive": [
                    "pileup", 
                    "sf_ele_reco", 
                    "sf_ele_id", 
                    "sf_mu_id", 
                    "sf_mu_iso", 
                    # "sf_btag"
                    ],
                "bycategory" : {
                }
            },
            "bysample": {
            }
        },
        "shape": {
            "common":{
                "inclusive": []
            }
        }
    },
    
    variables = {
        "deltaRbb_min" : HistConf(
            [Axis(coll="events", field="deltaRbb_min", bins=50, start=0, stop=5,
                  label="$\Delta R_{bb}$", overflow=True, underflow=True)]
        ),
        **count_hist(name="nJets", coll="JetGood",bins=10, start=2, stop=12),
        **count_hist(name="nBJets", coll="BJetGood",bins=14, start=0, stop=14),
        **count_hist(name="nLeptons", coll="LeptonGood",bins=3, start=0, stop=3),
        **lepton_hists(coll="LeptonGood",pos=0),
        **lepton_hists(coll="LeptonGood",pos=1),
        **jet_hists(name="bjet",coll="BJetGood", pos=0),
        **jet_hists(name="bjet",coll="BJetGood", pos=1),
        **jet_hists(name="bjet",coll="BJetGood", pos=2),
        **jet_hists(name="bjet",coll="BJetGood", pos=3),
        **jet_hists(name="bjet",coll="BJetGood", pos=4),
        **met_hists(coll="MET"),
    },
)

run_options = {
        "executor"       : "dask/slurm",
        "env"            : "conda",
        "workers"        : 1,
        "scaleout"       : 200,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
        "queue"          : "standard",
        "walltime"       : "12:00:00",
        "mem_per_worker" : "6GB", # GB
        "disk_per_worker" : "1GB", # GB
        "exclusive"      : False,
        "chunk"          : 500000,
        "retries"        : 50,
        "treereduction"  : 5,
        "adapt"          : False,
    }


if "dask"  in run_options["executor"]:
    import cloudpickle
    cloudpickle.register_pickle_by_value(workflow)
    cloudpickle.register_pickle_by_value(custom_cut_functions)
    cloudpickle.register_pickle_by_value(custom_cuts)
