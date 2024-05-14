from collections.abc import Iterable
import awkward as ak
import numpy as np
import correctionlib
from pocket_coffea.lib.triggers import get_trigger_mask
from pocket_coffea.lib.cut_functions import get_JetVetoMap_Mask
from pocket_coffea.lib.cut_definition import Cut


def dileptonic(events, params, year, processor_params, sample, isMC, **kwargs):
    # testing(processor_params, events, year)
    lepIds = abs(ak.mask(events.LeptonGood.pdgId, events.nLeptonGood == 2))

    is_em = ak.where(lepIds[:,0] + lepIds[:,1] == 24, True, False) #11 + 13 = 24
    is_em = is_em & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["MuonEG", "SingleMuon", "SingleEle"])
    
    is_mm = ak.where(lepIds[:,0] + lepIds[:,1] == 26, True, False) #13 + 13 = 26
    is_mm = is_mm & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["DoubleMuon", "SingleMuon"])
    
    is_ee = ak.where(lepIds[:,0] + lepIds[:,1] == 22, True, False) #11 + 11 = 22
    is_ee = is_ee & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["DoubleEle", "SingleEle"])
    
    if year in ["2022_preEE", "2022_postEE", "2023_preBPIX", "2023_postBPIX"]:
        mask_jetVetoMap = get_JetVetoMap_Mask(events, params, year, processor_params, sample, isMC, **kwargs)
    else:
        mask_jetVetoMap = True
    mask = (
        (events.nLeptonGood == 2)
        & (ak.firsts(events.LeptonGood.pt) >= params["pt_leading_lepton"])
        & (ak.mask(events.LeptonGood.pt>=params["pt_subleading_lepton"], ak.num(events.LeptonGood.pt)>=2)[:,1])
        & (ak.sum(events.LeptonGood.charge, axis=1) == 0)
        & (events.nJetGood >= params["njet"])
        & (events.nBJetGood >= params["nbjet"])
        & (events.PuppiMET.pt > params["met"])
        & mask_jetVetoMap
        & (is_em  
            | (is_mm & (
                (events.ll.mass > params["dy_window"]["stop"]) 
                | ((events.ll.mass > params["m_ee_mumu_min"])
                    & (events.ll.mass < params["dy_window"]["start"]))
            ))
            | (is_ee & (
                (events.ll.mass > params["dy_window"]["stop"]) 
                | ((events.ll.mass > params["m_ee_mumu_min"])
                    & (events.ll.mass < params["dy_window"]["start"]))
            ))
        )
    )
    return ak.where(ak.is_none(mask), False, mask)


    
def custom_dilepton(events, params, year, processor_params, **kwargs):
	if params["selec"] == "ee":
		SF = ((events.nMuonGood == 0) & (events.nElectronGood == 2))
		OS = events.ll.charge == 0

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF
		)
		return ak.where(ak.is_none(mask), False, mask)
	elif params["selec"] == "em":
		SF = ((events.nMuonGood == 1) & (events.nElectronGood == 1))
		OS = events.ll.charge == 0

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF
		)
		return ak.where(ak.is_none(mask), False, mask)
	elif params["selec"] == "mm":
		SF = ((events.nMuonGood == 2) & (events.nElectronGood == 0))
		OS = events.ll.charge == 0

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF
		)
		return ak.where(ak.is_none(mask), False, mask)
	else:
		raise Exception("selection name not valid") 



def ttB_masks(events, params, processor_params, year, isMC, **kwargs):
	""" Categorization of ttbar events according to genTtbar ID, see reference here:
	 https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher#Event_categorization_example_1 """ 
	
	genTtbarId = events["genTtbarId"]
	if params["ttBId"] == "ttB":
		return (((abs(genTtbarId) % 100) == 51)
			| ((abs(genTtbarId) % 100) == 52)
			| ((abs(genTtbarId) % 100) == 53)
			| ((abs(genTtbarId) % 100) == 54)
			| ((abs(genTtbarId) % 100) == 55))
	elif params["ttBId"] == "ttC":
		return (((abs(genTtbarId) % 100) == 41)
		    | ((abs(genTtbarId) % 100) == 42)
			| ((abs(genTtbarId) % 100) == 43)
			| ((abs(genTtbarId) % 100) == 44)
			| ((abs(genTtbarId) % 100) == 45))
	elif params["ttBId"] == "ttLF":
		return ~(((abs(genTtbarId) % 100) == 41)
			| ((abs(genTtbarId) % 100) == 42)
			| ((abs(genTtbarId) % 100) == 43)
			| ((abs(genTtbarId) % 100) == 44)
			| ((abs(genTtbarId) % 100) == 45)
			| ((abs(genTtbarId) % 100) == 51)
			| ((abs(genTtbarId) % 100) == 52)
			| ((abs(genTtbarId) % 100) == 53)
			| ((abs(genTtbarId) % 100) == 54)
			| ((abs(genTtbarId) % 100) == 55))
