import awkward as ak
from coffea.util import load
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import os
import argparse
from pocket_coffea.utils.plot_utils import Shape

Shapes = [

]


var_argNames = {
	"nJets" : "nJets",
	"deltaR" : "deltaRbb_min",
	"nBJets" : "nBJets",
	"nLeptons" : "nLeptons",
	"Lpt1" : "LeptonGood_pt_1",
	"Lpt2" : "LeptonGood_pt_2",
	"Leta1" : "LeptonGood_eta_1",
	"Leta2" : "LeptonGood_eta_2",
	"Lphi1" : "LeptonGood_phi_1",
	"Lphi2" : "LeptonGood_phi_2",
	"bjet_pt_1" : "bjet_pt_1",
	"bjet_pt_2" : "bjet_pt_2",
	"bjet_pt_3" : "bjet_pt_3",
	"bjet_pt_4" : "bjet_pt_4",
	"bjet_eta_1": "bjet_eta_1",
	"bjet_eta_2": "bjet_eta_2",
	"bjet_eta_3": "bjet_eta_3",
	"bjet_eta_4": "bjet_eta_4",
	"bjet_phi_1": 'bjet_phi_1',
	"bjet_phi_2": 'bjet_phi_2',
	"bjet_phi_3": 'bjet_phi_3',
	"bjet_phi_4": 'bjet_phi_4',
	"met_pt" : "MET_pt"
}

parser = argparse.ArgumentParser()
parser.add_argument("--var", "-v", default="all", help=f"variable to be plotted \n options: {var_argNames.keys()} \n default: \"all\"" )
parser.add_argument("--logscale", "-log", default=False, help="apply logscale, default: False" )
parser.add_argument("--Output", "-o", default=False, help="choose directory with output file" )
args = parser.parse_args()
var_args = [args.var]
log = args.logscale
Output = args.Output
print("creating plots for", var_args[0], "with logscale", str(log))


if os.path.isdir(Output):
	try:
	   os.makedirs(Output + "/my_plots")
	except FileExistsError:
	   pass

	try:
	   os.makedirs(Output + "/my_plots/3j3b")
	   os.makedirs(Output + "/my_plots/baseline")
	   
	except FileExistsError:
	   pass

def sum_over_hists(file):
	var = load(file)["variables"]
	outdict = {}
	for H in var.keys():
		concdict = {}
		samples = var[H]
		outdict.update({f"{str(H)}": 0})
		for s in samples.keys():
			datasets = samples[s]
			concdict.update({f"{str(s)}": 0})
			for d in datasets.keys():
				hist_obj = datasets[d]
				#print(hist_obj)
				if concdict[str(s)] == 0:
					concdict.update({f"{str(s)}": hist_obj})
				else:
					concdict[str(s)] += hist_obj
			outdict.update({f"{str(H)}": concdict})
	return outdict
	
filename = Output + "/output_all.coffea"
output = sum_over_hists(filename)
# ~ print(output["nJets"].keys())


def plotter(var,output,log):
	mcSamps = {
		"ttH": {"fancyName": r"$t\overline{t}H$",
				"color": "blue",
				"sampName": ["TTH_Hto2B"], 
				"hists": {},
				"var": {}},
		"VV":  {"fancyName": r"$VV$",
				"color": "darkgreen",
				"sampName": ["WWto2L2Nu","WZto3LNu"], 
				# "sampName": ["WWto2L2Nu","ZZto2L2Nu","WZto3LNu"], 
				"hists": {},
				"var": {}},
		"ttV": {"fancyName": r"$t\overline{t}V$",
				"color": "blue",
				"sampName": ["TTZ-ZtoQQ","TTNuNu","TTLNu","TTLL"], 
				"hists": {},
				"var": {}},
		"ttB": {"fancyName": r"$t\overline{t}B$",
				"color": "maroon",
				"sampName": ["TTTo2L2Nu__ttB","TTToLNu2Q__ttB"], 
				"hists": {},
				"var": {}},
		"ttC": {"fancyName": r"$t\overline{t}C$",
				"color": "firebrick",
				"sampName": ["TTTo2L2Nu__ttC","TTToLNu2Q__ttC"], 
				"hists": {},
				"var": {}},
		"ttLF": {"fancyName": r"$t\overline{t}LF$",
				"color": "salmon",
				"sampName": ["TTTo2L2Nu__ttLF","TTToLNu2Q__ttLF"], 
				"hists": {},
				"var": {}},
		"t":  {"fancyName": r"$Single Top$",
				"color": "blueviolet",
				"sampName": ["TWto2L2Nu"], 
				"hists": {},
				"var": {}},
		"DY": {"fancyName": r"DY",
				"color": "yellow",
				"sampName": ["DYJetsToLL"], 
				"hists": {},
				"var": {}},
	}
	
	
	
	betterNames = {
	"nJets": r"$N_{JetGood}$",
	"deltaRbb_min": r"min($\Delta R(bb)$)",
	"nBJets": r"$N_{BJetGood}$",
	"nLeptons" : r"$N_{LeptonsGood}$",
	'LeptonGood_eta_1' : r"$\eta^{l_L}$",
	'LeptonGood_pt_1' : r"$P_T^{l_L}~[GeV]$",
	'LeptonGood_phi_1': r"$\Phi^{l_L}$" ,
	'LeptonGood_pdgId_1' : 'LeptonGood_pdgId_1',
	'LeptonGood_eta_2' : r"$\eta^{l_T}$",
	'LeptonGood_pt_2'  : r"$P_T^{l_T}~[GeV]$",
	'LeptonGood_phi_2' : r"$\Phi^{l_T}$",
	'LeptonGood_pdgId_2' : 'LeptonGood_pdgId_2',
	"bjet_pt_1" : r"$p_T^{b_1}$",
	"bjet_pt_2" : r"$p_T^{b_2}$",
	"bjet_pt_3" : r"$p_T^{b_3}$",
	"bjet_pt_4" : r"$p_T^{b_4}$",
	"bjet_eta_1" : r"$\eta^{b_1}$",
	"bjet_eta_2" : r"$\eta^{b_2}$",
	"bjet_eta_3" : r"$\eta^{b_3}$",
	"bjet_eta_4" : r"$\eta^{b_4}$",
	"bjet_phi_1" : r"$\Phi^{b_1}$",
	"bjet_phi_2" : r"$\Phi^{b_2}$",
	"bjet_phi_3" : r"$\Phi^{b_3}$",
	"bjet_phi_4" : r"$\Phi^{b_4}$",
	"MET_pt" : r"$m_{ll}$"
	}
	
	
	data_keys = [
		'DATA_MuonEG__rmOverlap',
		'DATA_DoubleMuon__rmOverlap',
		'DATA_EGamma__rmOverlap'
	]
	
	data = 0
	dataHist = 0
	# ~ print(output.keys())
	data = output[var]
	for key in data_keys:
		dataHist += data[key]
	dataArr = dataHist.to_numpy(flow=True)
	for entry in dataArr[0]:
		entry[1]+=entry[0]
		entry[-2]+=entry[-1]
	dataCounts = {
		"3j3b" : dataArr[0][0][1:-1],
		"baseline" : dataArr[0][1][1:-1]
	}
	dataBins = dataArr[-1][1:-1]
	# ~ print(dataArr)
	#join MC Sets
	
	keyOrder = ["VV","ttV","t","DY","ttLF","ttC","ttB"]
	allSampsOrder = ["VV","ttV","t","DY","ttLF","ttC","ttB","ttH"]
	MClabels = [mcSamps[key]["fancyName"] for key in keyOrder]
	
	# ~ MC Hists to Numpy
	for key in keyOrder:
		MC_Hist = 0
		for samps in mcSamps[key]["sampName"]:
			MC_Hist += data[samps]
		MCArr = MC_Hist.to_numpy(flow=True)
		for entry in MCArr[0]:
			entry[0][1]+=entry[0][0]
			entry[0][-2]+=entry[0][-1]
		mcSamps[key].update({"hists" : {
			"3j3b" : MCArr[0][0][0][1:-1],
			"baseline" : MCArr[0][1][0][1:-1]}
			})
	
	mccolors = [mcSamps[key]["color"] for key in keyOrder]
	
	
	ttH_data = data["TTH_Hto2B"].to_numpy(flow=True)
	for entry in ttH_data[0]:
		entry[0][1]+=entry[0][0]
		entry[0][-2]+=entry[0][-1]
	mcSamps["ttH"].update({"hists" : {
			"3j3b" : ttH_data[0][0][0][1:-1],
			"baseline" : ttH_data[0][1][0][1:-1]}
			})
	
	shifts = ['pileupDown', 'pileupUp', 'sf_ele_idDown', 'sf_ele_idUp', 'sf_ele_recoDown', 'sf_ele_recoUp']
	cat = ["3j3b","baseline"]
	
	# ~ mcStat: overflow???
	for key in allSampsOrder:
		MC_Hist = 0
		for samps in mcSamps[key]["sampName"]:
			MC_Hist += data[samps]
		MCvar = MC_Hist.variances(flow=True)
		for i in range(2):
			MCvar[i][0][1] += MCvar[i][0][0]
			MCvar[i][0][-2] += MCvar[i][0][-1]
			# ~ print("sus", MCvar[1][0][1:-1])
		mcSamps[key]["var"].update({"mcstat" : {
			"3j3b" : {"up" : MCvar[0][0][1:-1], "down" : MCvar[0][0][1:-1]},
			"baseline" : {"up" : MCvar[1][0][1:-1], "down" : MCvar[1][0][1:-1]}}
			})


	hep.style.use("CMS")
	for key in dataCounts.keys():
		sys_var = []
		Up = []
		Down = []
		for samp in keyOrder:
			j=0
			if key == "baseline": j=1
			elif key == "3j3b" : j=0
			for samps in mcSamps[samp]["sampName"]:
				Hist= data[samps]
				nom = Hist.to_numpy(flow=True)[0][j][0]
				for i in range(1,len(shifts),2):
					Up.append(Hist.to_numpy(flow=True)[0][j][i+1]-nom)
					Down.append(Hist.to_numpy(flow=True)[0][j][i]-nom)
		all_up = []
		all_down = []
		for i in range(len(Up)):
			all_up.append(np.where(Up[i]>0,Up[i],0))
			all_up.append(np.where(Down[i]>0,Down[i],0))
			all_down.append(np.where(Up[i]<0,Up[i],0))
			all_down.append(np.where(Down[i]<0,Down[i],0))
		sys_var_up = np.zeros(len(dataBins)+1)
		sys_var_down = np.zeros(len(dataBins)+1)
		for i in all_up:
			sys_var_up += i**2
		for i in all_down:
			sys_var_down += i**2

		sys_var_up[-2] += sys_var_up[-1]
		sys_var_up[1] += sys_var_up[0]
		sys_var_down[-2] += sys_var_down[-1]
		sys_var_down[1] += sys_var_down[0]
		sys_var_up = sys_var_up[1:-1]
		sys_var_down = sys_var_down[1:-1]
		sys_unc_up = np.sqrt(sys_var_up)
		sys_unc_down = np.sqrt(sys_var_down)
		# ~ print("sysunc")
		# ~ print(sys_unc_up)
		# ~ print(sys_unc_down)
			
		
		var_up = np.zeros(len(dataBins)-1)
		var_down = np.zeros(len(dataBins)-1)
		for samp in mcSamps.keys():
			for unc in mcSamps[samp]["var"].keys():
				var_up += mcSamps[samp]["var"][unc][key]["up"]
				var_down += mcSamps[samp]["var"][unc][key]["down"]
		
		unc_up = np.sqrt(var_up+sys_var_up)
		unc_down = np.sqrt(var_down+sys_var_down)
		
		# ~ band = np.array([var_up,var_down])
		
		symmetric = []
		for i in range(len(var_down)):
			symmetric.append(np.sqrt(var_down[i]))
		symmetric = np.array(symmetric)
		# ~ print(band.shape)
		
		fig, ax = plt.subplots(2,1,height_ratios=[3,1], figsize=(10,10),sharex=True)
		hep.cms.label("Private Work", data=True, lumi=26.67,year=2022, ax=ax[0],com=13.6)
		nJet_sim_arr = []
		for samp in keyOrder:
			nJet_sim_arr.append((mcSamps[samp]["hists"][key], dataBins))
		hep.histplot((dataCounts[key],dataBins), label="Data", histtype="errorbar",color="black",ax=ax[0],yerr=np.sqrt(dataCounts[key]))
		hep.histplot(nJet_sim_arr, stack=True, histtype="fill", label=MClabels, color=mccolors, ax=ax[0])
		hep.histplot((mcSamps["ttH"]["hists"][key],dataBins), histtype="step", label=r"$t\overline{t}H$", color="red", ax=ax[0])

		
		MC_all = []
		for samps in keyOrder:
			MC_all.append(mcSamps[samps]["hists"][key])
		MC_summed = sum(MC_all)
		hep.histplot((dataCounts[key]/MC_summed, dataBins), yerr=np.sqrt(dataCounts[key])/np.abs(MC_summed), histtype="errorbar", color="black", ax=ax[1])
		# ~ hep.histplot((MC_summed,dataBins), histtype="band",color="b", label=r"uncertainties",yerr=symmetric, ax=ax[0])
		# ~ hep.histplot((np.ones(len(dataBins)-1),dataBins), histtype="band",color="b", label=r"MC_stat",yerr=dataCounts[key]/MC_summed**2*symmetric, ax=ax[1])
		binsize = dataBins[1]-dataBins[0]
		errorbins = dataBins[:-1]+binsize/2
		ax[0].fill_between(x=errorbins,step="mid",facecolor="none",y1=(MC_summed-unc_down),y2=(MC_summed+unc_up),edgecolor="darkblue",linewidth=0.0,hatch="///",label="full MCerr")
		ax[1].fill_between(x=errorbins,step="mid",facecolor="darkblue",y1=(MC_summed-unc_down)/MC_summed,y2=(MC_summed+unc_up)/MC_summed,edgecolor="gold",linewidth=0.0,label="full MCerr")
		ax[1].fill_between(x=errorbins,step="mid",facecolor="gold",y1=(MC_summed-sys_unc_down)/MC_summed,y2=(MC_summed+sys_unc_up)/MC_summed,edgecolor="yellow",linewidth=0.0,label="MCsys")
		# ~ ax[1].fill_between(x=errorbins,step="mid",facecolor="darkgreen",y1=(MC_summed-sys_unc_down)/MC_summed,y2=(MC_summed+sys_unc_up)/MC_summed,hatch="///",edgecolor="yellow",linewidth=0.0)
		# ~ print("MCsummed",MC_summed)
		
		ax[0].set_ylabel("Counts")
		ax[1].set_ylabel("Data/MC")
		ax[1].axhline(y=1, xmin=0, xmax=12, color="red",linestyle="--")
		ax[0].legend(prop={"size":20},ncol=3)
		if log:
			ax[0].set_yscale("log")
		ax[1].set_xlabel(f"{betterNames[var]}",size=20)
		ax[0].set_ylim(bottom = 1)
		ax[1].set_ylim(top = 1.5)
		ax[1].set_ylim(bottom = 0.5)
		plt.tight_layout(pad=0.1)
		hep.plot.yscale_legend(ax=ax[0])
		ax[1].legend(loc="best",prop={'size': 18})
		plt.savefig(f"{Output}/my_plots/{key}/{var}_{key}.pdf")
		plt.close("all")
		


if "all" in var_args:
	for i in var_argNames.keys():
		plotter(var_argNames[i],output,log)
else:
	for i in var_argNames.keys():
		if i in var_args:
			plotter(var_argNames[i],output,log)

	
