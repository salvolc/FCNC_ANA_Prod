import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os

untruth = False






if untruth:
	###################################################################################################
	########################################Photonen###################################################
	###################################################################################################
	if "plots" not in os.listdir():
		os.mkdir("plots")
	os.chdir("plots")
	if "top" not in os.listdir():
		os.mkdir("top")
	if "photon" not in os.listdir():
		os.mkdir("photon")
	if "WBoson" not in os.listdir():
		os.mkdir("WBoson")
	if "bJets" not in os.listdir():
		os.mkdir("bJets")
	os.chdir("../")

	lower_range = 0
	upper_range = 1000
	high_bin = 128
	low_bin = 64

	ev1 = np.genfromtxt("data/tua_LH_Photon_PT.txt")
	ev2 = np.genfromtxt("data/tua_RH_Photon_PT.txt")
	ev3 = np.genfromtxt("data/tca_LH_Photon_PT.txt")
	ev4 = np.genfromtxt("data/tca_RH_Photon_PT.txt")


	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev1[ev1!=0],label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.title("Photon PT")
	plt.legend(loc="best")
	plt.savefig("plots/photon/tua_LH_Photon_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev2[ev2!=0],label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.title("Photon PT")
	plt.legend(loc="best")
	plt.savefig("plots/photon/tua_RH_Photon_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev3[ev3!=0],label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γT}$")
	plt.ylabel("N")
	plt.title("Photon PT")
	plt.legend(loc="best")
	plt.savefig("plots/photon/tca_LH_Photon_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev4[ev4!=0],label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.title("Photon PT")
	plt.legend(loc="best")
	plt.savefig("plots/photon/tca_RH_Photon_PT.pdf",bbox_inches='tight')
	plt.close()



	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev1,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(ev2,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title("Photon PT Right vs. Left Handed")
	plt.savefig("plots/photon/tua_PPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev3,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(ev4,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title("Photon PT Right vs. Left Handed")
	plt.savefig("plots/photon/tca_PPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()



	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev1,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(ev3,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"Photon PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/photon/LH_PPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(ev2,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(ev4,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{γ_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"Photon PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/photon/RH_PPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()


	###################################################################################################
	########################################BJets######################################################
	###################################################################################################

	bJ_tua_LH = np.genfromtxt("data/tua_LH_bJet_PT.txt")
	bJ_tua_RH = np.genfromtxt("data/tua_RH_bJet_PT.txt")
	bJ_tca_LH = np.genfromtxt("data/tca_LH_bJet_PT.txt")
	bJ_tca_RH = np.genfromtxt("data/tca_RH_bJet_PT.txt")

	bJ_tua_LH = bJ_tua_LH[bJ_tua_LH!=0]
	bJ_tua_RH = bJ_tua_RH[bJ_tua_RH!=0]
	bJ_tca_LH = bJ_tca_LH[bJ_tca_LH!=0]
	bJ_tca_RH = bJ_tca_RH[bJ_tca_RH!=0]


	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(bJ_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(bJ_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"bJet PT Right vs. Left Handed")
	plt.savefig("plots/bJets/tua_bJPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(bJ_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(bJ_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"bJet PT Right vs. Left Handed")
	plt.savefig("plots/bJets/tca_bJPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(bJ_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(bJ_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"bJet PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/bJets/LH_bJPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(bJ_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(bJ_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"bJet PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/bJets/RH_bJPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()


	###################################################################################################
	########################################W Boson Elektron###########################################
	###################################################################################################

	WEPT_tua_LH = np.genfromtxt("data/tua_LH_WE_PT.txt")
	WEPT_tua_RH = np.genfromtxt("data/tua_RH_WE_PT.txt")
	WEPT_tca_LH = np.genfromtxt("data/tca_LH_WE_PT.txt")
	WEPT_tca_RH = np.genfromtxt("data/tca_RH_WE_PT.txt")

	WEPT_tua_LH = WEPT_tua_LH[WEPT_tua_LH!=0]
	WEPT_tua_RH = WEPT_tua_RH[WEPT_tua_RH!=0]
	WEPT_tca_LH = WEPT_tca_LH[WEPT_tca_LH!=0]
	WEPT_tca_RH = WEPT_tca_RH[WEPT_tca_RH!=0]

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.legend(loc="best")
	plt.title(r"$W_e$-PT")
	plt.savefig("plots/WBoson/tua_LH_WE_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.legend(loc="best")
	plt.title(r"$W_e$-PT")
	plt.savefig("plots/WBoson/tua_RH_WE_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.legend(loc="best")
	plt.title(r"$W_e$-PT")
	plt.savefig("plots/WBoson/tca_LH_WE_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.title(r"$W_e$-PT")
	plt.legend(loc="best")
	plt.savefig("plots/WBoson/tca_RH_WE_PT.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tua_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WEPT_tua_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$W_e$ PT Right vs. Left Handed")
	plt.savefig("plots/WBoson/tua_WE_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WEPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{We_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$W_e$ PT Right vs. Left Handed")
	plt.savefig("plots/WBoson/tca_WE_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WEPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$W_e$ PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/LH_WEPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WEPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WEPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{b_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$W_e$ PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/RH_WEPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	###################################################################################################
	########################################W Boson PT#################################################
	###################################################################################################

	WPT_tua_LH = np.genfromtxt("data/tua_LH_W_PT.txt")
	WPT_tua_RH = np.genfromtxt("data/tua_RH_W_PT.txt")
	WPT_tca_LH = np.genfromtxt("data/tca_LH_W_PT.txt")
	WPT_tca_RH = np.genfromtxt("data/tca_RH_W_PT.txt")
	WPT_tua_LH = WPT_tua_LH[WPT_tua_LH!=0]
	WPT_tua_RH = WPT_tua_RH[WPT_tua_RH!=0]
	WPT_tca_LH = WPT_tca_LH[WPT_tca_LH!=0]
	WPT_tca_RH = WPT_tca_RH[WPT_tca_RH!=0]

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tuγ$ LH")
	plt.savefig("plots/WBoson/WPT_tua_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tuγ$ RH")
	plt.savefig("plots/WBoson/WPT_tua_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tcγ$ LH")
	plt.savefig("plots/WBoson/WPT_tca_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tcγ$ RH")
	plt.savefig("plots/WBoson/WPT_tca_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT Right vs. Left Handed")
	plt.savefig("plots/WBoson/tua_WPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT Right vs. Left Handed")
	plt.savefig("plots/WBoson/tca_WPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/LH_WPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W PT $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/RH_WPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	###################################################################################################
	########################################W Boson Mass###############################################
	###################################################################################################

	WM_tua_LH = np.genfromtxt("data/tua_LH_W_M.txt")
	WM_tua_RH = np.genfromtxt("data/tua_RH_W_M.txt")
	WM_tca_LH = np.genfromtxt("data/tca_LH_W_M.txt")
	WM_tca_RH = np.genfromtxt("data/tca_RH_W_M.txt")
	WM_tua_LH = WM_tua_LH[WM_tua_LH!=0]
	WM_tua_RH = WM_tua_RH[WM_tua_RH!=0]
	WM_tca_LH = WM_tca_LH[WM_tca_LH!=0]
	WM_tca_RH = WM_tca_RH[WM_tca_RH!=0]

	lower_range = 0
	upper_range = 200
	low_bin = 40
	mW = 80.387

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mW,color="r",alpha=0.6,label=r"$m_W$")
	plt.hist(WM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tuγ$ LH")
	plt.savefig("plots/WBoson/WM_tua_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mW,color="r",alpha=0.6,label=r"$m_W$")
	plt.hist(WM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tuγ$ RH")
	plt.savefig("plots/WBoson/WM_tua_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mW,color="r",alpha=0.6,label=r"$m_W$")
	plt.hist(WM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tcγ$ LH")
	plt.savefig("plots/WBoson/WM_tca_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mW,color="r",alpha=0.6,label=r"$m_W$")
	plt.hist(WM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tcγ$ RH")
	plt.savefig("plots/WBoson/WM_tca_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M Right vs. Left Handed")
	plt.savefig("plots/WBoson/tua_WM_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M Right vs. Left Handed")
	plt.savefig("plots/WBoson/tca_WM_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/LH_WM_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(WM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(WM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_{W_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"W M $tuγ$ vs. $tcγ$")
	plt.savefig("plots/WBoson/RH_WM_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	###################################################################################################
	########################################Top Quark Mass#############################################
	###################################################################################################

	tM_tua_LH = np.genfromtxt("data/tua_LH_t_M.txt")
	tM_tua_RH = np.genfromtxt("data/tua_RH_t_M.txt")
	tM_tca_LH = np.genfromtxt("data/tca_LH_t_M.txt")
	tM_tca_RH = np.genfromtxt("data/tca_RH_t_M.txt")
	tM_tua_LH = tM_tua_LH[tM_tua_LH!=0]
	tM_tua_RH = tM_tua_RH[tM_tua_RH!=0]
	tM_tca_LH = tM_tca_LH[tM_tca_LH!=0]
	tM_tca_RH = tM_tca_RH[tM_tca_RH!=0]


	lower_range = 50
	upper_range = 500
	low_bin = 64
	mt = 173.5

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mt,color="r",alpha=0.6,label=r"$m_t$")
	plt.hist(tM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tuγ$ LH")
	plt.savefig("plots/top/tM_tua_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mt,color="r",alpha=0.6,label=r"$m_t$")
	plt.hist(tM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tuγ$ RH")
	plt.savefig("plots/top/tM_tua_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mt,color="r",alpha=0.6,label=r"$m_t$")
	plt.hist(tM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tcγ$ LH")
	plt.savefig("plots/top/tM_tca_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.axvline(x=mt,color="r",alpha=0.6,label=r"$m_t$")
	plt.hist(tM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tcγ$ RH")
	plt.savefig("plots/top/tM_tca_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ Right vs. Left Handed")
	plt.savefig("plots/top/tua_tM_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ Right vs. Left Handed")
	plt.savefig("plots/top/tca_tM_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tM_tua_LH,label=r"$tuγ$-LH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tM_tca_LH,label=r"$tcγ$-LH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tuγ$ vs. $tcγ$")
	plt.savefig("plots/top/LH_tM_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tM_tua_RH,label=r"$tuγ$-RH",bins=low_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tM_tca_RH,label=r"$tcγ$-RH",bins=low_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$m_t$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$m_t$ $tuγ$ vs. $tcγ$")
	plt.savefig("plots/top/RH_tM_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()


	###################################################################################################
	########################################top Quark PT###############################################
	###################################################################################################

	tPT_tua_LH = np.genfromtxt("data/tua_LH_t_PT.txt")
	tPT_tua_RH = np.genfromtxt("data/tua_RH_t_PT.txt")
	tPT_tca_LH = np.genfromtxt("data/tca_LH_t_PT.txt")
	tPT_tca_RH = np.genfromtxt("data/tca_RH_t_PT.txt")
	tPT_tua_LH = tPT_tua_LH[tPT_tua_LH!=0]
	tPT_tua_RH = tPT_tua_RH[tPT_tua_RH!=0]
	tPT_tca_LH = tPT_tca_LH[tPT_tca_LH!=0]
	tPT_tca_RH = tPT_tca_RH[tPT_tca_RH!=0]

	lower_range = 0
	upper_range = 1000
	mt = 173.5

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tuγ$ LH")
	plt.savefig("plots/top/tPT_tua_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tuγ$ RH")
	plt.savefig("plots/top/tPT_tua_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tcγ$ LH")
	plt.savefig("plots/top/tPT_tca_LH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tcγ$ RH")
	plt.savefig("plots/top/tPT_tca_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ Right vs. Left Handed")
	plt.savefig("plots/top/tua_tPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tca_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ Right vs. Left Handed")
	plt.savefig("plots/top/tca_tPT_LH_vs_RH.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tua_LH,label=r"$tuγ$-LH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tPT_tca_LH,label=r"$tcγ$-LH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tuγ$ vs. $tcγ$")
	plt.savefig("plots/top/LH_tPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()

	plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
	plt.hist(tPT_tua_RH,label=r"$tuγ$-RH",bins=high_bin,lw=0.5,color="blue",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.hist(tPT_tca_RH,label=r"$tcγ$-RH",bins=high_bin,lw=0.5,color="red",fill=False,histtype='step',normed=1,range=(lower_range,upper_range))
	plt.xlabel(r"$P_{t_T}$")
	plt.ylabel("N")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	plt.title(r"$p_{T_t}$ $tuγ$ vs. $tcγ$")
	plt.savefig("plots/top/RH_tPT_tua_vs_tca.pdf",bbox_inches='tight')
	plt.close()


###################################################################################################
########################################TRUTH PLOTS################################################
###################################################################################################

samples = ["tua_LH","tua_RH","tca_LH","tca_RH"]
vari = ["PT","Eta","Phi","M"]
part = ["WBoson","TopQuark","Photon","BQuark","UQuark"]
vgle = ["LH vs RH" , "$tuγ$ vs. $tcγ$"]
utyp = ["tua","tca"]
hand = ["LH","RH"]

PREFIX = "/home/salv/Dokumente/Masterarbeit/MadGraph/FCNC_NLO_Analysis/"

datas = os.listdir(PREFIX+"data")

ntruth = 0
for i in range(len(datas)):
	if("truth" in datas[i]):
		ntruth += 1


for par in part:
	if(par+"_truth" not in os.listdir(PREFIX+"plots/")):
		os.chdir(PREFIX+"plots/")
		os.mkdir(par+"_truth/")
		os.chdir(PREFIX)

topmassup=180
topmassdown=165
bmassup=5.2
bmassdown=3.2
umassup=1
umassdown=0
wmassup=95
wmassdown=65
default_range_up=1000
default_range_down=0

default_nbins = 128
mass_nbins = 64

ratiox=10
ratioy=4

for sam in samples:
	for par in part:
		for var in vari:
			nbin = default_nbins
			lower_range=default_range_down
			upper_range=default_range_up

			if(par == "Photon" and var=="M"):
				continue
			if(par == "Photon" and var=="PT"):
				upper_range=default_range_up
				lower_range=default_range_down
				nbin=default_nbins
			if (var == "Eta" or var == "Phi"):
				lower_range=-5
				upper_range=5
				nbin = 32
			if(var=="M"):
				nbin=mass_nbins
				if("W" in par):
					lower_range=wmassdown
					upper_range=wmassup
				if("Top" in par):
					lower_range=topmassdown
					upper_range=topmassup
				if("BQuark" in par):
					lower_range=bmassdown
					upper_range=bmassup
				if("UQuark" in par):
					lower_range=umassdown
					upper_range=umassup

			ev = np.genfromtxt("data/"+sam+"_"+par+"_"+var+"_truth.txt")
			ev = ev[ev!=0]
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			plt.hist(ev,label=sam[0:2]+r"$γ$"+sam[3:],bins=nbin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
			plt.xlabel(r"$"+var+"("+par+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.legend(loc="best")
			plt.title(var+"("+par+") "+sam[0:2]+r"$γ$"+sam[3:])
			plt.savefig("plots/"+par+"_truth/"+sam+"_"+par+"_"+var+".pdf",bbox_inches='tight')
			plt.close()
		

# tcatua Plots
for ha in hand:
	for par in part:
		for var in vari:
			nbin = default_nbins
			lower_range=default_range_down
			upper_range=default_range_up

			if(par == "Photon" and var=="M"):
				continue
			if(par == "Photon" and var=="PT"):
				upper_range=default_range_up
				lower_range=default_range_down
				nbin=default_nbins
			if (var == "Eta" or var == "Phi"):
				lower_range=-5
				upper_range=5
				nbin = 32
			if(var=="M"):
				nbin=mass_nbins
				if("W" in par):
					lower_range=wmassdown
					upper_range=wmassup
				if("Top" in par):
					lower_range=topmassdown
					upper_range=topmassup
				if("BQuark" in par):
					lower_range=bmassdown
					upper_range=bmassup
				if("UQuark" in par):
					lower_range=umassdown
					upper_range=umassup
			ev1 = np.genfromtxt("data/"+"tca_"+ha+"_"+par+"_"+var+"_truth.txt")
			ev1 = ev1[ev1!=0]
			ev2 = np.genfromtxt("data/"+"tua_"+ha+"_"+par+"_"+var+"_truth.txt")
			ev2 = ev2[ev2!=0]
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			plt.hist(ev1,label=r"$tcγ$-"+ha,bins=nbin,lw=0.5,color="red",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
			plt.hist(ev2,label=r"$tuγ$-"+ha,bins=nbin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))

			plt.xlabel(r"$"+var+"("+par+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.legend(loc="best")
			plt.title(var+"("+par+") "+ha+r" $tcγ$ vs. $tuγ$")
			plt.savefig("plots/"+par+"_truth/"+ha+"_"+par+"_"+var+"_tua_vs_tca.pdf",bbox_inches='tight')
			plt.close()

# LHRH Plots
for ut in utyp:
	for par in part:
		for var in vari:
			nbin = default_nbins
			lower_range=default_range_down
			upper_range=default_range_up

			if(par == "Photon" and var=="M"):
				continue
			if(par == "Photon" and var=="PT"):
				upper_range=default_range_up
				lower_range=default_range_down
				nbin=default_nbins
			if (var == "Eta" or var == "Phi"):
				lower_range=-5
				upper_range=5
				nbin = 32
			if(var=="M"):
				nbin=mass_nbins
				if("W" in par):
					lower_range=wmassdown
					upper_range=wmassup
				if("Top" in par):
					lower_range=topmassdown
					upper_range=topmassup
				if("BQuark" in par):
					lower_range=bmassdown
					upper_range=bmassup
				if("UQuark" in par):
					lower_range=umassdown
					upper_range=umassup

			ev1 = np.genfromtxt("data/"+ut+"_LH_"+par+"_"+var+"_truth.txt")
			ev1 = ev1[ev1!=0]
			ev2 = np.genfromtxt("data/"+ut+"_RH_"+par+"_"+var+"_truth.txt")
			ev2 = ev2[ev2!=0]
			plt.figure(num=None, figsize=(ratiox,ratioy), dpi=80, facecolor='w', edgecolor='k')
			plt.hist(ev1,label=ut[0:2]+r"$γ$"+ut[3:]+"_LH",bins=nbin,lw=0.5,color="blue",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))
			plt.hist(ev2,label=ut[0:2]+r"$γ$"+ut[3:]+"_RH",bins=nbin,lw=0.5,color="red",fill=False,histtype='step',normed=False,range=(lower_range,upper_range))

			plt.xlabel(r"$"+var+"("+par+")"+r"$")
			plt.ylabel("N")
			plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			plt.legend(loc="best")
			plt.title(var+"("+par+") "+ut[0:2]+r"$γ$"+ut[3:]+r" LH vs. RH")
			plt.savefig("plots/"+par+"_truth/"+ut+"_"+par+"_"+var+"_LH_vs_RH.pdf",bbox_inches='tight')
			plt.close()
		



