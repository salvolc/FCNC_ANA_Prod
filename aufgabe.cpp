#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "RConfig.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h" 
#include "TMultiGraph.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TSystem.h"
#include <iostream>
#include "TLorentzVector.h"
#include "math.h"
#include "ExRootClasses.h"
#include "ExRootTreeReader.h"
#include "TLeaf.h"
#include "classes/DelphesClasses.h"


#include "/home/salv/eigen/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;

int get_nevents(string fileName);
MatrixXd get_eventdisplay(string fileName, int event);
void ladebalken(int i, int max);
void speichere(std::string name, MatrixXd data);
void speichere(std::string name, VectorXd data);
VectorXd get_event_MET(string fileName, int event);
TLorentzVector get_event_MET_neutrino(string fileName, int event);
MatrixXd get_eventdisplay_particle(string fileName, int event, int PID);
int get_numberOfPhotons(string fileName);
int get_numberOfJets(string fileName);
int get_numberOfbJets(string fileName);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[]);
int get_nTruth(string fileName, int PID);
VectorXd vertex_match(TClonesArray* TCP);


int error=0;

int main(int argc, char const *argv[])
{

	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/MG5_aMC_v2_6_0/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/Delphes-3.4.1/libDelphes.so");

	
	string fileName = "samples/tua_LH_PY8_ATLASDELPHES_50000.root";

	int nEvents = get_nevents(fileName);

	MatrixXd ev = get_eventdisplay(fileName,1);
	speichere("einevent",ev);

	string fileNames[4];
	fileNames[0] = "samples/tua_LH_PY8_ATLASDELPHES_50000.root";
	fileNames[1] = "samples/tua_RH_PY8_ATLASDELPHES_50000.root";
	fileNames[2] = "samples/tca_LH_PY8_ATLASDELPHES_50000.root";
	fileNames[3] = "samples/tca_RH_PY8_ATLASDELPHES_50000.root";

	string filePPTNames[4];filePPTNames[0] = "data/tua_LH_Photon_PT";filePPTNames[1] = "data/tua_RH_Photon_PT";filePPTNames[2] = "data/tca_LH_Photon_PT";filePPTNames[3] = "data/tca_RH_Photon_PT";
	string filebJetPTNames[4];	filebJetPTNames[0] = "data/tua_LH_bJet_PT";filebJetPTNames[1] = "data/tua_RH_bJet_PT";filebJetPTNames[2] = "data/tca_LH_bJet_PT";filebJetPTNames[3] = "data/tca_RH_bJet_PT";
	string fileWPTENames[4];	fileWPTENames[0] = "data/tua_LH_WE_PT";fileWPTENames[1] = "data/tua_RH_WE_PT";fileWPTENames[2] = "data/tca_LH_WE_PT";fileWPTENames[3] = "data/tca_RH_WE_PT";
	string fileWPTNames[4];		fileWPTNames[0] = "data/tua_LH_W_PT";fileWPTNames[1] = "data/tua_RH_W_PT";fileWPTNames[2] = "data/tca_LH_W_PT";fileWPTNames[3] = "data/tca_RH_W_PT";
	string fileWMNames[4];		fileWMNames[0] = "data/tua_LH_W_M";fileWMNames[1] = "data/tua_RH_W_M";fileWMNames[2] = "data/tca_LH_W_M";fileWMNames[3] = "data/tca_RH_W_M";
	string filetPTNames[4];		filetPTNames[0] = "data/tua_LH_t_PT";filetPTNames[1] = "data/tua_RH_t_PT";filetPTNames[2] = "data/tca_LH_t_PT";filetPTNames[3] = "data/tca_RH_t_PT";
	string filetMNames[4];		filetMNames[0] = "data/tua_LH_t_M";filetMNames[1] = "data/tua_RH_t_M";filetMNames[2] = "data/tca_LH_t_M";filetMNames[3] = "data/tca_RH_t_M";
	

	string fileTPPTNames[4];fileTPPTNames[0] = "data/tua_LH_Photon_PT_truth";fileTPPTNames[1] = "data/tua_RH_Photon_PT_truth";fileTPPTNames[2] = "data/tca_LH_Photon_PT_truth";fileTPPTNames[3] = "data/tca_RH_Photon_PT_truth";
	string fileTPEtaNames[4];fileTPEtaNames[0] = "data/tua_LH_Photon_Eta_truth";fileTPEtaNames[1] = "data/tua_RH_Photon_Eta_truth";fileTPEtaNames[2] = "data/tca_LH_Photon_Eta_truth";fileTPEtaNames[3] = "data/tca_RH_Photon_Eta_truth";
	string fileTPPhiNames[4];fileTPPhiNames[0] = "data/tua_LH_Photon_Phi_truth";fileTPPhiNames[1] = "data/tua_RH_Photon_Phi_truth";fileTPPhiNames[2] = "data/tca_LH_Photon_Phi_truth";fileTPPhiNames[3] = "data/tca_RH_Photon_Phi_truth";
	

	string fileTWPTNames[4];fileTWPTNames[0] = "data/tua_LH_WBoson_PT_truth";fileTWPTNames[1] = "data/tua_RH_WBoson_PT_truth";fileTWPTNames[2] = "data/tca_LH_WBoson_PT_truth";fileTWPTNames[3] = "data/tca_RH_WBoson_PT_truth";
	string fileTWEtaNames[4];fileTWEtaNames[0] = "data/tua_LH_WBoson_Eta_truth";fileTWEtaNames[1] = "data/tua_RH_WBoson_Eta_truth";fileTWEtaNames[2] = "data/tca_LH_WBoson_Eta_truth";fileTWEtaNames[3] = "data/tca_RH_WBoson_Eta_truth";
	string fileTWPhiNames[4];fileTWPhiNames[0] = "data/tua_LH_WBoson_Phi_truth";fileTWPhiNames[1] = "data/tua_RH_WBoson_Phi_truth";fileTWPhiNames[2] = "data/tca_LH_WBoson_Phi_truth";fileTWPhiNames[3] = "data/tca_RH_WBoson_Phi_truth";
	string fileTWMNames[4];fileTWMNames[0] = "data/tua_LH_WBoson_M_truth";fileTWMNames[1] = "data/tua_RH_WBoson_M_truth";fileTWMNames[2] = "data/tca_LH_WBoson_M_truth";fileTWMNames[3] = "data/tca_RH_WBoson_M_truth";
	

	string fileTtPTNames[4];fileTtPTNames[0] = "data/tua_LH_TopQuark_PT_truth";fileTtPTNames[1] = "data/tua_RH_TopQuark_PT_truth";fileTtPTNames[2] = "data/tca_LH_TopQuark_PT_truth";fileTtPTNames[3] = "data/tca_RH_TopQuark_PT_truth";
	string fileTtEtaNames[4];fileTtEtaNames[0] = "data/tua_LH_TopQuark_Eta_truth";fileTtEtaNames[1] = "data/tua_RH_TopQuark_Eta_truth";fileTtEtaNames[2] = "data/tca_LH_TopQuark_Eta_truth";fileTtEtaNames[3] = "data/tca_RH_TopQuark_Eta_truth";
	string fileTtPhiNames[4];fileTtPhiNames[0] = "data/tua_LH_TopQuark_Phi_truth";fileTtPhiNames[1] = "data/tua_RH_TopQuark_Phi_truth";fileTtPhiNames[2] = "data/tca_LH_TopQuark_Phi_truth";fileTtPhiNames[3] = "data/tca_RH_TopQuark_Phi_truth";
	string fileTtMNames[4];fileTtMNames[0] = "data/tua_LH_TopQuark_M_truth";fileTtMNames[1] = "data/tua_RH_TopQuark_M_truth";fileTtMNames[2] = "data/tca_LH_TopQuark_M_truth";fileTtMNames[3] = "data/tca_RH_TopQuark_M_truth";


	string fileTbPTNames[4];fileTbPTNames[0] = "data/tua_LH_BQuark_PT_truth";fileTbPTNames[1] = "data/tua_RH_BQuark_PT_truth";fileTbPTNames[2] = "data/tca_LH_BQuark_PT_truth";fileTbPTNames[3] = "data/tca_RH_BQuark_PT_truth";
	string fileTbEtaNames[4];fileTbEtaNames[0] = "data/tua_LH_BQuark_Eta_truth";fileTbEtaNames[1] = "data/tua_RH_BQuark_Eta_truth";fileTbEtaNames[2] = "data/tca_LH_BQuark_Eta_truth";fileTbEtaNames[3] = "data/tca_RH_BQuark_Eta_truth";
	string fileTbPhiNames[4];fileTbPhiNames[0] = "data/tua_LH_BQuark_Phi_truth";fileTbPhiNames[1] = "data/tua_RH_BQuark_Phi_truth";fileTbPhiNames[2] = "data/tca_LH_BQuark_Phi_truth";fileTbPhiNames[3] = "data/tca_RH_BQuark_Phi_truth";
	string fileTbMNames[4];fileTbMNames[0] = "data/tua_LH_BQuark_M_truth";fileTbMNames[1] = "data/tua_RH_BQuark_M_truth";fileTbMNames[2] = "data/tca_LH_BQuark_M_truth";fileTbMNames[3] = "data/tca_RH_BQuark_M_truth";


	string fileTuPTNames[4];fileTuPTNames[0] = "data/tua_LH_UQuark_PT_truth";fileTuPTNames[1] = "data/tua_RH_UQuark_PT_truth";fileTuPTNames[2] = "data/tca_LH_UQuark_PT_truth";fileTuPTNames[3] = "data/tca_RH_UQuark_PT_truth";
	string fileTuEtaNames[4];fileTuEtaNames[0] = "data/tua_LH_UQuark_Eta_truth";fileTuEtaNames[1] = "data/tua_RH_UQuark_Eta_truth";fileTuEtaNames[2] = "data/tca_LH_UQuark_Eta_truth";fileTuEtaNames[3] = "data/tca_RH_UQuark_Eta_truth";
	string fileTuPhiNames[4];fileTuPhiNames[0] = "data/tua_LH_UQuark_Phi_truth";fileTuPhiNames[1] = "data/tua_RH_UQuark_Phi_truth";fileTuPhiNames[2] = "data/tca_LH_UQuark_Phi_truth";fileTuPhiNames[3] = "data/tca_RH_UQuark_Phi_truth";
	string fileTuMNames[4];fileTuMNames[0] = "data/tua_LH_UQuark_M_truth";fileTuMNames[1] = "data/tua_RH_UQuark_M_truth";fileTuMNames[2] = "data/tca_LH_UQuark_M_truth";fileTuMNames[3] = "data/tca_RH_UQuark_M_truth";



	for (int iFile = 0; iFile < 4; ++iFile)
	{
		
		TFile* file = new TFile(fileNames[iFile].c_str(),"READ");
		TTree* tree = (TTree*)file->Get("Delphes");


		TBranch *bP 		= tree->GetBranch("Particle");
		TBranch *bPhoton	= tree->GetBranch("Photon");
		TBranch *bJet   	= tree->GetBranch("Jet");
		TBranch *bElectron	= tree->GetBranch("Electron");
		TBranch *bMuon		= tree->GetBranch("Muon");
		TBranch *bMET		= tree->GetBranch("MissingET");



		TClonesArray *TCP 		= 0;
		TClonesArray *TCPhoton 	= 0;
		TClonesArray *TCJet 	= 0;
		TClonesArray *TCElectron= 0;
		TClonesArray *TCMuon 	= 0;
		TClonesArray *TCMET 	= 0;

		bP->SetAddress(&TCP);
		bPhoton->SetAddress(&TCPhoton);
		bJet->SetAddress(&TCJet);
		bElectron->SetAddress(&TCElectron);
		bMuon->SetAddress(&TCMuon);
		bMET->SetAddress(&TCMET);


		bP->GetEntry(0);
		bPhoton->GetEntry(0);
		bJet->GetEntry(0);
		bElectron->GetEntry(0);
		bMuon->GetEntry(0);
		bMET->GetEntry(0);

		bool truth = true;



		// bool untruth = true;
		// if (untruth)
		// {
		// 	int NPhotons = 0;
		// 	int numberOfPhotons = get_numberOfPhotons(fileNames[iFile].c_str());
		// 	VectorXd PhotonPT = VectorXd::Zero(numberOfPhotons);

		// 	int NbJets = 0;
		// 	int numberOfbJets = get_numberOfbJets(fileNames[iFile].c_str());
		// 	VectorXd bJet_PT = VectorXd::Zero(numberOfbJets);

		// 	VectorXd WE_PT = VectorXd::Zero(nEvents);

		// 	VectorXd W_PT = VectorXd::Zero(nEvents);
		// 	VectorXd W_M = VectorXd::Zero(nEvents);

		// 	VectorXd t_PT = VectorXd::Zero(nEvents);
		// 	VectorXd t_M = VectorXd::Zero(nEvents);


		// 	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		// 	{

		// 		bP->GetEntry(iEvent);
		// 		bPhoton->GetEntry(iEvent);
		// 		bJet->GetEntry(iEvent);
		// 		bElectron->GetEntry(iEvent);
		// 		bMuon->GetEntry(iEvent);
		// 		bMET->GetEntry(iEvent);

		// 		//####################################################
		// 		//##############Photon Zeugs##########################
		// 		//####################################################

		// 		int nPhotons = TCPhoton->GetEntries();
		// 		for (int iPhoton = 0; iPhoton < nPhotons; ++iPhoton)
		// 		{
		// 			Photon *P_Photon = (Photon*)TCPhoton->At(iPhoton);
		// 			PhotonPT(NPhotons) = P_Photon->PT;
		// 			NPhotons++;
		// 		}

		// 		//####################################################
		// 		//##############BJet Zeugs############################
		// 		//####################################################

		// 		TLorentzVector TV_bJet;
		// 		int nbtags = 0;
		// 		int nJets = TCJet->GetEntries();
		// 		for (int iJets = 0; iJets < nJets; ++iJets)
		// 		{
		// 		 	Jet	*P_Jet 	= (Jet*)TCJet->At(iJets);
		// 		 	if(P_Jet->BTag == 1){
		// 		 		nbtags++;
		// 		 		if(TV_bJet.Pt() < P_Jet->PT)
		// 		 			TV_bJet.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
		// 		 	}
		// 		}

		// 		if(nbtags == 1){
		// 			bJet_PT(NbJets) = TV_bJet.Pt();
		// 			NbJets++;
		// 		}else{
		// 			bJet_PT(NbJets) = 0;
		// 		}

		// 		//####################################################
		// 		//##############W Boson Reco Zeugs E##################
		// 		//####################################################
		// 		//Electronen
		// 		TLorentzVector TV_MET;
		// 		MissingET* P_MET = (MissingET*)TCMET->At(0);
		// 		TV_MET.SetPtEtaPhiE(P_MET->MET,P_MET->Eta,P_MET->Phi,P_MET->MET);


		// 		int nElectrons = TCElectron->GetEntries();

		// 		TLorentzVector TV_EElectron;
		// 		for (int iElectron = 0; iElectron < nElectrons; ++iElectron)
		// 		{
		// 			Electron* P_Electron = (Electron*)TCElectron->At(iElectron);
		// 			if (P_Electron->PT > TV_EElectron.Pt())
		// 			{
		// 				TV_EElectron.SetPtEtaPhiE(P_Electron->PT,P_Electron->Eta,P_Electron->Phi,P_Electron->PT);
		// 			}
		// 		}

		// 		WE_PT(iEvent) = (TV_EElectron+TV_MET).Pt();

		// 		//Richtig Reco
		// 		int nMuons = TCMuon->GetEntries();
		// 		double mW = 80.387;
		// 		TLorentzVector TV_W;
		// 		if(nJets >= 3 && nJets <= 6 && nMuons+nElectrons==0 && nbtags == 1)
		// 		{
		// 			TLorentzVector TV_qJet[5];
		// 			int nqjets = 0;
		// 			for (int iJets = 0; iJets < nJets; ++iJets)
		// 			{
		// 				Jet* P_Jet = (Jet*)TCJet->At(iJets);
		// 				if (P_Jet->BTag != 1)
		// 				{
		// 					TV_qJet[nqjets].SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
		// 					nqjets++;
		// 				}
		// 			}
		// 			if(nqjets+1 != nJets){cout << "ERROR: Irgenetwas stimmt mit deinen Jets nicht!!" << endl; continue;}
		// 			TV_W = permute_to_mass_reco(nqjets,TV_qJet,mW);
		// 			W_M(iEvent)	= TV_W.M();
		// 			W_PT(iEvent)= TV_W.Pt();
		// 		}
				
		// 		//####################################################
		// 		//##################Top Reco Zeugs E##################
		// 		//####################################################

		// 		TLorentzVector TV_top;
		// 		if(nJets >= 3 && nJets <= 6 && nMuons+nElectrons==0 && nbtags == 1)
		// 		{
		// 			TV_top = TV_W+TV_bJet;
		// 			t_M(iEvent)=TV_top.M();
		// 			t_PT(iEvent)=TV_top.Pt();
		// 		}
		// 	}
		// 	speichere(filePPTNames[iFile],PhotonPT);
		// 	speichere(filebJetPTNames[iFile],bJet_PT);
		// 	speichere(fileWPTENames[iFile],WE_PT);
		// 	speichere(fileWMNames[iFile],W_M);
		// 	speichere(fileWPTNames[iFile],W_PT);
		// 	speichere(filetMNames[iFile],t_M);
		// 	speichere(filetPTNames[iFile],t_PT);
		// }


		//####################################################
		//##################Truth Stuff#######################
		//####################################################



		VectorXd VTPhotonPT = VectorXd::Zero(nEvents);VectorXd VTPhotonEta = VectorXd::Zero(nEvents);VectorXd VTPhotonPhi = VectorXd::Zero(nEvents);
		VectorXd VTWBosonPT = VectorXd::Zero(nEvents);VectorXd VTWBosonEta = VectorXd::Zero(nEvents);VectorXd VTWBosonPhi = VectorXd::Zero(nEvents);VectorXd VTWBosonM = VectorXd::Zero(nEvents);
		VectorXd VTTopQuarkPT = VectorXd::Zero(nEvents);VectorXd VTTopQuarkEta = VectorXd::Zero(nEvents);VectorXd VTTopQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTTopQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTUQuarkPT = VectorXd::Zero(nEvents);VectorXd VTUQuarkEta = VectorXd::Zero(nEvents);VectorXd VTUQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTUQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTBQuarkPT = VectorXd::Zero(nEvents);VectorXd VTBQuarkEta = VectorXd::Zero(nEvents);VectorXd VTBQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTBQuarkM = VectorXd::Zero(nEvents);

		// for (int iEvent = 0; iEvent < 1 && truth; ++iEvent)
		// {
		// 	bP->GetEntry(iEvent);
		// 	GenParticle* P_Particle = (GenParticle*)TCP->At(0);
		// 	cout << P_Particle->PT << endl;

		// }


		for (int iEvent = 0; iEvent < nEvents && truth; ++iEvent)
		{
			bP->GetEntry(iEvent);

			int nParticles = TCP->GetEntries();

			VectorXd i_matched = vertex_match(TCP);
			if(i_matched(6)==0){continue;}


			GenParticle* P_Particle = (GenParticle*)TCP->At(i_matched(0));
			VTPhotonPT(iEvent) = P_Particle->PT;
			VTPhotonEta(iEvent) = P_Particle->Eta;
			VTPhotonPhi(iEvent) = P_Particle->Phi;

			P_Particle = (GenParticle*)TCP->At(i_matched(1));
			VTTopQuarkPT(iEvent) = P_Particle->PT;
			VTTopQuarkEta(iEvent) = P_Particle->Eta;
			VTTopQuarkPhi(iEvent) = P_Particle->Phi;
			VTTopQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(2));
			VTUQuarkPT(iEvent) = P_Particle->PT;
			VTUQuarkEta(iEvent) = P_Particle->Eta;
			VTUQuarkPhi(iEvent) = P_Particle->Phi;
			VTUQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(4));
			VTBQuarkPT(iEvent) = P_Particle->PT;
			VTBQuarkEta(iEvent) = P_Particle->Eta;
			VTBQuarkPhi(iEvent) = P_Particle->Phi;
			VTBQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(5));
			VTWBosonPT(iEvent) = P_Particle->PT;
			VTWBosonEta(iEvent) = P_Particle->Eta;
			VTWBosonPhi(iEvent) = P_Particle->Phi;
			VTWBosonM(iEvent) = P_Particle->Mass;
			if (iEvent%1000==0){ladebalken(iEvent,nEvents);}
		}


		speichere(fileTPPTNames[iFile],VTPhotonPT);
		speichere(fileTPEtaNames[iFile],VTPhotonEta);
		speichere(fileTPPhiNames[iFile],VTPhotonPhi);

		speichere(fileTWPTNames[iFile],VTWBosonPT);
		speichere(fileTWEtaNames[iFile],VTWBosonEta);
		speichere(fileTWPhiNames[iFile],VTWBosonPhi);
		speichere(fileTWMNames[iFile],VTWBosonM);

		speichere(fileTtPTNames[iFile],VTTopQuarkPT);
		speichere(fileTtEtaNames[iFile],VTTopQuarkEta);
		speichere(fileTtPhiNames[iFile],VTTopQuarkPhi);
		speichere(fileTtMNames[iFile],VTTopQuarkM);

		speichere(fileTbPTNames[iFile],VTBQuarkPT);
		speichere(fileTbEtaNames[iFile],VTBQuarkEta);
		speichere(fileTbPhiNames[iFile],VTBQuarkPhi);
		speichere(fileTbMNames[iFile],VTBQuarkM);

		speichere(fileTuPTNames[iFile],VTUQuarkPT);
		speichere(fileTuEtaNames[iFile],VTUQuarkEta);
		speichere(fileTuPhiNames[iFile],VTUQuarkPhi);
		speichere(fileTuMNames[iFile],VTUQuarkM);


		// int matched =0;
		// int not_matched = 0;
		// //for (int iEvent = 0; iEvent < 1; ++iEvent)
		// //{
		// 	int iEvent = 15950;
		// 	bP->GetEntry(iEvent);
		// 	VectorXd v = vertex_match(TCP);
		// 	cout << v << endl;
		// 	//ladebalken(iEvent+nEvents*iFile,nEvents*4);
		// 	if (v(4)==1){matched++;}
		// 	if (v(4)==0){not_matched++;}
		// 	// if (iEvent%10000 == 0)
		// 	// {
		// 	// 	cout << v << endl;
		// 	// }
		// //}

		// cout << endl;
		// cout << matched << "   " << not_matched << endl;

		file->Close();
	}

	return 0;

}




MatrixXd get_eventdisplay(string fileName, int event){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");

	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();
	MatrixXd display = MatrixXd::Zero(numberOfParticles,9);


    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{

		GenParticle *P = (GenParticle*)TCP->At(ipart);

		display(ipart,0)=P->PID;
		display(ipart,1)=P->Px;
		display(ipart,2)=P->Py;
		display(ipart,3)=P->Pz;
		display(ipart,4)=P->PT;
		display(ipart,5)=P->Mass;
		display(ipart,6)=P->E;
		display(ipart,7)=P->Eta;
		display(ipart,8)=P->Phi;
	}

	file->Close();
	return display;
}

MatrixXd get_eventdisplay_particle(string fileName, int event, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	int nPIDs = 0;
	for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID)nPIDs++;
	}
	if(nPIDs == 0){cout << "Kein Teilchen gefunden!" << endl; return MatrixXd::Zero(1,1);}
	MatrixXd display = MatrixXd::Zero(nPIDs,9);
	nPIDs = 0;
    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID){
		display(nPIDs,0)=P->PID;
		display(nPIDs,1)=P->Px;
		display(nPIDs,2)=P->Py;
		display(nPIDs,3)=P->Pz;
		display(nPIDs,4)=P->PT;
		display(nPIDs,5)=P->Mass;
		display(nPIDs,6)=P->E;
		display(nPIDs,7)=P->Eta;
		display(nPIDs,8)=P->Phi;
		nPIDs++;
		}
	}

	file->Close();
	return display;
}

int get_nTruth(string fileName, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(0);  
	int nEvents = branchP->GetEntries();

	int nPIDs = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		branchP->GetEntry(iEvent);  
		int nParticles = TCP->GetEntries();
		for (int ipart = 0; ipart < nParticles; ++ipart)
		{
		   	GenParticle *P = (GenParticle*)TCP->At(ipart);
			if(abs(P->PID) == PID)nPIDs++;
		}
	}

	file->Close();
	return nPIDs;
}


int get_nevents(string fileName){
	TFile* h_file = new TFile(fileName.c_str(),"READ");
	TTree* h_tree = (TTree*)h_file->Get("Delphes");
	int numberOfEntries = h_tree->GetEntries();
	h_file->Close();
	return numberOfEntries;

}


VectorXd get_event_MET(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	VectorXd MET = VectorXd::Zero(3);

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		MET(0)=MET(0)-P->Px;MET(1)=MET(1)-P->Py;
	}
	MET(2)=sqrt(pow(MET(0),2)+pow(MET(1),2));

	file->Close();
	return MET;
}


TLorentzVector get_event_MET_neutrino(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	TLorentzVector MET;

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		if (P->PID == 12 || P->PID == -12 || P->PID == 14 || P->PID == -14 || P->PID == 18 || P->PID == -18)
		{
			TLorentzVector nu(P->Px,P->Py,P->Pz,P->E);
			MET = MET+nu;
		}
	}
	MET.SetE(MET.Vect().Mag());

	file->Close();
	return MET;
}



int get_numberOfPhotons(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bPhoton = tree->GetBranch("Photon");
	TClonesArray* TCPhoton = 0;
	bPhoton->SetAddress(&TCPhoton);
	int nEvents = bPhoton->GetEntries();

	int nPhotons = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bPhoton->GetEntry(iEvent);
		nPhotons += TCPhoton->GetEntries();
	}
	file->Close();
	return nPhotons;
}

int get_numberOfJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();

	int nJet = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		nJet += TCJet->GetEntries();
	}
	file->Close();
	return nJet;
}

int get_numberOfbJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();
	int nbJet = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		int nJet = TCJet->GetEntries();
		for (int iJet = 0; iJet < nJet; ++iJet)
		{
			Jet* P_Jet = (Jet*)TCJet->At(iJet);
			if(P_Jet->BTag == 1)nbJet++;
		}
	}
	file->Close();
	return nbJet;
}


void ladebalken(int i, int max){
	double progress = (1.*i)/(1.*max)*100;
	#pragma omp critical
	std::cout << "\rSchritt " << i << " von " << max << " geschafft! " << "Also " << progress << " %";
	if(i == max-1 || i==max){std::cout << "\rSchritt " << max << " von " << max << " geschafft. Fertig!" << std::endl;}
	
	return;
}


void speichere(std::string name, MatrixXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}



void speichere(std::string name, VectorXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}




TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass)
{
	if(nVec == 1)
	{
		return vecs[0];
	}
	else if(nVec == 2)
	{
		return vecs[0]+vecs[1];
	}
	else if(nVec == 3)
	{
		TLorentzVector perms[3];
		double perm_m[3];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();
		perms[2] = vecs[1]+vecs[2];perm_m[2]=perms[2].M();
		if (fabs(perm_m[0]-mass) > fabs(perm_m[1]))
		{
			if (fabs(perm_m[1]-mass) > fabs(perm_m[2]))
			{
				return perms[2];
			}else
			{
				return perms[1];
			}
		}
		if (fabs(perm_m[0]-mass) > fabs(perm_m[2]))
		{
			return perms[2];
		}
		return perms[0];
	}
	else if(nVec == 4)
	{
		TLorentzVector perms[6];
		double perm_m[6];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[1]+vecs[2];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[3];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[2]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		for (int i = 0; i < 6; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%6]) && (perm_m[i] < perm_m[(i+2)%6]) && (perm_m[i] < perm_m[(i+3)%6]) && (perm_m[i] < perm_m[(i+4)%6]) && (perm_m[i] < perm_m[(i+5)%6]))
			{
				return perms[i];
			}
		}
	}
	else if(nVec == 5)
	{
		TLorentzVector perms[10];
		double perm_m[10];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[0]+vecs[4];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[2];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[1]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		perms[6] = vecs[1]+vecs[4];perm_m[6]=perms[6].M();perm_m[6]=fabs(perm_m[6]-mass);
		perms[7] = vecs[2]+vecs[3];perm_m[7]=perms[7].M();perm_m[7]=fabs(perm_m[7]-mass);
		perms[8] = vecs[2]+vecs[4];perm_m[8]=perms[8].M();perm_m[8]=fabs(perm_m[8]-mass);
		perms[9] = vecs[3]+vecs[4];perm_m[9]=perms[9].M();perm_m[9]=fabs(perm_m[9]-mass);
		for (int i = 0; i < 10; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%10]) && (perm_m[i] < perm_m[(i+2)%10]) && (perm_m[i] < perm_m[(i+3)%10]) && (perm_m[i] < perm_m[(i+4)%10]) && (perm_m[i] < perm_m[(i+5)%10]) && (perm_m[i] < perm_m[(i+6)%10]) && (perm_m[i] < perm_m[(i+7)%10]) && (perm_m[i] < perm_m[(i+8)%10]) && (perm_m[i] < perm_m[(i+9)%10]))
			{
				return perms[i];
			}
		}
	}
	else{
		cout << "ERROR: Not supported number of Jets" << endl;
		return vecs[0];
	}
}




VectorXd vertex_match(TClonesArray* TCP){
	int nPart = TCP->GetEntries();

	bool matched1=false;
	bool matched2=false;

	MatrixXd particle_infos(nPart,7);

	int nu = 0;int i_nu = 0;
	int nt = 0;int i_nt = 0; 
	int ng = 0;int i_ng = 0;
	int na = 0;int i_na = 0;

	for (int iPart = 0; iPart < nPart; ++iPart)
	{

		GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);
		GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(iPart);

		particle_infos(iPart,0) = P_Particle->PID;
		if(abs(P_Particle->M1) == 1){particle_infos(iPart,1) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M1));
			particle_infos(iPart,1) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->M2) == 1){particle_infos(iPart,2) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M2));
			particle_infos(iPart,2) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D1) == 1){particle_infos(iPart,3) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D1));
			particle_infos(iPart,3) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D2) == 1){particle_infos(iPart,4) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D2));
			particle_infos(iPart,4) = P_ParticleHelp->PID;
		}
		particle_infos(iPart,5) = iPart;
		particle_infos(iPart,6) = abs(P_Particle->PT);

		if(abs(particle_infos(iPart,0)) == 6)nt++;
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4)nu++;
		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21)ng++;
		if(particle_infos(iPart,0) == 22)na++;
	}

	MatrixXd u_infos = MatrixXd::Zero(nu,7);
	MatrixXd a_infos = MatrixXd::Zero(na,7);
	MatrixXd g_infos = MatrixXd::Zero(ng,7);
	MatrixXd t_infos = MatrixXd::Zero(nt,7);

	for (int iPart = 0; iPart < nPart; ++iPart)
	{
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
		if(abs(particle_infos(iPart,0)) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
		if(particle_infos(iPart,0) == 22){a_infos.row(i_na) = particle_infos.row(iPart);i_na++;}
	}


	//cout << a_infos << " " << na << "\n" << g_infos << " " << ng << "\n"  << u_infos << " " << nu << "\n"  << t_infos << " " << nt << "\n" ;


	VectorXd i_vertex_matched = VectorXd::Zero(7);
	bool matched_t = false;
	bool matched_b = false;
	GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(0);

	for (int iu = 0; iu < i_nu; ++iu)
	{
		if (u_infos(iu,3) == 22 && abs(u_infos(iu,4))==6)
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(1) = P_ParticleHelp->D2;
			i_vertex_matched(0) = P_ParticleHelp->D1;
			matched_t=true;
		}
		if(u_infos(iu,4) == 22 && abs(u_infos(iu,3))==6)
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(1) = P_ParticleHelp->D1;
			i_vertex_matched(0) = P_ParticleHelp->D2;
			matched_t=true;
		}
	}
	for (int it = 0; it < nt; ++it)
	{
		if ( abs(t_infos(it,3))==5 && abs(t_infos(it,4))==24 )   
		{
			P_ParticleHelp = (GenParticle*)TCP->At(t_infos(it,5));
			i_vertex_matched(4) = P_ParticleHelp->D1;
			i_vertex_matched(5) = P_ParticleHelp->D2;
			matched_b = true;
		}
		if ( abs(t_infos(it,4))==5 && abs(t_infos(it,3))==24 )  
		{
			i_vertex_matched(1) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(t_infos(it,5));
			i_vertex_matched(4) = P_ParticleHelp->D2;
			i_vertex_matched(5) = P_ParticleHelp->D1;
			matched_b = true;
		}
	}

	if (matched_t && matched_b){i_vertex_matched(6)=1;}else{i_vertex_matched(6)=0;}

	return i_vertex_matched;
}


TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[])
{
	if(error < 2){error++;
	cout << "INFO: No Mass was given choosing highest PT" << endl;}
	if (nVec < 2)
	{	
		if(error < 2){error++;
		cout << "INFO: Only one Vector given -> Outout = Input" << endl;}
		return vecs[0];
	}

	TLorentzVector reco[2];
	TLorentzVector swap;
	reco[0] = vecs[0];
	reco[1] = vecs[1];
	for (int i = 2; i < nVec; ++i)
	{
		if (reco[0].Pt() < reco[1].Pt())
		{
			swap = reco[1]; reco[1]=reco[0]; reco[0]=swap;
		}
		if (vecs[i].Pt() > reco[1].Pt())
		{
			reco[1] = vecs[i];
		}
	}
	swap = reco[0]+reco[1];
	return swap;
}




// bP->GetEntry(iEvent);  
// 		bMuon->GetEntry(iEvent);  
// 		bMet->GetEntry(iEvent);  
// 		bJet->GetEntry(iEvent);  
	
		
// 		Long64_t nParticles = TCP->GetEntries();
// 		Long64_t nMuons 	= TCMuon->GetEntries();
// 		Long64_t nJets 		= TCJet->GetEntries();

// 		GenParticle *P 		= (GenParticle*)TCP->At(0);
// 		MissingET 	*P_MET 	= (MissingET*)TCMet->At(0);
// 		Muon 		*P_Muon = (Muon*)TCMuon->At(0);
// 		Jet 		*P_Jet 	= (Jet*)TCJet->At(0);

// 		TLorentzVector TV_Met;
// 		TLorentzVector TV_Muon;
// 		TLorentzVector TV_bJet1;
// 		TLorentzVector TV_bJet2;
// 		TLorentzVector TV_qJet1;
// 		TLorentzVector TV_qJet2;


// 		if(nJets != 4)continue;

// 		int nbtags=0;
// 		int nq=0;
// 		for (int iJets = 0; iJets < nJets; ++iJets)
// 		{
// 			Jet	*P_Jet 	= (Jet*)TCJet->At(iJets);
// 			if (P_Jet->BTag == 1)nbtags++;
// 			else nq++;

// 			if(P_Jet->BTag == 1 && nbtags == 1)
// 				TV_bJet1.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
// 			if(P_Jet->BTag == 1 && nbtags == 2)
// 				TV_bJet2.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
// 			if(P_Jet->BTag == 0 && nq == 1)
// 				TV_qJet1.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
// 			if(P_Jet->BTag == 0 && nq == 2)
// 				TV_qJet2.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
			
// 		}
// 		if(nbtags != 2)continue;

// 		TLorentzVector TV_W = TV_qJet1+TV_qJet2;
// 		mW(sel) = TV_W.M();

// 		double mT=173;
// 		double RmT1=0,RmT2=0;

// 		RmT1 = (TV_W+TV_bJet1).M();
// 		RmT2 = (TV_W+TV_bJet2).M();

// 		if (fabs(RmT1 - mT)<fabs(RmT2-mT)){
// 			mT1(sel)=RmT1;
// 			mT2(sel)=RmT2;
// 		}
// 		else{
// 			mT1(sel)=RmT2;
// 			mT2(sel)=RmT1;
// 		}

// 		sel++;
// 	}

// 	speichere("QQ_INV_MASS",mW);
// 	speichere("WB1_INV_MASS",mT1);
// 	speichere("WB2_INV_MASS",mT2);



// TLorentzVector TV_Lepton;
			// TLorentzVector TV_W;
			// int nMuons = TCMuon->GetEntries();
			// if(nMuons+nElectrons == 1){
			// 	TV_Lepton = TV_EElectron;
			// 	for (int iMuon = 0; iMuon < nMuons; ++iMuon)
			// 	{
			// 		Muon* P_Muon = (Muon*)TCMuon->At(iMuon);
			// 		if (P_Muon->PT > TV_Lepton.Pt())
			// 		{
			// 			TV_Lepton.SetPtEtaPhiE(P_Muon->PT,P_Muon->Eta,P_Muon->Phi,P_Muon->PT);
			// 		}
			// 	}
			// 	TV_W = TV_Lepton+TV_MET;

			// 	W_M(iEvent) = TV_W.M();
			// 	W_PT(iEvent) = TV_W.Pt();
			//}


// for (int iElectron = 0; iElectron < nElectrons; ++iElectron)
// 				{
// 					Electron* P_Electron = (Electron*)TCElectron->At(iElectron);
// 					TV_Current_Lepton.SetPtEtaPhiE(P_Electron->PT,P_Electron->Eta,P_Electron->Phi,P_Electron->PT);
// 					if (fabs((TV_Lepton+TV_MET).M()-mW) > fabs((TV_Current_Lepton+TV_MET).M()-mW) && (TV_Current_Lepton+TV_MET).M() > 0)
// 					{
// 						TV_Lepton = TV_Current_Lepton;
// 					}
// 				}
// 				for (int iMuon = 0; iMuon < nMuons; ++iMuon)
// 				{
// 					Muon* P_Muon = (Muon*)TCMuon->At(iMuon);
// 					TV_Current_Lepton.SetPtEtaPhiE(P_Muon->PT,P_Muon->Eta,P_Muon->Phi,P_Muon->PT);
// 					if (fabs((TV_Lepton+TV_MET).M()-mW) > fabs((TV_Current_Lepton+TV_MET).M()-mW) && (TV_Current_Lepton+TV_MET).M() > 0)
// 					{
// 						TV_Lepton = TV_Current_Lepton;
// 					}
// 				}
// 				TV_W = TV_Lepton+TV_MET;

// 				W_M(iEvent) = TV_W.M();
// 				W_PT(iEvent) = TV_W.Pt();

















// VectorXd vertex_match(TClonesArray* TCP){
// 	int npart_sel = 20;
// 	int nPart = TCP->GetEntries();

// 	bool matched1=false;
// 	bool matched2=false;

// 	MatrixXd particle_infos(npart_sel,6);

// 	int nu = 0;int i_nu = 0;
// 	int nt = 0; int i_nt = 0; 
// 	int ng = 0;int i_ng = 0;
// 	int na = 0;int i_na = 0;
// 	for (int iPart = 0; iPart < npart_sel; ++iPart)
// 	{
// 		GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);

// 		particle_infos(iPart,0) = abs(P_Particle->PID);
// 		particle_infos(iPart,1) = abs(P_Particle->M1);
// 		particle_infos(iPart,2) = abs(P_Particle->M2);
// 		particle_infos(iPart,3) = abs(P_Particle->D1);
// 		particle_infos(iPart,4) = abs(P_Particle->D2);
// 		particle_infos(iPart,5) = iPart;

// 		if(particle_infos(iPart,0) == 2)nu++;
// 		if(particle_infos(iPart,0) == 6)nt++;
// 		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21)ng++;
// 		if(particle_infos(iPart,0) == 22)na++;
// 	}

// 	MatrixXd u_infos = MatrixXd::Zero(nu,6);
// 	MatrixXd a_infos = MatrixXd::Zero(na,6);
// 	MatrixXd g_infos = MatrixXd::Zero(ng,6);
// 	MatrixXd t_infos = MatrixXd::Zero(nt,6);

// 	for (int iPart = 0; iPart < npart_sel; ++iPart)
// 	{
// 		if(particle_infos(iPart,0) == 2){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
// 		if(particle_infos(iPart,0) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
// 		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
// 		if(particle_infos(iPart,0) == 22){a_infos.row(i_na) = particle_infos.row(iPart);i_na++;}
// 	}


// 	//cout << a_infos << "\n" << g_infos << "\n"  << u_infos << "\n"  << t_infos << "\n" ;
// 	VectorXd i_vertex_matched = VectorXd::Zero(5);
// 	for (int it = 0; it < nt; ++it)
// 	{
// 		for (int ia = 0; ia < na; ++ia)
// 		{
// 			if (a_infos(ia,1)==t_infos(it,1) && a_infos(ia,2)==t_infos(it,2))
// 			{
// 				i_vertex_matched(0)=t_infos(it,5);
// 				i_vertex_matched(1)=a_infos(ia,5);
// 				matched1=true;
// 			}
// 		}
// 	}
// 	for (int ig = 0; ig < ng; ++ig)
// 	{
// 		for (int iu = 0; iu < nu; ++iu)
// 		{
// 			if (g_infos(ig,3)==u_infos(iu,3) && g_infos(ig,4)==u_infos(iu,4))
// 			{
// 				i_vertex_matched(2)=u_infos(iu,5);
// 				i_vertex_matched(3)=g_infos(ig,5);
// 				matched2=true;
// 			}
// 		}
// 	}
// 	if (matched1 && matched2){i_vertex_matched(4)=1;}else{i_vertex_matched(4)=0;}

// 	return i_vertex_matched;
// }