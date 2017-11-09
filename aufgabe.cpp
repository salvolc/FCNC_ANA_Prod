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

	string filePPTNames[4];
	filePPTNames[0] = "data/tua_LH_Photon_PT";
	filePPTNames[1] = "data/tua_RH_Photon_PT";
	filePPTNames[2] = "data/tca_LH_Photon_PT";
	filePPTNames[3] = "data/tca_RH_Photon_PT";

	string filebJetPTNames[4];	filebJetPTNames[0] = "data/tua_LH_bJet_PT";filebJetPTNames[1] = "data/tua_RH_bJet_PT";filebJetPTNames[2] = "data/tca_LH_bJet_PT";filebJetPTNames[3] = "data/tca_RH_bJet_PT";
	string fileWPTENames[4];	fileWPTENames[0] = "data/tua_LH_WE_PT";fileWPTENames[1] = "data/tua_RH_WE_PT";fileWPTENames[2] = "data/tca_LH_WE_PT";fileWPTENames[3] = "data/tca_RH_WE_PT";
	string fileWPTNames[4];		fileWPTNames[0] = "data/tua_LH_W_PT";fileWPTNames[1] = "data/tua_RH_W_PT";fileWPTNames[2] = "data/tca_LH_W_PT";fileWPTNames[3] = "data/tca_RH_W_PT";
	string fileWMNames[4];		fileWMNames[0] = "data/tua_LH_W_M";fileWMNames[1] = "data/tua_RH_W_M";fileWMNames[2] = "data/tca_LH_W_M";fileWMNames[3] = "data/tca_RH_W_M";
	string filetPTNames[4];		filetPTNames[0] = "data/tua_LH_t_PT";filetPTNames[1] = "data/tua_RH_t_PT";filetPTNames[2] = "data/tca_LH_t_PT";filetPTNames[3] = "data/tca_RH_t_PT";
	string filetMNames[4];		filetMNames[0] = "data/tua_LH_t_M";filetMNames[1] = "data/tua_RH_t_M";filetMNames[2] = "data/tca_LH_t_M";filetMNames[3] = "data/tca_RH_t_M";
	

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



		int NPhotons = 0;
		int numberOfPhotons = get_numberOfPhotons(fileNames[iFile].c_str());
		VectorXd PhotonPT = VectorXd::Zero(numberOfPhotons);

		int NbJets = 0;
		int numberOfbJets = get_numberOfbJets(fileNames[iFile].c_str());
		VectorXd bJet_PT = VectorXd::Zero(numberOfbJets);

		VectorXd WE_PT = VectorXd::Zero(nEvents);

		VectorXd W_PT = VectorXd::Zero(nEvents);
		VectorXd W_M = VectorXd::Zero(nEvents);

		VectorXd t_PT = VectorXd::Zero(nEvents);
		VectorXd t_M = VectorXd::Zero(nEvents);


		for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		{

			bP->GetEntry(iEvent);
			bPhoton->GetEntry(iEvent);
			bJet->GetEntry(iEvent);
			bElectron->GetEntry(iEvent);
			bMuon->GetEntry(iEvent);
			bMET->GetEntry(iEvent);

			//####################################################
			//##############Photon Zeugs##########################
			//####################################################

			int nPhotons = TCPhoton->GetEntries();
			for (int iPhoton = 0; iPhoton < nPhotons; ++iPhoton)
			{
				Photon *P_Photon = (Photon*)TCPhoton->At(iPhoton);
				PhotonPT(NPhotons) = P_Photon->PT;
				NPhotons++;
			}

			//####################################################
			//##############BJet Zeugs############################
			//####################################################

			TLorentzVector TV_bJet;
			int nbtags = 0;
			int nJets = TCJet->GetEntries();
			for (int iJets = 0; iJets < nJets; ++iJets)
			{
			 	Jet	*P_Jet 	= (Jet*)TCJet->At(iJets);
			 	if(P_Jet->BTag == 1){
			 		nbtags++;
			 		if(TV_bJet.Pt() < P_Jet->PT)
			 			TV_bJet.SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
			 	}
			}

			if(nbtags == 1){
				bJet_PT(NbJets) = TV_bJet.Pt();
				NbJets++;
			}else{
				bJet_PT(NbJets) = 0;
			}

			//####################################################
			//##############W Boson Reco Zeugs E##################
			//####################################################
			//Electronen
			TLorentzVector TV_MET;
			MissingET* P_MET = (MissingET*)TCMET->At(0);
			TV_MET.SetPtEtaPhiE(P_MET->MET,P_MET->Eta,P_MET->Phi,P_MET->MET);


			int nElectrons = TCElectron->GetEntries();

			TLorentzVector TV_EElectron;
			for (int iElectron = 0; iElectron < nElectrons; ++iElectron)
			{
				Electron* P_Electron = (Electron*)TCElectron->At(iElectron);
				if (P_Electron->PT > TV_EElectron.Pt())
				{
					TV_EElectron.SetPtEtaPhiE(P_Electron->PT,P_Electron->Eta,P_Electron->Phi,P_Electron->PT);
				}
			}

			WE_PT(iEvent) = (TV_EElectron+TV_MET).Pt();

			//Richtig Reco
			int nMuons = TCMuon->GetEntries();
			double mW = 80.387;
			TLorentzVector TV_W;
			if(nJets >= 3 && nJets <= 6 && nMuons+nElectrons==0 && nbtags == 1)
			{
				TLorentzVector TV_qJet[5];
				int nqjets = 0;
				for (int iJets = 0; iJets < nJets; ++iJets)
				{
					Jet* P_Jet = (Jet*)TCJet->At(iJets);
					if (P_Jet->BTag != 1)
					{
						TV_qJet[nqjets].SetPtEtaPhiM(P_Jet->PT,P_Jet->Eta,P_Jet->Phi,P_Jet->Mass);
						nqjets++;
					}
				}
				if(nqjets+1 != nJets){cout << "ERROR: Irgenetwas stimmt mit deinen Jets nicht!!" << endl; continue;}
				TV_W = permute_to_mass_reco(nqjets,TV_qJet,mW);
				W_M(iEvent)	= TV_W.M();
				W_PT(iEvent)= TV_W.Pt();
			}
			
			//####################################################
			//##################Top Reco Zeugs E##################
			//####################################################

			TLorentzVector TV_top;
			if(nJets >= 3 && nJets <= 6 && nMuons+nElectrons==0 && nbtags == 1)
			{
				TV_top = TV_W+TV_bJet;
				t_M(iEvent)=TV_top.M();
				t_PT(iEvent)=TV_top.Pt();
			}







		}

		file->Close();

		speichere(filePPTNames[iFile],PhotonPT);
		speichere(filebJetPTNames[iFile],bJet_PT);
		speichere(fileWPTENames[iFile],WE_PT);
		speichere(fileWMNames[iFile],W_M);
		speichere(fileWPTNames[iFile],W_PT);
		speichere(filetMNames[iFile],t_M);
		speichere(filetPTNames[iFile],t_PT);
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