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

#include "functions.hpp"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;

VectorXd vertex_match(TClonesArray* TCP);


int main(int argc, char const *argv[])
{
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/MG5_aMC_v2_6_0/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MadGraph/Delphes-3.4.1/libDelphes.so");


	string 	fileNames[4];
			fileNames[0]="samples/tua_LH_PY8_ATLASDELPHES_50000.root";
			fileNames[1]="samples/tua_RH_PY8_ATLASDELPHES_50000.root";
			fileNames[2]="samples/tca_LH_PY8_ATLASDELPHES_50000.root";
			fileNames[3]="samples/tca_RH_PY8_ATLASDELPHES_50000.root";

	TFile* file = new TFile(fileNames[0].c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");


	TBranch *bP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	bP->SetAddress(&TCP);

	int nEvents = get_nevents(fileNames[2].c_str());	
	VectorXd i_matched = VectorXd::Zero(4);

	VectorXd PPT = VectorXd::Zero(nEvents);

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bool matched1 = false;
		bool matched2 = false;
		bP->GetEntry(iEvent);

		i_matched = vertex_match(TCP);

		if(i_matched(4)==0)continue;

		GenParticle* P_Particle = (GenParticle*)TCP->At(i_matched(1));
		PPT(iEvent) = P_Particle->PT;

	}
	cout << i_matched << endl;
	speichere("PPT",PPT);
	return 0;
}





VectorXd vertex_match(TClonesArray* TCP){
	int npart_sel = 20;
	int nPart = TCP->GetEntries();

	bool matched1=false;
	bool matched2=false;

	MatrixXd particle_infos(npart_sel,6);

	int nu = 0;int i_nu = 0;
	int nt = 0; int i_nt = 0; 
	int ng = 0;int i_ng = 0;
	int na = 0;int i_na = 0;
	for (int iPart = 0; iPart < npart_sel; ++iPart)
	{
		GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);

		particle_infos(iPart,0) = abs(P_Particle->PID);
		particle_infos(iPart,1) = abs(P_Particle->M1);
		particle_infos(iPart,2) = abs(P_Particle->M2);
		particle_infos(iPart,3) = abs(P_Particle->D1);
		particle_infos(iPart,4) = abs(P_Particle->D2);
		particle_infos(iPart,5) = iPart;

		if(particle_infos(iPart,0) == 2)nu++;
		if(particle_infos(iPart,0) == 6)nt++;
		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21)ng++;
		if(particle_infos(iPart,0) == 22)na++;
	}

	MatrixXd u_infos = MatrixXd::Zero(nu,6);
	MatrixXd a_infos = MatrixXd::Zero(na,6);
	MatrixXd g_infos = MatrixXd::Zero(ng,6);
	MatrixXd t_infos = MatrixXd::Zero(nt,6);

	for (int iPart = 0; iPart < npart_sel; ++iPart)
	{
		if(particle_infos(iPart,0) == 2){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
		if(particle_infos(iPart,0) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
		if(particle_infos(iPart,0) == 22){a_infos.row(i_na) = particle_infos.row(iPart);i_na++;}
	}


	//cout << a_infos << "\n" << g_infos << "\n"  << u_infos << "\n"  << t_infos << "\n" ;
	VectorXd i_vertex_matched = VectorXd::Zero(5);
	for (int it = 0; it < nt; ++it)
	{
		for (int ia = 0; ia < na; ++ia)
		{
			if (a_infos(ia,1)==t_infos(it,1) && a_infos(ia,2)==t_infos(it,2))
			{
				i_vertex_matched(0)=t_infos(it,5);
				i_vertex_matched(1)=a_infos(ia,5);
				matched1=true;
			}
		}
	}
	for (int ig = 0; ig < ng; ++ig)
	{
		for (int iu = 0; iu < nu; ++iu)
		{
			if (g_infos(ig,3)==u_infos(iu,3) && g_infos(ig,4)==u_infos(iu,4))
			{
				i_vertex_matched(2)=u_infos(iu,5);
				i_vertex_matched(3)=g_infos(ig,5);
				matched2=true;
			}
		}
	}
	if (matched1 && matched2){i_vertex_matched(4)=1;}else{i_vertex_matched(4)=0;}

	return i_vertex_matched;
}








		// int nPart = TCP->GetEntries();


		// MatrixXd particle_infos(npart_sel,6);


		// int nu = 0;int i_nu = 0;
		// int nt = 0; int i_nt = 0; 
		// int ng = 0;int i_ng = 0;
		// int na = 0;int i_na = 0;
		// for (int iPart = 0; iPart < npart_sel; ++iPart)
		// {
		// 	GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);

		// 	particle_infos(iPart,0) = abs(P_Particle->PID);
		// 	particle_infos(iPart,1) = abs(P_Particle->M1);
		// 	particle_infos(iPart,2) = abs(P_Particle->M2);
		// 	particle_infos(iPart,3) = abs(P_Particle->D1);
		// 	particle_infos(iPart,4) = abs(P_Particle->D2);
		// 	particle_infos(iPart,5) = iPart;

		// 	if(particle_infos(iPart,0) == 2)nu++;
		// 	if(particle_infos(iPart,0) == 6)nt++;
		// 	if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21)ng++;
		// 	if(particle_infos(iPart,0) == 22)na++;
		// }

		// MatrixXd u_infos = MatrixXd::Zero(nu,6);
		// MatrixXd a_infos = MatrixXd::Zero(na,6);
		// MatrixXd g_infos = MatrixXd::Zero(ng,6);
		// MatrixXd t_infos = MatrixXd::Zero(nt,6);

		// for (int iPart = 0; iPart < npart_sel; ++iPart)
		// {
		// 	if(particle_infos(iPart,0) == 2){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
		// 	if(particle_infos(iPart,0) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		// 	if(particle_infos(iPart,0) == 9 || particle_infos(iPart,0) == 21){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
		// 	if(particle_infos(iPart,0) == 22){a_infos.row(i_na) = particle_infos.row(iPart);i_na++;}
		// }


		// //cout << a_infos << "\n" << g_infos << "\n"  << u_infos << "\n"  << t_infos << "\n" ;
		// VectorXd i_vertex_matched = VectorXd::Zero(4);
		// for (int it = 0; it < nt; ++it)
		// {
		// 	for (int ia = 0; ia < na; ++ia)
		// 	{
		// 		if (a_infos(ia,1)==t_infos(it,1) && a_infos(ia,2)==t_infos(it,2))
		// 		{
		// 			i_vertex_matched(0)=t_infos(it,5);
		// 			i_vertex_matched(1)=a_infos(ia,5);
		// 			matched1=true;
		// 		}
		// 	}
		// }
		// for (int ig = 0; ig < ng; ++ig)
		// {
		// 	for (int iu = 0; iu < nu; ++iu)
		// 	{
		// 		if (g_infos(ig,3)==u_infos(iu,3) && g_infos(ig,4)==u_infos(iu,4))
		// 		{
		// 			i_vertex_matched(2)=u_infos(iu,5);
		// 			i_vertex_matched(3)=g_infos(ig,5);
		// 			matched2=true;
		// 		}
		// 	}
		// }
		// if (matched1 && matched2){matched++;}









// if (abs(P_Particle->PID) == 6 && P_Particle->IsPU == 0)
// {
// 	cout << "Top Quark : \n";
// 	cout << "Mother 1: " << P_Particle->M1 << endl;
// 	cout << "Mother 2: " << P_Particle->M2 << "\n" << endl;
// 	//cout << "Daugther 1: " << P_Particle->D1 << endl;
// 	//cout << "Daugther 2: " << P_Particle->D2 << "\n" <<endl;
// }
// if (abs(P_Particle->PID) == 2 && P_Particle->IsPU == 0)
// {
// 	cout << "Up Quark : \n";
// 	cout << "Mother 1: " << P_Particle->M1 << endl;
// 	cout << "Mother 2: " << P_Particle->M2 << endl;
//  	cout << "Daugther 1: " << P_Particle->D1 << endl;
//  	cout << "Daugther 2: " << P_Particle->D2 << "\n" <<endl;
// }
// if (abs(P_Particle->PID) == 22 && P_Particle->IsPU == 0)
// {
// 	cout << "Photon : \n";
// 	cout << "Mother 1: " << P_Particle->M1 << endl;
// 	cout << "Mother 2: " << P_Particle->M2 << "\n" << endl;
// 	// cout << "Daugther 1: " << P_Particle->D1 << endl;
// 	// cout << "Daugther 2: " << P_Particle->D2 << "\n" <<endl;
// }
// if ((abs(P_Particle->PID) == 9 || abs(P_Particle->PID) == 21) && P_Particle->IsPU == 0)
// {
// 	cout << "Gluon : \n";
// 	cout << "Mother 1: " << P_Particle->M1 << endl;
// 	cout << "Mother 2: " << P_Particle->M2 << endl;
// 	cout << "Daugther 1: " << P_Particle->D1 << endl;
// 	cout << "Daugther 2: " << P_Particle->D2 << "\n" <<endl;
// }
// int PPID = 	abs(P_Particle->PID);
// int MID = 	abs(P_Particle->M1);
// int MIID = 	abs(P_Particle->M2);
// int DID = 	abs(P_Particle->D1);
// int DIID = 	abs(P_Particle->D2);

// if(PPID == 22)j++;
// if(PPID==2 && (DIID==22||DID==6) && (DID==22||DIID==6)){
// 	matched++;
// }
// 	if(PPID==2 && (DIID==22||DID==8) && (DID==22||DIID==8)){
// 		matched++;
// 	}
// 	if(PPID==4 && (DIID==22||DID==8) && (DID==22||DIID==8)){
// 		matched++;
// 	}
// 	if(PPID==2 && (DIID==22||DID==6) && (DID==22||DIID==6)){
// 	matched++;
// }