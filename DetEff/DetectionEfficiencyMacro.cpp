/////////////////////////////////////////////////////////////////////
//                                                                 //
//  DetectionEfficiencyMacro.cpp                                   //
//                                                                 //
//  Jenna Chisholm                                                 //
//  August 19, 2020                                                //
//                                                                 //
//  Last Updated: August 20, 2020                                  //
//                                                                 //
//  Reads in simulated data from EvGen and the same data after it  //
//  is passed through A2Geant4 and Ant, and calculates a detector  //
//  efficiency.                                                    //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "physics.h"
#include "math.h"
using namespace std;


// A function that calculates the detector efficiency at a given incident photon energy and writes the histogram
// Input: incident photon energy (currently accepts: 150,175,200,225,250,275,300,325,350)
void GetDetEff(int ke)
{
	// Open files (quit program if they cannot be opened)
	
	// Input file (made by EvGen)
	TFile *evgenFile = TFile::Open(Form("~/opt/EvGen/out/5cm/pi0_he4_%d_in.root",ke));
	if(!evgenFile->IsOpen())
	{
		cout << "Error: EvGen file could not be opened.";
		exit(-1);
	}

	// Output file (from A2Geant4 and Ant)
	TFile *antFile = TFile::Open(Form("MC_for_Edet/he4pi0_%d.root",ke));
	if(!antFile->IsOpen())
	{
		cout << "Error: Ant file could not be opened.";
		exit(-1);
	}

	// Get histograms
	TH1D *hThetaIn = (TH1D*)evgenFile->Get("h4");
	TH3D *h3DOut = (TH3D*)antFile->Get("He4Pi0/h3D_MEpi0");

	// Project the Ant histogram so it shows just as a function of theta
        TH1D *hThetaOut = h3DOut->ProjectionY("hThetaOut");

	// Make the detection efficiency histogram
	hThetaIn->Rebin(10);  // if I want bigger bins I'll need to change Ant
	TH1D *hDetEff = (TH1D*)hThetaOut->Clone("hDetEff");
	hDetEff->Divide(hThetaIn);

	// Draw the detection efficiency histogram
	TCanvas *c1 = new TCanvas("c1","",200,10,1000,500);
	c1->Divide(3,1);
	c1->cd(1);
	hThetaIn->Draw();
	c1->cd(2);
	hThetaOut->Draw();
	c1->cd(3);
	hDetEff->SetTitle(Form("Detection Efficiency at %d MeV",ke));
	hDetEff->GetYaxis()->SetTitle("Detection Efficiency");
	hDetEff->Draw();

	// Write and save detector efficiency
	TFile *resultFile = new TFile(Form("DetEff_%d.root",ke),"RECREATE");
	hThetaIn->Write();
	hThetaOut->Write();
	hDetEff->Write();
	resultFile->Close();

	evgenFile->Close();
	antFile->Close();

}

