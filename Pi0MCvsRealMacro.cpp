/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Pi0MCvsRealMacro.cpp                                           //
//                                                                 //
//  Jenna Chisholm                                                 //
//  August 6, 2020                                                 //
//                                                                 //
//  Last Updated: August 20, 2020                                  //
//                                                                 //
//  Reads in simulated Pi0 missing energy histograms for coherent, //
//  QF P, and QF N, and adds them together, scaling as best as it  //
//  can to match the real data Pi0 missing energy histogram.       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "physics.h"
#include "math.h"
using namespace std;


// A function which adds the three MC histograms together, scales them, and
// compares them to the real data.
void Compare()
{
	// Open files (quit program if they cannot be opened)
	
	// Coherent (elastic) file
	TFile *cFile = TFile::Open("MC/G4He_4HePi0_300_result3.root");
	if(!cFile->IsOpen())
	{
		cout << "Error: coherent file could not be opened.";
		exit(-1);
	}

	// Quasi-free from proton file
	TFile *pFile = TFile::Open("MC/G4He_Pi0P3H_QF_300_result3.root");
	if(!pFile->IsOpen())
	{
		cout << "Error: P QF file could not be opened.";
		exit(-1);
	}

	// Quasi-free from neutron file
	TFile *nFile = TFile::Open("MC/G4He_Pi0N3He_QF_300_result3.root");
	if(!nFile->IsOpen())
	{
		cout << "Error: N QF file could not be opened.";
		exit(-1);
	}

	// Real data file
	TFile *dFile = TFile::Open("pi0_295-305MeV_result.root");
	if(!dFile->IsOpen())
	{
		cout << "Error: data file could not be opened.";
		exit(-1);
	}

	// Get Pi0 missing energy histograms
	TH1F *hCoherent = (TH1F*)cFile->Get("He4Pi0/h_MEpi0");
	TH1F *hPQF = (TH1F*)pFile->Get("He4Pi0/h_MEpi0");
	TH1F *hNQF = (TH1F*)nFile->Get("He4Pi0/h_MEpi0");	
	TH1F *hReal = (TH1F*)dFile->Get("He4Pi0/h_MMpi0_2");

	// Add together the three simulated files and scale them to match the real data
	TH1F *hTotal = (TH1F*) hCoherent->Clone("hTotal");
	//hTotal->Scale(0.75);  // this was based on eying the histograms
	hTotal->Add(hPQF);
	hTotal->Add(hNQF);
	//hTotal->Scale(0.45);  // this was based on eying the histograms
	
	// Create a canvas and draw the real and simulated data together
	TCanvas *c1 = new TCanvas("c1","",200,10,1000,750);
	hTotal->SetLineColor(kRed);
	hTotal->Draw();
	hReal->SetLineColor(kBlue);
	hReal->Draw("SAME");

	// Create a legend for the histogram
	TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
	legend->AddEntry(hTotal, "MC Data","l");
	legend->AddEntry(hReal, "Real Data","l");
	legend->Draw();


}

