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
	TFile *cFile = TFile::Open("MC/MC_4HePi0_300_Final.root");
	if(!cFile->IsOpen())
	{
		cout << "Error: coherent file could not be opened.";
		exit(-1);
	}

	// Quasi-free from proton file
	TFile *pFile = TFile::Open("MC/MC_Pi0P3H_QF_300_Final.root");
	if(!pFile->IsOpen())
	{
		cout << "Error: P QF file could not be opened.";
		exit(-1);
	}

	// Quasi-free from neutron file
	TFile *nFile = TFile::Open("MC/MC_Pi0N3He_QF_300_Final.root");
	if(!nFile->IsOpen())
	{
		cout << "Error: N QF file could not be opened.";
		exit(-1);
	}

	// Real data file
	TFile *dFile = TFile::Open("He4Pi0_295-305_FinalB.root");
	if(!dFile->IsOpen())
	{
		cout << "Error: data file could not be opened.";
		exit(-1);
	}

	// Get Pi0 missing energy histograms
	TH1F *hCoherent = (TH1F*)cFile->Get("He4Pi0/h_MEpi0");
	TH1F *hPQF = (TH1F*)pFile->Get("He4Pi0/h_MEpi0");
	TH1F *hNQF = (TH1F*)nFile->Get("He4Pi0/h_MEpi0");	
	TH1F *hReal = (TH1F*)dFile->Get("He4Pi0/h_MEpi0");

	// Create scaled versions of the MC data
	TH1F* hCoherent_Scaled = (TH1F*)hCoherent->Clone("hCoherent_Scaled");
	TH1F* hPQF_Scaled = (TH1F*)hPQF->Clone("hPQF_Scaled");
	TH1F* hNQF_Scaled = (TH1F*)hNQF->Clone("hNQF_Scaled");
	hCoherent_Scaled->Scale(0.98);
	hPQF_Scaled->Scale(1.28);
	hNQF_Scaled->Scale(0.8);

	// Add together the three scaled simulated files
	TH1F *hTotal = (TH1F*) hCoherent_Scaled->Clone("hTotal");
	hTotal->Add(hPQF_Scaled);
	hTotal->Add(hNQF_Scaled);

	// Chuck the error bars
	hTotal->Sumw2(0);
	hReal->Sumw2(0);
	hCoherent_Scaled->Sumw2(0);
	hPQF_Scaled->Sumw2(0);
	hNQF_Scaled->Sumw2(0);
	
	// Create a canvas and draw the real and simulated data together
	TCanvas *c1 = new TCanvas("c1","",200,10,1000,750);
	hReal->SetTitle("#pi^{0} Missing Energy in CM Frame");
	hReal->SetLineColor(kAzure);
        hReal->SetLineWidth(2);
        hReal->Draw();
	hTotal->SetLineColor(kSpring);
	hTotal->SetLineWidth(2);
	hTotal->Draw("SAME");
	hCoherent_Scaled->SetLineColor(kOrange-3);
	hCoherent_Scaled->SetLineWidth(2);
	hCoherent_Scaled->Draw("SAME");
	hPQF_Scaled->SetLineColor(kPink-9);
	hPQF_Scaled->SetLineWidth(2);
	hPQF_Scaled->Draw("SAME");
	hNQF_Scaled->SetLineColor(kCyan+1);
	hNQF_Scaled->SetLineWidth(2);
	hNQF_Scaled->Draw("SAME");

	// Create a legend for the histogram
	TLegend *legend = new TLegend(0.1,0.65,0.42,0.9);
	legend->AddEntry(hReal, "Real Data","l");
	legend->AddEntry(hTotal, "MC Data","l");
	legend->AddEntry(hCoherent_Scaled, "MC Coherent #pi^{0} Production", "l");
	legend->AddEntry(hPQF_Scaled, "MC QF #pi^{0} Production on p", "l");
	legend->AddEntry(hNQF_Scaled, "MC QF #pi^{0} Production on n", "l");
	legend->Draw();


}

