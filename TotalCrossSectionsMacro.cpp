///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TotalCrossSectionsMacro.cpp                                              //
//                                                                           //
//  Jenna Chisholm                                                           //
//  August 12, 2020                                                          //
//                                                                           //
//  Last Updated: November 17, 2020                                         //
//                                                                           //
//  Reads in full target and empty target results and divides their Tagger   //
//  counts to get the scaling factors for empty target subtraction. Then     //
//  uses the scaling factors to do an empty target subtraction on the pi0    //
//  missing energy histogram for Tagger channels corresponding to energies   //
//  150MeV to 350MeV. It fits each of those subtracted histograms with a     //
//  double gaussian, and then integrates under just the gaussian from the    //
//  elastic pi0 production peak (centered at 0). This integrated number is   //
//  the counts of elastic pi0 events and is plotted in a histogram as a      //
//  function of Tagger channel/photon energy. From this histogram and other  //
//  quantities, the total cross section is calculated.                       //          
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// This file does that total cross section calculations properly, but probably could be cleaned up


#include "physics.h"
#include "math.h"
using namespace std;



// ----------------------- TAGGER HISTOGRAM CONVERSION ----------------------------- //


// A method written by Phil which takes in a histogram that is a function of Tagger
// channel and converts it to a function of photon energy based on a given
// Tagger conversion file for the beamtime
// Input: conversion file name, histogram to be converted, Output: converted histogram
TH1D* TaggChanToEnergy(TString sFile, TH1D *hCrossSec)
{
    TTree *tTagg = new TTree("tTagg","tTagg");            // TTree for conversion?
    tTagg->ReadFile(sFile, "Channel:Energy");             // read conversion file into tree
    int iTagg = tTagg->Draw("Channel:Energy","","goff");  // get the number of channels/energies?
    double *dTaggCh = tTagg->GetV1();                     // first column in file/tree is channel
    double *dTaggEn = tTagg->GetV2();                     // second column in file/tree is energy

    // I think we're checking the order of things
    bool bRevrCh = (dTaggCh[1] < dTaggCh [0]);   // check if channels are in descending order
    bool bRevrEn = (dTaggEn[1] < dTaggEn[0]);    // check if energies are in descending order
    bool bSumw2 = (hCrossSec->GetSumw2N() > 0);  // sum of weights(???) positive?

    double *dTaggBn;                  // for lowest value of tagger bins
    dTaggBn = new double[iTagg+1];    // set it to be an array of length # of bins

    // So we're going through and getting the lowest value/edge in each bin?
    for (int i=0; i<=iTagg; i++)
    {
        if (bRevrEn) // if energy is in descending order, we'll need to go backwards through indices (i.e. iTagg-i)
        {
            if (i==0) dTaggBn[0] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
            else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
            else dTaggBn[i] = 0.5*(dTaggEn[iTagg-1-i] + dTaggEn[iTagg-i]);
        }
        else
        {
          if (i==0) dTaggBn[0] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else dTaggBn[i] = 0.5*(dTaggEn[i-1] + dTaggEn[i]);  // low edge bin energy is the average of this energy and the previous energy
        }
    }

    // Plot
    TH1D *hTaggCS = new TH1D("hTaggCS", "Pi0 Production Total Cross Section from 150MeV to 350MeV", iTagg, dTaggBn);

    for (int i=0; i<iTagg; i++)
    {
        if (bRevrCh == bRevrEn) // if both channels and energies are both descending order or both are in ascending order
        {
            hTaggCS->SetBinContent(i+1, hCrossSec->GetBinContent(i+1));
            if (bSumw2) hTaggCS->SetBinError(i+1, hCrossSec->GetBinError(i+1));
        }
        else // basically if channels are ascending but energies are descending or vice versa
        {
            hTaggCS->SetBinContent(i+1, hCrossSec->GetBinContent(iTagg-i));
            if (bSumw2) hTaggCS->SetBinError(i+1, hCrossSec->GetBinError(iTagg-i));
        }
    }

    return hTaggCS;
}


// ------------------ TOTAL DETECTION EFFICIENCY VALUE AT CHANNEL/ENERGY ------------------------- //


// A function that calculates the detection efficiency from the detection efficiency files
// Input: incident photon energy, output: total detection efficiency for that photon energy
double GetTotDetEff(int ke)
{
    // Open file
    TFile *edetFile = TFile::Open(Form("DetEff/DetEff_%d.root",ke));
    if(!edetFile->IsOpen())
    {
            cout << "Error: detection efficiency file could not be opened.";
            exit(-1);
    }
    TH1D* hEdet = (TH1D*)edetFile->Get("hDetEff");

    double Edet = hEdet->Integral()/180.0;

    //double in = hThetaIn->Integral();
    //double out = hThetaOut->Integral();
    //double Edet = out/in;

    return Edet;
}


// ---------------------- SCALING FACTOR FOR EMPTY TARGET SUBTRACTION --------------------------- //


// A function which draws the tagger counts histograms for full and empty target
// data and makes a histogram for the scaling ratio
// Input: tagger channel, output: scaling factor for that channel
double GetScalingFactor(int ch, TFile *ftFile, TFile *etFile)
{
        // Get tagger count histograms
        TH1F *hFull = (TH1F*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1F *hEmpty = (TH1F*)etFile->Get("He4Pi0/h_ScalarCounts");

	// Divide the full target data by empty target data to get a ratio
	TH1F *hRatio = (TH1F*) hFull->Clone("hRatio");
	hRatio->Divide(hEmpty);
	
        double factor = hRatio->GetBinContent(ch);

	return factor;
}


// ------------------------ COUNTS (YIELD) AT CHANNEL/ENERGY ----------------------------- //


// A function that subtracts the empty target data from the full target data
// Input: tagger channel you want to look at, output: counts in elastic peak
// Currently: these histograms are being plotted and I don't know why or how
int GetCounts(int ch, TFile *ftFile, TFile *etFile)
{
        // LETS CONVERT BEFOREHAND LIKE IN THE DXS FILE

	// Get 3D, 1 particle event ME Histograms
        TH3D *hFull_ME_3D = (TH3D*)ftFile->Get("He4Pi0/h3D_MEpi0");
        TH3D *hEmpty_ME_3D = (TH3D*)etFile->Get("He4Pi0/h3D_MEpi0");
	
        // Project 3D histograms at input Tagger channel
        TH1D *hFull_ME_projx = hFull_ME_3D->ProjectionX(Form("hME_projx_ch%d",ch),0,180,ch-1,ch);
        TH1D *hEmpty_ME_projx = hEmpty_ME_3D->ProjectionX(Form("hEmpty_ME_projx_ch%d",ch),0,180,ch-1,ch);

	// Get scaling factor
        double scale = GetScalingFactor(ch, ftFile, etFile);

	// Subtract empty from full with scaling
        TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone(Form("hSubtracted_ch%d",ch));
        hSubtracted->Add(hEmpty_ME_projx, -scale); 

	// Create a scaled version of the empty target histogram
        TH1D *hEmpty_ME_scaled = (TH1D*) hEmpty_ME_projx->Clone(Form("hEmpty_ME_scaled_ch%d",ch));
	hEmpty_ME_scaled->Scale(scale);
	
	// Create a fit for the subtracted data (gaussian for QF + gaussian for elastic)	
	TF1 *fit = new TF1("fit", "gaus(0)+gaus(3)", -200, 200);

	// Set fit parameters
	fit->SetParNames("BG Constant", "BG Mean", "BG Sigma", "Peak Constant", "Peak Mean", "Peak Sigma");
	fit->SetParameter("BG Mean", -30);
	fit->SetParameter("BG Sigma", 20);
	fit->SetParameter("Peak Mean", 0);
	fit->SetParameter("Peak Sigma", 5);

	// Fit hSubtracted and save the parameters
	double param[6];
        hSubtracted->Fit("fit","Q");
	fit->GetParameters(&param[0]);
	
	// Create a fit to mimic that which was fit to the elastic peak and integrate it
	TF1 *fPeak = new TF1("fit", "gaus", -200, 200);
	fPeak->SetParameters(param[3],param[4],param[5]);
	int counts = fPeak->Integral(-200,200);

	return counts;
}


// ----------------------- MAIN FUNCTION ------------------------- //


// The main function that calculates the elastic counts at several photon energies and plots them,
// and then calculates the total cross section
void PlotCS()
{
        // --------------------------- Open Files --------------------------- //

        // Full target file
        TFile *ftFile = TFile::Open("Full_Target_Result_2.root");
        if(!ftFile->IsOpen()){cout << "Error: full target file could not be opened."; exit(-1);}

        // Empty target file
        TFile *etFile = TFile::Open("Empty_Target_Result_2.root");
        if(!etFile->IsOpen()){cout << "Error: empty target file could not be opened.";exit(-1);}

        // Tagging Efficiency File
        TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
        if(!tagFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened.";exit(-1);}


        // ------------------- Calculate Target Thickness ------------------- //

        // Calculate target thickness
        const double length = 0.05;                     // 5cm??? not sure
        const double density = 124640;                  // g/m3??? found on the ELOG, may be approx value
        const double NA = 6.022*TMath::Power(10,23);    // Avogadro's number
        const double mm = 4.002602;                     // should be molar mass of helium-4
        const double tt = length*density*NA/mm;


        // ------------- Get the Yield of Elastic Pi0 Production ------------- //

        // Create a histogram to fill
        TH1F *hElastic = new TH1F("hElastic", "Counts of Elastic Pi0 Production",368,0,368);
        hElastic->GetXaxis()->SetTitle("Tagger Channel");   // want to change to photon energy
	hElastic->GetYaxis()->SetTitle("Counts");	

	// Go through the tagger channels and fill the histogram
        for (int ch=76; ch<=267; ch++)
	{
                int counts = GetCounts(ch,ftFile,etFile);  // get the pi0 yield for the channel
		for (int i=0; i<=counts; i++)              
		{
                        hElastic->Fill(ch);   // fill the yield histogram for every count pi0 we counted at that channel
                }
	}


        // ----------------- Get Tagging Efficiency Histogram ----------------- //

        TH1D* hTaggEff = (TH1D*)tagFile->Get("makeTaggEff/hist000");


        // ------------------------ Get Tagger Scalars ------------------------ //

        // Also we'll need the tagger scalars histogram
        TH1F* hScalars = (TH1F*)ftFile->Get("He4Pi0/h_ScalarCounts");




        // ----------------- Calculate the Total Cross Section ----------------- //

        // Calculate the total cross section histogram (=yield/Ne*taggeff*tt*edet)
        TH1D *hCrossSec = (TH1D*) hElastic->Clone("hCrossSec");
        hCrossSec->Divide(hElastic,hScalars,1,tt);                            // divide yield by tagg scalars (Ne) and by target thickness
        hCrossSec->Divide(hCrossSec,hTaggEff,TMath::Power(10,34),1);       // also divide by tagging efficiency and convert to microbarns


        // ----------------------- Draw the Histogram(s) ----------------------- //

        // Create a canvas
        TCanvas *c1 = new TCanvas("c1","",200,10,1000,750);

        // Need to clear first or the pi0 missing energy will plot underneath
        c1->Clear();
        c1->Divide(2,1);
	
        // Draw the total cross section histogram as function of Tagger channel
        c1->cd(1);
        hCrossSec->SetTitle("Pi0 Production Total Cross Section from 150MeV to 350MeV w/out Edet");
        hCrossSec->GetYaxis()->SetTitle("Total Cross Section (#mub)");
	hCrossSec->GetXaxis()->SetRangeUser(76,267);   // this range is about 150MeV-350MeV
	hCrossSec->Draw();

        // Draw the total cross section histogram as a function of incident photon energy
        TH1D *hCrossSec2 = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt", hCrossSec);
        hCrossSec2->GetYaxis()->SetTitle("Total Cross Section w/out Edet (#mub)");
        hCrossSec2->GetXaxis()->SetTitle("Incident Photon Energy (MeV)");
        c1->cd(2);
        hCrossSec2->Draw();


        // --------------------- Get Detection Efficiency --------------------- //


        // Make a histogram for the proper total cross section points
        TH1D *hCrossSecProper = new TH1D("hCrossSecProper", "Total Cross Sections",450,0,450);
        hCrossSecProper->GetXaxis()->SetTitle("Incident Photon Energy (MeV)");
        hCrossSecProper->GetYaxis()->SetTitle("Total Cross Section (#mub)");

        // Fill the histogram with total cross section divided by Edet for the energies we have
        for(int k=150; k<=350; k+=25)
        {
            double Edet = GetTotDetEff(k);
            int bin = hCrossSec2->GetXaxis()->FindBin(k);
            double xs = hCrossSec2->GetBinContent(bin)/Edet;
            hCrossSecProper->SetBinContent(k,xs);
            double error = hCrossSec2->GetBinError(bin);
            hCrossSecProper->SetBinError(k,error);
        }

        // Create a canvas
        TCanvas *c2 = new TCanvas("c2","",200,10,750,750);
        c2->cd();
        hCrossSecProper->Draw("E0");




        cout << Form("For ke = 200: \n"
                     "Counts: %f \n"
                     "Tagger Scalars: %f \n"
                     "Tagg Eff: %f \n"
                     "Target Thickness: %f \n"
                     "Det Eff: %f \n"
                     "Total Cross Sec: %f \n",
                     hElastic->GetBinContent(hElastic->GetXaxis()->FindBin(112)),
                     hScalars->GetBinContent(hScalars->GetXaxis()->FindBin(112)),
                     hTaggEff->GetBinContent(hTaggEff->GetXaxis()->FindBin(112)),
                     tt,
                     GetTotDetEff(200),
                     hCrossSecProper->Integral(195,205)
                     );



}
