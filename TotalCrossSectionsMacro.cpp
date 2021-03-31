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


// A method written by Phil which takes in a 1D histogram that is a function of Tagger
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


// Another method written by Phil which takes in a 3D histogram that is a function of Tagger
// channel and converts it to a function of photon energy based on a given
// Tagger conversion file for the beamtime
// Input: conversion file name, histogram to be converted, axis that Tagger channel is on, Output: converted histogram
TH3* Tagg_Chan_To_Ener(TString sFile, TH3 *hChan, TString sAxis="X")
{
  TTree *tTagg = new TTree("tTagg", "tTagg");
  tTagg->ReadFile(sFile, "Channel:Energy");
  Int_t iTagg = tTagg->Draw("Channel:Energy", "", "goff");
  Double_t *dTaggCh = tTagg->GetV1();
  Double_t *dTaggEn = tTagg->GetV2();

  Bool_t bRevrCh = (dTaggCh[1] < dTaggCh[0]);
  Bool_t bRevrEn = (dTaggEn[1] < dTaggEn[0]);
  Bool_t bSumw2 = (hChan->GetSumw2N() > 0);

  Double_t *dTaggBn;
  dTaggBn = new Double_t[iTagg+1];

  for (Int_t i=0; i<=iTagg; i++)
  {
      if (bRevrEn)
      {
          if (i==0) dTaggBn[0] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else dTaggBn[i] = 0.5*(dTaggEn[iTagg-1-i] + dTaggEn[iTagg-i]);
      }
      else
      {
          if (i==0) dTaggBn[0] = (dTaggEn[0] + 0.5*(dTaggEn[0] - dTaggEn[1]));
          else if (i==iTagg) dTaggBn[iTagg] = (dTaggEn[iTagg-1] + 0.5*(dTaggEn[iTagg-1] - dTaggEn[iTagg-2]));
          else dTaggBn[i] = 0.5*(dTaggEn[i-1] + dTaggEn[i]);
      }
  }

  TH3 *hEner = (TH3*)hChan->Clone("hEner");
  hEner->Reset();

  TAxis *axis;
  if (sAxis == "x" || sAxis == "X") axis = hEner->GetXaxis();
  else if (sAxis == "y" || sAxis == "Y") axis = hEner->GetYaxis();
  else if (sAxis == "z" || sAxis == "Z") axis = hEner->GetZaxis();

  if (axis->GetNbins() != iTagg) cout << "WARNING - Different number of bins between config file and histogram!" << endl << "This could lead to unexpected results!" << endl;
  axis->Set(iTagg, dTaggBn);

  Int_t iOld, iNew;
  for (Int_t iX=1; iX<=hChan->GetNbinsX(); iX++)
  {
      for (Int_t iY=1; iY<=hChan->GetNbinsY(); iY++)
      {
          for (Int_t iZ=1; iZ<=hChan->GetNbinsZ(); iZ++)
          {
              iOld = hChan->GetBin(iX, iY, iZ);
              if (bRevrCh == bRevrEn) iNew = hEner->GetBin(iX, iY, iZ);
              else
              {
                  if (sAxis == "x" || sAxis == "X") iNew = hEner->GetBin(iTagg-iX+1, iY, iZ);
                  else if (sAxis == "y" || sAxis == "Y") iNew = hEner->GetBin(iX, iTagg-iY+1, iZ);
                  else if (sAxis == "z" || sAxis == "Z") iNew = hEner->GetBin(iX, iY, iTagg-iZ+1);
              }
              hEner->SetBinContent(iNew, hChan->GetBinContent(iOld));
              if (bSumw2) hEner->SetBinError(iNew, hChan->GetBinError(iOld));
          }
      }
  }

  return hEner;
}


// ------------------ TOTAL DETECTION EFFICIENCY VALUE AT CHANNEL/ENERGY ------------------------- //


// A function that calculates the detection efficiency from the detection efficiency files
// Input: incident photon energy, output: total detection efficiency for that photon energy
double GetTotDetEff(int ke)
{
    // Open file
    TFile *edetFile = TFile::Open(Form("DetEff/DetEff_%d_FinalB.root",ke));
    if(!edetFile->IsOpen())
    {
            cout << "Error: detection efficiency file could not be opened.";
            exit(-1);
    }
    TH1D* hEdet = (TH1D*)edetFile->Get("hDetEff");

    double Edet = hEdet->Integral()/18.0;  // need to divide by the number of bins to normalize it (?)

    //double in = hThetaIn->Integral();
    //double out = hThetaOut->Integral();
    //double Edet = out/in;

    return Edet;
}


// ---------------------- SCALING FACTOR FOR EMPTY TARGET SUBTRACTION --------------------------- //


// A function which draws the tagger counts histograms for full and empty target
// data and makes a histogram for the scaling ratio
// Input: tagger channel, output: scaling factor for that channel
double GetScalingFactor(int ke, TFile *ftFile, TFile *etFile)
{
        // Get tagger count histograms
        TH1D *hFull = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1D *hEmpty = (TH1D*)etFile->Get("He4Pi0/h_ScalarCounts");

	// Divide the full target data by empty target data to get a ratio
        TH1D *hRatio = (TH1D*) hFull->Clone("hRatio");
	hRatio->Divide(hEmpty);

        TCanvas *ctemp = new TCanvas("ctemp","",200,10,750,750);
        hRatio->SetTitle("Empty Target Subtraction Scaling Factor");
        hRatio->GetYaxis()->SetTitle("Scaling Factor");
        hRatio->GetYaxis()->SetRangeUser(1.055,1.085);
        hRatio->SetLineWidth(2);
        hRatio->Draw("HIST");

        TH1D *hRatio_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hRatio);
	
        double factor = hRatio_conv->GetBinContent(hRatio_conv->GetXaxis()->FindBin(ke));


        // For plotting purposes only:
//        TCanvas *ctemp = new TCanvas("ctemp","",200,10,750,750);

//        //hSubtracted->SetTitle("#pi^{0} Missing Energy in CM Frame");
//        //hSubtracted->GetYaxis()->SetTitle("#");
//        ctemp->cd();
//        hFull->SetLineColor(kAzure-1);
//        hFull->SetFillColor(kAzure-2);
//        hFull->Draw("HIST");
//        hEmpty->SetLineColor(kMagenta);
//        hEmpty->SetFillColor(kMagenta-9);
//        hEmpty->Draw("SAME HIST");
//        gPad->RedrawAxis();

//        // Create a legend for the histogram
//        TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
//        legend->AddEntry(hFull, "Full Target Data","l");
//        legend->AddEntry(hEmpty,"Empty Target Data","l");
//        legend->Draw();



	return factor;
}





// TEMPORARY
//void PlotSub(int ke, TH3* hFull_ME_3D_conv, TH3* hEmpty_ME_3D_conv, double scale)
//{
//    int kebin = hFull_ME_3D_conv->GetZaxis()->FindBin(ke);

//    // Project 3D histograms at input Tagger channel
//    TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke));
//    TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke));

//    // Subtract empty from full with scaling
//    TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone(Form("hSubtracted_%dMeV",ke));
//    hSubtracted->Add(hEmpty_ME_projx, -scale);

//    // Create a scaled version of the empty target histogram
//    TH1D *hEmpty_ME_scaled = (TH1D*) hEmpty_ME_projx->Clone(Form("hEmpty_ME_scaled_%dMeV",ke));
//    hEmpty_ME_scaled->Scale(scale);

//    TCanvas *ctemp = new TCanvas("ctemp","",200,10,750,750);
//    hFull_ME_projx->Sumw2(0);
//    hFull_ME_projx->SetLineColor(kAzure-1);
//    hFull_ME_projx->SetFillColor(kAzure-2);
//    hFull_ME_projx->SetTitle("#pi^{0} Missing Energy in CM Frame");
//    hFull_ME_projx->GetYaxis()->SetTitle("#");
//    hFull_ME_projx->GetXaxis()->SetRangeUser(-200,200);
//    hFull_ME_projx->Draw();
//    hSubtracted->Sumw2(0);
//    hSubtracted->SetLineColor(kOrange-3);
//    hSubtracted->SetFillColor(kOrange);
//    hSubtracted->Draw("SAME");
//    hEmpty_ME_scaled->Sumw2(0);
//    hEmpty_ME_scaled->SetLineColor(kCyan);
//    hEmpty_ME_scaled->SetFillColor(kCyan);
//    hEmpty_ME_scaled->Draw("SAME");
//    hEmpty_ME_projx->Sumw2(0);
//    hEmpty_ME_projx->SetLineColor(kMagenta);
//    hEmpty_ME_projx->SetFillColor(kMagenta-9);
//    hEmpty_ME_projx->Draw("SAME");


//    // Create a legend for the histogram
//    TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
//    legend->AddEntry(hFull_ME_projx, "Full Target Data","l");
//    legend->AddEntry(hSubtracted,"Subtracted Data","l");
//    legend->AddEntry(hEmpty_ME_scaled,"Scaled Empty Target Data","l");
//    legend->AddEntry(hEmpty_ME_projx, "Empty Target Data","l");
//    legend->Draw();
//}







// ------------------------ COUNTS (YIELD) AT CHANNEL/ENERGY ----------------------------- //


// A function that subtracts the empty target data from the full target data
// Input: tagger channel you want to look at, output: counts in elastic peak
// Currently: these histograms are being plotted and I don't know why or how
int GetCounts(int ke, TH3* hFull_ME_3D_conv, TH3* hEmpty_ME_3D_conv, double scale)
{
        // Get bin numbers for desired ranges
        int kebin_l = hFull_ME_3D_conv->GetZaxis()->FindBin(ke-5);
        int kebin_h = hFull_ME_3D_conv->GetZaxis()->FindBin(ke+5);
        int NbinsY = hFull_ME_3D_conv->GetYaxis()->GetNbins();

        // Project 3D histograms at input energy
        TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX("hME_projx");
        TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX("hEmpty_ME_projx");
	//TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke),1,NbinsY,kebin_l,kebin_h);
        //TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke),1,NbinsY,kebin_l,kebin_h);

	// Subtract empty from full with scaling
        TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone("hSubtracted");
        hSubtracted->Add(hEmpty_ME_projx, -scale); 
	
	// Create a fit for the subtracted data (gaussian for QF + gaussian for elastic)	
	TF1 *fit = new TF1("fit", "gaus(0)+gaus(3)", -200, 200);

	// Set fit parameters
	fit->SetParNames("BG Constant", "BG Mean", "BG Sigma", "Peak Constant", "Peak Mean", "Peak Sigma");
        fit->SetParameter("BG Mean", -40);
	fit->SetParameter("BG Sigma", 20);
	fit->SetParameter("Peak Mean", 0);
	fit->SetParameter("Peak Sigma", 5);

	// Fit hSubtracted and save the parameters
	double param[6];
        hSubtracted->Fit("fit","Q");
	fit->GetParameters(&param[0]);

	// Create a fit to mimic that which was fit to the elastic peak and integrate it
        TF1 *fPeak = new TF1("fPeak", "gaus", -200, 200);
        double mean, sigma, constant;
        if (param[4] > param[1])   // make sure you're choosing the right Gaussian
        {
            mean = param[4];
            sigma = param[5];
            constant = param[3];
        }
        else
           {
            mean = param[1];
            sigma = param[2];
            constant = param[0];
        }
        fPeak->SetParameters(constant,mean,sigma);
	int counts = fPeak->Integral(-200,200);

//	// For plotting purposes only:
//	TCanvas *ctemp = new TCanvas("ctemp","",200,10,750,750);
//	hSubtracted->SetTitle("#pi^{0} Missing Energy in CM Frame");
//	hSubtracted->GetYaxis()->SetTitle("#");
//	hSubtracted->SetLineColor(kAzure);
//	hSubtracted->SetLineWidth(2);
//	hSubtracted->Draw("HIST");
//	TF1 *fOther = new TF1("fOther", "gaus",-200,200);
//	fOther->SetParameters(param[0],param[1],param[2]);
//	fit->SetLineWidth(2);
//	fit->SetLineColor(kSpring);
//	fit->Draw("SAME");
//	fPeak->SetLineWidth(2);
//        fPeak->SetLineColor(kOrange-3);
//	fPeak->Draw("SAME");
//	fOther->SetLineWidth(2);
//	fOther->SetLineColor(kPink-9);
//	fOther->Draw("SAME");

//        // Create a legend for the histogram
//        TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
//        legend->AddEntry(hSubtracted, "Data","l");
//	legend->AddEntry(fit,"Data Fit","l");
//	legend->AddEntry(fPeak,"Coherent #pi^{0} Production");
//	legend->AddEntry(fOther,"Incoherent #pi^{0} Production");
//	legend->Draw();


        if(counts < 0 || mean < -15)
        {
                cout << Form("Bad -> energy: %d, counts: %d, mean: %f",ke,counts,mean) << endl;
                counts = 0;
        }

	return counts;
}


// ----------------------- MAIN FUNCTION ------------------------- //


//The main function that calculates the elastic counts at several photon energies and plots them,
//and then calculates the total cross section
void PlotCS()
{
        // --------------------------- Open Files --------------------------- //

        // Full target file
        TFile *ftFile = TFile::Open("He4Pi0_Full_FinalB.root");
        if(!ftFile->IsOpen()){cout << "Error: full target file could not be opened."; exit(-1);}

        // Empty target file
        TFile *etFile = TFile::Open("He4Pi0_Empty_FinalB.root");
        if(!etFile->IsOpen()){cout << "Error: empty target file could not be opened.";exit(-1);}

        // Tagging Efficiency File
        TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
        if(!tagFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened.";exit(-1);}


        // ------------------- Calculate Target Thickness ------------------- //

        // Calculate target thickness
        const double length = 0.053;                     // 5.3cm??? not sure
        const double density = 124640;                  // g/m3??? found on the ELOG, may be approx value
        const double NA = 6.022*TMath::Power(10,23);    // Avogadro's number
        const double mm = 4.002602;                     // should be molar mass of helium-4
        const double tt = length*density*NA/mm;

        // -------------------- Get 3D Histograms and Convert -------------------- //

        // Get 3D, 1 particle event ME Histograms
        TH3D *hFull_ME_3D = (TH3D*)ftFile->Get("He4Pi0/h3D_MEpi0");
        TH3D *hEmpty_ME_3D = (TH3D*)etFile->Get("He4Pi0/h3D_MEpi0");

        // Convert the 3D histograms to be a function of photon energy instead of Tagger channel
        TH3 *hFull_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hFull_ME_3D,"Z");
        TH3 *hEmpty_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hEmpty_ME_3D,"Z");

        // -------------- Get Tagging Efficiency and Tagger Scalars  Histograms -------------- //

        // Get tagging efficiency histogram
        TH1D* hTaggEff = (TH1D*)tagFile->Get("makeTaggEff/hist000");
        TH1D* hTaggEff_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hTaggEff);

        // Also we'll need the tagger scalars histogram
        TH1D* hScalars = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1D* hScalars_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hScalars);

        // ------------------------ Calculate Total Cross Sections ------------------------ //

        // Arrays to save values in
        const int n = 9;
        double energy[n];
        double ke_err[n];

        int yield[n];
        double taggeff[n];
        double scalars[n];
        double deteff[n];

        double xs[n];
        double xs_err[n];

        // Go through each energy value
        int i = 0;
        for (int ke=150; ke<=350; ke=ke+25)
        {
            // Save the energy value
            energy[i] = ke;
            ke_err[i] = 5;

            // Get scaling factor
            double scale = GetScalingFactor(ke, ftFile, etFile);

            // Get yield
            yield[i] = GetCounts(ke,hFull_ME_3D_conv,hEmpty_ME_3D_conv, scale);

            // Get tagging efficiency
            taggeff[i] = hTaggEff_conv->GetBinContent(hTaggEff_conv->GetXaxis()->FindBin(ke));

            // Get tagger scalars
            scalars[i] = hScalars_conv->GetBinContent(hScalars_conv->GetXaxis()->FindBin(ke));

            // Get detection efficiency
            deteff[i] = GetTotDetEff(ke);

            // Calculate total cross section
            xs[i] = (yield[i]*TMath::Power(10,34))/(scalars[i]*taggeff[i]*deteff[i]*tt);

            // Calculate error in total cross section
            xs_err[i] =xs[i]/TMath::Sqrt(yield[i]);
            if (yield[i]==0){xs_err[i]=0;}  // zero error if there's no yield

            i++;
        }

        // ------------------------ Plot Total Cross Sections ------------------------ //

        // Make canvas and graph
        TCanvas *c1 = new TCanvas("c1","",200,10,750,750);
        auto *gCrossSec = new TGraphErrors(n,energy,xs,ke_err,xs_err);
        gCrossSec->SetTitle("Total Cross Sections");
        gCrossSec->GetXaxis()->SetTitle("Incidient Photon Energy (MeV)");
        gCrossSec->GetYaxis()->SetTitle("Total Cross Section (#mub)");

        // Draw
        c1->cd();
        gCrossSec->SetLineColor(kAzure-6);
        gCrossSec->Draw("AP Z");


        cout << Form("For ke = 200: \n"
                     "Counts: %d \n"
                     "Tagger Scalars: %f \n"
                     "Tagg Eff: %f \n"
                     "Target Thickness: %f \n"
                     "Det Eff: %f \n"
                     "Total Cross Sec: %f \n",
                     yield[3],
                     scalars[3],
                     taggeff[3],
                     tt,
                     deteff[3],
                     xs[3]
                     );



        // Plot Subtraction
        //double scale = GetScalingFactor(300,ftFile,etFile);
        //PlotSub(300,hFull_ME_3D_conv,hEmpty_ME_3D_conv, scale);

        GetScalingFactor(275, ftFile, etFile);



}
