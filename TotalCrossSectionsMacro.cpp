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
TH1D* Tagg_Chan_To_Energy1D(TString sFile, TH1 *hChan, TString sAxis="X")
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

  TH1D *hEner = (TH1D*)hChan->Clone("hEner");
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
double GetScalingFactor(int ke, TFile *ftFile, TFile *etFile)
{
        // Get tagger count histograms
        TH1D *hFull = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1D *hEmpty = (TH1D*)etFile->Get("He4Pi0/h_ScalarCounts");

	// Divide the full target data by empty target data to get a ratio
        TH1D *hRatio = (TH1D*) hFull->Clone("hRatio");
	hRatio->Divide(hEmpty);

        TH1D *hRatio_conv = Tagg_Chan_To_Energy1D("Tagger_Conversions_by_k_2019_06.txt",hRatio,"X");
	
        double factor = hRatio_conv->GetBinContent(ke);

	return factor;
}


// ------------------------ COUNTS (YIELD) AT CHANNEL/ENERGY ----------------------------- //


// A function that subtracts the empty target data from the full target data
// Input: tagger channel you want to look at, output: counts in elastic peak
// Currently: these histograms are being plotted and I don't know why or how
int GetCounts(int ke, TH3* hFull_ME_3D_conv, TH3* hEmpty_ME_3D_conv, double scale)
{

        // Project 3D histograms at input Tagger channel
        TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke),0,180,ke-1,ke);
        TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke),0,180,ke-1,ke);

	// Subtract empty from full with scaling
        TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone(Form("hSubtracted_%dMeV",ke));
        hSubtracted->Add(hEmpty_ME_projx, -scale); 

	// Create a scaled version of the empty target histogram
        //TH1D *hEmpty_ME_scaled = (TH1D*) hEmpty_ME_projx->Clone(Form("hEmpty_ME_scaled_%dMeV",ke));
        //hEmpty_ME_scaled->Scale(scale);
	
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
        if (param[4] > param[1])
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

	// For plotting purposes only:
//        TF1 *fOther = new TF1("fOther", "gaus",-200,200);
//        fOther->SetParameters(param[0],param[1],param[2]);
//        fPeak->SetLineColor(4);
//	fPeak->Draw("SAME");
//	fOther->SetLineColor(6);
//	fOther->Draw("SAME");

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
        TFile *ftFile = TFile::Open("Full_Target_Result_2.root");
        if(!ftFile->IsOpen()){cout << "Error: full target file could not be opened."; exit(-1);}

        // Empty target file
        TFile *etFile = TFile::Open("Empty_Target_Result_2.root");
        if(!etFile->IsOpen()){cout << "Error: empty target file could not be opened.";exit(-1);}

        // Tagging Efficiency File
        TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
        if(!tagFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened.";exit(-1);}

        // -------------------- Get 3D Histograms and Convert -------------------- //

        // Get 3D, 1 particle event ME Histograms
        TH3D *hFull_ME_3D = (TH3D*)ftFile->Get("He4Pi0/h3D_MEpi0");
        TH3D *hEmpty_ME_3D = (TH3D*)etFile->Get("He4Pi0/h3D_MEpi0");

        // Convert the 3D histograms to be a function of photon energy instead of Tagger channel
        TH3 *hFull_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hFull_ME_3D,"Z");
        TH3 *hEmpty_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hEmpty_ME_3D,"Z");

        // ------------------- Calculate Target Thickness ------------------- //

        // Calculate target thickness
        const double length = 0.05;                     // 5cm??? not sure
        const double density = 124640;                  // g/m3??? found on the ELOG, may be approx value
        const double NA = 6.022*TMath::Power(10,23);    // Avogadro's number
        const double mm = 4.002602;                     // should be molar mass of helium-4
        const double tt = length*density*NA/mm;


        // ------------- Get the Yield of Elastic Pi0 Production ------------- //


        // Create a histogram to fill
        TH1F *hElastic = new TH1F("hElastic", "Counts of Elastic Pi0 Production",200,150,350);
        hElastic->GetXaxis()->SetTitle("Incident Photon Energy (Pi0)");   // want to change to photon energy
        hElastic->GetYaxis()->SetTitle("Counts");

        // Go through the tagger channels and fill the histogram
        for (int ke=150; ke<=350; ke++)
        {
                // Get scaling factor
                double scale = GetScalingFactor(ke, ftFile, etFile);

                // Get yield
                int counts = GetCounts(ke,hFull_ME_3D_conv,hEmpty_ME_3D_conv, scale);  // get the pi0 yield for the channel
                for (int i=0; i<=counts; i++)
                {
                        hElastic->Fill(ke);   // fill the yield histogram for every count pi0 we counted at that channel
                }
        }


        // ----------------- Get Tagging Efficiency Histogram ----------------- //

        TH1D* hTaggEff = (TH1D*)tagFile->Get("makeTaggEff/hist000");
        TH1D* hTaggEff_conv = Tagg_Chan_To_Energy1D("Tagger_Conversions_by_k_2019_06.txt",hTaggEff,"X");

        double m;
        TH1D* hTaggEff_final = new TH1D("hTaggEff","Tagging Efficiency",200,150,350);
        for(int ke=150; ke<=350; ke++)
        {
            m = hTaggEff_conv->GetBinContent(hTaggEff_conv->GetXaxis()->FindBin(ke));
            hTaggEff_final->Fill(ke,m);
            m = 0;
        }



        // ------------------------ Get Tagger Scalars ------------------------ //

        // Also we'll need the tagger scalars histogram
        TH1D* hScalars = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1D* hScalars_conv = Tagg_Chan_To_Energy1D("Tagger_Conversions_by_k_2019_06.txt",hScalars,"X");

        long int n;
        TH1D* hScalars_final = new TH1D("hScalars", "Scalar Counts",200,150,350);
        for(int ke=150; ke<=350; ke++)
        {
            n = hScalars_conv->GetBinContent(hScalars_conv->GetXaxis()->FindBin(ke));
            hScalars_final->Fill(ke,n);
            n = 0;
        }




        // ----------------- Calculate the Total Cross Section ----------------- //

        // Calculate the total cross section histogram (=yield/Ne*taggeff*tt*edet)
        TH1D *hCrossSec = (TH1D*) hElastic->Clone("hCrossSec");
        hCrossSec->Divide(hElastic,hScalars_final,1,tt);                            // divide yield by tagg scalars (Ne) and by target thickness
        hCrossSec->Divide(hCrossSec,hTaggEff_final,TMath::Power(10,34),1);       // also divide by tagging efficiency and convert to microbarns


        // ----------------------- Draw the Histogram(s) ----------------------- //

        // Create a canvas
        TCanvas *c1 = new TCanvas("c1","",200,10,1000,1000);

        // Temp:
        c1->Clear();
//        c1->Divide(2,2);
//        c1->cd(1);
//        hTaggEff_conv->Draw();
//        c1->cd(2);
//        hTaggEff_final->Draw("HIST");
//        c1->cd(3);
//        hElastic->Draw();
//        c1->cd(4);
        c1->cd();
        hCrossSec->Draw("HIST");


//        // Need to clear first or the pi0 missing energy will plot underneath
//        c1->Clear();
//        c1->Divide(2,1);
	
//        // Draw the total cross section histogram as function of Tagger channel
//        c1->cd(1);
//        hCrossSec->SetTitle("Pi0 Production Total Cross Section from 150MeV to 350MeV w/out Edet");
//        hCrossSec->GetYaxis()->SetTitle("Total Cross Section (#mub)");
//	hCrossSec->GetXaxis()->SetRangeUser(76,267);   // this range is about 150MeV-350MeV
//	hCrossSec->Draw();

//        // Draw the total cross section histogram as a function of incident photon energy
//        TH1D *hCrossSec2 = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt", hCrossSec);
//        hCrossSec2->GetYaxis()->SetTitle("Total Cross Section w/out Edet (#mub)");
//        hCrossSec2->GetXaxis()->SetTitle("Incident Photon Energy (MeV)");
//        c1->cd(2);
//        hCrossSec2->Draw();


        // --------------------- Get Detection Efficiency --------------------- //


        // Make a histogram for the proper total cross section points
//        TH1D *hCrossSecProper = new TH1D("hCrossSecProper", "Total Cross Sections",450,0,450);
//        hCrossSecProper->GetXaxis()->SetTitle("Incident Photon Energy (MeV)");
//        hCrossSecProper->GetYaxis()->SetTitle("Total Cross Section (#mub)");

//        // Fill the histogram with total cross section divided by Edet for the energies we have
//        for(int k=150; k<=350; k+=25)
//        {
//            double Edet = GetTotDetEff(k);
//            int bin = hCrossSec->GetXaxis()->FindBin(k);
//            double xs = hCrossSec->GetBinContent(bin)/Edet;
//            hCrossSecProper->SetBinContent(k,xs);
//            //double error = hCrossSec2->GetBinError(bin);
//            //hCrossSecProper->SetBinError(k,error);
//        }

        TH1D *hCrossSecProper = (TH1D*) hCrossSec->Clone("hCrossSecProper");
        hCrossSecProper->SetTitle("Total Cross Sections");
        hCrossSecProper->GetXaxis()->SetTitle("Incident Photon Energy (MeV)");
        hCrossSecProper->GetYaxis()->SetTitle("Total Cross Section (#mub)");
        hCrossSecProper->Scale(1/GetTotDetEff(200));

        // Create a canvas
        TCanvas *c2 = new TCanvas("c2","",200,10,750,750);
        c2->cd();
        hCrossSecProper->Draw("HIST");




        cout << Form("For ke = 200: \n"
                     "Counts: %f \n"
                     "Tagger Scalars: %f \n"
                     "Tagg Eff: %f \n"
                     "Target Thickness: %f \n"
                     "Det Eff: %f \n"
                     "Total Cross Sec: %f \n",
                     hElastic->GetBinContent(hElastic->GetXaxis()->FindBin(200)),
                     hScalars_final->GetBinContent(hScalars->GetXaxis()->FindBin(200)),
                     hTaggEff_final->GetBinContent(hTaggEff->GetXaxis()->FindBin(200)),
                     tt,
                     GetTotDetEff(200),
                     hCrossSecProper->Integral(200,201)
                     );



}
