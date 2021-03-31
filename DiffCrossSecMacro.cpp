///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  DiffCrossSecMacro.cpp                                                    //
//                                                                           //
//  Jenna Chisholm                                                           //
//  October 7, 2020                                                          //
//                                                                           //
//  Last Updated: March 19, 2021                                             //
//                                                                           //
//  Calculates differential cross sections for neutral pion production at    //
//  150, 175, 200, 225, 250, 275, 300, 325, and 350MeV.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


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


// A function that gets the detection efficiency from the detection efficiency files
// Input: incident photon energy, output: total detection efficiency for that photon energy
TH1D* GetTotDetEff(int ke)
{
    // Open file
    TFile *edetFile = TFile::Open(Form("DetEff/DetEff_%d_FinalB.root",ke));
    if(!edetFile->IsOpen())
    {
            cout << "Error: detection efficiency file could not be opened.";
            exit(-1);
    }

    // Get detection efficiency histogram
    TH1D* hEdet = (TH1D*)edetFile->Get("hDetEff");

    return hEdet;
}



// ---------------------- SCALING FACTOR FOR EMPTY TARGET SUBTRACTION --------------------------- //


// A function which draws the tagger counts histograms for full and empty target
// data and makes a histogram for the scaling ratio
// Input: photon energy, output: scaling factor for that channel
double GetScalingFactor(int ke, TFile *ftFile, TFile *etFile)
{
        // Get tagger count histograms
        TH1D *hFull = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
        TH1D *hEmpty = (TH1D*)etFile->Get("He4Pi0/h_ScalarCounts");

	// Divide the full target data by empty target data to get a ratio
        TH1D *hRatio = (TH1D*) hFull->Clone("hRatio");
	hRatio->Divide(hEmpty);

        // Convert to a function of photon energy
        TH1D *hRatio_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hRatio);
	
        // Get factor for this energy
        double factor = hRatio_conv->GetBinContent(hRatio_conv->GetXaxis()->FindBin(ke));

	return factor;
}


// ------------------------ COUNTS (YIELD) AT ENERGY AND ANGLE ----------------------------- //


// A function that subtracts the scaled empty target data from the full target data,
// fits the peak, integrates the coherent peak, and returns a yield count
// Input: photon energy you want to look at, angle you want, converted 3D full and empty target histograms,
// scaling factor for this energy, output: counts in elastic peak
// Currently: these histograms are being plotted and I don't know why or how
int GetCounts(int ke, int th, TH3* hFull_ME_3D_conv, TH3* hEmpty_ME_3D_conv, double scale)
{
        // Get bin numbers for desired ranges
        int thbin = hFull_ME_3D_conv->GetYaxis()->FindBin(th);
        int kebin_l = hFull_ME_3D_conv->GetZaxis()->FindBin(ke-5);
        int kebin_h = hFull_ME_3D_conv->GetZaxis()->FindBin(ke+5);

        // Project 3D histograms at input energy range and theta
        /*TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke),th,th+1,ke-1,ke+1);
        TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke),th,th+1,ke-1,ke+1);*/  // maybe these should be th, th+1 ...
        TH1D *hFull_ME_projx = hFull_ME_3D_conv->ProjectionX(Form("hME_projx_%dMeV",ke),thbin,thbin,kebin_l,kebin_h);
        TH1D *hEmpty_ME_projx = hEmpty_ME_3D_conv->ProjectionX(Form("hEmpty_ME_projx_%dMeV",ke),thbin,thbin,kebin_l,kebin_h);

	// Subtract empty from full with scaling
        TH1D *hSubtracted = (TH1D*) hFull_ME_projx->Clone(Form("hSubtracted_%dMeV",ke));
        hSubtracted->Add(hEmpty_ME_projx, -scale); 

        //Create a scaled version of the empty target histogram -- for plotting only
        TH1D *hEmpty_ME_scaled = (TH1D*) hEmpty_ME_projx->Clone(Form("hEmpty_ME_scaled_%dMeV",ke));
	hEmpty_ME_scaled->Scale(scale);
	
	// Create a fit for the subtracted data (gaussian for QF + gaussian for elastic)	
	TF1 *fit = new TF1("fit", "gaus(0)+gaus(3)", -200, 200);

	// Set fit parameters
	fit->SetParNames("BG Constant", "BG Mean", "BG Sigma", "Peak Constant", "Peak Mean", "Peak Sigma");
        fit->SetParameter("BG Mean", -50);
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
        if (param[4] > param[1]) // make sure you're choosing the right Gaussian
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
//	TF1 *fOther = new TF1("fOther", "gaus",-200,200);
//	fOther->SetParameters(param[3],param[4],param[5]);
//	hSubtracted->Fit("fPeak","Q","+");
//	hSubtracted->Fit("fOther","Q","+");

        if(counts < 0 || mean < -15)
        {
                cout << Form("Bad -> energy: %d, theta: %d, counts: %d, mean: %f",ke,th,counts,mean) << endl;
                counts = 0;
        }

	return counts;
}


// ----------------------- CALCULATE DIFFERENTIAL CROSS SECTION ------------------------- //


// The main function that calculates the elastic counts at several photon energies and plots them,
// and then calculates the total cross section
TH1D* GetDXS(int ke=200)
//void GetDXS(int ke=200)
{
    // --------------------------- Open Files --------------------------- //

    // Full target file
    //TFile *ftFile = TFile::Open("Full_Target_Result_2.root");
    TFile *ftFile = TFile::Open("He4Pi0_Full_FinalB.root");
    if(!ftFile->IsOpen()){cout << "Error: full target file could not be opened."; exit(-1);}

    // Empty target file
   // TFile *etFile = TFile::Open("Empty_Target_Result_2.root");
    TFile *etFile = TFile::Open("He4Pi0_Empty_FinalB.root");
    if(!etFile->IsOpen()){cout << "Error: empty target file could not be opened.";exit(-1);}

    // Tagging Efficiency File
    TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
    if(!tagFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened.";exit(-1);}


	// TEMPORARY:
	
	//TH1F *h_MEpi0 = (TH1F*)ftFile->Get("He4Pi0/h_MEpi0");
	//TH1F *h_MMpi0 = (TH1F*)ftFile->Get("He4Pi0/h_MMpi0");


	// Create a canvas
    	//TCanvas *cTEMP1 = new TCanvas("cTEMP1","",200,10,750,750);
	//TCanvas *cTEMP2 = new TCanvas("cTEMP2","",200,10,750,750);

	//cTEMP1->cd();
	//h_MEpi0->SetLineColor(1);
	//h_MEpi0->SetFillColor(9);
	//h_MMpi0->SetTitle("Pi0 Invariant Mass");
	//h_MEpi0->Draw("hist");
	//cTEMP2->cd();
	//h_MMpi0->SetLineColor(1);
	//h_MMpi0->Draw("pfc");


    // ------------------- Calculate Target Thickness ------------------- //

    // Calculate target thickness
    const double length = 0.053;                     // ~5.3cm?
    const double density = 124640;                  // g/m3, found on the ELOG, corresponds to 1018mbar pressure?
    const double NA = 6.022*TMath::Power(10,23);    // Avogadro's number (/mol)
    const double mm = 4.002602;                     // should be molar mass of helium-4 (g/mol)
    const double tt = length*density*NA/mm;

    // -------------------- Get 3D Histograms and Convert -------------------- //

    // Get 3D, 1 particle event ME Histograms
    TH3D *hFull_ME_3D = (TH3D*)ftFile->Get("He4Pi0/h3D_MEpi0");
    TH3D *hEmpty_ME_3D = (TH3D*)etFile->Get("He4Pi0/h3D_MEpi0");

    // Convert the 3D histograms to be a function of photon energy instead of Tagger channel
    TH3 *hFull_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hFull_ME_3D,"Z");
    TH3 *hEmpty_ME_3D_conv = Tagg_Chan_To_Ener("Tagger_Conversions_by_k_2019_06.txt",hEmpty_ME_3D,"Z");


    // ------------- Get the Yield of Elastic Pi0 Production ------------- //

    // Create a histogram to fill
    TH1D *hElastic = new TH1D("hElastic", Form("Counts of Elastic Pi0 Production for KE=%d MeV",ke),18,0,180);
    hElastic->GetXaxis()->SetTitle("#theta");
    hElastic->GetYaxis()->SetTitle("Counts");


    // Get scaling factor
    double scale = GetScalingFactor(ke, ftFile, etFile);

    int yield[18];

    // Go through the theta and fill the histogram
    for (int th=0; th<18; th++) // only 18 theta bins made by Ant
    {

            int counts = GetCounts(ke,th*10,hFull_ME_3D_conv,hEmpty_ME_3D_conv, scale);  // get the pi0 yield for the energy (+/-5 MeV) and theta
            for (int i=0; i<=counts; i++)
            {
                hElastic->Fill(th*10);   // fill the yield histogram for every count pi0 we counted at that channel
            }
            yield[th] = counts;  // Save the yield for error calculations later
    }
    //hElastic->Sumw2();





    // ----------------- Get Tagging Efficiency Histogram ----------------- //

    TH1D* hTaggEff = (TH1D*)tagFile->Get("makeTaggEff/hist000");
    TH1D* hTaggEff_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hTaggEff);
    double taggeff = hTaggEff_conv->GetBinContent(hTaggEff_conv->GetXaxis()->FindBin(ke));



    // ------------------------ Get Tagger Scalars ------------------------ //

    // Also we'll need the tagger scalars histogram
    TH1D* hScalars = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
    TH1D* hScalars_conv = TaggChanToEnergy("Tagger_Conversions_by_k_2019_06.txt",hScalars);
    double scalars = hScalars_conv->GetBinContent(hScalars_conv->GetXaxis()->FindBin(ke));


    // ------------------------- Get Solid Angle ------------------------- // = 2pi int sintheta dtheta
    TH1D *hSolidAng = new TH1D("hSolidAng", "Solid Angle",18,0,180);
    hSolidAng->GetXaxis()->SetTitle("#theta (degrees)");
    hSolidAng->GetYaxis()->SetTitle("Solid Angle");
    TF1 *sinth = new TF1("sinth","sin(x)");

    double ang, content;
    for (int th=0; th<180; th=th+10) // only 18 theta bins made by Ant
    {
        ang = 2*kPI*(sinth->Integral(th*kD2R,(th+10)*kD2R)); // need to convert to radians first
        hSolidAng->Fill(th,ang);
    }
    hSolidAng->Sumw2(0);

    TCanvas *cA = new TCanvas("cA","",200,10,850,750);
    hSolidAng->SetLineWidth(2);
    hSolidAng->Draw();




    // --------------------- Get Detection Efficiency --------------------- //

    TH1D* hDetEff = GetTotDetEff(ke);
    hDetEff->Sumw2(0);


    // ------------------ Get Differential Cross Section ------------------ //


    // Make a histogram for the differential cross section
    TH1D *hDXS = new TH1D("hDXS", Form("Differential Cross Section at %d MeV #pm 5 MeV",ke),18,0,180);
    hDXS->GetXaxis()->SetTitle("#theta (degrees)");
    //hDXS->GetXaxis()->SetTitleSize(3);
    hDXS->GetYaxis()->SetTitle("Differential Cross Section (#mub)");
    //hDXS->GetYaxis()->SetTitleSize(0.4);
    hDXS->Divide(hElastic,hSolidAng,TMath::Power(10,34),scalars*taggeff*tt); // yield*(conversion factor)/(solidangle*scalars*taggeff*thickness) //need deteff
    hDXS->Divide(hDetEff);

    // Calculate error bars for differential cross section as \delta(dxs) = dxs/sqrt(Y)
    double dxs=0;
    double dxserr=0;
    hDXS->Sumw2(0);
    for (int thbin=1; thbin<18;thbin++)
    {
        dxs = hDXS->GetBinContent(thbin);
        dxserr = dxs/TMath::Sqrt(yield[thbin]);
	cout << "yield: " << yield[thbin] << ", error: " << dxserr << endl;
        if (yield[thbin]==0){hDXS->SetBinError(thbin,0);}  // Need to set error to zero if there is no data
        if (yield[thbin] != 0){hDXS->SetBinError(thbin,dxserr);}
    }



        //    cout << Form("For %dMeV: \n"
        //                 "Tagger Scalars: %f \n"
        //                 "Tagg Eff: %f \n"
        //                 "Target Thickness: %f \n \n",
        //                 ke, scalars, taggeff, tt);

        //    for (int th=0; th<180; th=th+10)
        //    {
        //    cout << Form("For theta = %d: \n"
        //                 "Counts: %f \n"
        //                 "Det Eff: %f \n"
        //                 "Solid Angle: %f \n \n",
        //                 th,
        //                 hElastic->GetBinContent(hElastic->GetXaxis()->FindBin(th)),
        //                 hDetEff->GetBinContent(hDetEff->GetXaxis()->FindBin(th)),
        //                 hSolidAng->GetBinContent(hSolidAng->GetXaxis()->FindBin(th)));
        //    }

        //    cout << Form("Total Cross Section should be: %f \n", hDXS->Integral());



        //    // Plot differential cross section, and pieces
        //    TCanvas *c3 = new TCanvas("c3","",200,10,750,1000);
        //    c3->Divide(2,3);
        //    c3->cd(1);
        //    hElastic->Draw();
        //    c3->cd(2);
        //    hSolidAng->Draw();
        //    c3->cd(3);
        //    hDetEff->Draw();
        //    c3->cd(4);
        //    hTaggEff->Draw();
        //    c3->cd(5);
        //    hScalars->Draw();
        //    c3->cd(6);
        //    //hDXS->SetMarkerStyle(20);
        //    //hDXS->SetMarkerSize(0.5);
        //    hDXS->Draw();  //"E1"



    return hDXS;
    //return hElastic;
    //return hDetEff;


}




// ----------------------- MAIN FUNCTION ------------------------- //


// The main function that calculates the differential cross sections at several photon energies and plots them
void PlotDXS()
{
    // Calculate the nine differential cross sections
    TH1D* hDiffCrossSec[9];
    int ke=150;
    for(int i=0; i<=8; i++)
    {
        hDiffCrossSec[i] = GetDXS(ke);
	hDiffCrossSec[i]->SetLineColor(kAzure-6);
	hDiffCrossSec[i]->SetLineWidth(1);
	//hDiffCrossSec[i]->SetMarkerStyle(8);
	//hDiffCrossSec[i]->SetMarkerSize(0.5);
	//hDiffCrossSec[i]->SetMarkerColor(kAzure);
        ke=ke+25;
    }

    // Create a canvas and draw the differential cross sections
    TCanvas *c2 = new TCanvas("c2","",200,10,850,750);
    c2->Divide(3,3);
    c2->cd(1);
    hDiffCrossSec[0]->Draw();  //"pfc"
    c2->cd(2);
    hDiffCrossSec[1]->Draw();
    c2->cd(3);
    hDiffCrossSec[2]->Draw();
    c2->cd(4);
    hDiffCrossSec[3]->Draw();
    c2->cd(5);
    hDiffCrossSec[4]->Draw();
    c2->cd(6);
    hDiffCrossSec[5]->Draw();
    c2->cd(7);
    hDiffCrossSec[6]->Draw();
    c2->cd(8);
    hDiffCrossSec[7]->Draw();
    c2->cd(9);
    hDiffCrossSec[8]->Draw();

//    TH1D *hDiffCrossSec = GetDXS(225);
    TCanvas *c3 = new TCanvas("c3","",200,10,750,750);
    c3->cd();
    hDiffCrossSec[5]->Draw();
}






