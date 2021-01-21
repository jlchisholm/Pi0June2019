///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  OldCrossSecMacro.cpp                                                     //
//                                                                           //
//  Jenna Chisholm                                                           //
//  August 12, 2020                                                          //
//                                                                           //
//  Last Updated: August 20, 2020                                            //
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


//NOTE: This is the old version of the macro that uses tagger channel instead of photon energy.
// It works, but we need to use the other macro that uses photon energy to get further in this analysis


// To do:
//   --> Run tagging efficiency and use that to finish total cross section (done, for one file anyways)
//   --> Get detection efficiency and bring that in
//   --> Also work on converting tagger channels to photon energies (think I can do so within Ant, will need to redo files)
//   --> Then try doing a 2D projection and get the differential cross sections
//   --> Finally, try this all with the Compton data, which will be a bit more difficult



#include "physics.h"
#include "math.h"
using namespace std;


// A function that gets the tagger conversion values from the file
// Input: array to store the tagger conversion values in
// Note: values will be stored as: channel, TDC, scaler, ADC, electron energy, photon energy
void GetTaggerConversion(double (&params)[368][6])
{
	// Open the conversion files but exit if it won't open
	ifstream inFile("Tagger_Conversions_2019_06.txt");
	if(!inFile){cerr << "Cannot open conversions file.\n"; return 1;}
	
	// Some variables we'll need
	string line;
	double ke;

	// Read the values and save them to the array
	for(int i=0; getline(inFile,line);i++)
	{
			istringstream iss(line);
			iss >> params[i][0] >> params[i][1] >> params[i][2] >> params[i][3] >> params[i][4]; 	
			ke = 450 - params[i][4];  // photon energy = beam energy - electron energy
			params[i][5] = ke;
	}

}

// A function that calculates the detection efficiency from the detection efficiency files
// Input: incident photon energy, output: total detection efficiency for that photon energy
double GetTotDetEff(int ke)
{
    // Open file
    TFile *edetFile = TFile::Open(Form("DetEff_%d.root",ke));
    if(!edetFile->IsOpen())
    {
            cout << "Error: detection efficiency file could not be opened.";
            exit(-1);
    }
    TH1F* hThetaIn = (TH1F*)edetFile->Get("h4");
    TH1F* hThetaOut = (TH1F*)edetFile->Get("hThetaOut");
    double in = hThetaIn->Integral();
    double out = hThetaOut->Integral();
    double Edet = out/in;

    edetFile->Close();

    return Edet;
}


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
	
	// Create a canvas and draw the histograms
/*	TCanvas *c1 = new TCanvas("c1","",200,10,1000,750);
	c1->Divide(1,2);
	c1->cd(1);
	hFull->SetLineColor(kBlue);
	hEmpty->SetLineColor(kRed);
	hFull->Draw();
	hEmpty->Draw("SAME");
	c1->cd(2);
	hRatio->GetYaxis()->SetTitle("Scaling Factor");
        hRatio->Draw();
*/
        double factor = hRatio->GetBinContent(ch);

	return factor;

}

// A function that subtracts the empty target data from the full target data
// Input: tagger channel you want to look at, output: counts in elastic peak
int GetCounts(int ch, TFile *ftFile, TFile *etFile)
{
	// Open files or quit if there is an error

/*	// Full target file
        TFile *ftFile = TFile::Open("Full_Target_Result_2.root");
	if(!ftFile->IsOpen())
	{
		cout << "Error: full target file could not be opened.";
		exit(-1);
	}

	// Empty target file
        TFile *etFile = TFile::Open("Empty_Target_Result_2.root");
	if(!etFile->IsOpen())
	{
		cout << "Error: empty target file could not be opened.";
		exit(-1);
	}
*/
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

/*	// Create canvas and draw histograms
        TCanvas *c1 = new TCanvas("c1","",200,10,750,750);
	c1->Divide(1,2);
	c1->cd(1);
	hFull_ME_projx->SetLineColor(kBlue);
	hEmpty_ME_projx->SetLineColor(kRed);
	hEmpty_ME_scaled->SetLineColor(kGreen);
	hFull_ME_projx->Draw();
	hEmpty_ME_projx->Draw("SAME");
	hEmpty_ME_scaled->Draw("SAME");

	// Create a legend for the first histogram
        TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
        legend->AddEntry(hFull_ME_projx,"Full Target Data","l");
        legend->AddEntry(hEmpty_ME_projx,"Empty Target Data","l");
	legend->AddEntry(hEmpty_ME_scaled,"Scaled Empty Target Data","l");
        legend->Draw();

	// Draw the empty subtracted histogram
	c1->cd(2);
	hSubtracted->Draw();
*/
	
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
	hSubtracted->Fit("fit");
//	c1->Update();
	fit->GetParameters(&param[0]);
	
	// Create a fit to mimic that which was fit to the elastic peak and integrate it
	TF1 *fPeak = new TF1("fit", "gaus", -200, 200);
	fPeak->SetParameters(param[3],param[4],param[5]);
	int counts = fPeak->Integral(-200,200);

	return counts;
}


// The main function that calculates the elastic counts at several photon energies and plots them,
// and then calculates the total cross section
void PlotCounts()

{
	// Full target file
        TFile *ftFile = TFile::Open("Full_Target_Result_2.root");
        if(!ftFile->IsOpen())
        {
                cout << "Error: full target file could not be opened.";
                exit(-1);
        }

        // Empty target file
        TFile *etFile = TFile::Open("Empty_Target_Result_2.root");
        if(!etFile->IsOpen())
        {
                cout << "Error: empty target file could not be opened.";
                exit(-1);
        }


	// Create a canvas and histogram to fill
	TCanvas *c3 = new TCanvas("c3","",200,10,750,750);
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


	// Need to clear first or the pi0 missing energy will plot underneath
	c3->Clear();
	
	// Draw the pi0 yield histogram
        c3->Divide(1,2);
        c3->cd(1);
        hElastic->GetXaxis()->SetRangeUser(76,267);
        hElastic->Draw();

	

	// now we get the cross section bois!!!
	
	// Calculate target thickness
	const double length = 0.05;                     // 5cm??? not sure
	const double density = 124640;                  // g/m3??? found on the ELOG, may be approx value
	const double NA = 6.022*TMath::Power(10,23);    // Avogadro's number
	const double mm = 4.002602;                     // should be molar mass of helium-4
	const double tt = length*density*NA/mm;

	// Open the tagging efficiency file and get the histogram
	TFile *tagFile = TFile::Open("~/TaggingEfficienciesJune2019/TaggEff-MOELLER-20.root");
	if(!tagFile->IsOpen())
	{
		cout << "Error: tagging efficiency file could not be opened.";
		exit(-1);
	}
        TH1D* hTaggEff = (TH1D*)tagFile->Get("makeTaggEff/hist000");

        // Get the conversion of tagger channels to incident photon energies
        double params[368][6];
        GetTaggerConversion(params);

        // Rebin tagging efficiency to be by photon energy (I think this is correct)
/*        TH1D* hTaggEffConvert = new TH1D("hTaggEffConvert","Tagging Efficiency",900,0,450);
        int lastbin=0;
        for(int ch=0; ch<=367; ch++)
        {
            double energy = params[ch][5];
            int energybin = hTaggEffConvert->GetXaxis()->FindBin(energy);
            double content = hTaggEff->GetBinContent(ch);
            double errorbar = hTaggEff->GetBinError(ch);
            hTaggEffConvert->SetBinContent(energybin, content);
            hTaggEffConvert->SetBinError(energybin,errorbar);
        }


        c3->cd(1);
        hTaggEff->Draw();
        c3->cd(2);
        hTaggEffConvert->Draw();
*/




        // Also we'll need the tagger scalars histogram
        TH1F* hScalars = (TH1F*)ftFile->Get("He4Pi0/h_ScalarCounts");



        // Get the detection efficiency
        double Edet = GetTotDetEff(300); // we'll use 300MeV edet for now





        // Calculate the total cross section histogram (=yield/Ne*taggeff*tt*edet) currently in nuclei/m2
	TH1F *hCrossSec = (TH1F*) hElastic->Clone("hCrossSec");
        hCrossSec->Divide(hElastic,hScalars,1,tt);
        hCrossSec->Divide(hCrossSec,hTaggEff,TMath::Power(10,34),Edet);
	// hCrossSec->Divide(hEdet);




	
	// Draw the total cross section histogram
        c3->cd(2);
	hCrossSec->SetTitle("Pi0 Production Total Cross Section from 150MeV to 350MeV");
        hCrossSec->GetYaxis()->SetTitle("Total Cross Section (#mub)");
	hCrossSec->GetXaxis()->SetRangeUser(76,267);   // this range is about 150MeV-350MeV
	hCrossSec->Draw();




}
