///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  PrettyHistogramsMacro.cpp                                                //
//                                                                           //
//  Jenna Chisholm                                                           //
//  March 26, 2021                                                           //
//                                                                           //
//  Last Updated: March 26, 2021                                             //
//                                                                           //
//  Takes histograms from analysis and makes them look all nice.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "physics.h"
#include "math.h"
using namespace std;



void PlotPretty(int bg1, int p1, int p2, int bg2)
{


    // Create canvas
    TCanvas *c1 = new TCanvas("c1","",200,10,1000,750);

    // --------------- PROMPT-RANDOM REGIONS --------------- //

//    // Open prompt-random windows file
//    TFile *prFile = TFile::Open("Windows_Final.root");
//    if(!prFile->IsOpen()){cout << "Error: windows file could not be opened."; exit(-1);}

//    // Get the prompt-random windows histogram
//    TH1D* hPR = (TH1D*)prFile->Get("GetPromptRandomWindows/h_TaggerTimeCBSubtraction");
//    TH1D* hPR_shade = (TH1D*)prFile->Get("GetPromptRandomWindows/h_TaggerTimeCBSubtraction");


    //c1->Divide(2,1);
    //c1->cd(1);
    //hPR->Draw();
    //c1->cd(2);
//    hPR_shade->Draw();

//    // Need all this to draw the sections, but I am getting some weird lines :(
//    Double_t bm = gPad->GetBottomMargin();
//    Double_t lm = gPad->GetLeftMargin();
//    Double_t rm = gPad->GetRightMargin();
//    Double_t tm = gPad->GetTopMargin();
//    Double_t x1 = hPR_shade->GetXaxis()->GetXmin();
//    Double_t y1 = hPR_shade->GetYaxis()->GetXmin();
//    Double_t x2 = hPR_shade->GetXaxis()->GetXmax();
//    Double_t y2 = hPR_shade->GetYaxis()->GetXmax();
//    TPad *null = new TPad("null","null",0,0,1,1);
//    null->SetFillStyle(0);
//    null->SetFrameFillStyle(0);
//    null->SetLineColor(0);
//    null->Draw("al");
//    null->cd();
//    null->Range(x1-(x2-x1)*(lm/(1-rm-lm)),
//                y1-(y2-y1)*(bm/(1-tm-bm)),
//                x2+(x2-x1)*(rm/(1-rm-lm)),
//                y2+(y2-y1)*(tm/(1-tm-bm)));

//    // Shade various sections

//    int style = 3002;

//    TBox *bR1 = new TBox(-100,0,bg1,1);
//    bR1->SetFillColor(kAzure+10);
//    bR1->SetFillStyle(style);
//    bR1->Draw();

//    TBox *b1 = new TBox(bg1,0,p1,1);
//    b1->SetFillColor(kGray+3);
//    b1->SetFillStyle(style);
//    b1->Draw();

//    TBox *b = new TBox(p1,0,p2,1);
//    b->SetFillColor(kOrange);
//    b->SetFillStyle(style);
//    b->Draw();

//    TBox *b2 = new TBox(p2,0,bg2,1);
//    b2->SetFillColor(kGray+3);
//    b2->SetFillStyle(style);
//    b2->Draw();

//    TBox *bR2 = new TBox(bg2,0,100,1);
//    bR2->SetFillColor(kAzure+10);
//    bR2->SetFillStyle(style);
//    bR2->Draw();

//    TMarker *m =new TMarker(0.0,0.0,2);
//    m->SetMarkerSize(2);
//    m->Draw();
//    null->Update();

//    hPR_shade->DrawClone("cont3 same");



    // --------------- STUFF --------------- //

    // Open run file
    TFile *ftFile = TFile::Open("He4Pi0_Full_Final.root");
    if(!ftFile->IsOpen()){cout << "Error: run file could not be opened."; exit(-1);}
    TFile *ftFileB = TFile::Open("He4Pi0_Full_FinalB.root");
    if(!ftFileB->IsOpen()){cout << "Error: run file could not be opened."; exit(-1);}

    // Get the histograms
//    TH1D* hScalars = (TH1D*)ftFile->Get("He4Pi0/h_ScalarCounts");
//    TH1D* hMM2 = (TH1D*)ftFile->Get("He4Pi0/h_MM2");
//    TH1D* hMMpi0 = (TH1D*)ftFile->Get("He4Pi0/h_MMpi0");
//    TH1D* hMEpi0 = (TH1D*)ftFile->Get("He4Pi0/h_MEpi0");

    TH3D* hMEpi0_3D = (TH3D*)ftFile->Get("He4Pi0/h3D_MEpi0");
    TH3D* hMEpi0_3D_B = (TH3D*)ftFileB->Get("He4Pi0/h3D_MEpi0");


    TH1D *hFull_ME_projz = hMEpi0_3D->ProjectionZ("hME_projz");
    TH1D *hFull_ME_projz_B = hMEpi0_3D_B->ProjectionZ("hME_projz_B");

    hFull_ME_projz_B->SetTitle("#pi^{0} Events");
    hFull_ME_projz_B->SetLineColor(kAzure+2);
    hFull_ME_projz_B->SetFillColor(kAzure+10);
    hFull_ME_projz_B->Draw("HIST");
    hFull_ME_projz->SetLineColor(kOrange-3);
    hFull_ME_projz->SetFillColor(kOrange);
    hFull_ME_projz->Draw("HIST SAME");
    hFull_ME_projz_B->GetYaxis()->SetRangeUser(0,14000);
    gPad->RedrawAxis();

    // Create a legend for the histogram
    TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
    legend->AddEntry(hFull_ME_projz_B, "#splitline{Events taken from}{prompt range (-6,30)}","l");
    legend->AddEntry(hFull_ME_projz, "#splitline{Events taken from}{prompt range (-6,8)}","l");
    legend->Draw();




//    hMM2->SetTitle("#pi^{0} Invariant Mass");
//    hMM2->GetXaxis()->SetRangeUser(0,260);
//    hMM2->GetXaxis()->SetTitle("Invariant Mass");
//    hMM2->SetLineWidth(2);
//    hMM2->Sumw2(0);
//    hMM2->Draw();

//    hMMpi0->SetTitle("#pi^{0} Invariant Mass");
//    hMMpi0->GetXaxis()->SetRangeUser(0,260);
//    hMMpi0->GetXaxis()->SetTitle("Invariant Mass");
//    hMMpi0->SetLineWidth(2);
//    hMMpi0->Sumw2(0);
//    hMMpi0->Draw();



    //hMM2->SetTitle("#pi^{0} Invariant Mass");
    //hMM2->GetXaxis()->SetTitle("Invariant Mass");
    //hMM2->Sumw2(0);
    //hMM2->GetXaxis()->SetRangeUser(0,260);
    //hMM2->SetLineColor(kAzure-1);
    //hMM2->SetFillColor(kAzure-2);
    //hMM2->Draw();
    //hMMpi0->Sumw2(0);
    //hMMpi0->SetTitle("#pi^{0} Invariant Mass");
    //hMMpi0->GetXaxis()->SetTitle("Invariant Mass (MeV)");
    //hMMpi0->GetXaxis()->SetRangeUser(0,260);
    //hMMpi0->GetYaxis()->SetRangeUser(0,250000);
    //hMMpi0->SetLineColor(kOrange-3);
    //hMMpi0->SetFillColor(kOrange);
    //hMMpi0->SetLineColor(kAzure-1);
    //hMMpi0->SetFillColor(kAzure-2);
    //hMMpi0->Draw("SAME");
    //gPad->RedrawAxis();

//    // Create a legend for the histogram
    //TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
    //legend->AddEntry(hMM2, "All Two-Particle-Events","l");
    //legend->AddEntry(hMMpi0, "Two-Photon-Events","l");
    //legend->Draw();

    //hMEpi0->SetTitle("#pi^{0} Missing Energy in CM Frame");
    //hMEpi0->Sumw2(0);
    //hMEpi0->GetXaxis()->SetRangeUser(-200,200);
    //hMEpi0->GetYaxis()->SetRangeUser(0,100000);
    //hMEpi0->SetLineColor(kAzure-1);
    //hMEpi0->SetFillColor(kAzure-2);
    //hMEpi0->Draw();
    //gPad->RedrawAxis();


        // --------------- TAGGING EFFICIENCY --------------- //

    // Open run file
//    TFile *teFile = TFile::Open("He4Pi0_Full_FinalB.root");
//    if(!teFile->IsOpen()){cout << "Error: tagging efficiency file could not be opened."; exit(-1);}



        // --------------- DETECTION EFFICIENCY --------------- //

    // Open run file
    //TFile *deFile = TFile::Open("DetEff/DetEff_275_FinalB.root");
    //if(!deFile->IsOpen()){cout << "Error: detection efficiency file could not be opened."; exit(-1);}

    // Get the histograms
    //TH1D* hIn = (TH1D*)deFile->Get("h4");
    //TH1D* hOut = (TH1D*)deFile->Get("hThetaOut");
    //TH1D* hDetEff = (TH1D*)deFile->Get("hDetEff");

    //hIn->SetLineWidth(2);
    //hIn->SetTitle("Scattered #theta at 275 MeV");
    //hIn->GetXaxis()->SetTitle("#theta (degrees)");
    //hIn->GetYaxis()->SetTitle("#");
    //hIn->Draw();

    //hOut->SetLineWidth(2);
    //hOut->SetTitle("Scattered #theta at 275 MeV");
    //hOut->GetXaxis()->SetTitle("#theta (degrees)");
    //hOut->GetYaxis()->SetTitle("#");
    //hOut->Draw();

    //hDetEff->SetLineWidth(2);
    //hDetEff->SetLineColor(kAzure-1);
    //hDetEff->SetFillColor(kAzure-2);
   // hDetEff->Draw();
	    
    //hIn->SetLineColor(kMagenta);
    //hIn->SetFillColor(kMagenta-9);
    //hIn->Draw();
    //hOut->SetLineColor(kCyan+1);
    //hOut->SetFillColor(kCyan);
    //hOut->Draw("SAME");
    //gPad->RedrawAxis();

	// Create a legend for the histogram
    //TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
    //legend->AddEntry(hIn, "Input Simulated Events","l");
    //legend->AddEntry(hOut, "Output Simulated Events","l");
    //legend->Draw();























}
