/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2020-07
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
***********************************************/

//Macro to fit a function to the invariant mass distribution of two gamma quanta
//and to determine the calibration factor for SEC

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TInterpreter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TArrow.h>
#include <TObjArray.h>
#include <vector>
#include <TMinuit.h>
#include <Riostream.h>

//POLYNOMIAL background function
Double_t Begin_reject_point = 0.08;
Double_t End_reject_point = 0.19;

Int_t pol_deg = 7;

Bool_t reject;

Double_t Pol_bkgnd(Double_t *x, Double_t *par) {
    if (reject && x[0] > Begin_reject_point && x[0] < End_reject_point) {
        TF1::RejectPoint();
        return 0;
    }

    Double_t fpol = 0;
    for(Int_t i=0; i<pol_deg+1; i++) {
        fpol += par[0+i]*pow(x[0],pol_deg-i);
    };
    return fpol;
}

//Novosibirsk signal function
Double_t Novosibirsk_sig(Double_t *x, Double_t *par) {
    Double_t lambda = sinh(par[1]*sqrt(log(4)))*pow(par[2]*sqrt(log(4)),-1);
    Double_t arg = 0;

    if(1+lambda*(x[0]-par[0])>0) arg=-0.5*pow(log(1+lambda*(x[0]-par[0])),2)*pow(par[1],-2)+pow(par[1],2);
    Double_t nov = 0;
    if(1+lambda*(x[0]-par[0])>0) nov = par[3]*exp(arg);
    return nov;
}

//TOTAL function
Double_t Total_Function(Double_t *x, Double_t *par) {

    Double_t lambda = sinh(par[1]*sqrt(log(4)))*pow(par[2]*sqrt(log(4)),-1);
    Double_t arg = 0;

    if(1+lambda*(x[0]-par[0])>0) arg=-0.5*pow(log(1+lambda*(x[0]-par[0])),2)*pow(par[1],-2)+pow(par[1],2);
    Double_t nov = 0;
    if(1+lambda*(x[0]-par[0])>0) nov = par[3]*exp(arg);

    Double_t fpol = 0;
    for(Int_t i=0; i<pol_deg+1; i++) {
        fpol += par[4+i]*pow(x[0],pol_deg-i);
    };

    Double_t out = 0;
    out = nov+fpol;
    return out;
}

////////////////////////////////////////////////////////////////////////////////////

void InvariantMass_fit() {

    TFile *myFile[3];

    myFile[0] = new TFile("input/DATA-SEC_calibration.root","READ");
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");

    TH1D* hInvariantMass[3];
    TH1D* hInvariantMass_fit = new TH1D();
    TH1D* hInvariantMass_signal = new TH1D();

    hInvariantMass[0] = (TH1D*)myFile[0]->Get("Histograms/InvariantMass/hIM_pion_lev2");
    hInvariantMass[1] = (TH1D*)myFile[1]->Get("Histograms/DATA_lev2_cut0/hIM_pion_lev2_cut0");
    hInvariantMass[2] = (TH1D*)myFile[2]->Get("Histograms/DATA_lev2_cut0/hIM_pion_lev2_cut0");

    hInvariantMass[0]->Rebin(3);
    hInvariantMass[1]->Rebin(2);
    hInvariantMass[2]->Rebin(2);

    hInvariantMass_fit = (TH1D*)hInvariantMass[0]->Clone("hInvariantMass_fit");
    hInvariantMass_signal = (TH1D*)hInvariantMass[0]->Clone("hInvariantMass_signal");

////////////////////////////////////////////////////////////////////////////////////

    Double_t Start_Reject_Bin;
    Double_t BinContent;
    Double_t BinCenterX;

    Int_t Begin_peak_bin;
    Int_t End_peak_bin;

    Int_t Max_peak_bin;
    Double_t Max_peak_bin_center;

    Int_t Max_histo_bin = hInvariantMass_fit->GetXaxis()->GetLast();

    //Get maximum for signal
    for(Int_t k=1; k<Max_histo_bin+1; k++) {

        Start_Reject_Bin = hInvariantMass_fit->GetXaxis()->FindBin(Begin_reject_point);

        BinContent = hInvariantMass_fit->GetBinContent(k);
        BinCenterX = hInvariantMass_fit->GetXaxis()->GetBinCenter(k);

        Begin_peak_bin = hInvariantMass_fit->GetXaxis()->FindBin(Begin_reject_point);
        End_peak_bin = hInvariantMass_fit->GetXaxis()->FindBin(End_reject_point);

        //maximum of signal peak
        hInvariantMass_fit->GetXaxis()->SetRange(Begin_peak_bin,End_peak_bin);
        Max_peak_bin = hInvariantMass_fit->GetMaximumBin();
        Max_peak_bin_center = hInvariantMass_fit->GetXaxis()->GetBinCenter(Max_peak_bin);

    }

    hInvariantMass_fit->GetXaxis()->SetRange(1,Max_histo_bin);

/////////////////////////////////////PARAMETERS/////////////////////////////////////

    //number of parameters
    Int_t bkgnd_param = pol_deg+1;
    Int_t signal_param = 4;
    Int_t total_param = bkgnd_param + signal_param;

/////////////////////////////////////BACKGROUND/////////////////////////////////////

    Double_t Begin_Fit_bkgnd = 0.02;
    Double_t End_Fit_bkgnd = 0.25;
    Double_t Begin_Draw_bkgnd = 0.02;
    Double_t End_Draw_bkgnd = 0.25;

    TF1 *bkgndFcn = new TF1("bkgndFcn",Pol_bkgnd,Begin_Draw_bkgnd,End_Draw_bkgnd,bkgnd_param);

    reject = kTRUE;

    hInvariantMass_fit->Fit(bkgndFcn,"","",Begin_Fit_bkgnd,End_Fit_bkgnd);

    reject=kFALSE;

///////////////////////////////////////SIGNAL///////////////////////////////////////

    Double_t Begin_Fit_signal = Max_peak_bin_center - 0.003;
    Double_t End_Fit_signal = Max_peak_bin_center + 0.003;
    Double_t Begin_Draw_signal = Max_peak_bin_center - 0.003;
    Double_t End_Draw_signal = Max_peak_bin_center + 0.003;

    TF1 *signalFcn = new TF1("signalFcn",Novosibirsk_sig,Begin_Draw_signal,End_Draw_signal,signal_param);

    //set some parameters to start
    signalFcn->SetParName(0, "pos");
    signalFcn->SetParName(1, "width");
    signalFcn->SetParName(2, "tail");
    signalFcn->SetParName(3, "height");

    signalFcn->SetParameter(0,0.139);
    signalFcn->SetParameter(1,0.8);
    signalFcn->SetParameter(2,0.02);
    signalFcn->SetParameter(3,80000.);

    hInvariantMass_fit->Fit(signalFcn,"","",Begin_Fit_signal,End_Fit_signal);

///////////////////////////////////TOTAL FUNCTION///////////////////////////////////

    Double_t Begin_Fit = 0.02;
    Double_t End_Fit = 0.25;
    Double_t Begin_Draw = 0.01;
    Double_t End_Draw = 0.25;

    TF1 *totalFcn = new TF1("totalFcn",Total_Function,Begin_Draw,End_Draw,total_param);

    //set parameters from previous fits

    for(Int_t k=0; k<signal_param; k++){
        totalFcn->SetParameter(k,signalFcn->GetParameter(k));
    }

    for(Int_t k=signal_param; k<total_param; k++){
        totalFcn->SetParameter(k,bkgndFcn->GetParameter(k-signal_param));
    }

    //fit
    hInvariantMass_fit->Fit(totalFcn,"","",Begin_Fit,End_Fit);

    //background
    TF1 *bkgndFcn_end = new TF1("bkgndFcn_end",Pol_bkgnd,Begin_Draw_bkgnd,End_Draw_bkgnd,bkgnd_param);

    for(Int_t k=0; k<bkgnd_param; k++) {
        bkgndFcn_end->SetParameter(k,totalFcn->GetParameter(k+signal_param));
    }

//new histogram of 2 gamma quanta invariant mass after background rejection

    hInvariantMass_signal->Add(bkgndFcn_end,-1);

    cout<<"Max_peak_bin_center: "<<Max_peak_bin_center<<endl;
    cout<<"calibration factor: "<<0.13497/Max_peak_bin_center<<endl;

////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");

    TCanvas* MyCanvas01 = new TCanvas; //create Canvas

    Double_t Ymax01 = 1.2*hInvariantMass[0]->GetMaximum();
    Double_t scale010 = (hInvariantMass[0]->GetMaximum())/(hInvariantMass[1]->GetMaximum());
    Double_t scale011 = (hInvariantMass[0]->GetMaximum())/(hInvariantMass[2]->GetMaximum());

    //hInvariantMass[0]->SetTitle("#pi^{0} #rightarrow #gamma #gamma");
    hInvariantMass[0]->SetMarkerStyle(2);
    hInvariantMass[0]->SetMarkerSize(0.7);
    hInvariantMass[0]->GetXaxis()->SetTitle("m_{#gamma_{1}#gamma_{2}}, GeV/c^{2}");
    hInvariantMass[0]->GetXaxis()->SetTitleSize(0.06);
    hInvariantMass[0]->GetXaxis()->SetTitleOffset(1.0);
    hInvariantMass[0]->GetXaxis()->SetLabelSize(0.05);
    hInvariantMass[0]->GetYaxis()->SetTitle("counts");
    hInvariantMass[0]->GetYaxis()->SetTitleSize(0.06);
    hInvariantMass[0]->GetYaxis()->SetTitleOffset(1.);
    hInvariantMass[0]->GetYaxis()->SetLabelSize(0.05);
    hInvariantMass[0]->GetYaxis()->SetRangeUser(0.,Ymax01);

    hInvariantMass[0]->Draw("E1");

    hInvariantMass[1]->SetMarkerStyle(2);
    hInvariantMass[1]->SetMarkerSize(0.7);
    hInvariantMass[1]->SetLineColor(kGreen);
    hInvariantMass[1]->Scale(scale010);
    //hInvariantMass[1]->DrawCopy("same C");

    hInvariantMass[2]->SetMarkerStyle(2);
    hInvariantMass[2]->SetMarkerSize(0.7);
    hInvariantMass[2]->SetLineColor(kYellow);
    hInvariantMass[2]->Scale(scale011);
    //hInvariantMass[2]->DrawCopy("same C");

    bkgndFcn_end->SetLineColor(kMagenta+2);
    bkgndFcn_end->Draw("same C");

    totalFcn->SetLineColor(kCyan-3);
    totalFcn->Draw("same C");

    hInvariantMass_signal->SetLineColor(kOrange+1);
    hInvariantMass_signal->SetLineWidth(2);
    hInvariantMass_signal->GetXaxis()->SetRangeUser(0.05,0.2);
    hInvariantMass_signal->Draw("same C");

    //pi0 mass
    TLine* line010 = new TLine(0.13497,0,0.13497,Ymax01);

    //peak position from gausian fit
    TLine* line011 = new TLine(Max_peak_bin_center,0,Max_peak_bin_center,Ymax01);

    line010->SetLineColor(2);
    line010->SetLineStyle(1);
    line010->SetLineWidth(1);

    line011->SetLineColor(1);
    line011->SetLineStyle(2);
    line011->SetLineWidth(1);

    line010->Draw("same");
    line011->Draw("same");

    //legend
    TLegend *MyLegend01 = new TLegend(0.485, 0.540, 0.885, 0.885);
    MyLegend01->SetFillStyle(1); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0); MyLegend01->SetTextSize(0.04);
    MyLegend01->AddEntry( hInvariantMass[0], "experimental points", "pe");
    MyLegend01->AddEntry( totalFcn, "fitting function", "l");
    MyLegend01->AddEntry( bkgndFcn_end, "7th degree polynomial", "l");
    MyLegend01->AddEntry( hInvariantMass_signal, "signal", "l");
    MyLegend01->AddEntry( line011, Form("peak m_{#gamma_{1}#gamma_{2}} = %g GeV/c^{2}",Max_peak_bin_center), "l");
    MyLegend01->AddEntry( line010, "m_{#pi^{0}} = 0.1349 GeV/c^{2}", "l");
    //MyLegend01->AddEntry( hInvariantMass[1], "WMC: pd #rightarrow (^{3}He#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    //MyLegend01->AddEntry( hInvariantMass[2], "WMC: pd #rightarrow dp#pi^{0}", "l");

    MyLegend01->Draw("same");

    MyCanvas01->Print("InvariantMass_fit.png","png");
    MyCanvas01->Print("InvariantMass_fit.eps","eps");

    TCanvas* MyCanvas02 = new TCanvas; //create Canvas

    hInvariantMass[0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hInvariantMass[0]->Draw("E1");

    //hInvariantMass[1]->DrawCopy("same C");
    //hInvariantMass[2]->DrawCopy("same C");

    bkgndFcn_end->Draw("same C");
    totalFcn->Draw("same C");

    hInvariantMass_signal->Draw("same C");

    line010->Draw("same");
    line011->Draw("same");

    //legend
    TLegend *MyLegend02 = new TLegend(0.485, 0.540, 0.885, 0.885);
    MyLegend02->SetFillStyle(1001); MyLegend02->SetFillColor(19); MyLegend02->SetLineColor(1); MyLegend02->SetTextSize(0.04); MyLegend02->SetBorderSize(5);
    MyLegend02->AddEntry( hInvariantMass[0], "dane eksperymentalne", "pe");
    MyLegend02->AddEntry( totalFcn, "funkcja dopasowania", "l");
    MyLegend02->AddEntry( bkgndFcn_end, "\\hbox{wielomian tła}", "l");
    MyLegend02->AddEntry( hInvariantMass_signal, "\\hbox{sygnał}", "l");
    MyLegend02->AddEntry( line011, Form("pik m_{#gamma_{1}#gamma_{2}} = %g GeV/c^{2}",Max_peak_bin_center), "l");
    MyLegend02->AddEntry( line010, "m_{#pi^{0}} = 0.1349 GeV/c^{2}", "l");
    //MyLegend02->AddEntry( hInvariantMass[1], "WMC: pd #rightarrow (^{3}He#eta)_{zwiazany} #rightarrow dp#pi^{0}", "l");
    //MyLegend02->AddEntry( hInvariantMass[2], "WMC: pd #rightarrow dp#pi^{0}", "l");

    MyLegend02->Draw();

    MyCanvas02->Print("InvariantMass_fit_pl.png","png");
    MyCanvas02->Print("InvariantMass_fit_pl.eps","eps");

}
