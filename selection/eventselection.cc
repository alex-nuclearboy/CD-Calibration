/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2019-06
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-07
***********************************************/

//Selection criteria for calibration of the electromagnetic calorimeter

#include "Wasa.hh"
#include "CDataManager.hh"
#include "CHistoManager.hh"
#include "CParameterManager.hh"
#include "CLog.hh"
#include "CConst.hh"
#include "EmsEvent.hh"

#include "WHitBank.hh"
#include "WHitScint.hh"
#include "WVertex.hh"
#include "WVertexBank.hh"
#include "WCluster.hh"
#include "WClusterBank.hh"
#include "WClusterChamb.hh"
#include <WClusterFinder.hh>
#include "WTrack.hh"
#include "WTrackBank.hh"
#include "CDTracksSimple.hh"
#include "FDFTHTracks.hh"

#include "TString.h"
#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>

#include "eventselection.hh"

ClassImp(eventselection);

eventselection::eventselection() {}

eventselection::eventselection(const char * name):CAnalysisModule(name) {

    fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default"));
    //FD table
    kFTH1_old = fDetectorTable->GetDet(CConst::kFTH)->Get1stPlane();
    kFRH1_old = fDetectorTable->GetDet(CConst::kFRH)->Get1stPlane();
    //kFVH1_old = fDetectorTable->GetDet(CConst::kFVH)->Get1stPlane();
    kFWC1_old = fDetectorTable->GetDet(CConst::kFWC)->Get1stPlane();
    //printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i,%i \n",kFTH1_old,kFRH1_old,kFVH1_old,kFWC1_old);
    printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i \n",kFTH1_old,kFRH1_old,kFWC1_old);

    //CD table
    //Get PlaneNumbers of first and last Planes of PS and SE
    kPSfirst_old = fDetectorTable->GetDet(CConst::kPSB)->Get1stPlane();     //141
    kPSlast_old = fDetectorTable->GetDet(CConst::kPSF)->Get1stPlane();      //143
    kSEfirst_old = fDetectorTable->GetDet(CConst::kSEB)->Get1stPlane();     //151
    kSElast_old = fDetectorTable->GetDet(CConst::kSEF)->Get1stPlane() + fDetectorTable->GetDet(CConst::kSEF)->GetWasaPlanes()-1;  //174
    gScreen<<kPSfirst_old<<"\t"<<kPSlast_old<<"\t"<<kSEfirst_old<<"\t"<<kSElast_old<<CLog::endl;

    fCDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));    //CD
    if(fCDTrackFinder!=0) fCDTrackBank = fCDTrackFinder->GetTrackBank();

    fFDTrackFinder = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));          //FD
    if(fFDTrackFinder!=0) fFDTrackBank = fFDTrackFinder->GetTrackBank();

    WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
    fMCTrackBank  = MCTrf->GetTrackBank();
    fMCVertexBank = MCTrf->GetVertexBank();

    fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));    //WMC Event header
    fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));                    //DATA Event Header

    //change Edep to Ekin (WasaParameters)
    //fFDEdep2Ekin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","default"));    //"default" is for protons

    SetupSpectra(name);

}

////////////////////////////////////////////////////////////////////////////////////////////////

eventselection::~eventselection() {}

void eventselection::ProcessEvent() {    //01//

    if (fProcessed) return;
    fProcessed = kTRUE;

/////////////////////////////////////////ANALYSIS START/////////////////////////////////////////

    ////PARTICLE MASSES////
    const Double_t m_pi0 = 0.13497;     //neutral pion mass     [GeV]

//////////////////////////////////////RECONSTRUCTED EVENTS//////////////////////////////////////

    ///////TRIGGER///////

    if (!(fHeader->TriggerNumSet(10)))  return; //trigger #10

    //
    Int_t NumNeutTrackCD = fCDTrackBank->GetEntries(11);    //Neutral Tracks in CD have Type 11
    Int_t NumCharTrackCD = fCDTrackBank->GetEntries(12);    //Charged Tracks in CD have Type 12
    Int_t NumCharTrackFD = fFDTrackBank->GetEntries(2);     //Charged Tracks in FD have type 2

    hNeutralTracksCD[0]->Fill(NumNeutTrackCD);
    hChargedTracksCD[0]->Fill(NumCharTrackCD);
    hChargedTracksFD[0]->Fill(NumCharTrackFD);

////////////////////////////////////////////////////////////////////////////////////////////////

    TLorentzVector Gamma1;
    TLorentzVector Gamma2;
    WTrack *TrackGamma1;
    WTrack *TrackGamma2;

    Double_t best_delta = 99999;

    Double_t InvariantMass;

    if(NumNeutTrackCD >= 2) {

        //Loop over tracks in CD to look the best gamma pair
        for (Int_t l = 0; l < fCDTrackBank->GetEntries(); l++) {

            WTrack *track1 = fCDTrackBank->GetTrack(l);

            if (track1->Type() == kCDN) {

                for (Int_t k = l + 1; k < fCDTrackBank->GetEntries(); k++) {

                    WTrack *track2 = fCDTrackBank->GetTrack(k);

                    if (track2->Type() == kCDN) {

                        Double_t Mom1 = track1->Momentum();
                        Double_t Mom2 = track2->Momentum();
                        Double_t Theta1 = track1->Theta();
                        Double_t Theta2 = track2->Theta();
                        Double_t Phi1 = track1->Phi();
                        Double_t Phi2 = track2->Phi();

                        TVector3 vec_1;
                        vec_1.SetMagThetaPhi(Mom1,Theta1,Phi1);
                        TLorentzVector P_1;
                        P_1.SetVectM(vec_1,0.);

                        TVector3 vec_2;
                        vec_2.SetMagThetaPhi(Mom2,Theta2,Phi2);
                        TLorentzVector P_2;
                        P_2.SetVectM(vec_2,0.);

                        Double_t inv_m = (P_1 + P_2).M();
                        Double_t delta = TMath::Abs(inv_m - m_pi0);

                        if (delta < best_delta) {

                            best_delta = delta;

                            Gamma1 = (Mom1 >= Mom2)?P_1:P_2;    //higher energy gamma has number one
                            Gamma2 = (Mom1 >= Mom2)?P_2:P_1;    //lower energy gamma has number two
                            TrackGamma1 = (Mom1 >= Mom2)?track1:track2;
                            TrackGamma2 = (Mom1 >= Mom2)?track2:track1;

                        }
                    }
                }
            }
        }

        //gamma quanta//
        Double_t p_g_lab[2];
        Double_t Theta_g_lab[2];
        Double_t Phi_g_lab[2];

        TVector3 vec_g_lab[2];
        TLorentzVector P_g[2];

        p_g_lab[0] = TrackGamma1->Momentum();
        Theta_g_lab[0] = TrackGamma1->Theta();
        Phi_g_lab[0] = TrackGamma1->Phi();

        vec_g_lab[0].SetMagThetaPhi(p_g_lab[0],Theta_g_lab[0],Phi_g_lab[0]);
        P_g[0].SetVectM(vec_g_lab[0],0.);

        p_g_lab[1] = TrackGamma2->Momentum();
        Theta_g_lab[1] = TrackGamma2->Theta();
        Phi_g_lab[1] = TrackGamma2->Phi();

        vec_g_lab[1].SetMagThetaPhi(p_g_lab[1],Theta_g_lab[1],Phi_g_lab[1]);
        P_g[1].SetVectM(vec_g_lab[1],0.);

        //Invariant and Missing Masses in CM//
        InvariantMass = (P_g[0] + P_g[1]).M();

        hIM_pion[0]->Fill(InvariantMass);

    }

////////////////////////////////////////////////////////////////////////////////////////////////

    if ( (NumCharTrackFD >= 1) && (NumNeutTrackCD == 2) ) {

        hNeutralTracksCD[1]->Fill(NumNeutTrackCD);
        hChargedTracksCD[1]->Fill(NumCharTrackCD);
        hChargedTracksFD[1]->Fill(NumCharTrackFD);

        hIM_pion[1]->Fill(InvariantMass);

    }

    return;

}   //01//

////////////////////////////////////////////////////////////////////////////////////////////////

void eventselection::SetupSpectra(const char * lpath) {   //02//

    TString hpath1 = "Tracks";
    TString hpath2 = "InvariantMass";

    for (Int_t i = 0; i < 2; i++) {

        hNeutralTracksCD[i] = new TH1F(Form("hNeutralTracksCD_lev%d",i),"",11,-0.5,10.5);
        hNeutralTracksCD[i]->GetXaxis()->SetTitle("Neutral Tracks in CD");
        hNeutralTracksCD[i]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hNeutralTracksCD[i],hpath1);

        hChargedTracksCD[i] = new TH1F(Form("hChargedTracksCD_lev%d",i),"",11,-0.5,10.5);
        hChargedTracksCD[i]->GetXaxis()->SetTitle("Charged Tracks in CD");
        hChargedTracksCD[i]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hChargedTracksCD[i],hpath1);

        hChargedTracksFD[i] = new TH1F(Form("hChargedTracksFD_lev%d",i),"",11,-0.5,10.5);
        hChargedTracksFD[i]->GetXaxis()->SetTitle("Charged Tracks in FD");
        hChargedTracksFD[i]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hChargedTracksFD[i],hpath1);

    }

    for (Int_t j = 0; j < 2; j++) {
        hIM_pion[j] = new TH1F(Form("hIM_pion_lev%d",j),"",1000,0.,0.4);
        hIM_pion[j]->GetXaxis()->SetTitle("m_{#gamma_{1}#gamma_{2}} [GeV]");
        hIM_pion[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hIM_pion[j],hpath2);
    }

    return;

}   //02//

void eventselection::Clear(Option_t *option){
    fProcessed = kFALSE;
    return;
}

void eventselection::Print(Option_t *option){
    return;
}

void eventselection::UserCommand(CCommand * command){
    return;
}
