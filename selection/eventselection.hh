#ifndef _EVENTSELECTION_HH_
#define _EVENTSELECTION_HH_

#include "CAnalysisModule.hh"
#include "REventWmcHeader.hh"
#include "REventHeader.hh"
#include "WTrackBank.hh"
#include "WVertexBank.hh"
#include "WHitBank.hh"
#include "WClusterBank.hh"
#include "WHitScint.hh"
#include "CCardWDET.hh"
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <WClusterFinder.hh>
#include <FDETClusters.hh>
#include <FDFTHTracks.hh>
#include "CDTracksSimple.hh"
#include "SECGeo.hh"
#include "SECCluster.hh"
#include "FDEdep2Ekin.hh"
#include <fstream>
#include "TCutG.h"

class eventselection : public CAnalysisModule {

public:
    eventselection();
    explicit eventselection(const char * name);
    virtual ~eventselection();

    virtual void    ProcessEvent();
    virtual void    Clear(Option_t *option = "");
    virtual void    Print(Option_t *option = "");
    virtual void    UserCommand(CCommand * command);

private:
    REventWmcHeader *fEventHeader;      //WMC Event header (contains event weight)
    REventHeader    *fHeader;           //DATA Event Header

    WTrackBank      *fFDTrackBank;      //FD Trackbank
    WTrackBank      *fCDTrackBank;      //CD Trackbank

    WTrackBank      *fMCTrackBank;      //MC Track bank (contains true MC tracks (from event generator)

    FDFTHTracks     *fFDTrackFinder;    //FD Trackfinder
    CDTracksSimple  *fCDTrackFinder;    //CD Trackfinder

    //MCTrackFinder   *fMCTrackFinder;

    //MC vertex bank. Each MC event should have one vertex,
    //which contains all emitted particles as generated by event generator
    WVertexBank     *fMCVertexBank;

    //Used to extract detector plane numbers, which are used eg in FD tracks
    CCardWDET       *fDetectorTable;

    SECGeo          *fGeometry;
    SECCluster      *fEMinCluster;

    //FDEdep2Ekin   *fFDEdep2Ekin;      //change Edep to Ekin (WasaParameters)

    //track type, particle type and FD planes new definition
    enum TrackTypes{kFDN=1,kFDC=2,kCDN=11,kCDC=12};  //FD neutral, FD charged, CD neutral, CD charged

    enum ParticleTypes{
            kDummy=0,kGamma=1,kPositron=2,kElectron=3,kPi0=7,kPiPlus=8,kPiMinus=9,
            kNeutron=13,kProton=14,kEta=17,kDeuteron=45,kTriton=46,kHe3=49,
            kDDummy=50
    };

    enum ForwardDetectorPlanes{
            kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,
            kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9
    };

    //WASA plane number for first(last) plane of these detectors
    //(take into acconunt at the beginning of the analysis program)
    Int_t   kFTH1_old, kFRH1_old, kFVH1_old, kFWC1_old, kPSfirst_old, kPSlast_old, kSEfirst_old, kSElast_old;

protected:
    void    SetupSpectra(const char * path);

    TH1F    *hNeutralTracksCD[4],*hChargedTracksCD[4],*hChargedTracksFD[4];
    TH1F    *hIM_pion[4];

    ClassDef(eventselection,0);

};

#endif