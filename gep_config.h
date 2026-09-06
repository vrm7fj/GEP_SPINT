#ifndef GEP_CONFIG_H
#define GEP_CONFIG_H

#include "TCut.h"
#include "TMath.h"

// ============================================================
// Constants
// ============================================================

const double Mp = 0.938272;       // Proton mass [GeV]
const double mu_p = 2.793;        // Proton magnetic moment
const double kappa_p = mu_p - 1.0;
const double rad2deg = 180.0 / TMath::Pi();
const int nmod_ft = 14;
const int nmod_fpp = 32; 

const int MAXHIT = 1000;

// ============================================================
// Input ROOT files
// ============================================================

const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP1/LH2/CFFOFF_DEFAULT_CUTS/rootfiles/gep5_*";
//const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/CFFOFF_DEFAULT_CUTS/rootfiles/gep5_fullreplay_*";

// ============================================================
// Global cut
// ============================================================

//TCut globalcut = "sbs.tr.n>0&&(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.gemFT.track.chi2ndf[0]<200&&sbs.gemFPP.track.chi2ndf[0]<200&&sbs.hcal.nblk>1";
TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.gemFT.track.chi2ndf[0]<200&&sbs.gemFPP.track.chi2ndf[0]<200&&sbs.hcal.nblk>1";

#endif
