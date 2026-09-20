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

const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/FPPA_F2B_WITH_JUNE8_9_DONE/rootfiles/gep5_*";
//const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/CFFOFF_DEFAULT_CUTS/rootfiles/gep5_fullreplay_*";

// ============================================================
// Global cut
// ============================================================

//TCut globalcut = "sbs.tr.n>0&&(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.gemFT.track.chi2ndf[0]<200&&sbs.gemFPP.track.chi2ndf[0]<200&&sbs.hcal.nblk>1";
TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&sbs.gemFPP.track.sclose[0]<0.003&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

#endif
