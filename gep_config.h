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

// ============================================================
// Polarimeter reconstruction cut (theta_FPP / DOCA / zclose /
// dxp,dyp / FT-FPP correlations), ported from polarimeter_recon.C.
// Stricter than globalcut above (adds heep.dt_ADC/heep.dpp and an
// FT-track-to-HCAL-position-match cut) and is evaluated at the FPP
// best-track index rather than a fixed index 0.
// ============================================================

TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.0034)<3*1.46&&abs(heep.dt_ADC[0]-0.2495)<3*1.271&&abs(heep.dpp[0]-0.0115)<3*0.01169&&sqrt(pow((heep.dxECAL-0.01011+0.025*earm.ecal.x)/0.01235,2)+pow((heep.dyECAL+0.005365-(0.005507+0.006248*earm.ecal.x+0.01591*pow(earm.ecal.x,2)+0.01066*pow(earm.ecal.x,3)))/0.01627,2))<=3.5&&sqrt(pow((sbs.gemFT.track.y[0]+sbs.gemFT.track.yp[0]*6.7-sbs.hcal.y[0]+0.002928)/0.05065,2)+pow((sbs.gemFT.track.x[0]+sbs.gemFT.track.xp[0]*6.7-sbs.hcal.x[0]-0.1926)/0.05886,2))<50.5";

// Secondary histogram-fill cuts / windows used downstream (theta range,
// zclose target window, DOCA cut) -- same numeric values as
// polarimeter_recon.C.
const double fpp_theta_min    = 0.6;             // deg
const double fpp_theta_max    = 9.0;             // deg
const double fpp_zclose_mean  = 1.5;             // m
const double fpp_zclose_sigma = (0.55 / 2.0) * 10.0; // m
const double fpp_sclose_cut   = 0.005;           // m

#endif
