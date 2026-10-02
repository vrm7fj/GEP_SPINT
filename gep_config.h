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

//const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/FPPA_F2B_WITH_JUNE8_4_SBSNEW/rootfiles/gep5_*";
//const char *rootfile = "/cache/halla/sbs/prod/GEP_REPLAYS/GEP3/LH2/June8_2025/gep5_fullreplay_

const bool use_runlist = false;   // flip to false to go back to wildcard

const char *rootfile_wildcard1 = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/FPPA_F2B_WITH_JUNE8_9_DONE/rootfiles/gep5_fullreplay_*";
//const char *rootfile_wildcard1 = "/volatile/halla/sbs/vidura/parsed-rootfiles/TODAY/*";

const char *rootdir = "/cache/halla/sbs/prod/GEP_REPLAYS/GEP3/LH2/June8_2025/";

const int runlist[] = {3628,3629,3635,3637,3639,3640,3641,3642,3643,3650,3654,3655,
                        3657,3658,3659,3661,3662,3664,3665,3667,3671,3672,3674,3675,
                        3676,3792,3793,3796,3799,3800,3803,3805,3816,3819,3821,3822,3825};
const int nruns = sizeof(runlist)/sizeof(runlist[0]);

// ============================================================
// Global cut
// ============================================================

//TCut globalcut = "sbs.tr.n>0&&(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.gemFT.track.chi2ndf[0]<200&&sbs.gemFPP.track.chi2ndf[0]<200&&sbs.hcal.nblk>1";
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&sbs.gemFPP.track.sclose[0]<0.003&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//Improved Final
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.1613)<3*1.3&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sbs.gemFPP.track.sclose[0]<0.0025&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

 TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)";

//June8
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(heep.dt_ADC[0]-10)<2.5*3.7&&abs(heep.dpp[0]-0.005)<2.5*0.019&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&sbs.gemFPP.track.sclose[0]<0.025";

//Jun8 with new timing
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-10)<3*3.7&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//MFFO
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.hcal.nblk[0]>1&&(abs(sbs.tr.vz[0]+0.23)<0.01577*3||abs(sbs.tr.vz[0]+0.1196)<0.01376*3||abs(sbs.tr.vz[0]+0.06531)<0.01276*3||abs(sbs.tr.vz[0]-0.0432)<0.01131*3)&&sqrt(pow(((sbs.gemFT.track.y[0]+6.7*sbs.gemFT.track.yp[0]-sbs.hcal.y[0]-0.0001234)/0.04903),2)+pow(((sbs.gemFT.track.x[0]+6.7*sbs.gemFT.track.xp[0]-sbs.hcal.x[0]-0.1986)/0.06138),2))<20";

// ============================================================
// Polarimeter reconstruction cut (theta_FPP / DOCA / zclose /
// dxp,dyp / FT-FPP correlations), ported from polarimeter_recon.C.
// Stricter than globalcut above (adds heep.dt_ADC/heep.dpp and an
// FT-track-to-HCAL-position-match cut) and is evaluated at the FPP
// best-track index rather than a fixed index 0.
// ============================================================

//Improved Final
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.1613)<3*1.3&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)";

//June8
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(heep.dt_ADC[0]-10)<2.5*3.7&&abs(heep.dpp[0]-0.005)<2.5*0.019&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7";

//Jun8 with new timing
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-10)<3*3.7&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//MFFO
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.hcal.nblk[0]>1&&(abs(sbs.tr.vz[0]+0.23)<0.01577*3||abs(sbs.tr.vz[0]+0.1196)<0.01376*3||abs(sbs.tr.vz[0]+0.06531)<0.01276*3||abs(sbs.tr.vz[0]-0.0432)<0.01131*3)&&sqrt(pow(((sbs.gemFT.track.y[0]+6.7*sbs.gemFT.track.yp[0]-sbs.hcal.y[0]-0.0001234)/0.04903),2)+pow(((sbs.gemFT.track.x[0]+6.7*sbs.gemFT.track.xp[0]-sbs.hcal.x[0]-0.1986)/0.06138),2))<20";


// Secondary histogram-fill cuts / windows used downstream (theta range,
// zclose target window, DOCA cut) -- same numeric values as
// polarimeter_recon.C.
const double fpp_theta_min    = 1.0;             // deg
const double fpp_theta_max    = 100.0;             // deg
const double fpp_zclose_mean  = 1.5;             // m
const double fpp_zclose_sigma = (0.55 / 2.0) * 10.0; // m
const double fpp_sclose_cut   = 0.0025;           // m

#endif
