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

// ------------------------------------------------------------
// Module -> layer mapping, from SBS-replay DB
// (DB/20250401/db_sbs.gemFT.dat and db_sbs.gemFPP.dat,
//  "sbs.gemXX.mN.layer = L"; same in 20240409 and 20250101).
//   FT : m0..m5 -> L0..L5 (one module each), m6-m9 -> L6, m10-m13 -> L7
//   FPP: 4 modules per layer, m(4L)..m(4L+3) -> L
// ------------------------------------------------------------
// Displayed y-range (mm) of the 2D residual maps (module-wise,
// layer-wise and comparison pages). Display only: histograms are
// still booked over -2..2 mm and fits use the full range.
const double resid_plot_min = -1.25;
const double resid_plot_max =  1.25;

const int nlayer_ft  = 8;
const int nlayer_fpp = 8;

const int layer_of_mod_ft[nmod_ft] = {
  0, 1, 2, 3, 4, 5,
  6, 6, 6, 6,
  7, 7, 7, 7
};

const int layer_of_mod_fpp[nmod_fpp] = {
  0, 0, 0, 0,   1, 1, 1, 1,   2, 2, 2, 2,   3, 3, 3, 3,
  4, 4, 4, 4,   5, 5, 5, 5,   6, 6, 6, 6,   7, 7, 7, 7
};

// Returns -1 for a module index outside the mapping.
inline int LayerOfModuleFT(int mod)  { return (mod >= 0 && mod < nmod_ft)  ? layer_of_mod_ft[mod]  : -1; }
inline int LayerOfModuleFPP(int mod) { return (mod >= 0 && mod < nmod_fpp) ? layer_of_mod_fpp[mod] : -1; }

const int MAXHIT = 1000;

// ============================================================
// Input ROOT files
// ============================================================

//const char *rootfile = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/FPPA_F2B_WITH_JUNE8_4_SBSNEW/rootfiles/gep5_*";
//const char *rootfile = "/cache/halla/sbs/prod/GEP_REPLAYS/GEP3/LH2/June8_2025/gep5_fullreplay_

const bool use_runlist = false;   // flip to false to go back to wildcard

//const char *rootfile_wildcard1 = "/volatile/halla/sbs/adr/gep_replayed/GEP3/cfoil_zerofield_May23_iter3/rootfiles/gep5_fullreplay_*";
const char *rootfile_wildcard1 = "/volatile/halla/sbs/vidura/parsed-rootfiles/DNP_J8/DNP_J8.root";
//const char *rootfile_wildcard1 = "/volatile/halla/sbs/adr/gep_replayed/GEP3/cfoil_zerofield_May23_iter3/rootfiles/gep5_fullreplay_*";

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

//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(sbs.tr.vz[0]+0.10)<3*0.02212&&sbs.gemFPP.track.sclose[0]<0.005";

//June8
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(heep.dt_ADC[0]-10)<2.5*3.7&&abs(heep.dpp[0]-0.005)<2.5*0.019&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&sbs.gemFPP.track.sclose[0]<0.05";

//June8 Better - **
TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&(sbs.gemFT.track.chi2ndf[0]<15&&sbs.gemFPP.track.chi2ndf[0]<15)&&abs(heep.dpp[0]-0.01586)<2.5*0.01539&&abs(heep.dt_ADC[0]-10)<2.5*3.66&&sbs.gemFPP.track.sclose[0]<0.005";

//Jun8 with new timing
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-10)<3*3.7&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//MFFO
//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.hcal.nblk[0]>1&&(abs(sbs.tr.vz[0]+0.23)<0.01577*3||abs(sbs.tr.vz[0]+0.1196)<0.01376*3||abs(sbs.tr.vz[0]+0.06531)<0.01276*3||abs(sbs.tr.vz[0]-0.0432)<0.01131*3)&&sqrt(pow(((sbs.gemFT.track.y[0]+6.7*sbs.gemFT.track.yp[0]-sbs.hcal.y[0]-0.0001234)/0.04903),2)+pow(((sbs.gemFT.track.x[0]+6.7*sbs.gemFT.track.xp[0]-sbs.hcal.x[0]-0.1986)/0.06138),2))<20";

//TCut globalcut = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.hcal.nblk[0]>1&&(abs(sbs.tr.vz[0]+0.23)<0.01577*3||abs(sbs.tr.vz[0]+0.1196)<0.01376*3||abs(sbs.tr.vz[0]+0.06531)<0.01276*3||abs(sbs.tr.vz[0]-0.0432)<0.01131*3)&&sbs.gemFPP.track.sclose[0]<0.005";

// ============================================================
// Polarimeter reconstruction cut (theta_FPP / DOCA / zclose /
// dxp,dyp / FT-FPP correlations), ported from polarimeter_recon.C.
// Stricter than globalcut above (adds heep.dt_ADC/heep.dpp and an
// FT-track-to-HCAL-position-match cut) and is evaluated at the FPP
// best-track index rather than a fixed index 0.
// ============================================================

//Improved Final
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.1613)<3*1.3&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(sbs.tr.vz[0]+0.10)<3*0.02212&&sbs.gemFPP.track.sclose[0]<0.005";

//June8
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(heep.dt_ADC[0]-10)<2.5*3.7&&abs(heep.dpp[0]-0.005)<2.5*0.019&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7";

//June8 Better - **
TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.013,2)+pow((heep.dyECAL-(0.0008+0.0007474*earm.ecal.x+0.01815*pow(earm.ecal.x,2)+0.005745*pow(earm.ecal.x,3)))/0.01506,2))<=3.5&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&(sbs.gemFT.track.chi2ndf[0]<15&&sbs.gemFPP.track.chi2ndf[0]<15)&&abs(heep.dpp[0]-0.01586)<2.5*0.01539&&abs(heep.dt_ADC[0]-10)<2.5*3.66";

//Jun8 with new timing
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-10)<3*3.7&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

//MFFO
//TCut globalcut_thetafpp = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&sbs.hcal.nblk[0]>1&&(abs(sbs.tr.vz[0]+0.23)<0.01577*3||abs(sbs.tr.vz[0]+0.1196)<0.01376*3||abs(sbs.tr.vz[0]+0.06531)<0.01276*3||abs(sbs.tr.vz[0]-0.0432)<0.01131*3)&&sbs.gemFPP.track.sclose[0]<0.005";


// ============================================================
// Second input set (for comparison). Every page is produced for
// set 1 and again for set 2, followed by side-by-side U/V residual
// comparison pages. Each files string may hold several files or
// wildcards separated by spaces or commas. Leave rootfile_set2
// empty to run on set 1 only. Files and labels can also be passed
// directly:
//   root -l -b -q 'gep_physics.C("set1/*.root", "set2/*.root", "label1", "label2")'
// Set 1 defaults to rootfile_wildcard1 / runlist above and uses
// globalcut / globalcut_thetafpp. Set 2 uses its own cuts below
// (default: same as set 1 -- replace with a different TCut string
// to compare selections).
// ============================================================
const char *rootfile_set2 = "/volatile/halla/sbs/vidura/parsed-rootfiles/DNP_DONE/gep5_fullreplay_*";
//const char *rootfile_set2 = "/volatile/halla/sbs/adr/gep_replayed/GEP3/mult_foil_optics_2/rootfiles/gep5_fullreplay_*";
const char *label_set1    = "Before Alignment";
const char *label_set2    = "After Alignment";

// Kinematic tag appended to the plot titles on the polarimeter
// comparison pages (ROOT TLatex; set to "" to leave titles as they are)
const char *q2_label      = "Q^{2} = 11 GeV^{2}";

//TCut globalcut_set2          =  "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(sbs.tr.vz[0]+0.10)<3*0.02212&&sbs.gemFPP.track.sclose[0]<0.005";
//TCut globalcut_thetafpp_set2 =  "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&abs(sbs.tr.vz[0]+0.10)<3*0.02212&&sbs.gemFPP.track.sclose[0]<0.005";


TCut globalcut_set2          = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.1613)<3*1.3&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sbs.gemFPP.track.sclose[0]<0.0025&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";
TCut globalcut_thetafpp_set2 = "(sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.ngoodhits[0]>2)&&(sbs.gemFPP.track.nhits[0]>4||sbs.gemFPP.track.ngoodhits[0]>2)&&earm.ecal.nblk>2&&sbs.hcal.nblk>1&&abs(sbs.tr.vz+0.1)<0.15&&(heep.ecalo[0]/heep.eprime_eth[0])>0.7&&abs(heep.dt_ADC[0]-0.1613)<3*1.3&&abs(heep.dpp[0]-0.005589)<3*0.02291&&sbs.gemFPP.track.sclose[0]<0.0025&&sqrt(pow((heep.dxECAL-0.01063+0.025*earm.ecal.x)/0.0125,2)+pow((heep.dyECAL+0.004986-(0.004374+0.004684*earm.ecal.x+0.01549*pow(earm.ecal.x,2)+0.009088*pow(earm.ecal.x,3)))/0.01705,2))<=3.5";

// Secondary histogram-fill cuts / windows used downstream (theta range,
// zclose target window, DOCA cut) -- same numeric values as
// polarimeter_recon.C.
const double fpp_theta_min    = 0.85;             // deg
const double fpp_theta_max    = 9.0;             // deg
const double fpp_zclose_mean  = 1.5;             // m
const double fpp_zclose_sigma = (0.55 / 2.0) * 10.0; // m
const double fpp_sclose_cut   = 0.0025;           // m

#endif
