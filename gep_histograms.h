#ifndef GEP_HISTOGRAMS_H
#define GEP_HISTOGRAMS_H

#include "TH1D.h"
#include "TH2D.h"
#include "gep_fill_vectors.h"
#include "gep_config.h"

struct GEPHistograms {
  TH1D *h_dx;
  TH1D *h_dy;
  TH1D *h_dxp;
  TH1D *h_dyp;
  TH1D *h_vz;
  TH1D *h_vx;
  TH1D *h_vy;
  TH2D *h_vxvy;
  TH2D *h_eresidu_FT_module;
  TH2D *h_eresidv_FT_module;
  TH2D *h_eresidu_FPP_module;
  TH2D *h_eresidv_FPP_module;

  // Layer-wise residual histograms: x-axis = layer index, y-axis = residual.
  // Modules are combined into layers via layer_of_mod_ft/fpp (gep_config.h).
  TH2D *h_eresidu_FT_layer;
  TH2D *h_eresidv_FT_layer;
  TH2D *h_eresidu_FPP_layer;
  TH2D *h_eresidv_FPP_layer;

  // Overall (all-modules-combined) residual distributions
  TH1D *h_eresidu_FT;
  TH1D *h_eresidv_FT;
  TH1D *h_eresidu_FPP;
  TH1D *h_eresidv_FPP;
};

GEPHistograms CreateHistograms() {
  GEPHistograms hist;
  hist.h_dx  = new TH1D("h_dx",  "Track #DeltaX;x_{FT} - x_{FPP} (m);Counts", 200, -0.1, 0.1);
  hist.h_dy  = new TH1D("h_dy",  "Track #DeltaY;y_{FT} - y_{FPP} (m);Counts", 200, -0.1, 0.1);
  hist.h_dxp = new TH1D("h_dxp", "Track #DeltaX';x'_{FT} - x'_{FPP};Counts", 200, -0.05, 0.05);
  hist.h_dyp = new TH1D("h_dyp", "Track #DeltaY';y'_{FT} - y'_{FPP};Counts", 200, -0.05, 0.05);

  //Target vertex
  hist.h_vz = new TH1D("h_vz", "Target vz; vz (m); Counts", 100, -0.5, 0.5);
  hist.h_vx = new TH1D("h_vx", "Target vx; vx (#mum); Counts", 100, 24.7, 25.7);
  hist.h_vy = new TH1D("h_vy", "Target vy; vy (#mum); Counts", 100, -52.8, -52.0);
  hist.h_vxvy = new TH2D("h_vxvy", "Target vx vs vy; vx (#mum);vy (#mum)", 100, 24.7, 25.7, 100, -52.8, -52.0);

  //Module-wise residual histograms: x-axis = module index, y-axis = residual
  hist.h_eresidu_FT_module  = new TH2D("h_eresidu_FT_module",  "eresidu FT; module; eresidu",  nmod_ft,  -0.5, nmod_ft-0.5,  200, -2, 2);
  hist.h_eresidv_FT_module  = new TH2D("h_eresidv_FT_module",  "eresidv FT; module; eresidv",  nmod_ft,  -0.5, nmod_ft-0.5,  200, -2, 2);
  hist.h_eresidu_FPP_module = new TH2D("h_eresidu_FPP_module", "eresidu FPP; module; eresidu", nmod_fpp, -0.5, nmod_fpp-0.5, 200, -2, 2);
  hist.h_eresidv_FPP_module = new TH2D("h_eresidv_FPP_module", "eresidv FPP; module; eresidv", nmod_fpp, -0.5, nmod_fpp-0.5, 200, -2, 2);

  //Layer-wise residual histograms: x-axis = layer index, y-axis = residual
  hist.h_eresidu_FT_layer  = new TH2D("h_eresidu_FT_layer",  "FT layer-wise eresidu; layer; eresidu (mm)",  nlayer_ft,  -0.5, nlayer_ft-0.5,  200, -2, 2);
  hist.h_eresidv_FT_layer  = new TH2D("h_eresidv_FT_layer",  "FT layer-wise eresidv; layer; eresidv (mm)",  nlayer_ft,  -0.5, nlayer_ft-0.5,  200, -2, 2);
  hist.h_eresidu_FPP_layer = new TH2D("h_eresidu_FPP_layer", "FPP layer-wise eresidu; layer; eresidu (mm)", nlayer_fpp, -0.5, nlayer_fpp-0.5, 200, -2, 2);
  hist.h_eresidv_FPP_layer = new TH2D("h_eresidv_FPP_layer", "FPP layer-wise eresidv; layer; eresidv (mm)", nlayer_fpp, -0.5, nlayer_fpp-0.5, 200, -2, 2);

  //Overall (all-modules-combined) residual distributions
  hist.h_eresidu_FT  = new TH1D("h_eresidu_FT",  "Overall eresidu FT; eresidu (mm); Counts",  120, -2, 2);
  hist.h_eresidv_FT  = new TH1D("h_eresidv_FT",  "Overall eresidv FT; eresidv (mm); Counts",  120, -2, 2);
  hist.h_eresidu_FPP = new TH1D("h_eresidu_FPP", "Overall eresidu FPP; eresidu (mm); Counts", 120, -2, 2);
  hist.h_eresidv_FPP = new TH1D("h_eresidv_FPP", "Overall eresidv FPP; eresidv (mm); Counts", 120, -2, 2);

  return hist;
}

void FillHistograms(const GEPData &data, GEPHistograms &hist) {
  for (size_t i = 0; i < data.ft_x.size(); i++) {
    hist.h_dx->Fill(data.ft_x[i] - data.fpp_x[i]);
    hist.h_dy->Fill(data.ft_y[i] - data.fpp_y[i]);
    hist.h_dxp->Fill(data.ft_xp[i] - data.fpp_xp[i]);
    hist.h_dyp->Fill(data.ft_yp[i] - data.fpp_yp[i]);
  }

  for (size_t i = 0; i < data.t_vz.size(); i++) {
    hist.h_vz->Fill(data.t_vz[i]);
    hist.h_vx->Fill(data.t_vx[i]*1e6);
    hist.h_vy->Fill(data.t_vy[i]*1e6);
    hist.h_vxvy->Fill(data.t_vx[i]*1e6, data.t_vy[i]*1e6);
  }

  for (size_t i = 0; i < data.mod_ft.size(); i++) {
    hist.h_eresidu_FT_module->Fill(data.mod_ft[i], data.eresidu_ft[i]);
    hist.h_eresidv_FT_module->Fill(data.mod_ft[i], data.eresidv_ft[i]);
    hist.h_eresidu_FT->Fill(data.eresidu_ft[i]);
    hist.h_eresidv_FT->Fill(data.eresidv_ft[i]);

    int layer = LayerOfModuleFT(int(data.mod_ft[i]));
    if (layer >= 0) {
      hist.h_eresidu_FT_layer->Fill(layer, data.eresidu_ft[i]);
      hist.h_eresidv_FT_layer->Fill(layer, data.eresidv_ft[i]);
    }
  }

  for (size_t i = 0; i < data.mod_fpp.size(); i++) {
    hist.h_eresidu_FPP_module->Fill(data.mod_fpp[i], data.eresidu_fpp[i]);
    hist.h_eresidv_FPP_module->Fill(data.mod_fpp[i], data.eresidv_fpp[i]);
    hist.h_eresidu_FPP->Fill(data.eresidu_fpp[i]);
    hist.h_eresidv_FPP->Fill(data.eresidv_fpp[i]);

    int layer = LayerOfModuleFPP(int(data.mod_fpp[i]));
    if (layer >= 0) {
      hist.h_eresidu_FPP_layer->Fill(layer, data.eresidu_fpp[i]);
      hist.h_eresidv_FPP_layer->Fill(layer, data.eresidv_fpp[i]);
    }
  }
}

// ============================================================
// Polarimeter reconstruction histograms (theta_FPP, DOCA, zclose,
// dxp/dyp, FT/FPP position-slope correlations, chi2/ndf), ported
// from polarimeter_recon.C. Fill colors use a red/blue theme in
// place of the source script's yellow/green.
// ============================================================

struct PolarimeterHistograms {
  TH1D *h_theta_fpp;
  TH1D *h_doca;

  TH1D *h_zclose_all;
  TH1D *h_zclose_sAng;
  TH1D *h_zclose_lAng;

  TH1D *h_dxp;
  TH1D *h_dyp;

  TH2D *h_dxpdyp;
  TH2D *h_dxpdyp_allth;

  TH2D *h_theta_vs_zclose;

  TH2D *h_xxp_ft;
  TH2D *h_xyp_ft;
  TH2D *h_yxp_ft;
  TH2D *h_yyp_ft;

  TH2D *h_xxp_fpp;
  TH2D *h_xyp_fpp;
  TH2D *h_yxp_fpp;
  TH2D *h_yyp_fpp;

  TH1D *h_chi2_ft;
  TH1D *h_chi2_fpp;
};

PolarimeterHistograms CreatePolarimeterHistograms() {
  PolarimeterHistograms h;

  h.h_theta_fpp = new TH1D("h_theta_fpp", "#theta_{FPP};#theta_{FPP} [deg];Counts", 100, 0, 10);
  h.h_doca      = new TH1D("h_doca", "DOCA;DOCA [cm];Counts", 50, 0, 2);

  h.h_zclose_all  = new TH1D("h_zclose_all",  "z_{close} all;z_{close} [m];Counts", 60, 0, 3.5);
  h.h_zclose_sAng = new TH1D("h_zclose_sAng", "z_{close} small angle;z_{close} [m];Counts", 60, 0, 3.5);
  h.h_zclose_lAng = new TH1D("h_zclose_lAng", "z_{close} large angle;z_{close} [m];Counts", 60, 0, 3.5);

  h.h_dxp = new TH1D("h_dxp_pol", "dxp;xp_{FT} - xp_{FPP} [deg];Counts", 400, -10, 11);
  h.h_dyp = new TH1D("h_dyp_pol", "dyp;yp_{FT} - yp_{FPP} [deg];Counts", 400, -10, 11);

  h.h_dxpdyp       = new TH2D("h_dxpdyp",       "dxp vs dyp;xp_{FT} - xp_{FPP} [deg];yp_{FT} - yp_{FPP} [deg]", 100, -5, 5, 100, -5, 5);
  h.h_dxpdyp_allth = new TH2D("h_dxpdyp_allth", "dxp vs dyp (No cut on #theta);xp_{FT} - xp_{FPP} [deg];yp_{FT} - yp_{FPP} [deg]", 100, -5, 5, 100, -5, 5);

  h.h_theta_vs_zclose = new TH2D("h_theta_vs_zclose", "#theta_{FPP} vs z_{close};z_{close} [m];#theta_{FPP} [deg]", 200, 0, 3.5, 200, 0, 10);

  h.h_xxp_ft = new TH2D("h_xxp_ft", "x vs xp FT;xp_{FT} [deg];x_{FT} [m]", 200, -16.0, 8.0, 200, -0.8, 0.8);
  h.h_xyp_ft = new TH2D("h_xyp_ft", "x vs yp FT;yp_{FT} [deg];x_{FT} [m]", 200, -4.0, 4.0, 200, -0.8, 0.8);
  h.h_yxp_ft = new TH2D("h_yxp_ft", "y vs xp FT;xp_{FT} [deg];y_{FT} [m]", 200, -16.0, 8.0, 200, -0.4, 0.4);
  h.h_yyp_ft = new TH2D("h_yyp_ft", "y vs yp FT;yp_{FT} [deg];y_{FT} [m]", 200, -4.0, 4.0, 200, -0.4, 0.4);

  h.h_xxp_fpp = new TH2D("h_xxp_fpp", "x vs xp FPP;xp_{FPP} [deg];x_{FPP} [m]", 200, -20.0, 10.0, 200, -1, 1);
  h.h_xyp_fpp = new TH2D("h_xyp_fpp", "x vs yp FPP;yp_{FPP} [deg];x_{FPP} [m]", 200, -8.0, 8.0, 200, -1.0, 1.0);
  h.h_yxp_fpp = new TH2D("h_yxp_fpp", "y vs xp FPP;xp_{FPP} [deg];y_{FPP} [m]", 200, -20.0, 10.0, 200, -0.6, 0.6);
  h.h_yyp_fpp = new TH2D("h_yyp_fpp", "y vs yp FPP;yp_{FPP} [deg];y_{FPP} [m]", 200, -8.0, 8.0, 200, -0.6, 0.6);

  h.h_chi2_ft  = new TH1D("h_chi2_ft_pol",  "chi2/ndf FT; chi2/ndf; Counts",  100, 0, 100);
  h.h_chi2_fpp = new TH1D("h_chi2_fpp_pol", "chi2/ndf FPP; chi2/ndf; Counts", 100, 0, 100);

  h.h_dxpdyp->SetContour(100);
  h.h_dxpdyp_allth->SetContour(100);
  h.h_theta_vs_zclose->SetContour(100);
  h.h_xxp_ft->SetContour(100);
  h.h_xyp_ft->SetContour(100);
  h.h_yxp_ft->SetContour(100);
  h.h_yyp_ft->SetContour(100);
  h.h_xxp_fpp->SetContour(100);
  h.h_xyp_fpp->SetContour(100);
  h.h_yxp_fpp->SetContour(100);
  h.h_yyp_fpp->SetContour(100);

  // Red/blue fill theme (source script used yellow/green)
  h.h_theta_fpp->SetLineColor(kBlack);
  h.h_theta_fpp->SetLineWidth(1);
  h.h_theta_fpp->SetFillColorAlpha(kBlue, 0.30);

  h.h_doca->SetLineColor(kBlack);
  h.h_doca->SetLineWidth(1);
  h.h_doca->SetMarkerColor(kBlack);
  h.h_doca->SetMarkerStyle(20);
  h.h_doca->SetMarkerSize(0.55);
  h.h_doca->SetFillStyle(0);

  h.h_dxp->SetLineColor(kBlack);
  h.h_dxp->SetLineWidth(1);
  h.h_dxp->SetFillColorAlpha(kRed, 0.30);

  h.h_dyp->SetLineColor(kBlack);
  h.h_dyp->SetLineWidth(1);
  h.h_dyp->SetFillColorAlpha(kBlue, 0.30);

  h.h_zclose_all->SetLineColor(kBlack);
  h.h_zclose_lAng->SetLineColor(kBlack);
  h.h_zclose_sAng->SetLineColor(kBlack);
  h.h_zclose_all->SetLineWidth(1);
  h.h_zclose_lAng->SetLineWidth(1);
  h.h_zclose_sAng->SetLineWidth(1);
  h.h_zclose_lAng->SetFillColor(kRed);
  h.h_zclose_sAng->SetFillColor(kBlue);

  return h;
}

void FillPolarimeterHistograms(const PolarimeterData &data, PolarimeterHistograms &h) {

  size_t n = data.theta_fpp_deg.size();

  for (size_t i = 0; i < n; i++) {

    double theta      = data.theta_fpp_deg[i];
    double sclose_cm   = data.sclose_m[i] * 100.0;
    double zclose      = data.zclose_m[i];

    h.h_theta_fpp->Fill(theta);
    h.h_doca->Fill(sclose_cm);

    bool close_track   = (data.sclose_m[i] < fpp_sclose_cut);
    bool near_target_z = (TMath::Abs(zclose - fpp_zclose_mean) < fpp_zclose_sigma);

    if (true) {
      h.h_zclose_all->Fill(zclose);

      if (theta > fpp_theta_min && theta < fpp_theta_max) {
        h.h_zclose_lAng->Fill(zclose);
      } else if (theta <= fpp_theta_min) {
        h.h_zclose_sAng->Fill(zclose);
      }

      h.h_theta_vs_zclose->Fill(zclose, theta);
    }

    if (true /*close_track && near_target_z*/) {
      h.h_dxp->Fill(data.dxp_deg[i]);
      h.h_dyp->Fill(data.dyp_deg[i]);
      h.h_dxpdyp->Fill(data.dxp_deg[i], data.dyp_deg[i]);
      h.h_dxpdyp_allth->Fill(data.dxp_deg[i], data.dyp_deg[i]);
    }

    if (true /*theta > fpp_theta_min && theta < fpp_theta_max && close_track && near_target_z*/) {
      h.h_xxp_ft->Fill(data.ft_xp_deg[i], data.ft_x[i]);
      h.h_xyp_ft->Fill(data.ft_yp_deg[i], data.ft_x[i]);
      h.h_yxp_ft->Fill(data.ft_xp_deg[i], data.ft_y[i]);
      h.h_yyp_ft->Fill(data.ft_yp_deg[i], data.ft_y[i]);

      h.h_xxp_fpp->Fill(data.fpp_xp_deg[i], data.fpp_x[i]);
      h.h_xyp_fpp->Fill(data.fpp_yp_deg[i], data.fpp_x[i]);
      h.h_yxp_fpp->Fill(data.fpp_xp_deg[i], data.fpp_y[i]);
      h.h_yyp_fpp->Fill(data.fpp_yp_deg[i], data.fpp_y[i]);
    }

    if (close_track) {
      h.h_chi2_ft->Fill(data.chi2_ft[i]);
      h.h_chi2_fpp->Fill(data.chi2_fpp[i]);
    }
  }
}

#endif
