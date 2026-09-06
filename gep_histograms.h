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
  hist.h_eresidu_FT_module  = new TH2D("h_eresidu_FT_module",  "eresidu FT; module; eresidu",  nmod_ft,  -0.5, nmod_ft-0.5,  100, -2, 2);
  hist.h_eresidv_FT_module  = new TH2D("h_eresidv_FT_module",  "eresidv FT; module; eresidv",  nmod_ft,  -0.5, nmod_ft-0.5,  100, -2, 2);
  hist.h_eresidu_FPP_module = new TH2D("h_eresidu_FPP_module", "eresidu FPP; module; eresidu", nmod_fpp, -0.5, nmod_fpp-0.5, 100, -2, 2);
  hist.h_eresidv_FPP_module = new TH2D("h_eresidv_FPP_module", "eresidv FPP; module; eresidv", nmod_fpp, -0.5, nmod_fpp-0.5, 100, -2, 2);

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
  }

  for (size_t i = 0; i < data.mod_fpp.size(); i++) {
    hist.h_eresidu_FPP_module->Fill(data.mod_fpp[i], data.eresidu_fpp[i]);
    hist.h_eresidv_FPP_module->Fill(data.mod_fpp[i], data.eresidv_fpp[i]);
  }
}

#endif
