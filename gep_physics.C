#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TString.h"
#include "TCut.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TROOT.h"
#include "TColor.h"
#include "TExec.h"
#include <algorithm>
#include <vector>
#include <cmath>

#include "gep_config.h"
#include "gep_fill_vectors.h"
#include "gep_histograms.h"
#include "gep_fit.h"
#include "gep_plot_utils.h"


// ============================================================
// Everything produced for one input set (files + cuts)
// ============================================================

struct SetResults {
  TString label;
  TString suffix;     // appended to histogram/function names ("" for set 1)

  GEPHistograms         hist;
  PolarimeterHistograms polhist;

  GEPFitResult fit_dx, fit_dy, fit_dxp, fit_dyp, fit_vz, fit_eresidu_FT, fit_eresidv_FT, fit_eresidu_FPP, fit_eresidv_FPP;

  TGraphErrors *g_eresidu_FT, *g_eresidv_FT, *g_eresidu_FPP, *g_eresidv_FPP, *g_eresidu_FT_layer, *g_eresidv_FT_layer, *g_eresidu_FPP_layer, *g_eresidv_FPP_layer;
};


// Fill, fit and normalize everything for one chain with its own cuts.
SetResults AnalyzeSet(TChain *C, TString label, TString suffix,
                      TCut cut, TCut cut_thetafpp) {

  SetResults r;
  r.label  = label;
  r.suffix = suffix;

  std::cout << "==== " << label << ": " << C->GetNtrees() << " file(s)" << std::endl;
  std::cout << "     cut:          " << cut.GetTitle() << std::endl;
  std::cout << "     cut_thetafpp: " << cut_thetafpp.GetTitle() << std::endl;

  // ----------------------------------------------------------
  // Fill vectors and histograms
  // ----------------------------------------------------------

  GEPData data;
  FillVectors(C, data, cut);

  r.hist = CreateHistograms(suffix);
  FillHistograms(data, r.hist);

  // Polarimeter reconstruction (separate cut, FPP-besttrack-indexed
  // pass over the same chain)
  PolarimeterData poldata;
  FillPolarimeterVectors(C, poldata, cut_thetafpp);

  r.polhist = CreatePolarimeterHistograms(suffix);
  FillPolarimeterHistograms(poldata, r.polhist);

  // ----------------------------------------------------------
  // Fits (on raw counts), per-module / per-layer mean +/- sigma
  // graphs, then per-column normalization for display
  // ----------------------------------------------------------

  r.fit_dx = FitPeak(r.hist.h_dx);
  r.fit_dy = FitPeak(r.hist.h_dy);
  r.fit_dxp = FitPeak(r.hist.h_dxp);
  r.fit_dyp = FitPeak(r.hist.h_dyp);
  r.fit_vz = FitPeak(r.hist.h_vz);

  // ==========================================================
  // Fit each module's residual distribution (mean +/- sigma).
  // Done on the raw (pre-normalization) 2D histograms so the fits
  // see real counts. Modules with no hits come back with
  // .valid == false and are skipped automatically -- nothing fails.
  // ==========================================================

  std::vector<GEPFitResult> fits_eresidu_FT = FitModuleColumns(r.hist.h_eresidu_FT_module,  nmod_ft);
  std::vector<GEPFitResult> fits_eresidv_FT = FitModuleColumns(r.hist.h_eresidv_FT_module,  nmod_ft);
  std::vector<GEPFitResult> fits_eresidu_FPP = FitModuleColumns(r.hist.h_eresidu_FPP_module, nmod_fpp);
  std::vector<GEPFitResult> fits_eresidv_FPP = FitModuleColumns(r.hist.h_eresidv_FPP_module, nmod_fpp);

  // Same per-column fit, but one column per layer (modules combined
  // into layers via layer_of_mod_ft/fpp in gep_config.h).
  std::vector<GEPFitResult> fits_eresidu_FT_layer = FitModuleColumns(r.hist.h_eresidu_FT_layer,  nlayer_ft);
  std::vector<GEPFitResult> fits_eresidv_FT_layer = FitModuleColumns(r.hist.h_eresidv_FT_layer,  nlayer_ft);
  std::vector<GEPFitResult> fits_eresidu_FPP_layer = FitModuleColumns(r.hist.h_eresidu_FPP_layer, nlayer_fpp);
  std::vector<GEPFitResult> fits_eresidv_FPP_layer = FitModuleColumns(r.hist.h_eresidv_FPP_layer, nlayer_fpp);

  // Fit the overall (all-modules-combined) U/V residual distributions
  r.fit_eresidu_FT = FitPeak(r.hist.h_eresidu_FT);
  r.fit_eresidv_FT = FitPeak(r.hist.h_eresidv_FT);
  r.fit_eresidu_FPP = FitPeak(r.hist.h_eresidu_FPP);
  r.fit_eresidv_FPP = FitPeak(r.hist.h_eresidv_FPP);

  r.g_eresidu_FT = BuildResidualGraph(fits_eresidu_FT);
  r.g_eresidv_FT = BuildResidualGraph(fits_eresidv_FT);
  r.g_eresidu_FPP = BuildResidualGraph(fits_eresidu_FPP);
  r.g_eresidv_FPP = BuildResidualGraph(fits_eresidv_FPP);

  r.g_eresidu_FT_layer = BuildResidualGraph(fits_eresidu_FT_layer);
  r.g_eresidv_FT_layer = BuildResidualGraph(fits_eresidv_FT_layer);
  r.g_eresidu_FPP_layer = BuildResidualGraph(fits_eresidu_FPP_layer);
  r.g_eresidv_FPP_layer = BuildResidualGraph(fits_eresidv_FPP_layer);

  // ==========================================================
  // Normalize module-wise 2D residual histograms for display
  // (each module column scaled to its own peak bin). Purely
  // cosmetic -- done after fitting so it doesn't affect the fits.
  // ==========================================================

  NormalizeModuleColumns(r.hist.h_eresidu_FT_module);
  NormalizeModuleColumns(r.hist.h_eresidv_FT_module);
  NormalizeModuleColumns(r.hist.h_eresidu_FPP_module);
  NormalizeModuleColumns(r.hist.h_eresidv_FPP_module);

  NormalizeModuleColumns(r.hist.h_eresidu_FT_layer);
  NormalizeModuleColumns(r.hist.h_eresidv_FT_layer);
  NormalizeModuleColumns(r.hist.h_eresidu_FPP_layer);
  NormalizeModuleColumns(r.hist.h_eresidv_FPP_layer);

  return r;
}


void WriteSet(const SetResults &r) {
  r.hist.h_dx->Write();
  r.hist.h_dy->Write();
  r.hist.h_dxp->Write();
  r.hist.h_dyp->Write();
  r.hist.h_vz->Write();
  r.hist.h_vx->Write();
  r.hist.h_vy->Write();
  r.hist.h_vxvy->Write();
  r.hist.h_eresidu_FT_module->Write();
  r.hist.h_eresidv_FT_module->Write();
  r.hist.h_eresidu_FPP_module->Write();
  r.hist.h_eresidv_FPP_module->Write();
  r.hist.h_eresidu_FT_layer->Write();
  r.hist.h_eresidv_FT_layer->Write();
  r.hist.h_eresidu_FPP_layer->Write();
  r.hist.h_eresidv_FPP_layer->Write();
  r.hist.h_eresidu_FT->Write();
  r.hist.h_eresidv_FT->Write();
  r.hist.h_eresidu_FPP->Write();
  r.hist.h_eresidv_FPP->Write();
}


// Small set label in the top-left corner of the canvas
void StampSetLabel(TCanvas *c, const TString &label) {
  c->cd();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.022);
  t.SetTextAlign(13);
  t.SetTextColor(kGray + 3);
  t.DrawLatex(0.003, 0.997, label);
}


// All per-set pages (alignment, target, module-wise residuals,
// layer-wise residuals, polarimeter kinematics, FT/FPP correlations,
// chi2/ndf), each stamped with the set label.
void DrawSetPages(TCanvas *c1, SetResults &r, const char *pdfname) {

  // Displayed residual range for all 2D residual maps (fits and
  // normalization were already done on the full range)
  TH2D *resid_maps[8] = { r.hist.h_eresidu_FT_module, r.hist.h_eresidv_FT_module,
                          r.hist.h_eresidu_FPP_module, r.hist.h_eresidv_FPP_module,
                          r.hist.h_eresidu_FT_layer, r.hist.h_eresidv_FT_layer,
                          r.hist.h_eresidu_FPP_layer, r.hist.h_eresidv_FPP_layer };
  for (TH2D *h : resid_maps) h->GetYaxis()->SetRangeUser(resid_plot_min, resid_plot_max);

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  //---------- Check Alignment

  c1->Clear();
  c1->SetCanvasSize(900, 700);
  c1->SetLeftMargin(gStyle->GetPadLeftMargin());
  c1->SetRightMargin(gStyle->GetPadRightMargin());
  c1->SetBottomMargin(gStyle->GetPadBottomMargin());
  c1->SetTopMargin(gStyle->GetPadTopMargin());
  c1->Divide(2,2);

  c1->cd(1);
  r.hist.h_dx->Draw();

  DrawTextBox({
  "Cuts:",
  "FT: N_{hits} > 4 || N_{goodhits} > 2",
  "FPP: N_{hits} > 4 || N_{goodhits} > 2",
  "FT: #chi^{2}/ndf < 200",
  "FPP: #chi^{2}/ndf < 200",
  "HCAL: N_{blk} > 1",
  "TAR: |vz+0.09494|<3*0.02332"
  }, 0.12, 0.60, 0.42, 0.88);

  c1->cd(2);
  r.hist.h_dy->Draw();

  c1->cd(3);
  r.hist.h_dxp->Draw();

  c1->cd(4);
  r.hist.h_dyp->Draw();

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);

  //---------- Check Target

  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1);
  r.hist.h_vz->Draw();

  DrawTextBox({
  "Cuts:",
  "FT: N_{hits} > 4 || N_{goodhits} > 2",
  "FPP: N_{hits} > 4 || N_{goodhits} > 2",
  "FT: #chi^{2}/ndf < 200",
  "FPP: #chi^{2}/ndf < 200",
  "HCAL: N_{blk} > 1"
  }, 0.12, 0.60, 0.42, 0.88);

  c1->cd(2);
  r.hist.h_vx->Draw();

  c1->cd(3);
  r.hist.h_vy->Draw();

  c1->cd(4);
  r.hist.h_vxvy->Draw();

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);

  //---------- Module-wise residuals: FT u/v and FPP u/v, plus overall
  //---------- (all-module) U/V residual distributions with fits, all
  //---------- on one page. Row 1 = FT (module map U, module map V,
  //---------- overall U, overall V), row 2 = same layout for FPP.

  c1->Clear();
  c1->SetCanvasSize(1800, 900);
  c1->Divide(4,2);

  gStyle->SetPalette(kRainBow);

  // Custom, larger TPaveText fit-stats boxes (like the polarimeter
  // kinematics page) replace the default ROOT stat/fit box on this
  // page's 1D overall-residual plots. The 2D module-map (COLZ) plots
  // get no stats box at all, so OptStat/OptFit stay off for the whole
  // page and the 1D pads draw their box manually.
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  c1->cd(1);
  r.hist.h_eresidu_FT_module->SetMinimum(0);
  r.hist.h_eresidu_FT_module->SetMaximum(1);
  r.hist.h_eresidu_FT_module->SetStats(0);
  r.hist.h_eresidu_FT_module->Draw("COLZ");
  r.g_eresidu_FT->Draw("P SAME");

  c1->cd(2);
  r.hist.h_eresidv_FT_module->SetMinimum(0);
  r.hist.h_eresidv_FT_module->SetMaximum(1);
  r.hist.h_eresidv_FT_module->SetStats(0);
  r.hist.h_eresidv_FT_module->Draw("COLZ");
  r.g_eresidv_FT->Draw("P SAME");

  c1->cd(3);
  r.hist.h_eresidu_FT->SetStats(0);
  r.hist.h_eresidu_FT->Draw();
  if (r.fit_eresidu_FT.valid) {
    MakeFitStatsBoxFromResult(r.hist.h_eresidu_FT, r.fit_eresidu_FT, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(4);
  r.hist.h_eresidv_FT->SetStats(0);
  r.hist.h_eresidv_FT->Draw();
  if (r.fit_eresidv_FT.valid) {
    MakeFitStatsBoxFromResult(r.hist.h_eresidv_FT, r.fit_eresidv_FT, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(5);
  r.hist.h_eresidu_FPP_module->SetMinimum(0);
  r.hist.h_eresidu_FPP_module->SetMaximum(1);
  r.hist.h_eresidu_FPP_module->SetStats(0);
  r.hist.h_eresidu_FPP_module->Draw("COLZ");
  r.g_eresidu_FPP->Draw("P SAME");

  c1->cd(6);
  r.hist.h_eresidv_FPP_module->SetMinimum(0);
  r.hist.h_eresidv_FPP_module->SetMaximum(1);
  r.hist.h_eresidv_FPP_module->SetStats(0);
  r.hist.h_eresidv_FPP_module->Draw("COLZ");
  r.g_eresidv_FPP->Draw("P SAME");

  c1->cd(7);
  r.hist.h_eresidu_FPP->SetStats(0);
  r.hist.h_eresidu_FPP->Draw();
  if (r.fit_eresidu_FPP.valid) {
    MakeFitStatsBoxFromResult(r.hist.h_eresidu_FPP, r.fit_eresidu_FPP, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(8);
  r.hist.h_eresidv_FPP->SetStats(0);
  r.hist.h_eresidv_FPP->Draw();
  if (r.fit_eresidv_FPP.valid) {
    MakeFitStatsBoxFromResult(r.hist.h_eresidv_FPP, r.fit_eresidv_FPP, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);

  //---------- Layer-wise residuals: one 2D plot per page
  //---------- (FT u, FT v, FPP u, FPP v). Same per-column
  //---------- normalization and mean +/- sigma overlay as the
  //---------- module-wise maps above.

  {
    struct LayerPage { TH2D *h; TGraphErrors *g; const char *modmap; };
    LayerPage layer_pages[4] = {
      { r.hist.h_eresidu_FT_layer,  r.g_eresidu_FT_layer,  "FT: L0-L5 = m0-m5, L6 = m6-m9, L7 = m10-m13" },
      { r.hist.h_eresidv_FT_layer,  r.g_eresidv_FT_layer,  "FT: L0-L5 = m0-m5, L6 = m6-m9, L7 = m10-m13" },
      { r.hist.h_eresidu_FPP_layer, r.g_eresidu_FPP_layer, "FPP: L_{n} = m_{4n} - m_{4n+3}" },
      { r.hist.h_eresidv_FPP_layer, r.g_eresidv_FPP_layer, "FPP: L_{n} = m_{4n} - m_{4n+3}" }
    };

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kRainBow);

    for (int ip = 0; ip < 4; ip++) {
      c1->Clear();
      c1->SetCanvasSize(1200, 900);
      c1->cd();
      c1->SetLeftMargin(0.11);
      c1->SetRightMargin(0.14);
      c1->SetBottomMargin(0.11);
      c1->SetTopMargin(0.08);

      TH2D *h = layer_pages[ip].h;
      h->SetMinimum(0);
      h->SetMaximum(1);
      h->SetStats(0);
      h->Draw("COLZ");
      layer_pages[ip].g->Draw("P SAME");

      DrawTextBox({ layer_pages[ip].modmap }, 0.13, 0.86, 0.60, 0.91);

      StampSetLabel(c1, r.label);
  c1->Print(pdfname);
    }

    // Restore canvas margins for the pages that follow
    c1->SetLeftMargin(gStyle->GetPadLeftMargin());
    c1->SetRightMargin(gStyle->GetPadRightMargin());
    c1->SetBottomMargin(gStyle->GetPadBottomMargin());
    c1->SetTopMargin(gStyle->GetPadTopMargin());
  }

  //---------- Polarimeter kinematics: theta_FPP, DOCA, z_close,
  //---------- theta vs z_close, dxp, dyp, dxp-vs-dyp. Ported from
  //---------- polarimeter_recon.C, keeping that script's own
  //---------- OptStat(0)/manual-stats-box convention and a golden
  //---------- -ratio canvas, with a red/blue fill theme in place of
  //---------- its original yellow/green.

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const int pol_canvas_width  = 2400 * 2;
  const int pol_canvas_height = int(pol_canvas_width / ((1.0 + TMath::Sqrt(5.0)) / 2.0));

  c1->Clear();
  c1->SetCanvasSize(pol_canvas_width, pol_canvas_height);
  c1->Divide(4, 2, 0.001, 0.001);

  for (int ipad = 1; ipad <= 8; ipad++) {
    c1->cd(ipad);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.035);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.075);
  }

  c1->cd(1);
  gPad->SetLogy();

  TLine *lfpp_theta1 = new TLine(fpp_theta_min, 0, fpp_theta_min, r.polhist.h_theta_fpp->GetMaximum());
  TLine *lfpp_theta2 = new TLine(fpp_theta_max, 0, fpp_theta_max, r.polhist.h_theta_fpp->GetMaximum());
  lfpp_theta1->SetLineWidth(1);
  lfpp_theta2->SetLineWidth(1);
  lfpp_theta1->SetLineColor(kRed+1);
  lfpp_theta2->SetLineColor(kRed+1);

  r.polhist.h_theta_fpp->Draw("hist");
  //lfpp_theta1->Draw("SAME");
  //lfpp_theta2->Draw("SAME");

  MakeStatsBox(r.polhist.h_theta_fpp, 0.62, 0.73, 0.91, 0.90)->Draw();

  c1->cd(2);
  gPad->SetLogy(0);
  r.polhist.h_doca->Draw("hist");

  double doca_fit_min = 0.0;
  double doca_fit_max = 0.2;
  TF1 *f_halfgaus = new TF1(Form("f_halfgaus%s", r.suffix.Data()), "[0]*exp(-0.5*x*x/([1]*[1]))", doca_fit_min, doca_fit_max);
  f_halfgaus->SetParNames("A", "#sigma");
  f_halfgaus->SetParameters(r.polhist.h_doca->GetMaximum(), 0.10);
  f_halfgaus->SetParLimits(0, 1e-6, 1e9);
  f_halfgaus->SetParLimits(1, 1e-5, 5.0);
  f_halfgaus->SetLineColor(kBlue+2);
  f_halfgaus->SetLineWidth(1);
  r.polhist.h_doca->Fit(f_halfgaus, "RQ0");

  double doca_ymax_hist = r.polhist.h_doca->GetMaximum();
  double doca_ymax_fit  = f_halfgaus->GetMaximum(doca_fit_min, doca_fit_max);
  r.polhist.h_doca->SetMaximum(1.2 * (doca_ymax_fit > doca_ymax_hist ? doca_ymax_fit : doca_ymax_hist));

  r.polhist.h_doca->Draw("hist");
  f_halfgaus->Draw("SAME");

  TLine *lfpp_sclose = new TLine(0.5, 0, 0.5, r.polhist.h_doca->GetMaximum());
  lfpp_sclose->SetLineWidth(1);
  lfpp_sclose->SetLineColor(kRed+1);
  //lfpp_sclose->Draw("SAME");

  TPaveText *pt_doca = new TPaveText(0.37, 0.61, 0.81, 0.88, "NDC");
  pt_doca->SetFillColorAlpha(kWhite, 0.88);
  pt_doca->SetBorderSize(0);
  pt_doca->SetTextFont(42);
  pt_doca->SetTextSize(0.05);
  pt_doca->SetTextAlign(12);
  pt_doca->SetMargin(0.06);
  pt_doca->AddText(Form("Range: %.3f - %.3f cm", doca_fit_min, doca_fit_max));
  pt_doca->AddText(Form("A = %.3g #pm %.3g", f_halfgaus->GetParameter(0), f_halfgaus->GetParError(0)));
  pt_doca->AddText(Form("#sigma = %.5f #pm %.5f cm", f_halfgaus->GetParameter(1), f_halfgaus->GetParError(1)));
  if (f_halfgaus->GetNDF() > 0) {
    pt_doca->AddText(Form("#chi^{2}/NDF = %.2f / %d = %.2f",
                          f_halfgaus->GetChisquare(), f_halfgaus->GetNDF(),
                          f_halfgaus->GetChisquare() / f_halfgaus->GetNDF()));
  }
  pt_doca->AddText(Form("N = %.0f", r.polhist.h_doca->GetEntries()));
  pt_doca->Draw("SAME");

  c1->cd(3);
  r.polhist.h_zclose_all->Draw("hist");
  r.polhist.h_zclose_sAng->Draw("hist same");
  r.polhist.h_zclose_lAng->Draw("hist same");

  TLine *lmin_zclose = new TLine(fpp_zclose_mean - fpp_zclose_sigma, 0,
                                 fpp_zclose_mean - fpp_zclose_sigma, r.polhist.h_zclose_all->GetMaximum());
  TLine *lmax_zclose = new TLine(fpp_zclose_mean + fpp_zclose_sigma, 0,
                                 fpp_zclose_mean + fpp_zclose_sigma, r.polhist.h_zclose_all->GetMaximum());
  lmin_zclose->SetLineWidth(1);
  lmax_zclose->SetLineWidth(1);
  lmin_zclose->SetLineColor(kRed+1);
  lmax_zclose->SetLineColor(kRed+1);
  lmin_zclose->Draw("SAME");
  lmax_zclose->Draw("SAME");

  TLegend *leg_zclose = new TLegend(0.59, 0.68, 0.93, 0.92);
  leg_zclose->SetBorderSize(0);
  leg_zclose->SetFillColorAlpha(kWhite, 0.88);
  leg_zclose->SetTextSize(0.05);
  leg_zclose->AddEntry(r.polhist.h_zclose_all, "all", "l");
  leg_zclose->AddEntry(r.polhist.h_zclose_lAng, Form("#theta_{FPP} > %.2f", fpp_theta_min), "lf");
  leg_zclose->AddEntry(r.polhist.h_zclose_sAng, Form("#theta_{FPP} <= %.2f", fpp_theta_min), "lf");
  leg_zclose->Draw();

  //MakeStatsBox(r.polhist.h_zclose_all, 0.59, 0.47, 0.93, 0.65)->Draw();

  c1->cd(4);
  gPad->SetRightMargin(0.16);
  r.polhist.h_theta_vs_zclose->SetStats(0);
  r.polhist.h_theta_vs_zclose->Draw("COLZ");

  //------------------------------------------------------

  c1->cd(5);
  r.polhist.h_dxp->GetXaxis()->SetRangeUser(-2.0, 4.0);
  r.polhist.h_dxp->Draw("hist");

  const double peak_fit_half_width = 0.2; // degrees
  double dxp_fit_mean = r.polhist.h_dxp->GetXaxis()->GetBinCenter(r.polhist.h_dxp->GetMaximumBin());
  double dxp_fit_rms  = r.polhist.h_dxp->GetRMS();
  TF1 *fgaus_dxp = new TF1(Form("fgaus_dxp%s", r.suffix.Data()), "gaus",
                           dxp_fit_mean - peak_fit_half_width,
                           dxp_fit_mean + 1.2*peak_fit_half_width);
  fgaus_dxp->SetParameters(r.polhist.h_dxp->GetMaximum(), dxp_fit_mean, 0.4);
  fgaus_dxp->SetLineColor(kRed+2);
  fgaus_dxp->SetLineWidth(1);

  bool dxp_fit_ok = r.polhist.h_dxp->GetEntries() > 3 && dxp_fit_rms > 0.0;
  if (dxp_fit_ok) {
    r.polhist.h_dxp->Fit(fgaus_dxp, "RQ0");
    fgaus_dxp->Draw("SAME");
  }

  TLine *ldxp = new TLine(0, 0, 0, r.polhist.h_dxp->GetMaximum());
  ldxp->SetLineWidth(1);
  ldxp->SetLineColor(kRed);
  ldxp->Draw("same");

  if (dxp_fit_ok) {
    MakeFitStatsBox(r.polhist.h_dxp, fgaus_dxp, 0.53, 0.58, 0.94, 0.91)->Draw();
  } else {
    MakeStatsBox(r.polhist.h_dxp, 0.62, 0.73, 0.91, 0.90)->Draw();
  }

  //------------------------------------------------------

  c1->cd(6);
  r.polhist.h_dyp->GetXaxis()->SetRangeUser(-2.0, 4.0);
  r.polhist.h_dyp->Draw("hist");

  double dyp_fit_mean = r.polhist.h_dyp->GetXaxis()->GetBinCenter(r.polhist.h_dyp->GetMaximumBin());
  double dyp_fit_rms  = r.polhist.h_dyp->GetRMS();
  TF1 *fgaus_dyp = new TF1(Form("fgaus_dyp%s", r.suffix.Data()), "gaus",
                           dyp_fit_mean - 1.2*peak_fit_half_width,
                           dyp_fit_mean + peak_fit_half_width);
  fgaus_dyp->SetParameters(r.polhist.h_dyp->GetMaximum(), dyp_fit_mean, 0.4);
  fgaus_dyp->SetLineColor(kBlue+2);
  fgaus_dyp->SetLineWidth(1);

  bool dyp_fit_ok = r.polhist.h_dyp->GetEntries() > 3 && dyp_fit_rms > 0.0;
  if (dyp_fit_ok) {
    r.polhist.h_dyp->Fit(fgaus_dyp, "RQ0");
    fgaus_dyp->Draw("SAME");
  }

  TLine *ldyp = new TLine(0, 0, 0, r.polhist.h_dyp->GetMaximum());
  ldyp->SetLineWidth(1);
  ldyp->SetLineColor(kRed);
  ldyp->Draw("same");

  if (dyp_fit_ok) {
    MakeFitStatsBox(r.polhist.h_dyp, fgaus_dyp, 0.53, 0.58, 0.94, 0.91)->Draw();
  } else {
    MakeStatsBox(r.polhist.h_dyp, 0.62, 0.73, 0.91, 0.90)->Draw();
  }

  //------------------------------------------------------

  c1->cd(7);
  gPad->SetRightMargin(0.16);
  r.polhist.h_dxpdyp->SetStats(0);
  r.polhist.h_dxpdyp->Draw("COLZ");

  {
    double xmin = r.polhist.h_dxpdyp->GetXaxis()->GetXmin();
    double xmax = r.polhist.h_dxpdyp->GetXaxis()->GetXmax();
    double ymin = r.polhist.h_dxpdyp->GetYaxis()->GetXmin();
    double ymax = r.polhist.h_dxpdyp->GetYaxis()->GetXmax();

    TLine *vertical   = new TLine(0.0, ymin, 0.0, ymax);
    TLine *horizontal = new TLine(xmin, 0.0, xmax, 0.0);
    vertical->SetLineColor(kRed);
    vertical->SetLineWidth(1);
    vertical->SetLineStyle(2);
    horizontal->SetLineColor(kRed);
    horizontal->SetLineWidth(1);
    horizontal->SetLineStyle(2);
    vertical->Draw("SAME");
    horizontal->Draw("SAME");
  }

  c1->cd(8); // left blank -- matches polarimeter_recon.C, where
             // hdxpdyp_allth is filled but never drawn.

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);

  //---------- FT/FPP position-vs-slope correlation maps

  c1->Clear();
  c1->Divide(4, 2, 0.001, 0.001);

  for (int ipad = 1; ipad <= 8; ipad++) {
    c1->cd(ipad);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.16);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.075);
  }

  c1->cd(1);
  r.polhist.h_xxp_ft->SetStats(0);
  r.polhist.h_xxp_ft->Draw("COLZ");

  c1->cd(2);
  r.polhist.h_xyp_ft->SetStats(0);
  r.polhist.h_xyp_ft->Draw("COLZ");

  c1->cd(3);
  r.polhist.h_yxp_ft->SetStats(0);
  r.polhist.h_yxp_ft->Draw("COLZ");

  c1->cd(4);
  r.polhist.h_yyp_ft->SetStats(0);
  r.polhist.h_yyp_ft->Draw("COLZ");

  c1->cd(5);
  r.polhist.h_xxp_fpp->SetStats(0);
  r.polhist.h_xxp_fpp->Draw("COLZ");

  c1->cd(6);
  r.polhist.h_xyp_fpp->SetStats(0);
  r.polhist.h_xyp_fpp->Draw("COLZ");

  c1->cd(7);
  r.polhist.h_yxp_fpp->SetStats(0);
  r.polhist.h_yxp_fpp->Draw("COLZ");

  c1->cd(8);
  r.polhist.h_yyp_fpp->SetStats(0);
  r.polhist.h_yyp_fpp->Draw("COLZ");

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);

  //---------- FT/FPP chi2/ndf

  c1->Clear();
  c1->SetCanvasSize(2000, 900);
  c1->Divide(2, 1);

  gStyle->SetOptStat(1111);

  c1->cd(1);
  r.polhist.h_chi2_ft->SetStats(1);
  r.polhist.h_chi2_ft->Draw();

  c1->cd(2);
  r.polhist.h_chi2_fpp->SetStats(1);
  r.polhist.h_chi2_fpp->Draw();

  StampSetLabel(c1, r.label);
  c1->Print(pdfname);
}


// ============================================================
// Polarimeter comparison pages (16:9, 1920x1080): rows = sets
// (set 1 on top, set 2 below), columns = the plots of one row of the
// per-set polarimeter page. No stats boxes; only the z_close legend.
// ============================================================

struct RowLayout {
  double label_h, gap, row_h;
  double rowy1[2], rowy2[2];
};

RowLayout SetupRowComparisonCanvas(TCanvas *c, const char *label1, const char *label2) {
  RowLayout L;
  L.label_h  = 0.040;
  L.gap      = 0.014;
  L.row_h    = (1.0 - 2.0 * L.label_h - L.gap) / 2.0;
  L.rowy2[0] = 1.0 - L.label_h;
  L.rowy1[0] = L.rowy2[0] - L.row_h;
  L.rowy2[1] = L.rowy1[0] - L.gap - L.label_h;
  L.rowy1[1] = 0.0;   // exact: (1 - ...) - row_h can round to -1e-17,
                      // and TPad rejects any edge outside [0,1]

  c->Clear();
  c->SetCanvasSize(1920, 1080);
  c->SetFillColor(kWhite);
  c->cd();

  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.028);
  t.SetTextAlign(22);
  t.DrawLatex(0.5, L.rowy2[0] + 0.5 * L.label_h, label1);
  t.DrawLatex(0.5, L.rowy2[1] + 0.5 * L.label_h, label2);

  // Divider between the two sets
  double ydiv = L.rowy1[0] - 0.5 * L.gap;
  TLine *divider = new TLine(0.0, ydiv, 1.0, ydiv);
  divider->SetNDC();
  divider->SetLineColor(kGray + 2);
  divider->SetLineWidth(3);
  divider->Draw();

  return L;
}

TPad *MakeRowPad(TCanvas *c, const TString &name, int icol, int ncol,
                 const RowLayout &L, int iset, bool colz) {
  c->cd();
  double w = 1.0 / ncol;
  auto clamp01 = [](double v) { return std::min(1.0, std::max(0.0, v)); };
  TPad *pad = new TPad(name, "", clamp01(icol * w), clamp01(L.rowy1[iset]),
                       clamp01((icol + 1) * w), clamp01(L.rowy2[iset]));
  pad->SetLeftMargin(0.14);
  pad->SetRightMargin(colz ? 0.14 : 0.04);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->Draw();
  pad->cd();
  return pad;
}


// Opacity of the set-1 (top) row on the polarimeter comparison pages
const double polcmp_set1_alpha = 0.90;

// Scale line/fill/marker opacity of a 1D histogram
void FadeHist(TH1 *h, double alpha) {
  auto a = [&](Color_t idx) {
    TColor *col = gROOT->GetColor(idx);
    return (col ? col->GetAlpha() : 1.0) * alpha;
  };
  h->SetLineColorAlpha(h->GetLineColor(), a(h->GetLineColor()));
  if (h->GetFillColor() != 0 && h->GetFillStyle() != 0)   // leave hollow histograms hollow
    h->SetFillColorAlpha(h->GetFillColor(), a(h->GetFillColor()));
  h->SetMarkerColorAlpha(h->GetMarkerColor(), a(h->GetMarkerColor()));
}

// Draw a 2D map with COLZ in kRainBow at the given opacity, using a
// per-pad TExec so pads with different opacities coexist on a page
// (axes-only first pass so no opaque copy sits underneath).
void DrawColzAlpha(TH2 *h, double alpha, const TString &uid) {
  TExec *ex = new TExec("palexec_" + uid,
                        TString::Format("gStyle->SetPalette(kRainBow, 0, %.3f);", alpha));
  h->SetContour(99);
  h->Draw("AXIS");
  ex->Draw();
  h->Draw("COLZ SAME");
  h->Draw("AXIS SAME");
}

// Small borderless box with a few lines of text (pad NDC)
void DrawValueBox(const std::vector<TString> &lines, double x1, double y1, double x2, double y2) {
  TPaveText *box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillColorAlpha(kWhite, 0.85);
  box->SetBorderSize(0);
  box->SetTextFont(42);
  box->SetTextSize(0.045);
  box->SetTextAlign(12);
  for (const auto &l : lines) box->AddText(l);
  box->Draw();
}

// Gaussian fit to the core of a dxp/dyp peak (same windows as the
// per-set polarimeter page). Not drawn; returns false if not fit.
bool FitPeakCore(TH1D *h, double lo_w, double hi_w, const TString &name,
                 double &mu, double &sigma) {
  if (h->GetEntries() <= 3 || h->GetRMS() <= 0) return false;
  double m = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  TF1 *f = new TF1(name, "gaus", m - lo_w, m + hi_w);
  f->SetParameters(h->GetMaximum(), m, 0.4);
  h->Fit(f, "RQ0");
  mu = f->GetParameter(1);
  sigma = std::fabs(f->GetParameter(2));
  return true;
}

// Clone with no stats box and no attached fit functions
template <class H>
H *CleanClone(H *src, const TString &name) {
  H *h = (H *)src->Clone(name);
  h->SetStats(0);
  if (q2_label && q2_label[0] != '\0')
    h->SetTitle(TString::Format("%s, %s", src->GetTitle(), q2_label));
  h->GetListOfFunctions()->Clear();
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleSize(0.050);
  h->GetYaxis()->SetTitleSize(0.050);
  h->GetYaxis()->SetTitleOffset(1.35);
  return h;
}

// Page 1: theta_FPP | DOCA | z_close | theta_FPP vs z_close
void DrawPolarimeterComparisonRow1(TCanvas *c, SetResults *sets[2],
                                   const char *label1, const char *label2,
                                   const char *pdfname) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(kRainBow);

  RowLayout L = SetupRowComparisonCanvas(c, label1, label2);

  for (int is = 0; is < 2; is++) {
    PolarimeterHistograms &ph = sets[is]->polhist;
    TString sfx = TString::Format("_cmpP1_%d", is);

    // theta_FPP (log y)
    TPad *p1 = MakeRowPad(c, "pth" + sfx, 0, 4, L, is, false);
    p1->SetLogy();
    const double alpha = (is == 0) ? polcmp_set1_alpha : 1.0;

    TH1D *hth = CleanClone(ph.h_theta_fpp, "hth" + sfx);
    if (alpha < 1.0) FadeHist(hth, alpha);
    hth->Draw("hist");

    // DOCA with the half-Gaussian curve, no stats box
    MakeRowPad(c, "pdoca" + sfx, 1, 4, L, is, false);
    TH1D *hdoca = CleanClone(ph.h_doca, "hdoca" + sfx);
    if (alpha < 1.0) FadeHist(hdoca, alpha);
    TF1 *fdoca = new TF1("fdoca" + sfx, "[0]*exp(-0.5*x*x/([1]*[1]))", 0.0, 0.2);
    fdoca->SetParameters(hdoca->GetMaximum(), 0.10);
    fdoca->SetParLimits(0, 1e-6, 1e9);
    fdoca->SetParLimits(1, 1e-5, 5.0);
    fdoca->SetLineColorAlpha(kBlue + 2, alpha);
    fdoca->SetLineWidth(3);
    hdoca->SetLineWidth(2);
    if (hdoca->GetEntries() > 3) hdoca->Fit(fdoca, "RQ0");
    double ymax = std::max(hdoca->GetMaximum(), fdoca->GetMaximum(0.0, 0.2));
    hdoca->SetMaximum(1.2 * ymax);
    hdoca->Draw("hist");
    if (hdoca->GetEntries() > 3) {
      fdoca->Draw("SAME");
      DrawValueBox({ TString::Format("Mean = %.4f cm", hdoca->GetMean()),
                     TString::Format("#sigma = %.4f cm", fdoca->GetParameter(1)) },
                   0.50, 0.74, 0.94, 0.90);
    }

    // z_close all / small / large angle, with legend
    MakeRowPad(c, "pz" + sfx, 2, 4, L, is, false);
    TH1D *hzall = CleanClone(ph.h_zclose_all,  "hzall" + sfx);
    TH1D *hzs   = CleanClone(ph.h_zclose_sAng, "hzs"   + sfx);
    TH1D *hzl   = CleanClone(ph.h_zclose_lAng, "hzl"   + sfx);
    if (alpha < 1.0) { FadeHist(hzall, alpha); FadeHist(hzs, alpha); FadeHist(hzl, alpha); }
    hzl->SetFillColorAlpha(kRed,  0.75);   // large-angle fill
    hzs->SetFillColorAlpha(kBlue, 0.75);   // small-angle fill
    hzall->Draw("hist");
    hzs->Draw("hist same");
    hzl->Draw("hist same");
    TLegend *leg = new TLegend(0.58, 0.68, 0.95, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite, 0.88);
    leg->SetTextSize(0.045);
    leg->AddEntry(hzall, "all", "l");
    leg->AddEntry(hzl, Form("#theta_{FPP} > %.2f", fpp_theta_min), "lf");
    leg->AddEntry(hzs, Form("#theta_{FPP} <= %.2f", fpp_theta_min), "lf");
    leg->Draw();

    // theta_FPP vs z_close
    MakeRowPad(c, "ptz" + sfx, 3, 4, L, is, true);
    TH2D *htz = CleanClone(ph.h_theta_vs_zclose, "htz" + sfx);
    DrawColzAlpha(htz, alpha, "htz" + sfx);
  }

  c->cd();
  c->Update();
  c->Print(pdfname);
  gStyle->SetPalette(kRainBow);
}

// Page 2: dxp | dyp | dxp vs dyp -- no stats, no fits, ranges centred
// on zero (-3..3 deg; both axes of the 2D plot).
void DrawPolarimeterComparisonRow2(TCanvas *c, SetResults *sets[2],
                                   const char *label1, const char *label2,
                                   const char *pdfname) {
  const double amin = -3.0, amax = 3.0;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(kRainBow);

  RowLayout L = SetupRowComparisonCanvas(c, label1, label2);

  for (int is = 0; is < 2; is++) {
    PolarimeterHistograms &ph = sets[is]->polhist;
    TString sfx = TString::Format("_cmpP2_%d", is);

    const double alpha = (is == 0) ? polcmp_set1_alpha : 1.0;

    TH1D *h1[2] = { CleanClone(ph.h_dxp, "hdxp" + sfx), CleanClone(ph.h_dyp, "hdyp" + sfx) };
    for (int k = 0; k < 2; k++) {
      MakeRowPad(c, TString::Format("pd%d", k) + sfx, k, 3, L, is, false);
      h1[k]->GetXaxis()->SetRangeUser(amin, amax);
      if (alpha < 1.0) FadeHist(h1[k], alpha);
      h1[k]->Draw("hist");
      gPad->Update();
      TLine *l0 = new TLine(0, 0, 0, h1[k]->GetMaximum());
      l0->SetLineWidth(1);
      l0->SetLineColor(kRed);
      l0->Draw("same");

      double mu = 0, sg = 0;
      bool ok = (k == 0) ? FitPeakCore(h1[k], 0.2, 0.24, TString::Format("fcore%d", k) + sfx, mu, sg)
                         : FitPeakCore(h1[k], 0.24, 0.2, TString::Format("fcore%d", k) + sfx, mu, sg);
      if (!ok) { mu = h1[k]->GetMean(); sg = h1[k]->GetStdDev(); }
      DrawValueBox({ TString::Format("Mean = %.3f deg", mu),
                     TString::Format("#sigma = %.3f deg", sg) },
                   0.58, 0.76, 0.95, 0.90);
    }

    MakeRowPad(c, "pdd" + sfx, 2, 3, L, is, true);
    TH2D *h2 = CleanClone(ph.h_dxpdyp, "hdxpdyp" + sfx);
    h2->GetXaxis()->SetRangeUser(amin, amax);
    h2->GetYaxis()->SetRangeUser(amin, amax);
    DrawColzAlpha(h2, alpha, "hdxpdyp" + sfx);
    TLine *v = new TLine(0.0, amin, 0.0, amax);
    TLine *hz = new TLine(amin, 0.0, amax, 0.0);
    for (TLine *ln : { v, hz }) {
      ln->SetLineColor(kRed);
      ln->SetLineWidth(1);
      ln->SetLineStyle(2);
      ln->Draw("SAME");
    }
  }

  c->cd();
  c->Update();
  c->Print(pdfname);
  gStyle->SetPalette(kRainBow);
}


// Overall (all-modules-combined) 1D residuals:
// FT U | FT V | FPP U | FPP V, set 1 (top row) vs set 2 (bottom row).
// Mean and sigma from the per-set FitPeak results; the Gaussian is
// drawn over mean +- 2 sigma (the FitPeak refit window).
void DrawResidual1DComparisonPage(TCanvas *c, SetResults *sets[2],
                                  const char *label1, const char *label2,
                                  const char *pdfname) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  RowLayout L = SetupRowComparisonCanvas(c, label1, label2);

  for (int is = 0; is < 2; is++) {
    SetResults &r = *sets[is];
    TString sfx = TString::Format("_cmpR1D_%d", is);
    const double alpha = (is == 0) ? polcmp_set1_alpha : 1.0;

    TH1D *src[4] = { r.hist.h_eresidu_FT,  r.hist.h_eresidv_FT,
                     r.hist.h_eresidu_FPP, r.hist.h_eresidv_FPP };
    GEPFitResult *fit[4] = { &r.fit_eresidu_FT,  &r.fit_eresidv_FT,
                             &r.fit_eresidu_FPP, &r.fit_eresidv_FPP };

    for (int k = 0; k < 4; k++) {
      MakeRowPad(c, TString::Format("pres%d", k) + sfx, k, 4, L, is, false);
      TH1D *h = CleanClone(src[k], TString::Format("hres%d", k) + sfx);
      h->SetLineColor(kBlack);
      h->SetLineWidth(2);
      h->SetFillStyle(1001);
      h->SetFillColorAlpha(kAzure - 9, 0.55);   // light blue, semi-transparent
      if (alpha < 1.0) FadeHist(h, alpha);

      const GEPFitResult &f = *fit[k];
      TF1 *g = nullptr;
      if (f.valid && f.sigma > 0) {
        g = new TF1(TString::Format("gres%d", k) + sfx, "gaus",
                    f.mean - 2.0 * f.sigma, f.mean + 2.0 * f.sigma);
        g->SetParameters(f.amplitude, f.mean, f.sigma);
        g->SetLineColorAlpha(kBlue + 2, alpha);
        g->SetLineWidth(3);
        h->SetMaximum(1.2 * std::max(h->GetMaximum(), g->GetMaximum()));
      } else {
        h->SetMaximum(1.2 * h->GetMaximum());
      }
      h->Draw("hist");
      if (g) g->Draw("SAME");

      double mu = (f.valid) ? f.mean  : h->GetMean();
      double sg = (f.valid) ? f.sigma : h->GetStdDev();
      DrawValueBox({ TString::Format("Mean = %.4f mm", mu),
                     TString::Format("#sigma = %.4f mm", sg) },
                   0.56, 0.76, 0.95, 0.90);
    }
  }

  c->cd();
  c->Update();
  c->Print(pdfname);
}


// files1 / files2: space- or comma-separated files or wildcards.
// Empty files1 -> runlist / rootfile_wildcard1 from gep_config.h.
// Empty files2 -> rootfile_set2 from gep_config.h; if that is also
// empty, only set 1 is analysed.
// Set 1 uses globalcut / globalcut_thetafpp, set 2 uses
// globalcut_set2 / globalcut_thetafpp_set2 (gep_config.h).
void gep_physics(TString files1 = "", TString files2 = "",
                 TString label1 = "", TString label2 = "") {

  if (files2 == "") files2 = rootfile_set2;
  if (label1 == "") label1 = label_set1;
  if (label2 == "") label2 = label_set2;
  const bool do_compare = (files2 != "");

  const char *pdfname = "gep_physics_output.pdf";

  // ==========================================================
  // Chains
  // ==========================================================

  TChain *C = new TChain("T");

  if (files1 != "") {
    AddFilesToChain(C, files1);
  } else if(use_runlist){
    for(int i=0; i<nruns; i++){
        int nf = C->Add(Form("%sgep5_fullreplay_%d*.root", rootdir, runlist[i]));
        if(nf==0) cout << "Warning: no files found for run " << runlist[i] << endl;
    }
  } else {
    C->Add(rootfile_wildcard1);
  }

  TChain *C2 = nullptr;
  if (do_compare) {
    C2 = new TChain("T");
    AddFilesToChain(C2, files2);
  }

  // ==========================================================
  // Analyse each set with its own cuts
  // ==========================================================

  SetResults set1 = AnalyzeSet(C, label1, "", globalcut, globalcut_thetafpp);

  SetResults set2;
  if (do_compare) {
    set2 = AnalyzeSet(C2, label2, "_set2", globalcut_set2, globalcut_thetafpp_set2);
  }

  // ==========================================================
  // ROOT output (set 2 histograms carry a _set2 suffix)
  // ==========================================================

  TFile *fout = new TFile( "gep_physics_output.root", "RECREATE" );
  WriteSet(set1);
  if (do_compare) WriteSet(set2);
  fout->Close();

  // ==========================================================
  // PDF output: all set-1 pages, all set-2 pages, then the
  // side-by-side U/V residual comparison pages (module-wise,
  // layer-wise; 1920x1080, no stats).
  // ==========================================================

  TCanvas *c1 = new TCanvas( "c1", "GEp Physics", 900, 700 );
  c1->Print(Form("%s[", pdfname));

  DrawSetPages(c1, set1, pdfname);

  if (do_compare) {
    DrawSetPages(c1, set2, pdfname);

    TH2D *mod1[4] = { set1.hist.h_eresidu_FT_module,  set1.hist.h_eresidv_FT_module,
                      set1.hist.h_eresidu_FPP_module, set1.hist.h_eresidv_FPP_module };
    TH2D *mod2[4] = { set2.hist.h_eresidu_FT_module,  set2.hist.h_eresidv_FT_module,
                      set2.hist.h_eresidu_FPP_module, set2.hist.h_eresidv_FPP_module };
    DrawResidualComparisonPage(c1, mod1, mod2, label1, label2, pdfname);

    TH2D *lay1[4] = { set1.hist.h_eresidu_FT_layer,  set1.hist.h_eresidv_FT_layer,
                      set1.hist.h_eresidu_FPP_layer, set1.hist.h_eresidv_FPP_layer };
    TH2D *lay2[4] = { set2.hist.h_eresidu_FT_layer,  set2.hist.h_eresidv_FT_layer,
                      set2.hist.h_eresidu_FPP_layer, set2.hist.h_eresidv_FPP_layer };
    DrawResidualComparisonPage(c1, lay1, lay2, label1, label2, pdfname);

    // Same maps again, 2 rows (U, V) x 4 columns (FT/FPP set 1 | FT/FPP set 2)
    DrawResidualComparisonPage2x4(c1, mod1, mod2, label1, label2, pdfname);
    DrawResidualComparisonPage2x4(c1, lay1, lay2, label1, label2, pdfname);

    // Overall 1D residuals, set 1 (top row) vs set 2 (bottom row)
    SetResults *both[2] = { &set1, &set2 };
    DrawResidual1DComparisonPage(c1, both, label1, label2, pdfname);

    // Polarimeter plots, set 1 (top row) vs set 2 (bottom row)
    DrawPolarimeterComparisonRow1(c1, both, label1, label2, pdfname);
    DrawPolarimeterComparisonRow2(c1, both, label1, label2, pdfname);
  }

  c1->Print(Form("%s]", pdfname));

  std::cout
    << "Output written to:"
    << std::endl
    << "  gep_physics_output.root"
    << std::endl
    << "  " << pdfname
    << std::endl;
}
