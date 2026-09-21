#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TMath.h"

#include "gep_config.h"
#include "gep_fill_vectors.h"
#include "gep_histograms.h"
#include "gep_fit.h"
#include "gep_plot_utils.h"


void gep_physics() {

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  // ==========================================================
  // Create chain
  // ==========================================================

  TChain *C = new TChain("T");

  if(use_runlist){
    for(int i=0; i<nruns; i++){
        int nf = C->Add(Form("%sgep5_fullreplay_%d*.root", rootdir, runlist[i]));
        if(nf==0) cout << "Warning: no files found for run " << runlist[i] << endl;
    }
  } else {
    C->Add(rootfile_wildcard1);
    C->Add(rootfile_wildcard2);
  }

  // ==========================================================
  // Fill vectors
  // ==========================================================

  GEPData data;

  FillVectors(C, data);

  // ==========================================================
  // Create and fill histograms
  // ==========================================================

  GEPHistograms hist = CreateHistograms();

  FillHistograms(data, hist);

  // ==========================================================
  // Fill polarimeter reconstruction vectors/histograms
  // (ported from polarimeter_recon.C -- separate cut, separate
  // FPP-besttrack-indexed pass over the same chain)
  // ==========================================================

  PolarimeterData poldata;

  FillPolarimeterVectors(C, poldata);

  PolarimeterHistograms polhist = CreatePolarimeterHistograms();

  FillPolarimeterHistograms(poldata, polhist);

  // ==========================================================
  // Fit histograms
  // ==========================================================

  GEPFitResult fit_dx  = FitPeak(hist.h_dx);
  GEPFitResult fit_dy  = FitPeak(hist.h_dy);
  GEPFitResult fit_dxp = FitPeak(hist.h_dxp);
  GEPFitResult fit_dyp = FitPeak(hist.h_dyp);
  GEPFitResult fit_vz = FitPeak(hist.h_vz);

  // ==========================================================
  // Fit each module's residual distribution (mean +/- sigma).
  // Done on the raw (pre-normalization) 2D histograms so the fits
  // see real counts. Modules with no hits come back with
  // .valid == false and are skipped automatically -- nothing fails.
  // ==========================================================

  std::vector<GEPFitResult> fits_eresidu_FT  = FitModuleColumns(hist.h_eresidu_FT_module,  nmod_ft);
  std::vector<GEPFitResult> fits_eresidv_FT  = FitModuleColumns(hist.h_eresidv_FT_module,  nmod_ft);
  std::vector<GEPFitResult> fits_eresidu_FPP = FitModuleColumns(hist.h_eresidu_FPP_module, nmod_fpp);
  std::vector<GEPFitResult> fits_eresidv_FPP = FitModuleColumns(hist.h_eresidv_FPP_module, nmod_fpp);

  // Fit the overall (all-modules-combined) U/V residual distributions
  GEPFitResult fit_eresidu_FT  = FitPeak(hist.h_eresidu_FT);
  GEPFitResult fit_eresidv_FT  = FitPeak(hist.h_eresidv_FT);
  GEPFitResult fit_eresidu_FPP = FitPeak(hist.h_eresidu_FPP);
  GEPFitResult fit_eresidv_FPP = FitPeak(hist.h_eresidv_FPP);

  TGraphErrors *g_eresidu_FT  = BuildResidualGraph(fits_eresidu_FT);
  TGraphErrors *g_eresidv_FT  = BuildResidualGraph(fits_eresidv_FT);
  TGraphErrors *g_eresidu_FPP = BuildResidualGraph(fits_eresidu_FPP);
  TGraphErrors *g_eresidv_FPP = BuildResidualGraph(fits_eresidv_FPP);

  // ==========================================================
  // Normalize module-wise 2D residual histograms for display
  // (each module column scaled to its own peak bin). Purely
  // cosmetic -- done after fitting so it doesn't affect the fits.
  // ==========================================================

  NormalizeModuleColumns(hist.h_eresidu_FT_module);
  NormalizeModuleColumns(hist.h_eresidv_FT_module);
  NormalizeModuleColumns(hist.h_eresidu_FPP_module);
  NormalizeModuleColumns(hist.h_eresidv_FPP_module);

  // ==========================================================
  // ROOT output
  // ==========================================================

  TFile *fout = new TFile( "gep_physics_output.root", "RECREATE" );

  hist.h_dx->Write();
  hist.h_dy->Write();
  hist.h_dxp->Write();
  hist.h_dyp->Write();
  hist.h_vz->Write();
  hist.h_vx->Write();
  hist.h_vy->Write();
  hist.h_vxvy->Write();
  hist.h_eresidu_FT_module->Write();
  hist.h_eresidv_FT_module->Write();
  hist.h_eresidu_FPP_module->Write();
  hist.h_eresidv_FPP_module->Write();
  hist.h_eresidu_FT->Write();
  hist.h_eresidv_FT->Write();
  hist.h_eresidu_FPP->Write();
  hist.h_eresidv_FPP->Write();

  fout->Close();

  // ==========================================================
  // PDF output
  // ==========================================================

  //---------- Check Alignment

  TCanvas *c1 = new TCanvas( "c1", "GEp Physics", 900, 700 );

  c1->Divide(2,2);

  c1->Print("gep_physics_output.pdf[");

  c1->cd(1);
  hist.h_dx->Draw();

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
  hist.h_dy->Draw();

  c1->cd(3);
  hist.h_dxp->Draw();

  c1->cd(4);
  hist.h_dyp->Draw();

  c1->Print("gep_physics_output.pdf");

  //---------- Check Target

  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1);
  hist.h_vz->Draw();

  DrawTextBox({
  "Cuts:",
  "FT: N_{hits} > 4 || N_{goodhits} > 2",
  "FPP: N_{hits} > 4 || N_{goodhits} > 2",
  "FT: #chi^{2}/ndf < 200",
  "FPP: #chi^{2}/ndf < 200",
  "HCAL: N_{blk} > 1"
  }, 0.12, 0.60, 0.42, 0.88);

  c1->cd(2);
  hist.h_vx->Draw();

  c1->cd(3);
  hist.h_vy->Draw();

  c1->cd(4);
  hist.h_vxvy->Draw();

  c1->Print("gep_physics_output.pdf");

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
  hist.h_eresidu_FT_module->SetMinimum(0);
  hist.h_eresidu_FT_module->SetMaximum(1);
  hist.h_eresidu_FT_module->SetStats(0);
  hist.h_eresidu_FT_module->Draw("COLZ");
  g_eresidu_FT->Draw("P SAME");

  c1->cd(2);
  hist.h_eresidv_FT_module->SetMinimum(0);
  hist.h_eresidv_FT_module->SetMaximum(1);
  hist.h_eresidv_FT_module->SetStats(0);
  hist.h_eresidv_FT_module->Draw("COLZ");
  g_eresidv_FT->Draw("P SAME");

  c1->cd(3);
  hist.h_eresidu_FT->SetStats(0);
  hist.h_eresidu_FT->Draw();
  if (fit_eresidu_FT.valid) {
    MakeFitStatsBoxFromResult(hist.h_eresidu_FT, fit_eresidu_FT, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(4);
  hist.h_eresidv_FT->SetStats(0);
  hist.h_eresidv_FT->Draw();
  if (fit_eresidv_FT.valid) {
    MakeFitStatsBoxFromResult(hist.h_eresidv_FT, fit_eresidv_FT, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(5);
  hist.h_eresidu_FPP_module->SetMinimum(0);
  hist.h_eresidu_FPP_module->SetMaximum(1);
  hist.h_eresidu_FPP_module->SetStats(0);
  hist.h_eresidu_FPP_module->Draw("COLZ");
  g_eresidu_FPP->Draw("P SAME");

  c1->cd(6);
  hist.h_eresidv_FPP_module->SetMinimum(0);
  hist.h_eresidv_FPP_module->SetMaximum(1);
  hist.h_eresidv_FPP_module->SetStats(0);
  hist.h_eresidv_FPP_module->Draw("COLZ");
  g_eresidv_FPP->Draw("P SAME");

  c1->cd(7);
  hist.h_eresidu_FPP->SetStats(0);
  hist.h_eresidu_FPP->Draw();
  if (fit_eresidu_FPP.valid) {
    MakeFitStatsBoxFromResult(hist.h_eresidu_FPP, fit_eresidu_FPP, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->cd(8);
  hist.h_eresidv_FPP->SetStats(0);
  hist.h_eresidv_FPP->Draw();
  if (fit_eresidv_FPP.valid) {
    MakeFitStatsBoxFromResult(hist.h_eresidv_FPP, fit_eresidv_FPP, 0.55, 0.60, 0.94, 0.90)->Draw();
  }

  c1->Print("gep_physics_output.pdf");

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

  TLine *lfpp_theta1 = new TLine(fpp_theta_min, 0, fpp_theta_min, polhist.h_theta_fpp->GetMaximum());
  TLine *lfpp_theta2 = new TLine(fpp_theta_max, 0, fpp_theta_max, polhist.h_theta_fpp->GetMaximum());
  lfpp_theta1->SetLineWidth(1);
  lfpp_theta2->SetLineWidth(1);
  lfpp_theta1->SetLineColor(kRed+1);
  lfpp_theta2->SetLineColor(kRed+1);

  polhist.h_theta_fpp->Draw("hist");
  lfpp_theta1->Draw("SAME");
  lfpp_theta2->Draw("SAME");

  MakeStatsBox(polhist.h_theta_fpp, 0.62, 0.73, 0.91, 0.90)->Draw();

  c1->cd(2);
  gPad->SetLogy(0);
  polhist.h_doca->Draw("E1 P");

  double doca_fit_min = 0.0;
  double doca_fit_max = 0.2;
  TF1 *f_halfgaus = new TF1("f_halfgaus", "[0]*exp(-0.5*x*x/([1]*[1]))", doca_fit_min, doca_fit_max);
  f_halfgaus->SetParNames("A", "#sigma");
  f_halfgaus->SetParameters(polhist.h_doca->GetMaximum(), 0.10);
  f_halfgaus->SetParLimits(0, 1e-6, 1e9);
  f_halfgaus->SetParLimits(1, 1e-5, 5.0);
  f_halfgaus->SetLineColor(kBlue+2);
  f_halfgaus->SetLineWidth(1);
  polhist.h_doca->Fit(f_halfgaus, "RQ0");

  double doca_ymax_hist = polhist.h_doca->GetMaximum();
  double doca_ymax_fit  = f_halfgaus->GetMaximum(doca_fit_min, doca_fit_max);
  polhist.h_doca->SetMaximum(1.2 * (doca_ymax_fit > doca_ymax_hist ? doca_ymax_fit : doca_ymax_hist));

  polhist.h_doca->Draw("E1 P");
  f_halfgaus->Draw("SAME");

  TLine *lfpp_sclose = new TLine(0.5, 0, 0.5, polhist.h_doca->GetMaximum());
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
  pt_doca->AddText(Form("N = %.0f", polhist.h_doca->GetEntries()));
  pt_doca->Draw("SAME");

  c1->cd(3);
  polhist.h_zclose_all->Draw("hist");
  polhist.h_zclose_sAng->Draw("hist same");
  polhist.h_zclose_lAng->Draw("hist same");

  TLine *lmin_zclose = new TLine(fpp_zclose_mean - fpp_zclose_sigma, 0,
                                 fpp_zclose_mean - fpp_zclose_sigma, polhist.h_zclose_all->GetMaximum());
  TLine *lmax_zclose = new TLine(fpp_zclose_mean + fpp_zclose_sigma, 0,
                                 fpp_zclose_mean + fpp_zclose_sigma, polhist.h_zclose_all->GetMaximum());
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
  leg_zclose->AddEntry(polhist.h_zclose_all, "all", "l");
  leg_zclose->AddEntry(polhist.h_zclose_lAng, Form("#theta_{FPP} > %.2f", fpp_theta_min), "lf");
  leg_zclose->AddEntry(polhist.h_zclose_sAng, Form("#theta_{FPP} <= %.2f", fpp_theta_min), "lf");
  leg_zclose->Draw();

  MakeStatsBox(polhist.h_zclose_all, 0.59, 0.47, 0.93, 0.65)->Draw();

  c1->cd(4);
  gPad->SetRightMargin(0.16);
  polhist.h_theta_vs_zclose->SetStats(0);
  polhist.h_theta_vs_zclose->Draw("COLZ");

  //------------------------------------------------------

  c1->cd(5);
  polhist.h_dxp->GetXaxis()->SetRangeUser(-2.0, 4.0);
  polhist.h_dxp->Draw("hist");

  const double peak_fit_half_width = 0.2; // degrees
  double dxp_fit_mean = polhist.h_dxp->GetXaxis()->GetBinCenter(polhist.h_dxp->GetMaximumBin());
  double dxp_fit_rms  = polhist.h_dxp->GetRMS();
  TF1 *fgaus_dxp = new TF1("fgaus_dxp", "gaus",
                           dxp_fit_mean - peak_fit_half_width,
                           dxp_fit_mean + 1.2*peak_fit_half_width);
  fgaus_dxp->SetParameters(polhist.h_dxp->GetMaximum(), dxp_fit_mean, 0.4);
  fgaus_dxp->SetLineColor(kRed+2);
  fgaus_dxp->SetLineWidth(1);

  bool dxp_fit_ok = polhist.h_dxp->GetEntries() > 3 && dxp_fit_rms > 0.0;
  if (dxp_fit_ok) {
    polhist.h_dxp->Fit(fgaus_dxp, "RQ0");
    fgaus_dxp->Draw("SAME");
  }

  TLine *ldxp = new TLine(0, 0, 0, polhist.h_dxp->GetMaximum());
  ldxp->SetLineWidth(1);
  ldxp->SetLineColor(kRed);
  ldxp->Draw("same");

  if (dxp_fit_ok) {
    MakeFitStatsBox(polhist.h_dxp, fgaus_dxp, 0.53, 0.58, 0.94, 0.91)->Draw();
  } else {
    MakeStatsBox(polhist.h_dxp, 0.62, 0.73, 0.91, 0.90)->Draw();
  }

  //------------------------------------------------------

  c1->cd(6);
  polhist.h_dyp->GetXaxis()->SetRangeUser(-2.0, 4.0);
  polhist.h_dyp->Draw("hist");

  double dyp_fit_mean = polhist.h_dyp->GetXaxis()->GetBinCenter(polhist.h_dyp->GetMaximumBin());
  double dyp_fit_rms  = polhist.h_dyp->GetRMS();
  TF1 *fgaus_dyp = new TF1("fgaus_dyp", "gaus",
                           dyp_fit_mean - 1.2*peak_fit_half_width,
                           dyp_fit_mean + peak_fit_half_width);
  fgaus_dyp->SetParameters(polhist.h_dyp->GetMaximum(), dyp_fit_mean, 0.4);
  fgaus_dyp->SetLineColor(kBlue+2);
  fgaus_dyp->SetLineWidth(1);

  bool dyp_fit_ok = polhist.h_dyp->GetEntries() > 3 && dyp_fit_rms > 0.0;
  if (dyp_fit_ok) {
    polhist.h_dyp->Fit(fgaus_dyp, "RQ0");
    fgaus_dyp->Draw("SAME");
  }

  TLine *ldyp = new TLine(0, 0, 0, polhist.h_dyp->GetMaximum());
  ldyp->SetLineWidth(1);
  ldyp->SetLineColor(kRed);
  ldyp->Draw("same");

  if (dyp_fit_ok) {
    MakeFitStatsBox(polhist.h_dyp, fgaus_dyp, 0.53, 0.58, 0.94, 0.91)->Draw();
  } else {
    MakeStatsBox(polhist.h_dyp, 0.62, 0.73, 0.91, 0.90)->Draw();
  }

  //------------------------------------------------------

  c1->cd(7);
  gPad->SetRightMargin(0.16);
  polhist.h_dxpdyp->SetStats(0);
  polhist.h_dxpdyp->Draw("COLZ");

  {
    double xmin = polhist.h_dxpdyp->GetXaxis()->GetXmin();
    double xmax = polhist.h_dxpdyp->GetXaxis()->GetXmax();
    double ymin = polhist.h_dxpdyp->GetYaxis()->GetXmin();
    double ymax = polhist.h_dxpdyp->GetYaxis()->GetXmax();

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

  c1->Print("gep_physics_output.pdf");

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
  polhist.h_xxp_ft->SetStats(0);
  polhist.h_xxp_ft->Draw("COLZ");

  c1->cd(2);
  polhist.h_xyp_ft->SetStats(0);
  polhist.h_xyp_ft->Draw("COLZ");

  c1->cd(3);
  polhist.h_yxp_ft->SetStats(0);
  polhist.h_yxp_ft->Draw("COLZ");

  c1->cd(4);
  polhist.h_yyp_ft->SetStats(0);
  polhist.h_yyp_ft->Draw("COLZ");

  c1->cd(5);
  polhist.h_xxp_fpp->SetStats(0);
  polhist.h_xxp_fpp->Draw("COLZ");

  c1->cd(6);
  polhist.h_xyp_fpp->SetStats(0);
  polhist.h_xyp_fpp->Draw("COLZ");

  c1->cd(7);
  polhist.h_yxp_fpp->SetStats(0);
  polhist.h_yxp_fpp->Draw("COLZ");

  c1->cd(8);
  polhist.h_yyp_fpp->SetStats(0);
  polhist.h_yyp_fpp->Draw("COLZ");

  c1->Print("gep_physics_output.pdf");

  //---------- FT/FPP chi2/ndf

  c1->Clear();
  c1->SetCanvasSize(2000, 900);
  c1->Divide(2, 1);

  gStyle->SetOptStat(1111);

  c1->cd(1);
  polhist.h_chi2_ft->SetStats(1);
  polhist.h_chi2_ft->Draw();

  c1->cd(2);
  polhist.h_chi2_fpp->SetStats(1);
  polhist.h_chi2_fpp->Draw();

  c1->Print("gep_physics_output.pdf");

  c1->Print("gep_physics_output.pdf]");

  std::cout
    << "Output written to:"
    << std::endl
    << "  gep_physics_output.root"
    << std::endl
    << "  gep_physics_output.pdf"
    << std::endl;
}
