#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"

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

  C->Add(rootfile);

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

  //---------- Module-wise residuals: FT u/v and FPP u/v, all on one page

  c1->Clear();
  c1->Divide(2,2);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kMint);

  c1->cd(1);
  hist.h_eresidu_FT_module->SetMinimum(0);
  hist.h_eresidu_FT_module->SetMaximum(1);
  hist.h_eresidu_FT_module->Draw("COLZ");
  g_eresidu_FT->Draw("P SAME");

  c1->cd(2);
  hist.h_eresidv_FT_module->SetMinimum(0);
  hist.h_eresidv_FT_module->SetMaximum(1);
  hist.h_eresidv_FT_module->Draw("COLZ");
  g_eresidv_FT->Draw("P SAME");

  c1->cd(3);
  hist.h_eresidu_FPP_module->SetMinimum(0);
  hist.h_eresidu_FPP_module->SetMaximum(1);
  hist.h_eresidu_FPP_module->Draw("COLZ");
  g_eresidu_FPP->Draw("P SAME");

  c1->cd(4);
  hist.h_eresidv_FPP_module->SetMinimum(0);
  hist.h_eresidv_FPP_module->SetMaximum(1);
  hist.h_eresidv_FPP_module->Draw("COLZ");
  g_eresidv_FPP->Draw("P SAME");

  c1->Print("gep_physics_output.pdf");

  gStyle->SetOptStat(1111);

  c1->Print("gep_physics_output.pdf]");

  std::cout
    << "Output written to:"
    << std::endl
    << "  gep_physics_output.root"
    << std::endl
    << "  gep_physics_output.pdf"
    << std::endl;
}
