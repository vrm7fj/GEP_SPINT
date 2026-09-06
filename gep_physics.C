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
 
  c1->Print("gep_physics_output.pdf]");

  std::cout
    << "Output written to:"
    << std::endl
    << "  gep_physics_output.root"
    << std::endl
    << "  gep_physics_output.pdf"
    << std::endl;
}
