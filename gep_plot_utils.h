#ifndef GEP_PLOT_UTILS_H
#define GEP_PLOT_UTILS_H

#include "TPaveText.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <vector>
#include <string>

#include "gep_fit.h"

void DrawTextBox(const std::vector<std::string> &text, double x1, double y1, double x2, double y2) {
  TPaveText *box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillStyle(0);
  box->SetBorderSize(0);
  box->SetTextAlign(12);
  box->SetTextSize(0.025);
  box->SetTextFont(42);
  for (const auto &line : text) box->AddText(line.c_str());
  box->Draw();
}

// Normalize each column (fixed module, i.e. x-bin) of a module-wise 2D
// residual histogram to its own peak bin content. This makes every
// module's residual shape visible on the same 0-1 color scale,
// regardless of how many hits that module happened to get. Modules
// with no entries (peak == 0) are simply left as-is (all zero).
void NormalizeModuleColumns(TH2D *h) {
  if (!h) return;

  int nbinsX = h->GetNbinsX();
  int nbinsY = h->GetNbinsY();

  for (int ix = 1; ix <= nbinsX; ix++) {

    TH1D *proj = h->ProjectionY("_py_normtmp", ix, ix);
    double maxVal = proj->GetMaximum();
    delete proj;

    if (maxVal <= 0) continue;

    for (int iy = 1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix, iy);
      h->SetBinContent(ix, iy, val / maxVal);
    }
  }
}

// Turn a vector of per-module fit results (from FitModuleColumns) into a
// TGraphErrors of mean +/- sigma vs module index, ready to overlay on the
// corresponding 2D histogram with Draw("P SAME"). Modules that were never
// fit (empty, .valid == false) are skipped, not plotted as zero.
TGraphErrors *BuildResidualGraph(const std::vector<GEPFitResult> &fits) {

  TGraphErrors *g = new TGraphErrors();

  int ipoint = 0;
  for (size_t imod = 0; imod < fits.size(); imod++) {

    if (!fits[imod].valid) continue;

    g->SetPoint(ipoint, (double)imod, fits[imod].mean);
    g->SetPointError(ipoint, 0, fits[imod].sigma);
    ipoint++;
  }

  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);

  return g;
}

// White-translucent, borderless stats box (N / Mean / RMS), ported from
// polarimeter_recon.C -- used for the polarimeter kinematics page, which
// keeps that script's own OptStat(0)/manual-TPaveText convention rather
// than this file's usual built-in OptStat/OptFit boxes.
TPaveText *MakeStatsBox(TH1 *hist, double x1, double y1, double x2, double y2) {
  TPaveText *box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillColorAlpha(kWhite, 0.88);
  box->SetBorderSize(0);
  box->SetLineColor(kBlack);
  box->SetLineWidth(1);
  box->SetTextColor(kBlack);
  box->SetTextFont(42);
  box->SetTextSize(0.05);
  box->SetTextAlign(12);
  box->SetMargin(0.08);
  box->AddText(Form("N = %.0f", hist->GetEntries()));
  box->AddText(Form("Mean = %.3g", hist->GetMean()));
  box->AddText(Form("RMS = %.3g", hist->GetRMS()));
  return box;
}

// Same, plus fit mu/sigma/chi2-ndf -- for a fitted 1D histogram.
TPaveText *MakeFitStatsBox(TH1 *hist, TF1 *fit, double x1, double y1, double x2, double y2) {
  TPaveText *box = MakeStatsBox(hist, x1, y1, x2, y2);
  box->AddText(Form("Fit #mu = %.3g #pm %.2g deg", fit->GetParameter(1), fit->GetParError(1)));
  box->AddText(Form("Fit #sigma = %.3g #pm %.2g deg", fit->GetParameter(2), fit->GetParError(2)));
  if (fit->GetNDF() > 0) {
    box->AddText(Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                      fit->GetChisquare(), fit->GetNDF(),
                      fit->GetChisquare() / fit->GetNDF()));
  }
  return box;
}

#endif
