#ifndef GEP_PLOT_UTILS_H
#define GEP_PLOT_UTILS_H

#include "TPaveText.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TExec.h"
#include "TList.h"
#include "TPaletteAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TChain.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStyle.h"
#include <iostream>
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
  g->SetLineWidth(1);

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

// Same as MakeFitStatsBox, but sourced from an already-computed
// GEPFitResult (mean/sigma/chi2/ndf) instead of a live TF1* -- used on
// pages where the fit was done earlier via FitPeak() and only the
// resulting numbers are needed for display. Same large, borderless,
// translucent-white style as the polarimeter kinematics page (page 4)
// in place of the small default ROOT stat/fit box.
TPaveText *MakeFitStatsBoxFromResult(TH1 *hist, const GEPFitResult &fit, double x1, double y1, double x2, double y2, const char *unit = "mm") {
  TPaveText *box = MakeStatsBox(hist, x1, y1, x2, y2);
  box->AddText(Form("Fit #mu = %.3g #pm %.2g %s", fit.mean, fit.mean_error, unit));
  box->AddText(Form("Fit #sigma = %.3g #pm %.2g %s", fit.sigma, fit.sigma_error, unit));
  if (fit.ndf > 0) {
    box->AddText(Form("#chi^{2}/NDF = %.2f/%d = %.2f",
                      fit.chi2, fit.ndf, fit.chi2 / fit.ndf));
  }
  return box;
}

// Add every file / wildcard in a space- or comma-separated list to a
// chain. Returns the number of files added.
int AddFilesToChain(TChain *C, const char *filelist) {
  int ntot = 0;
  TObjArray *tok = TString(filelist).Tokenize(" ,");
  for (int i = 0; i < tok->GetEntries(); i++) {
    TString pattern = ((TObjString *)tok->At(i))->GetString();
    int nf = C->Add(pattern);
    if (nf == 0) std::cout << "Warning: no files found for " << pattern << std::endl;
    ntot += nf;
  }
  delete tok;
  return ntot;
}

// One 16:9 page (1920x1080, for a slide), 4 rows x 2 columns:
// rows = FT U, FT V, FPP U, FPP V; left column = set 1, right column
// = set 2. Pads are placed by hand with tight margins; the set label
// is written once above each column and a gap + vertical divider
// separates the two sets. Plain COLZ maps on a common 0-1
// (per-column-normalized) z scale -- no stats, no fits. Clones are
// drawn so the original histograms keep their titles/styles.
void DrawResidualComparisonPage(TCanvas *c, TH2D *set1[4], TH2D *set2[4],
                                const char *label1, const char *label2,
                                const char *pdfname) {
  const char *rowname[4] = { "FT U-plane", "FT V-plane", "FPP U-plane", "FPP V-plane" };

  // Layout in canvas NDC
  const double header_h = 0.045;               // strip for column labels
  const double gap      = 0.016;               // separation between sets
  const double col_w    = (1.0 - gap) / 2.0;
  const double row_h    = (1.0 - header_h) / 4.0;
  const double colx[2]  = { 0.0, col_w + gap };

  // Pad margins (fractions of each pad)
  const double lm = 0.050, rm = 0.060, tm = 0.035, bm = 0.155;

  c->Clear();
  c->SetCanvasSize(1920, 1080);
  c->SetFillColor(kWhite);
  c->cd();

  // Column headers
  TLatex header;
  header.SetNDC();
  header.SetTextFont(62);
  header.SetTextSize(0.030);
  header.SetTextAlign(22);
  header.DrawLatex(colx[0] + 0.5 * col_w, 1.0 - 0.5 * header_h, label1);
  header.DrawLatex(colx[1] + 0.5 * col_w, 1.0 - 0.5 * header_h, label2);

  // Divider between the two sets
  TLine *divider = new TLine(0.5, 0.0, 0.5, 1.0);
  divider->SetNDC();
  divider->SetLineColor(kGray + 2);
  divider->SetLineWidth(3);
  divider->Draw();

  for (int irow = 0; irow < 4; irow++) {
    for (int icol = 0; icol < 2; icol++) {
      TH2D *src = (icol == 0) ? set1[irow] : set2[irow];
      if (!src) continue;

      double y2 = 1.0 - header_h - irow * row_h;
      double y1 = y2 - row_h;

      c->cd();
      TPad *pad = new TPad(Form("pcmp_%s_%d%d", src->GetName(), irow, icol), "",
                           colx[icol], y1, colx[icol] + col_w, y2);
      pad->SetLeftMargin(lm);
      pad->SetRightMargin(rm);
      pad->SetTopMargin(tm);
      pad->SetBottomMargin(bm);
      pad->SetFrameLineWidth(1);
      pad->Draw();
      pad->cd();

      // Per-pad palette: set 1 in kRainBow at 75% opacity (lighter,
      // over the white pad), set 2 in the usual opaque kRainBow. A TExec in each pad
      // switches the global palette right before that pad's histogram
      // (and its colour bar) is painted, so both coexist on one page.
      TExec *pal_exec = new TExec(Form("palexec_%d%d", irow, icol),
                                  (icol == 0)
                                    ? "gStyle->SetPalette(kRainBow, 0, 0.75);"
                                    : "gStyle->SetPalette(kRainBow);");

      TH2D *h = (TH2D *)src->Clone(Form("%s_cmp_%d%d", src->GetName(), irow, icol));
      h->SetTitle("");
      h->SetStats(0);
      h->SetMinimum(0);
      h->SetMaximum(1);

      h->GetYaxis()->SetTitle("");
      h->GetYaxis()->SetNdivisions(505);
      h->GetYaxis()->SetLabelSize(0.085);
      h->GetYaxis()->SetLabelOffset(0.006);

      h->GetXaxis()->SetTitleSize(0.085);
      h->GetXaxis()->SetTitleOffset(0.80);
      h->GetXaxis()->SetLabelSize(0.085);
      h->GetXaxis()->SetLabelOffset(0.008);
      h->GetXaxis()->SetTickLength(0.04);

      // Axes-only first draw (sets up the frame; TH1::Draw without
      // "same" clears the pad), then the TExec, then the colour map in
      // that palette, then the axes again on top of the cells. The
      // first pass must not paint the cells, otherwise an opaque copy
      // sits underneath and hides any transparency.
      h->SetContour(99);
      h->Draw("AXIS");
      pal_exec->Draw();
      h->Draw("COLZ SAME");
      h->Draw("AXIS SAME");
      pad->Update();

      // Narrow palette tucked into the right margin
      TPaletteAxis *pal = (TPaletteAxis *)h->GetListOfFunctions()->FindObject("palette");
      if (pal) {
        pal->SetX1NDC(1.0 - rm + 0.006);
        pal->SetX2NDC(1.0 - rm + 0.020);
        pal->SetY1NDC(bm);
        pal->SetY2NDC(1.0 - tm);
        pal->SetLabelSize(0.07);
        pal->GetAxis()->SetNdivisions(2);
      }

      // Row label inside the plot (top-left)
      TPaveText *tag = new TPaveText(lm + 0.008, 1.0 - tm - 0.20, lm + 0.34, 1.0 - tm - 0.03, "NDC");
      tag->SetFillColorAlpha(kWhite, 0.85);
      tag->SetBorderSize(0);
      tag->SetTextFont(42);
      tag->SetTextSize(0.10);
      tag->SetTextAlign(12);
      tag->AddText(Form("%s residual (mm)", rowname[irow]));
      tag->Draw();

      pad->Modified();
    }
  }

  c->cd();
  c->Update();
  c->Print(pdfname);

  // Leave the global palette as the rest of the macro expects
  gStyle->SetPalette(kRainBow);
}

#endif
