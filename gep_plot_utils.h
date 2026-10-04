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
#include <algorithm>
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

// ------------------------------------------------------------
// Shared pieces for the 16:9 (1920x1080) set-1 vs set-2 residual
// comparison pages. Plain COL maps (no colour bar) on the common 0-1
// (per-column-normalized) z scale -- no stats, no fits. Set 1 is drawn
// in kRainBow at 80% opacity, set 2 in opaque kRainBow. Clones are
// drawn so the original histograms keep their titles/styles.
// ------------------------------------------------------------

struct ResidPadStyle {
  double lm, rm, tm, bm;        // pad margins
  double label_size;            // axis label size
  double xtitle_size, xtitle_offset;
  double xtick, ytick;          // tick lengths
  double tag_x2, tag_h;         // row-label box right edge / height (pad NDC)
  double tag_text;              // row-label text size
};

void DrawResidualPad(TCanvas *c, TH2D *src, TString tagtext, bool is_set1,
                     double x1, double y1, double x2, double y2,
                     const ResidPadStyle &st, TString uid) {
  if (!src) return;

  // TPad rejects edges outside [0,1]; guard against round-off
  auto clamp01 = [](double v) { return std::min(1.0, std::max(0.0, v)); };
  x1 = clamp01(x1); y1 = clamp01(y1); x2 = clamp01(x2); y2 = clamp01(y2);

  c->cd();
  TPad *pad = new TPad("pcmp_" + uid, "", x1, y1, x2, y2);
  pad->SetLeftMargin(st.lm);
  pad->SetRightMargin(st.rm);
  pad->SetTopMargin(st.tm);
  pad->SetBottomMargin(st.bm);
  pad->SetFrameLineWidth(1);
  pad->Draw();
  pad->cd();

  // Per-pad palette via TExec, so both palettes coexist on one page
  TExec *pal_exec = new TExec("palexec_" + uid,
                              is_set1 ? "gStyle->SetPalette(kRainBow, 0, 0.80);"
                                      : "gStyle->SetPalette(kRainBow);");

  TH2D *h = (TH2D *)src->Clone(TString(src->GetName()) + "_cmp_" + uid);
  h->SetTitle("");
  h->SetStats(0);
  h->SetMinimum(0);
  h->SetMaximum(1);

  h->GetYaxis()->SetTitle("");
  h->GetYaxis()->SetNdivisions(505);
  h->GetYaxis()->SetLabelSize(st.label_size);
  h->GetYaxis()->SetLabelOffset(0.006);
  h->GetYaxis()->SetTickLength(st.ytick);

  h->GetXaxis()->SetTitleSize(st.xtitle_size);
  h->GetXaxis()->SetTitleOffset(st.xtitle_offset);
  h->GetXaxis()->SetLabelSize(st.label_size);
  h->GetXaxis()->SetLabelOffset(0.008);
  h->GetXaxis()->SetTickLength(st.xtick);

  // Axes-only first draw (sets up the frame; TH1::Draw without "same"
  // clears the pad), then the TExec, then the colour map in that
  // palette, then the axes again on top of the cells. The first pass
  // must not paint the cells, otherwise an opaque copy sits underneath
  // and hides the transparency.
  h->SetContour(99);
  h->Draw("AXIS");
  pal_exec->Draw();
  h->Draw("COL SAME");
  h->Draw("AXIS SAME");
  pad->Update();

  // Row label inside the plot (top-left)
  TPaveText *tag = new TPaveText(st.lm + 0.008, 1.0 - st.tm - 0.03 - st.tag_h,
                                 st.tag_x2, 1.0 - st.tm - 0.03, "NDC");
  tag->SetFillColorAlpha(kWhite, 0.85);
  tag->SetBorderSize(0);
  tag->SetTextFont(42);
  tag->SetTextSize(st.tag_text);
  tag->SetTextAlign(12);
  tag->AddText(tagtext);
  tag->Draw();

  pad->Modified();
}

// Canvas setup, set labels above each half and a divider between them
void SetupComparisonCanvas(TCanvas *c, double header_h,
                           double xcen1, double xcen2,
                           const char *label1, const char *label2) {
  c->Clear();
  c->SetCanvasSize(1920, 1080);
  c->SetFillColor(kWhite);
  c->cd();

  TLatex header;
  header.SetNDC();
  header.SetTextFont(62);
  header.SetTextSize(0.030);
  header.SetTextAlign(22);
  header.DrawLatex(xcen1, 1.0 - 0.5 * header_h, label1);
  header.DrawLatex(xcen2, 1.0 - 0.5 * header_h, label2);

  TLine *divider = new TLine(0.5, 0.0, 0.5, 1.0);
  divider->SetNDC();
  divider->SetLineColor(kGray + 2);
  divider->SetLineWidth(3);
  divider->Draw();
}

void FinishComparisonPage(TCanvas *c, const char *pdfname) {
  c->cd();
  c->Update();
  c->Print(pdfname);
  gStyle->SetPalette(kRainBow);   // leave the global palette as the macro expects
}

// Arrays set1[4] / set2[4] are ordered { FT U, FT V, FPP U, FPP V }.

// Page layout A: 4 rows x 2 columns.
//   rows = FT U, FT V, FPP U, FPP V; left column = set 1, right = set 2.
void DrawResidualComparisonPage(TCanvas *c, TH2D *set1[4], TH2D *set2[4],
                                const char *label1, const char *label2,
                                const char *pdfname) {
  const char *rowname[4] = { "FT U-plane", "FT V-plane", "FPP U-plane", "FPP V-plane" };

  const double header_h = 0.045;
  const double gap      = 0.016;
  const double col_w    = (1.0 - gap) / 2.0;
  const double row_h    = (1.0 - header_h) / 4.0;
  const double colx[2]  = { 0.0, col_w + gap };

  // Wide, short pads (~950 x 260 px)
  ResidPadStyle st = { 0.050, 0.015, 0.035, 0.155,
                       0.085, 0.085, 0.80, 0.04, 0.03,
                       0.39, 0.17, 0.10 };

  SetupComparisonCanvas(c, header_h, colx[0] + 0.5 * col_w, colx[1] + 0.5 * col_w, label1, label2);

  for (int irow = 0; irow < 4; irow++) {
    for (int icol = 0; icol < 2; icol++) {
      TH2D *src = (icol == 0) ? set1[irow] : set2[irow];
      double y2 = 1.0 - header_h - irow * row_h;
      TString tagtext = TString::Format("%s residual (mm)", rowname[irow]);
      TString uid     = TString::Format("A%s_%d%d", src ? src->GetName() : "x", irow, icol);
      DrawResidualPad(c, src, tagtext, icol == 0,
                      colx[icol], y2 - row_h, colx[icol] + col_w, y2, st, uid);
    }
  }

  FinishComparisonPage(c, pdfname);
}

// Page layout B: 2 rows x 4 columns.
//   columns = FT set 1, FPP set 1 | FT set 2, FPP set 2
//   rows    = U-plane, V-plane
void DrawResidualComparisonPage2x4(TCanvas *c, TH2D *set1[4], TH2D *set2[4],
                                   const char *label1, const char *label2,
                                   const char *pdfname) {
  const char *detname[2]   = { "FT", "FPP" };
  const char *planename[2] = { "U-plane", "V-plane" };

  const double header_h = 0.050;
  const double gap      = 0.016;
  const double col_w    = (1.0 - gap) / 4.0;
  const double row_h    = (1.0 - header_h) / 2.0;
  const double colx[4]  = { 0.0, col_w, 2.0 * col_w + gap, 3.0 * col_w + gap };

  // Near-square pads (~470 x 510 px)
  ResidPadStyle st = { 0.110, 0.020, 0.015, 0.100,
                       0.048, 0.050, 0.95, 0.02, 0.02,
                       0.82, 0.075, 0.048 };

  SetupComparisonCanvas(c, header_h, colx[0] + col_w, colx[2] + col_w, label1, label2);

  for (int iplane = 0; iplane < 2; iplane++) {
    for (int icol = 0; icol < 4; icol++) {
      int iset = icol / 2;          // 0 = set 1, 1 = set 2
      int idet = icol % 2;          // 0 = FT, 1 = FPP
      TH2D *src = (iset == 0 ? set1 : set2)[2 * idet + iplane];
      double y2 = 1.0 - header_h - iplane * row_h;
      TString tagtext = TString::Format("%s %s residual (mm)", detname[idet], planename[iplane]);
      TString uid     = TString::Format("B%s_%d%d", src ? src->GetName() : "x", iplane, icol);
      DrawResidualPad(c, src, tagtext, iset == 0,
                      colx[icol], y2 - row_h, colx[icol] + col_w, y2, st, uid);
    }
  }

  FinishComparisonPage(c, pdfname);
}

#endif
