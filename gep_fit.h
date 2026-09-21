#ifndef GEP_FIT_H
#define GEP_FIT_H

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <vector>

struct GEPFitResult {
  double peak;
  double fwhm;
  double mean;
  double sigma;
  double amplitude;
  double mean_error;
  double sigma_error;
  double chi2;
  int    ndf;
  bool   valid;   // false if the histogram had no data and was skipped
};

GEPFitResult FitPeak(TH1D *hist, double refit_nsigma = 2.0) {

  GEPFitResult result = {};
  result.valid = false;

  // Guard against empty (or missing) histograms, e.g. a module that
  // had no hits in this run. Nothing to fit — leave result.valid
  // false and return without touching TF1/Fit at all.
  if (!hist || hist->Integral() <= 0) {
    std::cout << (hist ? hist->GetName() : "(null histogram)")
              << ": no entries, skipping fit." << std::endl;
    return result;
  }

  int maxBin = hist->GetMaximumBin();
  double peakX = hist->GetBinCenter(maxBin);
  double maxContent = hist->GetBinContent(maxBin);
  double halfMax = 0.5 * maxContent;

  int leftBin = maxBin;
  while (leftBin > 1 && hist->GetBinContent(leftBin) > halfMax) leftBin--;

  int rightBin = maxBin;
  while (rightBin < hist->GetNbinsX() && hist->GetBinContent(rightBin) > halfMax) rightBin++;

  double leftX = hist->GetBinCenter(leftBin);
  double rightX = hist->GetBinCenter(rightBin);
  double fwhm = rightX - leftX;

  // Degenerate case (e.g. only one or two bins populated): fall back to
  // a minimal window around the peak so the fit range is never zero
  // width, instead of letting TF1/Fit choke on an empty range.
  if (fwhm <= 0) {
    fwhm = 2.0 * hist->GetBinWidth(maxBin);
    leftX = peakX - fwhm / 2.0;
    rightX = peakX + fwhm / 2.0;
  }

  double sigmaEstimate = fwhm / 2.35482;

  result.peak = peakX;
  result.fwhm = fwhm;

  std::string fit1name = std::string(hist->GetName()) + "_fit1";
  TF1 *fit1 = new TF1(fit1name.c_str(), "gaus", leftX, rightX);
  fit1->SetParameters(maxContent, peakX, sigmaEstimate);
  hist->Fit(fit1, "RQ0");

  double mean1 = fit1->GetParameter(1);
  double sigma1 = std::abs(fit1->GetParameter(2));
  double amp1 = fit1->GetParameter(0);

  if (sigma1 <= 0) sigma1 = sigmaEstimate > 0 ? sigmaEstimate : hist->GetBinWidth(maxBin);

  double refitMin = mean1 - refit_nsigma * sigma1;
  double refitMax = mean1 + refit_nsigma * sigma1;

  std::string fit2name = std::string(hist->GetName()) + "_fit";
  TF1 *fit2 = new TF1(fit2name.c_str(), "gaus", refitMin, refitMax);
  fit2->SetParameters(amp1, mean1, sigma1);

  hist->Fit(fit2, "RQ");

  result.amplitude = fit2->GetParameter(0);
  result.mean = fit2->GetParameter(1);
  result.sigma = std::abs(fit2->GetParameter(2));
  result.mean_error = fit2->GetParError(1);
  result.sigma_error = fit2->GetParError(2);
  result.chi2 = fit2->GetChisquare();
  result.ndf = fit2->GetNDF();
  result.valid = true;

  std::cout << hist->GetName() << std::endl;
  std::cout << "  Peak estimate = " << peakX << std::endl;
  std::cout << "  FWHM estimate = " << fwhm << std::endl;
  std::cout << "  Fit mean      = " << result.mean << " +/- " << result.mean_error << std::endl;
  std::cout << "  Fit sigma     = " << result.sigma << " +/- " << result.sigma_error << std::endl;
  std::cout << "  Fit FWHM      = " << 2.35482 * result.sigma << std::endl;

  delete fit1;

  return result;
}

// Fit the residual distribution of every module (column) of a module-wise
// 2D residual histogram (x = module index, y = residual). Modules with no
// entries come back with .valid == false and are simply skipped — nothing
// throws or crashes when a module happens to be empty.
std::vector<GEPFitResult> FitModuleColumns(TH2D *h2d, int nmodules) {

  std::vector<GEPFitResult> results(nmodules);

  if (!h2d) return results;

  for (int imod = 0; imod < nmodules; imod++) {

    int ix = imod + 1; // TH2 bins are 1-indexed; module 0 -> bin 1

    std::string pname = std::string(h2d->GetName()) + "_proj_mod" + std::to_string(imod);
    TH1D *proj = h2d->ProjectionY(pname.c_str(), ix, ix);

    results[imod] = FitPeak(proj);

    delete proj;
  }

  return results;
}

#endif
