#ifndef GEP_FIT_H
#define GEP_FIT_H

#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <string>

struct GEPFitResult {
  double peak;
  double fwhm;
  double mean;
  double sigma;
  double amplitude;
  double mean_error;
  double sigma_error;
};

GEPFitResult FitPeak(TH1D *hist, double refit_nsigma = 2.0) {

  GEPFitResult result = {};

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

  std::cout << hist->GetName() << std::endl;
  std::cout << "  Peak estimate = " << peakX << std::endl;
  std::cout << "  FWHM estimate = " << fwhm << std::endl;
  std::cout << "  Fit mean      = " << result.mean << " +/- " << result.mean_error << std::endl;
  std::cout << "  Fit sigma     = " << result.sigma << " +/- " << result.sigma_error << std::endl;
  std::cout << "  Fit FWHM      = " << 2.35482 * result.sigma << std::endl;

  delete fit1;

  return result;
}

#endif
