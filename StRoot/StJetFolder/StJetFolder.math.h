// 'StJetFolder.math.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates various mathematical routines.
//
// Last updated: 05.27.2018


#pragma once

using namespace std;



TH1D* StJetFolder::CalculateRatio(const TH1D *hA, const TH1D *hB, const Char_t *rName) {

  // check denominator and numerator dimension
  const Int_t    nA = hA -> GetNbinsX();
  const Int_t    nB = hB -> GetNbinsX();
  const Double_t a1 = hA -> GetBinLowEdge(1);
  const Double_t a2 = hA -> GetBinLowEdge(nA + 1);
  const Double_t b1 = hB -> GetBinLowEdge(1);
  const Double_t b2 = hB -> GetBinLowEdge(nB + 1);

  Bool_t hasSameDimension = ((nA == nB) && (a1 == b1) && (a2 == b2));
  if (!hasSameDimension) {
    PrintError(12);
    assert(hasSameDimension);
  }


  // initialize ratio histogram
  TH1D *hR = (TH1D*) hA -> Clone();
  hR -> SetName(rName);
  hR -> Reset("ICE");

  // calculate ratio
  Int_t nRpts = 0;
  for (Int_t i = 1; i < nA + 1; i++) {

    const Double_t xA = hA -> GetBinCenter(i);
    const Double_t yA = hA -> GetBinContent(i);
    const Double_t yB = hB -> GetBinContent(i);
    const Double_t eA = hA -> GetBinError(i);
    const Double_t eB = hB -> GetBinError(i);
    if ((yA <= 0.) || (yB <= 0.)) continue;

    const Double_t r   = yA / yB;
    const Double_t rA  = eA / yA;
    const Double_t rB  = eB / yB;
    const Double_t eR2 = pow(r, 2) * (pow(rA, 2) + pow(rB, 2));
    const Double_t eR  = sqrt(eR2);

    const Int_t iR = hR -> FindBin(xA);
    hR -> SetBinContent(iR, r);
    hR -> SetBinError(iR, eR);
    nRpts++;

  }

  hR -> SetEntries(nRpts);
  return hR;

}  // end 'CalculateRatio(TH1D*, TH1D*, TH1D*)'


Double_t StJetFolder::Smear(const Double_t yP) {

  TH1D     *hSmear;
  Int_t    iPrior;
  Int_t    nSmear;
  Double_t xS;

  iPrior = _hResponse -> GetYaxis() -> FindBin(yP);
  hSmear = (TH1D*) _hResponse -> ProjectionX("hSmear", iPrior, iPrior);
  nSmear = hSmear -> GetEntries();
  if (nSmear < 1)
    xS = -1000.;
  else
    xS = hSmear -> GetRandom();

  if (xS > Umax)
    xS = -1000.;
  return xS;

}  // end 'Smear(Double_t)'


Double_t StJetFolder::CalculateChi2(const TH1D *hA, TH1D *hB) {

  // determine where to start and stop comparing
  const Int_t    aMin = hA -> FindFirstBinAbove(0.);
  const Int_t    aMax = hA -> FindLastBinAbove(0.);
  const Int_t    bMin = hB -> FindFirstBinAbove(0.);
  const Int_t    bMax = hB -> FindLastBinAbove(0.);
  const Double_t iMin = TMath::Max(aMin, bMin);
  const Double_t iMax = TMath::Min(aMax, bMax);
  const Double_t xMin = hA -> GetBinCenter(iMin);
  const Double_t xMax = hA -> GetBinCenter(iMax);


  // calculate chi2
  Int_t    nChi = 0;
  Int_t    nA   = hA -> GetNbinsX();
  Double_t chi2 = 0;
  for (Int_t i = 0; i < nA+1; ++i) {

    const Double_t xA  = hA -> GetBinCenter(i);
    const Double_t yA  = hA -> GetBinContent(i);
    const Double_t yLo = hA -> GetBinError(i);
    const Double_t yHi = yLo;
    if (xA < xMin)
      continue;
    if (xA > xMax)
      continue;

    const Int_t    j  = hB -> FindBin(xA);
    const Double_t yB = hB -> GetBinContent(j);
    const Double_t eB = hB -> GetBinError(j);
    if ((yA < 0) || (yB < 0))
      continue;

    Double_t eUse = 0.;
    if (yB > yA)
      eUse = yHi;
    else
      eUse = yLo;

    if ((eUse > 0.) && (eB > 0.)) {
      const Double_t dErr = sqrt(eB*eB + eUse*eUse);
      const Double_t num  = pow(yA - yB, 2.);
      const Double_t den  = pow(dErr, 2.);
      const Double_t c2   = num / den;

      chi2 += c2;
      ++nChi;
    }

  }
  if (nChi > 0) chi2 /= (Double_t) nChi;
  return chi2;

}  // end 'CalculateChi2(TH1D*, TH1D*)'

// End ------------------------------------------------------------------------
 
