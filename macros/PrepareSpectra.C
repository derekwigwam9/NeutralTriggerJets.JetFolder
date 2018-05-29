// 'PrepareSpectra.C'
// Derek Anderson
//
// Use this to rebin and/or normalize the response matrix (if necessary)
// as well as calculate efficiencies, etc. before doing the unfolding.
//
// NOTE1: if 'doEffCalc' is equal to 'false', then script will look for
// a histogram with the name 'eName' in the file 'eFile'.  Otherwise, it
// will look for histograms with the names 'eNameP' and 'eNameS' in the
// file 'eName' to calculate the efficiency.
//
// NOTE2: it's assumed that the response matrix will have the same binning
// as the prior and smeared spectra (in the y and x axes respectively)
// since the matrix was probably trained on those 2 spectra.


#include <stdio>
#include <cassert>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"

using namespace std;


// filepaths and namecycles
static const TString pFile("input/match.r03a02rm1.root");
static const TString sFile("input/match.r03a02rm1.root");
static const TString mFile("input/pp200corr.r03a02rm1chrg.d15m1y2017.plots.root");
static const TString rFile("input/match.r03a02rm1.root");
static const TString oFile("forTesting.d17m2y2017.root");
static const TString pName("GeantJets/hJetPtCorrG");
static const TString sName("MatchJets/hJetPtCorrM");
static const TString mName("Gam/hJetPtCorrG");
static const TString rName("hResponsePtc");

// for efficiency
static const Bool_t  doEffCalc = true;
static const TString eFile("input/match.r03a02rm1.root");
static const TString eName("hEfficiencyPt");
static const TString eNameP("EventInfo/hGeantPtCorr");
static const TString eNameS("EventInfo/hMatchPtCorr");



void PrepareSpectra() {

  // lower verbosity
  gErrorIgnoreLevel = kError;


  cout << "\n  Preparing spectra..." << endl;
  TFile *fPrior   = new TFile(pFile.Data(), "read");
  TFile *fSmear   = new TFile(sFile.Data(), "read");
  TFile *fMeasure = new TFile(mFile.Data(), "read");
  TFile *fMatrix  = new TFile(rFile.Data(), "read");
  TFile *fOut     = new TFile(oFile.Data(), "recreate");


  cout << "    Grabbing spectra..." << endl;
  TH1D *hPrior   = (TH1D*) fPrior   -> Get(pName.Data());
  TH1D *hSmear   = (TH1D*) fSmear   -> Get(sName.Data());
  TH1D *hMeasure = (TH1D*) fMeasure -> Get(mName.Data());

  // determine prior and smeared binning
  const Int_t    nP   = hPrior -> GetNbinsX();
  const Int_t    nS   = hSmear -> GetNbinsX();
  const Double_t p1   = hPrior -> GetBinLowEdge(1);
  const Double_t p2   = hPrior -> GetBinLowEdge(nP + 1); 
  const Double_t s1   = hSmear -> GetBinLowEdge(1);
  const Double_t s2   = hSmear -> GetBinLowEdge(nS + 1);
  const Double_t pBin = (p2 - p1) / nP;
  const Double_t sBin = (s2 - s1) / nS;
  cout << "    Prior, Smeared, and measured binning: \n"
       << "      Prior -- nBin = " << nP << ", width = " << pBin
       << ", (x1, x2) = (" << p1 << ", " << p2 << ")\n"
       << "      Smear -- nBin = " << nS << ", width = " << sBin
       << ", (x1, x2) = (" << s1 << ", " << s2 << ")\n"
       << "    === For simplicity, binning should be same! ==="
       << endl;

  Bool_t binsAreSame = ((pBin == sBin) && (p1 == s1) && (p2 == s2));
  if (!binsAreSame) assert(binsAreSame);

  // determine measured binning
  const Int_t    nM   = hMeasure -> GetNbinsX();
  const Double_t m1   = hMeasure -> GetBinLowEdge(1);
  const Double_t m2   = hMeasure -> GetBinLowEdge(nM + 1);
  const Double_t mBin = (m2 - m1) / nM;
  cout << "    Measured binning: \n"
       << "      Measure -- nBin = " << nM << ", width = " << mBin
       << ", (x1, x2) = (" << m1 << ", " << m2 << ")\n"
       << endl;


  cout << "    Grabbing response matrix..." << endl;
  TH2D *hResponse = (TH2D*) fMatrix -> Get(rName.Data());

  // determine response binning
  const Int_t    nX   = hResponse -> GetNbinsX();
  const Int_t    nY   = hResponse -> GetNbinsY();
  const Double_t x1   = hResponse -> GetXaxis() -> GetBinLowEdge(1);
  const Double_t x2   = hResponse -> GetXaxis() -> GetBinLowEdge(nX + 1);
  const Double_t y1   = hResponse -> GetYaxis() -> GetBinLowEdge(1);
  const Double_t y2   = hResponse -> GetYaxis() -> GetBinLowEdge(nY + 1);
  const Double_t xBin = (x2 - x1) / nX;
  const Double_t yBin = (y2 - y1) / nY;
  cout << "    Response binning: \n"
       << "      X -- nBin = " << nX << ", width = " << xBin
       << ", (x1, x2) = (" << x1 << ", " << x2 << ")\n"
       << "      Y -- nBin = " << nY << ", width = " << yBin
       << ", (y1, y2) = (" << y1 << ", " << y2 << ")\n"
       << "    === For simplicity, binning should be same! ==="
       << endl;

  binsAreSame = ((xBin == yBin) && (x1 == y1) && (x2 == y2));
  if (binsAreSame) assert(binsAreSame);


  // determine how to rebin
  const Double_t newBin = TMath::Max(mBin, xBin);
  const Double_t oldBin = TMath::Min(mBin, xBin);
  const Double_t newX1  = TMath::Max(m1, x1);
  const Double_t newX2  = TMath::Min(m2, x2);
  const Int_t    newNX  = (Int_t) ((newX2 - newX1) / newBin);
  const Int_t    nGroup = (Int_t) (newBin / oldBin);
  cout << "    New binning:\n"
       << "      nBin = " << newNX << ", width = " << newBin << ", (x1, x2) = (" << newX1 << ", " << newX2 << ")\n"
       << "      nGroup = " << nGroup
       << endl;

  // determine what to rebin
  Bool_t rebinMeasured = false;
  Bool_t rebinResponse = false;
  if (mBin > xBin)
    rebinResponse = true;
  else
    rebinMeasured = true;

  // do rebinning
  if (rebinResponse) {
    hPrior    -> Rebin(nGroup);
    hSmear    -> Rebin(nGroup);
    hResponse -> Rebin2D(nGroup, nGroup);
  }
  else if (rebinMeasured) {
    hMeasure  -> Rebin(nGroup);
  }

  // check if axis needs to be truncated / shifted
  if ((p1 < newX1) || (p2 > newX2)) {

    TH1D *hNewP = new TH1D("hPrior", "Prior, re-binned", newNX, newX1, newX2);
    TH1D *hNewS = new TH1D("hSmeared", "Smeared, re-binned", newNX, newX1, newX2);
    for (Int_t i = 1; i <= newNX; i++) {
      const Double_t xP = hNewP -> GetBinCenter(i);
      const Double_t xS = hNewS -> GetBinCenter(i);
      const Int_t    iP = hPrior -> FindBin(xP);
      const Int_t    iS = hSmear -> FindBin(xS);
      const Double_t yP = hPrior -> GetBinContent(iP);
      const Double_t yS = hSmear -> GetBinContent(iS);
      const Double_t eP = hPrior -> GetBinError(iP);
      const Double_t eS = hSmear -> GetBinError(iS);
      hNewP -> SetBinContent(i, yP);
      hNewS -> SetBinContent(i, yS);
      hNewP -> SetBinError(i, eP);
      hNewS -> SetBinError(i, eS);
    }
    hPrior = hNewP;
    hSmear = hNewS;

  }
  if ((m1 < newX1) || (m2 > newX2)) {

    TH1D *hNewM = new TH1D("hMeasured", "Measured, re-binned", newNX, newX1, newX2);
    for (Int_t i = 1; i <= newNX; i++) {
      const Double_t xM = hNewM -> GetBinCenter(i);
      const Int_t    iM = hMeasure -> FindBin(xM);
      const Double_t yM = hMeasure -> GetBinContent(iM);
      const Double_t eM = hMeasure -> GetBinError(iM);
      hNewM -> SetBinContent(i, yM);
      hNewM -> SetBinError(i, eM);
    }
    hMeasure = hNewM;

  }

  // set user ranges
  hPrior    -> GetXaxis() -> SetRangeUser(newX1, newX2);
  hSmear    -> GetXaxis() -> SetRangeUser(newX1, newX2);
  hMeasure  -> GetXaxis() -> SetRangeUser(newX1, newX2);
  hResponse -> GetXaxis() -> SetRangeUser(newX1, newX2);
  hResponse -> GetYaxis() -> SetRangeUser(newX1, newX2);
  cout << "    Rebinned spectra!" << endl;


  cout << "    Normalizing response matrix..." << endl;
  const Int_t nNewX = hResponse -> GetNbinsX();
  const Int_t nNewY = hResponse -> GetNbinsY();
  for (Int_t j = 1; j <= nNewY; j++) {

    // project prior (y) bin and normalize
    TH1D *hX = hResponse -> ProjectionX("", j, j);
    const Double_t norm = hX -> Integral();
    if (norm == 0)
      continue;
    else
      hX -> Scale(1. / norm);

    // update response matrix's content
    for (Int_t i = 1; i <= nNewX; i++) {
      const Double_t cnt = hX -> GetBinContent(i);
      const Double_t err = hX -> GetBinError(i);
      hResponse -> SetBinContent(i, j, cnt);
      hResponse -> SetBinError(i, j, err);
    }

  }
  cout << "    Response matrix normalized!" << endl;


  TH1D  *hEff;
  TH1D  *hEffP;
  TH1D  *hEffS;
  TFile *fEff = new TFile(eFile.Data(), "read");
  if (doEffCalc) {

    cout << "    Calculating efficiency..." << endl;
    hEffP = (TH1D*) fEff -> Get(eNameP.Data());
    hEffS = (TH1D*) fEff -> Get(eNameS.Data());

    // rebin if necessary
    const Int_t    nEP   = hEffP -> GetNbinsX();
    const Int_t    nES   = hEffS -> GetNbinsX();
    const Double_t ep1   = hEffP -> GetBinLowEdge(1);
    const Double_t ep2   = hEffP -> GetBinLowEdge(nEP + 1);
    const Double_t es1   = hEffS -> GetBinLowEdge(1);
    const Double_t es2   = hEffS -> GetBinLowEdge(nES + 1);
    const Double_t epBin = (ep2 - ep1) / nEP;
    const Double_t esBin = (es2 - es1) / nES;
    cout << "    Efficiency binning: \n"
         << "      Prior -- nBin = " << nEP << ", width = " << epBin
         << ", (x1, x2) = (" << ep1 << ", " << ep2 << ")\n"
         << "      Smear -- nBin = " << nES << ", width = " << esBin
         << ", (x1, x2) = (" << es1 << ", " << es2 << ")\n"
         << "    === Prior and Smear binning should be same!  ===\n"
         << "    === Also should be same as response binning! ==="
         << endl;

    Bool_t psBinsAreSame  = ((epBin == esBin) && (ep1 == es1) && (ep2 == es2));
    Bool_t rxBinsAreSame  = ((esBin == xBin) && (es1 == x1) && (es2 == x2));
    Bool_t ryBinsAreSame  = ((epBin == yBin) && (ep1 == y1) && (ep2 == y2));
    Bool_t effBinsAreSame = (psBinsAreSame && rxBinsAreSame && ryBinsAreSame);
    if (effBinsAreSame) assert(effBinsAreSame);

    // rebin spectra if necessary
    if (rebinResponse) {

      hEffP -> Rebin(nGroup);
      hEffS -> Rebin(nGroup);
      if ((ep1 < newX1) || (ep2 > newX2)) {
        TH1D *hNewEP = new TH1D("hPrior", "Prior, re-binned", newNX, newX1, newX2);
        TH1D *hNewES = new TH1D("hSmeared", "Smeared, re-binned", newNX, newX1, newX2);
        for (Int_t i = 1; i <= newNX; i++) {
          const Double_t xEP = hNewP -> GetBinCenter(i);
          const Double_t xES = hNewS -> GetBinCenter(i);
          const Int_t    iEP = hPrior -> FindBin(xP);
          const Int_t    iES = hSmear -> FindBin(xS);
          const Double_t yEP = hPrior -> GetBinContent(iP);
          const Double_t yES = hSmear -> GetBinContent(iS);
          const Double_t eEP = hPrior -> GetBinError(iP);
          const Double_t eES = hSmear -> GetBinError(iS);
          hNewEP -> SetBinContent(i, yEP);
          hNewES -> SetBinContent(i, yES);
          hNewEP -> SetBinError(i, eEP);
          hNewES -> SetBinError(i, eES);
        }
        hEffP = hNewEP;
        hEffS = hNewES;
      }

    }
    hEffP -> GetXaxis() -> SetRangeUser(newX1, newX2);
    hEffS -> GetXaxis() -> SetRangeUser(newX1, newX2);

    // calculate efficiency
    hEff = (TH1D*) hEffP -> Clone();
    hEff -> Divide(hEffS, hEffP, 1., 1.);
    hEff -> SetTitle("Reconstruction efficiency, #epsilon = N_{matched}/N_{pythia}");
    cout << "    Efficiency calculated!" << endl;

  }  // end eff calc
  else {

    cout << "    Grabbing efficiency..." << endl;
    hEff = (TH1D*) fEff -> Get(eName.Data());
    const Int_t    nE   = hEff -> GetNbinsX();
    const Double_t e1   = hEff -> GetBinLowEdge(1);
    const Double_t e2   = hEff -> GetBinLowEdge(nE + 1);
    const Double_t eBin = (e2 - e1) / nE;
    cout << "    Efficiency binning:\n"
         << "      Eff. -- nBin = " << nE << ", width = " << eBin
         << ", (x1, x2) = (" << e1 << ", " << e2 << ")"
         << "    === Should be same as response binning! ==="
         << endl;

    Bool_t resBinsAreSame = ((eBin == yBin) && (e1 == y1) && (e2 == y2));
    if (!resBinsAreSame) assert(resBinsAreSame);

    // rebin if necessary
    if (rebinResponse) {

      hEff -> Rebin(nGroup);
      if ((e1 < newX1) || (e2 > newX2)) {
        TH1D *hNewE = new TH1D("hEfficiency", "Effciency, re-binned", newNX, newX1, newX2);
        for (Int_t i = 1; i <= newNX; i++) {
          const Double_t xE = hNewE -> GetBinCenter(i);
          const Int_t    iE = hEff  -> FindBin(xP);
          const Double_t yE = hEff  -> GetBinContent(iP);
          const Double_t eE = hEff  -> GetBinError(iP);
          hNewE -> SetBinContent(i, yE);
          hNewE -> SetBinError(i, eE);
        }
        hEff = hNewE;
      }

    }
    hEff -> GetXaxis() -> SetRangeUser(newX1, newX2);

  }  // end else statement


  // set names and titles
  hPrior    -> SetName("hPrior");
  hSmear    -> SetName("hSmeared");
  hMeasure  -> SetName("hMeasured");
  hEff      -> SetName("hEfficiency");
  hResponse -> SetName("hResponse");

  cout << "    Saving spectra..." << endl;
  fOut      -> cd();
  hPrior    -> Write();
  hSmear    -> Write();
  hMeasure  -> Write();
  hEff      -> Write();
  hResponse -> Write();
  fOut      -> Close();


  cout << "  Finished preparations!\n" << endl;

}

// End ------------------------------------------------------------------------
