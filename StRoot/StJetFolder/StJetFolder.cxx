// 'StJetFolder.cxx'
// Derek Anderson
// 02.14.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// contains the 'Init()', 'Unfold()', 'Backfold()', and 'Finish()'
// routines.
//
// Last updated: 02.17.2017


#define StJetFolder_cxx

// user includes
#include "StJetFolder.h"
#include "StJetFolder.io.h"
#include "StJetFolder.sys.h"
#include "StJetFolder.math.h"
#include "StJetFolder.plot.h"

ClassImp(StJetFolder)

using namespace std;



void StJetFolder::Init() {

  Bool_t inputOK = CheckFlags();
  if (!inputOK) assert(inputOK);

  //_response = new RooUnfoldResponse(_hSmeared, _hPrior, _hResponse);
  _response = new RooUnfoldResponse(0, 0, _hResponse);
  if (_response) {
    PrintInfo(4);
    _flag[10] = true;
  }
  else {
    PrintError(11);
    assert(_response);
  }

}  // end 'Init()'


void StJetFolder::Unfold() {

  PrintInfo(5);

  RooUnfoldBayes    *bay;
  RooUnfoldSvd      *svd;
  RooUnfoldBinByBin *bin;
  RooUnfoldTUnfold  *tun;
  RooUnfoldInvert   *inv;
  switch (_method) {
    case 0:
      _hUnfolded = (TH1D*) _hMeasured -> Clone("hUnfolded");
      break;
    case 1:
      bay        = new RooUnfoldBayes(_response, _hMeasured, _kReg);
      _hUnfolded = (TH1D*) bay -> Hreco();
      break;
    case 2:
      svd        = new RooUnfoldSvd(_response, _hMeasured, _kReg, _nToy);
      _hUnfolded = (TH1D*) svd -> Hreco();
      break;
    case 3:
      bin        = new RooUnfoldBinByBin(_response, _hMeasured);
      _hUnfolded = (TH1D*) bin -> Hreco();
      break;
    case 4:
      tun        = new RooUnfoldTUnfold(_response, _hMeasured, TUnfold::kRegModeDerivative);
      _hUnfolded = (TH1D*) tun -> Hreco();
      break;
    case 5:
      inv        = new RooUnfoldInvert(_response, _hMeasured);
      _hUnfolded = (TH1D*) inv -> Hreco();
      break;
  }

  // TEST
  _hUnfolded -> Divide(_hEfficiency);

  // make sure unfolded didn't exceed max bin
  Int_t nU = _hUnfolded -> GetNbinsX();
  for (Int_t i = 1; i < nU + 1; i++) {
    Double_t uBin = _hUnfolded -> GetBinLowEdge(i);
    if (uBin > Umax) {
      _hUnfolded -> SetBinContent(i, 0.);
      _hUnfolded -> SetBinError(i, 0.);
    }
  }

  PrintInfo(6);

}  // end 'Unfold()'


void StJetFolder::Backfold(Double_t &chi2) {

  PrintInfo(7);

  if (_method == 0) {
    _hBackfolded = (TH1D*) _hMeasured -> Clone("hBackfolded");
    return;
  }

  Int_t    nU  = _hUnfolded -> GetNbinsX();
  Int_t    nM  = _hMeasured -> GetNbinsX();
  Double_t u1  = _hUnfolded -> GetBinLowEdge(1);
  Double_t u2  = _hUnfolded -> GetBinLowEdge(nU + 1);
  Double_t m1  = _hMeasured -> GetBinLowEdge(1);
  Double_t m2  = _hMeasured -> GetBinLowEdge(nM + 1);
  _hNormalize  = new TH1D("hNormalize", "For normalizng backfolded spectrum", nU, u1, u2);
  _hBackfolded = new TH1D("hBackfolded", "Backfolded spectrum", nM, m1, m2);
  _hNormalize  -> Sumw2();
  _hBackfolded -> Sumw2();


  // monte-carlo loop
  Double_t u = 0.;
  Double_t b = 0.;
  for (Int_t i = 0; i < _nMC; i++) {
    u = _hUnfolded -> GetRandom();
    b = Smear(u);
    _hNormalize -> Fill(u);
    if (b > -1000.) _hBackfolded -> Fill(b);
  }

  // normalize backfolded spectrum / apply efficiency
  Double_t iU = _hUnfolded  -> Integral();
  Double_t iN = _hNormalize -> Integral();
  if (iU > 0.) {
    Double_t scale = iU / iN;
    _hBackfolded -> Scale(scale);
  }
  _hBackfolded -> Multiply(_hEfficiency);
  _chi2 = CalculateChi2(_hMeasured, _hBackfolded);
  chi2  = _chi2;

  PrintInfo(8);

}  // end 'Backfold()'


void StJetFolder::Finish() {

  // calculate ratios
  _hBackVsMeasRatio  = CalculateRatio(_hBackfolded, _hMeasured, "hBackVsMeasRatio");
  _hPriVsUnfoldRatio = CalculateRatio(_hPrior, _hUnfolded, "hPriVsUnfoldRatio");
  _hSmearVsMeasRatio = CalculateRatio(_hSmeared, _hMeasured, "hSmearVsMeasRatio");
  PrintInfo(9);

  // set names
  _hPrior            -> SetName("hPrior");
  _hSmeared          -> SetName("hSmeared");
  _hMeasured         -> SetName("hMeasured");
  _hUnfolded         -> SetName("hUnfolded");
  _hEfficiency       -> SetName("hEfficiency");
  _hResponse         -> SetName("hResponse");


  CreateLabel();
  CreatePlots();
  PrintInfo(11);

  // save and close file
  _fOut              -> cd();
  _hPrior            -> Write();
  _hSmeared          -> Write();
  _hMeasured         -> Write();
  _hUnfolded         -> Write();
  _hNormalize        -> Write();
  _hBackfolded       -> Write();
  _hBackVsMeasRatio  -> Write();
  _hPriVsUnfoldRatio -> Write();
  _hSmearVsMeasRatio -> Write();
  _hEfficiency       -> Write();
  _hResponse         -> Write();
  _label             -> Write();
  _fOut              -> Close();
  PrintInfo(12);

}  // end 'Finish()'

// End ------------------------------------------------------------------------
