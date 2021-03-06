// 'StJetFolder.cxx'
// Derek Anderson
// 07.15.2018
//
// This class handles the unfolding of a provided spectrum.  This file
// contains the 'Init()', 'Unfold()', 'Backfold()', and 'Finish()'
// routines.   Pearson Coefficient calculation adapted from Rhagav K.
// Elayavalli.
//
// Last updated: 10.03.2018


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

  // initialize response
  if (_differentPrior) {
    InitializePriors();
    _response = new RooUnfoldResponse(0, 0, _hResponseDiff);
  }
  else
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


void StJetFolder::Unfold(Double_t &chi2unfold) {

  PrintInfo(5);

  // do unfolding
  RooUnfoldBayes    *bay;
  RooUnfoldSvd      *svd;
  RooUnfoldBinByBin *bin;
  RooUnfoldTUnfold  *tun;
  RooUnfoldInvert   *inv;
  RooUnfoldErrors   *err;
  TMatrixD          *cov;
  switch (_method) {
    case 0:
      _hUnfolded = (TH1D*) _hMeasured -> Clone("hUnfolded");
      break;
    case 1:
      bay        = new RooUnfoldBayes(_response, _hMeasured, _kReg);
      err        = new RooUnfoldErrors(_nToy, bay);
      cov        = (TMatrixD*) bay -> Ereco().Clone();
      _hUnfolded = (TH1D*)     bay -> Hreco();
      break;
    case 2:
      svd        = new RooUnfoldSvd(_response, _hMeasured, _kReg, _nToy);
      err        = new RooUnfoldErrors(_nToy, svd);
      cov        = (TMatrixD*) svd -> Ereco().Clone();
      _hUnfolded = (TH1D*)     svd -> Hreco();
      break;
    case 3:
      bin        = new RooUnfoldBinByBin(_response, _hMeasured);
      err        = new RooUnfoldErrors(_nToy, bin);
      cov        = (TMatrixD*) bin -> Ereco().Clone();
      _hUnfolded = (TH1D*)     bin -> Hreco();
      break;
    case 4:
      tun        = new RooUnfoldTUnfold(_response, _hMeasured, TUnfold::kRegModeDerivative);
      err        = new RooUnfoldErrors(_nToy, tun);
      cov        = (TMatrixD*) tun -> Ereco().Clone();
      _hUnfolded = (TH1D*)     tun -> Hreco();
      break;
    case 5:
      inv        = new RooUnfoldInvert(_response, _hMeasured);
      err        = new RooUnfoldErrors(_nToy, inv);
      cov        = (TMatrixD*) inv -> Ereco().Clone();
      _hUnfolded = (TH1D*)     inv -> Hreco();
      break;
  }

  // correct for efficiency, calculate pearson coef's, grab errors, and grab D vector
  if (_differentPrior)
    _hUnfolded -> Divide(_hEfficiencyDiff);
  else
    _hUnfolded -> Divide(_hEfficiency);


  _hPearson  = GetPearsonCoefficient(cov, _pearsonDebug, "hPearson");
  switch (_method) {
    case 0:
      _hUnfoldErrors = (TH1D*) _hUnfolded -> Clone();
      _hSVvector     = (TH1D*) _hUnfolded -> Clone();
      _hDvector      = (TH1D*) _hUnfolded -> Clone();
      break;
    case 1:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) bay -> Hreco();
      _hDvector      = (TH1D*) bay -> Hreco();
      break;
    case 2:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) svd -> Impl() -> GetSV();
      _hDvector      = (TH1D*) svd -> Impl() -> GetD();
      break;
    case 3:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) bin -> Hreco();
      _hDvector      = (TH1D*) bin -> Hreco();
      break;
    case 4:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) tun -> Hreco();
      _hDvector      = (TH1D*) tun -> Hreco();
      break;
    case 5:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) inv -> Hreco();
      _hDvector      = (TH1D*) inv -> Hreco();
      break;
  }


  // make sure unfolded didn't exceed max bin
  Int_t nU = _hUnfolded -> GetNbinsX();
  for (Int_t i = 1; i < nU + 1; i++) {
    Double_t uBin = _hUnfolded -> GetBinLowEdge(i);
    if (uBin > _uMax) {
      _hUnfolded -> SetBinContent(i, 0.);
      _hUnfolded -> SetBinError(i, 0.);
    }
  }

  // calculate chi2
  _chi2unfold = CalculateChi2(_hPrior, _hUnfolded);
  chi2unfold  = _chi2unfold;

  PrintInfo(6);

}  // end 'Unfold(Double_t)'


void StJetFolder::Backfold(Double_t &chi2backfold) {

  PrintInfo(7);

  if (_method == 0) {
    _hBackfolded = (TH1D*) _hMeasured -> Clone("hBackfolded");
    return;
  }

  _hNormalize  = (TH1D*) _hUnfolded -> Clone();
  _hBackfolded = (TH1D*) _hMeasured -> Clone();
  _hNormalize  -> SetNameTitle("hNormalize", "For normalizing backfolded spectrum");
  _hBackfolded -> SetNameTitle("hBackfolded", "Backfolded spectrum");
  _hNormalize  -> Reset("ICE");
  _hBackfolded -> Reset("ICE");


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

  // apply efficiency
  _hBackfolded -> Multiply(_hEfficiency);
  _chi2backfold = CalculateChi2(_hMeasured, _hBackfolded);
  chi2backfold  = _chi2backfold;

  PrintInfo(8);

}  // end 'Backfold(Double_t)'


void StJetFolder::Finish() {

  // calculate ratios
  _hBackVsMeasRatio   = CalculateRatio(_hBackfolded, _hMeasured, "hBackVsMeasRatio");
  _hUnfoldVsPriRatio  = CalculateRatio(_hUnfolded, _hPrior, "hUnfoldVsPriRatio");
  _hSmearVsMeasRatio  = CalculateRatio(_hSmeared, _hMeasured, "hSmearVsMeasRatio");
  _hUnfoldVsMeasRatio = CalculateRatio(_hUnfolded, _hMeasured, "hUnfoldVsMeasRatio");
  _hSmearVsPriRatio   = CalculateRatio(_hSmeared, _hPrior, "hSmearVsPriRatio");
  PrintInfo(9);


  // set names
  _hPrior        -> SetName("hPrior");
  _hSmeared      -> SetName("hSmeared");
  _hMeasured     -> SetName("hMeasured");
  _hUnfolded     -> SetName("hUnfolded");
  _hDvector      -> SetName("hDvector");
  _hSVvector     -> SetName("hSVvector");
  _hUnfoldErrors -> SetName("hUnfoldErrors");
  _hEfficiency   -> SetName("hEfficiency");
  _hPearson      -> SetName("hPearson");
  _hResponse     -> SetName("hResponse");
  if (_differentPrior) {
    _hEfficiencyDiff -> SetName("hEfficiencyDiff");
    _hResponseDiff   -> SetName("hResponseDiff");
  }

  CreateLabel();
  CreatePlots();
  PrintInfo(11);


  // save and close file
  _fOut               -> cd();
  _hPrior             -> Write();
  _hSmeared           -> Write();
  _hMeasured          -> Write();
  _hUnfolded          -> Write();
  _hNormalize         -> Write();
  _hBackfolded        -> Write();
  _hBackVsMeasRatio   -> Write();
  _hUnfoldVsPriRatio  -> Write();
  _hSmearVsMeasRatio  -> Write();
  _hUnfoldVsMeasRatio -> Write();
  _hPearson           -> Write();
  _hDvector           -> Write();
  _hSVvector          -> Write();
  _hUnfoldErrors      -> Write();
  _hEfficiency        -> Write();
  _hResponse          -> Write();
  if (_differentPrior) {
    _hResponseDiff    -> Write();
    _hEfficiencyDiff  -> Write();
  }
  _fOut               -> Close();
  PrintInfo(12);

}  // end 'Finish()'

// End ------------------------------------------------------------------------
