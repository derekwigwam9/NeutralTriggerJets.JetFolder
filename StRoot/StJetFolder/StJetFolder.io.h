// 'StJetFolder.io.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates I/O routines.
//
// Last updated: 08.23.2018


#pragma once

using namespace std;



void StJetFolder::SetPrior(const Char_t *pFile, const Char_t *pName) {

  TFile *fPrior = (TFile*) gROOT -> GetListOfFiles() -> FindObject(pFile);
  if (!fPrior || !fPrior->IsOpen()) {
    fPrior = new TFile(pFile);
  }

  TH1D *hPrior;
  if (fPrior) {
    hPrior = (TH1D*) fPrior -> Get(pName);
  }


  if (hPrior) {
    _hPrior  = (TH1D*) hPrior -> Clone();
    _flag[0] = true;
  }
  else {
    PrintError(0);
    assert(hPrior);
  }

}  // end 'SetPrior(Char_t*, Char_t*)'


void StJetFolder::SetSmeared(const Char_t *sFile, const Char_t *sName) {

  TFile *fSmeared = (TFile*) gROOT -> GetListOfFiles() -> FindObject(sFile);
  if (!fSmeared || !fSmeared->IsOpen()) {
    fSmeared = new TFile(sFile);
  }

  TH1D *hSmeared;
  if (fSmeared) {
    hSmeared = (TH1D*) fSmeared -> Get(sName);
  }


  if (hSmeared) {
    _hSmeared = (TH1D*) hSmeared -> Clone();
    _flag[1]  = true;
  }
  else {
    PrintError(1);
    assert(hSmeared);
  }

}  // end 'SetSmeared(Char_t*, Char_t*)'


void StJetFolder::SetMeasured(const Char_t *mFile, const Char_t *mName) {

  TFile *fMeasured = (TFile*) gROOT -> GetListOfFiles() -> FindObject(mFile);
  if (!fMeasured || !fMeasured->IsOpen()) {
    fMeasured = new TFile(mFile);
  }

  TH1D *hMeasured;
  if (fMeasured) {
    hMeasured = (TH1D*) fMeasured -> Get(mName);
  }


  if (hMeasured) {
    _hMeasured = (TH1D*) hMeasured -> Clone();
    _flag[2]   = true;
  }
  else {
    PrintError(2);
    assert(hMeasured);
  }

}  // end 'SetMeasured(Char_t*, Char_t*)'


void StJetFolder::SetResponse(const Char_t *rFile, const Char_t *rName) {

  TFile *fResponse = (TFile*) gROOT -> GetListOfFiles() -> FindObject(rFile);
  if (!fResponse || !fResponse->IsOpen()) {
    fResponse = new TFile(rFile);
  }

  TH2D *hResponse;
  if (fResponse) {
    hResponse = (TH2D*) fResponse -> Get(rName);
  }


  if (hResponse) {
    _hResponse = (TH2D*) hResponse -> Clone();
    _flag[3]   = true;
  }
  else {
    PrintError(3);
    assert(hResponse);
  }

}  // end 'SetResponse(Char_t*, Char_t*)'


void StJetFolder::SetEfficiency(const Char_t *eFile, const Char_t *eName, const Bool_t doSmoothing, const Bool_t removeErrors) {

  TFile *fEfficiency = (TFile*) gROOT -> GetListOfFiles() -> FindObject(eFile);
  if (!fEfficiency || !fEfficiency->IsOpen()) {
    fEfficiency = new TFile(eFile);
  }

  TH1D *hEfficiency;
  if (fEfficiency) {
    hEfficiency = (TH1D*) fEfficiency -> Get(eName);
  }


  if (hEfficiency) {
    _hEfficiency = (TH1D*) hEfficiency -> Clone();
    _flag[4]     = true;
  }
  else {
    PrintError(4);
    assert(hEfficiency);
  }

  // smooth efficiency at high pT
  const Float_t fitGuess(0.87);
  const Float_t fitRange[2] = {10., 30.};
  if (doSmoothing) {
    TF1 *fFit = new TF1("fFit", "[0]", fitRange[0], fitRange[1]);
    fFit -> SetParameter(0, fitGuess);

    _hEfficiency -> Fit(fFit, "RQ0");
    if (_hEfficiency -> GetFunction("fFit")) {
      const UInt_t nBins  = _hEfficiency -> GetNbinsX();
      const UInt_t iStart = _hEfficiency -> FindBin(fitRange[0]);
      for (UInt_t iBin = iStart; iBin < (nBins + 1); iBin++) {
        const Float_t binVal = _hEfficiency -> GetBinContent(iBin);
        const Float_t binErr = _hEfficiency -> GetBinError(iBin);
        const Float_t relErr = binErr / binVal;
        const Float_t newVal = fFit -> GetParameter(0);
        const Float_t newErr = relErr * newVal;
        if (binVal > 0.) {
          _hEfficiency -> SetBinContent(iBin, newVal);
          _hEfficiency -> SetBinError(iBin, newErr);
        }  // end if (val > 0)
      }  // end bin loop
    }  // end if (function)
  }  // end smoothing

  // remove errors
  if (removeErrors) {
    const UInt_t nBins = _hEfficiency -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
      _hEfficiency -> SetBinError(iBin, 0.);
    }  // end bin loop
  }  // end removing bin errors

}  // end 'SetEfficiency(Char_t*, Char_t*)'


void StJetFolder::SetEventInfo(const Int_t beam, const Double_t energy) {

  const Int_t nDecE = 1;


  TString bTxt("");
  TString eStr("");
  TString eTxt("");
  switch (beam) {
    case 0:
      bTxt = "pp collisions, #sqrt{s} = ";
      break;
    case 1:
      bTxt = "AuAu collisions, #sqrt{s_{NN}} = ";
      break;
  }

  eStr += energy;
  ResizeString(eStr, nDecE);
  eTxt.Append(eStr);
  eTxt.Append(" GeV");

  // combine strings
  TString evnt(bTxt);
  evnt.Append(eTxt);
  _sEvnt = new TString(evnt);


  if (_sEvnt) {
    _flag[5] = true;
  }
  else {
    PrintError(5);
    assert(_sEvnt);
  }
  

}  // end 'SetEventInfo(Int_t, Double_t)'


void StJetFolder::SetTriggerInfo(const Int_t trigger, const Double_t eTmin, const Double_t eTmax, const Double_t hMax) {

  // no. of decimals
  const Int_t nDecE = 1;


  TString tTxt("");
  TString eStr("");
  TString eTxt("");
  switch (trigger) {
    case 0:
      tTxt = "#gamma^{dir} trigger, ";
      break;
    case 1:
      tTxt = "#gamma^{rich} trigger, ";
      break;
    case 2:
      tTxt = "#pi^{0} trigger, ";
      break;
    case 3:
      tTxt = "h^{#pm} trigger, ";
      break;
  }
  _trigger = trigger;

  eStr += eTmin;
  ResizeString(eStr, nDecE);
  eTxt.Append("E_{T}^{trg} #in (");
  eTxt.Append(eStr);

  eStr  = "";
  eStr += eTmax;
  ResizeString(eStr, nDecE);
  eTxt.Append(", ");
  eTxt.Append(eStr);
  eTxt.Append(") GeV");

  // combine strings
  TString trig(tTxt);
  trig.Append(eTxt);
  _sTrig = new TString(trig);


  if (_sTrig) {
    _flag[6] = true;
  }
  else {
    PrintError(6);
    assert(_sTrig);
  }

}  // end 'SetTriggerInfo(Int_t, Double_t, Double_t, Double_t, Double_t)'


void StJetFolder::SetJetInfo(const Int_t type, const Int_t nRM, const Double_t rJet, const Double_t aMin, const Double_t pTmin) {

  const Int_t nDecR = 1;
  const Int_t nDecA = 2;
  const Int_t nDecP = 1;


  TString tTxt("");
  TString nTxt("");
  TString rStr("");
  TString rTxt("");
  TString aStr("");
  TString aTxt("");
  TString pStr("");
  TString pTxt("");
  switch (type) {
    case 0:
      tTxt = "#bf{charged jets}";
      break;
    case 1:
      tTxt = "#bf{full jets};";
      break;
  }

  nTxt += "N_{rm} = ";
  nTxt += nRM;

  rStr += rJet;
  ResizeString(rStr, nDecR);
  rTxt.Append("R = ");
  rTxt.Append(rStr);

  aStr += aMin;
  ResizeString(aStr, nDecA);
  aTxt.Append("A_{jet} > ");
  aTxt.Append(aStr);
  aTxt.Append(", ");

  pStr += pTmin;
  ResizeString(pStr, nDecP);
  pTxt.Append("p_{T}^{cst} > ");
  pTxt.Append(pStr);
  pTxt.Append(", ");

  // combine strings
  _sJet1 = new TString("anti-k_{T}, ");
  _sJet2 = new TString(aTxt);
  _sJet3 = new TString(tTxt);
  _sJet1 -> Append(rTxt);
  _sJet2 -> Append(pTxt);
  _sJet2 -> Append(nTxt);


  Bool_t jetInfoIsSet = (_sJet1 && _sJet2 && _sJet3);
  if (jetInfoIsSet) {
    _flag[7] = true;
  }
  else {
    PrintError(7);
    assert(jetInfoIsSet);
  }

}  // end 'SetJetInfo(Int_t, Int_t, Double_t, Double_t, Double_t)'


void StJetFolder::SetPriorParameters(const Int_t prior, const Double_t bPrior, const Double_t mPrior, const Double_t nPrior, const Double_t tPrior) {

  _prior  = prior;
  _bPrior = bPrior;
  _mPrior = mPrior;
  _nPrior = nPrior;
  _tPrior = tPrior;
  if (_prior > 0)
    _differentPrior = true;
  else
    _differentPrior = false;

  // create prior functions
  const Double_t nBinsP = _hPrior -> GetNbinsX();
  const Double_t stopP  = _hPrior -> GetBinLowEdge(nBinsP + 1);
  const Double_t startP = XminPrior;
  _fLevy        = new TF1("fLevy", StJetFolder::Levy, startP, stopP, 4);
  _fTsallis     = new TF1("fTsallis", StJetFolder::Tsallis, startP, stopP, 3);
  _fExponential = new TF1("fExponential", StJetFolder::Exponential, startP, stopP, 2);
  _fPowerLaw    = new TF1("fPowerLaw", StJetFolder::PowerLaw, startP, stopP, 2);
  _fLevy        -> SetParameters(_bPrior, _mPrior, _nPrior, _tPrior);
  _fTsallis     -> SetParameters(_bPrior, _nPrior, _tPrior);
  _fExponential -> SetParameters(_bPrior, _tPrior);
  _fPowerLaw    -> SetParameters(_bPrior, _tPrior);


  _flag[8] = true;

}  // end 'SetPriorParameters(Int_t, Double_t, Double_t, Double_t, Double_t)'


void StJetFolder::SetUnfoldParameters(const Int_t method, const Int_t kReg, const Int_t nMC, const Int_t nToy, const Double_t uMax, const Double_t bMax) {

  _method = method;
  _kReg   = kReg;
  _nMC    = nMC;
  _nToy   = nToy;
  _uMax   = uMax;
  _bMax   = bMax;


  _flag[9] = true;

}  // end 'SetUnfoldParameters(Int_t, Int_t, Int_t)'

// End ------------------------------------------------------------------------

