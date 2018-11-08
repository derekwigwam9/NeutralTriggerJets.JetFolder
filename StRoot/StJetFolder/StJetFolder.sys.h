// 'StJetFolder.sys.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates various internal routines (e.g. printing error messages).
//
// Last updated: 07.15.2018


#pragma once

using namespace std;



void StJetFolder::PrintInfo(const Int_t code) {

  switch (code) {
    case 0:
      cout << "\n  Folder created!\n"
           << "    Writing to '" << _fOut -> GetName() << "'"
           << endl;
      break;
    case 1:
      cout << "    Spectra grabbed..." << endl;
      break;
    case 2:
      cout << "    Info set..." << endl;
      break;
    case 3:
      cout << "    Prior parameters set...\n"
           << "      b = " << _bPrior << ", m = " << _mPrior << "\n"
           << "      n = " << _nPrior << ", t = " << _tPrior
           << endl;
      cout << "    Unfolding parameters set...\n"
           << "      method = " << _method << ", k = " << _kReg << "\n"
           << "      nMc = " << _nMC << ", nToy = " << _nToy << "\n"
           << "      uMax = " << _uMax << ", bMax = " << _bMax
           << endl;
      break;
    case 4:
      cout << "    Folder initialized..." << endl;
      break;
    case 5:
      cout << "    Unfolding..." << endl;
      break;
    case 6:
      cout << "    Unfolding finished!\n"
           << "      Chi2 (unfold) = " << _chi2unfold
           << endl;
      break;
    case 7:
      cout << "    Backfolding..." << endl;
      break;
    case 8:
      cout << "    Backfolding finished!\n"
           << "      Chi2 (backfold) = " << _chi2backfold
           << endl;
      break;
    case 9:
      cout << "    Ratios calculated!" << endl;
      break;
    case 10:
      cout << "    Creating plots..." << endl;
      break;
    case 11:
      cout << "    Plots created; saving..." << endl;
      break;
    case 12:
      cout << "  Folding finished!\n" << endl;
      break;
  }

}  // end 'PrintInfo(Int_t)'


void StJetFolder::PrintError(const Int_t code) {

  switch (code) {
    case 0:
      cerr << "PANIC: couldn't grab prior spectrum!" << endl;
      break;
    case 1:
      cerr << "PANIC: couldn't grab smeared spectrum!" << endl;
      break;
    case 2:
      cerr << "PANIC: couldn't grab measured spectrum!" << endl;
      break;
    case 3:
      cerr << "PANIC: couldn't grab response matrix!" << endl;
      break;
    case 4:
      cerr << "PANIC: couldn't grab effiency!" << endl;
      break;
    case 5:
      cerr << "PANIC: event info couldn't be set!" << endl;
      break;
    case 6:
      cerr << "PANIC: trigger info couldn't be set!" << endl;
      break;
    case 7:
      cerr << "PANIC: jet info couldn't be set!" << endl;
      break;
    case 8:
      cerr << "PANIC: at least one specturm not set!" << endl;
      break;
    case 9:
      cerr << "PANIC: some info not set!" << endl;
      break;
    case 10:
      cerr << "PANIC: some parameters not set!" << endl;
      break;
    case 11:
      cerr << "PANIC: folder couldn't be initialized!" << endl;
      break;
    case 12:
      cerr << "PANIC: trying to take ratio of 2 histograms with different dimensions!" << endl;
      break;
  }

}  // end 'PrintInfo(Int_t)'


void StJetFolder::InitializePriors() {

  // for normalization
  TH1D *hSmearNorm = (TH1D*) _hMeasured -> Clone();
  TH1D *hAfterEff  = (TH1D*) _hMeasured -> Clone();
  hSmearNorm -> Reset("ICE");
  hAfterEff  -> Divide(_hEfficiency);

  // for efficiency
  TH1D *hParEffDif = (TH1D*) _hPrior   -> Clone();
  TH1D *hDetEffDif = (TH1D*) _hSmeared -> Clone();
  hParEffDif -> Reset("ICE");
  hDetEffDif -> Reset("ICE");

  // initialize response and efficiency
  _hEfficiencyDiff = (TH1D*) _hEfficiency -> Clone();
  _hResponseDiff   = (TH2D*) _hResponse   -> Clone();
  _hEfficiencyDiff -> Reset("ICE");
  _hResponseDiff   -> Reset("ICE");


  // create particle-level prior
  switch (_prior) {
    case 1:
      _hPrior -> Reset("ICE");
      for (UInt_t iMC = 0; iMC < _nMC; iMC++) {
        const Double_t l = _fLevy -> GetRandom(XminPrior, _uMax);
        _hPrior -> Fill(l);
      }
      break;
    case 2:
      _hPrior -> Reset("ICE");
      for (UInt_t iMC = 0; iMC < _nMC; iMC++) {
        const Double_t t = _fTsallis -> GetRandom(XminPrior, _uMax);
        _hPrior -> Fill(t);
      }
      break;
    case 3:
      _hPrior -> Reset("ICE");
      for (UInt_t iMC = 0; iMC < _nMC; iMC++) { 
        const Double_t e = _fExponential -> GetRandom(XminPrior, _uMax);
        _hPrior -> Fill(e);
      }
      break;
    case 4:
      _hPrior -> Reset("ICE");
      for (UInt_t iMC = 0; iMC < _nMC; iMC++) {
        const Double_t p = _fPowerLaw -> GetRandom(XminPrior, _uMax);
        _hPrior -> Fill(p);
      }
      break;
  }

  const Double_t iPar   = hAfterEff -> Integral();
  const Double_t iParMC = _hPrior   -> Integral();
  const Double_t scaleP = iPar / iParMC;
  if (iParMC > 0.) _hPrior -> Scale(scaleP);

  // normalize by bin width
  const UInt_t nBinsP = _hPrior -> GetNbinsX();
  for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
    const Double_t pVal = _hPrior -> GetBinContent(iBinP);
    const Double_t pErr = _hPrior -> GetBinError(iBinP);
    const Double_t pBin = _hPrior -> GetBinWidth(iBinP);
    _hPrior -> SetBinContent(iBinP, pVal / pBin);
    _hPrior -> SetBinError(iBinP, pErr / pBin);
  }


  // create detector-level prior and response
  _hSmeared -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMC; iMC++) {
    const Double_t p = _hPrior -> GetRandom();
    const Double_t s = Smear(p);
    const UInt_t   e = ApplyEff(p);
    if ((s > -1000.) && (e == 1)) {
      _hSmeared      -> Fill(s);
      _hResponseDiff -> Fill(s, p);
      hDetEffDif     -> Fill(p);
    }
    hSmearNorm -> Fill(p);
    hParEffDif -> Fill(p);
  }
  _hEfficiencyDiff -> Divide(hDetEffDif, hParEffDif, 1., 1.);

  const Double_t iNorm  = _hPrior    -> Integral();
  const Double_t iDetMC = hSmearNorm -> Integral();
  const Double_t scaleS = iNorm / iDetMC;
  if (iNorm > 0.) _hSmeared -> Scale(scaleS);

}  // end 'InitializePriors()'


Bool_t StJetFolder::CheckFlags() {

  // check spectra
  Bool_t spectraOK = true;
  for (Int_t i = 0; i < 5; i++) {
    if (!_flag[i]) {
      spectraOK = false;
      PrintError(8);
      break;
    }
  }
  if (spectraOK)
    PrintInfo(1);

  // check info
  Bool_t infoOK = true;
  for (Int_t i = 5; i < 8; i++) {
    if (!_flag[i]) {
      infoOK = false;
      PrintError(9);
      break;
    }
  }
  if (infoOK)
    PrintInfo(2);

  // check parameters
  Bool_t parametersOK = true;
  for (Int_t i = 8; i < 10; i++) {
    if (!_flag[i]) {
      parametersOK = false;
      PrintError(10);
      break;
    }
  }
  if (parametersOK)
    PrintInfo(3);


  Bool_t inputOK = (spectraOK && infoOK && parametersOK);
  return inputOK;

}  // end 'CheckFlags()'

// End ------------------------------------------------------------------------
