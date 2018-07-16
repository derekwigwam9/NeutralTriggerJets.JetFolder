// 'StJetFolder.h'
// Derek Anderson
// 02.14.2017
//
// This class handles the unfolding of a provided spectrum.  Please see below
// for definition of names:
//
//   prior    -- prior spectrum, used to seed unfolding algorithm.
//   smear    -- smeared prior spectrum, the prior spectrum with
//               efficiencies, smearing, etc. applied.
//   measure  -- measured spectrum, to be unfolded.
//   unfold   -- unfolded measured spectrum.
//   backfold -- unfolded spectrum with efficiencies, smearing,
//               etc. applied.
//
// Last updated: 07.15.2018


#ifndef StJetFolder_h
#define StJetFolder_h

#include <cmath>
#include <cassert>
#include <iostream>
// ROOT includes
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TPaveText.h"
// RooUnfold includes
#include "../RooUnfold/RooUnfoldResponse.h"
#include "../RooUnfold/RooUnfoldBayes.h"
#include "../RooUnfold/RooUnfoldSvd.h"
#include "../RooUnfold/RooUnfoldBinByBin.h"
#include "../RooUnfold/RooUnfoldTUnfold.h"
#include "../RooUnfold/RooUnfoldInvert.h"

using namespace std;


// global constants
const Int_t    Nflag = 11;
const Double_t Umax  = 100.;
const Double_t Mpion = 0.140;



class StJetFolder {

public:

  StJetFolder(const Char_t *oFile);
  virtual ~StJetFolder();

  // public methods ('StJetFolder.io.h')
  void SetPrior(const Char_t *pFile, const Char_t *pName);
  void SetSmeared(const Char_t *sFile, const Char_t *sName);
  void SetMeasured(const Char_t *mFile, const Char_t *mName);
  void SetResponse(const Char_t *rFile, const Char_t *rName);
  void SetEfficiency(const Char_t *eFile, const Char_t *eName);
  void SetEventInfo(const Int_t beam, const Double_t energy);
  void SetTriggerInfo(const Int_t trigger, const Double_t eTmin, const Double_t eTmax, const Double_t hMax);
  void SetJetInfo(const Int_t type, const Int_t nRM, const Double_t rJet, const Double_t aMin, const Double_t pTmin);
  void SetPriorParameters(const Int_t prior, const Double_t bPrior, const Double_t mPrior, const Double_t nPrior, const Double_t tPrior);
  void SetUnfoldParameters(const Int_t method, const Int_t kReg, const Int_t nMC, const Int_t nToy);
  // public methods ('StJetFolder.cxx')
  void Init();
  void Unfold(Double_t &chi2unfold);
  void Backfold(Double_t &chi2backfold);
  void Finish();

  // static public methods ('StJetFolder.math.h')
  static Double_t Levy(const Double_t *x, const Double_t *p);
  static Double_t Tsallis(const Double_t *x, const Double_t *p);
  static Double_t Exponential(const Double_t *x, const Double_t *p);


private:

  // atomic members
  Int_t     _type;
  Int_t     _prior;
  Int_t     _method;
  Int_t     _kReg;
  Int_t     _nMC;
  Int_t     _nToy;
  Bool_t    _flag[Nflag];
  Double_t  _bPrior;
  Double_t  _mPrior;
  Double_t  _nPrior;
  Double_t  _tPrior;
  Double_t  _chi2unfold;
  Double_t  _chi2backfold;
  // ROOT members
  TF1       *_fLevy;
  TF1       *_fTsallis;
  TF1       *_fExponential;
  TH1D      *_hPrior;
  TH1D      *_hSmeared;
  TH1D      *_hMeasured;  
  TH1D      *_hUnfolded;  
  TH1D      *_hBackfolded;
  TH1D      *_hEfficiency;
  TH1D      *_hNormalize;
  TH1D      *_hBackVsMeasRatio;
  TH1D      *_hUnfoldVsPriRatio;
  TH1D      *_hSmearVsMeasRatio;
  TH2D      *_hResponse;
  TFile     *_fOut;
  TString   *_sEvnt;
  TString   *_sTrig;
  TString   *_sJet1;
  TString   *_sJet2;
  TString   *_sJet3;
  TPaveText *_label;
  // RooUnfold members
  RooUnfoldResponse *_response;

  // private methods ('StJetFolder.sys.h')
  void     PrintInfo(const Int_t code);
  void     PrintError(const Int_t code);
  void     InitializePriors();
  Bool_t   CheckFlags();
  // private methods ('StJetFolder.plot.h')
  void     CreateLabel();
  void     CreatePlots();
  void     ResizeString(TString &str, const Int_t nDec);
  void     DrawHistogram(TH1 *h, const Char_t *option, const Int_t mColor, const Int_t lColor, const Int_t fColor, const Int_t mStyle, const Int_t lStyle, const Int_t fStyle, const Double_t mSize);
  TString* CreateTitle();
  // private methods ('StJetFolder.math.h')
  TH1D*    CalculateRatio(const TH1D *hA, const TH1D *hB, const Char_t *rName);
  Double_t Smear(const Double_t yP);
  Double_t CalculateChi2(const TH1D *hA, TH1D *hB);


  ClassDef(StJetFolder, 1)

};



#endif
#ifdef StJetFolder_cxx

StJetFolder::StJetFolder(const Char_t *oFile) {

  _fOut = new TFile(oFile, "recreate");
  for (Int_t i = 0; i < Nflag; i++) {
    _flag[i] = false;
  }
  PrintInfo(0);

}  // end 'StJetFolder(Char_t*)'


StJetFolder::~StJetFolder() {

}  // end '~StJetFolder()'

#endif

// End ------------------------------------------------------------------------
