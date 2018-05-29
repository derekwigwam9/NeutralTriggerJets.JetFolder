// 'Debug.C'
// Derek Anderson
// 02.14.2017
//
// Simple macro for testing StJetFolder

#include <TSystem>
#include <iostream>
#include "TFile.h"
#include "TString.h"

using namespace std;

class StJetFolder;

// filepaths
static const TString  oFile("debug.root");
static const TString  pFile("forTesting.root");
static const TString  sFile("forTesting.root");
static const TString  mFile("forTesting.root");
static const TString  rFile("forTesting.root");
static const TString  eFile("forTesting.root");
static const TString  pName("hPri");
static const TString  sName("hSmeS");
static const TString  mName("hSmeS");
static const TString  rName("hResS");
static const TString  eName("hEffS");
// event info
static const Int_t    beam    = 0;
static const Double_t energy  = 200.;
// trigger info
static const Int_t    trig    = 0;
static const Double_t eTmin   = 9.;
static const Double_t eTmax   = 30.;
static const Double_t hTrgMax = 0.9;
// jet info
static const Int_t    type    = 0;
static const Int_t    nRM     = 1;
static const Double_t rJet    = 0.3;
static const Double_t aMin    = 0.2;
static const Double_t pTmin   = 0.2;
// prior parameters
static const Int_t    prior   = 0;
static const Double_t bPrior  = 0.1;
static const Double_t mPrior  = 0.14;
static const Double_t nPrior  = 5.6;
static const Double_t tPrior  = 0.4;
// unfolding parameters
static const Int_t    method  = 3;
static const Int_t    kReg    = 1;
static const Int_t    nMC     = 10000;
static const Int_t    nToy    = 10;



void Debug() {

  gSystem -> Load("../../RooUnfold/libRooUnfold.so");
  gSystem -> Load("StJetFolder");

  StJetFolder folder(oFile.Data());
  // set spectra
  folder.SetPrior(pFile.Data(), pName.Data());
  folder.SetSmeared(sFile.Data(), sName.Data());
  folder.SetMeasured(mFile.Data(), mName.Data());
  folder.SetResponse(rFile.Data(), rName.Data());
  folder.SetEfficiency(eFile.Data(), eName.Data());
  // set info and parameters
  folder.SetEventInfo(beam, energy);
  folder.SetTriggerInfo(trig, eTmin, eTmax, hTrgMax);
  folder.SetJetInfo(type, nRM, rJet, aMin, pTmin);
  folder.SetPriorParameters(prior, bPrior, mPrior, nPrior, tPrior);
  folder.SetUnfoldParameters(method, kReg, nMC, nToy);
  // preform unfolding
  folder.Init();
  folder.Unfold();
  folder.Backfold();
  folder.Finish();

}

// End ------------------------------------------------------------------------
