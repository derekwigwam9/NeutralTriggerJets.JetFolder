// 'FoldJets.C'
// Derek Anderson
// 06.16.2016
// 
// This macro performs the unfolding / backfolding of a provided jet spectrum
// using the 'StJetFolder' class.
//
// NOTE: if 'sName' is not provided, then the code will assume the response
// matrix and the smeared spectrum have already been calculated in the file
// 'sIn' and are named 'hResponse' and 'hSmeared' respectively. The argument
// of 'Unfold()' specifies the namecycle of the response matrix contained in
// 'sIn'; if no argument is provided, the code looks for a TH2D by the name
// of 'hResponse'. 
//
// Last updated: 09.28.2016

#include <TSystem>
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TPaveText.h"

using namespace std;


class StJetFolder;

const Bool_t   pdsf = true;
const Bool_t   rcas = false;
// input and output files
const TString  tIn("input/ExponentEfficiencyFixed.root");
const TString  sIn("");
const TString  mIn("input/ExponentEfficiencyFixed.root");
const TString  aIn("input/ExponentEfficiencyFixed.root");
const TString  out("Debug.root");
// input namecycles
const TString  tName("hExp");
const TString  sName("test");
const TString  mName("hExpP");
const TString  aName("hExp");
// jet parameters
const Int_t    nRM   = 1;
const Double_t rJet  = 0.5;
const Double_t aMin  = 0.65;
const Double_t pTmin = 0.2;
const Double_t pTmax = 20.;
const Double_t eTmin = 9.;
const Double_t eTmax = 100.;
// unfolding parameters
const Int_t    method = 1;
const Int_t    k      = 2;
const Int_t    nMC    = 1000000;
const Int_t    nToy   = 10;
// prior parameters
const Int_t    prior  = 3;
const Double_t bPrior = 0.1;
const Double_t nPrior = 5.8;
const Double_t tPrior = 5.;
const Double_t mPrior = 0.140;


void FoldJets() {

  if (pdsf) {
    gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjet.so");
    gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjettools.so");
  }
  if (rcas) {
    gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjet.so");
    gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjettools.so");
  }
  gSystem -> Load("../../RooUnfold/libRooUnfold.so");
  gSystem -> Load("StJetFolder");
  
  // lower verbosity
  gErrorIgnoreLevel = kError;

  TPaveText *label = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC NB");
  label -> AddText("pp collisions, #sqrt{s_{NN}}=200 GeV");
  label -> AddText("#gamma^{dir} trigger, p_{T}^{trg}>9 GeV/c");
  label -> AddText("R=0.7, A_{jet}>1.2, p_{T}^{cst}>0.2 GeV/c");
  label -> AddText("Charge jets");
  label -> SetFillColor(kWhite);

  Double_t    chi2 = 0.;
  StJetFolder f(method, prior);
  f.Init(out, tIn, sIn, mIn);
  f.SetNamecycles(tName, sName, mName);
  f.SetJetParameters(nRM, rJet, aMin, pTmin, pTmax, eTmin, eTmax);
  f.SetUnfoldParameters(k, nMC, nToy, bPrior, nPrior, tPrior);
  f.SetLabel(label);
  f.SetAbsoluteTruth(aIn, aName);
  chi2 = f.Unfold();
  f.Finish();

}

// End ------------------------------------------------------------------------
