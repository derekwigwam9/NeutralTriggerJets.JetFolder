// 'MatchJets.C'
// Derek Anderson
// 07.13.2016
// 
// This macro performs just the jet-matching calculation used to calculate
// the response matrix in unfolding, which can be very time-intensive.
// Here, tName and mName are used solely to define the binning of the
// response matrix / matrices.
//
// NOTE1: The argument of 'Unfold()' specifies the
// namecycle of the response matrix contained in 'sIn'; if no argument
// is provided, the code looks for a TH2D by the name of 'hResponse'.
//
// NOTE2: 'prior' controls what source to use as the prior spectrum.
// Use 0 for Pythia, 1 for a Levy distribution, and 2 for a Tsallis
// distribution.
//
// Last updated: 09.28.2016

#include <TSystem>
#include <iostream>
#include "TString.h"
#include "TPaveText.h"

using namespace std;

class StJetFolder;

const Bool_t   pdsf = true;
const Bool_t   rcas = false;
// input and output files
const TString  tIn("../JetRun/output/Pythia20p.r05a065rm1.gCharged.root");
const TString  sIn("../JetRun/input/Pythia20p.gMerged.root");
const TString  mIn("../JetRun/output/Pythia20p.r05a065rm1.gCharged.root");
const TString  out("Pythia20r.r05a065rm1.gBetterBins.root");
// input namecycles
const TString  tName("QA/hPtCorr");
const TString  sName("ParTree");
const TString  mName("QA/hPtCorr");
// jet parameters
const Int_t    nRM   = 1;
const Double_t rJet  = 0.5;
const Double_t aMin  = 0.65;
const Double_t pTmin = 0.2;
const Double_t pTmax = 20.;
const Double_t eTmin = 9.;
const Double_t eTmax = 20.;
// unfolding parameters
const Int_t    k      = 1;
const Int_t    nMC    = 1000;
const Int_t    nToy   = 10.;
// prior function parameters
const Int_t    prior  = 0;
const Double_t bPrior = 0.1;
const Double_t nPrior = 6.5;
const Double_t tPrior = 1.05;
const Double_t mPrior = 0.14;


void MatchJets() {

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
  gErrorIgnoreLevel = kSysError;


  TPaveText *label = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC NB");
  label -> AddText("pp collisions, #sqrt{s_{NN}}=200 GeV");
  label -> AddText("#gamma^{dir} trigger, p_{T}^{trg}>9 GeV/c");
  label -> AddText("R=0.5, A_{jet}>0.2, p_{T}^{cst}>0.2 GeV/c");
  label -> AddText("Charge jets");
  label -> SetFillColor(kWhite);

  Double_t    chi2(0.);
  StJetFolder f(0, prior);
  f.Init(out, tIn, sIn, mIn);
  f.SetNamecycles(tName, sName, mName);
  f.SetJetParameters(nRM, rJet, aMin, pTmin, pTmax, eTmin, eTmax);
  f.SetUnfoldParameters(k, nMC, nToy, bPrior, nPrior, tPrior, mPrior);
  f.SetLabel(label);
  chi2 = f.Unfold();
  f.Finish();

}

// End ------------------------------------------------------------------------
