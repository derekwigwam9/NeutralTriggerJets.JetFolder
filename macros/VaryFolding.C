// 'VaryFolding.C'
// Derek Anderson
// 08.08.2016
//
// This macro automates the unfolding of a measured jet spectrum for
// numerous different parameters.
//
// Last updated: 09.28.2016

#include <iostream>
#include "TSystem.h"
#include "TString.h"

using namespace std;

class StJetFolder;


// input and output files
const TString  tIn("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/UnfoldingMacro/input/ChargedParticlePt.root");
const TString  sIn("");
const TString  mIn("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/UnfoldingMacro/input/ChargedParticlePt.root");
const TString  aIn("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/UnfoldingMacro/input/ChargedParticlePt.root");
const TString  oPath("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/UnfoldingMacro");
const TString  oName("Debug");
// input namecycles
const TString  tName("hPtPar");
const TString  sName("tpc");
const TString  mName("hPtDet");
const TString  aName("hPtPar");
// jet parameters
const Int_t    nRM    = 1;
const Double_t rJet   = 0.3;
const Double_t aMin   = 0.2;
const Double_t pTmin  = 0.2;
const Double_t pTmax  = 20.;
const Double_t eTmin  = 9.;
const Double_t eTmax  = 100.;
// unfolding parameters
const Int_t    nToy   = 10.;
const Int_t    nMC    = 1000000;
// levy function parameters;
const Int_t    prior  = 1;
const Double_t bPrior = 0.1;
const Double_t mPrior = 0.14;



void VaryFolding(const Int_t p=0, const Int_t m=1, const Int_t k=1, const Double_t n=5.8, const Double_t t=0.4, const TString outDir=oPath) {

  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjet.so");
  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjettools.so");
  gSystem -> Load("/global/project/projectdirs/star/pwg/starjetc/dmawxc/RooUnfold/libRooUnfold.so");
  gSystem -> Load("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/UnfoldingMacro/.sl64_gcc482/lib/libStJetFolder.so");

  // lower verbosity
  gErrorIgnoreLevel = kError;


  // create output name
  Int_t Ntxt = (Int_t) (n * 10.);
  Int_t Ttxt = (Int_t) (t * 10.);
  TString out(outDir);
  out += "/";
  out += oName;
  out += ".p";
  out += p;
  out += "m";
  out += m;
  out += "k";
  out += k;
  out += "n";
  out += Ntxt;
  out += "t";
  out += Ttxt;
  out += ".root";


  // unfolding
  Double_t chi2;
  StJetFolder f(m, p);
  f.Init(out, tIn, sIn, mIn);
  f.SetNamecycles(tName, sName, mName);
  f.SetJetParameters(nRM, rJet, aMin, pTmin, pTmax, eTmin, eTmax);
  f.SetUnfoldParameters(k, nMC, nToy, bPrior, n, t, mPrior);
  f.SetAbsoluteTruth(aIn, aName);
  chi2 = f.Unfold();
  f.Finish();

}

// End ------------------------------------------------------------------------
