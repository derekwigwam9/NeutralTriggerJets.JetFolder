// 'SumErrors.C'
// Derek Anderson
// 09.21.2018
//
// Sum up some histograms and their
// errors why don'cha?  There should
// be NSys + 2 histograms: NSys
// histograms with systematic errors,
// 1 histogram with statistical errors,
// and the last one will be the sum
// of the two.


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;


// global constants
static const UInt_t NSys(2);
static const UInt_t NBin(23);



void SumErrors() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Summing errors..." << endl;

  // files
  const TString sOut("summedSystematics.et911pt0230vz55pi0.r02a005rm1chrg.d16m11y2018.root");
  const TString sStat("sysEt911r02/systematics.priorAndReg.et911vz55pt0230pi0.r02a005rm1chrg.d13m11y2018.root");
  const TString sSys[NSys] = {"sysEt911r02/systematics.priorAndReg.et911vz55pt0230pi0.r02a005rm1chrg.d13m11y2018.root", "sysEt911r02/systematics.effAndRes.et911vz55pt0230pi0.r02a005rm1chrg.d13m11y2018.root"};

  // input histograms
  const TString sHistStat("hDefault");
  const TString sHistSys[NSys] = {"hTotal", "hTotal"};


  // open files
  TFile *fOut  = new TFile(sOut.Data(), "recreate");
  TFile *fStat = new TFile(sStat.Data(), "read");
  if (!fStat) {
    cerr << "PANIC: couldn't open statistics file!" << endl;
    return;
  }

  TFile *fSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] = new TFile(sSys[iSys].Data(), "read");
    if (!fSys[iSys]) {
      cerr << "PANIC: couldn't systematic file " << iSys << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;


  // grab histograms
  TH1D *hStat = (TH1D*) fStat -> Get(sHistStat.Data());
  if (!hStat) {
    cerr << "PANIC: couldn't grab statistics histogram!" << endl;
    return;
  }

  TH1D *hSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] = (TH1D*) fSys[iSys] -> Get(sHistSys[iSys].Data());
    if (!hSys[iSys]) {
      cerr << "PANIC: couldn't grab systematic histogram " << iSys << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;


  // sum histograms
  TH1D *hStatOut = (TH1D*) hStat -> Clone();
  TH1D *hSysOut  = (TH1D*) hStat -> Clone();
  TH1D *hTotOut  = (TH1D*) hStat -> Clone();
  hStatOut -> SetName("hStatistics");
  hStatOut -> Reset("ICE");
  hSysOut  -> SetName("hSystematics");
  hSysOut  -> Reset("ICE");
  hTotOut  -> SetName("hTotal");
  hTotOut  -> Reset("ICE");

  Double_t valStat[NBin];
  Double_t errStat[NBin];
  Double_t perSys2[NBin];
  Double_t perTot2[NBin];

  const UInt_t iStart = hStat -> FindFirstBinAbove(0.);
  const UInt_t iLast  = hStat -> FindLastBinAbove(0.);
  const UInt_t iStop  = hStat -> GetNbinsX();
  for (UInt_t iBin = 0; iBin < NBin; iBin++) {
    const UInt_t iVal   = iBin + iStart;
    const Bool_t isLast = (iVal > iLast);
    const Bool_t isDone = (iVal > iStop);
    if (isLast || isDone) break;

    const Double_t val  = hStat -> GetBinContent(iVal);
    const Double_t stat = hStat -> GetBinError(iVal);
    const Double_t per  = stat / val;
    valStat[iBin] = val;
    errStat[iBin] = stat;
    perSys2[iBin] = 0.;
    perTot2[iBin] = per * per;
  }
  cout << "    Initialized errors." << endl;

  // loop over systematics
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    const UInt_t iCheck = hSys[iSys] -> FindFirstBinAbove(0.);
    if (iStart != iCheck) {
      cerr << "WARNING: looks like bin numbering is off in systematic no. " << iSys << "!\n"
           << "         iStart = " << iStart << ", iCheck = " << iCheck
           << endl;
    }
    for (UInt_t iBin = 0; iBin < NBin; iBin++) {
      const UInt_t iVal   = iBin + iStart;
      const Bool_t isLast = (iVal > iLast);
      const Bool_t isDone = (iVal > iStop);
      if (isLast || isDone) break;

      const Double_t valAdd = hSys[iSys] -> GetBinContent(iBin + iStart);
      const Double_t sysAdd = hSys[iSys] -> GetBinError(iBin + iStart);
      const Double_t perAdd = sysAdd / valAdd;
      if (valAdd > 0.) {
        perSys2[iBin] += perAdd * perAdd;
        perTot2[iBin] += perAdd * perAdd;
      }
    }  // end bin loop
  }  // end systematic loop
  cout << "    Summed errors." << endl;


  Double_t perSys[NBin];
  Double_t errSys[NBin];
  Double_t perTot[NBin];
  Double_t errTot[NBin];
  for (UInt_t iBin = 0; iBin < NBin; iBin++) {
    perSys[iBin] = TMath::Sqrt(perSys2[iBin]);
    perTot[iBin] = TMath::Sqrt(perTot2[iBin]);
    errSys[iBin] = perSys[iBin] * valStat[iBin];
    errTot[iBin] = perTot[iBin] * valStat[iBin];
    hStatOut -> SetBinContent(iBin + iStart, valStat[iBin]);
    hStatOut -> SetBinError(iBin + iStart, errStat[iBin]);
    hSysOut  -> SetBinContent(iBin + iStart, valStat[iBin]);
    hSysOut  -> SetBinError(iBin + iStart, errSys[iBin]);
    hTotOut  -> SetBinContent(iBin + iStart, valStat[iBin]);
    hTotOut  -> SetBinError(iBin + iStart, errTot[iBin]);
  }
  cout << "    Calculated errors." << endl;

  // write histograms
  fOut     -> cd();
  hStatOut -> Write();
  hSysOut  -> Write();
  hTotOut  -> Write();
  fOut     -> Close();
  fStat    -> cd();
  fStat    -> Close();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] -> cd();
    fSys[iSys] -> Close();
  }
  cout << "  DONE." << endl;

}

// End ------------------------------------------------------------------------
