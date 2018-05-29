// 'MergeData.C'
// Derek Anderson
// 10.24.2016
//
// This macro merges numerous TTree's into a TChain.

#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;


void MergeData() {

  const Int_t   nFiles = 20;
  const TString oName("input/Pythia23d.gMerged.root");
  const TString tName("DetTree");
  const TString cName("DetTree");
  const TString fPath("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Ana_nutralTgr_Jet/PythiaData/");
  const TString fPrefix("Pythia23.g");
  const TString fSuffix(".root");


  cout << "\n  Merging files..." << endl;

  TFile  *oFile = new TFile(oName, "recreate");
  TChain *chain = new TChain(tName, cName);
  for (Int_t i = 1; i <= nFiles; i++) {

    if (i == 2) continue;

    TString fName(fPath);
    fName += fPrefix;
    fName += i;
    fName += fSuffix;
    chain -> Add(fName);

    cout << "    File '" << fName << "' merged." << endl;

  }


  cout << "  Files merged!\n" << endl;

  // draw a few plots to check
  TCanvas *cTotalMult = new TCanvas("cTotalMult", "Total multiplicity", 200, 10, 700, 500);
  cTotalMult -> SetGrid(0, 0);
  chain      -> Draw("Events_refmult");
  cTotalMult -> Write();
  cTotalMult -> Close();

  TCanvas *cTrackPt = new TCanvas("cTrackPt", "Track pT", 200, 10, 700, 500);
  cTrackPt -> SetGrid(0, 0);
  cTrackPt -> SetLogy(1);
  chain    -> Draw("pTracks_pT");
  cTrackPt -> Write();
  cTrackPt -> Close(); 

  oFile -> cd();
  chain -> Write();
  oFile -> Close();

}
