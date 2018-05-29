// 'MakeMCinput.C'
// Derek Anderson
// 08.29.2016
//
// Use this to produce the input histograms for the 'MonteCarloMatchJets'
// function of the class 'StJetFolder'. If you don't do it beforehand,
// 'StJetFolder' will call this macro while running to produce the input.
//
// Last update: 08.30.2016

#include <TSystem>
#include <iostream>
#include "TString.h"

using namespace std;

class StMonteCarloInputMaker;


// min. trigger eT and max trigger eta
const Double_t eTmin = 9.;
const Double_t hMax  = 1.;
// input, output files
const TString iFile("Pythia19_GammaJet2.root");
const TString oFile("MonteCarloInput3.root");


void MakeMcInput() {

  gSystem -> Load("StMonteCarloInputMaker");

  // lower verbosity
  gErrorIgnoreLevel = kError;


  StMonteCarloInputMaker mc(iFile, oFile);
  mc.Make(eTmin, hMax);

}

// End ------------------------------------------------------------------------
