// 'DoUnfolding.C'
// Derek Anderson
// 12.17.2016
// 
// This macro performs the unfolding / backfolding of a provided jet spectrum
// using the 'StJetFolder' class. Parameters:
// 
//   mIn    -- name of root-file containing measured spectrum
//   rIn    -- name of root-file containing the response matrix
//             and reconstruction efficiency
//   out    -- name of output root-file
//   mName  -- name of measured spectrum's histogram
//   rName  -- name of response matrix's histogram
//   eName  -- name of reconstruction effiency's profile
//   method -- which unfolding algorithm to be used:
//               0 = no unfolding
//               1 = bayesian
//               2 = SVD
//   kReg   -- regularization parameter (1~5 is usually
//             good, anything more will be too sensitive
//             to statistical fluctuations)
//   nMC    -- no. of iterations of backfolding (~100000 is 
//             usually sufficient)
//   prior  -- which type of prior to be used:
//               0 = pythia
//               1 = Levy
//               2 = Tsallis
//               3 = Exponential
//   nPrior -- adjusts how fast Levy fncn. drops off
//             (~5.8 typically produces good results)
//   tPrior -- adjusts slope of Levy fncn. (~0.4
//             typically produces good results)

#include <TSystem>
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TPaveText.h"

using namespace std;


class StJetFolder;


// input and output files
const TString  pFile("input/pp200r9embed.unfoldingVariableBins.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  sFile("input/pp200r9embed.unfoldingVariableBins.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  mFile("input/pp200r9embed.unfoldingVariableBins.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  eFile("input/pp200r9embed.unfoldingVariableBins.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  rFile("input/pp200r9embed.unfoldingVariableBins.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  oFile("pp200r9embed.closureTestVariableBinsRFFwFF.et9vz55.r03a02rm1chrg.dr03q15");
// input namecycles
const TString  pName("hSumParRFF");
const TString  sName("hSumDetRFF");
const TString  mName("hSumDetRFF");
const TString  eName("hEfficiencyFF");
const TString  rName("hResponseFF");

// unfolding parameters (to loop over)
const Int_t nM  = 2;
const Int_t M[] = {1, 2};
const Int_t nK  = 40;
const Int_t K[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
// prior parameters (to loop over)
const Int_t    P   = 0;
const Int_t    nN  = 1;
const Double_t N[] = {5.8};
const Double_t nT  = 1;
const Double_t T[] = {0.4};


// jet parameters (won't impact unfolding)
const Int_t    nRM     = 1;
const Int_t    type    = 0;
const Double_t rJet    = 0.3;
const Double_t aMin    = 0.2;
const Double_t pTmin   = 0.2;
const Double_t eTmin   = 9.;
const Double_t eTmax   = 20.;
const Double_t hTrgMax = 0.9;

// these don't need to be changed
const Int_t    beam   = 0;       // 0 = "pp", 1 = "AuAu"
const Int_t    trig   = 2;       // 0 = "gamma-dir", 1 = "gamma-rich", 2 = "pi0"
const Int_t    nToy   = 10;      // used to calculate covariances
const Int_t    nMC    = 100000;  // number of MC iterations for backfolding
const Double_t energy = 200.;    // sqrt(s)
const Double_t bPrior = 0.1;     // normalization of prior
const Double_t mPrior = 0.140;   // m-parameter of prior


void DoUnfolding() {

  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjet.so");
  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjettools.so");
  gSystem -> Load("../../RooUnfold/libRooUnfold.so");
  gSystem -> Load("StJetFolder");
  
  // lower verbosity
  gErrorIgnoreLevel = kError;


  TDatime start;
  cout << "\nStarting folding: " << start.AsString() << "\n" << endl;

  // for recording chi2
  TString sBay("hBayes");
  TString sSVD("hSVD");
  TString sChiX("k_{reg}");
  TString sChiY("#chi^{2}(backfold, measured)");
  TString sChi2(oFile);
  sChi2 += ".performance.root";

  TFile *fChi2 = new TFile(sChi2.Data(), "recreate");
  TH1D  *hBay  = new TH1D(sBay.Data(), "Bayesian Algorithm", nK, K[0], K[nK - 1]);
  TH1D  *hSVD  = new TH1D(sSVD.Data(), "SVD Algorithm", nK, K[0], K[nK - 1]);
  hBay -> GetXaxis() -> SetName(sChiX.Data());
  hSVD -> GetXaxis() -> SetName(sChiY.Data());
  hBay -> Sumw2();
  hSVD -> Sumw2();


  // loop over parameters
  Double_t chi2best = 999.;
  TString  bestFile;
  for (Int_t n = 0; n < nN; n++) {
    for (Int_t t = 0; t < nT; t++) {
      for (Int_t m = 0; m < nM; m++) {
        for (Int_t k = 0; k < nK; k++) {

          const Int_t    prior  = P;
          const Int_t    method = M[m];
          const Int_t    kReg   = K[k];
          const Double_t nPrior = N[n];
          const Double_t tPrior = T[t];

          // create output name
          Int_t   Ntxt = (Int_t) (n * 10.);
          Int_t   Ttxt = (Int_t) (t * 10.);
          TString output(oFile);
          output += ".p";
          output += prior;
          output += "m";
          output += method;
          output += "k";
          output += k;
          output += "n";
          output += Ntxt;
          output += "t";
          output += Ttxt;
          output += ".root";

          // create folder
          Double_t    chi2 = 0.;
          StJetFolder f(output.Data());
          // set spectra
          f.SetPrior(pFile.Data(), pName.Data());
          f.SetSmeared(sFile.Data(), sName.Data());
          f.SetMeasured(mFile.Data(), mName.Data());
          f.SetResponse(rFile.Data(), rName.Data());
          f.SetEfficiency(eFile.Data(), eName.Data());
          // set info and parameters
          f.SetEventInfo(beam, energy);
          f.SetTriggerInfo(trig, eTmin, eTmax, hTrgMax);
          f.SetJetInfo(type, nRM, rJet, aMin, pTmin);
          f.SetPriorParameters(prior, bPrior, mPrior, nPrior, tPrior);
          f.SetUnfoldParameters(method, kReg, nMC, nToy);
          // do unfolding
          f.Init();
          f.Unfold();
          f.Backfold(chi2);
          f.Finish();

          if (chi2 < chi2best) {
            chi2best = chi2;
            bestFile = output;
          }

          // record chi2
          UInt_t iReg(0);
          switch (method) {
            case 1:
              iReg = hBay -> FindBin(kReg);
              hBay -> SetBinContent(iReg, chi2);
              hBay -> SetBinError(iReg, 0.);
              break;
            case 2:
              iReg = hSVD -> FindBin(kReg);
              hSVD -> SetBinContent(iReg, chi2);
              hSVD -> SetBinError(iReg, 0.);
              break;
          }

        }  // end k loop
      }  // end method loop
    }  // end tPrior loop
  }  // end nPrior loop


  // save chi2
  fChi2 -> cd();
  hBay  -> Write();
  hSVD  -> Write();
  fChi2 -> Close();

  // announce winner
  TDatime end;
  cout << "\nFinished folding! " << end.AsString() << "\n"
       << "  Best chi2 = " << chi2best << "\n"
       << "  Best file = " << bestFile << "\n"
       << endl;

}

// End ------------------------------------------------------------------------
