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
const TString  pFile("input/pp200r9embed.unfoldingPtBin1.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  sFile("input/pp200r9embed.unfoldingPtBin1.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  mFile("input/pp200r9.forUnfoldingPtBin1.et920vz55.r03a02rm1chrg.root");
const TString  eFile("input/pp200r9embed.unfoldingPtBin1.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  rFile("input/pp200r9embed.unfoldingPtBin1.et9vz55.r03a02rm1chrg.dr03q15.root");
const TString  oFile("pp200r9.unfolding.et9vz55.r03a02rm1chrg.dr03q15");
// input namecycles
const TString  pName("hSumParAll");
const TString  sName("hSumDetAll");
const TString  mName("Pi0/hJetPtCorrP");
const TString  eName("hEfficiencyAll");
const TString  rName("hResponseAll");

// unfolding parameters (to loop over)
const Int_t nM  = 2;
const Int_t M[] = {1, 2};
const Int_t nK  = 30;
const Int_t K[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
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
  TString sBayU("hBayUnfold");
  TString sBayB("hBayBackfold");
  TString sSvdU("hSvdUnfold");
  TString sSvdB("hSvdBackfold");
  TString sChiX("k_{reg}");
  TString sChiYU("#chi^{2}(unfold, prior)");
  TString sChiYB("#chi^{2}(backfold, measured)");
  TString sChi2(oFile);
  sChi2 += ".performance.root";

  TFile *fChi2 = new TFile(sChi2.Data(), "recreate");
  TH1D  *hBayUnfold    = new TH1D(sBayU.Data(), "Bayesian Algorithm", nK, K[0], K[nK - 1]);
  TH1D  *hBayBackfold  = new TH1D(sBayB.Data(), "Bayesian Algorithm", nK, K[0], K[nK - 1]);
  TH1D  *hSvdUnfold    = new TH1D(sSvdU.Data(), "SVD Algorithm", nK, K[0], K[nK - 1]);
  TH1D  *hSvdBackfold  = new TH1D(sSvdB.Data(), "SVD Algorithm", nK, K[0], K[nK - 1]);
  hBayUnfold   -> GetXaxis() -> SetName(sChiX.Data());
  hBayUnfold   -> GetYaxis() -> SetName(sChiYU.Data());
  hBayBackfold -> GetXaxis() -> SetName(sChiX.Data());
  hBayBackfold -> GetYaxis() -> SetName(sChiYB.Data());
  hSvdUnfold   -> GetXaxis() -> SetName(sChiX.Data());
  hSvdUnfold   -> GetYaxis() -> SetName(sChiYU.Data());
  hSvdBackfold -> GetXaxis() -> SetName(sChiX.Data());
  hSvdBackfold -> GetYaxis() -> SetName(sChiYB.Data());
  hBayUnfold   -> Sumw2();
  hBayBackfold -> Sumw2();
  hSvdUnfold   -> Sumw2();
  hSvdBackfold -> Sumw2();


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
          Double_t    chi2u = 0.;
          Double_t    chi2b = 0.;
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
          f.Unfold(chi2u);
          f.Backfold(chi2b);
          f.Finish();

          const Double_t merit = TMath::Abs(chi2b - 1);
          const Double_t best  = TMath::Abs(chi2best - 1);
          if (merit < best) {
            chi2best = chi2b;
            bestFile = output;
          }

          // record chi2
          UInt_t iReg(0);
          switch (method) {
            case 1:
              iReg = hBayBackfold -> FindBin(kReg);
              hBayUnfold   -> SetBinContent(iReg, chi2u);
              hBayUnfold   -> SetBinError(iReg, 0.);
              hBayBackfold -> SetBinContent(iReg, chi2b);
              hBayBackfold -> SetBinError(iReg, 0.);
              break;
            case 2:
              iReg = hSvdBackfold -> FindBin(kReg);
              hSvdUnfold   -> SetBinContent(iReg, chi2u);
              hSvdUnfold   -> SetBinError(iReg, 0.);
              hSvdBackfold -> SetBinContent(iReg, chi2b);
              hSvdBackfold -> SetBinError(iReg, 0.);
              break;
          }

        }  // end k loop
      }  // end method loop
    }  // end tPrior loop
  }  // end nPrior loop


  // save chi2
  fChi2        -> cd();
  hBayUnfold   -> Write();
  hBayBackfold -> Write();
  hSvdUnfold   -> Write();
  hSvdBackfold -> Write();
  fChi2        -> Close();

  // announce winner
  TDatime end;
  cout << "\nFinished folding! " << end.AsString() << "\n"
       << "  Best chi2 = " << chi2best << "\n"
       << "  Best file = " << bestFile << "\n"
       << endl;

}

// End ------------------------------------------------------------------------
