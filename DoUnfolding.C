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
//               3 = bin-by-bin
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
#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;


class StJetFolder;


// input and output files
static const TString pFile("input/pp200r9embed.pTbinRes.et920vz55.r02a005rm1chrg.dr02q015185.root");
static const TString sFile("input/pp200r9embed.pTbinRes.et920vz55.r02a005rm1chrg.dr02q015185.root");
static const TString mFile("input/pp200r9.pTbinRes.et911vz55.r02a005rm1chrg.d16m8y2018.root");
static const TString eFile("input/pp200py8.defaultResponse.pTbinRes.et920pi0.r02a005rm1chrg.dr02q015185.root");
static const TString rFile("input/pp200py8.defaultResponse.pTbinRes.et920pi0.r02a005rm1chrg.dr02q015185.root");
static const TString oFile("pp200r9.unfoldSystematicCheck_WithSystematics.et911vz55pi0.r02a005rm1chrg");
// input namecycles
static const TString pName("hSumParAll");
static const TString sName("hSumDetAll");
static const TString mName("Pi0/hJetPtCorrP");
static const TString eName("hEfficiency");
static const TString rName("hResponse");
// unfolding parameters (to loop over)
static const Int_t nM  = 1;
static const Int_t M[] = {1};
static const Int_t nK  = 1;
static const Int_t K[] = {2};
// prior parameters (to loop over)
static const Int_t    nP  = 1;
static const Int_t    nN  = 1;
static const Int_t    nT  = 1;
static const Int_t    P[] = {0};
static const Double_t N[] = {5.8};
static const Double_t T[] = {0.4};


// trigger and jet parameters (for plot labels)
static const Int_t   beam   = 0;        // 0 = "pp", 1 = "AuAu"
static const Int_t   trig   = 2;        // 0 = "gamma-dir", 1 = "gamma-rich", 2 = "pi0"
static const Int_t   type   = 0;        // 0 = "charged jets", 1 = "full jets"
static const Float_t energy = 200.;     // sqrt(s)
static const Float_t eTmin  = 9.;
static const Float_t eTmax  = 11.;
static const Float_t rJet   = 0.2;

// jet parameters (won't impact unfolding)
static const Int_t    nRM     = 1;
static const Double_t aMin    = 0.05;
static const Double_t pTmin   = 0.2;
static const Double_t pTmaxU  = 47.;
static const Double_t pTmaxB  = 38.;
static const Double_t hTrgMax = 0.9;

// these don't need to be changed
static const Int_t    nToy     = 10;      // used to calculate covariances
static const Int_t    nMC      = 100000;  // number of MC iterations for backfolding
static const Bool_t   smooth   = true;    // smooth efficiency at high pT
static const Bool_t   noErrors = true;    // remove errors on efficiency
static const Double_t bPrior   = 0.1;     // normalization of prior
static const Double_t mPrior   = 0.140;   // m-parameter of prior


void DoUnfolding() {

  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjet.so");
  gSystem -> Load("/common/star/star64/opt/star/sl64_gcc447/lib/libfastjettools.so");
  gSystem -> Load("../../RooUnfold/libRooUnfold.so");
  gSystem -> Load("StJetFolder");
  
  // lower verbosity
  gErrorIgnoreLevel = kError;

  TDatime start;
  cout << "\nStarting folding: " << start.AsString() << "\n" << endl;


  // create output stream
  const TString sStream(oFile.Data());
  sStream += ".bestFiles.list";

  ofstream bestFiles(sStream.Data());
  if (!bestFiles) {
    cerr << "PANIC: couldn't open output stream!" << endl;
    return;
  }


  // prior loops
  Double_t chi2bestest = 999.;
  TString  bestestFile;
  for (Int_t p = 0; p < nP; p++) {
    for (Int_t n = 0; n < nN; n++) {
      for (Int_t t = 0; t < nT; t++) {

        // don't double count priors...
        const Bool_t isPyth   = (p == 0);
        const Bool_t isExpo   = (p == 3);
        const Bool_t isFirstN = (n == 0);
        const Bool_t isFirstT = (t == 0);
        const Bool_t isFirst  = (isFirstN || isFirstT);
        if (isPyth && !isFirst) continue;
        if (isExpo && !isFirst) continue;

        // for file names
        const Double_t nPrior = N[n];
        const Double_t tPrior = T[t];
        const Float_t  Ntxt = N[n] * 10.;
        const Float_t  Ttxt = T[t] * 10.;

        // for recording chi2
        const TString sBayU("hBayUnfold");
        const TString sBayB("hBayBackfold");
        const TString sSvdU("hSvdUnfold");
        const TString sSvdB("hSvdBackfold");
        const TString sChiX("k_{reg}");
        const TString sChiYU("#chi^{2}(unfold, prior)");
        const TString sChiYB("#chi^{2}(backfold, measured)");

        // create performance file
        TString sChi2(oFile.Data());
        sChi2 += ".p";
        sChi2 += P[p];
        sChi2 += "n";
        sChi2 += Ntxt;
        sChi2 += "t";
        sChi2 += Ttxt;
        sChi2 += ".performance.root";

        TFile *fChi2 = new TFile(sChi2.Data(), "recreate");
        TH1D  *hBayUnfold    = new TH1D(sBayU.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hBayBackfold  = new TH1D(sBayB.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hSvdUnfold    = new TH1D(sSvdU.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hSvdBackfold  = new TH1D(sSvdB.Data(), "", nK, K[0], K[nK - 1] + 1);
        hBayUnfold   -> SetTitleFont(42);
        hBayUnfold   -> GetXaxis() -> SetTitle(sChiX.Data());
        hBayUnfold   -> GetXaxis() -> SetTitleOffset(1.);
        hBayUnfold   -> GetXaxis() -> SetTitleFont(42);
        hBayUnfold   -> GetXaxis() -> SetLabelFont(42);
        hBayUnfold   -> GetYaxis() -> SetTitle(sChiYU.Data());
        hBayUnfold   -> GetYaxis() -> SetTitleFont(42);
        hBayUnfold   -> GetYaxis() -> SetLabelFont(42);
        hBayBackfold -> SetTitleFont(42);
        hBayBackfold -> GetXaxis() -> SetTitle(sChiX.Data());
        hBayBackfold -> GetXaxis() -> SetTitleOffset(1.);
        hBayBackfold -> GetXaxis() -> SetTitleFont(42);
        hBayBackfold -> GetXaxis() -> SetLabelFont(42);
        hBayBackfold -> GetYaxis() -> SetTitle(sChiYB.Data());
        hBayBackfold -> GetYaxis() -> SetTitleFont(42);
        hBayBackfold -> GetYaxis() -> SetLabelFont(42);
        hSvdUnfold   -> SetTitleFont(42);
        hSvdUnfold   -> GetXaxis() -> SetTitle(sChiX.Data());
        hSvdUnfold   -> GetXaxis() -> SetTitleOffset(1.);
        hSvdUnfold   -> GetXaxis() -> SetTitleFont(42);
        hSvdUnfold   -> GetXaxis() -> SetLabelFont(42);
        hSvdUnfold   -> GetYaxis() -> SetTitle(sChiYU.Data());
        hSvdUnfold   -> GetYaxis() -> SetTitleFont(42);
        hSvdUnfold   -> GetYaxis() -> SetLabelFont(42);
        hSvdBackfold -> SetTitleFont(42);
        hSvdBackfold -> GetXaxis() -> SetTitle(sChiX.Data());
        hSvdBackfold -> GetXaxis() -> SetTitleOffset(1.);
        hSvdBackfold -> GetXaxis() -> SetTitleFont(42);
        hSvdBackfold -> GetXaxis() -> SetLabelFont(42);
        hSvdBackfold -> GetYaxis() -> SetTitle(sChiYB.Data());
        hSvdBackfold -> GetYaxis() -> SetTitleFont(42);
        hSvdBackfold -> GetYaxis() -> SetLabelFont(42);
        hBayUnfold   -> Sumw2();
        hBayBackfold -> Sumw2();
        hSvdUnfold   -> Sumw2();
        hSvdBackfold -> Sumw2();

        // method and k loop
        Double_t chi2best = 999.;
        TString  bestFile;
        for (Int_t m = 0; m < nM; m++) {
          for (Int_t k = 0; k < nK; k++) {

            const Int_t prior  = P[p];
            const Int_t method = M[m];
            const Int_t kReg   = K[k];

            // don't double count bin-by-bin corrections
            const Bool_t isBinByBin = (method == 3);
            const Bool_t isFirstK   = (k == 0);
            if (isBinByBin && !isFirstK) continue;

            // skip unreasonable kReg
            const Bool_t isBay      = (method == 1);
            const Bool_t isSVD      = (method == 2);
            const Bool_t isGoodBayK = (kReg < 6);
            const Bool_t isGoodSvdK = ((kReg > 5) && (kReg < 12));
            //if (isBay && !isGoodBayK) continue;
            //if (isSVD && !isGoodSvdK) continue;

            // create output name
            TString output(oFile);
            output += ".p";
            output += prior;
            output += "m";
            output += method;
            output += "k";
            output += kReg;
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
            f.SetEfficiency(eFile.Data(), eName.Data(), smooth, noErrors);
            // set info and parameters
            f.SetEventInfo(beam, energy);
            f.SetTriggerInfo(trig, eTmin, eTmax, hTrgMax);
            f.SetJetInfo(type, nRM, rJet, aMin, pTmin);
            f.SetPriorParameters(prior, bPrior, mPrior, nPrior, tPrior);
            f.SetUnfoldParameters(method, kReg, nMC, nToy, pTmaxU, pTmaxB);
            // do unfolding
            f.Init();
            f.Unfold(chi2u);
            f.Backfold(chi2b);
            f.Finish();

            const Double_t merit   = TMath::Abs(chi2b - 1);
            const Double_t best    = TMath::Abs(chi2best - 1);
            const Double_t bestest = TMath::Abs(chi2bestest - 1);
            if (merit < best) {
              chi2best = chi2b;
              bestFile = output;
            }
            if (merit < bestest) {
              chi2bestest = chi2b;
              bestestFile = output;
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
              default:
                break;
            }

          }  // end k loop
        }  // end method loop


        // make performance plots
        const UInt_t cBay(810);
        const UInt_t cSvd(860);
        hBayUnfold   -> SetLineColor(cBay);
        hBayUnfold   -> SetMarkerColor(cBay);
        hBayBackfold -> SetLineColor(cBay);
        hBayBackfold -> SetMarkerColor(cBay);
        hSvdUnfold   -> SetLineColor(cSvd);
        hSvdUnfold   -> SetMarkerColor(cSvd);
        hSvdBackfold -> SetLineColor(cSvd);
        hSvdBackfold -> SetMarkerColor(cSvd);

        TLegend *lUnfold   = new TLegend(0.1, 0.1, 0.3, 0.3);
        TLegend *lBackfold = new TLegend(0.1, 0.1, 0.3, 0.3);
        lUnfold   -> SetFillColor(0);
        lUnfold   -> SetLineColor(0);
        lUnfold   -> SetTextFont(42);
        lUnfold   -> SetTextAlign(12);
        lUnfold   -> AddEntry(hBayUnfold, "Bayes.");
        lUnfold   -> AddEntry(hSvdUnfold, "SVD");
        lBackfold -> SetFillColor(0);
        lBackfold -> SetLineColor(0);
        lBackfold -> SetTextFont(42);
        lBackfold -> SetTextAlign(12);
        lBackfold -> AddEntry(hBayBackfold, "Bayes.");
        lBackfold -> AddEntry(hSvdBackfold, "SVD");

        TLine *lOne = new TLine(K[0], 1, K[nK - 1] + 1, 1);
        lOne  -> SetLineColor(1);
        lOne  -> SetLineStyle(2);
        fChi2 -> cd();

        TCanvas *cUnfold   = new TCanvas("cUnfold", "", 750, 500);
        TCanvas *cBackfold = new TCanvas("cBackfold", "", 750, 500);
        cUnfold      -> SetGrid(0, 0);
        cBackfold    -> SetGrid(0, 0);
        cUnfold      -> cd();
        hBayUnfold   -> Draw();
        hSvdUnfold   -> Draw("same");
        lUnfold      -> Draw();
        lOne         -> Draw();
        cBackfold    -> cd();
        hBayBackfold -> Draw();
        hSvdBackfold -> Draw("same");
        lBackfold    -> Draw();
        lOne         -> Draw();
        cUnfold      -> Write();
        cUnfold      -> Close();
        cBackfold    -> Write();
        cBackfold    -> Close();

        // save chi2
        fChi2        -> cd();
        hBayUnfold   -> Write();
        hBayBackfold -> Write();
        hSvdUnfold   -> Write();
        hSvdBackfold -> Write();
        fChi2        -> Close();

        // announce winner
        TDatime endPrior;
        cout << "\nFinished folding prior! " << endPrior.AsString() << "\n"
             << "  Best chi2 = " << chi2best << "\n"
             << "  Best file = " << bestFile << "\n"
             << endl;

        // stream winner
        bestFiles << bestFile.Data();
        bestFiles << endl;

      }  // end tPrior loop
    }  // end nPrior loop
  }  // end prior loop


  // announce biggest winner
  TDatime end;
  cout << "\nFinished all folding! " << end.AsString() << "\n"
       << "  Bestest chi2 = " << chi2bestest << "\n"
       << "  Bestest file = " << bestestFile << "\n"
       << endl;
}

// End ------------------------------------------------------------------------
