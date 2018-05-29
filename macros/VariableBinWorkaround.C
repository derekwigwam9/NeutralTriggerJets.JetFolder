// 'VariableBinWorkaround.C'
// Derek Anderson
// 05.28.2018
//
// Use this to turn a variably-binned
// histogram into one with uniform
// bins by assigning a number to each
// bin or vice versa.

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"

using namespace std;


// global constants
static const UInt_t   NBinsX(37);
static const UInt_t   NBinsY(37);
static const Double_t BinsX[NBinsX] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 22.5, 25., 27.5, 30., 35., 40., 50., 60., 80.};
static const Double_t BinsY[NBinsY] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 22.5, 25., 27.5, 30., 35., 40., 50., 60., 80.};




void VariableBinWorkaround() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning variable bin workaround..." << endl;

  // io parameters
  const TString sIn("p200r9embed.closureTestALLwALL.PtBinVar.et9vz55.r03a02rm1chrg.dr03q15.p0m2k7n0t0.root");
  const TString sOut("PriorVsUnfoldRatio.afterUnfolding.d28m5y2018.root");
  const TString sInHist("hPriVsUnfoldRatio");
  const TString sInName("hPriVsUnfoldRatio_Uniform");
  const TString sOutHist("hPriVsUnfoldRatio_Variable");

  // misc parameters
  const UInt_t nBinsHistX(NBinsX - 1);
  const UInt_t nBinsHistY(NBinsY - 1);
  const Bool_t inputIsTwoDim(false);
  const Bool_t inputIsVariable(false);


  // open files
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open a file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;


  // grab histogram
  TH1D *hIn1D;
  TH2D *hIn2D;
  if (inputIsTwoDim) {
    hIn2D = (TH2D*) fIn -> Get(sInHist.Data());
    if (!hIn2D) {
      cerr << "PANIC: couldn't grab (2D) input histogram." << endl;
    }
    cout << "    Grabbed (2D) input histogram." << endl;
  }
  else {
    hIn1D = (TH1D*) fIn -> Get(sInHist.Data());
    if (!hIn1D) {
      cerr << "PANIC: couldn't grab (1D) input histogram." << endl;
    }
    cout << "    Grabbed (1D) input histogram." << endl;
  }


  // do workaround
  TH1D *hOut1D;
  TH2D *hOut2D;
  if (inputIsVariable) {
    if (inputIsTwoDim) {
      hOut2D = new TH2D(sOutHist.Data(), "", nBinsHistX, 0., nBinsHistX, nBinsHistY, 0., nBinsHistY);
      hOut2D -> Sumw2();
      for (UInt_t iBinX = 1; iBinX < NBinsX; iBinX++) {
        for (UInt_t iBinY = 1; iBinY < NBinsY; iBinY++) {
          const Double_t val = hIn2D -> GetBinContent(iBinX, iBinY);
          const Double_t err = hIn2D -> GetBinError(iBinX, iBinY);
          hOut2D -> SetBinContent(iBinX, iBinY, val);
          hOut2D -> SetBinError(iBinX, iBinY, err);
        }  // end y-bin loop
      }  // end x-bin loop
    }
    else {
      hOut1D = new TH1D(sOutHist.Data(), "", nBinsHistX, 0., nBinsHistX);
      hOut1D -> Sumw2();
      for (UInt_t iBinX = 1; iBinX < NBinsX; iBinX++) {
        const Double_t val = hIn1D -> GetBinContent(iBinX);
        const Double_t err = hIn1D -> GetBinError(iBinX);
        hOut1D -> SetBinContent(iBinX, val);
        hOut1D -> SetBinError(iBinX, err);
      }  // end x-bin loop
    }
  }  // end if(isVariable)
  else {
    if (inputIsTwoDim) {
      hOut2D = new TH2D(sOutHist.Data(), "", nBinsHistX, BinsX, nBinsHistY, BinsY);
      hOut2D -> Sumw2();
      for (UInt_t iBinX = 1; iBinX < NBinsX; iBinX++) {
        for (UInt_t iBinY = 1; iBinY < NBinsY; iBinY++) {
          const Double_t val = hIn2D -> GetBinContent(iBinX, iBinY);
          const Double_t err = hIn2D -> GetBinError(iBinX, iBinY);
          hOut2D -> SetBinContent(iBinX, iBinY, val);
          hOut2D -> SetBinError(iBinX, iBinY, err);
        }  // end y-bin loop
      }  // end x-bin loop
    }
    else {
      hOut1D = new TH1D(sOutHist.Data(), "", nBinsHistX, BinsX);
      hOut1D -> Sumw2();
      for (UInt_t iBinX = 1; iBinX < NBinsX; iBinX++) {
        const Double_t val = hIn1D -> GetBinContent(iBinX);
        const Double_t err = hIn1D -> GetBinError(iBinX);
        hOut1D -> SetBinContent(iBinX, val);
        hOut1D -> SetBinError(iBinX, err);
      }  // end x-bin loop
    }
  }
  cout << "    Did workaround." << endl;


  // write histograms
  fOut -> cd();
  if (inputIsTwoDim) {
    hIn2D  -> SetName(sInName.Data());
    hIn2D  -> Write();
    hOut2D -> Write();
  }
  else {
    hIn1D  -> SetName(sInName.Data());
    hIn1D  -> Write();
    hOut1D -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Workaround finished!\n" << endl;

}

// End ------------------------------------------------------------------------
