// 'LocateCorrectDataPoints.C'
// Derek Anderson
// 09.26.2018
//
// Use this to determine the
// correct location of data
// points in wide bins.


#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TGraphAsymmErrors.h"

using namespace std;

// global constants
const UInt_t NIter(1);
const UInt_t NFunc(3);
const UInt_t NParam(2);



void LocateCorrectDataPoints() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Locating data points..." << endl;

  // input and graph constants
  const TString sIn("sysFinished/summedSystematics.et911pt0215vz55pi0.r02a005rm1chrg.d26m9y2018.root");
  const TString sOut("shiftedDataPoints.nIter1.et911vz55pi0.r02a005rm1chrg.d27m9y2018.root");
  const TString sTxt("shiftedDataPoints.nIter1.et911vz55pi0.r02a005rm1chrg.d27m9y2018.txt");
  const TString sHist("hStatistics");
  const TString sGraph[NIter + 1] = {"gIter0", "gIter1"};
  const UInt_t  cGraph[NIter + 1] = {858, 898};
  const UInt_t  mGraph[NIter + 1] = {24, 25};
  const UInt_t  fGraph[NIter + 1] = {0, 0};

  // function constants
  const Float_t xStart[NFunc] = {0.2, 3., 23.};
  const Float_t xStop[NFunc]  = {3., 23., 47.};
  const Float_t guess[NParam] = {-3., 3.};
  const TString sName[NFunc]  = {"fExp023", "fExp323", "fExp2347"};
  const TString sFunc[NFunc]  = {"expo", "expo", "expo"};


  // open files and histograms
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  TH1D *hIn = (TH1D*) fIn -> Get(sHist.Data());
  if (!hIn) {
    cerr << "PANIC: couldn't grab input histogram!" << endl;
    return;
  }
  cout << "    Grabbed histograms." << endl;


  // make initial graph
  TGraphAsymmErrors *gIn = new TGraphAsymmErrors(hIn);
  gIn -> SetName("gInput");
  cout << "    Made initial graph." << endl;

  // make functions
  TF1 *fFit[NFunc];
  for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
    fFit[iFunc] = new TF1(sName[iFunc].Data(), sFunc[iFunc].Data(), xStart[iFunc], xStop[iFunc]);
    for (UInt_t iPar = 0; iPar < NParam; iPar++) {
      fFit[iFunc] -> SetParameter(iPar, guess[iPar]);
    }
  }
  const UInt_t nPts = gIn -> GetN();
  cout << "    Made functions.\n"
       << "    Calculating new points..."
       << endl;


  // determine correct points
  TVectorD xValNew(nPts);
  TVectorD yValNew(nPts);
  TVectorD eXloNew(nPts);
  TVectorD eXhiNew(nPts);
  TVectorD eYloNew(nPts);
  TVectorD eYhiNew(nPts);

  TGraph *gIter[NIter + 1];
  gIter[0] = (TGraphAsymmErrors*) gIn -> Clone();
  gIter[0] -> SetName(sGraph[0].Data());
  gIter[0] -> SetLineColor(cGraph[0]);
  gIter[0] -> SetFillColor(cGraph[0]);
  gIter[0] -> SetFillStyle(fGraph[0]);
  gIter[0] -> SetMarkerColor(cGraph[0]);
  gIter[0] -> SetMarkerStyle(mGraph[0]);
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {

    // fit graph
    cout << "      Iteration " << iIter << "..." << endl;
    for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
      gIter[iIter] -> Fit(sName[iFunc].Data(), "+RN0");
    }

    // new point calculation
    for (UInt_t iPt = 0; iPt < nPts; iPt++) {
      Double_t x(0.);
      Double_t y(0.);
      gIter[iIter] -> GetPoint(iPt, x, y);

      // determine relevant fit
      UInt_t iFit(0);
      Bool_t isFitted(false);
      for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
        const Bool_t isInFitRange = ((x >= xStart[iFunc]) && (x < xStop[iFunc]));
        if (isInFitRange) {
          iFit     = iFunc;
          isFitted = true;
          break;
        }
      }

      // produce means
      if (isFitted) {
        const Double_t xOld = x;
        const Double_t yOld = y;
        const Double_t eXlo = gIter[iIter] -> GetErrorXlow(iPt);
        const Double_t eXhi = gIter[iIter] -> GetErrorXhigh(iPt);
        const Double_t eYlo = gIter[iIter] -> GetErrorYlow(iPt);
        const Double_t eYhi = gIter[iIter] -> GetErrorYhigh(iPt);
        const Double_t xLo  = xOld - eXlo;
        const Double_t xHi  = xOld + eXhi;
        const Double_t mean = fFit[iFit] -> Mean(xLo, xHi);
        const Double_t diff = TMath::Abs(mean - xOld);

        // set new points
        const Double_t xNew   = mean;
        const Double_t eLoNew = xNew - xLo;
        const Double_t eHiNew = xHi - xNew;
        xValNew[iPt] = xNew;
        yValNew[iPt] = yOld;
        eXloNew[iPt] = eLoNew;
        eXhiNew[iPt] = eHiNew;
        eYloNew[iPt] = eYlo;
        eYhiNew[iPt] = eYhi;
      }
      else {
        xValNew[iPt] = x;
        yValNew[iPt] = y;
        eXloNew[iPt] = gIter[iIter] -> GetErrorXlow(iPt);
        eXhiNew[iPt] = gIter[iIter] -> GetErrorXhigh(iPt);
        eYloNew[iPt] = gIter[iIter] -> GetErrorYlow(iPt);
        eYhiNew[iPt] = gIter[iIter] -> GetErrorYhigh(iPt);
      }
    }  // end point loop

    // create new graph
    gIter[iIter + 1] = new TGraphAsymmErrors(xValNew, yValNew, eXloNew, eXhiNew, eYloNew, eYhiNew);
    gIter[iIter + 1] -> SetName(sGraph[iIter + 1].Data());
    gIter[iIter + 1] -> SetLineColor(cGraph[iIter + 1]);
    gIter[iIter + 1] -> SetFillColor(cGraph[iIter + 1]);
    gIter[iIter + 1] -> SetFillStyle(fGraph[iIter + 1]);
    gIter[iIter + 1] -> SetMarkerColor(cGraph[iIter + 1]);
    gIter[iIter + 1] -> SetMarkerStyle(mGraph[iIter + 1]);

    // reset vectors
    xValNew.Zero();
    yValNew.Zero();
    eXloNew.Zero();
    eXhiNew.Zero();
    eYloNew.Zero();
    eYhiNew.Zero();

  }
  cout << "    Location calculation finished." << endl;


  // make legend
  const UInt_t  cLeg(0);
  const UInt_t  txt(42);
  const UInt_t  aln(12);
  const Float_t aLeg(0.7);
  const Float_t bLeg(0.9);
  const TString sLeg("iteration");

  TLegend *leg = new TLegend(aLeg, aLeg, bLeg, bLeg);
  leg -> SetFillColor(cLeg);
  leg -> SetLineColor(cLeg);
  leg -> SetTextFont(txt);
  leg -> SetTextAlign(aln);
  for (UInt_t iIter = 0; iIter < (NIter + 1); iIter++) {
    TString sLabel(sLeg.Data());
    sLabel += " ";
    sLabel += iIter;
    leg -> AddEntry(gIter[iIter], sLabel.Data());
  }
  cout << "    Made legend." << endl;


  // make canvas
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  logy(1);
  const TString sCan("cCalc");

  TCanvas *cCalc = new TCanvas(sCan.Data(), "", width, height);
  cCalc    -> SetGrid(grid, grid);
  cCalc    -> SetLogy(logy);
  fOut     -> cd();
  cCalc    -> cd();
  gIter[0] -> Draw("ALP");
  for (UInt_t iIter = 1; iIter < (NIter + 1); iIter++) {
    gIter[iIter] -> Draw("SAME LP");
  }
  leg   -> Draw();
  cCalc -> Write();
  cCalc -> Close();
  cout << "    Made plot." << endl;


  // write new points to text file
  ofstream oNew(sTxt.Data());
  for (UInt_t iPt = 0; iPt < nPts; iPt++) {
    Double_t x(0.);
    Double_t y(0.);
    gIter[NIter] -> GetPoint(iPt, x, y);
    if (iPt == 0) oNew << "x y xErrLo xErrHi yErrLo yErrHi" << endl;

    oNew << x;
    oNew << " ";
    oNew << y;
    oNew << " ";
    oNew << gIter[NIter] -> GetErrorXlow(iPt);
    oNew << " ";
    oNew << gIter[NIter] -> GetErrorXhigh(iPt);
    oNew << " ";
    oNew << gIter[NIter] -> GetErrorYlow(iPt);
    oNew << " ";
    oNew << gIter[NIter] -> GetErrorYhigh(iPt);
    oNew << endl;
  }
  cout << "    Streamed to output file." << endl;


  // write output and close files
  fOut -> cd();
  gIn  -> Write();
  for (UInt_t iIter = 0; iIter < (NIter + 1); iIter++) {
    gIter[iIter] -> Write();
  }
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Located data points.\n" << endl;

}

// End ------------------------------------------------------------------------
