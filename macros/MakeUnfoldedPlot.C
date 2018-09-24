// 'MakeUnfoldedPlot.C'
// Derek Anderson
// 07.16.2018
//
// Use this to plot an unfolded distribution
// vs. pythia (or otherwise) distribution.
//
// NOTE: the gamma-subtraction assumes that
//       the histogram specified by 'sInU'
//       is the gamma-rich distribution.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t NPad(2);
static const UInt_t NPlot(2);
static const UInt_t NHist(3);
static const UInt_t NRebin(2);
static const UInt_t NRebinVar(16);
static const Bool_t DoRebin(false);
static const Bool_t DoVariableRebin(false);
static const Bool_t DoGammaRichPlot(true);
static const Bool_t DoGammaSubtraction(false);



void MakeUnfoldedPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // io parameters
  const TString sOut("summedErrorsvsPythia.sysVsStat.et911vz55pi0.r02a005rm1chrg.d23m9y2018.root");
  const TString sInU("sysFinished/summedSystematics.et911vz55pi0.r02a005rm1chrg.d21m9y2018.root");
  const TString sInP("input/pp200py8.defaultResponse.pTbinRes.et911pi0.r02a005rm1chrg.dr02q015185.root");
  const TString sHistU("hSystematics");
  const TString sHistP("hParticle");

  // plot parameters
  const TString sTitle("");
  const TString sNameU("hSystematics");
  const TString sNameP("hPythia");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rhoA^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("unfold / particle");
  const TString sLabelU("unfolded [sys.]");
  const TString sLabelP("pythia");

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T}, R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // subtraction parameters
  const UInt_t   fFilGR(3017);
  const UInt_t   fLinGR(1);
  const UInt_t   fWidGR(1);
  const UInt_t   fMarGR(4);
  const UInt_t   fColGR(879);
  const UInt_t   fColP[NHist] = {859, 921, 859};
  const UInt_t   fColG[NHist] = {899, 921, 899};
  const TString  sFileGam("sysFinished/summedSystematics.et911vz55pi0.r02a005rm1chrg.d21m9y2018.root");
  const TString  sFilePi0("sysEt11/pp200r9.forSys.et1115vz55pi0.r02a005rm1chrg.p0m1k3n58t4.root");
  const TString  sHistGam("hStatistics");
  const TString  sHistPi0("hUnfolded");
  const TString  sNameGR("hStatistics");
  const TString  sLabelGR("unfolded [stat.]");
  const Double_t gammaPurity(0.52027 - 0.0356863);

  // misc parameters
  const Double_t plotRange[NPlot]    = {-1., 37.};
  const Double_t varRebin[NRebinVar] = {0., 1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 18., 21., 24., 27., 30.};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInU = new TFile(sInU.Data(), "read");
  TFile *fInP = new TFile(sInP.Data(), "read");
  if (!fOut || !fInU || !fInP) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInU = " << fInU << ", fInP" << fInP
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hUnfold = (TH1D*) fInU -> Get(sHistU.Data());
  TH1D *hPythia = (TH1D*) fInP -> Get(sHistP.Data());
  if (!hUnfold || !hPythia) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hUnfold = " << hUnfold << ", hPythia = " << hPythia
         << endl;
    return;
  }
  hUnfold -> SetName(sNameU.Data());
  hPythia -> SetName(sNameP.Data());
  cout << "    Grabbed histograms." << endl;


  // grab gamma-rich histogram (if need be)
  TH1D *hGamRich;
  if (DoGammaRichPlot) {
    TFile *fGam = new TFile(sFileGam.Data(), "read");
    if (!fGam) {
      cerr << "PANIC: couldn't open gamma-rich file!" << endl;
      return;
    }
    hGamRich = (TH1D*) fGam -> Get(sHistGam.Data());
    if (!hGamRich) {
      cerr << "PANIC: couldn't grab gamma-rich histogram!" << endl;
      return;
    }
    hGamRich -> SetName(sNameGR.Data());
    cout << "    Grabbed gamma-rich histogram." << endl;
  }


  // do gamma-subtraction (if need be)
  if (DoGammaSubtraction) {
    TFile *fPi0 = new TFile(sFilePi0.Data(), "read");
    if (!fPi0) {
      cerr << "PANIC: couldn't open pi0 file!" << endl;
      return;
    }
    TH1D  *hPi0 = (TH1D*) fPi0    -> Get(sHistPi0.Data());
    TH1D  *hGam = (TH1D*) hUnfold -> Clone();
    if (!hPi0 || !hGam) {
      cerr << "PANIC: couldn't grab pi0 or gamma histogram!\n"
           << "       hPi0 = " << hPi0 << ", hGam = " << hGam
           << endl;
      return;
    }
    hPi0    -> Scale(gammaPurity);
    hUnfold -> Add(hGam, hPi0, 1., -1.);
    hUnfold -> Scale(1. / (1. - gammaPurity));
    fPi0    -> Close();
    cout << "    Did gamma subtraction." << endl;
  }
  fOut -> cd();


  // rebin (if need be)
  hUnfold -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  hPythia -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  if (DoRebin) {
    hUnfold -> Rebin(NRebin);
    hPythia -> Rebin(NRebin);

    const UInt_t nBinU = hUnfold -> GetNbinsX();
    const UInt_t nBinP = hPythia -> GetNbinsX();
    for (UInt_t iBinU = 1; iBinU <( nBinU + 1); iBinU++) {
      const Double_t valU = hUnfold -> GetBinContent(iBinU);
      const Double_t errU = hUnfold -> GetBinError(iBinU);
      const Double_t sizU = hUnfold -> GetBinWidth(iBinU);
      hUnfold -> SetBinContent(iBinU, valU / sizU);
      hUnfold -> SetBinError(iBinU, errU / sizU);
    }
    for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
      const Double_t valP = hPythia -> GetBinContent(iBinP);
      const Double_t errP = hPythia -> GetBinError(iBinP);
      const Double_t sizP = hPythia -> GetBinWidth(iBinP);
      hPythia -> SetBinContent(iBinP, valP / sizP);
      hPythia -> SetBinError(iBinP, errP / sizP);
    }
    cout << "    Rebinned histograms (uniform)." << endl;
  }
  if (DoVariableRebin) {
    hUnfold = (TH1D*) hUnfold -> Rebin(NRebinVar - 1, sNameU.Data(), varRebin);
    hPythia = (TH1D*) hPythia -> Rebin(NRebinVar - 1, sNameP.Data(), varRebin);

    const UInt_t nBinU = hUnfold -> GetNbinsX();
    const UInt_t nBinP = hPythia -> GetNbinsX();
    for (UInt_t iBinU = 1; iBinU < (nBinU + 1); iBinU++) {
      const Double_t valU = hUnfold -> GetBinContent(iBinU);
      const Double_t errU = hUnfold -> GetBinError(iBinU);
      const Double_t sizU = hUnfold -> GetBinWidth(iBinU);
      hUnfold -> SetBinContent(iBinU, valU / sizU);
      hUnfold -> SetBinError(iBinU, errU / sizU);
    }
    for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
      const Double_t valP = hPythia -> GetBinContent(iBinP);
      const Double_t errP = hPythia -> GetBinError(iBinP);
      const Double_t sizP = hPythia -> GetBinWidth(iBinP);
      hPythia -> SetBinContent(iBinP, valP / sizP);
      hPythia -> SetBinError(iBinP, errP / sizP);
    }
    cout << "    Rebinned histograms (variable)." << endl;
  }

  // calculate ratio(s)
  TH1D *hRatio = (TH1D*) hUnfold -> Clone();
  hRatio -> SetName("hRatio");

  const UInt_t goodDiv = hRatio  -> Divide(hUnfold, hPythia, 1., 1.);
  const UInt_t nBinR   = hUnfold -> GetNbinsX();
  if (goodDiv == 0) {
    for (UInt_t iBinR = 0; iBinR < nBinR; iBinR++) {
      const UInt_t   iBinP = hPythia -> FindBin(hUnfold -> GetBinCenter(iBinR));
      const Double_t num   = hUnfold -> GetBinContent(iBinR);
      const Double_t den   = hPythia -> GetBinContent(iBinP);
      const Double_t nErr  = hUnfold -> GetBinError(iBinR);
      const Double_t dErr  = hUnfold -> GetBinError(iBinP);
      const Double_t nRel  = nErr / num;
      const Double_t dRel  = dErr / den;
      const Double_t ratio = num / den;
      const Double_t error = ratio * TMath::Sqrt((nRel * nRel) + (dRel * dRel));
      if (den > 0.) {
        hRatio -> SetBinContent(iBinR, ratio);
        hRatio -> SetBinError(iBinR, error);
      }
      else {
        hRatio -> SetBinContent(iBinR, 0.);
        hRatio -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
  }
  hRatio -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);


  TH1D *hRatioGR;
  if (DoGammaRichPlot) {
    hRatioGR = (TH1D*) hGamRich -> Clone();
    hRatioGR -> SetName("hRatioStat");

    const UInt_t goodDivGR = hRatioGR -> Divide(hGamRich, hPythia, 1., 1.);
    const UInt_t nBinGR    = hGamRich -> GetNbinsX();
    if (goodDivGR == 0) {
      for (UInt_t iBinGR = 0; iBinGR < nBinGR; iBinGR++) {
        const UInt_t   iBinP = hPythia  -> FindBin(hGamRich -> GetBinCenter(iBinGR));
        const Double_t num   = hGamRich -> GetBinContent(iBinGR);
        const Double_t den   = hPythia  -> GetBinContent(iBinP);
        const Double_t nErr  = hGamRich -> GetBinError(iBinGR);
        const Double_t dErr  = hGamRich -> GetBinError(iBinP);
        const Double_t nRel  = nErr / num;
        const Double_t dRel  = dErr / den;
        const Double_t ratio = num / den;
        const Double_t error = ratio * TMath::Sqrt((nRel * nRel) + (dRel * dRel));
        if (den > 0.) {
          hRatioGR -> SetBinContent(iBinGR, ratio);
          hRatioGR -> SetBinError(iBinGR, error);
        }
        else {
          hRatioGR -> SetBinContent(iBinGR, 0.);
          hRatioGR -> SetBinError(iBinGR, 0.);
        }
      }  // end bin loop
    }
    hRatioGR -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Calculated ratio(s)." << endl;


  // set styles
  UInt_t fCol[NHist];
  if (DoGammaSubtraction) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColG[iHist];
    }
  }
  else {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColP[iHist];
    }
  }

  const UInt_t  fMar[NHist] = {4, 1, 4};
  const UInt_t  fFil[NHist] = {3018, 0, 3018};
  const UInt_t  fLin[NHist] = {1, 1, 1};
  const UInt_t  fWid[NHist] = {1, 2, 1};
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hUnfold -> SetMarkerColor(fCol[0]);
  hUnfold -> SetMarkerStyle(fMar[0]);
  hUnfold -> SetFillColor(fCol[0]);
  hUnfold -> SetFillStyle(fFil[0]);
  hUnfold -> SetLineColor(fCol[0]);
  hUnfold -> SetLineStyle(fLin[0]);
  hUnfold -> SetLineWidth(fWid[0]);
  hUnfold -> SetTitle(sTitle.Data());
  hUnfold -> SetTitleFont(fTxt);
  hUnfold -> GetXaxis() -> SetTitle(sTitleX.Data());
  hUnfold -> GetXaxis() -> SetTitleFont(fTxt);
  hUnfold -> GetXaxis() -> SetTitleSize(fTit[1]);
  hUnfold -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hUnfold -> GetXaxis() -> SetLabelFont(fTxt);
  hUnfold -> GetXaxis() -> SetLabelSize(fLab[1]);
  hUnfold -> GetXaxis() -> CenterTitle(fCnt);
  hUnfold -> GetYaxis() -> SetTitle(sTitleY.Data());
  hUnfold -> GetYaxis() -> SetTitleFont(fTxt);
  hUnfold -> GetYaxis() -> SetTitleSize(fTit[1]);
  hUnfold -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hUnfold -> GetYaxis() -> SetLabelFont(fTxt);
  hUnfold -> GetYaxis() -> SetLabelSize(fLab[1]);
  hUnfold -> GetYaxis() -> CenterTitle(fCnt);
  hPythia -> SetMarkerColor(fCol[1]);
  hPythia -> SetMarkerStyle(fMar[1]);
  hPythia -> SetFillColor(fCol[1]);
  hPythia -> SetFillStyle(fFil[1]);
  hPythia -> SetLineColor(fCol[1]);
  hPythia -> SetLineStyle(fLin[1]);
  hPythia -> SetLineWidth(fWid[1]);
  hPythia -> SetTitle(sTitle.Data());
  hPythia -> SetTitleFont(fTxt);
  hPythia -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPythia -> GetXaxis() -> SetTitleFont(fTxt);
  hPythia -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPythia -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPythia -> GetXaxis() -> SetLabelFont(fTxt);
  hPythia -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPythia -> GetXaxis() -> CenterTitle(fCnt);
  hPythia -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPythia -> GetYaxis() -> SetTitleFont(fTxt);
  hPythia -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPythia -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPythia -> GetYaxis() -> SetLabelFont(fTxt);
  hPythia -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPythia -> GetYaxis() -> CenterTitle(fCnt);
  hRatio  -> SetMarkerColor(fCol[2]);
  hRatio  -> SetMarkerStyle(fMar[2]);
  hRatio  -> SetFillColor(fCol[2]);
  hRatio  -> SetFillStyle(fFil[2]);
  hRatio  -> SetLineColor(fCol[2]);
  hRatio  -> SetLineStyle(fLin[2]);
  hRatio  -> SetLineWidth(fWid[2]);
  hRatio  -> SetTitle(sTitle.Data());
  hRatio  -> SetTitleFont(fTxt);
  hRatio  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio  -> GetXaxis() -> SetTitleFont(fTxt);
  hRatio  -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatio  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatio  -> GetXaxis() -> SetLabelFont(fTxt);
  hRatio  -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatio  -> GetXaxis() -> CenterTitle(fCnt);
  hRatio  -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatio  -> GetYaxis() -> SetTitleFont(fTxt);
  hRatio  -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatio  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatio  -> GetYaxis() -> SetLabelFont(fTxt);
  hRatio  -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatio  -> GetYaxis() -> CenterTitle(fCnt);
  if (DoGammaRichPlot) {
    hGamRich  -> SetMarkerColor(fColGR);
    hGamRich  -> SetMarkerStyle(fMarGR);
    hGamRich  -> SetFillColor(fColGR);
    hGamRich  -> SetFillStyle(fFilGR);
    hGamRich  -> SetLineColor(fColGR);
    hGamRich  -> SetLineStyle(fLinGR);
    hGamRich  -> SetLineWidth(fWidGR);
    hGamRich  -> SetTitle(sTitle.Data());
    hGamRich  -> SetTitleFont(fTxt);
    hGamRich  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hGamRich  -> GetXaxis() -> SetTitleFont(fTxt);
    hGamRich  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hGamRich  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hGamRich  -> GetXaxis() -> SetLabelFont(fTxt);
    hGamRich  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hGamRich  -> GetXaxis() -> CenterTitle(fCnt);
    hGamRich  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hGamRich  -> GetYaxis() -> SetTitleFont(fTxt);
    hGamRich  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hGamRich  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hGamRich  -> GetYaxis() -> SetLabelFont(fTxt);
    hGamRich  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hGamRich  -> GetYaxis() -> CenterTitle(fCnt);
    hRatioGR  -> SetMarkerColor(fColGR);
    hRatioGR  -> SetMarkerStyle(fMarGR);
    hRatioGR  -> SetFillColor(fColGR);
    hRatioGR  -> SetFillStyle(fFilGR);
    hRatioGR  -> SetLineColor(fColGR);
    hRatioGR  -> SetLineStyle(fLinGR);
    hRatioGR  -> SetLineWidth(fWidGR);
    hRatioGR  -> SetTitle(sTitle.Data());
    hRatioGR  -> SetTitleFont(fTxt);
    hRatioGR  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioGR  -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioGR  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatioGR  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatioGR  -> GetXaxis() -> SetLabelFont(fTxt);
    hRatioGR  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatioGR  -> GetXaxis() -> CenterTitle(fCnt);
    hRatioGR  -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatioGR  -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioGR  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatioGR  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatioGR  -> GetYaxis() -> SetLabelFont(fTxt);
    hRatioGR  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatioGR  -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;


  // make legend(s)
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.3};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hPythia, sLabelP.Data());
  leg -> AddEntry(hUnfold, sLabelU.Data());
  if (DoGammaRichPlot)
    leg -> AddEntry(hGamRich, sLabelGR.Data());

  TLegend *legR = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  legR -> SetFillColor(fColLe);
  legR -> SetFillStyle(fFilLe);
  legR -> SetLineColor(fColLe);
  legR -> SetLineStyle(fLinLe);
  legR -> SetTextFont(fTxt);
  legR -> SetTextAlign(fAln);
  if (DoGammaRichPlot) {
    legR -> AddEntry(hRatio, sLabelU.Data());
    legR -> AddEntry(hRatioGR, sLabelGR.Data());
  }
  cout << "    Made legend(s)." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NPlot * NPlot] = {0.3, 0.1, 0.5, 0.3};
  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(2);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NPlot * NPlot] = {plotRange[0], 1., plotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;


  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.05);
  const Float_t fMarginT1(0.);
  const Float_t fMarginT2(0.05);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.);
  const Float_t fPadXY1[NPlot * NPlot] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NPlot * NPlot] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot   -> SetGrid(fGrid, fGrid);
  cPlot   -> SetTicks(fTick, fTick);
  cPlot   -> SetBorderMode(fMode);
  cPlot   -> SetBorderSize(fBord);
  pPad1   -> SetGrid(fGrid, fGrid);
  pPad1   -> SetTicks(fTick, fTick);
  pPad1   -> SetLogx(fLogX);
  pPad1   -> SetLogy(fLogY);
  pPad1   -> SetBorderMode(fMode);
  pPad1   -> SetBorderSize(fBord);
  pPad1   -> SetFrameBorderMode(fFrame);
  pPad1   -> SetLeftMargin(fMarginL);
  pPad1   -> SetRightMargin(fMarginR);
  pPad1   -> SetTopMargin(fMarginT1);
  pPad1   -> SetBottomMargin(fMarginB1);
  pPad2   -> SetGrid(fGrid, fGrid);
  pPad2   -> SetTicks(fTick, fTick);
  pPad2   -> SetLogx(fLogX);
  pPad2   -> SetLogy(fLogY);
  pPad2   -> SetBorderMode(fMode);
  pPad2   -> SetBorderSize(fBord);
  pPad2   -> SetFrameBorderMode(fFrame);
  pPad2   -> SetLeftMargin(fMarginL);
  pPad2   -> SetRightMargin(fMarginR);
  pPad2   -> SetTopMargin(fMarginT2);
  pPad2   -> SetBottomMargin(fMarginB2);
  cPlot   -> cd();
  pPad1   -> Draw();
  pPad2   -> Draw();
  pPad1   -> cd();
  hRatio  -> Draw("E5");
  if (DoGammaRichPlot) {
    hRatioGR -> Draw("SAME E5");
    legR     -> Draw();
  }
  line    -> Draw();
  pPad2   -> cd();
  hUnfold -> Draw("E5");
  if (DoGammaRichPlot)
    hGamRich -> Draw("SAME E5");
  hPythia -> Draw("SAME HIST C");
  leg     -> Draw();
  txt     -> Draw();
  fOut    -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOut    -> cd();
  hUnfold -> Write();
  hPythia -> Write();
  hRatio  -> Write();
  if (DoGammaRichPlot)
    hGamRich -> Write();
  fOut    -> Close();
  fInU    -> cd();
  fInU    -> Close();
  fInP    -> cd();
  fInP    -> Close();
  if (DoGammaRichPlot) {
    fGam  -> cd();
    fGam  -> Close();
  }
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
