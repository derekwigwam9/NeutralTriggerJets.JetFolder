// 'CalculateSystematicError.C'
// Derek Anderson
// 07.17.2018
//
// Use this to calculate a systematic
// error wrt a defualt value.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TColor.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t NPad(2);
static const UInt_t NSys(2);
static const UInt_t NPlot(2);
static const UInt_t NRebin(2);
static const UInt_t NRebinVar(16);
static const Bool_t DoRebin(false);
static const Bool_t UseAverage(true);
static const Bool_t DoVariableRebin(false);
static const Bool_t DoGammaSubtraction(false);



void CalculateSystematicError() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // io parameters
  const TString sOut("test.root");
  //const TString sOut("systematics.varyingTrkEff.et911vz55gd.r02a005rm1chrg.d2m9y2018.root");
  const TString sInD("output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffDefault.et911vz55gam.r02a005rm1chrg.p0m2k11n58t4.root");
  const TString sHistD("hUnfolded");
  const TString sInS[NSys]   = {"output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffM5.et911vz55gam.r02a005rm1chrg.p0m2k9n58t4.root", "output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffP5.et911vz55gam.r02a005rm1chrg.p0m2k10n58t4.root"};
  const TString sHistS[NSys] = {"hUnfolded", "hUnfolded"};

  // plot parameters
  const TString sTitle("");
  const TString sNameD("hDefault");
  const TString sNameA("hAverage");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rhoA^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleRD("1 + [(var. - def.) / def.]");
  const TString sTitleRA("1 + [(var. - avg.) / avg.]");
  const TString sLabelD("default #epsilon_{trk}");
  const TString sLabelA("average");
  const TString sNameS[NSys]  = {"hSys_TrkEffMinus5", "hSys_TrkEffPlus5"};
  const TString sNameR[NSys]  = {"hDif_TrkEffMinus5", "hDif_TrkEffPlus5"};
  const TString sLabelS[NSys] = {"#epsilon_{trk} - 5%", "#epsilon_{trk} + 5%"};

  // subtraction parameters
  const TString  sDefPi0("output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffDefault.et911vz55pi0.r02a005rm1chrg.p0m2k11n58t4.root");
  const TString  sHistDefPI0("hUnfolded");
  const TString  sInPi0[NSys]   = {"output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffM5.et911vz55pi0.r02a005rm1chrg.p0m2k9n58t4.root", "output/AugustUnfolding2018/TrkEffSysR02/pp200r9.trkEffP5.et911vz55pi0.r02a005rm1chrg.p0m2k10n58t4.root"};
  const TString  sHistPi0[NSys] = {"hUnfolded", "hUnfolded"};
  const Double_t gammaPurity(0.569875);

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("#gamma^{dir} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T}, R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // misc parameters
  const UInt_t  fColS[NSys]          = {806, 886};
  const UInt_t  fFilS[NSys]          = {3345, 3354};
  const Double_t plotRange[NPlot]    = {-1., 37.};
  const Double_t varRebin[NRebinVar] = {0., 1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 18., 21., 24., 27., 30.};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInD = new TFile(sInD.Data(), "read");
  if (!fOut || !fInD) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInD = " << fInD
         << endl;
    return;
  }

  TFile *fInS[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fInS[iSys] = new TFile(sInS[iSys].Data(), "read");
    if (!fInS[iSys]) {
      cerr << "PANIC: couldn't open a file!\n"
           << "       fInS[" << iSys <<"] = " << fInS[iSys]
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hDefault = (TH1D*) fInD -> Get(sHistD.Data());
  if (!hDefault) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hDefault = " << hDefault
         << endl;
    return;
  }
  hDefault -> SetName(sNameD.Data());
  hDefault -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);

  TH1D *hSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] = (TH1D*) fInS[iSys] -> Get(sHistS[iSys].Data());
    if (!hSys[iSys]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hSys[" << iSys << "] = " << hSys[iSys]
           << endl;
      return;
    }
    hSys[iSys] -> SetName(sNameS[iSys].Data());
    hSys[iSys] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Grabbed histograms." << endl;


  // do gamma-subtraction (if need be)
  TFile *fPi0[NSys];
  TH1D  *hPi0[NSys];
  TH1D  *hGam[NSys];
  if (DoGammaSubtraction) {
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      fPi0[iSys] = new TFile(sInPi0[iSys].Data(), "read");
      if (!fPi0[iSys]) {
        cerr << "PANIC: couldn't open a pi0 file!\n"
             << "       fInPi0[" << iSys << "] = " << fInPi0[iSys]
             << endl;
        return;
      }
      hPi0[iSys] = (TH1D*) fPi0[iSys] -> Get(sHistPi0[iSys].Data());
      hGam[iSys] = (TH1D*) hSys[iSys] -> Clone();
      if (!hPi0[iSys] || !hGam[iSys]) {
        cerr << "PANIC: couldn't grab pi0 or gamma histogram!\n"
             << "       hPi0[" << iSys << "] = " << hPi0[iSys] << ", hGam[" << iSys << "] = " << hGam[iSys]
             << endl;
        return;
      }
      hPi0[iSys] -> Scale(1. - gammaPurity);
      hSys[iSys] -> Add(hGam[iSys], hPi0[iSys], 1., -1.);
      hSys[iSys] -> Scale(1. / gammaPurity);
      fPi0[iSys] -> Close();
    }  // end systematic loop
    cout << "    Did gamma subtraction." << endl;
  }


  // rebin (if need be)
  if (DoRebin) {
    // rebin default histogram
    hDefault -> Rebin(NRebin);

    const UInt_t nBinD = hDefault -> GetNbinsX();
    for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
      const Double_t valD = hDefault -> GetBinContent(iBinD);
      const Double_t errD = hDefault -> GetBinError(iBinD);
      const Double_t sizD = hDefault -> GetBinWidth(iBinD);
      hDefault -> SetBinContent(iBinD, valD / sizD);
      hDefault -> SetBinError(iBinD, errD / sizD);
    }

    // rebin systematic histograms
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      hSys[iSys] -> Rebin(NRebin);

      const UInt_t nBinS = hSys[iSys] -> GetNbinsX();
      for (UInt_t iBinS = 1; iBinS < (nBinS + 1); iBinS++) {
        const Double_t valS = hSys[iSys] -> GetBinContent(iBinS);
        const Double_t errS = hSys[iSys] -> GetBinError(iBinS);
        const Double_t sizS = hSys[iSys] -> GetBinWidth(iBinS);
        hSys[iSys] -> SetBinContent(iBinS, valS / sizS);
        hSys[iSys] -> SetBinError(iBinS, errS / sizS);
      }
    }
    cout << "    Rebinned histograms (uniform)." << endl;
  }
  if (DoVariableRebin) {
    // rebin default histogram
    hDefault = (TH1D*) hDefault -> Rebin(NRebinVar - 1, sNameD.Data(), varRebin);

    const UInt_t nBinD = hDefault -> GetNbinsX();
    for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
      const Double_t valD = hDefault -> GetBinContent(iBinD);
      const Double_t errD = hDefault -> GetBinError(iBinD);
      const Double_t sizD = hDefault -> GetBinWidth(iBinD);
      hDefault -> SetBinContent(iBinD, valD / sizD);
      hDefault -> SetBinError(iBinD, errD / sizD);
    }

    // rebin systematic histograms
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      hSys[iSys] = (TH1D*) hSys[iSys] -> Rebin(NRebinVar - 1, sNameS[iSys].Data(), varRebin);

      const UInt_t nBinS = hSys[iSys] -> GetNbinsX();
      for (UInt_t iBinS = 1; iBinS < (nBinS + 1); iBinS++) {
        const Double_t valS = hSys[iSys] -> GetBinContent(iBinS);
        const Double_t errS = hSys[iSys] -> GetBinError(iBinS);
        const Double_t sizS = hSys[iSys] -> GetBinWidth(iBinS);
        hSys[iSys] -> SetBinContent(iBinS, valS / sizS);
        hSys[iSys] -> SetBinError(iBinS, errS / sizS);
      }
    }
    cout << "    Rebinned histograms (variable)." << endl;
  }


  // calculate average
  TH1D *hAverage = (TH1D*) hDefault -> Clone();
  hAverage -> SetName(sNameA.Data());
  hAverage -> Reset("ICE");
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hAverage -> Add(hSys[iSys], 1.);
  }
  hAverage -> Scale(1. / (Double_t) NSys);
  cout << "    Calculated average." << endl;


  // calculate %-difference
  TH1D *hDif[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hDif[iSys] = (TH1D*) hDefault -> Clone();
    hDif[iSys] -> SetName(sNameR[iSys].Data());
    hDif[iSys] -> Reset("ICE");

    const UInt_t nBinR = hDif[iSys] -> GetNbinsX();
    for (UInt_t iBinR = 1; iBinR < (nBinR + 1); iBinR++) {
      const Double_t def  = hDefault   -> GetBinContent(iBinR);
      const Double_t avg  = hAverage   -> GetBinContent(iBinR);
      const Double_t sys  = hSys[iSys] -> GetBinContent(iBinR);
      const Double_t raw  = hSys[iSys] -> GetBinError(iBinR);
      const Double_t rel  = (raw / sys);
      const Double_t difD = (sys - def) / def;
      const Double_t difA = (sys - avg) / avg;

      Double_t dif = 0.;
      if (UseAverage)
        dif = difA;
      else
        dif = difD;

      const Double_t val = 1. + dif;
      const Double_t err = val * rel;
      if ((def > 0.) && (sys > 0.)) {
        hDif[iSys] -> SetBinContent(iBinR, val);
        hDif[iSys] -> SetBinError(iBinR, err);
      }
      else {
        hDif[iSys] -> SetBinContent(iBinR, 0.);
        hDif[iSys] -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
    hDif[iSys] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Calculated difference." << endl;


  // set styles
  TString sTitleR("");
  if (UseAverage)
    sTitleR = sTitleRA;
  else
    sTitleR = sTitleRD;

  const UInt_t  fColD(1);
  const UInt_t  fColA(1);
  const UInt_t  fMarD(33);
  const UInt_t  fMarA(8);
  const UInt_t  fMarS(4);
  const UInt_t  fFilD(0);
  const UInt_t  fFilA(0);
  const UInt_t  fLinD(1);
  const UInt_t  fLinA(1);
  const UInt_t  fLinS(1);
  const UInt_t  fWidD(1);
  const UInt_t  fWidA(1);
  const UInt_t  fWidS(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fAlphaD(0.);
  const Float_t fAlphaS(0.5);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hDefault -> SetMarkerColor(fColD);
  hDefault -> SetMarkerStyle(fMarD);
  hDefault -> SetFillColor(fColD);
  hDefault -> SetFillStyle(fFilD);
  hDefault -> SetLineColor(fColD);
  hDefault -> SetLineStyle(fLinD);
  hDefault -> SetLineWidth(fWidD);
  hDefault -> SetTitle(sTitle.Data());
  hDefault -> SetTitleFont(fTxt);
  hDefault -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDefault -> GetXaxis() -> SetTitleFont(fTxt);
  hDefault -> GetXaxis() -> SetTitleSize(fTit[1]);
  hDefault -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hDefault -> GetXaxis() -> SetLabelFont(fTxt);
  hDefault -> GetXaxis() -> SetLabelSize(fLab[1]);
  hDefault -> GetXaxis() -> CenterTitle(fCnt);
  hDefault -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDefault -> GetYaxis() -> SetTitleFont(fTxt);
  hDefault -> GetYaxis() -> SetTitleSize(fTit[1]);
  hDefault -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hDefault -> GetYaxis() -> SetLabelFont(fTxt);
  hDefault -> GetYaxis() -> SetLabelSize(fLab[1]);
  hDefault -> GetYaxis() -> CenterTitle(fCnt);
  hAverage -> SetMarkerColor(fColA);
  hAverage -> SetMarkerStyle(fMarA);
  hAverage -> SetFillColor(fColA);
  hAverage -> SetFillStyle(fFilA);
  hAverage -> SetLineColor(fColA);
  hAverage -> SetLineStyle(fLinA);
  hAverage -> SetLineWidth(fWidA);
  hAverage -> SetTitle(sTitle.Data());
  hAverage -> SetTitleFont(fTxt);
  hAverage -> GetXaxis() -> SetTitle(sTitleX.Data());
  hAverage -> GetXaxis() -> SetTitleFont(fTxt);
  hAverage -> GetXaxis() -> SetTitleSize(fTit[1]);
  hAverage -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hAverage -> GetXaxis() -> SetLabelFont(fTxt);
  hAverage -> GetXaxis() -> SetLabelSize(fLab[1]);
  hAverage -> GetXaxis() -> CenterTitle(fCnt);
  hAverage -> GetYaxis() -> SetTitle(sTitleY.Data());
  hAverage -> GetYaxis() -> SetTitleFont(fTxt);
  hAverage -> GetYaxis() -> SetTitleSize(fTit[1]);
  hAverage -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hAverage -> GetYaxis() -> SetLabelFont(fTxt);
  hAverage -> GetYaxis() -> SetLabelSize(fLab[1]);
  hAverage -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] -> SetMarkerColor(fColS[iSys]);
    hSys[iSys] -> SetMarkerStyle(fMarS);
    hSys[iSys] -> SetFillColor(fColS[iSys]);
    hSys[iSys] -> SetFillStyle(fFilS[iSys]);
    hSys[iSys] -> SetLineColor(fColS[iSys]);
    hSys[iSys] -> SetLineStyle(fLinS);
    hSys[iSys] -> SetLineWidth(fWidS);
    hSys[iSys] -> SetTitle(sTitle.Data());
    hSys[iSys] -> SetTitleFont(fTxt);
    hSys[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hSys[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hSys[iSys] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hSys[iSys] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hSys[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hSys[iSys] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hSys[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hSys[iSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hSys[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hSys[iSys] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hSys[iSys] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hSys[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hSys[iSys] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hSys[iSys] -> GetYaxis() -> CenterTitle(fCnt);
    hDif[iSys] -> SetMarkerColor(fColS[iSys]);
    hDif[iSys] -> SetMarkerStyle(fMarS);
    hDif[iSys] -> SetFillColor(fColS[iSys]);
    hDif[iSys] -> SetFillStyle(fFilS[iSys]);
    hDif[iSys] -> SetLineColor(fColS[iSys]);
    hDif[iSys] -> SetLineStyle(fLinS);
    hDif[iSys] -> SetLineWidth(fWidS);
    hDif[iSys] -> SetTitle(sTitle.Data());
    hDif[iSys] -> SetTitleFont(fTxt);
    hDif[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDif[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hDif[iSys] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDif[iSys] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDif[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hDif[iSys] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDif[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hDif[iSys] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDif[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hDif[iSys] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDif[iSys] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDif[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hDif[iSys] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDif[iSys] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.5};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hDefault, sLabelD.Data());
  if (UseAverage)
    leg -> AddEntry(hAverage, sLabelA.Data());
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    leg -> AddEntry(hSys[iSys], sLabelS[iSys].Data());
  }
  cout << "    Made legend." << endl;

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
  hDif[0] -> Draw("E2");
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hDif[iSys] -> Draw("E2 SAME");
  }
  line    -> Draw();
  pPad2   -> cd();
  hSys[0] -> Draw("E2");
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hSys[iSys] -> Draw("E2 SAME");
  }
  hDefault -> Draw("same");
  if (UseAverage)
    hAverage -> Draw("same");
  leg      -> Draw();
  txt      -> Draw();
  fOut     -> cd();
  cPlot    -> Write();
  cPlot    -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOut     -> cd();
  hDefault -> Write();
  hAverage -> Write();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] -> Write();
    hDif[iSys] -> Write();
  }
  fOut -> Close();
  fInD -> cd();
  fInD -> Close();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fInS[iSys] -> cd();
    fInS[iSys] -> Close();
  }
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
