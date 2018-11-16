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
static const UInt_t NSys(4);
static const UInt_t NPlot(2);
static const UInt_t NRebin(2);
static const UInt_t NRebinVar(16);
static const Bool_t DoRebin(false);
static const Bool_t UseAverage(false);
static const Bool_t DoVariableRebin(false);
static const Bool_t DoGammaSubtraction(false);



void CalculateSystematicError() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // io parameters
  const TString sOut("systematics.priorSys.et1115pt0230vz55pt0230pi0.r02a005rm1chrg.d16m11y2018.root");
  const TString sInD("sysEt1115r02/pp200r9.default.et1115vz55pt0230pi0.r02a005rm1chrg.p0m1k4n58t27.root");
  const TString sHistD("hUnfolded");
  const TString sInS[NSys]   = {"sysEt1115r02/pp200r9.pythiaSys.et1115vz55pi0.r02a005rm1chrg.p0m1k4n58t4.root", "sysEt1115r02/pp200r9.levySys1.et1115vz55pi0.r02a005rm1chrg.p1m1k4n58t4.root", "sysEt1115r02/pp200r9.powSysM.et1115vz55pi0.r02a005rm1chrg.p4m1k4n58t19.root", "sysEt1115r02/pp200r9.powSysP.et1115vz55pi0.r02a005rm1chrg.p4m1k4n58t23.root"};
  const TString sHistS[NSys] = {"hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded"};

  // general plot parameters
  const TString sTitle("");
  const TString sNameD("hDefault");
  const TString sNameA("hAverage");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rhoA^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sLabelD("default [pythia 6]");
  const TString sLabelDS("default [pythia 6]");
  const TString sLabelA("average");
  const TString sLabelAS("average");

  // variation parameters
  const TString sTitleRD("var. / def.");
  const TString sTitleRA("var. / avg.");
  const TString sNameV[NSys]  = {"hVarPy8", "hVarLevy", "hVarPowM", "hVarPowP"};
  const TString sNameR[NSys]  = {"hRatPy8", "hRatLevy", "hRatPowM", "hRatPowP"};
  const TString sLabelS[NSys] = {"pythia 8", "levy [n = 5.8, t = 0.4]", "power [b = 2.1 - 0.2]", "power [b = 2.1 + 0.2]"};
  const UInt_t  fColS[NSys]   = {810, 850, 870, 890};
  const UInt_t  fFilS[NSys]   = {3305, 3325, 3345, 3395};

  // systematic parameters
  const UInt_t  fColT(922);
  const UInt_t  fFilT(0);
  const TString sNameT("hTotal");
  const TString sNamePT("hPerTotal");
  const TString sTitleP("fractional error");
  const TString sLabelT("total systematic");
  const TString sNameS[NSys] = {"hSysPy8", "hSysLevy", "hSysPowM", "hSysPowP"};
  const TString sNameP[NSys] = {"hPerPy8", "hPerLevy", "hPerPowM", "hPerPowP"};

  // subtraction parameters
  const TString  sDefPi0("");
  const TString  sHistDefPi0("");
  const TString  sInPi0[NSys]   = {"", "", "", ""};
  const TString  sHistPi0[NSys] = {"", "", "", ""};
  const Double_t gammaPurity(0.520270);

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trigger, E_{T}^{trg} #in (11, 15) GeV");
  const TString sJet("anti-k_{T}, R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // misc parameters
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

  TH1D *hVar[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hVar[iSys] = (TH1D*) fInS[iSys] -> Get(sHistS[iSys].Data());
    if (!hVar[iSys]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hVar[" << iSys << "] = " << hVar[iSys]
           << endl;
      return;
    }
    hVar[iSys] -> SetName(sNameV[iSys].Data());
    hVar[iSys] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Grabbed histograms." << endl;


  // do gamma-subtraction (if need be)
  TFile *fDefPi0;
  TFile *fPi0[NSys];
  TH1D  *hDefPi0;
  TH1D  *hDefGam;
  TH1D  *hPi0[NSys];
  TH1D  *hGam[NSys];
  if (DoGammaSubtraction) {
    // subtract default distribtion
    fDefPi0 = new TFile(sDefPi0.Data(), "read");
    if (!fDefPi0) {
      cerr << "PANIC: couldn't open default pi0 file!" << endl;
      return;
    }
    hDefPi0 = (TH1D*) fDefPi0  -> Get(sHistDefPi0.Data());
    hDefGam = (TH1D*) hDefault -> Clone();
    if (!hDefPi0 || !hDefGam) {
      cerr << "PANIC: couldn't grab default pi0 or gamma histogram!\n"
           << "       hDefPi0 = " << hDefPi0 << ", hDefGam = " << hDefGam
           << endl;
      return;
    }
    hDefPi0  -> Scale(gammaPurity);
    hDefault -> Add(hDefGam, hDefPi0, 1., -1.);
    hDefault -> Scale(1. / (1. - gammaPurity));
    fDefPi0  -> Close();

    // subtract systematic distributions
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      fPi0[iSys] = new TFile(sInPi0[iSys].Data(), "read");
      if (!fPi0[iSys]) {
        cerr << "PANIC: couldn't open a pi0 file!\n"
             << "       fInPi0[" << iSys << "] = " << fInPi0[iSys]
             << endl;
        return;
      }
      hPi0[iSys] = (TH1D*) fPi0[iSys] -> Get(sHistPi0[iSys].Data());
      hGam[iSys] = (TH1D*) hVar[iSys] -> Clone();
      if (!hPi0[iSys] || !hGam[iSys]) {
        cerr << "PANIC: couldn't grab pi0 or gamma histogram!\n"
             << "       hPi0[" << iSys << "] = " << hPi0[iSys] << ", hGam[" << iSys << "] = " << hGam[iSys]
             << endl;
        return;
      }
      hPi0[iSys] -> Scale(gammaPurity);
      hVar[iSys] -> Add(hGam[iSys], hPi0[iSys], 1., -1.);
      hVar[iSys] -> Scale(1. / (1. - gammaPurity));
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
      hVar[iSys] -> Rebin(NRebin);

      const UInt_t nBinS = hVar[iSys] -> GetNbinsX();
      for (UInt_t iBinS = 1; iBinS < (nBinS + 1); iBinS++) {
        const Double_t valS = hVar[iSys] -> GetBinContent(iBinS);
        const Double_t errS = hVar[iSys] -> GetBinError(iBinS);
        const Double_t sizS = hVar[iSys] -> GetBinWidth(iBinS);
        hVar[iSys] -> SetBinContent(iBinS, valS / sizS);
        hVar[iSys] -> SetBinError(iBinS, errS / sizS);
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
      hVar[iSys] = (TH1D*) hVar[iSys] -> Rebin(NRebinVar - 1, sNameV[iSys].Data(), varRebin);

      const UInt_t nBinS = hVar[iSys] -> GetNbinsX();
      for (UInt_t iBinS = 1; iBinS < (nBinS + 1); iBinS++) {
        const Double_t valS = hVar[iSys] -> GetBinContent(iBinS);
        const Double_t errS = hVar[iSys] -> GetBinError(iBinS);
        const Double_t sizS = hVar[iSys] -> GetBinWidth(iBinS);
        hVar[iSys] -> SetBinContent(iBinS, valS / sizS);
        hVar[iSys] -> SetBinError(iBinS, errS / sizS);
      }
    }
    cout << "    Rebinned histograms (variable)." << endl;
  }


  // calculate average
  TH1D *hAverage = (TH1D*) hDefault -> Clone();
  hAverage -> SetName(sNameA.Data());
  hAverage -> Reset("ICE");
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hAverage -> Add(hVar[iSys], 1.);
  }
  hAverage -> Scale(1. / (Double_t) NSys);
  cout << "    Calculated average." << endl;


  // calculate ratio
  TH1D *hDiv[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hDiv[iSys] = (TH1D*) hDefault -> Clone();
    hDiv[iSys] -> SetName(sNameR[iSys].Data());
    hDiv[iSys] -> Reset("ICE");

    const UInt_t nBinR = hDiv[iSys] -> GetNbinsX();
    for (UInt_t iBinR = 1; iBinR < (nBinR + 1); iBinR++) {
      const Double_t def  = hDefault   -> GetBinContent(iBinR);
      const Double_t avg  = hAverage   -> GetBinContent(iBinR);
      const Double_t var  = hVar[iSys] -> GetBinContent(iBinR);
      const Double_t rawD = hDefault   -> GetBinError(iBinR);
      const Double_t rawA = hAverage   -> GetBinError(iBinR);
      const Double_t rawV = hVar[iSys] -> GetBinError(iBinR);
      const Double_t relD = (rawD / def);
      const Double_t relA = (rawA / avg);
      const Double_t relV = (rawV / var);
      const Double_t divD = (var / def);
      const Double_t divA = (var / avg);

      Double_t div = 0.;
      Double_t rel = 0.;
      if (UseAverage) {
        div = divA;
        rel = TMath::Sqrt(TMath::Abs((relV * relV) - (relA * relA)));
      }
      else {
        div = divD;
        rel = TMath::Sqrt(TMath::Abs((relV * relV) - (relD * relD)));
      }

      // check for nan's
      const Bool_t isGoodDef = ((TMath::Abs(def) > 0.) && !UseAverage);
      const Bool_t isGoodAvg = ((TMath::Abs(avg) > 0.) && UseAverage);
      const Bool_t isGoodVar = (TMath::Abs(var) > 0.);

      // assign bin values
      const Double_t val = div;
      const Double_t err = val * rel;
      if ((isGoodDef || isGoodAvg) && isGoodVar) {
        hDiv[iSys] -> SetBinContent(iBinR, val);
        hDiv[iSys] -> SetBinError(iBinR, err);
      }
      else {
        hDiv[iSys] -> SetBinContent(iBinR, 0.);
        hDiv[iSys] -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
    hDiv[iSys] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Calculated difference." << endl;


  // calculate systematics
  TH1D *hSys[NSys + 1];
  TH1D *hPer[NSys + 1];
  if (UseAverage) {
    hSys[NSys] = (TH1D*) hAverage -> Clone();
    hPer[NSys] = (TH1D*) hAverage -> Clone();
  }
  else {
    hSys[NSys] = (TH1D*) hDefault -> Clone();
    hPer[NSys] = (TH1D*) hDefault -> Clone();
  }
  hSys[NSys] -> SetName(sNameT.Data());
  hPer[NSys] -> SetName(sNamePT.Data());

  // clear total systematic histograms
  const UInt_t nBinTot = hSys[NSys] -> GetNbinsX();
  for (UInt_t iBinTot = 0; iBinTot < nBinTot; iBinTot++) {
    hPer[NSys] -> SetBinContent(iBinTot, 0.);
    hPer[NSys] -> SetBinError(iBinTot, 0.);
    hSys[NSys] -> SetBinError(iBinTot, 0.);
  }

  // loop over systematic variations
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {

    // define histograms
    if (UseAverage) {
      hSys[iSys] = (TH1D*) hAverage -> Clone();
      hPer[iSys] = (TH1D*) hAverage -> Clone();
    }
    else {
      hSys[iSys] = (TH1D*) hDefault -> Clone();
      hPer[iSys] = (TH1D*) hDefault -> Clone();
    }
    hSys[iSys] -> SetName(sNameS[iSys].Data());
    hPer[iSys] -> SetName(sNameP[iSys].Data());

    const UInt_t nBins = hSys[iSys] -> GetNbinsX();
    for (UInt_t iBin = 0; iBin < nBins; iBin++) {

      const Double_t defVal = hDefault   -> GetBinContent(iBin);
      const Double_t defErr = hDefault   -> GetBinError(iBin);
      const Double_t avgVal = hAverage   -> GetBinContent(iBin);
      const Double_t avgErr = hAverage   -> GetBinError(iBin);
      const Double_t varVal = hVar[iSys] -> GetBinContent(iBin);
      const Double_t varErr = hVar[iSys] -> GetBinError(iBin);

      Double_t valDif(0.);
      Double_t errDif(0.);
      Double_t sysErr(0.);
      if (UseAverage) {
        valDif = TMath::Abs(avgVal - varVal);
        errDif = TMath::Sqrt(TMath::Abs((avgErr * avgErr) - (varErr * varErr)));
      }
      else {
        valDif = TMath::Abs(defVal - varVal);
        errDif = TMath::Sqrt(TMath::Abs((defErr * defErr) - (varErr * varErr)));
      }
      if (valDif > errDif)
        sysErr = TMath::Sqrt((valDif * valDif) - (errDif * errDif));
      else
        sysErr = 0.;

      Double_t perErr(0.);
      Double_t totErr(0.);
      Double_t totPer(0.);
      totErr = hSys[NSys] -> GetBinError(iBin);
      totErr = TMath::Sqrt((totErr * totErr) + (sysErr * sysErr));
      if (sysErr > 0.) {
        totPer = totErr / (hSys[NSys] -> GetBinContent(iBin));
        if (UseAverage)
          perErr = sysErr / avgVal;
        else
          perErr = sysErr / defVal;
      }
      hPer[iSys] -> SetBinContent(iBin, perErr);
      hPer[NSys] -> SetBinContent(iBin, totPer);
      hPer[iSys] -> SetBinError(iBin, 0.);
      hPer[NSys] -> SetBinError(iBin, 0.);
      hSys[iSys] -> SetBinError(iBin, sysErr);
      hSys[NSys] -> SetBinError(iBin, totErr);

    }  // end bin loop
  }  // end systematic loop
  cout << "    Calculated systematic uncertainties." << endl;


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
  const UInt_t  fMarT(8);
  const UInt_t  fFilD(0);
  const UInt_t  fFilA(0);
  const UInt_t  fLinD(1);
  const UInt_t  fLinA(1);
  const UInt_t  fLinS(1);
  const UInt_t  fLinT(1);
  const UInt_t  fWidD(1);
  const UInt_t  fWidA(1);
  const UInt_t  fWidS(1);
  const UInt_t  fWidT(1);
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
    hVar[iSys] -> SetMarkerColor(fColS[iSys]);
    hVar[iSys] -> SetMarkerStyle(fMarS);
    hVar[iSys] -> SetFillColor(fColS[iSys]);
    hVar[iSys] -> SetFillStyle(fFilS[iSys]);
    hVar[iSys] -> SetLineColor(fColS[iSys]);
    hVar[iSys] -> SetLineStyle(fLinS);
    hVar[iSys] -> SetLineWidth(fWidS);
    hVar[iSys] -> SetTitle(sTitle.Data());
    hVar[iSys] -> SetTitleFont(fTxt);
    hVar[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVar[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hVar[iSys] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVar[iSys] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVar[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hVar[iSys] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVar[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hVar[iSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVar[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hVar[iSys] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVar[iSys] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVar[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hVar[iSys] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVar[iSys] -> GetYaxis() -> CenterTitle(fCnt);
    hDiv[iSys] -> SetMarkerColor(fColS[iSys]);
    hDiv[iSys] -> SetMarkerStyle(fMarS);
    hDiv[iSys] -> SetFillColor(fColS[iSys]);
    hDiv[iSys] -> SetFillStyle(fFilS[iSys]);
    hDiv[iSys] -> SetLineColor(fColS[iSys]);
    hDiv[iSys] -> SetLineStyle(fLinS);
    hDiv[iSys] -> SetLineWidth(fWidS);
    hDiv[iSys] -> SetTitle(sTitle.Data());
    hDiv[iSys] -> SetTitleFont(fTxt);
    hDiv[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDiv[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hDiv[iSys] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDiv[iSys] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDiv[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hDiv[iSys] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDiv[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hDiv[iSys] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDiv[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hDiv[iSys] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDiv[iSys] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDiv[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hDiv[iSys] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDiv[iSys] -> GetYaxis() -> CenterTitle(fCnt);
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
    hPer[iSys] -> SetMarkerColor(fColS[iSys]);
    hPer[iSys] -> SetMarkerStyle(fMarS);
    hPer[iSys] -> SetFillColor(fColS[iSys]);
    hPer[iSys] -> SetFillStyle(fFilS[iSys]);
    hPer[iSys] -> SetLineColor(fColS[iSys]);
    hPer[iSys] -> SetLineStyle(fLinS);
    hPer[iSys] -> SetLineWidth(fWidS);
    hPer[iSys] -> SetTitle(sTitle.Data());
    hPer[iSys] -> SetTitleFont(fTxt);
    hPer[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPer[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hPer[iSys] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hPer[iSys] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hPer[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hPer[iSys] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hPer[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hPer[iSys] -> GetYaxis() -> SetTitle(sTitleP.Data());
    hPer[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hPer[iSys] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hPer[iSys] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hPer[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hPer[iSys] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hPer[iSys] -> GetYaxis() -> CenterTitle(fCnt);
  }
  hSys[NSys] -> SetMarkerColor(fColT);
  hSys[NSys] -> SetMarkerStyle(fMarT);
  hSys[NSys] -> SetFillColor(fColT);
  hSys[NSys] -> SetFillStyle(fFilT);
  hSys[NSys] -> SetLineColor(fColT);
  hSys[NSys] -> SetLineStyle(fLinT);
  hSys[NSys] -> SetLineWidth(fWidT);
  hSys[NSys] -> SetTitle(sTitle.Data());
  hSys[NSys] -> SetTitleFont(fTxt);
  hSys[NSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSys[NSys] -> GetXaxis() -> SetTitleFont(fTxt);
  hSys[NSys] -> GetXaxis() -> SetTitleSize(fTit[1]);
  hSys[NSys] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hSys[NSys] -> GetXaxis() -> SetLabelFont(fTxt);
  hSys[NSys] -> GetXaxis() -> SetLabelSize(fLab[1]);
  hSys[NSys] -> GetXaxis() -> CenterTitle(fCnt);
  hSys[NSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSys[NSys] -> GetYaxis() -> SetTitleFont(fTxt);
  hSys[NSys] -> GetYaxis() -> SetTitleSize(fTit[1]);
  hSys[NSys] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hSys[NSys] -> GetYaxis() -> SetLabelFont(fTxt);
  hSys[NSys] -> GetYaxis() -> SetLabelSize(fLab[1]);
  hSys[NSys] -> GetYaxis() -> CenterTitle(fCnt);
  hPer[NSys] -> SetMarkerColor(fColT);
  hPer[NSys] -> SetMarkerStyle(fMarT);
  hPer[NSys] -> SetFillColor(fColT);
  hPer[NSys] -> SetFillStyle(fFilT);
  hPer[NSys] -> SetLineColor(fColT);
  hPer[NSys] -> SetLineStyle(fLinT);
  hPer[NSys] -> SetLineWidth(fWidT);
  hPer[NSys] -> SetTitle(sTitle.Data());
  hPer[NSys] -> SetTitleFont(fTxt);
  hPer[NSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPer[NSys] -> GetXaxis() -> SetTitleFont(fTxt);
  hPer[NSys] -> GetXaxis() -> SetTitleSize(fTit[0]);
  hPer[NSys] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hPer[NSys] -> GetXaxis() -> SetLabelFont(fTxt);
  hPer[NSys] -> GetXaxis() -> SetLabelSize(fLab[0]);
  hPer[NSys] -> GetXaxis() -> CenterTitle(fCnt);
  hPer[NSys] -> GetYaxis() -> SetTitle(sTitleP.Data());
  hPer[NSys] -> GetYaxis() -> SetTitleFont(fTxt);
  hPer[NSys] -> GetYaxis() -> SetTitleSize(fTit[0]);
  hPer[NSys] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hPer[NSys] -> GetYaxis() -> SetLabelFont(fTxt);
  hPer[NSys] -> GetYaxis() -> SetLabelSize(fLab[0]);
  hPer[NSys] -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legends
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.5};
  TLegend *lVar = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  TLegend *lSys = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  lVar -> SetFillColor(fColLe);
  lVar -> SetFillStyle(fFilLe);
  lVar -> SetLineColor(fColLe);
  lVar -> SetLineStyle(fLinLe);
  lVar -> SetTextFont(fTxt);
  lVar -> SetTextAlign(fAln);
  lSys -> SetFillColor(fColLe);
  lSys -> SetFillStyle(fFilLe);
  lSys -> SetLineColor(fColLe);
  lSys -> SetLineStyle(fLinLe);
  lSys -> SetTextFont(fTxt);
  lSys -> SetTextAlign(fAln);
  if (UseAverage) {
    lVar -> AddEntry(hAverage, sLabelA.Data());
    lSys -> AddEntry(hAverage, sLabelAS.Data());
  }
  else {
    lSys -> AddEntry(hDefault, sLabelDS.Data());
  }
  lVar -> AddEntry(hDefault, sLabelD.Data());
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    lVar -> AddEntry(hVar[iSys], sLabelS[iSys].Data());
    lSys -> AddEntry(hSys[iSys], sLabelS[iSys].Data());
  }
  lSys -> AddEntry(hSys[NSys], sLabelT.Data());
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

  TCanvas *cPlot1 = new TCanvas("cVariations", "", width, height);
  TPad    *pPad1  = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2  = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot1  -> SetGrid(fGrid, fGrid);
  cPlot1  -> SetTicks(fTick, fTick);
  cPlot1  -> SetBorderMode(fMode);
  cPlot1  -> SetBorderSize(fBord);
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
  cPlot1  -> cd();
  pPad1   -> Draw();
  pPad2   -> Draw();
  pPad1   -> cd();
  hDiv[0] -> Draw("E2");
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hDiv[iSys] -> Draw("E2 SAME");
  }
  line    -> Draw();
  pPad2   -> cd();
  hDefault -> Draw("E2");
  if (UseAverage)
    hAverage -> Draw("same");
  hVar[0] -> Draw("P hist same");
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hVar[iSys] -> Draw("P hist same");
  }
  lVar   -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cPlot1 -> Write();
  cPlot1 -> Close();

  TCanvas *cPlot2 = new TCanvas("cSystematics", "", width, height);
  TPad    *pPad3 = new TPad("pPad3", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad4 = new TPad("pPad4", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot2  -> SetGrid(fGrid, fGrid);
  cPlot2  -> SetTicks(fTick, fTick);
  cPlot2  -> SetBorderMode(fMode);
  cPlot2  -> SetBorderSize(fBord);
  pPad3   -> SetGrid(fGrid, fGrid);
  pPad3   -> SetTicks(fTick, fTick);
  pPad3   -> SetLogx(fLogX);
  pPad3   -> SetLogy(fLogY);
  pPad3   -> SetBorderMode(fMode);
  pPad3   -> SetBorderSize(fBord);
  pPad3   -> SetFrameBorderMode(fFrame);
  pPad3   -> SetLeftMargin(fMarginL);
  pPad3   -> SetRightMargin(fMarginR);
  pPad3   -> SetTopMargin(fMarginT1);
  pPad3   -> SetBottomMargin(fMarginB1);
  pPad4   -> SetGrid(fGrid, fGrid);
  pPad4   -> SetTicks(fTick, fTick);
  pPad4   -> SetLogx(fLogX);
  pPad4   -> SetLogy(fLogY);
  pPad4   -> SetBorderMode(fMode);
  pPad4   -> SetBorderSize(fBord);
  pPad4   -> SetFrameBorderMode(fFrame);
  pPad4   -> SetLeftMargin(fMarginL);
  pPad4   -> SetRightMargin(fMarginR);
  pPad4   -> SetTopMargin(fMarginT2);
  pPad4   -> SetBottomMargin(fMarginB2);
  cPlot2  -> cd();
  pPad3   -> Draw();
  pPad4   -> Draw();
  pPad3   -> cd();
  hPer[0] -> Draw();
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hPer[iSys] -> Draw("SAME");
  }
  hPer[NSys] -> Draw("SAME");
  line       -> Draw();
  pPad4      -> cd();
  hSys[0]    -> Draw("E2");
  for (UInt_t iSys = 1; iSys < NSys; iSys++) {
    hSys[iSys] -> Draw("E2 SAME");
  }
  hSys[NSys] -> Draw("E2 SAME");
  if (UseAverage)
    hAverage -> Draw("E2 SAME");
  else
    hDefault -> Draw("E2 SAME");
  lSys   -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cPlot2 -> Write();
  cPlot2 -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOut     -> cd();
  hDefault -> Write();
  hAverage -> Write();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hVar[iSys] -> Write();
    hDiv[iSys] -> Write();
    hSys[iSys] -> Write();
    hPer[iSys] -> Write();
  }
  hSys[NSys] -> Write();
  hPer[NSys] -> Write();
  fOut       -> Close();
  fInD       -> cd();
  fInD       -> Close();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fInS[iSys] -> cd();
    fInS[iSys] -> Close();
  }
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
