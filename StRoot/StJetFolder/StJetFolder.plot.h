// 'StJetFolder.plot.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates various routines associated with plotting results.
//
// Last updated: 02.17.2017


#pragma once

using namespace std;



void StJetFolder::CreateLabel() {

  _label = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC NB");
  _label -> SetLineColor(0);
  _label -> SetFillColor(0);
  _label -> SetFillStyle(0);
  _label -> SetTextColor(1);
  _label -> SetTextAlign(12);
  _label -> SetTextFont(42);
  _label -> AddText(_sEvnt -> Data());
  _label -> AddText(_sTrig -> Data());
  _label -> AddText(_sJet1 -> Data());
  _label -> AddText(_sJet2 -> Data());
  _label -> AddText(_sJet3 -> Data()); 

}  // end 'Createlabel()'


void StJetFolder::CreatePlots() {

  _fOut -> cd();
  PrintInfo(10);


  // histogram colors
  const Int_t    cM  = 810;
  const Int_t    cU  = 870;
  const Int_t    cB  = 800;
  const Int_t    cP  = 860;
  const Int_t    cS  = 860;
  const Int_t    cBM = 810;
  const Int_t    cUP = 860;
  const Int_t    cSM = 418;
  const Int_t    cL  = 1;
  // histogram marker styles
  const Int_t    mM = 29;
  const Int_t    mU = 24;
  const Int_t    mB = 24;
  const Int_t    mP = 29;
  const Int_t    mS = 24;
  const Int_t    mR = 7;
  // histogram line styles
  const Int_t    lM = 1;
  const Int_t    lU = 1;
  const Int_t    lB = 1;
  const Int_t    lP = 1;
  const Int_t    lS = 1;
  const Int_t    lR = 1;
  const Int_t    lL = 2;
  // histogram fill styles
  const Int_t    fM = 0;
  const Int_t    fU = 3017;
  const Int_t    fB = 3018;
  const Int_t    fP = 0;
  const Int_t    fS = 0;
  const Int_t    fR = 0;
  // histogram dimensions
  const Int_t    nYup = 1000;
  const Int_t    nYlo = 65;
  const Double_t x1   = 0.;
  const Double_t x2   = 60.;
  const Double_t y1up = 0.000005;
  const Double_t y2up = 0.5;
  const Double_t y1lo = 0.;
  const double_t y2lo = 6.5;
  // axis titles
  const TString sX("p_{T}^{jet} [GeV/c]");
  const TString sYup("(1/N_{trg}) dN_{jet}/(dp_{T}^{jet} d#eta^{jet})");
  const TString sYlo1("backfolded / measured");
  const TString sYlo2("prior / unfolded");
  const TString sYlo3("smeared / measured");


  // generate legends
  TLegend *lAll = new TLegend(0.3, 0.1, 0.5, 0.3);
  lAll -> SetFillColor(0);
  lAll -> SetFillStyle(0);
  lAll -> SetLineColor(0);
  lAll -> SetTextColor(1);
  lAll -> SetTextFont(42);
  lAll -> SetTextAlign(12);
  lAll -> AddEntry(_hMeasured, "Measured");
  lAll -> AddEntry(_hUnfolded, "Unfolded");
  lAll -> AddEntry(_hBackfolded, "Backfolded");
  lAll -> AddEntry(_hPrior, "Prior");

  TLegend *lUvB = new TLegend(0.3, 0.1, 0.5, 0.3);
  lUvB -> SetFillColor(0);
  lUvB -> SetFillStyle(0);
  lUvB -> SetLineColor(0);
  lUvB -> SetTextColor(1);
  lUvB -> SetTextFont(42);
  lUvB -> SetTextAlign(12);
  lUvB -> AddEntry(_hUnfoldVsPriRatio, "unfolded / prior");
  lUvB -> AddEntry(_hBackVsMeasRatio, "backfolded / measured");

  TLegend *lBvM = new TLegend(0.3, 0.1, 0.5, 0.3);
  lBvM -> SetFillColor(0);
  lBvM -> SetFillStyle(0);
  lBvM -> SetLineColor(0);
  lBvM -> SetTextColor(1);
  lBvM -> SetTextFont(42);
  lBvM -> SetTextAlign(12);
  lBvM -> AddEntry(_hMeasured, "Measured");
  lBvM -> AddEntry(_hBackfolded, "Backfolded");

  TLegend *lPvU = new TLegend(0.3, 0.1, 0.5, 0.3);
  lPvU -> SetFillColor(0);
  lPvU -> SetFillStyle(0);
  lPvU -> SetLineColor(0);
  lPvU -> SetTextColor(1);
  lPvU -> SetTextFont(42);
  lPvU -> SetTextAlign(12);
  lPvU -> AddEntry(_hUnfolded, "Unfolded");
  lPvU -> AddEntry(_hPrior, "Prior");

  TLegend *lSvM = new TLegend(0.3, 0.1, 0.5, 0.3);
  lSvM -> SetFillColor(0);
  lSvM -> SetFillStyle(0);
  lSvM -> SetLineColor(0);
  lSvM -> SetTextColor(1);
  lSvM -> SetTextFont(42);
  lSvM -> SetTextAlign(12);
  lSvM -> AddEntry(_hMeasured, "Measured");
  lSvM -> AddEntry(_hSmeared, "Smeared Prior");


  // create chi2 text
  TString x2txt("");
  TString x2str("");
  x2str += _chi2backfold;

  ResizeString(x2str, 2);
  x2txt.Append("#chi^{2}/dof = ");
  x2txt.Append(x2str);

  TPaveText *pChi2 = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pChi2 -> SetFillColor(0);
  pChi2 -> SetFillStyle(0);
  pChi2 -> SetLineColor(0);
  pChi2 -> SetTextColor(1);
  pChi2 -> SetTextFont(42);
  pChi2 -> SetTextAlign(12);
  pChi2 -> AddText(x2txt.Data());



  // create blank histograms
  const Int_t    nM   = _hMeasured -> GetNbinsX();
  const Int_t    x1m  = _hMeasured -> GetBinLowEdge(1);
  const Int_t    x2m  = _hMeasured -> GetBinLowEdge(nM + 1);
  const Double_t wM   = (x2m - x1m) / nM;
  const TString  *sUp = CreateTitle();

  Int_t  nX  = (Int_t) ((x2 - x1) * wM);
  TH2D  *hUp = new TH2D("hUpper", "", nX, x1, x2, nYup, y1up, y2up);
  hUp -> GetXaxis() -> SetTitle("");
  hUp -> GetXaxis() -> SetTitleFont(42);
  hUp -> GetXaxis() -> SetTitleSize(0.0);
  hUp -> GetXaxis() -> CenterTitle(true);
  hUp -> GetXaxis() -> SetLabelFont(42);
  hUp -> GetXaxis() -> SetLabelSize(0.0);
  hUp -> GetYaxis() -> SetTitle(sYup);
  hUp -> GetYaxis() -> SetTitleFont(42);
  hUp -> GetYaxis() -> SetTitleSize(0.05);
  hUp -> GetYaxis() -> CenterTitle(true);
  hUp -> GetYaxis() -> SetLabelFont(42);
  hUp -> GetYaxis() -> SetLabelSize(0.05);
  hUp -> SetTitle(sUp -> Data());
  hUp -> SetTitleFont(42);

  Int_t nD  = (Int_t) y2lo;
  TH2D *hLo = new TH2D("hLower", "", nX, x1, x2, nYlo, y1lo, y2lo);
  hLo -> GetXaxis() -> SetTitle(sX);
  hLo -> GetXaxis() -> SetTitleFont(42);
  hLo -> GetXaxis() -> SetTitleSize(0.0);
  hLo -> GetXaxis() -> CenterTitle(true);
  hLo -> GetXaxis() -> SetLabelFont(42);
  hLo -> GetXaxis() -> SetLabelSize(0.0);
  hLo -> GetYaxis() -> SetTitle(sYlo1);
  hLo -> GetYaxis() -> SetTitleFont(42);
  hLo -> GetYaxis() -> SetTitleSize(0.05);
  hLo -> GetYaxis() -> CenterTitle(true);
  hLo -> GetYaxis() -> SetLabelFont(42);
  hLo -> GetYaxis() -> SetLabelSize(0.05);
  hLo -> GetYaxis() -> SetNdivisions(nD, 5, 0);


  // create line
  TLine *lOne = new TLine(x1, 1., x2, 1.);
  lOne -> SetLineStyle(lL);
  lOne -> SetLineColor(cL);


  // draw plots
  TCanvas *cAll = new TCanvas("cAll", "All 4 spectra", 850, 850);
  TPad    *pLoA = new TPad("pLoA", "backfold vs. measured ratio", 0, 0, 1, 0.33);
  TPad    *pUpA = new TPad("pUpA", "all 4 spectra", 0, 0.33, 1, 1);
  // set plot options
  pLoA   -> SetFillStyle(4000);
  pLoA   -> SetFillColor(0);
  pLoA   -> SetBorderMode(0);
  pLoA   -> SetFrameBorderMode(0);
  pLoA   -> SetTopMargin(0);
  pLoA   -> SetGrid(0, 0);
  pLoA   -> SetTickx(1);
  pLoA   -> SetTicky(1);
  pUpA   -> SetFillStyle(4000);
  pUpA   -> SetFillColor(0);
  pUpA   -> SetBorderMode(0);
  pUpA   -> SetFrameBorderMode(0);
  pUpA   -> SetBottomMargin(0);
  pUpA   -> SetGrid(0, 0);
  pUpA   -> SetTickx(1);
  pUpA   -> SetTicky(1);
  pUpA   -> SetLogy(1);
  pLoA   -> Draw();
  pUpA   -> Draw();
  // draw histograms
  pLoA   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hUnfoldVsPriRatio, "PE2 same", cUP, cUP, cUP, mR, lR, fR, 1.);
  DrawHistogram(_hBackVsMeasRatio, "PE2 same", cBM, cBM, cBM, mR, lR, fR, 1.);
  lUvB   -> Draw();
  lOne   -> Draw();
  pChi2  -> Draw();
  pUpA   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  DrawHistogram(_hUnfolded, "PE6 same", cU, cU, cU, mU, lU, fU, 1.);
  DrawHistogram(_hBackfolded, "PE6 same", cB, cB, cB, mB, lB, fB, 1.);
  DrawHistogram(_hPrior, "PE2 same", cP, cP, cP, mP, lP, fP, 1.);
  lAll   -> Draw();
  _label -> Draw();
  cAll   -> Write();
  cAll   -> Close();

  TCanvas *cBvM = new TCanvas("cBvM", "Backfold vs. measured", 850, 850);
  TPad    *pLo1 = new TPad("pLo1", "ratio", 0, 0, 1, 0.33);
  TPad    *pUp1 = new TPad("pUp1", "backfold and measured", 0, 0.33, 1, 1);
  // set plot options
  pLo1   -> SetFillStyle(4000);
  pLo1   -> SetFillColor(0);
  pLo1   -> SetBorderMode(0);
  pLo1   -> SetFrameBorderMode(0);
  pLo1   -> SetTopMargin(0);
  pLo1   -> SetGrid(0, 0);
  pLo1   -> SetTickx(1);
  pLo1   -> SetTicky(1);
  pUp1   -> SetFillStyle(4000);
  pUp1   -> SetFillColor(0);
  pUp1   -> SetBorderMode(0);
  pUp1   -> SetFrameBorderMode(0);
  pUp1   -> SetBottomMargin(0);
  pUp1   -> SetGrid(0, 0);
  pUp1   -> SetTickx(1);
  pUp1   -> SetTicky(1);
  pUp1   -> SetLogy(1);
  pLo1   -> Draw();
  pUp1   -> Draw();
  // draw histograms
  pLo1   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hBackVsMeasRatio, "PE2 same", cBM, cBM, cBM, mR, lR, fR, 1.);
  lOne   -> Draw();
  pChi2  -> Draw();
  pUp1   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  DrawHistogram(_hBackfolded, "PE2 same", cB, cB, cB, mB, lB, fB, 1.);
  lBvM   -> Draw();
  _label -> Draw();
  cBvM   -> Write();
  cBvM   -> Close();

  TCanvas *cPvU = new TCanvas("cPvU", "Prior vs. unfolded", 850, 850);
  TPad    *pLo2 = new TPad("pLo2", "ratio", 0, 0, 1, 0.33);
  TPad    *pUp2 = new TPad("pUp2", "prior and unfolded", 0, 0.33, 1, 1);
  // set plot options
  pLo2   -> SetFillStyle(4000);
  pLo2   -> SetFillColor(0);
  pLo2   -> SetBorderMode(0);
  pLo2   -> SetFrameBorderMode(0);
  pLo2   -> SetTopMargin(0);
  pLo2   -> SetGrid(0, 0);
  pLo2   -> SetTickx(1);
  pLo2   -> SetTicky(1);
  pUp2   -> SetFillStyle(4000);
  pUp2   -> SetFillColor(0);
  pUp2   -> SetBorderMode(0);
  pUp2   -> SetFrameBorderMode(0);
  pUp2   -> SetBottomMargin(0);
  pUp2   -> SetGrid(0, 0);
  pUp2   -> SetTickx(1);
  pUp2   -> SetTicky(1);
  pUp2   -> SetLogy(1);
  pLo2   -> Draw();
  pUp2   -> Draw();
  // draw histograms
  pLo2   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hUnfoldVsPriRatio, "PE2 same", cUP, cUP, cUP, mR, lR, fR, 1.);
  lOne   -> Draw();
  pChi2  -> Draw();
  pUp2   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hUnfolded, "PE2 same", cU, cU, cU, mU, lU, fU, 1.);
  DrawHistogram(_hPrior, "PE2 same", cP, cP, cP, mP, lP, fP, 1.);
  lPvU   -> Draw();
  _label -> Draw();
  cPvU   -> Write();
  cPvU   -> Close();

  TCanvas *cSvM = new TCanvas("cSvM", "Smear vs. measured", 850, 850);
  TPad    *pLo3 = new TPad("pLo3", "ratio", 0, 0, 1, 0.33);
  TPad    *pUp3 = new TPad("pUp3", "smeared and measured", 0, 0.33, 1, 1);
  // set plot options
  pLo3   -> SetFillStyle(4000);
  pLo3   -> SetFillColor(0);
  pLo3   -> SetBorderMode(0);
  pLo3   -> SetFrameBorderMode(0);
  pLo3   -> SetTopMargin(0);
  pLo3   -> SetGrid(0, 0);
  pLo3   -> SetTickx(1);
  pLo3   -> SetTicky(1);
  pUp3   -> SetFillStyle(4000);
  pUp3   -> SetFillColor(0);
  pUp3   -> SetBorderMode(0);
  pUp3   -> SetFrameBorderMode(0);
  pUp3   -> SetBottomMargin(0);
  pUp3   -> SetGrid(0, 0);
  pUp3   -> SetTickx(1);
  pUp3   -> SetTicky(1);
  pUp3   -> SetLogy(1);
  pLo3   -> Draw();
  pUp3   -> Draw();
  // draw histograms
  pLo3   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hSmearVsMeasRatio, "PE2 same", cSM, cSM, cSM, mR, lR, fR, 1.);
  lOne   -> Draw();
  pChi2  -> Draw();
  pUp3   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  DrawHistogram(_hSmeared, "PE2 same", cS, cS, cS, mS, lS, fS, 1.);
  lSvM   -> Draw();
  _label -> Draw();
  cSvM   -> Write();
  cSvM   -> Close();


}  // end 'CreatePlots()'


void StJetFolder::ResizeString(TString &str, const Int_t nDec) {

  Int_t nDig = 0;
  Int_t size = 0;

  nDig = str.First(".");
  if (nDig > 0)
    size = (nDig + 1) + nDec;
  else
    size = str.Length();
  str.Resize(size);

}  // end 'ResizeString(TString&, Int_t)'


void StJetFolder::DrawHistogram(TH1 *h, const Char_t *option, const Int_t mColor, const Int_t lColor, const Int_t fColor, const Int_t mStyle, const Int_t lStyle, const Int_t fStyle, const Double_t mSize) {

  h -> SetMarkerColor(mColor);
  h -> SetLineColor(lColor);
  h -> SetFillColor(fColor);
  h -> SetMarkerStyle(mStyle);
  h -> SetLineStyle(lStyle);
  h -> SetFillStyle(fStyle);
  h -> SetMarkerSize(mSize);
  h -> Draw(option);

}  // end 'DrawHistogram(TH1*, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Char_t*)'


TString* StJetFolder::CreateTitle() {

  // determine prior
  TString pTxt("");
  switch (_prior) {
    case 0:
      pTxt.Append("Pythia");
      break;
    case 1:
      pTxt.Append("Levy");
      break;
    case 2:
      pTxt.Append("Tsallis");
      break;
    case 3:
      pTxt.Append("Exponential");
      break;
  }

  // determine method
  TString mTxt("");
  switch (_method) {
    case 0:
      mTxt.Append("No unfolding");
      break;
    case 1:
      mTxt.Append("Bayes.");
      break;
    case 2:
      mTxt.Append("SVD");
      break;
    case 3:
      mTxt.Append("Bin-by-bin");
      break;
    case 4:
      mTxt.Append("TUnfold");
      break;
    case 5:
      mTxt.Append("Matrix inversion");
      break;
  }

  // convert kReg to a string
  TString kTxt("k_{reg} = ");
  kTxt += _kReg;

  // convert nPrior to a string
  TString nTxt("");
  TString nStr("");
  nStr += _nPrior;
  ResizeString(nStr, 2);
  nTxt.Append("N_{prior} = ");
  nTxt.Append(nStr);

  // convert tPrior to a string
  TString tTxt("");
  TString tStr("");
  tStr += _tPrior;
  ResizeString(tStr, 2);
  tTxt.Append("T_{prior} = ");
  tTxt.Append(tStr);


  // create title
  TString title("");
  title.Append(mTxt);
  title.Append(", ");
  title.Append(kTxt);
  title.Append(": ");
  title.Append(pTxt);
  title.Append(", ");
  title.Append(nTxt);
  title.Append(", ");
  title.Append(tTxt);

  return new TString(title.Data());

}  // end 'CreateTitle()'

// End ------------------------------------------------------------------------
