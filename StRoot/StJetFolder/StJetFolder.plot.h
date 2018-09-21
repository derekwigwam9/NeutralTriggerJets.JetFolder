// 'StJetFolder.plot.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates various routines associated with plotting results.
//
// Last updated: 09.12.2018


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
  _label -> AddText(_sJet3 -> Data()); 

}  // end 'Createlabel()'


void StJetFolder::CreatePlots() {

  _fOut -> cd();
  PrintInfo(10);


  // determine colors
  UInt_t cUnfold(0);
  UInt_t cBackfold(0);
  switch (_trigger) {
    case 0:
      cUnfold   = 899;
      cBackfold = 896;
      break;
    case 1:
      cUnfold   = 879;
      cBackfold = 876;
      break;
    case 2:
      cUnfold   = 859;
      cBackfold = 856;
      break;
    default:
      cUnfold   = 870;
      cBackfold = 800;
      break;
  }

  // histogram colors
  const Int_t    cM  = 923;
  const Int_t    cU  = cUnfold;
  const Int_t    cB  = cBackfold;
  const Int_t    cP  = 921;
  const Int_t    cS  = cBackfold;
  const Int_t    cBM = cBackfold;
  const Int_t    cUP = cUnfold;
  const Int_t    cSM = cBackfold;
  const Int_t    cUM = cUnfold;
  const Int_t    cSP = cBackfold;
  const Int_t    cRP = 1;
  const Int_t    cL  = 923;
  // histogram marker styles
  const Int_t    mM  = 29;
  const Int_t    mU  = 24;
  const Int_t    mB  = 24;
  const Int_t    mP  = 20;
  const Int_t    mS  = 24;
  const Int_t    mR  = 24;
  const Int_t    mRP = 20;
  // histogram line styles
  const Int_t    lM  = 1;
  const Int_t    lU  = 1;
  const Int_t    lB  = 1;
  const Int_t    lP  = 1;
  const Int_t    lS  = 1;
  const Int_t    lR  = 1;
  const Int_t    lRP = 1;
  const Int_t    lL  = 2;
  // histogram fill styles
  const Int_t    fM  = 0;
  const Int_t    fU  = 3017;
  const Int_t    fB  = 3018;
  const Int_t    fP  = 0;
  const Int_t    fS  = 0;
  const Int_t    fRU = 3017;
  const Int_t    fRB = 3018;
  const Int_t    fRS = 0;
  const Int_t    fRP = 0;
  // histogram dimensions
  const Double_t x1   = -5.;
  const Double_t x2   = 47.;
  const Double_t y1up = 0.0000005;
  const Double_t y2up = 5.;
  const Double_t y1lo = 0.;
  const double_t y2lo = 3.33;
  // axis titles
  const TString sX("p_{T}^{jet} [GeV/c]");
  const TString sXres("p_{T}^{jet}(det.) [GeV/c]");
  const TString sXeff("p_{T}^{jet}(par.) [GeV/c]");
  const TString sYup("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sYlo("ratio");
  const TString sYres("p_{T}^{jet}(par.) [GeV/c]");
  const TString sYeff("#epsilon_{jet}(p_{T}^{jet})");


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

  TLegend *lUvM = new TLegend(0.3, 0.1, 0.5, 0.3);
  lUvM -> SetFillColor(0);
  lUvM -> SetFillStyle(0);
  lUvM -> SetLineColor(0);
  lUvM -> SetTextColor(1);
  lUvM -> SetTextFont(42);
  lUvM -> SetTextAlign(12);
  lUvM -> AddEntry(_hMeasured, "Measured");
  lUvM -> AddEntry(_hUnfolded, "Unfolded");

  TLegend *lSvP = new TLegend(0.3, 0.1, 0.5, 0.3);
  lSvP -> SetFillColor(0);
  lSvP -> SetFillStyle(0);
  lSvP -> SetLineColor(0);
  lSvP -> SetTextColor(1);
  lSvP -> SetTextFont(42);
  lSvP -> SetTextAlign(12);
  lSvP -> AddEntry(_hSmeared, "Smeared Prior");
  lSvP -> AddEntry(_hPrior, "Prior");


  // create chi2 texts
  TString x2txtB("");
  TString x2strB("");
  x2strB += _chi2backfold;

  TString x2txtU("");
  TString x2strU("");
  x2strU += _chi2unfold;

  ResizeString(x2strB, 2);
  x2txtB.Append("#chi^{2}/dof(backfold) = ");
  x2txtB.Append(x2strB);

  ResizeString(x2strU, 2);
  x2txtU.Append("#chi^{2}/dof(unfold) = ");
  x2txtU.Append(x2strU);

  TPaveText *pChi2B = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pChi2B -> SetFillColor(0);
  pChi2B -> SetFillStyle(0);
  pChi2B -> SetLineColor(0);
  pChi2B -> SetTextColor(1);
  pChi2B -> SetTextFont(42);
  pChi2B -> SetTextAlign(12);
  pChi2B -> AddText(x2txtB.Data());

  TPaveText *pChi2U = new TPaveText(0.7, 0.3, 0.9, 0.5, "NDC NB");
  pChi2U -> SetFillColor(0);
  pChi2U -> SetFillStyle(0);
  pChi2U -> SetLineColor(0);
  pChi2U -> SetTextColor(1);
  pChi2U -> SetTextFont(42);
  pChi2U -> SetTextAlign(12);
  pChi2U -> AddText(x2txtU.Data());



  // create blank histograms
  const Int_t    nM   = _hMeasured -> GetNbinsX();
  const Int_t    x1m  = _hMeasured -> GetBinLowEdge(1);
  const Int_t    x2m  = _hMeasured -> GetBinLowEdge(nM + 1);
  const Double_t wM   = (x2m - x1m) / nM;
  const TString  *sUp = CreateTitle();

  Int_t  nX  = (Int_t) ((x2 - x1) * wM);
  TH1D  *hUp = new TH1D("hUpper", "", nX, x1, x2);
  hUp -> GetXaxis() -> SetTitle("");
  hUp -> GetXaxis() -> SetTitleFont(42);
  hUp -> GetXaxis() -> SetTitleSize(0.04);
  hUp -> GetXaxis() -> CenterTitle(true);
  hUp -> GetXaxis() -> SetLabelFont(42);
  hUp -> GetXaxis() -> SetLabelSize(0.04);
  hUp -> GetXaxis() -> SetTitleOffset(1.);
  hUp -> GetYaxis() -> SetTitle(sYup);
  hUp -> GetYaxis() -> SetTitleFont(42);
  hUp -> GetYaxis() -> SetTitleSize(0.04);
  hUp -> GetYaxis() -> SetTitleOffset(1.3);
  hUp -> GetYaxis() -> CenterTitle(true);
  hUp -> GetYaxis() -> SetLabelFont(42);
  hUp -> GetYaxis() -> SetLabelSize(0.04);
  hUp -> GetYaxis() -> SetRangeUser(y1up, y2up);
  hUp -> SetTitle(sUp -> Data());
  hUp -> SetTitleFont(42);

  Int_t nD  = (Int_t) y2lo;
  TH1D *hLo = new TH1D("hLower", "", nX, x1, x2);
  hLo -> GetXaxis() -> SetTitle(sX);
  hLo -> GetXaxis() -> SetTitleFont(42);
  hLo -> GetXaxis() -> SetTitleSize(0.074);
  hLo -> GetXaxis() -> SetTitleOffset(1.1);
  hLo -> GetXaxis() -> CenterTitle(true);
  hLo -> GetXaxis() -> SetLabelFont(42);
  hLo -> GetXaxis() -> SetLabelSize(0.074);
  hLo -> GetYaxis() -> SetTitle(sYlo);
  hLo -> GetYaxis() -> SetTitleFont(42);
  hLo -> GetYaxis() -> SetTitleSize(0.074);
  hLo -> GetYaxis() -> SetTitleOffset(0.7);
  hLo -> GetYaxis() -> CenterTitle(true);
  hLo -> GetYaxis() -> SetLabelFont(42);
  hLo -> GetYaxis() -> SetLabelSize(0.074);
  hLo -> GetYaxis() -> SetNdivisions(nD, 5, 0);
  hLo -> GetYaxis() -> SetRangeUser(y1lo, y2lo);

  // set response and efficiency styles
  _hResponse -> GetXaxis() -> SetTitle(sXres);
  _hResponse -> GetXaxis() -> SetTitleFont(42);
  _hResponse -> GetXaxis() -> SetTitleSize(0.04);
  _hResponse -> GetXaxis() -> CenterTitle(true);
  _hResponse -> GetXaxis() -> SetLabelFont(42);
  _hResponse -> GetXaxis() -> SetLabelSize(0.04);
  _hResponse -> GetXaxis() -> SetTitleOffset(1.3);
  _hResponse -> GetXaxis() -> SetRangeUser(x1, x2);
  _hResponse -> GetYaxis() -> SetTitle(sYres);
  _hResponse -> GetYaxis() -> SetTitleFont(42);
  _hResponse -> GetYaxis() -> SetTitleSize(0.04);
  _hResponse -> GetYaxis() -> SetTitleOffset(1.3);
  _hResponse -> GetYaxis() -> CenterTitle(true);
  _hResponse -> GetYaxis() -> SetLabelFont(42);
  _hResponse -> GetYaxis() -> SetLabelSize(0.04);
  _hResponse -> GetYaxis() -> SetRangeUser(x1, x2);
  _hResponse -> GetZaxis() -> SetLabelFont(42);
  _hResponse -> GetZaxis() -> SetLabelSize(0.04);
  _hResponse -> SetTitleFont(42);

  _hEfficiency -> GetXaxis() -> SetTitle(sXeff);
  _hEfficiency -> GetXaxis() -> SetTitleFont(42);
  _hEfficiency -> GetXaxis() -> SetTitleSize(0.04);
  _hEfficiency -> GetXaxis() -> CenterTitle(true);
  _hEfficiency -> GetXaxis() -> SetLabelFont(42);
  _hEfficiency -> GetXaxis() -> SetLabelSize(0.04);
  _hEfficiency -> GetXaxis() -> SetTitleOffset(1.3);
  _hEfficiency -> GetXaxis() -> SetRangeUser(x1, x2);
  _hEfficiency -> GetYaxis() -> SetTitle(sYeff);
  _hEfficiency -> GetYaxis() -> SetTitleFont(42);
  _hEfficiency -> GetYaxis() -> SetTitleSize(0.04);
  _hEfficiency -> GetYaxis() -> SetTitleOffset(1.3);
  _hEfficiency -> GetYaxis() -> CenterTitle(true);
  _hEfficiency -> GetYaxis() -> SetLabelFont(42);
  _hEfficiency -> GetYaxis() -> SetLabelSize(0.04);
  _hEfficiency -> GetYaxis() -> SetRangeUser(y1lo, y2lo);
  _hEfficiency -> SetTitleFont(42);


  // create response profile
  TH1D *hResProfile = _hResponse -> ProfileX("hResProfile", 0, -1, "S");
  hResProfile -> GetXaxis() -> SetTitle(sXres);
  hResProfile -> GetXaxis() -> SetTitleFont(42);
  hResProfile -> GetXaxis() -> SetTitleSize(0.04);
  hResProfile -> GetXaxis() -> CenterTitle(true);
  hResProfile -> GetXaxis() -> SetLabelFont(42);
  hResProfile -> GetXaxis() -> SetLabelSize(0.04);
  hResProfile -> GetXaxis() -> SetTitleOffset(1.3);
  hResProfile -> GetXaxis() -> SetRangeUser(x1, x2);
  hResProfile -> GetYaxis() -> SetTitle(sYres);
  hResProfile -> GetYaxis() -> SetTitleFont(42);
  hResProfile -> GetYaxis() -> SetTitleSize(0.04);
  hResProfile -> GetYaxis() -> SetTitleOffset(1.3);
  hResProfile -> GetYaxis() -> CenterTitle(true);
  hResProfile -> GetYaxis() -> SetLabelFont(42);
  hResProfile -> GetYaxis() -> SetLabelSize(0.04);
  hResProfile -> GetYaxis() -> SetRangeUser(x1, x2);
  hResProfile -> SetTitleFont(42);


  // create line
  TLine *lOne = new TLine(x1, 1., x2, 1.);
  lOne -> SetLineStyle(lL);
  lOne -> SetLineColor(cL);


  // draw plots
  TCanvas *cAll = new TCanvas("cAll", "All 4 spectra", 750, 950);
  TPad    *pLoA = new TPad("pLoA", "backfold vs. measured ratio", 0, 0, 1, 0.35);
  TPad    *pUpA = new TPad("pUpA", "all 4 spectra", 0, 0.35, 1, 1);
  // set plot options
  pLoA   -> SetFillStyle(4000);
  pLoA   -> SetFillColor(0);
  pLoA   -> SetBorderMode(0);
  pLoA   -> SetFrameBorderMode(0);
  pLoA   -> SetLeftMargin(0.15);
  pLoA   -> SetRightMargin(0.05);
  pLoA   -> SetTopMargin(0.);
  pLoA   -> SetBottomMargin(0.25);
  pLoA   -> SetGrid(0, 0);
  pLoA   -> SetTickx(1);
  pLoA   -> SetTicky(1);
  pUpA   -> SetFillStyle(4000);
  pUpA   -> SetFillColor(0);
  pUpA   -> SetBorderMode(0);
  pUpA   -> SetFrameBorderMode(0);
  pUpA   -> SetLeftMargin(0.15);
  pUpA   -> SetRightMargin(0.05);
  pUpA   -> SetTopMargin(0.05);
  pUpA   -> SetBottomMargin(0.);
  pUpA   -> SetGrid(0, 0);
  pUpA   -> SetTickx(1);
  pUpA   -> SetTicky(1);
  pUpA   -> SetLogy(1);
  pLoA   -> Draw();
  pUpA   -> Draw();
  // draw histograms
  pLoA   -> cd();
  hLo    -> Draw();
  if (_method == 3) {
    DrawHistogram(_hUnfoldVsPriRatio, "PE2 same", cUP, cUP, cUP, mR, lR, fRU, 1.);
    DrawHistogram(_hBackVsMeasRatio, "PE2 same", cBM, cBM, cBM, mR, lR, fRB, 1.);
  }
  else {
    DrawHistogram(_hUnfoldVsPriRatio, "PE5 same", cUP, cUP, cUP, mR, lR, fRU, 1.);
    DrawHistogram(_hBackVsMeasRatio, "PE5 same", cBM, cBM, cBM, mR, lR, fRB, 1.);
  }
  lUvB   -> Draw();
  lOne   -> Draw();
  pChi2U -> Draw();
  pChi2B -> Draw();
  pUpA   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  DrawHistogram(_hPrior, "PE2 same", cP, cP, cP, mP, lP, fP, 1.);
  if (_method == 3) {
    DrawHistogram(_hUnfolded, "PE2 same", cU, cU, cU, mU, lU, fU, 1.);
    DrawHistogram(_hBackfolded, "PE2 same", cB, cB, cB, mB, lB, fB, 1.);
  }
  else {
    DrawHistogram(_hUnfolded, "PE5 same", cU, cU, cU, mU, lU, fU, 1.);
    DrawHistogram(_hBackfolded, "PE5 same", cB, cB, cB, mB, lB, fB, 1.);
  }
  lAll   -> Draw();
  _label -> Draw();
  cAll   -> Write();
  cAll   -> Close();

  TCanvas *cBvM = new TCanvas("cBvM", "Backfold vs. measured", 750, 950);
  TPad    *pLo1 = new TPad("pLo1", "ratio", 0, 0, 1, 0.35);
  TPad    *pUp1 = new TPad("pUp1", "backfold and measured", 0, 0.35, 1, 1);
  // set plot options
  pLo1   -> SetFillStyle(4000);
  pLo1   -> SetFillColor(0);
  pLo1   -> SetBorderMode(0);
  pLo1   -> SetFrameBorderMode(0);
  pLo1   -> SetLeftMargin(0.15);
  pLo1   -> SetRightMargin(0.05);
  pLo1   -> SetTopMargin(0.);
  pLo1   -> SetBottomMargin(0.25);
  pLo1   -> SetGrid(0, 0);
  pLo1   -> SetTickx(1);
  pLo1   -> SetTicky(1);
  pUp1   -> SetFillStyle(4000);
  pUp1   -> SetFillColor(0);
  pUp1   -> SetBorderMode(0);
  pUp1   -> SetFrameBorderMode(0);
  pUp1   -> SetLeftMargin(0.15);
  pUp1   -> SetRightMargin(0.05);
  pUp1   -> SetTopMargin(0.05);
  pUp1   -> SetBottomMargin(0.);
  pUp1   -> SetGrid(0, 0);
  pUp1   -> SetTickx(1);
  pUp1   -> SetTicky(1);
  pUp1   -> SetLogy(1);
  pLo1   -> Draw();
  pUp1   -> Draw();
  // draw histograms
  pLo1   -> cd();
  hLo    -> Draw();
  if (_method == 3)
    DrawHistogram(_hBackVsMeasRatio, "PE2 same", cBM, cBM, cBM, mR, lR, fRB, 1.);
  else
    DrawHistogram(_hBackVsMeasRatio, "PE5 same", cBM, cBM, cBM, mR, lR, fRB, 1.);
  lOne   -> Draw();
  pChi2B -> Draw();
  pUp1   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  if (_method == 3)
    DrawHistogram(_hBackfolded, "PE2 same", cB, cB, cB, mB, lB, fB, 1.);
  else
    DrawHistogram(_hBackfolded, "PE5 same", cB, cB, cB, mB, lB, fB, 1.);
  lBvM   -> Draw();
  _label -> Draw();
  cBvM   -> Write();
  cBvM   -> Close();

  TCanvas *cPvU = new TCanvas("cPvU", "Prior vs. unfolded", 750, 950);
  TPad    *pLo2 = new TPad("pLo2", "ratio", 0, 0, 1, 0.35);
  TPad    *pUp2 = new TPad("pUp2", "prior and unfolded", 0, 0.35, 1, 1);
  // set plot options
  pLo2   -> SetFillStyle(4000);
  pLo2   -> SetFillColor(0);
  pLo2   -> SetBorderMode(0);
  pLo2   -> SetFrameBorderMode(0);
  pLo2   -> SetLeftMargin(0.15);
  pLo2   -> SetRightMargin(0.05);
  pLo2   -> SetTopMargin(0.);
  pLo2   -> SetBottomMargin(0.25);
  pLo2   -> SetGrid(0, 0);
  pLo2   -> SetTickx(1);
  pLo2   -> SetTicky(1);
  pUp2   -> SetFillStyle(4000);
  pUp2   -> SetFillColor(0);
  pUp2   -> SetBorderMode(0);
  pUp2   -> SetFrameBorderMode(0);
  pUp2   -> SetLeftMargin(0.15);
  pUp2   -> SetRightMargin(0.05);
  pUp2   -> SetTopMargin(0.05);
  pUp2   -> SetBottomMargin(0.);
  pUp2   -> SetGrid(0, 0);
  pUp2   -> SetTickx(1);
  pUp2   -> SetTicky(1);
  pUp2   -> SetLogy(1);
  pLo2   -> Draw();
  pUp2   -> Draw();
  // draw histograms
  pLo2   -> cd();
  hLo    -> Draw();
  if (_method == 3)
    DrawHistogram(_hUnfoldVsPriRatio, "PE2 same", cUP, cUP, cUP, mR, lR, fRU, 1.);
  else
    DrawHistogram(_hUnfoldVsPriRatio, "PE5 same", cUP, cUP, cUP, mR, lR, fRU, 1.);
  lOne   -> Draw();
  pChi2U -> Draw();
  pUp2   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hPrior, "PE2 same", cP, cP, cP, mP, lP, fP, 1.);
  if (_method == 3)
    DrawHistogram(_hUnfolded, "PE2 same", cU, cU, cU, mU, lU, fU, 1.);
  else
    DrawHistogram(_hUnfolded, "PE5 same", cU, cU, cU, mU, lU, fU, 1.);
  lPvU   -> Draw();
  _label -> Draw();
  cPvU   -> Write();
  cPvU   -> Close();

  TCanvas *cSvM = new TCanvas("cSvM", "Smear vs. measured", 750, 950);
  TPad    *pLo3 = new TPad("pLo3", "ratio", 0, 0, 1, 0.35);
  TPad    *pUp3 = new TPad("pUp3", "smeared and measured", 0, 0.35, 1, 1);
  // set plot options
  pLo3   -> SetFillStyle(4000);
  pLo3   -> SetFillColor(0);
  pLo3   -> SetBorderMode(0);
  pLo3   -> SetFrameBorderMode(0);
  pLo3   -> SetLeftMargin(0.15);
  pLo3   -> SetRightMargin(0.05);
  pLo3   -> SetTopMargin(0.);
  pLo3   -> SetBottomMargin(0.25);
  pLo3   -> SetGrid(0, 0);
  pLo3   -> SetTickx(1);
  pLo3   -> SetTicky(1);
  pUp3   -> SetFillStyle(4000);
  pUp3   -> SetFillColor(0);
  pUp3   -> SetBorderMode(0);
  pUp3   -> SetFrameBorderMode(0);
  pUp3   -> SetLeftMargin(0.15);
  pUp3   -> SetRightMargin(0.05);
  pUp3   -> SetTopMargin(0.05);
  pUp3   -> SetBottomMargin(0.);
  pUp3   -> SetGrid(0, 0);
  pUp3   -> SetTickx(1);
  pUp3   -> SetTicky(1);
  pUp3   -> SetLogy(1);
  pLo3   -> Draw();
  pUp3   -> Draw();
  // draw histograms
  pLo3   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hSmearVsMeasRatio, "PE2 same", cSM, cSM, cSM, mR, lR, fRS, 1.);
  lOne   -> Draw();
  pUp3   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  DrawHistogram(_hSmeared, "PE2 same", cS, cS, cS, mS, lS, fS, 1.);
  lSvM   -> Draw();
  _label -> Draw();
  cSvM   -> Write();
  cSvM   -> Close();

  TCanvas *cUvM = new TCanvas("cUvM", "Unfolded vs. measured", 750, 950);
  TPad    *pLo4 = new TPad("pLo4", "ratio", 0, 0, 1, 0.35);
  TPad    *pUp4 = new TPad("pUp4", "unfolded and measured", 0, 0.35, 1, 1);
  // set plot options
  pLo4   -> SetFillStyle(4000);
  pLo4   -> SetFillColor(0);
  pLo4   -> SetBorderMode(0);
  pLo4   -> SetFrameBorderMode(0);
  pLo4   -> SetLeftMargin(0.15);
  pLo4   -> SetRightMargin(0.05);
  pLo4   -> SetTopMargin(0.);
  pLo4   -> SetBottomMargin(0.25);
  pLo4   -> SetGrid(0, 0);
  pLo4   -> SetTickx(1);
  pLo4   -> SetTicky(1);
  pUp4   -> SetFillStyle(4000);
  pUp4   -> SetFillColor(0);
  pUp4   -> SetBorderMode(0);
  pUp4   -> SetFrameBorderMode(0);
  pUp4   -> SetLeftMargin(0.15);
  pUp4   -> SetRightMargin(0.05);
  pUp4   -> SetTopMargin(0.05);
  pUp4   -> SetBottomMargin(0.);
  pUp4   -> SetGrid(0, 0);
  pUp4   -> SetTickx(1);
  pUp4   -> SetTicky(1);
  pUp4   -> SetLogy(1);
  pLo4   -> Draw();
  pUp4   -> Draw();
  // draw histograms
  pLo4   -> cd();
  hLo    -> Draw();
  if (_method == 3)
    DrawHistogram(_hUnfoldVsMeasRatio, "PE2 same", cUM, cUM, cUM, mR, lR, fRU, 1.);
  else
    DrawHistogram(_hUnfoldVsMeasRatio, "PE5 same", cUM, cUM, cUM, mR, lR, fRU, 1.);
  lOne   -> Draw();
  pUp4   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hMeasured, "PE2 same", cM, cM, cM, mM, lM, fM, 1.);
  if (_method == 3)
    DrawHistogram(_hUnfolded, "PE2 same", cU, cU, cU, mU, lU, fU, 1.);
  else
    DrawHistogram(_hUnfolded, "PE5 same", cU, cU, cU, mU, lU, fU, 1.);
  lUvM   -> Draw();
  _label -> Draw();
  cUvM   -> Write();
  cUvM   -> Close();

  TCanvas *cSvP = new TCanvas("cSvP", "Smeared vs. prior", 750, 950);
  TPad    *pLo5 = new TPad("pLo5", "ratio", 0, 0, 1, 0.35);
  TPad    *pUp5 = new TPad("pUp5", "smeared and prior", 0, 0.35, 1, 1);
  // set plot options
  pLo5   -> SetFillStyle(4000);
  pLo5   -> SetFillColor(0);
  pLo5   -> SetBorderMode(0);
  pLo5   -> SetFrameBorderMode(0);
  pLo5   -> SetLeftMargin(0.15);
  pLo5   -> SetRightMargin(0.05);
  pLo5   -> SetTopMargin(0.);
  pLo5   -> SetBottomMargin(0.25);
  pLo5   -> SetGrid(0, 0);
  pLo5   -> SetTickx(1);
  pLo5   -> SetTicky(1);
  pUp5   -> SetFillStyle(4000);
  pUp5   -> SetFillColor(0);
  pUp5   -> SetBorderMode(0);
  pUp5   -> SetFrameBorderMode(0);
  pUp5   -> SetLeftMargin(0.15);
  pUp5   -> SetRightMargin(0.05);
  pUp5   -> SetTopMargin(0.05);
  pUp5   -> SetBottomMargin(0.);
  pUp5   -> SetGrid(0, 0);
  pUp5   -> SetTickx(1);
  pUp5   -> SetTicky(1);
  pUp5   -> SetLogy(1);
  pLo5   -> Draw();
  pUp5   -> Draw();
  // draw histograms
  pLo5   -> cd();
  hLo    -> Draw();
  DrawHistogram(_hSmearVsPriRatio, "PE2 same", cSP, cSP, cSP, mR, lR, fRS, 1.);
  lOne   -> Draw();
  pUp5   -> cd();
  hUp    -> Draw();
  DrawHistogram(_hPrior, "PE2 same", cP, cP, cP, mP, lP, fP, 1.);
  DrawHistogram(_hSmeared, "PE2 same", cS, cS, cS, mS, lS, fS, 1.);
  lSvP   -> Draw();
  _label -> Draw();
  cSvP   -> Write();
  cSvP   -> Close();


  TCanvas *cResponse = new TCanvas("cResponse", "Efficiency and response matrix", 1500, 750);
  TPad    *pResponse = new TPad("pResponse", "response matrix", 0, 0, 0.5, 1);
  TPad    *pEfficiency = new TPad("pEfficiency", "efficiency", 0.5, 0., 1, 1);
  // set plot options
  pResponse   -> SetFillStyle(4000);
  pResponse   -> SetFillColor(0);
  pResponse   -> SetBorderMode(0);
  pResponse   -> SetFrameBorderMode(0);
  pResponse   -> SetLeftMargin(0.15);
  pResponse   -> SetRightMargin(0.15);
  pResponse   -> SetTopMargin(0.05);
  pResponse   -> SetBottomMargin(0.15);
  pResponse   -> SetGrid(0, 0);
  pResponse   -> SetTickx(1);
  pResponse   -> SetTicky(1);
  pResponse   -> SetLogz(1);
  pEfficiency -> SetFillStyle(4000);
  pEfficiency -> SetFillColor(0);
  pEfficiency -> SetBorderMode(0);
  pEfficiency -> SetFrameBorderMode(0);
  pEfficiency -> SetLeftMargin(0.15);
  pEfficiency -> SetRightMargin(0.15);
  pEfficiency -> SetTopMargin(0.05);
  pEfficiency -> SetBottomMargin(0.15);
  pEfficiency -> SetGrid(0, 0);
  pEfficiency -> SetTickx(1);
  pEfficiency -> SetTicky(1);
  pResponse   -> Draw();
  pEfficiency -> Draw();
  // draw histograms
  pResponse   -> cd();
  _hResponse  -> Draw("colz");
  DrawHistogram(hResProfile, "same", cRP, cRP, cRP, mRP, lRP, fRP, 1.);
  pEfficiency -> cd();
  DrawHistogram(_hEfficiency, "PE2", cM, cM, cM, mM, lM, fM, 1.);
  _label      -> Draw();
  cResponse   -> Write();
  cResponse   -> Close();


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
    case 4:
      pTxt.Append("Power");
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
