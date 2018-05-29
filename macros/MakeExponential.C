// 'MakeExponential.C'
// Derek Anderson
//
// Take an exponential, put it in a histogram.

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TRandom3.h"

using namespace std;


void MakeExponential() {

  const Int_t    nBin = 1100;
  const Double_t bin1 = -10.;
  const Double_t bin2 = 100.;
  const Double_t bin  = (bin2 - bin1) / nBin;
  TFile *oFile = new TFile("forTesting.root", "recreate");
  TH1D  *hPri  = new TH1D("hPri", "Prior: Exponential distribution, #tau = 3", nBin, bin1, bin2);
  TH1D  *hSmeE = new TH1D("hSmeE", "Smeared prior: flat efficiency applied", nBin, bin1, bin2);
  TH1D  *hSmeP = new TH1D("hSmeP", "Smeared prior: 'p' dependent efficiency applied", nBin, bin1, bin2);
  TH1D  *hSmeS = new TH1D("hSmeS", "Smeared prior: Exponential with TPC-esque smearing and eff.", nBin, bin1, bin2);
  TH1D  *hEffE = new TH1D("hEffE", "Efficiency: flat efficiency applied", nBin, bin1, bin2);
  TH1D  *hEffP = new TH1D("hEffP", "Efficiency: 'p' dependent efficiency applied", nBin, bin1, bin2);
  TH1D  *hEffS = new TH1D("hEffS", "Efficiency: Exponential with TPC-esque smearing and eff.", nBin, bin1, bin2);
  TH2D  *hResE = new TH2D("hResE", "Response: flat efficiency applied", nBin, bin1, bin2, nBin, bin1, bin2);
  TH2D  *hResP = new TH2D("hResP", "Response: 'p' dependent efficiency applied", nBin, bin1, bin2, nBin, bin1, bin2);
  TH2D  *hResS = new TH2D("hResS", "Response: Exponential with TPC-esque smearing and eff.", nBin, bin1, bin2, nBin, bin1, bin2);

  // set up errors
  hPri  -> Sumw2();
  hSmeE -> Sumw2();
  hSmeP -> Sumw2();
  hSmeS -> Sumw2();
  hEffE -> Sumw2();
  hEffP -> Sumw2();
  hEffS -> Sumw2();
  hResE -> Sumw2();
  hResP -> Sumw2();
  hResS -> Sumw2();


  // fill histograms
  const Int_t    n  = 10000000;
  const Double_t e1 = 0.4;
  const Double_t t  = 3.;
  for (Int_t i = 0; i < n; ++i) {

    Double_t p = gRandom -> Exp(t);
    hPri -> Fill(p);

    // apply flat eff.
    Double_t r1 = gRandom -> Uniform(0., 1.);
    if (r1 < e1) {
      hSmeE -> Fill(p);
      hResE -> Fill(p, p);
    }

    // apply p-dependent eff.
    Double_t e2 = exp(-1.*p/2.);
    Double_t r2 = gRandom -> Uniform(0., 1.);
    if (r2 < e2) {
      hSmeP -> Fill(p);
      hResP -> Fill(p, p);
    }

    // apply TPC-esque smearing
    Double_t r3 = (gRandom -> Gaus(0., (0.01 * p))) + (gRandom -> Gaus(0., 0.017));
    Double_t pS = p * (1.0 + r3);
    Double_t e3 = 0.91 * (1.0 - exp(-4.0 * pS));
    Double_t r4 = gRandom -> Uniform(0., 1.);
    if (r3 < e3) {
      hSmeS -> Fill(pS);
      hResS -> Fill(pS, p);
    }

  }  // end for loop


  // calculate effiencies
  hEffE -> Divide(hSmeE, hPri, 1., 1.);
  hEffP -> Divide(hSmeP, hPri, 1., 1.);
  hEffS -> Divide(hSmeS, hPri, 1., 1.);

  // normalize histograms
  hPri  -> Scale(1. / n);
  hPri  -> Scale(1. / bin);
  hSmeE -> Scale(1. / n);
  hSmeE -> Scale(1. / bin);
  hSmeP -> Scale(1. / n);
  hSmeP -> Scale(1. / bin);
  hSmeS -> Scale(1. / n);
  hSmeS -> Scale(1. / bin);

  const Double_t intE = hResE -> Integral();
  const Double_t intP = hResP -> Integral();
  const Double_t intS = hResS -> Integral();
  hResE -> Scale(1. / intE);
  hResP -> Scale(1. / intP);
  hResS -> Scale(1. / intS);


  // save and close files
  oFile -> cd();
  hPri  -> Write();
  hSmeE -> Write();
  hSmeP -> Write();
  hSmeS -> Write();
  hEffE -> Write();
  hEffP -> Write();
  hEffS -> Write();
  hResE -> Write();
  hResP -> Write();
  hResS -> Write();
  oFile -> Close();

}

// End ------------------------------------------------------------------------
