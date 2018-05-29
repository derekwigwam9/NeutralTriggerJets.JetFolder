// 'MakePtSpectrum.C'
// Derek Anderson
//
// Make a pT spectrum from particle and detector level tracks

#include <cmath>
#include <cassert>
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;



Double_t Levy(Double_t *x, Double_t *p) {

  Double_t tau = TMath::TwoPi();
  Double_t pT  = x[0];
  Double_t b   = p[0];
  Double_t n   = p[1];
  Double_t t   = p[2];
  Double_t m   = 0.14;
  Double_t mT  = sqrt(pT*pT + m*m);

  Double_t num = tau * b * pT;
  Double_t arg = 1 + (mT - m) / (n * t);
  Double_t den = pow(arg, n);
  Double_t lev = num / den;
  return lev;

}  // end 'Levy(Double_t*, Double_t*)'



void MakePtSpectrum() {

  // global constants
  const Int_t    nTmax = 200;
  const Double_t eTmin = 9.;
  const Double_t eTmax = 100.;
  const Double_t pTmin = 0.2;
  const Double_t pTmax = 20.;
  const Double_t hMax  = 1.;
  const Double_t pi    = TMath::Pi();

  // input / output filenames
  const TString  in("../PythiaData/Pythia20.g2.root");
  const TString  out("Check.root");


  // open output / input files
  TTree *ParTree;
  TTree *DetTree;
  TFile *oFile = new TFile(out, "recreate");
  TFile *iFile = (TFile*) gROOT -> GetListOfFiles() -> FindObject(in);
  if (!iFile) {
    iFile = new TFile(in);
    iFile -> GetObject("ParTree", ParTree);
    iFile -> GetObject("DetTree", DetTree);
  }
  else {
    cout << "PANIC: input file could not be opened!" << endl;
    assert(iFile);
  }


  // create histograms
  const Int_t    nBin  = 1100;
  const Double_t bin1  = -10.;
  const Double_t bin2  = 100.;
  const Double_t bin   = (bin2 - bin1) / nBin;
  TH1D *hPtPar = new TH1D("hPtPar", "p_{T}^{trk}(particle)", nBin, bin1, bin2);
  TH1D *hPtDet = new TH1D("hPtDet", "p_{T}^{trk}(detector)", nBin, bin1, bin2);
  hPtPar -> Sumw2();
  hPtDet -> Sumw2();


  // declare particle leaves
  Int_t    Pvents_num;
  Int_t    Pvents_Process;
  Int_t    Pvents_refmult;
  Int_t    Pvents_refPos;
  Int_t    Pvents_refNeg;
  Int_t    Pvents_runId;
  Double_t Pvents_MagF;
  Int_t    Pvents_nVertex;
  Float_t  Pvents_rankV;
  Double_t Pvents_pvpdz;
  Double_t Pvents_BBCz;
  Double_t Pvents_ZDCz;
  Double_t Pvents_primVx;
  Double_t Pvents_primVy;
  Double_t Pvents_primVz;
  Int_t    Pvents_TrigId;
  Int_t    Pvents_TrigStat;
  Int_t    Pvents_TrigIndex;
  Int_t    Pvents_PartonId_NS;
  Int_t    Pvents_PartonId_AS;
  Int_t    Pvents_PartonStat_NS;
  Int_t    Pvents_PartonStat_AS;
  Float_t  Pvents_PartonEta_NS;
  Float_t  Pvents_PartonEta_AS;
  Float_t  Pvents_PartonPhi_NS;
  Float_t  Pvents_PartonPhi_AS;
  Float_t  Pvents_PartonE_NS;
  Float_t  Pvents_PartonE_AS;
  Float_t  Pvents_PartonEt_NS;
  Float_t  Pvents_PartonEt_AS;
  Float_t  Pvents_tsp;
  Float_t  Pvents_iso;
  Int_t    Pvents_Twr_didT;
  Int_t    Pvents_Twr_adc11;
  Float_t  Pvents_Twr_eneT0;
  Float_t  Pvents_Twr_eT;
  Float_t  Pvents_Twr_ENET0;
  Float_t  Pvents_Twr_phT;
  Float_t  Pvents_Twr_PTower;
  Float_t  Pvents_Twr_pidTower;
  Float_t  Pvents_Twr_moduleT;
  Float_t  Pvents_Clust_EneT0;
  Float_t  Pvents_Clust_EneR;
  Float_t  Pvents_Clust_EneH;
  Float_t  Pvents_Clust_rapv1;
  Float_t  Pvents_Clust_etav1;
  Float_t  Pvents_Clust_phiv1;
  Float_t  Pvents_Estrp_en01;
  Float_t  Pvents_Estrp_en02;
  Float_t  Pvents_Estrp_en03;
  Float_t  Pvents_Estrp_en0;
  Float_t  Pvents_Estrp_en1;
  Float_t  Pvents_Estrp_en2;
  Float_t  Pvents_Estrp_en3;
  Float_t  Pvents_Estrp_en4;
  Float_t  Pvents_Estrp_en5;
  Float_t  Pvents_Estrp_en6;
  Float_t  Pvents_Estrp_en7;
  Float_t  Pvents_Estrp_en8;
  Float_t  Pvents_Estrp_en9;
  Float_t  Pvents_Estrp_en10;
  Float_t  Pvents_Estrp_en11;
  Float_t  Pvents_Estrp_en12;
  Float_t  Pvents_Estrp_en13;
  Float_t  Pvents_Estrp_en14;
  Float_t  Pvents_Estrp_en15;
  Int_t    Pvents_Twr_didE;
  Float_t  Pvents_Pstrip_enp01;
  Float_t  Pvents_Pstrip_enp02;
  Float_t  Pvents_Pstrip_enp03;
  Float_t  Pvents_Pstrip_enp0;
  Float_t  Pvents_Pstrip_enp1;
  Float_t  Pvents_Pstrip_enp2;
  Float_t  Pvents_Pstrip_enp3;
  Float_t  Pvents_Pstrip_enp4;
  Float_t  Pvents_Pstrip_enp5;
  Float_t  Pvents_Pstrip_enp6;
  Float_t  Pvents_Pstrip_enp7;
  Float_t  Pvents_Pstrip_enp8;
  Float_t  Pvents_Pstrip_enp9;
  Float_t  Pvents_Pstrip_enp10;
  Float_t  Pvents_Pstrip_enp11;
  Float_t  Pvents_Pstrip_enp12;
  Float_t  Pvents_Pstrip_enp13;
  Float_t  Pvents_Pstrip_enp14;
  Float_t  Pvents_Pstrip_enp15;
  Float_t  Pvents_clust_Ennq1;
  Float_t  Pvents_clust_Ennq20;
  Float_t  Pvents_clust_Ennq19;
  Float_t  Pvents_clust_Enpq1;
  Float_t  Pvents_clust_Enpq20;
  Float_t  Pvents_clust_Enpq19;
  Float_t  Pvents_clust_Enpq21;
  Int_t    Pvents_noOfprimaryTrks;
  Float_t  pTracks_pT[nTmax];
  Float_t  pTracks_px[nTmax];
  Float_t  pTracks_py[nTmax];
  Float_t  pTracks_pz[nTmax];
  Float_t  pTracks_Eta[nTmax];
  Float_t  pTracks_Phi[nTmax];
  Float_t  pTracks_dEdx[nTmax];
  Float_t  pTracks_chrg[nTmax];
  Float_t  pTracks_gdca[nTmax];
  Int_t    pTracks_Fp[nTmax];
  Int_t    pTracks_Ppo[nTmax];
  Float_t  pTracks_nSigPi[nTmax];
  Float_t  pTracks_nSigK[nTmax];
  Float_t  pTracks_nSigP[nTmax];
  Float_t  pTracks_nSigE[nTmax];

  // declare detector leaves
  Int_t    Dvents_num;
  Int_t    Dvents_Process;
  Int_t    Dvents_refmult;
  Int_t    Dvents_refPos;
  Int_t    Dvents_refNeg;
  Int_t    Dvents_runId;
  Double_t Dvents_MagF;
  Int_t    Dvents_nVertex;
  Float_t  Dvents_rankV;
  Double_t Dvents_pvpdz;
  Double_t Dvents_BBCz;
  Double_t Dvents_ZDCz;
  Double_t Dvents_primVx;
  Double_t Dvents_primVy;
  Double_t Dvents_primVz;
  Int_t    Dvents_TrigId;
  Int_t    Dvents_TrigStat;
  Int_t    Dvents_TrigIndex;
  Int_t    Dvents_PartonId_NS;
  Int_t    Dvents_PartonId_AS;
  Int_t    Dvents_PartonStat_NS;
  Int_t    Dvents_PartonStat_AS;
  Float_t  Dvents_PartonEta_NS;
  Float_t  Dvents_PartonEta_AS;
  Float_t  Dvents_PartonPhi_NS;
  Float_t  Dvents_PartonPhi_AS;
  Float_t  Dvents_PartonE_NS;
  Float_t  Dvents_PartonE_AS;
  Float_t  Dvents_PartonEt_NS;
  Float_t  Dvents_PartonEt_AS;
  Float_t  Dvents_tsp;
  Float_t  Dvents_iso;
  Int_t    Dvents_Twr_didT;
  Int_t    Dvents_Twr_adc11;
  Float_t  Dvents_Twr_eneT0;
  Float_t  Dvents_Twr_eT;
  Float_t  Dvents_Twr_ENET0;
  Float_t  Dvents_Twr_phT;
  Float_t  Dvents_Twr_PTower;
  Float_t  Dvents_Twr_pidTower;
  Float_t  Dvents_Twr_moduleT;
  Float_t  Dvents_Clust_EneT0;
  Float_t  Dvents_Clust_EneR;
  Float_t  Dvents_Clust_EneH;
  Float_t  Dvents_Clust_rapv1;
  Float_t  Dvents_Clust_etav1;
  Float_t  Dvents_Clust_phiv1;
  Float_t  Dvents_Estrp_en01;
  Float_t  Dvents_Estrp_en02;
  Float_t  Dvents_Estrp_en03;
  Float_t  Dvents_Estrp_en0;
  Float_t  Dvents_Estrp_en1;
  Float_t  Dvents_Estrp_en2;
  Float_t  Dvents_Estrp_en3;
  Float_t  Dvents_Estrp_en4;
  Float_t  Dvents_Estrp_en5;
  Float_t  Dvents_Estrp_en6;
  Float_t  Dvents_Estrp_en7;
  Float_t  Dvents_Estrp_en8;
  Float_t  Dvents_Estrp_en9;
  Float_t  Dvents_Estrp_en10;
  Float_t  Dvents_Estrp_en11;
  Float_t  Dvents_Estrp_en12;
  Float_t  Dvents_Estrp_en13;
  Float_t  Dvents_Estrp_en14;
  Float_t  Dvents_Estrp_en15;
  Int_t    Dvents_Twr_didE;
  Float_t  Dvents_Pstrip_enp01;
  Float_t  Dvents_Pstrip_enp02;
  Float_t  Dvents_Pstrip_enp03;
  Float_t  Dvents_Pstrip_enp0;
  Float_t  Dvents_Pstrip_enp1;
  Float_t  Dvents_Pstrip_enp2;
  Float_t  Dvents_Pstrip_enp3;
  Float_t  Dvents_Pstrip_enp4;
  Float_t  Dvents_Pstrip_enp5;
  Float_t  Dvents_Pstrip_enp6;
  Float_t  Dvents_Pstrip_enp7;
  Float_t  Dvents_Pstrip_enp8;
  Float_t  Dvents_Pstrip_enp9;
  Float_t  Dvents_Pstrip_enp10;
  Float_t  Dvents_Pstrip_enp11;
  Float_t  Dvents_Pstrip_enp12;
  Float_t  Dvents_Pstrip_enp13;
  Float_t  Dvents_Pstrip_enp14;
  Float_t  Dvents_Pstrip_enp15;
  Float_t  Dvents_clust_Ennq1;
  Float_t  Dvents_clust_Ennq20;
  Float_t  Dvents_clust_Ennq19;
  Float_t  Dvents_clust_Enpq1;
  Float_t  Dvents_clust_Enpq20;
  Float_t  Dvents_clust_Enpq19;
  Float_t  Dvents_clust_Enpq21;
  Int_t    Dvents_noOfprimaryTrks;
  Float_t  dTracks_pT[nTmax];
  Float_t  dTracks_px[nTmax];
  Float_t  dTracks_py[nTmax];
  Float_t  dTracks_pz[nTmax];
  Float_t  dTracks_Eta[nTmax];
  Float_t  dTracks_Phi[nTmax];
  Float_t  dTracks_dEdx[nTmax];
  Float_t  dTracks_chrg[nTmax];
  Float_t  dTracks_gdca[nTmax];
  Int_t    dTracks_Fp[nTmax];
  Int_t    dTracks_Ppo[nTmax];
  Float_t  dTracks_nSigPi[nTmax];
  Float_t  dTracks_nSigK[nTmax];
  Float_t  dTracks_nSigP[nTmax];
  Float_t  dTracks_nSigE[nTmax];

  // set particle branches
  ParTree -> SetBranchAddress("Events_num", &Pvents_num);
  ParTree -> SetBranchAddress("Events_Process", &Pvents_Process);
  ParTree -> SetBranchAddress("Events_refmult", &Pvents_refmult);
  ParTree -> SetBranchAddress("Events_refPos", &Pvents_refPos);
  ParTree -> SetBranchAddress("Events_refNeg", &Pvents_refNeg);
  ParTree -> SetBranchAddress("Events_runId", &Pvents_runId);
  ParTree -> SetBranchAddress("Events_MagF", &Pvents_MagF);
  ParTree -> SetBranchAddress("Events_nVertex", &Pvents_nVertex);
  ParTree -> SetBranchAddress("Events_rankV", &Pvents_rankV);
  ParTree -> SetBranchAddress("Events_pvpdz", &Pvents_pvpdz);
  ParTree -> SetBranchAddress("Events_BBCz", &Pvents_BBCz);
  ParTree -> SetBranchAddress("Events_ZDCz", &Pvents_ZDCz);
  ParTree -> SetBranchAddress("Events_primVx", &Pvents_primVx);
  ParTree -> SetBranchAddress("Events_primVy", &Pvents_primVy);
  ParTree -> SetBranchAddress("Events_primVz", &Pvents_primVz);
  ParTree -> SetBranchAddress("Events_TrigId", &Pvents_TrigId);
  ParTree -> SetBranchAddress("Events_TrigStat", &Pvents_TrigStat);
  ParTree -> SetBranchAddress("Events_TrigIndex", &Pvents_TrigIndex);
  ParTree -> SetBranchAddress("Events_PartonId_NS", &Pvents_PartonId_NS);
  ParTree -> SetBranchAddress("Events_PartonId_AS", &Pvents_PartonId_AS);
  ParTree -> SetBranchAddress("Events_PartonStat_NS", &Pvents_PartonStat_NS);
  ParTree -> SetBranchAddress("Events_PartonStat_AS", &Pvents_PartonStat_AS);
  ParTree -> SetBranchAddress("Events_PartonEta_NS", &Pvents_PartonEta_NS);
  ParTree -> SetBranchAddress("Events_PartonEta_AS", &Pvents_PartonEta_AS);
  ParTree -> SetBranchAddress("Events_PartonPhi_NS", &Pvents_PartonPhi_NS);
  ParTree -> SetBranchAddress("Events_PartonPhi_AS", &Pvents_PartonPhi_AS);
  ParTree -> SetBranchAddress("Events_PartonE_NS", &Pvents_PartonE_NS);
  ParTree -> SetBranchAddress("Events_PartonE_AS", &Pvents_PartonE_AS);
  ParTree -> SetBranchAddress("Events_PartonEt_NS", &Pvents_PartonEt_NS);
  ParTree -> SetBranchAddress("Events_PartonEt_AS", &Pvents_PartonEt_AS);
  ParTree -> SetBranchAddress("Events_tsp", &Pvents_tsp);
  ParTree -> SetBranchAddress("Events_iso", &Pvents_iso);
  ParTree -> SetBranchAddress("Events_Twr_didT", &Pvents_Twr_didT);
  ParTree -> SetBranchAddress("Events_Twr_adc11", &Pvents_Twr_adc11);
  ParTree -> SetBranchAddress("Events_Twr_eneT0", &Pvents_Twr_eneT0);
  ParTree -> SetBranchAddress("Events_Twr_eT", &Pvents_Twr_eT);
  ParTree -> SetBranchAddress("Events_Twr_ENET0", &Pvents_Twr_ENET0);
  ParTree -> SetBranchAddress("Events_Twr_phT", &Pvents_Twr_phT);
  ParTree -> SetBranchAddress("Events_Twr_PTower", &Pvents_Twr_PTower);
  ParTree -> SetBranchAddress("Events_Twr_pidTower", &Pvents_Twr_pidTower);
  ParTree -> SetBranchAddress("Events_Twr_moduleT", &Pvents_Twr_moduleT);
  ParTree -> SetBranchAddress("Events_Clust_EneT0", &Pvents_Clust_EneT0);
  ParTree -> SetBranchAddress("Events_Clust_EneR", &Pvents_Clust_EneR);
  ParTree -> SetBranchAddress("Events_Clust_EneH", &Pvents_Clust_EneH);
  ParTree -> SetBranchAddress("Events_Clust_rapv1", &Pvents_Clust_rapv1);
  ParTree -> SetBranchAddress("Events_Clust_etav1", &Pvents_Clust_etav1);
  ParTree -> SetBranchAddress("Events_Clust_phiv1", &Pvents_Clust_phiv1);
  ParTree -> SetBranchAddress("Events_Estrp_en01", &Pvents_Estrp_en01);
  ParTree -> SetBranchAddress("Events_Estrp_en02", &Pvents_Estrp_en02);
  ParTree -> SetBranchAddress("Events_Estrp_en03", &Pvents_Estrp_en03);
  ParTree -> SetBranchAddress("Events_Estrp_en0", &Pvents_Estrp_en0);
  ParTree -> SetBranchAddress("Events_Estrp_en1", &Pvents_Estrp_en1);
  ParTree -> SetBranchAddress("Events_Estrp_en2", &Pvents_Estrp_en2);
  ParTree -> SetBranchAddress("Events_Estrp_en3", &Pvents_Estrp_en3);
  ParTree -> SetBranchAddress("Events_Estrp_en4", &Pvents_Estrp_en4);
  ParTree -> SetBranchAddress("Events_Estrp_en5", &Pvents_Estrp_en5);
  ParTree -> SetBranchAddress("Events_Estrp_en6", &Pvents_Estrp_en6);
  ParTree -> SetBranchAddress("Events_Estrp_en7", &Pvents_Estrp_en7);
  ParTree -> SetBranchAddress("Events_Estrp_en8", &Pvents_Estrp_en8);
  ParTree -> SetBranchAddress("Events_Estrp_en9", &Pvents_Estrp_en9);
  ParTree -> SetBranchAddress("Events_Estrp_en10", &Pvents_Estrp_en10);
  ParTree -> SetBranchAddress("Events_Estrp_en11", &Pvents_Estrp_en11);
  ParTree -> SetBranchAddress("Events_Estrp_en12", &Pvents_Estrp_en12);
  ParTree -> SetBranchAddress("Events_Estrp_en13", &Pvents_Estrp_en13);
  ParTree -> SetBranchAddress("Events_Estrp_en14", &Pvents_Estrp_en14);
  ParTree -> SetBranchAddress("Events_Estrp_en15", &Pvents_Estrp_en15);
  ParTree -> SetBranchAddress("Events_Twr_didE", &Pvents_Twr_didE);
  ParTree -> SetBranchAddress("Events_Pstrip_enp01", &Pvents_Pstrip_enp01);
  ParTree -> SetBranchAddress("Events_Pstrip_enp02", &Pvents_Pstrip_enp02);
  ParTree -> SetBranchAddress("Events_Pstrip_enp03", &Pvents_Pstrip_enp03);
  ParTree -> SetBranchAddress("Events_Pstrip_enp0", &Pvents_Pstrip_enp0);
  ParTree -> SetBranchAddress("Events_Pstrip_enp1", &Pvents_Pstrip_enp1);
  ParTree -> SetBranchAddress("Events_Pstrip_enp2", &Pvents_Pstrip_enp2);
  ParTree -> SetBranchAddress("Events_Pstrip_enp3", &Pvents_Pstrip_enp3);
  ParTree -> SetBranchAddress("Events_Pstrip_enp4", &Pvents_Pstrip_enp4);
  ParTree -> SetBranchAddress("Events_Pstrip_enp5", &Pvents_Pstrip_enp5);
  ParTree -> SetBranchAddress("Events_Pstrip_enp6", &Pvents_Pstrip_enp6);
  ParTree -> SetBranchAddress("Events_Pstrip_enp7", &Pvents_Pstrip_enp7);
  ParTree -> SetBranchAddress("Events_Pstrip_enp8", &Pvents_Pstrip_enp8);
  ParTree -> SetBranchAddress("Events_Pstrip_enp9", &Pvents_Pstrip_enp9);
  ParTree -> SetBranchAddress("Events_Pstrip_enp10", &Pvents_Pstrip_enp10);
  ParTree -> SetBranchAddress("Events_Pstrip_enp11", &Pvents_Pstrip_enp11);
  ParTree -> SetBranchAddress("Events_Pstrip_enp12", &Pvents_Pstrip_enp12);
  ParTree -> SetBranchAddress("Events_Pstrip_enp13", &Pvents_Pstrip_enp13);
  ParTree -> SetBranchAddress("Events_Pstrip_enp14", &Pvents_Pstrip_enp14);
  ParTree -> SetBranchAddress("Events_Pstrip_enp15", &Pvents_Pstrip_enp15);
  ParTree -> SetBranchAddress("Events_clust_Ennq1", &Pvents_clust_Ennq1);
  ParTree -> SetBranchAddress("Events_clust_Ennq20", &Pvents_clust_Ennq20);
  ParTree -> SetBranchAddress("Events_clust_Ennq19", &Pvents_clust_Ennq19);
  ParTree -> SetBranchAddress("Events_clust_Enpq1", &Pvents_clust_Enpq1);
  ParTree -> SetBranchAddress("Events_clust_Enpq20", &Pvents_clust_Enpq20);
  ParTree -> SetBranchAddress("Events_clust_Enpq19", &Pvents_clust_Enpq19);
  ParTree -> SetBranchAddress("Events_clust_Enpq21", &Pvents_clust_Enpq21);
  ParTree -> SetBranchAddress("Events_noOfprimaryTrks", &Pvents_noOfprimaryTrks);
  ParTree -> SetBranchAddress("pTracks_pT", pTracks_pT);
  ParTree -> SetBranchAddress("pTracks_px", pTracks_px);
  ParTree -> SetBranchAddress("pTracks_py", pTracks_py);
  ParTree -> SetBranchAddress("pTracks_pz", pTracks_pz);
  ParTree -> SetBranchAddress("pTracks_Eta", pTracks_Eta);
  ParTree -> SetBranchAddress("pTracks_Phi", pTracks_Phi);
  ParTree -> SetBranchAddress("pTracks_dEdx", pTracks_dEdx);
  ParTree -> SetBranchAddress("pTracks_chrg", pTracks_chrg);
  ParTree -> SetBranchAddress("pTracks_gdca", pTracks_gdca);
  ParTree -> SetBranchAddress("pTracks_Fp", pTracks_Fp);
  ParTree -> SetBranchAddress("pTracks_Ppo", pTracks_Ppo);
  ParTree -> SetBranchAddress("pTracks_nSigPi", pTracks_nSigPi);
  ParTree -> SetBranchAddress("pTracks_nSigK", pTracks_nSigK);
  ParTree -> SetBranchAddress("pTracks_nSigP", pTracks_nSigP);
  ParTree -> SetBranchAddress("pTracks_nSigE", pTracks_nSigE);

  // set detector branches
  DetTree -> SetBranchAddress("Events_num", &Dvents_num);
  DetTree -> SetBranchAddress("Events_Process", &Dvents_Process);
  DetTree -> SetBranchAddress("Events_refmult", &Dvents_refmult);
  DetTree -> SetBranchAddress("Events_refPos", &Dvents_refPos);
  DetTree -> SetBranchAddress("Events_refNeg", &Dvents_refNeg);
  DetTree -> SetBranchAddress("Events_runId", &Dvents_runId);
  DetTree -> SetBranchAddress("Events_MagF", &Dvents_MagF);
  DetTree -> SetBranchAddress("Events_nVertex", &Dvents_nVertex);
  DetTree -> SetBranchAddress("Events_rankV", &Dvents_rankV);
  DetTree -> SetBranchAddress("Events_pvpdz", &Dvents_pvpdz);
  DetTree -> SetBranchAddress("Events_BBCz", &Dvents_BBCz);
  DetTree -> SetBranchAddress("Events_ZDCz", &Dvents_ZDCz);
  DetTree -> SetBranchAddress("Events_primVx", &Dvents_primVx);
  DetTree -> SetBranchAddress("Events_primVy", &Dvents_primVy);
  DetTree -> SetBranchAddress("Events_primVz", &Dvents_primVz);
  DetTree -> SetBranchAddress("Events_TrigId", &Dvents_TrigId);
  DetTree -> SetBranchAddress("Events_TrigStat", &Dvents_TrigStat);
  DetTree -> SetBranchAddress("Events_TrigIndex", &Dvents_TrigIndex);
  DetTree -> SetBranchAddress("Events_PartonId_NS", &Dvents_PartonId_NS);
  DetTree -> SetBranchAddress("Events_PartonId_AS", &Dvents_PartonId_AS);
  DetTree -> SetBranchAddress("Events_PartonStat_NS", &Dvents_PartonStat_NS);
  DetTree -> SetBranchAddress("Events_PartonStat_AS", &Dvents_PartonStat_AS);
  DetTree -> SetBranchAddress("Events_PartonEta_NS", &Dvents_PartonEta_NS);
  DetTree -> SetBranchAddress("Events_PartonEta_AS", &Dvents_PartonEta_AS);
  DetTree -> SetBranchAddress("Events_PartonPhi_NS", &Dvents_PartonPhi_NS);
  DetTree -> SetBranchAddress("Events_PartonPhi_AS", &Dvents_PartonPhi_AS);
  DetTree -> SetBranchAddress("Events_PartonE_NS", &Dvents_PartonE_NS);
  DetTree -> SetBranchAddress("Events_PartonE_AS", &Dvents_PartonE_AS);
  DetTree -> SetBranchAddress("Events_PartonEt_NS", &Dvents_PartonEt_NS);
  DetTree -> SetBranchAddress("Events_PartonEt_AS", &Dvents_PartonEt_AS);
  DetTree -> SetBranchAddress("Events_tsp", &Dvents_tsp);
  DetTree -> SetBranchAddress("Events_iso", &Dvents_iso);
  DetTree -> SetBranchAddress("Events_Twr_didT", &Dvents_Twr_didT);
  DetTree -> SetBranchAddress("Events_Twr_adc11", &Dvents_Twr_adc11);
  DetTree -> SetBranchAddress("Events_Twr_eneT0", &Dvents_Twr_eneT0);
  DetTree -> SetBranchAddress("Events_Twr_eT", &Dvents_Twr_eT);
  DetTree -> SetBranchAddress("Events_Twr_ENET0", &Dvents_Twr_ENET0);
  DetTree -> SetBranchAddress("Events_Twr_phT", &Dvents_Twr_phT);
  DetTree -> SetBranchAddress("Events_Twr_PTower", &Dvents_Twr_PTower);
  DetTree -> SetBranchAddress("Events_Twr_pidTower", &Dvents_Twr_pidTower);
  DetTree -> SetBranchAddress("Events_Twr_moduleT", &Dvents_Twr_moduleT);
  DetTree -> SetBranchAddress("Events_Clust_EneT0", &Dvents_Clust_EneT0);
  DetTree -> SetBranchAddress("Events_Clust_EneR", &Dvents_Clust_EneR);
  DetTree -> SetBranchAddress("Events_Clust_EneH", &Dvents_Clust_EneH);
  DetTree -> SetBranchAddress("Events_Clust_rapv1", &Dvents_Clust_rapv1);
  DetTree -> SetBranchAddress("Events_Clust_etav1", &Dvents_Clust_etav1);
  DetTree -> SetBranchAddress("Events_Clust_phiv1", &Dvents_Clust_phiv1);
  DetTree -> SetBranchAddress("Events_Estrp_en01", &Dvents_Estrp_en01);
  DetTree -> SetBranchAddress("Events_Estrp_en02", &Dvents_Estrp_en02);
  DetTree -> SetBranchAddress("Events_Estrp_en03", &Dvents_Estrp_en03);
  DetTree -> SetBranchAddress("Events_Estrp_en0", &Dvents_Estrp_en0);
  DetTree -> SetBranchAddress("Events_Estrp_en1", &Dvents_Estrp_en1);
  DetTree -> SetBranchAddress("Events_Estrp_en2", &Dvents_Estrp_en2);
  DetTree -> SetBranchAddress("Events_Estrp_en3", &Dvents_Estrp_en3);
  DetTree -> SetBranchAddress("Events_Estrp_en4", &Dvents_Estrp_en4);
  DetTree -> SetBranchAddress("Events_Estrp_en5", &Dvents_Estrp_en5);
  DetTree -> SetBranchAddress("Events_Estrp_en6", &Dvents_Estrp_en6);
  DetTree -> SetBranchAddress("Events_Estrp_en7", &Dvents_Estrp_en7);
  DetTree -> SetBranchAddress("Events_Estrp_en8", &Dvents_Estrp_en8);
  DetTree -> SetBranchAddress("Events_Estrp_en9", &Dvents_Estrp_en9);
  DetTree -> SetBranchAddress("Events_Estrp_en10", &Dvents_Estrp_en10);
  DetTree -> SetBranchAddress("Events_Estrp_en11", &Dvents_Estrp_en11);
  DetTree -> SetBranchAddress("Events_Estrp_en12", &Dvents_Estrp_en12);
  DetTree -> SetBranchAddress("Events_Estrp_en13", &Dvents_Estrp_en13);
  DetTree -> SetBranchAddress("Events_Estrp_en14", &Dvents_Estrp_en14);
  DetTree -> SetBranchAddress("Events_Estrp_en15", &Dvents_Estrp_en15);
  DetTree -> SetBranchAddress("Events_Twr_didE", &Dvents_Twr_didE);
  DetTree -> SetBranchAddress("Events_Pstrip_enp01", &Dvents_Pstrip_enp01);
  DetTree -> SetBranchAddress("Events_Pstrip_enp02", &Dvents_Pstrip_enp02);
  DetTree -> SetBranchAddress("Events_Pstrip_enp03", &Dvents_Pstrip_enp03);
  DetTree -> SetBranchAddress("Events_Pstrip_enp0", &Dvents_Pstrip_enp0);
  DetTree -> SetBranchAddress("Events_Pstrip_enp1", &Dvents_Pstrip_enp1);
  DetTree -> SetBranchAddress("Events_Pstrip_enp2", &Dvents_Pstrip_enp2);
  DetTree -> SetBranchAddress("Events_Pstrip_enp3", &Dvents_Pstrip_enp3);
  DetTree -> SetBranchAddress("Events_Pstrip_enp4", &Dvents_Pstrip_enp4);
  DetTree -> SetBranchAddress("Events_Pstrip_enp5", &Dvents_Pstrip_enp5);
  DetTree -> SetBranchAddress("Events_Pstrip_enp6", &Dvents_Pstrip_enp6);
  DetTree -> SetBranchAddress("Events_Pstrip_enp7", &Dvents_Pstrip_enp7);
  DetTree -> SetBranchAddress("Events_Pstrip_enp8", &Dvents_Pstrip_enp8);
  DetTree -> SetBranchAddress("Events_Pstrip_enp9", &Dvents_Pstrip_enp9);
  DetTree -> SetBranchAddress("Events_Pstrip_enp10", &Dvents_Pstrip_enp10);
  DetTree -> SetBranchAddress("Events_Pstrip_enp11", &Dvents_Pstrip_enp11);
  DetTree -> SetBranchAddress("Events_Pstrip_enp12", &Dvents_Pstrip_enp12);
  DetTree -> SetBranchAddress("Events_Pstrip_enp13", &Dvents_Pstrip_enp13);
  DetTree -> SetBranchAddress("Events_Pstrip_enp14", &Dvents_Pstrip_enp14);
  DetTree -> SetBranchAddress("Events_Pstrip_enp15", &Dvents_Pstrip_enp15);
  DetTree -> SetBranchAddress("Events_clust_Ennq1", &Dvents_clust_Ennq1);
  DetTree -> SetBranchAddress("Events_clust_Ennq20", &Dvents_clust_Ennq20);
  DetTree -> SetBranchAddress("Events_clust_Ennq19", &Dvents_clust_Ennq19);
  DetTree -> SetBranchAddress("Events_clust_Enpq1", &Dvents_clust_Enpq1);
  DetTree -> SetBranchAddress("Events_clust_Enpq20", &Dvents_clust_Enpq20);
  DetTree -> SetBranchAddress("Events_clust_Enpq19", &Dvents_clust_Enpq19);
  DetTree -> SetBranchAddress("Events_clust_Enpq21", &Dvents_clust_Enpq21);
  DetTree -> SetBranchAddress("Events_noOfprimaryTrks", &Dvents_noOfprimaryTrks);
  DetTree -> SetBranchAddress("pTracks_pT", dTracks_pT);
  DetTree -> SetBranchAddress("pTracks_px", dTracks_px);
  DetTree -> SetBranchAddress("pTracks_py", dTracks_py);
  DetTree -> SetBranchAddress("pTracks_pz", dTracks_pz);
  DetTree -> SetBranchAddress("pTracks_Eta", dTracks_Eta);
  DetTree -> SetBranchAddress("pTracks_Phi", dTracks_Phi);
  DetTree -> SetBranchAddress("pTracks_dEdx", dTracks_dEdx);
  DetTree -> SetBranchAddress("pTracks_chrg", dTracks_chrg);
  DetTree -> SetBranchAddress("pTracks_gdca", dTracks_gdca);
  DetTree -> SetBranchAddress("pTracks_Fp", dTracks_Fp);
  DetTree -> SetBranchAddress("pTracks_Ppo", dTracks_Ppo);
  DetTree -> SetBranchAddress("pTracks_nSigPi", dTracks_nSigPi);
  DetTree -> SetBranchAddress("pTracks_nSigK", dTracks_nSigK);
  DetTree -> SetBranchAddress("pTracks_nSigP", dTracks_nSigP);
  DetTree -> SetBranchAddress("pTracks_nSigE", dTracks_nSigE);


  Int_t nPvt = ParTree -> GetEntries();
  Int_t nDvt = DetTree -> GetEntries();
  assert(nPvt == nDvt);
  cout << "Processing " << nPvt << " events:" << endl;

  // event loop
  Int_t nPbyte = 0;
  Int_t nDbyte = 0;
  Int_t nTrg   = 0;
  for (Int_t i = 0;  i < nPvt; i++) {

     if ((i % 1000) == 0)
       cout << "  " << i << " events processed..." << endl;

     // load entries
     nPbyte += ParTree -> GetEntry(i);
     nDbyte += DetTree -> GetEntry(i);

     // trigger cuts
     Int_t    iTrg  = Pvents_TrigIndex;
     Double_t eTrg  = Pvents_Clust_EneT0;
     Double_t hTrg  = Pvents_Clust_etav1;
     Double_t tTrg  = 2. * atan(exp(-1. * hTrg));
     Double_t eTtrg = eTrg * sin(tTrg);
     if ((eTtrg < eTmin) || (eTtrg > eTmax))
       continue;
     if (abs(hTrg) > hMax)
       continue;
     ++nTrg;


     // particle track loop
     Int_t nPtrk = Pvents_noOfprimaryTrks;
     for (Int_t j = 0; j < nPtrk; j++) {

       // track cuts
       Int_t    iTrk  = pTracks_Ppo[j];
       Double_t pTtrk = pTracks_pT[j];
       Double_t hTrk  = pTracks_Eta[j];
       if (iTrk == iTrg)
         continue;
       if ((pTtrk < pTmin) || (pTtrk > pTmax))
         continue;
       if (abs(hTrk) < hMax)
         continue;
       hPtPar -> Fill(pTtrk);

     }

     // detector track loop
     Int_t nDtrk = Dvents_noOfprimaryTrks;
     for (Int_t j = 0; j < nDtrk; j++) {

       // track cuts
       Int_t    iTrk  = dTracks_Ppo[j];
       Double_t pTtrk = dTracks_pT[j];
       Double_t hTrk  = dTracks_Eta[j];
       if (iTrk == iTrg)
         continue;
       if ((pTtrk < pTmin) || (pTtrk > pTmax))
         continue;
       if (abs(hTrk) < hMax)
         continue;
       hPtDet -> Fill(pTtrk);

     }

  }  // end event loop

  cout << "Events processed: " << nTrg << " accepted." << endl;


  // normalize histograms
  hPtPar -> Scale(1./nTrg);
  hPtPar -> Scale(1./bin);
  hPtDet -> Scale(1./nTrg);
  hPtDet -> Scale(1./bin);


  // fit particle spectra
  Int_t    nP = 3;
  Double_t b  = 1.;
  Double_t n  = 5.8;
  Double_t t  = 0.4;
  //Double_t m  = 0.14; 
  TF1 *fLevy  = new TF1("fLevy", Levy, pTmin, pTmax, nP);
  fLevy  -> SetParameters(b, n, t);
  fLevy  -> SetParNames("B", "N", "T");
  fLevy  -> SetLineColor(kRed);
  fLevy  -> SetLineStyle(2);
  hPtPar -> Fit("fLevy");


  // close files
  oFile  -> cd();
  hPtPar -> Write();
  hPtDet -> Write();
  fLevy  -> Write();
  oFile  -> Close();
  iFile  -> cd();
  iFile  -> Close();

}  // end 'MakePtSpectrum()'

// End ------------------------------------------------------------------------
