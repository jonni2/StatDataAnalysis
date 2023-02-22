R__LOAD_LIBRARY($HOME/RooUnfold/libRooUnfold.so)
//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldEx()
{
  cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse response (40, -10.0, 10.0);
  
  // Train with a Breit-Wigner, mean 0.3 and width 2.5. or a Gaus (0.7, 3.)
  for (Int_t i= 0; i<100000; i++) {
    //Double_t xt= gRandom->Gaus (0.7, 3.);
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    Double_t x= smear (xt);
    if (x!=cutdummy)
      response.Fill (x, xt);
    else
      response.Miss (xt);
  }
  
  cout << "==================================== TEST =====================================" << endl;
  TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
    hTrue->Fill(xt);
    if (x!=cutdummy) hMeas->Fill(x);
  }

  cout << "==================================== UNFOLD ===================================" << endl;
RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D* hUnfold= (TH1D*) unfold.Hunfold();
  
  
  // plot the response matrix
  auto responseMatrix = response.Mresponse();
  new TCanvas("response matrix");
  responseMatrix.Draw("colz");

  // correctly define and plot its inverse
  auto responseInv = *static_cast<TMatrixD*>(responseMatrix.Clone());
  responseInv.Invert();
  new TCanvas("inverted matrix");
  responseInv.Draw("colz");

  auto model_meas_Rinv = static_cast<TH1D*>(hMeas->Clone("model_meas_Rinv"));
  model_meas_Rinv->Reset();
  for (int ibin =1; ibin < model_meas_Rinv->GetNbinsX()+1; ibin++){
    double matsum = 0;
    for (int jbin =1; jbin < model_meas_Rinv->GetNbinsX()+1; jbin++) {
        matsum += hMeas->GetBinContent(jbin)*responseInv[ibin-1][jbin-1];
        }
    model_meas_Rinv->SetBinContent(ibin,matsum);
  }

  auto c4 = new TCanvas("matrix inversion approach", "matrix inversion approach");
  c4->Divide(2,2);
  c4->cd(1); hTrue->Draw("hist");
  c4->cd(2); hMeas->Draw("pe same");
  model_meas_Rinv->SetLineColor(EColor::kOrange);
  c4->cd(3); model_meas_Rinv->Draw("histsame");
    
    // plot unfolded distribution
    c4->cd(4);
    model_meas_Rinv->Draw("hist");

  

  TCanvas* c1= new TCanvas("canvas","canvas");

  unfold.PrintTable (cout, hTrue);
  hUnfold->Draw();
  hMeas->Draw("SAME");
  hTrue->SetLineColor(8);
  hTrue->Draw("SAME");

  //c1->SaveAs("RooUnfoldExample.pdf");

}
#ifndef __CINT__
int main () { RooUnfoldEx(); return 0; }  // Main program when run stand-alone
#endif