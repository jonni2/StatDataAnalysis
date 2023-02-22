
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"

using namespace RooFit;

void ex_2022_2_2()
{
  // Observable definition
  RooRealVar x("x", "m_{#gamma#gamma}", 5080, 5555, "MeV/c^{2}");
  x.setBins(25); // taken from the paper
  
  // Load dataset from rarest_b0_decay.dat
  RooDataSet data = *RooDataSet::read("rarest_b0_decay.dat", x, "v");
  
  // Model creation
  RooWorkspace w;
  w.import(x);

  // Define Background of model
  w.factory("Exponential::bkg(x, r[-0.1, -10, 0])");
  
  // Define gaussian peak of B0 mass
  w.factory("Gaussian::g_B0(x, m_B0[5279, 5260, 5290], s_B0[50, 0.1, 100])");

  // Define gaussian peak of B0s mass
  w.factory("Gaussian::g_B0s(x, m_B0s[5366, 5350, 5380], s_B0s[10, 0.1, 100])");

  // Create composite model: bkg + B0 + B0s
  w.factory("SUM::model(N_bkg[50, 0, 100]*bkg, N_B0[20, 0, 100]*g_B0, N_B0s[10, 0, 50]*g_B0s)");
  
  // Fit model ti data
  w.pdf("model")->fitTo(data);

  // Plot model and data
  RooPlot* frame = x.frame(Title("B0 Invariant Mass"));
  data.plotOn(frame);
  w.pdf("model")->plotOn(frame, Components(*w.pdf("g_B0")), LineColor(kRed));
  w.pdf("model")->plotOn(frame, Components(*w.pdf("g_B0s")), LineColor(kGreen));
  w.pdf("model")->plotOn(frame, Components(*w.pdf("bkg")), LineColor(kBlack));
  w.pdf("model")->plotOn(frame);

  // Drawing all into a unique Canva
  TCanvas *c = new TCanvas("ex_2022_2_2", "ex_2022_2_2", 900, 900);
  c->Divide(1,3);
  c->cd(1);
  frame->Draw();

  // Residuals of LAST plotted hist and LAST plotted curve on frame without specifying arguments
  RooHist* h_res = frame->residHist();
  RooHist* h_pull = frame->pullHist();

  // New frames for residuals and pull distributions
  RooPlot* frame_res = x.frame(Title("Residuals distribution"));
  RooPlot* frame_pul = x.frame(Title("Pulls distribution"));
  
  // Adding distribution to frame
  frame_res->addPlotable(h_res, "P");
  frame_pul->addPlotable(h_pull, "P");

  // Drawing frames
  c->cd(2);
  frame_res->Draw();
  c->cd(3);
  frame_pul->Draw();

  // Saving into a .png file
  c->SaveAs("ex_2022_2_2.png");
}