
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

void ex_2022_2_3()
{
  // The reconstructed neutrino energy is the observable
  RooRealVar x("x", "reconstructed energy", 0.5, 14, "GeV");

  // Load dataset of MINOS observed
  RooDataSet data = *RooDataSet::read("minos_2013_data.dat", x, "v");

  // Load dataset of non oscillated
  RooDataSet mc_noosc = *RooDataSet::read("minos_2013_mc.dat", x, "v");
  
  // Creating an Histogram based function for non oscillating neutrinos
  RooDataSet* dd = (RooDataSet*) mc_noosc.reduce(RooArgSet(x));
  RooDataHist* dh_mc_noosc = dd->binnedClone();

  RooHistFunc func_noosc {"func_mc_noosc", "No oscillation", x, *dh_mc_noosc, 2};

  // MOdel PDF of the energy distribution
  RooRealVar mixing("mixing", "( sin(2 * theta) )^2", 0, 1);
  RooRealVar dm2("dm2", "dm2", 0.5, 14);
  RooRealVar distance("distance", "distance L", 730, "km");

  RooFormulaVar prob_osc("prob_osc", "1 - (mixing * (sin(1.267 * dm2 * distance / x))^2)", RooArgList(mixing, dm2, distance, x));
  RooGenericPdf model("model", "model", "@0*@1", RooArgSet(prob_osc, func_noosc));

  // Fitting model to data
  model.fitTo(data);

  // Plotting data and model on frame
  RooPlot* frame = x.frame(Title("Energy distribution for oscillated neutrinos"));
  data.plotOn(frame, LineColor(kRed));
  model.plotOn(frame);

  // Drawing and Printing
  TCanvas *c= new TCanvas("ex_2022_2_3","ex_2022_2_3"); 
  frame->Draw();

  // Save image as .png
  //c->Print("minos_data.png");

  // Costruct function object representing -log(L)
  RooAbsReal* nll = pdf.createNLL(data);
  
  // Create a Minuit interface object
  RooMinuit m(*nll);

  m.setVerbose(kTRUE);
  m.migrad();
  dm2.Print();
  mixing.Print();
  m.setVerbose(kFALSE);
 
  m.hesse();
 }
