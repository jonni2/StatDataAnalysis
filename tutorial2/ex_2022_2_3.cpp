
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
  RooRealVar y("y", "energy_noosc", 0.5, 14, "GeV");
  RooDataSet mc_noosc = *RooDataSet::read("minos_2013_mc.dat", y, "v");
  
  // Creating an Histogram based function for non oscillating neutrinos
  RooDataSet* dd = (RooDataSet*) mc_noosc.reduce(RooArgSet(y));
  RooDataHist* dh_mc_noosc = dd->binnedClone();

  RooHistFunc func_noosc {"func_mc_noosc", "No oscillation", y, *dh_mc_noosc, 2};

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
  TCanvas *c= new TCanvas("ex_2022_2_3_data","ex_2022_2_3_data"); 
  frame->Draw();

  // Print folowing results on file Results.txt
  std::ofstream outfile("ex_2022_2_3_Results.txt");

  // Costruct unbinned likelihood of model w.r.t. data 
  RooAbsReal* nll = model.createNLL(data);
  
  // Create a MINUIT interface object
  RooMinimizer m(*nll);
  
  // Activate verbose logging of MINUIT parameter
  m.setVerbose(kTRUE);

  // Call MIGRAD to minimize the likelihood
  m.migrad();

  // Print values of all parameters, that reflect values (and error estimates)
  // that are back propagated from MINUIT. [A]
  dm2.Print();
  mixing.Print();

  // Disable verbose logging
  m.setVerbose(kFALSE);
 
  // Run HESSE to calculate errors from d2L/dp2
  m.hesse();

  // Print value (and error) of (dm2, mixing) parameter, that reflects
  // value and error back propagated from MINUIT. [B]
  dm2.Print();
  mixing.Print();

  // Run MINOS on dm2 parameter only
  m.minos(dm2);

  // Print value (and error) of dm2 parameter, that reflects
  // value and error back propagated from MINUIT. [C]
  dm2.Print();

  // Save a snapshot of the fit result. This object contains the initial
  // fit parameters, the final fit parameters, the complete correlation
  // matrix, the EDM, the minimized FCN , the last MINUIT status code and
  // the number of times the RooFit function object has indicated evaluation
  // problems (e.g. zero probabilities during likelihood evaluation)
  RooFitResult *res = m.save();

  // Make contour plot of dm2 vs mixing at 1,2,3 sigma
  RooPlot *frame_contour = m.contour(mixing, dm2, 1, 2, 3);
  frame_contour->SetTitle("Minuit contour plot");

  // Drawing the plot
  TCanvas *c2= new TCanvas("ex_2022_2_3_likelihood","ex_2022_2_3_likelihood"); 
  frame_contour->Draw();
  
  // Saving snapshot of fit result and cloing outfile. [D]
  res->printTitle(outfile);
  res->printMultiline(outfile, 0, true);
  outfile.close();

  // Save image as .png
  c->Print("minos_data.png");
  c2->Print("minos_likelihood.png");
 }
