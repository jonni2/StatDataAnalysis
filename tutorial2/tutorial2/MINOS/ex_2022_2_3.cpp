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
  // PART 1: working with the model --------------------------

  // The reconstructed neutrino energy is the observable
  RooRealVar x("x", "reconstructed energy", 0.5, 14, "GeV");
  x.setBins(23);

  // Load dataset of MINOS observed
  RooDataSet data = *RooDataSet::read("minos_2013_data.dat", x, "v");

  // Load MC SIMULATED dataset
  RooDataSet mc_noosc = *RooDataSet::read("minos_2013_mc.dat", x, "v");
  
  // Creating an Histogram based function for non oscillating neutrinos
  RooDataSet* dd = (RooDataSet*) mc_noosc.reduce(RooArgSet(x));
  RooDataHist* dh_mc_noosc = dd->binnedClone();

  RooHistFunc func_noosc("func_mc_noosc", "No oscillation", x, *dh_mc_noosc, 2);

  // MOdel PDF of the energy distribution
  RooRealVar mixing("mixing", "sin^{2}(2#theta)", 0.95, 0., 1.);
  RooRealVar dm2("dm2", "|#Deltam^{2}|", 0.00232, 0.001, 0.004, "eV^{2}");
  RooRealVar L("L", "distance", 730, "km");

  RooFormulaVar prob_osc("prob_osc", "1 - (mixing * (sin(1.267 * dm2 * L / x))^2)", RooArgList(mixing, dm2, L, x));
  RooGenericPdf model("model", "model", "@0*@1", RooArgSet(prob_osc, func_noosc));

  // Fitting model to data
  model.fitTo(data);

  // Plotting data and model on frame
  RooPlot* frame = x.frame(Title("Energy distribution for #nu_{#mu}->#nu_{#mu} oscillated neutrinos"));
  data.plotOn(frame, LineColor(kRed));
  model.plotOn(frame);
  

  // Drawing and Printing
  TCanvas *c= new TCanvas("ex_2022_2_3_data","ex_2022_2_3_data"); 
  c->cd();
  frame->Draw();

  // Save image as .png
  c->Print("minos_data.png");

  // PART 2: working with LIKELIHOOD --------------------------

  // Save following results on file Results.txt
  std::ofstream file("ex_2022_2_3_Results.txt");

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
  file << "Migrad [A]\ndm2 (eV^2): ";
  dm2.printValue(file);
  file << "\nmixing: ";
  mixing.printValue(file);
  file << "\n\n";

  // Disable verbose logging
  m.setVerbose(kFALSE);
 
  // Run HESSE to calculate errors from d2L/dp2
  m.hesse();

  // Print value (and error) of (dm2, mixing) parameter, that reflects
  // value and error back propagated from MINUIT. [B]
  file << "Hesse [B]:\ndm2 (eV^2): ";
  dm2.printValue(file);
  file << "\nmixing: ";
  mixing.printValue(file);
  file << "\n\n";

  // Run MINOS on dm2 parameter only
  m.minos(dm2);

  // Print value (and error) of dm2 parameter, that reflects
  // value and error back propagated from MINUIT. [C]
  file << "Minos [C]:\ndm2 (eV^2): ";
  dm2.printValue(file);
  file << "\n\n";

  // Save a snapshot of the fit result. [D]
  RooFitResult *res = m.save();
  res->Print("v");
  res->printMultiline(file, 0, kTRUE);
  file.close();
 
  // Make contour plot of dm2 vs mixing at 1,2,3 sigma
  RooPlot *frame_contour = m.contour(mixing, dm2, 1, 2, 3);
  frame_contour->SetTitle("Minuit contour plot, |#Deltam^{2}| vs sin^{2}(2#theta)");

  // Drawing the plot
  TCanvas* c2= new TCanvas("ex_2022_2_3_likelihood","ex_2022_2_3_likelihood"); 
  c2->cd();
  frame_contour->Draw();

  // Save image as .png
  c2->Print("minos_likelihood.png");
 }