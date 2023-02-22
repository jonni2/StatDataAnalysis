
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

void ex_2022_2_1b()
{
  // Time as a variable
  RooRealVar x("x", "dI/dX", 9400, 29656, "ADU/cm");
  x.setBins(49); // taken from the paper
  
  // Define the pdfs for the convolution model
  RooWorkspace w;
  w.import(x);
  w.factory("Landau::lan(x, p0[9000, 1, 16000], p1[500, 1, 10000])");
  w.factory("Gaussian::gaus(x, m[4000, 1, 16000], s[1000, 1, 10000])");
    
  // Define the convolution model
  w.factory("FCONV::model(x, lan, gaus)");
    
  // Define the dataset of input data
  RooDataHist data("data", "data", x);
    
  // Read input data from file
  ifstream is("ariadne_g012.dat");
  double val, weight;
  while(!is.eof()) {
      is >> val >> weight;
      x.setVal(val);
      data.set(x, weight);
  }
  is.close();
    
  // Fit and save results
  RooFitResult* res = w.pdf("model")->fitTo(data, Save());
    
  // Print fit results on a .txt file
  std::ofstream outfile("ex_2022_2_1b_FitResults.txt");
  res->printTitle(outfile);
  res->printMultiline(outfile, 0, true);
  outfile.close();
    
  // Draw data and pdf on frame
  RooPlot* frame = x.frame(Title("Data fitted to Landau-Gaussian convolution"));
  data.plotOn(frame);
  w.pdf("model")->plotOn(frame, LineColor(kRed));

  // Physical (Landau) component in green
  // w.pdf("lan")->plotOn(frame, LineColor(kGreen), LineStyle(kDashed));
  // Detector resolution (Gaussian) component in blue
  // w.pdf("gaus")->plotOn(frame, LineColor(kBlue), LineStyle(kDashed));

  // Plot fit params in box
  w.pdf("model")->paramOn(frame, RooFit::Label("Fit parameters"));
    
  // Drawing and Printing
  TCanvas *c= new TCanvas("ex_2022_2_1b","ex_2022_2_1b"); 
  frame->Draw();

  // Save image as .png
  c->Print("ex2022_2_1b.png");
}