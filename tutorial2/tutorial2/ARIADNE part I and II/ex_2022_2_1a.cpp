
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

void ex_2022_2_1a()
{
  // Drift time as observable
  RooRealVar x("x", "Drift time", 460, 530, "#mus");
  x.setBins(70); 

  RooRealVar fs("fs", "signal strength", 0.5, -0.5, 1.5);
  
  // Composite model
  RooWorkspace w;
  w.import(x);
  w.factory("Uniform::bkg(x)");
  w.factory("Landau::lan(x, mean[450, 400, 500], gamma[20, 0, 100])");

  w.import(fs);
  w.factory("SUM:model(fs * lan, bkg)");
    
  RooDataHist data("data", "data", x);
    
  // Reading data from .dat file
  std::ifstream is("ariadne_g006_plus_400.dat");
  double val, weight;
   while(!is.eof()) {
        is >> val >> weight;
        x.setVal(val);
        data.set(x, weight);
   }
    
  is.close();
    
  // Save fit results on RooFitResult
  RooFitResult* res = w.pdf("model")->fitTo(data, Save());
    
  // Print fit results on file FitResults.txt
  std::ofstream outfile("ex_2022_2_1a_FitResults.txt");
  res->printTitle(outfile);
  res->printMultiline(outfile, 0, true);
  outfile.close();

  // Plotting datas and fitted model
  RooPlot* frame = x.frame(Title("Fitted drift times"));
  data.plotOn(frame);
  w.pdf("model")->plotOn(frame);
   
  // Signal (Landau) component in red
  w.pdf("model")->plotOn(frame, Components(*w.pdf("lan")), LineColor(kRed), LineStyle(kDashed));
  // Background (Uniform) component in black
  w.pdf("model")->plotOn(frame, Components(*w.pdf("bkg")), LineColor(kBlack), LineStyle(kDashed));
    
  // Plot fit params in box
  w.pdf("model")->paramOn(frame, Label("Fit parameters"));
  
  // Drawing and Printing
  TCanvas *c= new TCanvas("ex_2022_2_1a","ex_2022_2_1a"); 
  frame->Draw();

  // Save image as .png
  c->Print("ex_2022_2_1a.png");
}