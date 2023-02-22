#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"


using namespace RooFit;
using namespace RooStats;

void exam_project(){

   // Define invariant mass as the observable
   RooRealVar invMass("invMass", "M_{inv}", 60, 200, "GeV");
   invMass.setBins(10); // from the plot

   // make a simple signal model.
   RooRealVar mH("mH", "Higgs Mass", 120, 90, 160);
   RooRealVar sigma1("sigma1", "Width of Gaussian", 12., 2, 100);
   RooGaussian sigModel("sigModel", "Signal Model", invMass, mH, sigma1);
 
   // make zjj model.  Just like signal model
   RooRealVar mZ("mZ", "Z Mass", 91.2, 0, 100);
   RooRealVar sigma1_z("sigma1_z", "Width of Gaussian", 10., 6, 100);
   RooGaussian zjjModel("zjjModel", "Z+jets Model", invMass, mZ, sigma1_z);
    
   // make QCD model
   RooRealVar a0("a0", "a0", 0.26, -1, 1);
   RooRealVar a1("a1", "a1", -0.17596, -1, 1);
   RooRealVar a2("a2", "a2", 0.018437, -1, 1);
   //RooRealVar a3("a3", "a3", 0.02, -1, 1);
   RooChebychev qcdModel("qcdModel", "A  Polynomial for QCD", invMass, RooArgList(a0, a1, a2));
 
   // The shape is known
   a0.setConstant(true);
   a1.setConstant(true);
   a2.setConstant(true);
 
   // Spectrum components
   RooRealVar fsig("fsig", "fraction of signal events", .2, 0., 1);
   RooRealVar fzjj("fzjj", "fraction of zjj background events", .4, 0., 1);

   // use fsig as main parameter, so fix this.
   // fsig.setConstant(true); 
   
   // combined model
   // RooAddPdf model("model", "sig+zjj+qcd", RooArgList(sigModel, zjjModel, qcdModel), RooArgList(fsig, fzjj, fqcd)); // not correct
   RooAddPdf model("model", "sig+zjj+qcd", RooArgList(sigModel, zjjModel, qcdModel), RooArgList(fsig, fzjj));

   // Define a workspace and import model
   RooWorkspace wks("wks");
   wks.import(model);

   // Add a toy dataset
   RooDataSet *data = model.generate(invMass, 150);
   wks.import(*data, Rename("data"));
   
   // Define a Model Config for Higgs mass
   ModelConfig mc_mass("mc_mass");
   mc_mass.SetWorkspace(wks);
   mc_mass.SetPdf(*wks.pdf("model"));
   mc_mass.SetParametersOfInterest(*wks.var("mH"));
   mc_mass.SetObservables(*wks.var("invMass"));         
    //floating p are nuisance
    // wks.defineSet("nuisParams","fzjj,fqcd"); 
    wks.defineSet("nuisParams","fsig,fzjj");
   mc_mass.SetNuisanceParameters(*wks.set("nuisParams")); 
   mc_mass.SetSnapshot(mH);

   // Setting constantness
   // Z mass known
   wks.var("mZ")->setConstant(true);
   // Knowledge of resolution
   wks.var("sigma1_z")->setConstant(true);
   // This specific mass point for the signal will be tested
   // wks.var("mH")->setConstant(true);
   // Knowledge of resolution
   wks.var("sigma1")->setConstant(true);

   wks.import(mc_mass);

   ///////////////////////////////////////////////////////////////////////////////
   // PROFILE LIKELIHOOD
   ProfileLikelihoodCalculator plc;
   plc.SetData(*(wks.data("data")));
   plc.SetModel(mc_mass);
   plc.SetConfidenceLevel(0.68);    
   LikelihoodInterval* interval = plc.GetInterval();
   auto poi = static_cast<RooRealVar*>(mc_mass.GetParametersOfInterest()->first());
   double lowerLimit = interval->LowerLimit(*poi);
   double upperLimit = interval->UpperLimit(*poi);
   
   // Plotting result
   TCanvas* c1 = new TCanvas ("c1", "c1");
   LikelihoodIntervalPlot plot1(interval);
   plot1.Draw();
   
   // Print results HT
   cout << "\n>>>> RESULT PROFILE LIKELIHOOD : " << 0.68 * 100 << "% interval on " << poi->GetName() << " is : ["
         << interval->LowerLimit(*poi) << ", " << interval->UpperLimit(*poi) << "]\n " << endl;

   ///////////////////////////////////////////////////////////////////////////////
   // Make plots for the Alternate hypothesis
   fsig.setConstant(false);

   model.fitTo(*data);
   
   // plot sig candidates, full model, and individual components
   TCanvas* c2 = new TCanvas ("c2", "c2");
   RooPlot *frame2 = invMass.frame();
   data->plotOn(frame2);
   model.plotOn(frame2);
   model.plotOn(frame2, Components(sigModel), LineStyle(kDashed), LineColor(kRed));
   model.plotOn(frame2, Components(zjjModel), LineStyle(kDashed), LineColor(kBlack));
   model.plotOn(frame2, Components(qcdModel), LineStyle(kDashed), LineColor(kGreen));
 
   frame2->SetTitle("An example fit to the signal + background model");
   frame2->Draw();

   // Do Fit to the Null hypothesis
   fsig.setVal(0);          
   fsig.setConstant(true); 

   model.fitTo(*data, Save(kTRUE), Minos(kFALSE), Hesse(kFALSE), PrintLevel(-1));
 
   // plot signal candidates with background model and components
   TCanvas* c3 = new TCanvas ("c3", "c3");
   RooPlot *frame3 = invMass.frame();
   data->plotOn(frame3, DataError(RooAbsData::SumW2));
   model.plotOn(frame3);
   model.plotOn(frame3, Components(zjjModel), LineStyle(kDashed), LineColor(kBlack));
   model.plotOn(frame3, Components(qcdModel), LineStyle(kDashed), LineColor(kGreen));
 
   frame3->SetTitle("An example fit to the background-only model");
   frame3->Draw();
   }

