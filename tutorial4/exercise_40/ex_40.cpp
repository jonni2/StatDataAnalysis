#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

using namespace RooFit;
using namespace std;
using namespace RooStats;

void ex_40(){

    RooWorkspace w("w");

    // make Poisson model * Gaussian constraint
    w.factory("sum:nexp(s[3,0.1,15],b[1,0,10])");

    // Poisson of (n | s+b)
    w.factory("Poisson:pdf(nobs[0,50],nexp)");

    // Gaussian constraint to express uncertainty in b
    w.factory("Gaussian:constraint(b0[1,0,10],b,sigmab[0.2])");

    // the total model p.d.f will be the product of the two
    w.factory("PROD:model(pdf,constraint)");
    
    // w.var("b0")->setRange();
    w.var("b0")->setConstant(true);

    // Creating ModelConfig object
    ModelConfig mc("mc",&w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("s"));
    mc.SetObservables(*w.var("nobs"));
    mc.SetNuisanceParameters(*w.var("b"));

    // need now to set the global observable
    mc.SetGlobalObservables(*w.var("b0"));

    // this is needed for the hypothesis tests
    mc.SetSnapshot(*w.var("s"));

    // import model in the workspace 
    w.import(mc);

////////////////////////////////////////////////////////////////////////////////////////
    
    // make data set with the number of observed events
    RooDataSet data("data","data", *w.var("nobs"));
    w.var("nobs")->setVal(3);
    data.add(*w.var("nobs") );
    
    // import data set in workspace and save it in a file
    w.import(data);
    w.writeToFile("CountingModel.root", true);

    // Compute upper limits on s with 95% CL
    // PLH case
    ProfileLikelihoodCalculator plc(data, mc);
    plc.SetConfidenceLevel(0.95); // 95% interval
    LikelihoodInterval *interval = plc.GetInterval();

    RooRealVar *firstPOI = (RooRealVar *)mc.GetParametersOfInterest()->first();

    LikelihoodIntervalPlot plot(interval);
    plot.SetNPoints(50);
    plot.SetRange(interval->LowerLimit(*firstPOI), interval->UpperLimit(*firstPOI));

    plot.Draw();

   // FC case
   FeldmanCousins fc(data, mc);
   fc.SetConfidenceLevel(0.95); 
   fc.UseAdaptiveSampling(true); 
   fc.SetNBins(10);             
   fc.CreateConfBelt(true);      

   if (!mc.GetPdf()->canBeExtended()) {
      if (data.numEntries() == 1)
         fc.FluctuateNumDataEntries(false);
      else
         cout << "Not sure what to do about this model" << endl;
   }

   PointSetInterval *interval2 = fc.GetInterval();
   ConfidenceBelt *belt = fc.GetConfidenceBelt();

   cout << "\n>>>> RESULT PROFILE LIKELIHOOD : " << 0.95 * 100 << "% interval on " << firstPOI->GetName() << " is : ["
         << interval->LowerLimit(*firstPOI) << ", " << interval->UpperLimit(*firstPOI) << "]\n " << endl;

   cout << "\n>>>> RESULT NEYMANN-COUSINS : "<< "95% interval on " << firstPOI->GetName() << " is : [" << interval2->LowerLimit(*firstPOI) << ", "
        << interval2->UpperLimit(*firstPOI) << "] " << endl;
   }

