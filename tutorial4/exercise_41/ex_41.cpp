#include "RooDataSet.h"

#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooPlot.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void ex_41(){
    
    // observables (signal = nobs, control = noff)
    RooRealVar nobs{"nobs", "nobs", 0, 50};
    RooRealVar noff{"noff", "noff", 0, 500};

    RooRealVar s{"s", "s", 3, 0, 15}; // p.o.Interest
    RooRealVar b{"b", "b", 2, 0, 10}; // nuisance

    // total n. of exp. events for signal region model
    RooFormulaVar nexp{"nexp", "@0 + @1", RooArgSet{s, b}};

    RooWorkspace w("w");
    w.import(nobs);
    w.import(noff);
    w.import(s);
    w.import(b);
    w.import(nexp);

    // Poisson of (n|s+b)
    w.factory("Poisson::sb_pdf(nobs,nexp)");

    // Auxiliary measurement for background
    // tau is the ratio of means of noff and nobs. Tau = 100
    RooRealVar tau{"tau", "tau", 100};

    // total n. of exp. events for control region model
    RooFormulaVar nexp_off{"nexp_off", "@0 * @1", RooArgSet{tau, b}};
    
    w.import(tau);
    w.import(nexp_off);

    // control region model 
    w.factory("Poisson::b_pdf(noff, nexp_off)");

    // the total model PDF will be the product of the two
    w.factory("PROD:model(sb_pdf, b_pdf)");
    
    // Creating ModelConfig object
    ModelConfig mc("mc", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(s);
    mc.SetObservables(RooArgSet{nobs, noff});
    mc.SetNuisanceParameters(b);
    mc.SetSnapshot(s);

    // make data set with the number of observed events
    RooDataSet data("data", "data", RooArgSet{nobs, noff});
    nobs.setVal(3);
    noff.setVal(90);
    data.add(RooArgSet{nobs, noff});
    
    // import data set in workspace and save it in a file
    w.import(data);
    w.import(mc);
    w.Print();
    w.writeToFile("ex_41.root", true);

    ProfileLikelihoodCalculator plc(data, mc);
    plc.SetConfidenceLevel(0.95); // 95% interval
    LikelihoodInterval *interval = plc.GetInterval();

    RooRealVar *firstPOI = (RooRealVar *)mc.GetParametersOfInterest()->first();

    LikelihoodIntervalPlot plot(interval);
    plot.SetNPoints(50);
    plot.SetRange(interval->LowerLimit(*firstPOI), interval->UpperLimit(*firstPOI));

    plot.Draw();


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