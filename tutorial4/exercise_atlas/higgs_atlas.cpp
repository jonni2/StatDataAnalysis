//READ ONLY FROM LINE 445
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

#include <cassert>

using namespace std;
using namespace RooFit;
using namespace RooStats;

struct HypoTestOptions {

   bool noSystematics = false; // force all systematics to be off (i.e. set all nuisance parameters as constat
   double nToysRatio = 4;      // ratio Ntoys Null/ntoys ALT
   double poiValue = -1;       // change poi snapshot value for S+B model (needed for expected p0 values)
   int printLevel = 0;
   bool generateBinned = false;       // for binned generation
   bool useProof = false;             // use Proof
   bool enableDetailedOutput = false; // for detailed output
};

HypoTestOptions optHT;

void StandardHypoTestDemo(const char *infile = "", const char *workspaceName = "combined",
                          const char *modelSBName = "ModelConfig", const char *modelBName = "",
                          const char *dataName = "obsData", int calcType = 0, /* 0 freq 1 hybrid, 2 asymptotic */
                          int testStatType = 3, /* 0 LEP, 1 TeV, 2 LHC, 3 LHC - one sided*/
                          int ntoys = 5000, bool useNC = false, const char *nuisPriorName = 0)
{

   bool noSystematics = optHT.noSystematics;
   double nToysRatio = optHT.nToysRatio; // ratio Ntoys Null/ntoys ALT
   double poiValue = optHT.poiValue;     // change poi snapshot value for S+B model (needed for expected p0 values)
   int printLevel = optHT.printLevel;
   bool generateBinned = optHT.generateBinned; // for binned generation
   bool useProof = optHT.useProof;             // use Proof
   bool enableDetOutput = optHT.enableDetailedOutput;

   // Other Parameter to pass in tutorial
   // apart from standard for filename, ws, modelconfig and data

   // type = 0 Freq calculator
   // type = 1 Hybrid calculator
   // type = 2 Asymptotic calculator

   // testStatType = 0 LEP
   // = 1 Tevatron
   // = 2 Profile Likelihood
   // = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)

   // ntoys:         number of toys to use

   // useNumberCounting:  set to true when using number counting events

   // nuisPriorName:   name of prior for the nuisance. This is often expressed as constraint term in the global model
   // It is needed only when using the HybridCalculator (type=1)
   // If not given by default the prior pdf from ModelConfig is used.

   // extra options are available as global parameters of the macro. They major ones are:

   // generateBinned       generate binned data sets for toys (default is false) - be careful not to activate with
   // a too large (>=3) number of observables
   // nToyRatio            ratio of S+B/B toys (default is 2)
   // printLevel

   // disable - can cause some problems
   // ToyMCSampler::SetAlwaysUseMultiGen(true);

   SimpleLikelihoodRatioTestStat::SetAlwaysReuseNLL(true);
   ProfileLikelihoodTestStat::SetAlwaysReuseNLL(true);
   RatioOfProfiledLikelihoodsTestStat::SetAlwaysReuseNLL(true);

   // RooRandom::randomGenerator()->SetSeed(0);

   // to change minimizers
   // ~~~{.bash}
   // ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
   // ~~~

   // -------------------------------------------------------
   // First part is just to access a user-defined file
   // or create the standard example file if it doesn't exist
   const char *filename = "";
   if (!strcmp(infile, "")) {
      filename = "results/example_combined_GaussExample_model.root";
      bool fileExist = !gSystem->AccessPathName(filename); // note opposite return code
      // if file does not exists generate with histfactory
      if (!fileExist) {
#ifdef _WIN32
         cout << "HistFactory file cannot be generated on Windows - exit" << endl;
         return;
#endif
         // Normally this would be run on the command line
         cout << "will run standard hist2workspace example" << endl;
         gROOT->ProcessLine(".! prepareHistFactory .");
         gROOT->ProcessLine(".! hist2workspace config/example.xml");
         cout << "\n\n---------------------" << endl;
         cout << "Done creating example input" << endl;
         cout << "---------------------\n\n" << endl;
      }

   } else
      filename = infile;

   // Try to open the file
   TFile *file = TFile::Open(filename);

   // if input file was specified byt not found, quit
   if (!file) {
      cout << "StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl;
      return;
   }

   // -------------------------------------------------------
   // Tutorial starts here
   // -------------------------------------------------------

   // get the workspace out of the file
   RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);
   if (!w) {
      cout << "workspace not found" << endl;
      return;
   }
   w->Print();

   // get the modelConfig out of the file
   ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);

   // get the modelConfig out of the file
   RooAbsData *data = w->data(dataName);

   // make sure ingredients are found
   if (!data || !sbModel) {
      w->Print();
      cout << "data or ModelConfig was not found" << endl;
      return;
   }
   // make b model
   ModelConfig *bModel = (ModelConfig *)w->obj(modelBName);

   // case of no systematics
   // remove nuisance parameters from model
   if (noSystematics) {
      const RooArgSet *nuisPar = sbModel->GetNuisanceParameters();
      if (nuisPar && nuisPar->getSize() > 0) {
         std::cout << "StandardHypoTestInvDemo"
                   << "  -  Switch off all systematics by setting them constant to their initial values" << std::endl;
         RooStats::SetAllConstant(*nuisPar);
      }
      if (bModel) {
         const RooArgSet *bnuisPar = bModel->GetNuisanceParameters();
         if (bnuisPar)
            RooStats::SetAllConstant(*bnuisPar);
      }
   }

   if (!bModel) {
      Info("StandardHypoTestInvDemo", "The background model %s does not exist", modelBName);
      Info("StandardHypoTestInvDemo", "Copy it from ModelConfig %s and set POI to zero", modelSBName);
      bModel = (ModelConfig *)sbModel->Clone();
      bModel->SetName(TString(modelSBName) + TString("B_only"));
      RooRealVar *var = dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
      if (!var)
         return;
      double oldval = var->getVal();
      var->setVal(0);
      // bModel->SetSnapshot( RooArgSet(*var, *w->var("lumi"))  );
      bModel->SetSnapshot(RooArgSet(*var));
      var->setVal(oldval);
   }

   if (!sbModel->GetSnapshot() || poiValue > 0) {
      Info("StandardHypoTestDemo", "Model %s has no snapshot  - make one using model poi", modelSBName);
      RooRealVar *var = dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
      if (!var)
         return;
      double oldval = var->getVal();
      if (poiValue > 0)
         var->setVal(poiValue);
      // sbModel->SetSnapshot( RooArgSet(*var, *w->var("lumi") ) );
      sbModel->SetSnapshot(RooArgSet(*var));
      if (poiValue > 0)
         var->setVal(oldval);
      // sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
   }

   // part 1, hypothesis testing
   SimpleLikelihoodRatioTestStat *slrts = new SimpleLikelihoodRatioTestStat(*bModel->GetPdf(), *sbModel->GetPdf());
   // null parameters must includes snapshot of poi plus the nuisance values
   RooArgSet nullParams(*bModel->GetSnapshot());
   if (bModel->GetNuisanceParameters())
      nullParams.add(*bModel->GetNuisanceParameters());

   slrts->SetNullParameters(nullParams);
   RooArgSet altParams(*sbModel->GetSnapshot());
   if (sbModel->GetNuisanceParameters())
      altParams.add(*sbModel->GetNuisanceParameters());
   slrts->SetAltParameters(altParams);

   ProfileLikelihoodTestStat *profll = new ProfileLikelihoodTestStat(*bModel->GetPdf());

   RatioOfProfiledLikelihoodsTestStat *ropl =
      new RatioOfProfiledLikelihoodsTestStat(*bModel->GetPdf(), *sbModel->GetPdf(), sbModel->GetSnapshot());
   ropl->SetSubtractMLE(false);

   if (testStatType == 3)
      profll->SetOneSidedDiscovery(1);
   profll->SetPrintLevel(printLevel);

   if (enableDetOutput) {
      slrts->EnableDetailedOutput();
      profll->EnableDetailedOutput();
      ropl->EnableDetailedOutput();
   }

   /* profll.SetReuseNLL(mOptimize);*/
   /* slrts.SetReuseNLL(mOptimize);*/
   /* ropl.SetReuseNLL(mOptimize);*/

   AsymptoticCalculator::SetPrintLevel(printLevel);

   HypoTestCalculatorGeneric *hypoCalc = 0;
   // note here Null is B and Alt is S+B
   if (calcType == 0)
      hypoCalc = new FrequentistCalculator(*data, *sbModel, *bModel);
   else if (calcType == 1)
      hypoCalc = new HybridCalculator(*data, *sbModel, *bModel);
   else if (calcType == 2)
      hypoCalc = new AsymptoticCalculator(*data, *sbModel, *bModel);

   if (calcType == 0) {
      ((FrequentistCalculator *)hypoCalc)->SetToys(ntoys, ntoys / nToysRatio);
      if (enableDetOutput)
         ((FrequentistCalculator *)hypoCalc)->StoreFitInfo(true);
   }
   if (calcType == 1) {
      ((HybridCalculator *)hypoCalc)->SetToys(ntoys, ntoys / nToysRatio);
      // n. a. yetif (enableDetOutput) ((HybridCalculator*) hypoCalc)->StoreFitInfo(true);
   }
   if (calcType == 2) {
      if (testStatType == 3)
         ((AsymptoticCalculator *)hypoCalc)->SetOneSidedDiscovery(true);
      if (testStatType != 2 && testStatType != 3)
         Warning("StandardHypoTestDemo",
                 "Only the PL test statistic can be used with AsymptoticCalculator - use by default a two-sided PL");
   }

   // check for nuisance prior pdf in case of nuisance parameters
   if (calcType == 1 && (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters())) {
      RooAbsPdf *nuisPdf = 0;
      if (nuisPriorName)
         nuisPdf = w->pdf(nuisPriorName);
      // use prior defined first in bModel (then in SbModel)
      if (!nuisPdf) {
         Info("StandardHypoTestDemo",
              "No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the   model");
         if (bModel->GetPdf() && bModel->GetObservables())
            nuisPdf = RooStats::MakeNuisancePdf(*bModel, "nuisancePdf_bmodel");
         else
            nuisPdf = RooStats::MakeNuisancePdf(*sbModel, "nuisancePdf_sbmodel");
      }
      if (!nuisPdf) {
         if (bModel->GetPriorPdf()) {
            nuisPdf = bModel->GetPriorPdf();
            Info("StandardHypoTestDemo",
                 "No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",
                 nuisPdf->GetName());
         } else {
            Error("StandardHypoTestDemo", "Cannot run Hybrid calculator because no prior on the nuisance parameter is "
                                          "specified or can be derived");
            return;
         }
      }
      assert(nuisPdf);
      Info("StandardHypoTestDemo", "Using as nuisance Pdf ... ");
      nuisPdf->Print();

      const RooArgSet *nuisParams =
         (bModel->GetNuisanceParameters()) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
      RooArgSet *np = nuisPdf->getObservables(*nuisParams);
      if (np->getSize() == 0) {
         Warning("StandardHypoTestDemo",
                 "Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
      }
      delete np;

      ((HybridCalculator *)hypoCalc)->ForcePriorNuisanceAlt(*nuisPdf);
      ((HybridCalculator *)hypoCalc)->ForcePriorNuisanceNull(*nuisPdf);
   }

   /* hypoCalc->ForcePriorNuisanceAlt(*sbModel->GetPriorPdf());*/
   /* hypoCalc->ForcePriorNuisanceNull(*bModel->GetPriorPdf());*/

   ToyMCSampler *sampler = (ToyMCSampler *)hypoCalc->GetTestStatSampler();

   if (sampler && (calcType == 0 || calcType == 1)) {

      // look if pdf is number counting or extended
      if (sbModel->GetPdf()->canBeExtended()) {
         if (useNC)
            Warning("StandardHypoTestDemo", "Pdf is extended: but number counting flag is set: ignore it ");
      } else {
         // for not extended pdf
         if (!useNC) {
            int nEvents = data->numEntries();
            Info("StandardHypoTestDemo",
                 "Pdf is not extended: number of events to generate taken  from observed data set is %d", nEvents);
            sampler->SetNEventsPerToy(nEvents);
         } else {
            Info("StandardHypoTestDemo", "using a number counting pdf");
            sampler->SetNEventsPerToy(1);
         }
      }

      if (data->isWeighted() && !generateBinned) {
         Info("StandardHypoTestDemo", "Data set is weighted, nentries = %d and sum of weights = %8.1f but toy "
                                      "generation is unbinned - it would be faster to set generateBinned to true\n",
              data->numEntries(), data->sumEntries());
      }
      if (generateBinned)
         sampler->SetGenerateBinned(generateBinned);

      // use PROOF
      if (useProof) {
         ProofConfig pc(*w, 0, "", kFALSE);
         sampler->SetProofConfig(&pc); // enable proof
      }

      // set the test statistic
      if (testStatType == 0)
         sampler->SetTestStatistic(slrts);
      if (testStatType == 1)
         sampler->SetTestStatistic(ropl);
      if (testStatType == 2 || testStatType == 3)
         sampler->SetTestStatistic(profll);
   }

   HypoTestResult *htr = hypoCalc->GetHypoTest();
   htr->SetPValueIsRightTail(true);
   htr->SetBackgroundAsAlt(false);
   htr->Print(); // how to get meaningful CLs at this point?

   delete sampler;
   delete slrts;
   delete ropl;
   delete profll;

   if (calcType != 2) {
      HypoTestPlot *plot = new HypoTestPlot(*htr, 100);
      plot->SetLogYaxis(true);
      plot->Draw();
   } else {
      std::cout << "Asymptotic results " << std::endl;
   }

   // look at expected significances
   // found median of S+B distribution
   if (calcType != 2) {

      SamplingDistribution *altDist = htr->GetAltDistribution();
      HypoTestResult htExp("Expected Result");
      htExp.Append(htr);
      // find quantiles in alt (S+B) distribution
      double p[5];
      double q[5];
      for (int i = 0; i < 5; ++i) {
         double sig = -2 + i;
         p[i] = ROOT::Math::normal_cdf(sig, 1);
      }
      std::vector<double> values = altDist->GetSamplingDistribution();
      TMath::Quantiles(values.size(), 5, &values[0], q, p, false);

      for (int i = 0; i < 5; ++i) {
         htExp.SetTestStatisticData(q[i]);
         double sig = -2 + i;
         std::cout << " Expected p -value and significance at " << sig << " sigma = " << htExp.NullPValue()
                   << " significance " << htExp.Significance() << " sigma " << std::endl;
      }
   } else {
      // case of asymptotic calculator
      for (int i = 0; i < 5; ++i) {
         double sig = -2 + i;
         // sigma is inverted here
         double pval = AsymptoticCalculator::GetExpectedPValues(htr->NullPValue(), htr->AlternatePValue(), -sig, false);
         std::cout << " Expected p -value and significance at " << sig << " sigma = " << pval << " significance "
                   << ROOT::Math::normal_quantile_c(pval, 1) << " sigma " << std::endl;
      }
   }

   // write result in a file in case of toys
   bool writeResult = (calcType != 2);

   if (enableDetOutput) {
      writeResult = true;
      Info("StandardHypoTestDemo", "Detailed output will be written in output result file");
   }

   if (htr != NULL && writeResult) {

      // write to a file the results
      const char *calcTypeName = (calcType == 0) ? "Freq" : (calcType == 1) ? "Hybr" : "Asym";
      TString resultFileName = TString::Format("%s_HypoTest_ts%d_", calcTypeName, testStatType);
      // strip the / from the filename

      TString name = infile;
      name.Replace(0, name.Last('/') + 1, "");
      resultFileName += name;

      TFile *fileOut = new TFile(resultFileName, "RECREATE");
      htr->Write();
      Info("StandardHypoTestDemo", "HypoTestResult has been written in the file %s", resultFileName.Data());

      fileOut->Close();
   }
}
//HERE///////////////////////////////////////////////
void higgs_atlas(){
    using namespace RooFit;
    using namespace RooStats;

    // Define the observable
    RooRealVar x("x","x",110,135,"GeV");
    x.setBins(100); // from the plot

    // ASCII import the unbinned dataset
    RooDataSet data = *RooDataSet::read("Higgs.txt", x);
    
    // model for inv mass distribution
    // background model is a polinomial of degree 2
    RooRealVar a1("a1","a1",-160, -200, -100);
    RooRealVar a2("a2","a2", 2.7, 2, 4);
    RooPolynomial bModel("bModel", "bModel", x, RooArgList{a1, a2});

    // signal model
    RooRealVar mass("mass", "mass", 110, 150);
    RooRealVar width("width", "width", 1.744681); // 4.1/2.35
    RooRealVar alpha("alpha", "alpha", 0.6);
    RooRealVar n("n", "n", 20);
    RooCBShape sModel("sModel", "sModel", x, mass, width, alpha, n);

    // spectrum components
    RooRealVar nsignal("nsignal", "nsignal", 70, 0, 500);
    RooRealVar nbackground("nbackground", "nbackground", 70, 0, 500);
    
    // Extended Composite Model using LH formalism
    RooAddPdf model("model","model",RooArgSet(sModel,bModel),RooArgSet(nsignal,nbackground));
    
    // Fitting with maximum LH
    model.fitTo(data, Extended(kTRUE)); 
    
    // Plot data and model
    TCanvas* c1 = new TCanvas ("c1", "c1");
    RooPlot* rp = x.frame();
    data.plotOn(rp);
    model.plotOn(rp);
    rp->Draw();

    // Import all into a workspace
    RooWorkspace w("w");
    w.import(data);
    w.import(model);

    // Define a Model Config for signal + bgk model
    ModelConfig mc("mc", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("nsignal"));
    mc.SetObservables(*w.var("x"));         
     //floating p are nuisance
     w.defineSet("nuisParams","nbackground,a1,a2"); 
    mc.SetNuisanceParameters(*w.set("nuisParams")); 
    mc.SetSnapshot(nsignal);

    // Set constantness
    w.var("width")->setConstant(true); 
    w.var("alpha")->setConstant(true);
    w.var("n")->setConstant(true);
    // required for significance vs mass plot
    w.var("mass")->setConstant(true);

    w.import(mc);    

    // Model config for Higgs mass
    ModelConfig mc_mass("mc_mass", &w);
    mc_mass.SetPdf(*w.pdf("model"));
    mc_mass.SetParametersOfInterest(*w.var("mass"));
    mc_mass.SetObservables(*w.var("x"));
     // Nuisance
     w.defineSet("nuisParams2","nbackground,a1,a2,nsignal,width");
    mc_mass.SetNuisanceParameters(*w.set("nuisParams2"));
    mc.SetSnapshot(mass);

    w.import(mc_mass);

    // Save the workspace to a file
    TFile* output = new TFile("output.root", "RECREATE");
    w.writeToFile("output.root");

////////////// PART 2 //////////////////

    ProfileLikelihoodCalculator plc(data, mc_mass);
    plc.SetConfidenceLevel(0.68);
    LikelihoodInterval* interval = plc.GetInterval();
    auto poi = static_cast<RooRealVar*>(mc_mass.GetParametersOfInterest()->first());
    double lowerLimit = interval->LowerLimit(*poi);
    double upperLimit = interval->UpperLimit(*poi);

    TCanvas* c2 = new TCanvas ("c2", "c2");
    LikelihoodIntervalPlot plot2(interval);
    plot2.Draw();

////////////// PART 3 //////////////////

    FeldmanCousins fc(data, mc_mass);
    fc.SetConfidenceLevel(0.9);
    fc.SetNBins(10); // set how many points per parameter of interest to scan
    fc.CreateConfBelt(true);

    PointSetInterval *interval_fc = fc.GetInterval();
    ConfidenceBelt *belt = fc.GetConfidenceBelt();

    //StandardHypoTestDemo("output.root", "w","mc", "mc2", "data", 2,3);

    RooRealVar *firstPOI = (RooRealVar *)mc_mass.GetParametersOfInterest()->first();

   cout << "\n>>>> RESULT PROFILE LIKELIHOOD : " << 0.68 * 100 << "% interval on " << firstPOI->GetName() << " is : ["
         << interval->LowerLimit(*firstPOI) << ", " << interval->UpperLimit(*firstPOI) << "]\n " << endl;

   cout << "\n>>>> RESULT NEYMANN-COUSINS : "<< "90% interval on " << firstPOI->GetName() << " is : [" << interval_fc->LowerLimit(*firstPOI) << ", "
        << interval_fc->UpperLimit(*firstPOI) << "] " << endl;

}

////////////// PART 4 //////////////////

void test(){
    StandardHypoTestDemo("output.root", "w","mc", "mc2", "data", 2,3);
}


