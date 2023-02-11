// Test TMVA
// This macro is the BASIC Training TMVA macro which has ALL the NECESSARY commands and declarations to make TMVA work.
// To EXECUTE this macro just type from Shell: root TmvaTest.cpp
// ATTENTION: this macro will generate few FILES while executing, so try to execute it inside a SPECIFIC FOLDER

void TmvaTest() {
    
    // 1) Declare FACTORY
    auto outputFile = TFile::Open("TmvaTest.root", "RECREATE");
    
    TMVA::Factory factory("TMVAClassification", outputFile, "!V:ROC:!Correlations:!Silent:Color:!DrawProgressBar:AnalysisType=Classification");
    
    // 2) Declare DataLoader
    TMVA::DataLoader loader("dataset");
    
    loader.AddVariable("var1", 'F');
    loader.AddVariable("var2", 'F');
    //loader.AddVariable("var3 := a*b", 'F'); // We can use EXPRESSIONS
    
    // Set individual weights
    // loader->SetSignalWeightExpression("weight1*weight2");
    
    // 3) Setup DATASET for TRAINING
    auto inputFile = TFile::Open(
    "https://raw.githubusercontent.com/iml-wg/tmvatutorials/master/inputdata.root");
    
    TTree* tsignal;
    TTree* tbackground;
    
    // "Sig" and "Bkg" are the NAMES of the TTrees in the TFile (you can see them if you do file->ls())
    inputFile->GetObject("Sig", tsignal);
    inputFile->GetObject("Bkg", tbackground);
    
    TCut mycuts, mycutb;
    
    loader.AddSignalTree(tsignal, 1.0); // signal weight = 1
    loader.AddBackgroundTree(tbackground, 1.0); // background weight = 1
    
    // tsignal->ls();
    
    loader.PrepareTrainingAndTestTree(mycuts, mycutb,
    "nTrain_Signal=1000:nTrain_Background=1000:"
    "SplitMode=Random:NormMode=NumEvents:!V");
    
    // DECIDE which ALGORITHM to use for TRAINING PHASE:
    // HERE we're choosing BoostDecisionTree (BDT)
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    
    // Start TRAINING PHASE
    factory.TrainAllMethods();
    
    factory.TestAllMethods();
    factory.EvaluateAllMethods();
    
    // Plot ROC curve
    auto c = factory.GetROCCurve(&loader);
    c->Draw();
    
    // Open the Gui to check training (the name must be tha same of the outputFile)
    TMVA::TMVAGui("TmvaTest.root");
    
}
