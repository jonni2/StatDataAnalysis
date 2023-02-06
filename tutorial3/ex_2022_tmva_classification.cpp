// Exercise 2022.tmva - Tutorial 3
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it

void ex_2022_tmva() {
    
    // Create output Root File
    TFile* outputFile = TFile::Open("ex_2022_tmva.root", "RECREATE");
    
    // Declare TMVA Factory
    TMVA::Factory factory{"jonni2_TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification"};
    
    // Declare DataLoader
    TMVA::DataLoader loader("dataset");
    
    // loader.AddVariable("var1", 'F');
    // loader.AddVariable("var2", 'F');
    loader.AddVariable("DER_pt_ratio_lep_tau", 'F');
    loader.AddVariable("PRI_lep_phi");
    loader.AddVariable("DER_deltar_tau_lep");
    
    TFile* input_sigFile = TFile::Open("atlas-higgs-challenge-2014-v2-sig.root");
    TFile* input_bkgFile = TFile::Open("atlas-higgs-challenge-2014-v2-bkg.root");
    
    // Declare TTrees for training
    TTree* tsig;
    TTree* tbkg;
    
    // Put the input datasets in the trees
    input_sigFile->GetObject("tree_sig", tsig);
    tsig->Print();
    input_bkgFile->GetObject("tree_bkg", tbkg);
    
    // tsig->ls();
    
    loader.AddSignalTree(tsig, 1.0); // signal weight = 1
    loader.AddBackgroundTree(tbkg, 1.0); // background weight = 1
    
    TCut mycuts, mycutb;
    loader.PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=10000:nTrain_Background=20000:SplitMode=Random:NormMode=NumEvents:!V");
    
    //factory.BookMethod(&loader, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator");
    
    factory.BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    
    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();
    
    auto c = factory.GetROCCurve(&loader);
    
    c->Draw();
}
