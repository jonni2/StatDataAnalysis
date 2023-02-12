// Exercise 2022.tmva - Tutorial 3
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it

void ex_2022_tmva_classification() {
    
    // Create output Root File
    TString S_outfile("ex_2022_tmva.root");
    TFile* outputFile = TFile::Open(S_outfile, "RECREATE");
    
    // Declare TMVA Factory
    TMVA::Factory factory{"ex_2022_TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification"};
    
    // Declare DataLoader
    TMVA::DataLoader loader("dataset");
    
    // Add variables of the dataset
    // Best variables to use (variables that have more separated distributions in sig, bkg):
    // DER_mass_MMC
    // DER_mass_transverse_met_lep
    // DER_mass_vis
    // DER_pt_ratio_lep_tau
    // PRI_tau_pt
    // PRI_met
    loader.AddVariable("DER_mass_MMC", 'F');
    loader.AddVariable("DER_mass_transverse_met_lep", 'F');
    loader.AddVariable("DER_mass_vis", 'F');
    loader.AddVariable("DER_pt_ratio_lep_tau", 'F');
    loader.AddVariable("PRI_tau_pt", 'F');
    loader.AddVariable("PRI_met", 'F');
    
    TFile* input_sigFile = TFile::Open("atlas-higgs-challenge-2014-v2-sig.root");
    TFile* input_bkgFile = TFile::Open("atlas-higgs-challenge-2014-v2-bkg.root");
    
    // Declare TTrees for training
    TTree* tsig;
    TTree* tbkg;
    
    // Put the input datasets in the trees
    input_sigFile->GetObject("tree_sig", tsig);
    tsig->Print();
    input_bkgFile->GetObject("tree_bkg", tbkg);
    
    // List Tree content
    // tsig->ls();
    
    loader.AddSignalTree(tsig, 1.0); // signal weight = 1
    loader.AddBackgroundTree(tbkg, 1.0); // background weight = 1
    
    TCut mycuts, mycutb;
    loader.PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=10000:nTrain_Background=20000:SplitMode=Random:NormMode=NumEvents:!V");
    
    // -----------------  BOOKING METHODS  ----------------------
    
    // Cuts optimisation
    factory.BookMethod(&loader, TMVA::Types::kCuts, "ex2022_Cuts", "!H:!V:FitMethod=MC");
    
    // Multi Layer Perceptron
    factory.BookMethod(&loader, TMVA::Types::kMLP, "ex2022_MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:TestRate=5:HiddenLayers=N+5:!UseRegulator");
    
    // Boost Decision Tree
    factory.BookMethod(&loader, TMVA::Types::kBDT, "ex2022_BDT", "!V:NTrees=200:BoostType=AdaBoost");
    
    // Fisher
    factory.BookMethod(&loader, TMVA::Types::kFisher, "ex2022_Fisher", "H:!V:Fisher");
    
    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();
    
    TMVA::TMVAGui("ex_2022_tmva.root");
    
}

void ex_2022_tmva_application() {
    
    TString file("application_event.csv");
    TTree* t = new TTree("t", "application event");
    t->ReadFile("application_event.csv", "EventId/F:DER_mass_MMC:DER_mass_transverse_met_lep:DER_mass_vis:DER_pt_h:DER_deltaeta_jet_jet:DER_mass_jet_jet:DER_prodeta_jet_jet:DER_deltar_tau_lep:DER_pt_tot:DER_sum_pt:DER_pt_ratio_lep_tau:DER_met_phi_centrality:DER_lep_eta_centrality:PRI_tau_pt:PRI_tau_eta:PRI_tau_phi:PRI_lep_pt:PRI_lep_eta:PRI_lep_phi:PRI_met:PRI_met_phi:PRI_met_sumet:PRI_jet_num:PRI_jet_leading_pt:PRI_jet_leading_eta:PRI_jet_leading_phi:PRI_jet_subleading_pt:PRI_jet_subleading_eta:PRI_jet_subleading_phi:PRI_jet_all_pt:Weight:Label:KaggleSet:KaggleWeight");
    
    // Declare Reader
    TMVA::Reader reader("!Color:!Silent");
    
    // Add variables to reader
    float x1, x2, x3, x4, x5, x6;
    reader.AddVariable("DER_mass_MMC", &x1);
    reader.AddVariable("DER_mass_transverse_met_lep", &x2);
    reader.AddVariable("DER_mass_vis", &x3);
    reader.AddVariable("DER_pt_ratio_lep_tau", &x4);
    reader.AddVariable("PRI_tau_pt", &x5);
    reader.AddVariable("PRI_met", &x6);
    
    // Select weights file path - chosen MLP
    reader.BookMVA("ex2022_MLP", "dataset/weights/ex_2022_TMVAClassification_ex2022_MLP.weights.xml");
    
    // Select variables from TTree
    t->SetBranchAddress("DER_mass_MMC", &x1);
    t->SetBranchAddress("DER_mass_transverse_met_lep", &x2);
    t->SetBranchAddress("DER_mass_vis", &x3);
    t->SetBranchAddress("DER_pt_ratio_lep_tau", &x4);
    t->SetBranchAddress("PRI_tau_pt", &x5);
    t->SetBranchAddress("PRI_met", &x6);
    
    // Classifier response
    double tMLP = reader.EvaluateMVA("ex2022_MLP");
    double mvaErr = reader.GetMVAError();
    
    std::cout << "MLP response: " << tMLP << ", error: " << mvaErr << '\n';
    
}
