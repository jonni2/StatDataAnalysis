// File to read and visualize TTree. Useful to read data for TMVA

void readTree() {
    TString path("inputdata.root");
    
    TFile* file = new TFile(path, "READ");
    
    file->ls();
    TTree* t_sig = (TTree*)file->Get("Sig;1");
    t_sig->Print();
    
    TTree* t_bkg = (TTree*)file->Get("Bkg;1");
    t_bkg->Print();
    
    t_bkg->Draw("weight");
}
