// Exercise 12.2 tutorial1 - Alessandro Mancini
// alessandro.mancini16@studio.unibo.it

TCanvas c("c", "Plots");

// Function that returns a (pointer to) HISTOGRAM read from ROOT File
auto ReadData(char* path = "B0sInvariantMass.root") {
    TFile* file = new TFile(path, "READ");
    file->ls(); // Just to visualize the root file content
    TH1F* h = (TH1F*)file->Get("massaB0;1");
    return h;
}

// Function that creates Breit Wigner with Factory
void makeBreitWigner(RooWorkspace& w) {
    // With the following line we can AUTOMATICALLY IMPORT the RooRealVar x
    w.factory("BreitWigner::bw(x, mean_bw[5, 0, 6], gamma[0.1, 0.00001, 5])");
    // There's NO NEED to call w.import(bw);
}

// Function that creates Gaussian with Factory
void makeGaussian(RooWorkspace& w) {
    // We need ANOTHER name for the mean ("mean" is already used in the bw)
    w.factory("Gaussian::gauss(x, mean_gauss[5, 0, 6], sigma[0.1, 0.00001, 5])");
}

// Function that uses bw and gauss and FITS them to the data
void useModel(RooWorkspace& w) {
    
    // Extracting RooRealVar x from RooWorkspace
    auto x = w.var("x");
    RooPlot* frame_bw = x->frame(RooFit::Title("Breit Wigner fit"));
    RooPlot* frame_gauss = x->frame(RooFit::Title("Gaussian fit"));
    
    // Extracting data and models from RooWorkspace
    auto data = w.data("data");
    auto bw = w.pdf("bw");
    auto gauss = w.pdf("gauss");
    
    
    // Fitting models
    bw->fitTo(*data);
    gauss->fitTo(*data);
    
    // Drawing models
    data->plotOn(frame_bw);
    data->plotOn(frame_gauss);
    bw->plotOn(frame_bw);
    gauss->plotOn(frame_gauss, RooFit::LineColor(kRed));
    
    // Printing fit parameters
    w.var("mean_bw")->Print();
    w.var("mean_gauss")->Print();
        
    c.Divide(2,1);
    c.cd(1);
    frame_gauss->Draw();
    
    c.cd(2);
    frame_bw->Draw();
    
    // Save graph in png format
    c.Print("ex12_2.png");
    
    // Export RooWorkspace to a Root File for future use
    w.writeToFile("B0s_Fit.root");
}


void ex12_2() {
    
    TH1F* h_B0Mass = ReadData();
    
    RooRealVar x("x", "Inv. mass B0s", 0, 6, "GeV/c^{2}"); //B0 Inv-Mass (GeV/c^2)
    
    RooDataHist data("data", "B0 InvMass", x, RooFit::Import(*h_B0Mass));
    
    RooWorkspace w("w");
    w.import(x);
    w.import(data);
    
    makeBreitWigner(w);
    makeGaussian(w);
    
    useModel(w);
    
}
