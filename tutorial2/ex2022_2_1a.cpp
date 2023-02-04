// Exercise 2022.2.1 part 1 ARIADNE - Tutorial 2
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it

TCanvas c("c", "Plots");

void ex2022_2_1a() {
    RooRealVar x("x", "Drift time", 460, 530, "#mus");
    x.setBins(70);
    
    RooUniform bkg("bkg", "uniform component", x);
    
    RooRealVar mean("mean", "mean", 450, 400, 500);
    RooRealVar gamma("gamma", "gamma", 20, 0, 100);
    RooLandau lan("lan", "landau component", x, mean, gamma);
    
    RooRealVar fs("fs", "signal fraction", 0.5, -0.5, 1.5);
    
    RooAddPdf model("model", "model", RooArgList(lan, bkg), fs);
    
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
    RooFitResult* res = model.fitTo(data, RooFit::Save());
    
    // Print fit results on file FitResults.txt
    std::ofstream outfile("ex2022_2_1a_FitResults.txt");
    res->printTitle(outfile);
    res->printMultiline(outfile, 0, true);
    outfile.close();
    
    RooPlot* xframe = x.frame(RooFit::Title("Drift times with Landau fit"));
    data.plotOn(xframe);
    model.plotOn(xframe);
    
    // Signal (Landau) component in red
    model.plotOn(xframe, RooFit::Components("lan"), RooFit::LineColor(kRed));
    // Background (Uniform) component in black
    model.plotOn(xframe, RooFit::Components("bkg"), RooFit::LineColor(kBlack));
    
    
    // Plot fit params in box. We can use different arguments:
    // RooFit::Parameters(RooArgSet(desired params));
    // RooFit::Label("label");
    // RooFit::Layout(xmin, xmax, ymax);
    model.paramOn(xframe, RooFit::Label("Fit parameters"));
    
    xframe->Draw();
    c.Print("ex2022_2_1a.png");
}
