// Exercise 2022.2.2 rarest B0 decay - Tutorial 2
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it


void makeModel(RooWorkspace& w) {
    
    // Define Background of model
    w.factory("Exponential::bkg(x, r[-0.1, -10, 0])");
    
    // Define gaussian peak of B0 mass
    w.factory("Gaussian::g_B0(x, m_B0[5279, 5260, 5290], s_B0[50, 0, 100])");
    
    // Define gaussian peak of B0s mass
    w.factory("Gaussian::g_B0s(x, m_B0s[5366, 5350, 5380], s_B0s[10, 0, 100])");
    
    // Create composite model: bkg + B0 + B0s
    w.factory("SUM::model(N_bkg[50, 0, 100]*bkg, N_B0[20, 0, 100]*g_B0, N_B0s[10, 0, 50]*g_B0s)");
    
}

void useModel(RooWorkspace& w) {
    
    // Extract variable from workspace
    auto x = w.var("x");
    RooDataSet data = *RooDataSet::read("rarest_b0_decay.dat", *x, "v");
    
    // Extract model from workspace
    auto model = w.pdf("model");
    
    // Fit model to data
    model->fitTo(data);
    
    // Plot model and data
    RooPlot* frame = x->frame(RooFit::Title("B0 Invariant Mass"));
    data.plotOn(frame);
    model->plotOn(frame, RooFit::Components("g_B0"), RooFit::LineColor(kRed));
    model->plotOn(frame, RooFit::Components("g_B0s"), RooFit::LineColor(kGreen));
    model->plotOn(frame, RooFit::Components("bkg"), RooFit::LineColor(kBlack));
    model->plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1", "Fit");
    c1->cd();
    frame->Draw();
    // Save fit image as PNG
    c1->Print("ex2022_2_2_Fit.png");
    
    // Create residual and pull distributions
    // The following computes residuals of LAST plotted hist and LAST plotted curve on frame if we DON'T specify arguments
    RooHist* h_res = frame->residHist();
    RooHist* h_pull = frame->pullHist();
    
    // Plot residuals and pull distributions
    RooPlot* frame2 = x->frame(RooFit::Title("Residuals distribution"));
    // Add the res histo to the frame
    frame2->addPlotable(h_res, "P");
    
    TCanvas* c2 = new TCanvas("c2", "Residuals");
    c2->cd();
    frame2->Draw();
    // Save Residuals distrib image as PNG
    c2->Print("ex2022_2_2_Res.png");
    
    RooPlot* frame3 = x->frame(RooFit::Title("Pull distribution"));
    frame3->addPlotable(h_pull, "P");
    
    TCanvas* c3 = new TCanvas("c3", "Pull");
    c3->cd();
    frame3->Draw();
    // Save Pull distrib image as PNG
    c3->Print("ex2022_2_2_Pull.png");
    
}

void ex2022_2_2() {
    
    RooWorkspace w("w");
    
    // Define model variable
    RooRealVar x("x", "m(p#bar{p})", 5080, 5555, "MeV/c^{2}");
    x.setBins(25);
    w.import(x);
    
    
    makeModel(w);
    
    useModel(w);
    
}
