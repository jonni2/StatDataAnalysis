// Exercise 2022.2.3 (MINOS) composite model - Tutorial 2
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it

void ex2022_2_3() {

    // PART 1: working with the model
    
    // Variable of the model
    RooRealVar x("x", "Neutrino Energy", 0.5, 14, "GeV");
    x.setBins(27);
    
    // Load unbinned dataset of observed events by MINOS
    RooDataSet data = *RooDataSet::read("minos_2013_data.dat", x, "v");
    
    // RooRealVar noosc_E("noosc_E", "non-oscillated neutrino energy", 0.5, 14, "GeV");
    // Load SIMULATED (MC) unbinned dataset
    RooDataSet mc_noosc = *RooDataSet::read("minos_2013_mc.dat", x, "v");
    
    RooDataSet* dd = (RooDataSet*)mc_noosc.reduce(RooArgSet(x));
    RooDataHist* dh_mc_noosc = dd->binnedClone();
    
    RooHistFunc func_noosc("func_mc_noosc", "No oscillation", x, *dh_mc_noosc, 2);
    
    // Parameters of the model
    // mixing := (sin(2theta))^2
    RooRealVar mixing("mixing", "sin^{2}(2#theta)", 0.95, 0., 1.);
    // mass difference
    RooRealVar dm2("dm2", "|#Deltam^{2}|", 0.0023, 0., 0.0035, "eV^{2}");
    // L = 730 km
    RooRealVar L("L", "distance", 730, "km");
    RooFormulaVar prob_osc("prob_osc", "1-mixing*(sin(1.267*dm2*L/x))^2", RooArgSet(mixing, dm2, L, x));
    RooGenericPdf model("model", "model", "@0*@1", RooArgSet(prob_osc, func_noosc));
    
    
    RooPlot* frame = x.frame(RooFit::Title("MINOS, Energy distribution for #nu_{#mu}->#nu_{#mu} oscillated neutrinos"));
    data.plotOn(frame);
    
    model.fitTo(data);
    model.plotOn(frame);
    
    std::cout << "Parameters from Fit:\n";
    mixing.Print();
    dm2.Print();
    
    TCanvas* c = new TCanvas("c", "Neutrino oscillation");
    c->cd();
    frame->Draw();
    
    c->Print("minos_data.png");
    
    // PART 2: working with LIKELIHOOD
    
    std::cout << "\n \u001b[31m ---------- PART 2: working with likelihood ---------- \u001b[0m \n\n";
    
    // Create negative likelihood object
    RooAbsReal* nll = model.createNLL(data);
    
    // Create RooMinuit interface
    RooMinuit m(*nll);
    
    // SetVerbose option: show each step in Minuit algorithm
    m.setVerbose(kTRUE);
    
    // Call Migrad algorithm to find nll minimum
    m.migrad();
    
    // Printing fit parameters and errors
    dm2.Print();
    mixing.Print();
    
    // Disable Verbose option
    m.setVerbose(kFALSE);
    
    // Run Hesse
    m.hesse();
    
    // Print parameters and errors
    dm2.Print();
    mixing.Print();
    
    // Run Minos only on dm2
    m.minos(RooArgSet(dm2));
    
    // Print value of dm2
    dm2.Print();
    
    // Save results
    RooFitResult* result = m.save();
    
    // Print results
    result->Print("v");
    
    // Make contour plot
    RooPlot* frame_contout = m.contour(mixing, dm2, 1, 2, 3);
    frame_contout->SetTitle("Contour plot, |#Deltam^{2}| vs sin^{2}(2#theta)");
    
    TCanvas* c_contour = new TCanvas("c_contour", "Contour plot");
    c_contour->cd();
    frame_contout->Draw();
    c_contour->Print("minos_likelihood.png");
    
    
}
