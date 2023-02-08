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
    RooRealVar mixing("mixing", "(sin(2theta))^2", 0.95, 0., 1.);
    // mass difference
    RooRealVar dm2("dm2", "mass difference", 0.0023, 0., 10.);
    // L = 730 km
    RooRealVar L("L", "distance", 730, "km");
    RooFormulaVar prob_osc("prob_osc", "1-mixing*(sin(1.267*dm2*L/x))^2", RooArgSet(mixing, dm2, L, x));
    RooGenericPdf model("model", "model", "@0*@1", RooArgSet(prob_osc, func_noosc));
    
    
    RooPlot* frame = x.frame(RooFit::Title("MINOS #nu_{#mu}->#nu_{#mu} oscillation"));
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
    
    
    
}
