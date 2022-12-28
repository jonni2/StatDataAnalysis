// Exercise 16.1 tutorial1 - Alessandro Mancini
// alessandro.mancini16@studio.unibo.it

TCanvas c("c", "Plots");

void ex16() {
    //RooWorkspace w("w");
    
    // Definition of 8 indep. random vars
    RooRealVar x1("x1", "var1", 0, 1);
    RooRealVar x2("x2", "var1", 0, 1);
    RooRealVar x3("x3", "var1", 0, 1);
    RooRealVar x4("x4", "var1", 0, 1);
    RooRealVar x5("x5", "var1", 0, 1);
    RooRealVar x6("x6", "var1", 0, 1);
    RooRealVar x7("x7", "var1", 0, 1);
    RooRealVar x8("x8", "var1", 0, 1);
    
    // Definition of 8 uniform pdfs for the 8 random vars
    RooUniform u1("u1", "Uniform 1", x1);
    RooUniform u2("u2", "Uniform 1", x2);
    RooUniform u3("u3", "Uniform 1", x3);
    RooUniform u4("u4", "Uniform 1", x4);
    RooUniform u5("u5", "Uniform 1", x5);
    RooUniform u6("u6", "Uniform 1", x6);
    RooUniform u7("u7", "Uniform 1", x7);
    RooUniform u8("u8", "Uniform 1", x8);
    
    // Vector of all the variables
    RooArgSet x_set(x1,x2,x3,x4,x5,x6,x7,x8);
    
    // Vector of all the distributions
    RooArgSet model_set(u1,u2,u3,u4,u5,u6,u7,u8);
    
    // model_set.add(pdf9); //In case you want to do it with MORE THAN 9 vars (RooArgList has a LIMIT to 9 variables)
    
    // The distrib of the set of 8 vars is the PRODUCT of the 8 pdfs, since we assume that they are INDEPENDENT
    RooProdPdf model("model", "model", model_set);
    
    // Generate data for the set of 8 indep. vars
    RooDataSet* data = model.generate(x_set, 10000);
    
    // SUM variable of all the 8 variables
    RooFormulaVar fsum("xsum", "x1+x2+x3+x4+x5+x6+x7+x8-4", x_set);
    
    auto xsum = static_cast<RooRealVar*>(data->addColumn(fsum));
    
    // SUM of 2 of the 8 variables.
    RooFormulaVar fsum2("xsum2", "x1+x2-1", RooArgList(x1,x2));
    auto xsum2 = static_cast<RooRealVar*>(data->addColumn(fsum2));
    
    // Frame for the sum of all the 8 vars
    // We need to set bins, title and range because the xsum defined by data->addColumn is a RooFormulaVar, not a RooRealVar
    auto xframe = xsum->frame(RooFit::Bins(40), RooFit::Title("Sum of 8 variables"), RooFit::Range(-4, 4));
    
    // Frame for the sum of 2 of the 8 vars
    auto xframe2 = xsum2->frame(RooFit::Bins(40), RooFit::Title("Sum of 2 variables"), RooFit::Range(-2, 2));
    
    // Parameters for the gaussian fit for xsum
    RooRealVar mean("mean", "mean", 0, -1, 1);
    RooRealVar sigma("sigma", "sigma", 0.5, 0, 1);
    RooGaussian gaus("gaus", "gaussian", *xsum, mean, sigma);
    
    // Fit
    gaus.fitTo(*data);
    data->plotOn(xframe);
    gaus.plotOn(xframe);
    
    // RooFit automatically UNDERSTANDS what it has to plot according to the definition of the xframe
    data->plotOn(xframe2);
    
    // Drawing plots
    c.Divide(2,1);
    
    c.cd(1);
    
    xframe->Draw();
    
    c.cd(2);
    xframe2->Draw();
    
    c.Print("ex16.png");
}
