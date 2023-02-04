// Exercise 2022.2.1 part 2 ARIADNE - Tutorial 2
// Alessandro Mancini - alessandro.mancini16@studio.unibo.it

TCanvas c("c", "Plots");

void ex2022_2_1b() {
    RooWorkspace w("w");
    
    // Define the variable of the model
    RooRealVar x("x", "dI/dX", 9400, 29656, "ADU/cm");
    w.import(x);
    
    x.setBins(49);
    
    // Define the pdfs for the convolution model
    w.factory("Landau::lan(x, p0[9000, 1, 16000], p1[500, 1, 10000])");
    w.factory("Gaussian::gaus(x, m[4000, 1, 16000], s[1000, 1, 10000])");
    
    // Define the convolution model
    w.factory("FCONV::model(x, lan, gaus)");
    
    // Define the dataset of input data
    RooDataHist data("data", "data", x);
    
    // Read input data from file
    ifstream file("ariadne_g012.dat");
    double val, weight;
    while(!file.eof()) {
        file >> val >> weight;
        x.setVal(val);
        data.set(x, weight);
    }
    
    file.close();
    
    // Fit and save results
    RooFitResult* res = w.pdf("model")->fitTo(data, RooFit::Save());
    
    // Print fit results on a .txt file
    std::ofstream outfile("ex2022_2_1b_FitResults.txt");
    res->printTitle(outfile);
    res->printMultiline(outfile, 0, true);
    outfile.close();
    
    // Draw data and pdf on frame
    RooPlot* frame = x.frame(RooFit::Title("Landau-Gaussian convolution fit"));
    data.plotOn(frame);
    
    // Uncomment following to see also the pdfs of convolution
    // w.pdf("lan")->plotOn(frame, RooFit::LineColor(kGreen));
    // w.pdf("gaus")->plotOn(frame, RooFit::LineColor(kRed));
    w.pdf("model")->plotOn(frame, RooFit::LineColor(kBlue));
    
    w.pdf("model")->paramOn(frame, RooFit::Label("Fit parameters"));
    
    frame->Draw();
    
    // Save image as .png
    c.Print("ex2022_2_1b.png");
}
