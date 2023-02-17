// Exercise 2022 RooUnfold - Tutorial 4
// alessandro.mancini16@studio.unibo.it

const double cutdummy = -99999.0;

double smear(double xt) {
    double xeff = 0.3+(1.0-0.3)/20*(xt+10.0); // Efficiency
    double x = gRandom->Rndm();
    
    if(x > xeff) return cutdummy;
    
    // Smear (0.2) and bias (-2.5)
    double xsmear = gRandom->Gaus(-2.5, 0.2);
    return xt+xsmear;
}

void ex_2022_roounfold() {
    
    // Include the RooUnfold library
    gSystem->Load("$HOME/RooUnfold/build/libRooUnfold.so");
    
    RooUnfoldResponse response(40, -10.0, 10.0);
    
    // Train with a Breit-Wigner, mean 0.3 and width 2.5
    for(int i = 0; i != 100000; ++i) {
        double xt = gRandom->BreitWigner(0.3, 2.5);
        double x = smear(xt);
        if(x != cutdummy) {
            response.Fill(x, xt);
        }
        else {
            response.Miss(xt);
        }
    }
    
    TH1D* hTrue = new TH1D("hTrue", "Test Truth", 40, -10.0, 10.);
    TH1D* hMeas = new TH1D("hMeas", "Test Measured", 40, -10.0, 10.);
    
    // Test with Gaussian
    for(int i = 0;  i != 10000; ++i) {
        double xt = gRandom->Gaus(0.0, 2.0);
        double x = smear(xt);
        hTrue->Fill(xt);
        if(x != cutdummy) hMeas->Fill(x);
    }
    
    // Choose (uncomment) one of the following unfolding methods
    RooUnfoldBayes unfold(&response, hMeas, 4);
    // RooUnfoldSvd unfold(&response, hMeas, 20);
    // RooUnfoldIds unfold(&response, hMeas, 1);
    
    TH1D* hReco = (TH1D*)unfold.Hreco();
    
    TLegend* leg = new TLegend();
    leg->AddEntry(hTrue, "True histo");
    leg->AddEntry(hMeas, "Smeared (Measured) histo");
    leg->AddEntry(hReco, "Unfolded histo");
    
    TCanvas* c = new TCanvas("c", "RooUnfold Example");
    unfold.PrintTable(std::cout, hTrue);
    hReco->SetLineColor(kRed);
    // hReco->Reset();
    hReco->Draw();
    hMeas->Draw("SAME");
    hTrue->SetLineColor(kGreen);
    hTrue->Draw("SAME");
    leg->Draw("SAME");
    
    // Response matrix: (row, column) = (measured, truth)
    auto responseMatrix = response.Mresponse();
    
    TCanvas* c2 = new TCanvas("c2", "Response matrix");
    c2->Divide(2,1);
    c2->cd(1);
    responseMatrix.Draw("colz");
    
    // Inverse of responseMatrix
    auto responseInv = *static_cast<TMatrixD*>(responseMatrix.Clone());
    responseInv.Invert();
    
    // Plot inverse matrix
    // TCanvas* c3 = new TCanvas("c3", "Inverse response matrix");
    c2->cd(2);
    responseInv.Draw("colz");
    
    // Apply inverse of response matrix to hMeas
    auto hMeas_Rinv = static_cast<TH1D*>(hMeas->Clone("hMeas_Rinv"));
    hMeas_Rinv->Reset();
    for(int ibin = 1; ibin != hMeas_Rinv->GetNbinsX()+1; ++ibin) {
        double matsum = 0;
        for(int jbin = 1; jbin != hMeas_Rinv->GetNbinsX()+1; ++jbin) {
            // The following is the operation of multiplying every bin content (vector of bins) times the elements of the matrix responseInv (operation Sum(Rij*vj))
            matsum += hMeas->GetBinContent(jbin)*responseInv[ibin-1][jbin-1];
        }
        hMeas_Rinv->SetBinContent(ibin, matsum);
    }
    
    hMeas_Rinv->SetLineColor(kOrange);
    // leg->AddEntry(hMeas_Rinv, "R inv");
    
    TCanvas* c4 = new TCanvas("c4", "Matrix inversion approach");
    // c4->cd();
    // hTrue->Draw();
    // hMeas->Draw("SAME");
    // hReco->Draw("SAME");
    // hMeas_Rinv->Draw("SAME");
    c4->Divide(2, 2);
    c4->cd(1);  hTrue->Draw("hist");
    c4->cd(2);  hMeas->Draw("pe same");
    c4->cd(3);  hMeas_Rinv->Draw("histsame");
    c4->cd(4);  hReco->Draw("hist");
    // leg->Draw("SAME");
    
}
