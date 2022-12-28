// Exercise 2022_2 tutorial1 - Alessandro Mancini
// alessandro.mancini16@studio.unibo.it


TCanvas c{"c", "Plots"};

int ex_2022_2()
{
  // Impostazione del modello
  // ---------------------

  // Dichiara le variabili x,tau associando nome, titolo, valore
  // iniziale e intervallo. x e' l'osservabile, tau un parametri
    RooRealVar x("x", "x", 0, 10);
    RooRealVar tau("tau", "exp param", 3, 0, 5);
    RooFormulaVar rate("rate", "rate", "-1/tau", RooArgSet(tau));
    
  // Costruisci un modello utilizzando una p.d.f. esponenziale in termini di x, rate
    RooExponential exp("exp", "exponential pdf", x, rate);
    
  // Costruisci la vista (plot frame) per 'x'
    RooPlot* frame = x.frame(RooFit::Title("Exponential pdf"));
  
  // Visualizza il modello. Cambia un parametro
  // ------------------------------------------
    
  // Disegna il modello nella vista (i.e. in x)
    exp.plotOn(frame);
    
  // ToyMC: Generare 1000 eventi
  // ----------------------------
    
  // Genera un dataset (binned) di 1000 eventi x secondo il modello
    x.setBins(20); //bin width = 0.5
    RooDataHist* data = exp.generateBinned(x, 1000);
    
  // Costruisci una seconda vista (plot frame) per 'x'
  // e disegna sia il modello che il dataset nel in questa nuova vista
    RooPlot* frame2 = x.frame(RooFit::Title("PDF fitted to data"));
    
    data->plotOn(frame2);
    exp.plotOn(frame2);

  // Fittare il modello ai dati
  // --------------------------

  // Fittare modello ai dati
    exp.fitTo(*data);
    
  // Visualizzare (sulla shell) i valori dei parametri tau e rate
  // (ora dovrebbero rappresentare il risultato del fit)
    tau.Print();
    rate.Print();
    
  // Visualizza le due viste su una canvas
    c.Divide(2, 1);
    
    c.cd(1);
    frame->Draw();
    
    c.cd(2);
    frame2->Draw();
    
    c.Print("ex_2022_2.png");

    return 0;
}
