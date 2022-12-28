// Exercise 2022_1 tutorial1 - Alessandro Mancini
// alessandro.mancini16@studio.unibo.it

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"

using namespace RooFit;

TCanvas c1{"c1", "Plots"};

int ex_2022_1()
{
  // Impostazione del modello
  // ---------------------

  // Dichiara le variabili x,media,sigma associando nome, titolo, valore
  // iniziale e intervallo x e' l'osservabile, media e sigma sono parametri
    RooRealVar x("x", "x", -5, 5);
    RooRealVar mean("mean", "mean", 0, -3, 3);
    RooRealVar sigma("sigma", "sigma", 1, 0, 5);
    RooRealVar alpha("alpha", "alpha", 1.5, 0, 10);
    RooRealVar n("n", "n", 1.5, 0, 10);


  // Costruisci un modello utilizzando una p.d.f. crystal ball in termini di x,mean
  // and sigma
    RooCBShape model("model", "Crystal Ball", x, mean, sigma, alpha, n);

  // Costruisci la vista (plot frame) per 'x'
    RooPlot* frame = x.frame();

  // Visualizza il modello. Cambia un parametro
  // ------------------------------------------

  // Disegna il modello nella vista (i.e. in x)
    model.plotOn(frame);

  // Cambiare il valore del parametro sigma a 0.3  ( suggerimento: usare il metodo
  // setVal di RooRealVar)
    sigma.setVal(0.3);

  // Disegna il modello modificato nella stessa vista (i.e. in x)
    model.plotOn(frame);
    
  // ToyMC: Generare 10000 eventi
  // ----------------------------

  // Genera un dataset (unbinned) di 1000 eventi x secondo il modello
    RooDataSet* data = model.generate(x, 10000);

  // Costruisci una seconda vista (plot frame) per 'x'
  // e disegna sia il modello che il dataset nel in questa nuova vista
    RooPlot* frame2 = x.frame();
    data->plotOn(frame2);
    model.plotOn(frame2);
    
  // Fittare il modello ai dati
  // --------------------------

  // Fittare modello ai dati
    model.fitTo(*data); //We can FIT also AFTER THE PLOTTING!

  // Visualizzare (sulla shell) i valori dei parametri media e sigma
  // (ora dovrebbero rappresentare il risultato del fit)
    mean.Print();
    sigma.Print();
    

  // Visualizza le due viste su una canvas
    c1.Divide(1,2);
    c1.cd(1);
    frame->Draw();
    c1.cd(2);
    frame2->Draw();
    
    c1.Print("ex_2022_1.png");
    
    return 0;
}
