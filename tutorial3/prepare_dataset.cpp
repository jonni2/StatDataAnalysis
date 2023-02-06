// Macro written by professor to EXTRACT Root files from CSV files
//  allows to analyse data stored in TTrees (CSVs) with a high level interface.

void prepare_dataset()
{
  auto df = ROOT::RDF::MakeCsvDataFrame("atlas-higgs-challenge-2014-v2.csv");

  // Get the number of entries

  auto nevents = *df.Count(); 
  cout << "n. dati " << nevents << "\n"; 
  // auto h1 = df.Histo1D("totale_casi");

  // Returns the names of the available columns.

  auto colNames = df.GetColumnNames();
  for (auto &&colName : colNames)  {
    std::cout << colName << ", " ;  
  }
  std::cout << '\n';

  // auto colNames = df.GetColumnNames();
  // for (auto &&colName : colNames)  {
  //   std::cout << "  loader.AddVariable(\"" << colName << "\", 'F');  \t// "
  //     << colName << " (" << df.GetColumnType(colName) << ")\n";
  // }

  // Save selected columns to disk, in a new TTree in new file.
  // "Label" is an ENTRY COLUMN NAME in the HEADER of the DATASET!
  // The "LABEL" values can be (LOOK AT THE DATASET by UNZOOMING) either "b" (BACKGROUND) or "s" (SIGNAL). So this dataset contains ALREADY CLASSIFIED data in sig and bkg!! These data can be used to TRAIN TMVA!!
  df.Filter("Label == \"s\" ")
    .Snapshot("tree_sig", "atlas-higgs-challenge-2014-v2-sig.root");  

  df.Filter("Label == \"b\" ")
    .Snapshot("tree_bkg", "atlas-higgs-challenge-2014-v2-bkg.root");

}
