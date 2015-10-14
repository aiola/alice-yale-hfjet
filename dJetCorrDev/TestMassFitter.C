// root macro to test the class MassFitter

class MassFitter;

void TestMassFitter()
{
  if (gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  gROOT->LoadMacro("MassFitter.cxx+g");

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  MassFitter *fitter = new MassFitter("MyMassFitter", MassFitter::kGaus, MassFitter::kExpo, 1.7, 2);
  
  fitter->GetFitFunction()->SetParameter(0,    99.00);
  fitter->GetFitFunction()->SetParameter(1,    -2.00);
  fitter->GetFitFunction()->SetParameter(2,     1.00);
  fitter->GetFitFunction()->SetParameter(3,     1.87);
  fitter->GetFitFunction()->SetParameter(4,     0.01);

  fitter->DisableSig();
  
  TCanvas *c = new TCanvas();
  TH1* hist = new TH1F("hist", "hist", 100, 1.7, 2);
  hist->FillRandom(fitter->GetFitFunction()->GetName(),1000000);
  hist->Sumw2();
  hist->Draw();
  //hist->GetYaxis()->SetRangeUser(0,v*1.5);

  //fitter->GetFitFunction()->FixParameter(0, 1000000);
  fitter->GetFitFunction()->SetParameter(1, -5);

  fitter->Fit(hist, "0 E S");

  fitter->Draw("same");

  Double_t histInt = hist->Integral(hist->GetXaxis()->FindBin(1.7), hist->GetXaxis()->FindBin(2));
  Double_t fitInt = fitter->GetFitFunction()->Integral(1.7,2);

  Printf("hist integral = %.3f, fit integral = %.3f", histInt, fitInt);

  Double_t bkg = fitter->GetBackground();
  Double_t sig = fitter->GetSignal();

  Printf("bkg = %.3f, sig = %.3f", bkg, sig);
}
