// root macro to test the class MassFitter

class MassFitter;

void TestMassFitter()
{
  gROOT->LoadMacro("MassFitter.cxx+g");

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  MassFitter *fitter = new MassFitter("MyMassFitter", MassFitter::kGaus, MassFitter::kExpo, 1.7, 2);
  
  fitter->GetFitFunction()->SetParameter(0,    99.00);
  fitter->GetFitFunction()->SetParameter(1,    -2.00);
  fitter->GetFitFunction()->SetParameter(2,     1.00);
  fitter->GetFitFunction()->SetParameter(3,     1.87);
  fitter->GetFitFunction()->SetParameter(4,     0.01);

  fitter->DisableSig();

  Double_t v = fitter->GetFitFunction()->Eval(1.8);

  Printf("%.3f",v);
  
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

  Double_t bkgInt = fitter->GetBkgFunction()->Integral(1.7, 2);
  Double_t totInt = fitter->GetFitFunction()->Integral(1.7, 2);

  Printf("bkgInt = %.3f, totInt = %.3f", bkgInt, totInt);

  Double_t bkg = fitter->GetBackground();
  Double_t sig = fitter->GetSignal();

  Printf("bkg = %.3f, sig = %.3f", bkg, sig);
}
