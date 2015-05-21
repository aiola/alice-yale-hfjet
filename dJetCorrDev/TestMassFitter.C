// root macro to test the class MassFitter

class MassFitter;

void TestMassFitter()
{
  gROOT->LoadMacro("MassFitter.cxx+g");

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  MassFitter *fitter = new MassFitter("MyMassFitter", MassFitter::kGaus, MassFitter::kExpoPower, mpi, 0.18);

  fitter->GetFitFunction()->SetParameter(0,600);
  fitter->GetFitFunction()->SetParameter(1,0.01);
  fitter->GetFitFunction()->SetParameter(2,50);
  fitter->GetFitFunction()->SetParameter(3,0.145);
  fitter->GetFitFunction()->SetParameter(4,0.001);

  Double_t v = fitter->GetFitFunction()->Eval(0.145);
  Printf("%.3f",v);
  
  TCanvas *c = new TCanvas();
  TH1* hist = new TH1F("hist", "hist", 100, 0.12, 0.18);
  hist->Draw();
  hist->GetYaxis()->SetRangeUser(0,v*1.5);
  fitter->Draw("same");

  Double_t bkgInt = fitter->GetBkgFunction()->Integral(mpi, 0.18);
  Double_t totInt = fitter->GetFitFunction()->Integral(mpi, 0.18);

  Printf("bkgInt = %.3f, totInt = %.3f", bkgInt, totInt);
}
