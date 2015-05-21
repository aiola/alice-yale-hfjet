// root macro to test the class MassFitter

class MassFitter;

void TestMassFitter()
{
  gROOT->LoadMacro("MassFitter.cxx+g");

  Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  MassFitter *fitter = new MassFitter("MyMassFitter", MassFitter::kGaus, MassFitter::kExpoPower, 0.105, 0.185);
  
  fitter->GetFitFunction()->SetParameter(0,8.80000e+01);
  fitter->GetFitFunction()->SetParameter(1,2.34355e+01);
  fitter->GetFitFunction()->SetParameter(2,3.17321e+01);
  fitter->GetFitFunction()->SetParameter(3,1.45430e-01);
  fitter->GetFitFunction()->SetParameter(4,6.87373e-04);

  Double_t v = fitter->GetFitFunction()->Eval(0.14543);
  Printf("%.3f",v);
  
  TCanvas *c = new TCanvas();
  TH1* hist = new TH1F("hist", "hist", 100, 0.105, 0.185);
  hist->Draw();
  hist->GetYaxis()->SetRangeUser(0,v*1.5);
  fitter->Draw("same");

  Double_t bkgInt = fitter->GetBkgFunction()->Integral(mpi, 0.185);
  Double_t totInt = fitter->GetFitFunction()->Integral(mpi, 0.185);

  Printf("bkgInt = %.3f, totInt = %.3f", bkgInt, totInt);
}
