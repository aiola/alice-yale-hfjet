#include <TString.h>
#include <TH1F.h>
#include <THnBase.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>

using namespace std;

TString xsecHistName;
TString scalingHistName;
TString trialsHistName;

TFile* OpenFile(const char* localPath, Int_t pt);
Double_t GetTotalXsec(const char* dataset);
TH1* CreateTH1F(const char* name, const char* title);
void ProcessDir(TDirectoryFile *file, Double_t scaling, TDirectoryFile *outfile, const char* type, const char* exclude);
Bool_t ProcessObject(TObject *obj, Double_t scaling);

void ScaleResults(const char* localPath, Int_t ptHard, const char* dataset, const char* type="", const char* exclude="")
{
  TH1::AddDirectory(kFALSE);

  xsecHistName = Form("%s_xsec",dataset);
  scalingHistName = Form("%s_scaling",dataset);
  trialsHistName = Form("%s_trials",dataset);

  Double_t totalxsec = GetTotalXsec(dataset);

  TFile *xsec_file = 0;
  TH1 *newxsec = 0;
  TH1 *newscalingFactors = 0;
  TH1 *newtrials = 0;
  TString newxsecFileName(Form("%s/%s.newxsec.root", localPath, dataset));
  if (!gSystem->AccessPathName(newxsecFileName)) {
    xsec_file = TFile::Open(newxsecFileName);
    if (xsec_file && !xsec_file->IsZombie()) {    
      TH1* temp = dynamic_cast<TH1*>(xsec_file->Get(xsecHistName));
      if (temp)
	newxsec = static_cast<TH1*>(temp->Clone(xsecHistName));

      TH1* temp2 = dynamic_cast<TH1*>(xsec_file->Get(scalingHistName));
      if (temp2)
	newscalingFactors = static_cast<TH1*>(temp2->Clone(scalingHistName));

      TH1* temp3 = dynamic_cast<TH1*>(xsec_file->Get(trialsHistName));
      if (temp3)
	newtrials = static_cast<TH1*>(temp3->Clone(trialsHistName));

      xsec_file->Close();
      delete xsec_file;
    }
  }
  xsec_file = 0;

  TFile *file = OpenFile(localPath, ptHard);

  TList *keys = dynamic_cast<TList*>(file->GetListOfKeys()->Clone("mykeys"));
  if (!keys) {
    cout << "Could not get keys" << endl;
    return;
  }

  Int_t nlists = keys->GetEntries();

  Double_t xsec = 1;
  Double_t ntrials = 1;

  Bool_t xSecOk = kFALSE;
  const Int_t histbin = ptHard + 1;

  for (Int_t j = 0; j < nlists; j++) {
    TString listname(keys->At(j)->GetName());
    
    TObject *orig = file->Get(listname);

    TList *listOrig = dynamic_cast<TList*>(orig);
    if (listOrig) {
      Printf("Looking in list %s", listOrig->GetName());

      TH1* histntrials = dynamic_cast<TH1*>(orig->FindObject("fHistTrials")); 
      TH1* histxsec = dynamic_cast<TH1*>(orig->FindObject("fHistXsection"));

      if (histntrials && histntrials) {
	xsec = histxsec->GetBinContent(histbin);
	ntrials = histntrials->GetBinContent(histbin);
	xSecOk = kTRUE;
      }
      delete listOrig;
      listOrig = 0;
    }

    TDirectoryFile *dirOrig = dynamic_cast<TDirectoryFile*>(orig);
    if (dirOrig) {
      Printf("Looking in directory %s", dirOrig->GetName());
      TString listName(dirOrig->GetName());
      listName.ReplaceAll("PWG4_","chist0");
      Printf("Looking for list %s", listName.Data());
      TList *listOrig2 = dynamic_cast<TList*>(dirOrig->Get(listName));
      if (listOrig2) {
	Printf("Looking in list %s", listOrig2->GetName());
	TH1* histntrials = dynamic_cast<TH1*>(listOrig2->FindObject("fh1Trials")); 
	TH1* histxsec = dynamic_cast<TH1*>(listOrig2->FindObject("fh1Xsec"));

	if (histntrials && histntrials) {
	  xsec = histxsec->GetBinContent(1);
	  ntrials = histntrials->GetBinContent(1);
	  xSecOk = kTRUE;
	}

	delete listOrig2;
	listOrig2 = 0;
      }
    }

    if (xSecOk) break;
  }

  if (!xSecOk) { 
    Printf("I could not find any x section and trials info. Returning...");
    return;
  }

  Double_t scaling = xsec / ntrials * totalxsec; 

  if (!newxsec) newxsec = CreateTH1F(xsecHistName, "xsection");
  if (!newscalingFactors) newscalingFactors = CreateTH1F(scalingHistName, "scaling");
  if (!newtrials) newtrials = CreateTH1F(trialsHistName, "trials");

  Printf("Using pt hard bin: %s", newtrials->GetXaxis()->GetBinLabel(histbin));
  Printf("Scaling = %f, xsec = %f, trials = %f, totalxsec = %f", scaling, xsec, ntrials, totalxsec);

  newxsec->SetBinContent(histbin, xsec);
  newtrials->SetBinContent(histbin, ntrials);
  newscalingFactors->SetBinContent(histbin, scaling);

  if (newxsec) {
    xsec_file = TFile::Open(newxsecFileName,"RECREATE");
    if (!xsec_file || xsec_file->IsZombie()) {
      Printf("Could not store xsec file.");
    }
    else {
      xsec_file->cd();
      newxsec->Write();
      newscalingFactors->Write();
      newtrials->Write();
      xsec_file->Close();
      delete xsec_file;
    }
  }
  xsec_file = 0;

  TString fileoutName = Form("%s/%d/ScaledResults.root", localPath, ptHard);
    
  TFile *outfile = TFile::Open(fileoutName, "RECREATE");

  ProcessDir(file, scaling, outfile, type, exclude);

  file->Close();
  outfile->Close();

  //delete file; // for some reason this fails, to be investigated
  file = 0;
  
  delete outfile;
  outfile = 0;
}

void ProcessDir(TDirectoryFile *file, Double_t scaling, TDirectoryFile *outfile, const char* type, const char* exclude)
{
  Printf("Working on direcory %s...", file->GetName());

  TList *keys = dynamic_cast<TList*>(file->GetListOfKeys()->Clone("mykeys"));
  if (!keys) {
    cout << "Could not get keys" << endl;
    return;
  }
  
  Int_t nobjects = keys->GetEntries();
  TObject **objects = new TObject*[nobjects];
  Int_t nobjectsAcc = 0;

  for (Int_t j = 0; j < nobjects; j++) {
    TString objname(keys->At(j)->GetName());
    if (type && strcmp(type,"")!=0 && objname.BeginsWith("AliJetResponseMaker") && !objname.Contains(type)) {
      cout << "Skipping object " << objname.Data() << "..." << endl;
      continue;
    }
    if (exclude && strcmp(exclude,"")!=0 && objname.Contains(exclude)) {
      cout << "Skipping object " << objname.Data() << "..." << endl;
      continue;
    }
    cout << "Now working with object " << objname.Data() << endl;
    TObject *orig = file->Get(objname);

    TDirectoryFile *dirOrig = dynamic_cast<TDirectoryFile*>(orig);
    if (dirOrig) {
      /// create the outsubfile
      outfile->cd();
      TDirectoryFile *outsubfile = new TDirectoryFile(objname, objname);
      ProcessDir(dirOrig, scaling, outsubfile, type, exclude);

      // automatically added to the file

      continue;
    }

    TList *listOrig = dynamic_cast<TList*>(orig);
    if (listOrig) {
      
      objects[nobjectsAcc] = orig;
      nobjectsAcc++;

      for (Int_t k = 0; k < listOrig->GetEntries(); k++) {
	Bool_t success = ProcessObject(listOrig->At(k), scaling);
	if (!success) listOrig->RemoveAt(k);
      }

      continue;
    }

    Bool_t success = ProcessObject(orig, scaling);

    if (success) {
      objects[nobjectsAcc] = orig;
      nobjectsAcc++;
    }
  }     

  Printf("Done with this directory, now saving...");

  outfile->cd();
  for (Int_t j = 0; j < nobjectsAcc; j++) {
    objects[j]->Write(objects[j]->GetName(), TObject::kSingleKey);
  }

  Printf("Saving completed.");
}


Bool_t ProcessObject(TObject *obj, Double_t scaling)
{
  TH1* histo = dynamic_cast<TH1*>(obj);
  if (histo) {
    histo->Sumw2();
	
    TString histoName(histo->GetName());
    if (histoName != "fHistTrialsAfterSel" && histoName != "fHistEventsAfterSel" &&
	histoName != "fHistTrials" && histoName != "fHistXsection" &&
	histoName != "fHistEvents") {
      histo->Scale(scaling);
    }
    return kTRUE;
  }

  THnBase *hn = dynamic_cast<THnBase*>(obj);
  if (hn) {
    hn->Sumw2();
    for (Int_t ibin = 0; ibin < hn->GetNbins(); ibin++) {
      Double_t c = hn->GetBinContent(ibin);
      hn->SetBinContent(ibin, c*scaling);
      hn->SetBinError2(ibin, c*scaling*scaling);	    
    }
    //hn->Scale(scaling); // not safe: if Sumw2() was not called before filling the THn errors are not properly propagated
    return kTRUE;
  }

  Printf("Object %s not recognized!", obj->GetName());

  return kFALSE;
}


    
TFile* OpenFile(const char* localPath, Int_t pt)
{
  static TFile file;
  static Int_t ptopen  = -1;

  if (pt == ptopen) return &file;

  if (file.IsOpen()) file.Close("R");

  TString fullpath = Form("%s/%d/AnalysisResults.root", localPath, pt);

  new (&file) TFile(fullpath);

  if (file.IsOpen()) {
    ptopen = pt;
    cout << "File " << fullpath << " open" << endl;
  }
  else {
    cout << "Could not open file " << fullpath << endl;
  }

  return &file;
}

Double_t GetTotalXsec(const char* dataset)
{
  Double_t totalxsec = 1;

  TFile *xsec_file = 0;
  TString xsecFileName(Form("%s.xsec.root",dataset));
  if (!gSystem->AccessPathName(xsecFileName)) {
    xsec_file = TFile::Open(xsecFileName);
    if (!xsec_file || xsec_file->IsZombie()) {
      Printf("Unable to open x sec file. The spectra will not be scaled...");
    }
    else {
      Printf("File %s open...", xsec_file->GetName());
      TH1F *xsec = dynamic_cast<TH1F*>(xsec_file->Get(xsecHistName));
      if (xsec) {
	totalxsec = xsec->Integral();
      }
      else {
	Printf("Unable to get x sec histogram. The spectra will not be scaled...");
      }
      delete xsec;
      xsec = 0;
      xsec_file->Close();
      delete xsec_file;
    }
  }
  else {
    Printf("Unable to find x sec file. The spectra will not be scaled...");
  }
  xsec_file = 0;

  return totalxsec;
}
 
TH1* CreateTH1F(const char* name, const char* title) 
{
  TH1* result = new TH1F(name, title, 11, 0, 11);
  result->GetXaxis()->SetTitle("p_{T} hard bin");
  result->GetYaxis()->SetTitle(title);

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
    
  for (Int_t i = 1; i < 12; i++) 
    result->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));

  return result;
}
  
