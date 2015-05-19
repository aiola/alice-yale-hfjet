#include <fstream>
#include <TString.h>
#include <Riostream.h>
#include <algorithm>
#include <TArrayI.h>
#include <TObjString.h>
#include <TObjArray.h>

TArrayI* GetRunList(const char* listname);

void CompareRunLists(const char* list1 = "LHC10d_flat", const char* list2 = "LHC10d_nonflat")
{
  TArrayI* runlist1 = GetRunList(list1);
  TArrayI* runlist2 = GetRunList(list2);

  TArrayI* runlistOR = new TArrayI(runlist1->GetSize()+runlist2->GetSize());
  TArrayI* runlistAND = new TArrayI(TMath::Min(runlist1->GetSize(),runlist2->GetSize()));

  Int_t i1 = 0;
  Int_t i2 = 0;
  Int_t nAND = 0;
  Int_t nOR = 0;

  TString listOR;
  TString listAND;

  Printf("%s %s", list1, list2);
  
  while (i1 < runlist1->GetSize() && i2 < runlist2->GetSize()) {
    if (runlist1->At(i1) == runlist2->At(i2)) {
      Printf("%d %d", runlist1->At(i1), runlist2->At(i2));
      (*runlistAND)[nAND] = runlist1->At(i1);
      (*runlistOR)[nOR] = runlist1->At(i1);
      nAND++;
      nOR++;

      listOR += runlist1->At(i1);
      listOR += ", ";
      listAND += runlist1->At(i1);
      listAND += ", ";

      i1++;
      i2++;
      
      continue;
    }
    while (runlist1->At(i1) < runlist2->At(i2) && i1 < runlist1->GetSize()) {
      Printf("%d n/a", runlist1->At(i1));
      (*runlistOR)[nOR] = runlist1->At(i1);
      nOR++;

      listOR += runlist1->At(i1);
      listOR += ", ";

      i1++;
    }
    while (runlist2->At(i2) < runlist1->At(i1) && i2 < runlist2->GetSize()) {
      Printf("n/a    %d", runlist2->At(i2));
      (*runlistOR)[nOR] = runlist2->At(i2);
      nOR++;

      listOR += runlist2->At(i2);
      listOR += ", ";

      i2++;
    }
  }
  while (i1 < runlist1->GetSize()) {
    Printf("%d n/a", runlist1->At(i1));
    (*runlistOR)[nOR] = runlist1->At(i1);
    nOR++;

    listOR += runlist1->At(i1);
    listOR += ", ";

    i1++;
  }
  while (i2 < runlist2->GetSize()) {
    Printf("n/a    %d", runlist2->At(i2));
    (*runlistOR)[nOR] = runlist2->At(i2);
    nOR++;

    listOR += runlist2->At(i2);
    listOR += ", ";

    i2++;
  }
  
  listAND.Remove(listAND.Length()-2);
  listOR.Remove(listOR.Length()-2);

  Printf("Total runs in list %s: %d", list1, runlist1->GetSize());
  Printf("Total runs in list %s: %d", list2, runlist2->GetSize());
  Printf("Total shared runs: %d", nAND);
  Printf("Runs: %s", listAND.Data());
  Printf("Total runs: %d", nOR);
  Printf("Runs: %s", listOR.Data());
}

TArrayI* GetRunList(const char* listname)
{
  const Int_t NMAX = 10000;
  char buffer[NMAX];
  
  TString fname(listname);
  fname += ".txt";

  ifstream file;

  file.open(fname);
  file.read(buffer, NMAX);
  buffer[file.gcount()] = '\0';
  TString runliststr(buffer);
  file.close();

  TObjArray* runlist = runliststr.Tokenize(" ");
  TArrayI* runlistint = new TArrayI(runlist->GetEntriesFast());

  TIter next(runlist);
  TObjString* runstr = 0;
  Int_t i = 0;
  while ((runstr = static_cast<TObjString*>(next()))) {
    (*runlistint)[i] = runstr->GetString().Atoi();
    i++;
  }

  std::sort(runlistint->GetArray(), runlistint->GetArray()+runlistint->GetSize());

  return runlistint;
}
