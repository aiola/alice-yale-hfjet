#include <TFileMerger.h>
#include <Riostream.h>
#include <TGrid.h>
#include <TString.h>
#include <TObjArray.h>

#include <string>
#include <fstream>

Int_t MergeFiles(const char* output, const char* fileList, const char* skipList = "", const char* acceptList = "", Int_t n=2)
{
  TGrid::Connect("alien://");

  TFileMerger merger(kFALSE);
  merger.OutputFile(output);
  merger.SetMaxOpenedFiles(n);

  ifstream in(fileList);

  Int_t nFiles = 0;
  while (in.good()) {
    std::string f;

    in >> f;

    if (f.length() == 0) continue;

    Printf("Adding file %s", f.c_str());
    merger.AddFile(f.c_str());

    nFiles++;
  }

  in.close();

  Printf("Total number of files is %d", nFiles);

  UInt_t mode = TFileMerger::kAllIncremental;

  TObject* obj = 0;

  TString strSkipList(skipList);
  TObjArray *skipArray = strSkipList.Tokenize(" ");

  TString strAcceptList(acceptList);
  TObjArray *acceptArray = strAcceptList.Tokenize(" ");

  if (skipArray && skipArray->GetEntriesFast() > 0) {
    mode = mode | TFileMerger::kSkipListed;
    if (acceptArray && acceptArray->GetEntriesFast() > 0) {
      Printf("Accept list is being ignored!!!");
    }
    TIter nextSkip(skipArray);
    
    while ((obj = nextSkip())) {
      merger.AddObjectNames(obj->GetName());
    }
  }
  else if (acceptArray && acceptArray->GetEntriesFast() > 0) {
    mode = mode | TFileMerger::kOnlyListed;
    TIter nextAccept(acceptArray);
    while ((obj = nextAccept())) {
      merger.AddObjectNames(obj->GetName());
    }
  }

  delete skipArray; skipArray = 0;
  delete acceptArray; acceptArray = 0;

  merger.PrintFiles("");
  Int_t r = merger.PartialMerge(mode);

  if (!r) Printf("Merge error!");

  return r;
}
