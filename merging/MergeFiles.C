#include <TFileMerger.h>
#include <TString.h>
#include <TObjArray.h>

void MergeFiles(TString output, TString list, Int_t n = 2)
{
  TFileMerger merger(kFALSE);

  TObjArray *fileList = list.Tokenize(" ");
  TIter next(fileList);
  TObject *fileName = 0;

  while ((fileName = next())) {
    merger.AddFile(fileName->GetName());
  }

  delete fileList;
  fileList = 0;

  merger.AddObjectNames("AnalysisCutsDzero");
  merger.AddObjectNames("AnalysisCutsDStar");
  
  merger.OutputFile(output);
  merger.SetMaxOpenedFiles(n);
  merger.PrintFiles("");
  Bool_t r = merger.PartialMerge(TFileMerger::kAllIncremental | TFileMerger::kSkipListed);

  if (!r) Printf("Merge error!");
}
