const char *anatype = "AOD";

void lego_train()
{
// Analysis using AOD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// load base root libraries
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");

   gSystem->Load("libMinuit");

// Add aditional libraries
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");
   gSystem->Load("libTree");
   gSystem->Load("libMinuit");
   gSystem->Load("libProof");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
// Load analysis framework libraries
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libOADB");
   gSystem->Load("libCORRFW");

// include path
   TString intPath = gInterpreter->GetIncludePath();
   TObjArray *listpaths = intPath.Tokenize(" ");
   TIter nextpath(listpaths);
   TObjString *pname;
   while ((pname=(TObjString*)nextpath())) {
      TString current = pname->GetName();
      if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
      gSystem->AddIncludePath(current);
   }
   if (listpaths) delete listpaths;
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());


// analysis source to be compiled at runtime (if any)

// read the analysis manager from file
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("train.root");
   if (!mgr) return;
   mgr->SetDebugLevel(100);
   mgr->PrintStatus();
   gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
   TChain* chain = CreateAODChain("fileLists/files_LHC10c_pass4_AOD.txt", 2, 0, kFALSE , "AliAOD.VertexingHF.root");
   mgr->StartAnalysis("local", chain, 5000);
   timer.Stop();
   timer.Print();
}

