// root macro to use the class ProjectDJetCorr

class ProjectDJetCorr;

void runProjectDJetCorr()
{
  gROOT->LoadMacro("ProjectDJetCorr.cxx+g");

  ProjectDJetCorr* projDjet = new ProjectDJetCorr("train");
  projDjet->SetInputListName("D0", "Charged", "R040");
  projDjet->SetUseTest();
  projDjet->SetInputPath("~/");
  projDjet->SetOverwrite(kTRUE);

  projDjet->Run();
}
