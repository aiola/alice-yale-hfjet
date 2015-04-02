// root macro to run the DecaySimulator

class DecaySimulator;

void runDecaySimulator()
{
  gROOT->LoadMacro("DecaySimulator.cxx+g");

  DecaySimulator* simulator = new DecaySimulator();
  simulator->StartSimulation(10000);
  simulator->SaveTree("simulation.root");
}
