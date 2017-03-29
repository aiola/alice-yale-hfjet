# How to run the analysis

## Detector response

Start with generating the detector response and reconstruction efficiency:

`./runResponseAnalysis.py LHC15i2response_Train1073.yaml `

When this is done, it should be repeated using results from the previous to apply the reconstruction efficiency:

`./runResponseAnalysis.py LHC15i2response_Train1073_efficiency.yaml`

## Reflection templates

Generate reflection templates without reconstruction efficiency:

`./runDataAnalysis.py LHC15i2reflections_Train1073.yaml`

`cd ../rawYieldUnc`

`./GenerateReflTemplate.py ../DMesonJetAnalysis/LHC15i2reflections_Train1073.yaml`

and with reconstruction efficiency:

`./runDataAnalysis.py LHC15i2reflections_Train1073_efficiency.yaml`

`cd ../rawYieldUnc`

`./GenerateReflTemplate.py ../DMesonJetAnalysis/LHC15i2reflections_Train1073_efficiency.yaml`

## Reconstruct data

Analyze data without reconstruction efficiency:

`./runDataAnalysis.py LHC10analysis_Train883.yaml`

and with reconstruction efficiency:

`./runDataAnalysis.py LHC10analysis_Train883_efficiency.yaml`

## Raw yield uncertainty

Go in `rawYieldUnc` folder:

`cd ../rawYieldUnc`.

Then execute:

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train883.yaml --refl DoubleGaus`

`./ExtractDZeroJetRawYieldReflectionUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train883.yaml`

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train883_efficiency.yaml --refl DoubleGaus`

`./ExtractDZeroJetRawYieldReflectionUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train883_efficiency.yaml`

## Unfolding

`./runDataUnfolding.py LHC10_Train883_LHC15i2_Train1073_efficiency.yaml`

`./runDataUnfolding.py LHC10_Train883_LHC15i2_Train1073.yaml`

To run full systematics:

`./runDataUnfolding.py LHC10_Train883_LHC15i2_Train1073_efficiency.yaml --fd-syst --ry-syst`