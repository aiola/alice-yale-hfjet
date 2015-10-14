#! /bin/bash

echo "Root: $(which root)"
echo "AliRoot: $(which aliroot)"
echo "Alien: $(which alien-token-info)"

TOKEN=$(alien-token-info)
if [[ $TOKEN != *Token\ is\ still\ valid!* ]]
then
    echo $TOKEN
    echo "Alien token not available. Creating a token for you..."
    alien-token-init saiola
fi

#defining some paths

TrainNumber="${1}"

if [ -z "${TrainNumber}" ]
then
    echo "${0} [trainNumber] [overwrite=0] [year=2015] [dataset=LHC15i2b]  [trainName=Jets_EMC_pp_MC]"
    exit 0
fi

StrippedTrainNumber=${TrainNumber:0:3}

Overwrite="${2}"

if [ -z "${Overwrite}" ]
then
    Overwrite=0
fi

Year="${3}"

if [ -z "${Year}" ]
then
    Year="2015"
fi

Dataset="${4}"

if [ -z "${Dataset}" ]
then
    Dataset="LHC15i2b"
fi

TrainName="${5}"

if [ -z "${TrainName}" ]
then
    TrainName="Jets_EMC_pp_MC"
fi

if [ "${Dataset}" = "LHC15i2b" ]
then
    RunList=( "114786" "114798" "114918" "114920" "114924" "114930" "114931" "115186" "115193" "115310" "115318" "115322" "115328" "115335" "115345" "115393" "115399" "115401" "115414" "115521" "116079" "116081" "116102" "116288" "116402" "116403" "116562" "116571" "116574" "116643" "116645" "117048" "117050" "117052" "117053" "117059" "117060" "117063" "117092" "117099" "117109" "117112" "117116" "117220" "117222" )
    PtHardBinList=( "0" "1" "2" "3" "4" "5" "6" "7" "8" )
else
    echo "Dataset ${Dataset} is not configured for this script!"
    exit 0
fi

AlienPath="alien:///alice/sim/${Year}/${Dataset}"
LocalPath="${JETRESULTS}/${TrainName}_${StrippedTrainNumber}"

echo "Dataset: ${Dataset}"
echo "Train: ${TrainName}"
echo "Train n.: ${TrainNumber}"
echo "Run list: ${RunList[*]}"
echo "Alien path: ${AlienPath}"
echo "Local path: ${LocalPath}"
echo "Overwrite mode: ${Overwrite}"

if [ ! -d "${LocalPath}" ]
then
    echo "Creating directory ${LocalPath}"
    mkdir "${LocalPath}"
fi

cp ./ScaleResults.C $LocalPath/

if [ -e "./${Dataset}.xsec.root" ]
then
    cp ./${Dataset}.xsec.root ${LocalPath}/
fi

cd $LocalPath/

pwd

for PtHardBin in "${PtHardBinList[@]}"
do

  echo "Working on pT hard bin: ${PtHardBin}"

  FileList=""
  
  for Run in "${RunList[@]}"
  do
      if [ ! -d "./${PtHardBin}/${Run}" ]
      then
	  echo "Creating directory ./${PtHardBin}/${Run}"
	  mkdir -p ./${PtHardBin}/${Run}
      fi

      if [ -e "./${PtHardBin}/${Run}/AnalysisResults.root" ] && [ "${Overwrite}" -ge "4" ]
      then
	  echo "Deleting file: ./${PtHardBin}/${Run}/AnalysisResults.root"
	  rm "./${PtHardBin}/${Run}/AnalysisResults.root"
      fi

      if [ ! -e "./${PtHardBin}/${Run}/AnalysisResults.root" ]
      then
	  AlienFile="${AlienPath}/${Run}/${PtHardBin}/PWGJE/${TrainName}/${TrainNumber}/AnalysisResults.root"
	  echo "Copying from alien location: ${AlienFile}"
	  alien_cp ${AlienFile} ./${PtHardBin}/${Run}
      fi

      if [ -e "./${PtHardBin}/${Run}/AnalysisResults.root" ]
      then
	  FileList="${FileList} ./${PtHardBin}/${Run}/AnalysisResults.root"
      fi
  done

  if [ -e "./${PtHardBin}/AnalysisResults.root" ] && [ "${Overwrite}" -ge "3" ]
  then
      echo "Deleting file: ./${PtHardBin}/AnalysisResults.root"
      rm "./${PtHardBin}/AnalysisResults.root"
  fi

  if [ ! -e "./${PtHardBin}/AnalysisResults.root" ]
  then
      echo "Merging for pT hard bin: ${PtHardBin}"
      hadd "./${PtHardBin}/AnalysisResults.root" ${FileList}
  fi

  if [ -e "./${PtHardBin}/ScaledResults.root" ] && [ "${Overwrite}" -ge "2" ]
  then
      echo "Deleting file: ./${PtHardBin}/ScaledResults.root"
      rm "./${PtHardBin}/ScaledResults.root"
  fi
  
  if [ ! -e "./${PtHardBin}/ScaledResults.root" ] && [ -e "./${PtHardBin}/AnalysisResults.root" ]
  then
      root -l -b -q ScaleResults.C+g\(${PtHardBin},\"${Dataset}\"\)
  fi

  if [ -e "./${PtHardBin}/ScaledResults.root" ] && [ "${PtHardBin}" -gt "0" ]
  then
      FinalFileList="${FinalFileList} ./${PtHardBin}/ScaledResults.root"
  fi
done

if [ -e "./AnalysisResults.root" ] && [ "${Overwrite}" -ge "1" ]
then
    echo "Deleting file: ./AnalysisResults.root"
    rm "./AnalysisResults.root"
fi

if [ ! -e "./AnalysisResults.root" ]
then
    echo "Final merging..."
    hadd "./AnalysisResults.root" ${FinalFileList}
fi

echo "Done."

ls
