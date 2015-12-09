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
    
elif [ "${Dataset}" = "LHC15i2c" ]
then
    RunList=( "118506" "118507" "118512" "118518" "118556" "118558" "118560" "118561" "119159" "119161" "119163" "119841" "119842" "119844" "119845" "119846" "119849" "119853" "119856" "119859" "119862" "120067" "120069" "120072" "120073" "120076" "120079" "120244" "120503" "120504" "120505" "120616" "120617" "120671" "120741" "120750" "120758" "120820" "120821" "120822" "120823" "120824" "120825" "120829" "121039" "121040" )
    PtHardBinList=( "0" "1" "2" "3" "4" "5" "6" "7" "8" )

elif [ "${Dataset}" = "LHC15i2d" ]
then
    RunList=( "118506" "118507" "118512" "118518" "118556" "118558" "118560" "118561" "119159" "119161" "119163" "119841" "119842" "119844" "119845" "119846" "119849" "119853" "119856" "119859" "119862" "120067" "120069" "120072" "120073" "120076" "120079" "120244" "120503" "120504" "120505" "120616" "120617" "120671" "120741" "120750" "120758" "120820" "120821" "120822" "120823" "120824" "120825" "120829" "121039" "121040" )
    PtHardBinList=( "0" "1" "2" "3" "4" "5" "6" "7" "8" )

elif [ "${Dataset}" = "LHC15i2e" ]
then
    RunList=( "128366" "128452" "128486" "128494" "128495" "128498" "128503" "128504" "128505" "128506" "128582" "128590" "128592" "128594" "128596" "128605" "128609" "128611" "128615" "128621" "128677" "128678" "128777" "128778" "128819" "128820" "128823" "128824" "128833" "128834" "128835" "128836" "128843" "128850" "128853" "128855" "128913" "129042" "129512" "129513" "129514" "129515" "129516" "129519" "129520" "129521" "129523" "129524" "129525" "129527" "129528" "129536" "129540" "129586" "129587" "129599" "129639" "129641" "129647" "129650" "129651" "129652" "129653" "129659" "129666" "129723" "129725" "129726" "129729" "129734" "129735" "129736" "129738" "129742" "129744" "129959" "129960" "129961" "129962" "129966" "129983" "130149" "130151" "130157" "130158" "130168" "130172" "130178" "130342" "130343" "130354" "130356" "130358" "130360" "130375" "130479" "130480" "130481" "130517" "130519" "130520" "130524" "130526" "130601" "130608" "130609" "130620" "130621" "130623" "130628" "130696" "130704" "130793" "130795" "130798" "130799" "130802" "130803" "130804" "130834" "130840" "130842" "130844" "130847" "130848" "130850" )
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
cp ./MergeFiles.C $LocalPath/

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
      root -l -b -q MergeFiles.C+g\(\""./${PtHardBin}/AnalysisResults.root"\",\""${FileList}"\"\)
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
    root -l -b -q MergeFiles.C+g\(\""./AnalysisResults.root"\",\""${FinalFileList}"\"\)
fi

echo "Done."

ls
