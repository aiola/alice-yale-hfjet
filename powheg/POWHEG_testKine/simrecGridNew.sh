#!/bin/bash

echo "------------------ job starts ---------------------"
date +%s

export XROOTD_MAXWAIT=10
export GCLIENT_SERVER_LIST="pcapiserv01.cern.ch:10000|pcapiserv02.cern.ch:10000|pcapiserv04.cern.ch:10000|pcapiserv05.cern.ch:10000|pcapiserv06.cern.ch:10000"
export GCLIENT_COMMAND_MAXWAIT=60
export GCLIENT_COMMAND_RETRY=10
export GCLIENT_SERVER_RESELECT=4
export GCLIENT_SERVER_RECONNECT=2
export GCLIENT_RETRY_DAMPING=1.2
export GCLIENT_RETRY_SLEEPTIME=1
export GRID_TOKEN=OK

export PWHGPROC=$1
NEVT=$2

echo "Untar LHAPDF datasets..."
tar zxf lhapdf6.tgz
export PATH="./LHAPDF/6.0.6a0/bin:$PATH"
export LD_LIBRARY_PATH="./LHAPDF/6.0.6a0/lib:$LD_LIBRARY_PATH"
export LHAPDF_DATA_PATH=./LHAPDF/6.0.6a0/share/LHAPDF 

echo "Running $1 MC production on:"
uname -a

echo "Running POWHEG..."

#if [ "$PWHGPROC" = "dijet" ] 
#then
#echo "iseed $RANDOM" >> powheg$1.input
#else
#echo "randomseed $RANDOM" >> powheg$1.input
#fi
echo "iseed $RANDOM" >> powheg$1.input
echo "numevts $NEVT" >> powheg$1.input

cat powheg$1.input

cp powheg$1.input powheg.input
chmod u+x pwhg_main_$1
./pwhg_main_$1 < powheg.input &> pwhg.log

export CONFIG_SEED=$RANDOM
echo "Setting PYTHIA seed to $CONFIG_SEED"

NE="$(echo "scale=0; $NEVT-0.1*$NEVT" | bc)"
export PYNEVT="$(echo "($NE+0.5)/1" | bc)"

echo "Untar OCDB..."
tar zxf OCDB.tgz

echo "Running simulation..."
aliroot -b >& sim.log << EOF
.x simSimpleGrid.C($PYNEVT,"./ConfigSimpleGrid.C")
.q
EOF

echo "Done"
echo "...see results in the log file"

ls -l

echo "############### SCRATCH SIZE ####################"
du -sh

echo "------------------ job ends ----------------------"
date +%s

