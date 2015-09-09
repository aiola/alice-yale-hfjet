#!/bin/bash

free
echo _____________________________________________
echo "HOME IS $HOME"
ls $HOME
length=`echo $HOME |wc -c`
if   (( $length >= 100 )) ;
then
     echo "WARNING: The home directory $HOME is longer than 100 char"
     OLDHOME=$HOME
     NEWHOME="/tmp/alien_home_dir.${ALIEN_PROC_ID}"
     echo "CHANGING HOME TO $NEWHOME"
     ln -s "$HOME" "$NEWHOME"
     export HOME=$NEWHOME
fi
echo _____________________________________________

export PRODUCTION_METADATA="$ALIEN_JDL_LPMMETADATA"

if [ "$1" = "OCDB" ]; then
    echo "Setting env for generating OCDB.root"
    
    export OCDB_SNAPSHOT_CREATE="kTRUE"
    export OCDB_SNAPSHOT_FILENAME="OCDB.root"
    
    touch OCDB.generating.job
    
    shift
fi

if [ -f simrun.sh ]; then
    echo "Calling simrun.sh with $*"
    chmod a+x simrun.sh
    ./simrun.sh $*

    error=$?

    if [ $error -ne 0 ]; then
        echo "*! Command 'simrun.sh $*' exited with error code $error"
        echo "Command 'simrun.sh $*' exited with error code $error" > validation_error.message
        exit $error
    fi
else
    echo "Executing root.exe with $*"
    time root.exe -b -q -x $*

    error=$?

    if [ $error -ne 0 ]; then
        echo "*! Command 'root.exe -b -q -x $*' exited with error code $error"
        echo "Command 'root.exe -b -q -x $*' exited with error code $error" > validation_error.message
        exit 2
    fi
fi

if [ ! -f AliESDs.root ]; then
    echo "*! Could not find AliESDs.root, the simulation/reconstruction chain failed!"
    echo "Could not find AliESDs.root, the simulation/reconstruction chain failed!" > validation_error.message
    exit 2
fi

if [ -f QAtrainsim.C ]; then
    echo "QAresults.root" >> validation_extrafiles.list

    while [ ! -z "$1" ]; do
        if [ "$1" == "--run" ]; then
            runno="$2"
            break
        fi

        shift
    done

    echo "* Running the QA train for run $runno"
    time aliroot -b -q -x QAtrainsim.C\($runno\) &> qa.log

    error=$?

    for file in *.stat; do
        mv "$file" "$file.qa"
    done

    if [ $error -ne 0 ]; then
        echo "*! QA exited with error code $error"
        echo "QA exited with error code $error" > validation_error.message
        exit 100
    fi
fi

if [ -f AODtrainsim.C ]; then
    echo "AliAOD.root" >> validation_extrafiles.list

    rm -f outputs_valid &>/dev/null

    echo "* Running the FILTERING train..."
    time aliroot -b -q -x AODtrainsim.C\(\) &> aod.log

    error=$?

    for file in *.stat; do
        mv $file $file.aod
    done

    if [ $error -ne 0 ]; then
        echo "*! AOD filtering exited with error code $error"
        echo "AOD filtering exited with error code $error" > validation_error.message
        exit 200
    fi
fi

if [ ! -z "$NEWHOME" ]; then
    echo "DELETING $NEWHOME"
    export HOME=$OLDHOME
    rm -rf $NEWHOME
fi

exit 0
