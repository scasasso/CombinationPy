#!/bin/bash

TYPE=$1
MASS=$2
TOOL=$3
CPRIME=$4
BRNEW=$5
OPTIONS=$6

OPTIONS="$OPTIONS --setPhysicsModelParameters CMS_zz4l_csquared_BSM=${CPRIME},CMS_zz4l_brnew_BSM=${BRNEW} --freezeNuisances CMS_zz4l_csquared_BSM,CMS_zz4l_brnew_BSM"

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd` with options $TYPE $MASS $OPTIONS

eval `scram runtime -sh`

if [[ "$TYPE" == "ASCLS" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_ASCLS_BSM.sh $OPTIONS -l $MASS comb_hzz4l.root; fi;
    if [[ "$TOOL" == "lands" ]]; then bash make_ASCLS_lands.sh $OPTIONS $MASS comb_hzz4l.txt; fi;

elif [[ "$TYPE" == "PLP" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_PLC.sh $OPTIONS -P $MASS comb_hzz4l.root; fi;
    if [[ "$TOOL" == "lands" ]]; then bash make_PLC_lands.sh $OPTIONS -P $MASS comb_hzz4l.txt; fi;

elif [[ "$TYPE" == "PLPE" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_PLC.sh $OPTIONS --PE $MASS comb_hzz4l.root; fi;
    if [[ "$TOOL" == "lands" ]]; then bash make_PLC_lands.sh $OPTIONS --PE $MASS comb_hzz4l.txt; fi;

elif [[ "$TYPE" == "PLS" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_PLC.sh $OPTIONS -S $MASS comb_hzz4l.root; fi
    if [[ "$TOOL" == "lands" ]]; then bash make_PLC_lands.sh $OPTIONS -S $MASS comb_hzz4l.txt; fi;

elif [[ "$TYPE" == "PLSE" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_PLC.sh $OPTIONS --SE $MASS comb_hzz4l.root; fi;
    if [[ "$TOOL" == "lands" ]]; then bash make_PLC_lands.sh $OPTIONS --SE $MASS comb_hzz4l.txt; fi;


elif [[ "$TYPE" == "ML" ]]; then
    
    if [[ "$TOOL" == "combine" ]]; then bash make_ML.sh $OPTIONS $MASS comb_hzz4l.root; fi;
    if [[ "$TOOL" == "lands" ]]; then bash make_ML_lands.sh $OPTIONS $MASS comb_hzz4l.txt; fi;


else 
    echo "Unkown Type: $TYPE"
    exit;

fi



