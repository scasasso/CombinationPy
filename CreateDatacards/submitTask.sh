#!/bin/bash


verbosity=0

queue=1nw

taskName=$1
mass=$2
card1=$3
card2=$4


rndSeed=123456
boosFactor=1000

## Combine options
combineOpt="--X-rtd TMCSO_AdaptivePseudoAsimov --X-rtd TMCSO_AdaptivePseudoAsimov_Boost=${boostFactor} --run expected -s ${rndSeed}" #pseudo Asimov
#combineOpt="-t -1 -s 654321" #full Asimov
#combineOpt="-s 654321" #from observation
#combineOpt="-t 1000" #toys

#HMP cards or not?
#combineOpt="${combineOpt} --setPhysicsModelParameters CMS_zz4l_csquared_BSM=1.0,CMS_zz4l_brnew_BSM=0.0 --freezeNuisances CMS_zz4l_csquared_BSM,CMS_zz4l_brnew_BSM,width_mult"


wd=`pwd`
cmssw_dir=$CMSSW_BASE

name1=${card1/hzz4l_/}
name1=${name1/.txt}
name2=${card2/hzz4l_/}
name2=${name2/.txt}

comb1=hzz4l_${name1}_${name2}
comb2=hzz4l_${name2}_${name1}

if [[ -f ${taskName}.lsf ]]; then rm ${taskName}.lsf; else touch ${taskName}.lsf; fi

# Go to the directory and setup CMSSW
echo "echo \"Working directory is: ${wd}\"" >> ${taskName}.lsf
echo "echo \"CMSSW directory is: ${cmssw_dir}\"" >> ${taskName}.lsf
echo "cd ${cmssw_dir}/src;" >> ${taskName}.lsf
echo "eval \`scramv1 runtime -sh\`;" >> ${taskName}.lsf
echo "cd ${wd}" >> ${taskName}.lsf

# Combine the cards with the two orders
echo "combineCards.py ${card1} ${card2} > ${comb1}_${taskName}.txt;" >> ${taskName}.lsf 

# Produce combined workspaces
echo "text2workspace.py -b ${comb1}_${taskName}.txt;" >> ${taskName}.lsf 

echo "echo \"---> Combine will be run with the following options:\"" >> ${taskName}.lsf
echo "echo \"${combineOpt}\"" >> ${taskName}.lsf

# Run combine
echo "combine ${comb1}_${taskName}.root -M Asymptotic -n _hzz4l_${comb1/hzz4l_/}_${mass}_${taskName} -m ${mass} -v ${verbosity} --minosAlgo=stepping --minimizerStrategy=0 --minimizerTolerance=0.0001 --cminFallback Minuit2:0.01 --cminFallback Minuit:0.001 ${combineOpt} 2>&1 | tee log_hzz4l_${comb1}_${taskName}.txt;" >> ${taskName}.lsf  


# If 2nd card is specified run also the other combination
if ! [[ "${card2}" == "" ]]; then
    echo "combineCards.py ${card2} ${card1} > ${comb2}_${taskName}.txt;" >> ${taskName}.lsf 
    echo "text2workspace.py -b ${comb2}_${taskName}.txt;" >> ${taskName}.lsf 
    echo "combine ${comb2}_${taskName}.root -M Asymptotic -n _hzz4l_${comb2/hzz4l_/}_${mass}_${taskName} -m ${mass} -v ${verbosity} --minosAlgo=stepping --minimizerStrategy=0 --minimizerTolerance=0.0001 --cminFallback Minuit2:0.01 --cminFallback Minuit:0.001 ${combineOpt} 2>&1 | tee log_hzz4l_${comb2}_${taskName}.txt;" >> ${taskName}.lsf
fi

chmod u+x ${taskName}.lsf
bsub -q ${queue} -o $wd/lsflog_${taskName}.stdout -e $wd/lsflog_${taskName}.stderr ${taskName}.lsf




