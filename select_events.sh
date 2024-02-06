#!/usr/bin/env bash
filepath_tttt='/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer20UL18NanoAODv9/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2520000'
filepath_tttW='/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer20UL18NanoAODv9/TTTW_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2'
filepath_tttJ='/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer20UL18NanoAODv9/TTTJ_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2'

tttt_filenames[0]='1323716F-65CC-1B46-96FE-2E660766A235.root'
tttt_filenames[1]='2D36D54E-8CCB-E144-90C5-3312F0805C7F.root'
tttt_filenames[2]='3876D6A8-8D88-424D-8DE2-408C95A1AF78.root'
tttt_filenames[3]='5F6252A3-1634-404F-A493-2A36345FDA1B.root'
tttt_filenames[4]='6AAF523E-FCEE-ED42-B161-C6204614BDF8.root'
tttt_filenames[5]='9135C82E-46A9-0B4A-AE88-2149604E80F8.root'
tttt_filenames[6]='AED740DE-A0C3-9248-846F-3B6EA9B689EE.root'
N_tttt_files=7
tttJ_filenames[0]='2520000/DBB31367-6197-CC4F-AE0C-CED640EE24A2.root'
tttJ_filenames[1]='2530000/46BC5096-1D58-1A45-81F3-F3E18CD4695A.root'
tttJ_filenames[2]='260000/4668BDC5-735F-D14C-8B2B-5DC437CA5B25.root'
tttJ_filenames[3]='270000/11CCD6D3-B6C1-8C45-9B20-08BF279D6BA4.root'
N_tttJ_files=4
tttW_filenames[0]='2520000/9660E5EF-E305-7546-B26E-9A759DBC1976.root'
tttW_filenames[1]='2530000/07E85166-4C72-3B4E-9572-2D070F2EF385.root'
tttW_filenames[2]='260000/07350478-4397-054C-B6F8-9F6979FF7357.root'
tttW_filenames[3]='260000/ECACC29E-698C-EE42-A460-BAD444727AFA.root'
N_tttW_files=4
echo How many tttt files?
read num_tttt

if [ "$num_tttt" == 'all' ]; then
    num_tttt=7
    echo "Selecting $num_tttt tttt files"
fi
echo "Selecting $num_tttt tttt files"

if [ $num_tttt != 0 ]; then
    for ((i=0; i<$N_tttt_files; i++)); do
        echo $i
        python3 event_selector.py -f /pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer20UL18NanoAODv9/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2520000/${tttt_filenames[$i]} -p tttt
    done
fi
# echo Do ttt processes?
# read ttt
ttt='yes'
if [ "$ttt"=='yes' ]; then
    for ((i=0; i<$N_tttJ_files; i++)); do
        echo $i
        python3 event_selector.py -f $filepath_tttJ/${tttJ_filenames[$i]} -p tttJ
    done
    for ((i=0; i<$N_tttW_files; i++)); do
        echo $i
        python3 event_selector.py -f $filepath_tttW/${tttW_filenames[$i]} -p tttW
    done
fi
