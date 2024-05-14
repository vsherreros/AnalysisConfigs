rm Run3_MC_ttBkg.json Run3_MC_ttBkg_redirector.json
build_datasets.py \
    --cfg samples_Run3_ttBKG.json \
    -bs T1_DE_KIT_Disk \
    -bs T1_FR_CCIN2P3_Disk \
    -bs T1_IT_CNAF_Disk \
    -bs T1_IT_CNAF_Tape \
    -bs T1_RU_JINR_Disk \
    -bs T1_UK_RAL_Disk \
    -bs T1_US_FNAL_Disk \
    -bs T2_RU_INR \
    -bs T2_UA_KIPT \
    -bs T2_CH_CSCS \
    -bs T2_US_MIT \
    -bs T2_DE_DESY \
    -bs T2_UK_SGrid_RALPP \
    -bs T2_RU_IHEP \
    -bs T2_RU_ITEP \
    -bs T2_PL_Swierk
# build_datasets.py --cfg samples_Run3_ttBKG.json -ws T2_BE_UCL
# rm Run3_MC_ttBkg.json Run3_MC_ttBkg_redirector.json
# build_datasets.py --cfg samples_Run3_ttBKG.json -bs T1_FR_CCIN2P3_Disk -bs T1_IT_CNAF_Disk -bs T1_IT_CNAF_Tape -bs T1_RU_JINR_Disk -bs T1_UK_RAL_Disk -bs T1_US_FNAL_Disk -bs T2_RU_INR -bs T2_UA_KIPT -bs T2_CH_CSCS -bs T2_US_MIT -bs T2_DE_DESY -bs T2_RU_IHEP -bs T2_UK_London_IC -bs T2_UK_SGrid_RALPP -bs T2_RU_IHEP -bs T2_ES_CIEMAT -bs T2_US_Vanderbilt
