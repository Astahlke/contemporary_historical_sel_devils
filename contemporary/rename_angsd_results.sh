## Aug 7 2020
## Goal is to rename output of angsd and upload single files to OwnCloud then Dryad. 

dir=~/lfs2/devils/contemporary_sel/angsd_2019-01-18/results
cd $dir

find . -type f -exec sh -c 'for f do x=${f#./}; y="${x// /_}"; echo "cp ${x// /\ } for_dryad/${y////_}"; done' {} +

find . -type f -exec sh -c 'for f do x=${f#./}; y="${x// /_}"; eval "cp ${x// /\ } for_dryad/${y////_}"; done' {} +

