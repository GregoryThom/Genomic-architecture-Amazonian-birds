source activate selection2

cd /home/gthom/Dropbox/1_Chapman_fellowship/data/DiploShic
mkdir rawFVFiles
mv sims/*.fvec rawFVFiles

# apparently the maketrainingsets command cannot have extra underscores
rename 's/phleg_/phleg/g' rawFVFiles/*
rename 's/xipho_/xipho/g' rawFVFiles/*
rename 's/lipau_/lipau/g' rawFVFiles/*

mkdir training_phleg
mkdir training_xipho
mkdir training_lipau


python diploSHIC/diploSHIC.py makeTrainingSets rawFVFiles/phlegneutral.fvec \
rawFVFiles/phlegsoft rawFVFiles/phleghard 5 0,1,2,3,4,6,7,8,9,10 training_phleg

python diploSHIC/diploSHIC.py makeTrainingSets rawFVFiles/xiphoneutral.fvec \
rawFVFiles/xiphosoft rawFVFiles/xiphohard 5 0,1,2,3,4,6,7,8,9,10 training_xipho

python diploSHIC/diploSHIC.py makeTrainingSets rawFVFiles/lipauneutral.fvec \
rawFVFiles/lipausoft rawFVFiles/lipauhard 5 0,1,2,3,4,6,7,8,9,10 training_lipau