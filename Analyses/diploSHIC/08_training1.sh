
cd /home/gthom/Dropbox/1_Chapman_fellowship/data/DiploShic
conda deactivate
conda activate selection3

# run several training sessions (e.g., ~5) and choose the best training model for each species

for i in $(seq 1 5); do python diploSHIC/diploSHIC.py train training_phleg/ training_phleg/ phlegModel$i > phlegModel_run$i; done

for i in $(seq 1 5); do python diploSHIC/diploSHIC.py train training_xipho/ training_xipho/ xiphoModel$i > xiphoModel_run$i; done

for i in $(seq 1 5); do python diploSHIC/diploSHIC.py train training_lipau/ training_lipau/ lipauModel$i > lipauModel_run$i; done
