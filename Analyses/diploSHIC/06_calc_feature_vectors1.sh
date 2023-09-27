source activate selection2
cd ../
#cd /home/gthom/Dropbox/1_Chapman_fellowship/data/DiploShic

seq 0 10 | parallel -j 11 python diploSHIC/diploSHIC.py fvecSim diploid sims/phleg_soft_{}.txt.gz \
sims/phleg_soft_{}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11

seq 0 10 | parallel -j 11 python diploSHIC/diploSHIC.py fvecSim diploid sims/phleg_hard_{}.txt.gz \
sims/phleg_hard_{}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11

python diploSHIC/diploSHIC.py fvecSim diploid sims/phleg_neutral.txt.gz \
sims/phleg_neutral.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11