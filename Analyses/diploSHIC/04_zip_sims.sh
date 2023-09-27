cd sims

ls *.txt | parallel gzip {}


#for i in $( ls phleg*.txt ); do gzip $i; done