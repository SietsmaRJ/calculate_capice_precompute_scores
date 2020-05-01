This folder contains the output of CAPICE (Shuang et al., unpublished).
It is a collection of all pre-computed scores using CADD 1.4 as input.

Input CADD data:
/apps/software/CADD/v1.4-foss-2018b/data/prescored/GRCh37/incl_anno/whole_genome_SNVs_inclAnno.tsv.gz

Input model:
/groups/umcg-gcc/tmp01/umcg-rsietsma/model/xgb_weightedSample_randomsearch.pickle.dat

Output:
/groups/umcg-gcc/tmp01/umcg-rsietsma/output/

Output is:
 - a single logfile (if all goes correctly), logs every x minutes.
 - folders per chromosome (1-22, x, y) containing multiple txt files of pre-computed CAPICE scores.