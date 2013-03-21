SCRIPT_DIR=../Scripts
COUNTS_DIR=counts
DESEQ_DIR=deseq
PLOT_DIR=plot

REP_1=../rep1/counts
REP_2=../rep2/counts
REP_3=../rep3/counts
REP_4=../rep7/counts

# $(REP_1)/Control_HeLa_L1_replicate.txt.gz $(REP_2)/rep1_Control_HeLa_L1.txt.gz $(REP_2)/rep2_Control_HeLa_L1.txt.gz $(REP_3)/rep1_Control_HeLa_L1.txt.gz $(REP_3)/rep2_Control_HeLa_L1.txt.gz $(REP_3)/rep3_Control_HeLa_L1.txt.gz $(REP_3)/rep4_Control_HeLa_L1.txt.gz $(COUNTS_DIR)/.d $(REP_4)rep1_Control_HeLa_L1.txt.gz $(REP_4)rep2_Control_HeLa_L1.txt.gz

concatenate-all-replicates : $(COUNTS_DIR)/.d
	@ gunzip $(REP_1)/Control_HeLa_L1_replicate.txt.gz $(REP_2)/rep1_Control_HeLa_L1.txt.gz $(REP_2)/rep2_Control_HeLa_L1.txt.gz $(REP_3)/rep1_Control_HeLa_L1.txt.gz $(REP_3)/rep2_Control_HeLa_L1.txt.gz $(REP_3)/rep3_Control_HeLa_L1.txt.gz $(REP_3)/rep4_Control_HeLa_L1.txt.gz $(REP_4)/rep1_Control_HeLa_L1.txt.gz $(REP_4)/rep2_Control_HeLa_L1.txt.gz
	@ paste $(REP_1)/Control_HeLa_L1_replicate.txt $(REP_2)/rep1_Control_HeLa_L1.txt $(REP_2)/rep2_Control_HeLa_L1.txt $(REP_3)/rep1_Control_HeLa_L1.txt $(REP_3)/rep2_Control_HeLa_L1.txt $(REP_3)/rep3_Control_HeLa_L1.txt $(REP_3)/rep4_Control_HeLa_L1.txt $(REP_4)/rep1_Control_HeLa_L1.txt $(REP_4)/rep2_Control_HeLa_L1.txt | cut -f1,2,4,6,8,10,12,14,16,18- | sed '1 i\rep1\trep2\trep2\trep3\trep3\trep3\trep3\trep4\trep4' > $(COUNTS_DIR)/ALL_replicates.txt
	@ gzip $(REP_1)/Control_HeLa_L1_replicate.txt $(REP_2)/rep1_Control_HeLa_L1.txt $(REP_2)/rep2_Control_HeLa_L1.txt $(REP_3)/rep1_Control_HeLa_L1.txt $(REP_3)/rep2_Control_HeLa_L1.txt $(REP_3)/rep3_Control_HeLa_L1.txt $(REP_3)/rep4_Control_HeLa_L1.txt $(REP_4)/rep1_Control_HeLa_L1.txt $(REP_4)/rep2_Control_HeLa_L1.txt $(COUNTS_DIR)/ALL_replicates.txt
	
	
deseq-between-all-rep : $(COUNTS_DIR)/ALL_replicates.txt.gz $(DESEQ_DIR)/.d
	@ gunzip $(COUNTS_DIR)/ALL_replicates.txt.gz
	@ Rscript $(SCRIPT_DIR)/compare.all.rep.R $(COUNTS_DIR)/ALL_replicates.txt $(DESEQ_DIR)/
	@ gzip $(COUNTS_DIR)/ALL_replicates.txt $(DESEQ_DIR)/*.csv
	
plot-deseq : $(PLOT_DIR)/.d
	@ gunzip $(DESEQ_DIR)/DESeq_*_ALL.csv.gz
	@ Rscript $(SCRIPT_DIR)/analysis.deseq.rep.R $(DESEQ_DIR)/ $(PLOT_DIR)/
	@ gzip  $(DESEQ_DIR)/DESeq_*_ALL.csv

	
# Create a directory (use DIR/.d)
%/.d:
	@ mkdir -p $(@D)
	@ touch $@

