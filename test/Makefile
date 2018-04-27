
SKESA := skesa --gz --use_paired_end

.DEFAULT: all
.PHONY: all clean

WTFQ = ../wtfq

all : orig.skesa.fa wtfq.skesa.fa
	fa -e $^
	blastn -dust no -query lmo_0559.fasta -subject $(word 1,$^) | grep -A 5 '> '
	blastn -dust no -query lmo_0559.fasta -subject $(word 2,$^) | grep -A 5 '> '
	
orig.skesa.fa : R1.fq.gz R2.fq.gz
	$(SKESA) --fastq R1.fq.gz,R2.fq.gz > $@
	
wtfq.skesa.fa : R1.fq.gz R2.fq.gz SUPP.fq.gz
	$(SKESA) --fastq R1.fq.gz,R2.fq.gz --fastq SUPP.fq.gz > $@

SUPP.fq.gz : R1.fq.gz R2.fq.gz
	seqtk mergepe $^ | seqtk trimfq - | $(WTFQ) - | pigz > $@
	fq $@

clean:
	rm -fv orig.skesa.fa wtfq.skesa.fa SUPP.fq.gz
