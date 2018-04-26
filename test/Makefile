
WTFQ = ../wtfq

all: orig.skesa.fa wtfq.skesa.fa
	fa -e $^
#	seqtk fqchk $(word 1,$^) | head -n 3
#	seqtk fqchk $(word 2,$^) | head -n 3

orig.skesa.fa: R1.fq.gz R2.fq.gz
	skesa --gz --fastq R1.fq.gz,R2.fq.gz --use_paired_end > $@
	
wtfq.skesa.fa: D1.fq.gz D2.fq.gz
	skesa --gz --fastq R1.fq.gz,R2.fq.gz --use_paired_end > $@

PE.fq.gz: R1.fq.gz R2.fq.gz
	seqtk mergepe $^ | seqtk trimfq - | $(WTFQ) - | seqtk dropse - | pigz > $@

D1.fq.gz: PE.fq.gz
	seqtk seq -1 $< | pigz > $@

D2.fq.gz: PE.fq.gz
	seqtk seq -2 $< | pigz > $@
	