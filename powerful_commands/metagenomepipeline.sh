

##substitute rare gene of choice for intein


####set up psi blast to get pssms
####  *.fa = amino acid sequences of inteins
####   transcont.fasta = aminoacid sequence of metagenome (use emboss transeq on assembled reads)
for b in *.fa; do psiblast -query "$b" -db transcont.fasta -out ${b%.fa}.blast -outfmt 6 -num_iterations 5 -out_pssm ${b%.fa}.pssm -save_pssm_after_last_round; done 
# this takes a while, recommend doin git in  a script
#to get rid of converged sequences/pssm
for file in $(grep -l CONVERGED *.blast); do rm -i ${file%.blast}.pssm; done


###################### search an AA metagenome with all intein pssms ###################


######### 
######### *.pssm= pssm
#########  contigs.fastaa = nucleotide sequence from assembled reads

for b in *.pssm; do tblastn -in_pssm "$b" -db contigs.fasta -out "$b".search -outfmt 6 -evalue 1e-10; done

#use each pssm to blast a metagenome (-d) this command will only work for aa metagenomes
#this must be done in a directory containing all the pssms, an infile for each pssm, and the
#metagenome that the user intends to search which has been previously formated into a db,
#and a list of the names of each pssm in a file called pssm.list

###########extract all the hits from a metagenome using all PSSMs #######################
 ################ remove redundant sequences  hit in multiple PSSMs #####################


cat *.search |cut -f 2 >all.bout
#list all of the hits from a psiblast using all intein PSSMs

sort all.hits >sort.hits
#sort all of the hits alphabetically by contig name from metagenome file

uniq sort.hits >mg_name.hits
#extract only unique hits mg_name.hits is a list of contigs which were found by PSSMs
#can also do in excel, data -> remove duplicates


####need to remove \n in contigs so each sequence is all in one line
#linebreak inbetween multifasta, this is all one line
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' contigs.fasta  >> contigs.eol



##### mg.fas = fasta of metagenome assembly, usually contigs.eol from above
for filn in `cat all.hits`; do grep -A 1 $filn contigs.eol> $filn.seqfile;done
#makes a list of all of the contigs and searches for the corresponding fasta sequence in the metagenome .fas file

#### *.seqfile = fasta sequences of contigs identified by pssms
#### intein.db == multifasta of intein sequences containing all inteins from inteins.com (AA)
for filn in `cat all.hits`; do blastx -query $filn.seqfile -db /home/CAM/yfeng/inteindb/intein.db -outfmt 6 -out $filn.blast ;done
#new blast



for f in `cat all.hits`; do head -1 $f.blast >> all.tsbh; done
#works for allf
#extract tsbh for each contig, identifies hit as being BIL, HNT, HE, mid, or intein

cat *.tsbh > all.tsbh
#make one file of all the tsbh for all contigs 

for filn in *.seqfile; do blastx -query $filn -db /home/CAM/yfeng/inteindb/intein.db -outfmt "6 qaccver saccver qstart qend sstart send evalue qseq" -max_target_seqs 1 >> seqs.tsbh ;done

#sort all.tsbh by bitscore in excel or equivalent  (highlight all -> data -> sort )

# at this point the user needs to determine how many of the hits should be extracts I have 
#determined that bits of 100 or more are actually inteins anything under that cannot be reliably aligned
# to the identified intein hit, once the user has determined how many of the hits are above 100
#use tail to grab _ hits from the bottom of the list for instance if the last 20 hits are 
#above 100: tail -n 20 mg.int.sorted >mg.int.significant




############ NR blast is optional: used to identify gene/species


################ once hits are extracted they need to be blasted against the #############
################ nr db to identify the extein the intein may be inserted into ############

cut -f 1 mg.int.significant > contigs.list
#makes a list of the contigs in the list which have hits over 100
for filn in `cat contigs.list`; do blastall -p blastp -i $filn.seqfile -d nr -m 8 -o $filn.outfile; done
#blasts nr db with significant contigs

for filn in `cat contigs.list`; do grep -m 1 "gi" $filn.outfile > $filn.nr.tsbh;done
#extracts tsbh from nr db

for filn in `cat contigs.list`; do cut -f 2 $filn.nr.tsbh >$filn.nr.gi;done
#extracts gi # from the tsbh of each contig

for filn in *.fa; do blastx -query $filn -db /isg/shared/databases/blast/nr -outfmt 6 -out $filn.outfile -evalue 1e-10 -max_target_seqs 1 -num_threads 10; done

#extract sequences from nr which correspond to the contigs used to search nr 





##########################################################
#########################################################
########################################################
#####################################################


#mapping



#for mapping you need to prepare two sets of nucl sequences
#1) contigs as is: which contains extein/intein/extein
#2) contigs with intein removed: extein/extein
#sometimes the contig seqs might be very large so you have to cut them down to size
#attached is contigtrimmer.pl which does that. usage: perl contigtrimmer.pl all.tsbh
	#you may need to adjust or look at script
	#if you trim contigs then you need to blastx again to find location of the intein in the sequence
#to artificially remove inteins from the exteins to get the second set of sequences use extein.pl 
	#similar to contigtrimmer.pl, except it takes theflanking region instead of the intein
	#usage: perl extein.pl all.tsbh
	#again may need to adjust script
	
############After you have both sets of sequences you can map them to the reads#####


###do this for both sets of sequences,  make sure you change the names of output and input files####

##i would do all of this with a submission script unless the reads file is small####

###build bowtie index
bowtie2-build extonly.txt extonly
#extonly.txt = multifasta with second set of sequences from above
#extonly name of bowtie index (can be anything)



#map reads back to sequences
bowtie2 -x exind -1 /home/CAM/yfeng/metagenomes/lakevida/trimmedf.fastq -2 /home/CAM/yfeng/metagenomes/lakevida/trimmr.fastq -S exmap.sam -p 10 
#  -x: base name of the bowtie index
# 	-1 and -2 path to forward and reveerse reads, can also work for unpaired reads
#	-S sam file to output to
#	-p number of threads


#makes your sam file much smaller so its downloadable, and also convert sam to bam
samtools view -b -F 4 combmap.sam > combmapped.bam
#only thing you should chabge is the .sam and .bam file name
#these other parameters are required for filtering to work

#sort bam file by genome position so it actually makes sense
samtools sort exmapped.bam -o exsort.bam
#change .bam to output file name of previous step

#you can then load the files into IGV and see the mapped reads and coverage or use them in bedtools
#for Integrated genome viewer(IGV) you need to generate a bam index for your map file and a fasta index for your sequences


#bam index
samtools index -b name.bam

#fasta index
samtools faidx name.fasta 




