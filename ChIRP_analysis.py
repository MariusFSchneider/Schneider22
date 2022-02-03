
import subprocess
import os

# function for automatic alginment, sorting, removing duplicates and generationg of bdg.files

def Alignment(file_name):
''' Function to align and filter sequencing reads and to generated bedgraph file. Ref genome = mm10, aligner = bowtie 2.0
Values:
	file_name(str): file_name of file
Returns:
	files: aligned and filtered sequencing files
Raises:
	Nothing
''' 
	fastq = file_name + ".fastq"
	fastq1 = file_name + "_fw.fastq"
	fastq2 = file_name + "_rv.fastq"
	sam = file_name + ".sam"
	bam = file_name + ".bam"
	bam_bq = file_name + "_bq.bam"
	bam_sorted = file_name + "_sorted.bam"
	bam_rmdup = file_name + "_rmdup.bam"
	bdg_file = file_name + ".bdg"
	## check if single end or paired end
	if os.path.exists(fastq):
	# process file as single end
		print subprocess.call(("bowtie2 -q -x mm10 -U " + fastq + " -S " +sam + " -p 4"), shell=True)
	else:
	## process files as paired end
		print subprocess.call(("bowtie2 -q -x mm10 -1 " + fastq1 + " -2 " + fastq2 + " -S " +sam + " -p 4"), shell=True)
	print subprocess.call(("samtools view -S -b " + sam + " > "  + bam ), shell=True)
	print subprocess.call(("samtools view -bq " + bam + " > "  + bam_bq ), shell=True)
	print subprocess.call(("samtools sort " + bam_bq + " " + bam_sorted), shell=True)
	print subprocess.call(("samtools rmdup " + bam_sorted + " " + bam_rmdup), shell=True)
	print subprocess.call(("bamCoverage  -b " + bam_sorted + " -o " + bdg_file + " --numberOfProcessors 4 --binSize 10 --normalizeUsing CPM --outFileFormat bedgraph"), shell=True)

# function to call ChIP peaks
	
def CallPeaks(threat, input):

''' Function to call peaks using MACS14
Values:
	file_name(threat): file_name of threatment file
	file_name(input): file_name of control file
Returns:
	files: abed files: peaks and summits
Raises:
	Nothing
''' 

	threat_file = threat + "_rmdup.bam"
	input_file = input +"_rmdup.bam"
	name_output = threat 
	print subprocess.call(("macs14 -t" + threat_file + " -c " + input_file + " -g mm -n " +name_output), shell=True)


#loop through triplicates	
#call alignment function and cakk peaks function for input, even and odd sample
for i in range(3):

	input = "input_rep" + str(i+1)
	odd = "odd_rep"  + str(i+1)
	even = "even_rep" + str(i+1)
	twoCol = "twoCol_rep" + str(i+1)
	oneCol = "RUS_rep" + str(i+1) 
	
	Alignment(input)
	Alignment(odd)
	Alignment(even)

	CallPeaks(odd, input)
	CallPeaks(even, input)




#combine bedfiles ad call lowest peaks , no noramlization due to artifical peaks and loss of peaks 
	print subprocess.call(("bedtools unionbedg -i "+ even +"_rv.bdg "+ odd+ "_rv.bdg > "+ twoCol + ".bdg"), shell=True)
	print subprocess.call(("python2 takeLower.py "+ twoCol + ".bdg "+ oneCol+ ".bdg"), shell=True)


### find overlaps and remove genomic interval probes using R

print subprocess.call(("Rscript findOverlaps.s > peaks_report.txt"), shell=True)