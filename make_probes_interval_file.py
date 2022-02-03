import subprocess, csv

files = ['probe1','probe2','probe3','probe4','probe5','probe6','probe7','probe8',]
fimo_sites_file = "probes_intervals.txt"

for file in files:
	print subprocess.call(("rna2meme " +file +".fa"), shell=True)
	print subprocess.call(("fimo --o " +file + " " +file +".meme mm10.fa"), shell=True)
	infile = file + "/fimo.txt"
	with open(infile, "rb") as infile, open(fimo_sites_file, "a") as outfile:
		reader = csv.reader(infile)
		next(reader, None)  # skip the headers
		writer = csv.writer(outfile)
		for row in reader:
       # process each row
			writer.writerow(row)
