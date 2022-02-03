library(GenomicRanges)

setwd('C:/Users/mariu/Documents/newChRIP_analysis')
files_OI <- c( "even", "odd")

reps_OI_p <- c("_rep1_peakss.bed", "_rep2_peaks.bed", "_rep3_peaks.bed")


reps_OI_s <- c("_rep1_summits.bed", "_rep2_summits.bed", "_rep3_summits.bed")





fun_find_overlaps_summit <- function(files_list, rep_name){

file_even <- paste(files_list[1], rep_name, sep = "",collapse = NULL)
file_odd <- paste(files_list[2], rep_name, sep = "",collapse = NULL)
data_even <- read.delim(file_even, sep = "\t", header = FALSE)
data_odd <- read.delim(file_odd, sep = "\t", header = FALSE)
mean_odd <- mean(data_odd[,5])
data_odd <- data_odd[data_odd[,5]>= 0.5*mean_odd,]
mean_even <- mean(data_even[,5])
data_even <- data_even[data_even[,5]>= 0.5*mean_even,]
data_even <- data_even[1:3]
colnames(data_even) <- c("chr", "start", "end")
gr_even <- GRanges(data_even)
gr_even <- resize(gr_even, 250, fix = "center")
data_odd <- data_odd[1:3]
colnames(data_odd) <- c("chr", "start", "end")
gr_odd <- GRanges(data_odd)
gr_odd <- resize(gr_odd, 250, fix = "center")
gr_overlap1 <- subsetByOverlaps(gr_even, gr_odd)
gr_overlap2 <- subsetByOverlaps(gr_odd, gr_even)
gr_overlap <- c(gr_overlap1, gr_overlap2)
gr_overlap <- reduce(gr_overlap)
return(gr_overlap)
}




fun_find_overlaps_peaks <- function(files_list, rep_name){

file_even <- paste(files_list[1], rep_name, sep = "",collapse = NULL)
file_odd <- paste(files_list[2], rep_name, sep = "",collapse = NULL)
data_even <- read.delim(file_even, sep = "\t", header = FALSE)
data_odd <- read.delim(file_odd, sep = "\t", header = FALSE)
data_even <- data_even[1:3]
colnames(data_even) <- c("chr", "start", "end")
gr_even <- GRanges(data_even)
data_odd <- data_odd[1:3]
colnames(data_odd) <- c("chr", "start", "end")
gr_odd <- GRanges(data_odd)
gr_overlap <- subsetByOverlaps(gr_even, gr_odd)
return(gr_overlap)
}


### call every overlapping summits as background peaks file
gr_merge_summits <- fun_find_overlaps_summit(files_OI, reps_OI_s[1])
i = 2
while (i <= length(reps_OI_s)){
	gr_rep <- fun_find_overlaps_summit(files_OI, reps_OI_s[i])

	gr_merge_summits <-c(gr_merge_summits, gr_rep)
	i = i +1 

}
gr_merge_summits <- reduce(gr_merge_summits)


#### call overlapping peaks between all 3 data sets

gr_merge_peaks <- fun_find_overlaps_peaks(files_OI, reps_OI_p[1])
gr_merge_peaks <- subsetByOverlaps(gr_merge_peaks, gr_merge_summits)
print(c("peak in ",reps_OI_p[1]))
print(gr_merge_peaks)
print(c("--------------"))
i = 2
while (i <= length(reps_OI_p)){
	gr_rep <- fun_find_overlaps_peaks(files_OI, reps_OI_p[i])
	gr_rep <- subsetByOverlaps(gr_rep, gr_merge_summits)
	print(c("peak in ",reps_OI_p[i]))
	print(gr_rep)
	print(c("--------------"))
	gr_merge_peaks <- subsetByOverlaps(gr_merge_peaks, gr_rep)
	print(gr_merge_peaks)
	print(c("--------------"))
	i = i +1 

}

print("final peaks list):")

print(gr_merge_peaks)

dataDF <- data.frame(gr_merge_peaks)

write.table(dataDF, "overlapping_peaks.csv", sep = ";", quote = FALSE, row.names = FALSE)
write.table(dataDF, "overlapping_peaks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names= FALSE)


###remove probe intervalls


dataGR <- gr_merge_peaks
dataGR <-reduce(dataGR)


probes <-  read.delim('probes.bed', sep = '\t', header = FALSE)

colnames(probes) <- c("chr","start","end")

PROBES <- GRanges(probes)

PROBES <- reduce(PROBES)
print(c("input peaksremove interval"))
print(PROBES)
print(c("-----------------"))

hits_probes <- findOverlaps(dataGR, PROBES)
hits_probes_df <- data.frame(hits_probes)
remove_intervals <- hits_probes_df$queryHits



dataDF <- data.frame(dataGR)
kept_intervals <- seq(1,length(dataDF[,1]))
kept_intervals <- kept_intervals[! kept_intervals %in% remove_intervals]
dataDF <- dataDF[kept_intervals,]
selected_GR <- GRanges(dataDF)

print(c("-----------------"))
print("final_peak set")
print(selected_GR)

write.table(dataDF, "ol2_peaks_probes_intervals_removed.csv", sep = ";", quote = FALSE, row.names = FALSE)
write.table(dataDF, "ol2_peaks_probes_intervals_removed.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names= FALSE)




