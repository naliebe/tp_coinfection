library(tidyverse)
library(phylotools)
library(reshape2)

`%!in%` <- negate(`%in%`)

path <- "/Users/Nicole/Desktop/U002/TP0865_all_fastq/FINAL_renamed/merged/fastas/"
fastas <- list.files("/Users/Nicole/Desktop/U002/TP0865_all_fastq/FINAL_renamed/merged/fastas", pattern=".fasta")


df <- as.data.frame(matrix(nrow=0, ncol=4))
colnames(df) <- c("seq.name", "seq.text", "n", "merged_reads")

for (i in fastas) {  
  print(match(i, fastas))
  f <- phylotools::read.fasta(paste0(path, i))
  g <- f %>% group_by(seq.text) %>% summarize(n=n()) %>% arrange(desc(n))
  g$sample <- i
  g$merged_reads <- sum(g$n)
  df <- rbind(df, g)
}

write_csv(df, paste0(path, "FINAL_merged_reads.csv"))






#df <- read_csv("/Users/Nicole/Desktop/nicole_analysis_2x300/merged_reads.csv")

#remove anything that has Ns: 
df2 <- df %>% filter(!str_detect(seq.text, "N"))

#remove anything without <0.1% abundance (arbitrary, primarily to make calculations easier, rounding error)
df2$perc <- df2$n*100 / df2$merged_reads
df3 <- df2 %>% filter(perc > 0.1)

unique_vars <- df3 %>% group_by(seq.text) %>% summarize(n_samples=n())
unique_vars$name <- paste0("vars", str_pad(rownames(unique_vars), 4, side="left", pad=0))

df4 <- left_join(unique_vars, df3)


delim2 <- as.data.frame(str_split_fixed(df4$sample, "-", 4))
#ctrl <- df4 %>% filter(str_detect(sample, "CS") | str_detect(sample, "LS"))
#biop <- df4 %>% filter(sample %!in% ctrl$sample)

#biop_maxstrict_P1 <- biop %>% filter(str_detect(sample, "P1")) %>% filter(merge_strictness == "maxstrict") %>% filter(perc >1)
#biop_maxstrict_P2 <- biop %>% filter(str_detect(sample, "P2")) %>% filter(merge_strictness == "maxstrict") %>% filter(perc >1)

loose <- df4 %>% filter(perc >1)
#ctrl_loose <- ctrl %>% filter(merge_strictness == "loose") %>% filter(perc >1)

#biop_ctrl_total <- rbind(biop, ctrl)
#all <- rbind(biop_loose, ctrl_loose)

#all_loose2 <- biop_ctrl_total %>% filter(merge_strictness == "loose") %>% filter(name %in% all$name)

#now rescale to the actual amount of reads in final analsyis:
all_final_reads <- loose %>% group_by(sample) %>% summarize(final_reads=sum(n))
slim <- inner_join(loose, all_final_reads)
slim$final_percent <- slim$n*100 / slim$final_reads

A_slim <- slim %>% filter(str_detect(sample, "ampliconA"))
B_slim <- slim %>% filter(str_detect(sample, "ampliconB"))

#now cast: 
A_cast <- dcast(A_slim, seq.text + name ~ sample, value.var="final_percent")
B_cast <- dcast(B_slim, seq.text + name ~ sample, value.var="final_percent")


#exporting and doing the remaining calculations in excel/geneious.... 
write_csv(A_cast, paste0(path, "A_cast.csv"))
write_csv(B_cast, paste0(path, "B_cast.csv"))
write_csv(slim,paste0(path, "for_recombination_calculations.csv"))

write_csv(A_slim, paste0(path, "A_slim.csv"))
write_csv(B_slim, paste0(path, "B_slim.csv"))


#Now need to determine Number of raw and merged reads per sample - 

raw <- read_tsv("/Users/Nicole/Desktop/U002/TP0865_all_fastq/FINAL_renamed/x2reads.txt")
names <- read_tsv("/Users/Nicole/Desktop/U002/TP0865_all_fastq/FINAL_renamed/names.txt", col_names = FALSE)

reads <- cbind(names, raw)
colnames(reads) <- c("sample", "raw_reads")

reads$raw_reads <- reads$raw_reads / 4
reads$sample <- gsub("_R1.fastq.gz", "", reads$sample)

merged_reads <- slim %>% select(sample, merged_reads)
merged_reads <- merged_reads[!duplicated(merged_reads$sample),]
merged_reads$sample <- gsub("_merged_loose.fasta", "", merged_reads$sample)
reads <- left_join(reads, merged_reads)
reads$percent_merged <- reads$merged_reads*100 / reads$raw_reads

write_csv(reads, paste0(path, "merged_reads_FINAL.csv"))





