library(ggplot2)

path = "/mnt/project1/projects/POLIO/DDNS_protocol/DDNS_on_going_val/batch1_1_jan/PoP_output/reports/"
qc = read_excel(paste(path,"reports.xlsx",sep = ""),sheet = "General_QC")
s1_qc = read_excel(paste(path,"reports.xlsx",sep = ""),sheet = "Sabin1")
s2_qc = read_excel(paste(path,"reports.xlsx",sep = ""),sheet = "Sabin2")
s3_qc = read_excel(paste(path,"reports.xlsx",sep = ""),sheet = "Sabin3")


#### % mapped reads for each reference ######

qc[c('ddns_num', 'sample','dup')] <- str_split_fixed(qc$Sample, '_', 3)
qc$dup[grepl("POOL", qc$ddns_num, fixed = TRUE)] <- "POOL"
qc <- qc[order(qc$sample,decreasing=TRUE),]

melt_data <- melt(qc[c("sample","dup","%Sabin1 (from mapped_reads)","%Sabin2 (from mapped_reads)","%Sabin3 (from mapped_reads)")], id = c("sample","dup")) 
ggplot(melt_data, aes(x = reorder(factor(dup), value,sum ), y = value, fill = factor(variable),)) +
  geom_bar(stat="identity", width = 0.9) +
  labs(x = "Sample", y = "Mapped Reads", fill = "Reference", title = "%Mapped Reads to Reference by Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(. ~ sample) 


###### coverage cns5 ######
coverage = data.frame(s1_qc$Sample, s1_qc$coverage_cns5, s2_qc$coverage_cns5, s3_qc$coverage_cns5)
colnames(coverage) <- c("Sample" ,"Sabin1", "Sabin2", "Sabin3")

coverage[c('ddns_num', 'sample','dup')] <- str_split_fixed(coverage$Sample, '_', 3)
coverage$dup[grepl("POOL", coverage$ddns_num, fixed = TRUE)] <- "POOL"
coverage <- coverage[order(coverage$sample,decreasing=TRUE),]

melt_data <- melt(coverage[c("sample","dup","Sabin1","Sabin2","Sabin3")], id = c("sample","dup")) 
ggplot(melt_data, aes(x = reorder(factor(dup), value,sum ), y = value, fill = factor(dup),)) +
  geom_bar(stat="identity", width = 0.9) +
  labs(x = "Sample", y = "Coverage", fill = "Replicate", title = "%Coverage (CNS_5) to Reference by Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(variable ~ sample) 


#####mutations####