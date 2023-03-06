######### Jan28,2023
## launch rstudio under rosetta environment at terminal x86-64

###########


library(PureCN)
library(dplyr)



setwd("/Users/weiwu/Desktop/UCSF500_Bam/PureCN_results/PureCN_Rds")

# Set the directory containing the .rds files
#dir <- "/Users/weiwu/Desktop/UCSF500_Bam/PureCN_results/PureCN_Rds"


# Get a list of all .rds files in the directory
files <- list.files(getwd(), pattern = "*.rds", full.names = TRUE)
files

# Initialize an empty list to store the R objects
Fnames <-list()
objects <- list()
results <- list()
nalist <- list()
L2list <-list()
Tolist <-list()
# counter variable

# Use a for loop to read each file and store the R object in the list
for (file in files) {
  Fnames <- str_extract(file, "(?<=PureCN_Rds/).*(?=.rds)")
  objects <- readRDS(file)
 
  df1 <-predictSomatic(objects)
  
  df2 <-df1[!(df1$chr %in% "chrX"| df1$chr %in% "chrY"), ]
  df2$MCN <- (df2$ML.C-df2$ML.M.SEGMENT)
  
 
  save(df2,file = paste("/Users/weiwu/Desktop/UCSF500_Bam/PureCN_results/PureCN_Rds_CSV/",Fnames,"_CNV_Metrix.csv"))
  
  seg <- read.delim(paste(Fnames,"_dnacopy.seg",sep = ""), sep="\t", header = T, stringsAsFactors = F)
  CNVmatrix <- df2
 
  CNVmatrix2 <- CNVmatrix[,c("seg.id", "ML.C", "MCN", "ML.M.SEGMENT", "chr")]
  CNVmatrix2 <- unique(CNVmatrix2)
  check1 <- any(duplicated(CNVmatrix2$seg.id)) ## if TRUE - means there is a duplication - this should flag an error and a check
  if(check1){
    print("ERROR: duplicated seg id; please review this case")
  }
 
  ## set the seg.is for the seg file
  seg$seg.id <- 1:nrow(seg)
  
  ## now we merge the minor CN/MCN with the seg file using the unique seg file
  seg <- merge(seg, CNVmatrix2, by="seg.id", all.x=T, all.y=T)
  table(seg$seg.id) ## ensure no duplicated segIDs
  if(any(table(seg$seg.id)>1)){
    print("ERROR: duplicated seg id; please review this case")
    break
  }
 
  
  ## calculate seg sizes
  seg$segsize <- seg$loc.end-seg$loc.start
  
  ## calculate proportion of segments where we do not have minor/major allele data
  sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100
  
  ## calculate proportion of segments where major allele CN >= 2
  sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100
  
  ## final check to ensure we have total of all seg (NA, <2 MCN & >2MCN)==100
  sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100 + 
    sum(seg$segsize[seg$MCN<2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 + 
    sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 == 100
  
  ## is there GD?
  if((sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100) > 50) {
    isGD <- "TRUE"
  }else{
    isGD <- "FALSE"
  }
 
 results[[Fnames]] <- isGD
 #nalist[[Fnames]] <- sum(is.na(df2$MCN)) / length(df2$MCN)
 nalist[[Fnames]] <- sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100
 
 #L2list[[Fnames]] <- sum(df2$MCN >=2) /length(df2$MCN)
 L2list[[Fnames]] <- sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100
  
 Tolist[[Fnames]] <- sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100 + 
   sum(seg$segsize[seg$MCN<2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 + 
   sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 == 100
}

  results_df <- data.frame(stack(results))
  nalist_df <- data.frame(stack(nalist))
  L2list_df <- data.frame(stack(L2list))
  Tolist_df <- data.frame(stack(Tolist))

  
names(results_df) <-c("WGD","sampleid")
names(nalist_df) <-c("NA%","sampleid")
names(L2list_df) <-c("G2%","sampleid")
names(Tolist_df) <-c("Total_100","sampleid")

Final_result <- merge(results_df,L2list_df,by="sampleid")
Final_result <- merge(Final_result, nalist_df,by="sampleid")
Final_result <- merge(Final_result, Tolist_df,by="sampleid")

write.csv(Final_result,file = "/Users/weiwu/Desktop/UCSF500_Bam/PureCN_results/PureCN_Rds_CSV/UCSF500_EGFRmut_TP53_WGD_results.csv")



############################## script from Maise

#Can you please collate for us in a table for all cases:
  
# % genome with NA MCN
# % genome with MCN >=2
# % genome with MCN <2
  seg <- read.delim("GP4960_4992_dnacopy.seg", sep="\t", header = T, stringsAsFactors = F)
CNVmatrix <- read.delim("GP4960_4992_CNV_Metrix.csv", sep=",", header = T, row.names = 1,stringsAsFactors = F)

# there is a segmentation file, which only contains the total integer CN state for each segment
# CNVmatrix file has all the GL and somatic variants called with the minor and total copy number state
# we can merge these 2 files and pull out segment data to call WGD.

## subset for unique segid rows from the CNVmatrix file
CNVmatrix2 <- CNVmatrix[,c("seg.id", "ML.C", "MCN", "ML.M.SEGMENT", "chr")]
CNVmatrix2 <- unique(CNVmatrix2)
check1 <- any(duplicated(CNVmatrix2$seg.id)) ## if TRUE - means there is a duplication - this should flag an error and a check
if(check1){
  print("ERROR: duplicated seg id; please review this case")
}

## set the seg.is for the seg file
seg$seg.id <- 1:nrow(seg)

## now we merge the minor CN/MCN with the seg file using the unique seg file
seg <- merge(seg, CNVmatrix2, by="seg.id", all.x=T, all.y=T)
table(seg$seg.id) ## ensure no duplicated segIDs
if(any(table(seg$seg.id)>1)){
  print("ERROR: duplicated seg id; please review this case")
}
##now we check that the total CN from both files are consistent
table(seg$C-seg$ML.C, useNA="always")
if(any(!unique(seg$C-seg$ML.C)%in%c(NA,0))){
  print("ERROR: duplicated seg id; please review this case")
}

## calculate seg sizes
seg$segsize <- seg$loc.end-seg$loc.start

## calculate proportion of segments where we do not have minor/major allele data
sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100

## calculate proportion of segments where major allele CN >= 2
sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100

## final check to ensure we have total of all seg (NA, <2 MCN & >2MCN)==100
sum(seg$segsize[is.na(seg$MCN)])/sum(seg$segsize)*100 + 
  sum(seg$segsize[seg$MCN<2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 + 
  sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100 == 100

## is there GD?
if((sum(seg$segsize[seg$MCN>=2 & !is.na(seg$MCN)])/sum(seg$segsize) *100) > 50) {
  isGD <- "TRUE"
}else{
  isGD <- "FALSE"
}
is GD
