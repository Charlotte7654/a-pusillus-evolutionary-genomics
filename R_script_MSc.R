#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory")

library(BiocManager)
library(devtools)

# House admin ####

#Restar
.rs.restartR()

#Install 
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Generate citations
  #For R
  citation()
  
  #For R Studio
  library(rstudioapi)
  versionInfo()$citation
  citation("pkgname")

#Update R in RStudio (didn't work)
devtools::install_github("AndreaCirilloAC/updateR")
library(updateR)
updateR(admin_password = "Kikapoo98") #Didn't work
updateR() #Didn't work, I think because I need the intell processor version. 

#Update R manually 
https://cran.r-project.org/bin/macosx/


# Missing data 96 (DONE) ####

#package: poppr 
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

#made a copy of populations.snps.genepop from Stacks outputs and named it populations.snps.genepop96.gen"
poppr_genind<-import2genind("/Users/Charlotte/Documents/Stellenbosch_University/masters/Bioinformatics/populations_outputs_96/populations.snps.gen")
View(poppr_genind)
levels(poppr_genind$pop) <- c("K", "PP", "CC", "M", "S", "BL", "LJP", "TS", "SR") #Based on order in popmap_ref
View(poppr_genind)

#Heatmap of missing data 
poppr_geneclone<-as.genclone(poppr_genind)
missing_data_plot_96<-info_table(poppr_geneclone, plot = TRUE) #plot heat map
missing_data_plot_96

#Table of missing data
info_table<-info_table(poppr_geneclone, plot = FALSE)
write.table(info_table, "missing_data_96", quote=FALSE, sep="\t", row.names=TRUE)

#Missing data by individuals
genind_obj <- poppr_genind

# Calculate missing data for each individual
missing_individuals <- apply(genind_obj@tab, 1, function(x) sum(is.na(x)) / length(x))

# Convert to percentages
missing_individuals_percent <- missing_individuals * 100

# Create a data frame with the results
missing_individuals_df <- data.frame(
  Individual = rownames(genind_obj@tab),
  PercentMissing = missing_individuals_percent
)

# View the results
print(missing_individuals_df)
                  Individual PercentMissing
P01_A01_K4         P01_A01_K4      16.623214
P01_A02_K5         P01_A02_K5      13.278915
P01_A03_K6         P01_A03_K6      14.248079
P01_A04_K7         P01_A04_K7      22.644317
P01_A05_K8         P01_A05_K8      13.685690
P01_A06_K9         P01_A06_K9      12.607325
P01_A07_K10       P01_A07_K10      14.612539
P01_A08_K14       P01_A08_K14      13.212029
P01_A09_K15       P01_A09_K15      13.673405
P01_A10_K16       P01_A10_K16      15.047980
P01_A11_K17       P01_A11_K17      18.115180
P01_A12_K19       P01_A12_K19      17.601933
P01_B01_K45       P01_B01_K45      14.482862
P01_B02_K46       P01_B02_K46      12.573199
P01_B03_K48       P01_B03_K48      14.486957
P01_B04_PP1       P01_B04_PP1      15.770076
P01_B05_PP2       P01_B05_PP2      15.614464
P01_B06_PP3       P01_B06_PP3      13.674770
P01_B07_PP5       P01_B07_PP5      20.400224
P01_B08_PP6       P01_B08_PP6      13.562839
P01_B09_PP7       P01_B09_PP7      15.239083
P01_B10_PP8       P01_B10_PP8      12.169153
P01_B11_PP10     P01_B11_PP10      17.816241
P01_B12_PP11     P01_B12_PP11      15.441106
P01_C01_PP12     P01_C01_PP12      13.252979
P01_C02_PP13     P01_C02_PP13      15.688175
P01_C03_PP14     P01_C03_PP14      13.935489
P01_C04_PP15     P01_C04_PP15      13.906824
P01_C05_PP16     P01_C05_PP16      13.905459
P01_C06_PP17     P01_C06_PP17      17.214267
P01_C07_PP18     P01_C07_PP18      20.299212
P01_C08_PP21     P01_C08_PP21      14.469212
P01_C09_CC22     P01_C09_CC22      14.234429
P01_C10_CC23     P01_C10_CC23      15.143532
P01_C11_CC24     P01_C11_CC24      14.295854
P01_C12_CC25     P01_C12_CC25      12.821633
P01_D01_CC26     P01_D01_CC26      13.411321
P01_D02_CC28     P01_D02_CC28      14.181193
P01_D03_CC29     P01_D03_CC29      13.606519
P01_D04_CC30     P01_D04_CC30      15.230893
P01_D05_CC31     P01_D05_CC31      14.638474
P01_D06_CC32     P01_D06_CC32      14.820022
P01_D07_CC33     P01_D07_CC33      18.370439
P01_D08_CC34     P01_D08_CC34      15.551673
P01_D09_CC36     P01_D09_CC36      13.403131
P01_D10_CC37     P01_D10_CC37      13.666580
P01_D11_CC38     P01_D11_CC38      16.903043
P01_D12_CC39     P01_D12_CC39      18.961493
P01_E01_M3         P01_E01_M3      14.768151
P01_E02_M5         P01_E02_M5      10.692202
P01_E03_S4         P01_E03_S4      16.646419
P01_E04_S6         P01_E04_S6      15.793281
P01_E05_BLANK2 P01_E05_BLANK2      99.972700
P01_G02_BLANK1 P01_G02_BLANK1      97.807778
P01_E06_LJP1     P01_E06_LJP1      12.509043
P01_E07_LJP2     P01_E07_LJP2      15.349650
P01_E08_LJP8     P01_E08_LJP8      15.761886
P01_E09_LJP9     P01_E09_LJP9      12.964960
P01_E10_LJP12   P01_E10_LJP12      18.227112
P01_E11_LJP13   P01_E11_LJP13      14.456927
P01_E12_LJP14   P01_E12_LJP14      14.485592
P01_F01_LJP15   P01_F01_LJP15      14.065166
P01_F02_LJP16   P01_F02_LJP16      13.910919
P01_F03_LJP17   P01_F03_LJP17      14.102022
P01_F04_LJP18   P01_F04_LJP18      15.073916
P01_F05_LJP19   P01_F05_LJP19      13.831748
P01_F06_LJP20   P01_F06_LJP20      11.239575
P01_F07_LJP26   P01_F07_LJP26      14.526543
P01_F08_LJP28   P01_F08_LJP28       8.669242
P01_F09_LJP33   P01_F09_LJP33      17.206077
P01_F10_LJP35   P01_F10_LJP35      17.375340
P01_F11_LJP37   P01_F11_LJP37      34.409424
P01_F12_TS2       P01_F12_TS2      99.993175
P01_G01_TS3       P01_G01_TS3      14.749041
P01_G03_TS7       P01_G03_TS7      12.547264
P01_G04_TS9       P01_G04_TS9      14.976999
P01_G05_TS12     P01_G05_TS12      14.168908
P01_G06_TS16     P01_G06_TS16      16.352940
P01_G07_TS20     P01_G07_TS20      12.979975
P01_G08_TS29     P01_G08_TS29      15.472502
P01_G09_TS34     P01_G09_TS34      15.764616
P01_G10_TS51     P01_G10_TS51      17.405370
P01_G11_TS52     P01_G11_TS52      19.757299
P01_G12_SR1       P01_G12_SR1      15.465677
P01_H01_SR11     P01_H01_SR11      13.377196
P01_H02_SR16     P01_H02_SR16      11.803328
P01_H03_SR17     P01_H03_SR17      14.215318
P01_H04_SR18     P01_H04_SR18      13.505508
P01_H05_SR25     P01_H05_SR25      18.999713
P01_H06_SR26     P01_H06_SR26      12.083157
P01_H07_SR29     P01_H07_SR29      15.260924
P01_H08_SR32     P01_H08_SR32      13.998280
P01_H09_SR36     P01_H09_SR36      12.207374
P01_H10_SR38     P01_H10_SR38      13.598329
P01_H11_SR39     P01_H11_SR39      99.350250
P01_H12_SR40     P01_H12_SR40      14.380486

#Now missing data per population 
# Calculate missing data for each locus across all individuals within each population
missing_per_population <- apply(genind_obj@tab, 2, function(x) tapply(is.na(x), genind_obj@pop, mean))
View(missing_per_population)

# Convert to percentages
missing_per_population_percent <- missing_per_population * 100

# Calculate average missing data per population
average_missing_population <- colMeans(missing_per_population_percent, na.rm = TRUE)

# Create a data frame with the results
missing_population_df <- data.frame(
  Population = names(average_missing_population),
  PercentMissing = average_missing_population
)

#ok... still for each locus... 
# Calculate the average missing data per population
average_missing_data <- rowMeans(missing_per_population, na.rm = TRUE)

# Convert to a data frame for easier inspection
average_missing_data_df <- data.frame(
  Population = rownames(missing_per_population),
  AverageMissingData = average_missing_data
)

# View the results NB THIS IS PROPORTION 
print(average_missing_data_df)
Population AverageMissingData
K            K          0.1512624
PP          PP          0.1543292
CC          CC          0.1495251
M            M          0.1273018
S            S          0.1621985
BL          BL          0.9889024
LJP        LJP          0.1545362
TS          TS          0.2310619
SR          SR          0.2063427

# missing data 92 (DONE) ####

#package: poppr 
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

#made a copy of populations.snps.genepop from Stacks outputs and named it populations.snps.genepop.gen"
poppr_genind<-import2genind("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.snps.genepop.gen")
View(poppr_genind)
levels(poppr_genind$pop) <- c("K", "PP", "CC", "M", "S", "LJP", "TS", "SR") #Based on order in popmap_ref_92
View(poppr_genind)

#Heatmap of missing data 
poppr_geneclone<-as.genclone(poppr_genind)
missing_data_plot_92<-info_table(poppr_geneclone, plot = TRUE) #plot heat map
missing_data_plot_92

#Table of missing data (optional - I skipped because it takes an hour+ and is huge)
#info_table<-info_table(poppr_geneclone, plot = FALSE)
#write.table(info_table, "missing_data_92", quote=FALSE, sep="\t", row.names=TRUE)

#Missing data by individuals
genind_obj <- poppr_genind

# Calculate missing data for each individual
missing_individuals <- apply(genind_obj@tab, 1, function(x) sum(is.na(x)) / length(x))

# Convert to percentages
missing_individuals_percent <- missing_individuals * 100

# Create a data frame with the results
missing_individuals_df <- data.frame(
  Individual = rownames(genind_obj@tab),
  PercentMissing = missing_individuals_percent
)

# View the results
print(missing_individuals_df)
Individual PercentMissing
P01_A01_K4       P01_A01_K4       17.40322
P01_A02_K5       P01_A02_K5       14.00928
P01_A03_K6       P01_A03_K6       14.87638
P01_A04_K7       P01_A04_K7       23.55139
P01_A05_K8       P01_A05_K8       14.34442
P01_A06_K9       P01_A06_K9       13.27650
P01_A07_K10     P01_A07_K10       15.05459
P01_A08_K14     P01_A08_K14       14.01726
P01_A09_K15     P01_A09_K15       14.36171
P01_A10_K16     P01_A10_K16       15.81397
P01_A11_K17     P01_A11_K17       18.76637
P01_A12_K19     P01_A12_K19       18.38203
P01_B01_K45     P01_B01_K45       15.15833
P01_B02_K46     P01_B02_K46       13.18208
P01_B03_K48     P01_B03_K48       15.11444
P01_B04_PP1     P01_B04_PP1       16.57335
P01_B05_PP2     P01_B05_PP2       16.21029
P01_B06_PP3     P01_B06_PP3       14.29122
P01_B07_PP5     P01_B07_PP5       21.29986
P01_B08_PP6     P01_B08_PP6       14.26197
P01_B09_PP7     P01_B09_PP7       15.77806
P01_B10_PP8     P01_B10_PP8       12.80305
P01_B11_PP10   P01_B11_PP10       18.49241
P01_B12_PP11   P01_B12_PP11       16.02676
P01_C01_PP12   P01_C01_PP12       13.93215
P01_C02_PP13   P01_C02_PP13       16.34860
P01_C03_PP14   P01_C03_PP14       14.56120
P01_C04_PP15   P01_C04_PP15       14.58114
P01_C05_PP16   P01_C05_PP16       14.57051
P01_C06_PP17   P01_C06_PP17       17.97907
P01_C07_PP18   P01_C07_PP18       21.05382
P01_C08_PP21   P01_C08_PP21       15.07986
P01_C09_CC22   P01_C09_CC22       14.77132
P01_C10_CC23   P01_C10_CC23       15.92169
P01_C11_CC24   P01_C11_CC24       14.98544
P01_C12_CC25   P01_C12_CC25       13.54514
P01_D01_CC26   P01_D01_CC26       14.11302
P01_D02_CC28   P01_D02_CC28       14.88569
P01_D03_CC29   P01_D03_CC29       14.33511
P01_D04_CC30   P01_D04_CC30       16.01346
P01_D05_CC31   P01_D05_CC31       15.43495
P01_D06_CC32   P01_D06_CC32       15.47484
P01_D07_CC33   P01_D07_CC33       19.21189
P01_D08_CC34   P01_D08_CC34       16.30205
P01_D09_CC36   P01_D09_CC36       14.05716
P01_D10_CC37   P01_D10_CC37       14.31782
P01_D11_CC38   P01_D11_CC38       17.50296
P01_D12_CC39   P01_D12_CC39       19.64944
P01_E01_M3       P01_E01_M3       15.43628
P01_E02_M5       P01_E02_M5       11.38803
P01_E03_S4       P01_E03_S4       17.33805
P01_E04_S6       P01_E04_S6       16.54675
P01_E06_LJP1   P01_E06_LJP1       13.24458
P01_E07_LJP2   P01_E07_LJP2       16.00681
P01_E08_LJP8   P01_E08_LJP8       16.45366
P01_E09_LJP9   P01_E09_LJP9       13.71803
P01_E10_LJP12 P01_E10_LJP12       19.03236
P01_E11_LJP13 P01_E11_LJP13       15.20088
P01_E12_LJP14 P01_E12_LJP14       15.21418
P01_F01_LJP15 P01_F01_LJP15       14.78462
P01_F02_LJP16 P01_F02_LJP16       14.65695
P01_F03_LJP17 P01_F03_LJP17       14.90830
P01_F04_LJP18 P01_F04_LJP18       15.83791
P01_F05_LJP19 P01_F05_LJP19       14.64897
P01_F06_LJP20 P01_F06_LJP20       11.84951
P01_F07_LJP26 P01_F07_LJP26       15.31792
P01_F08_LJP28 P01_F08_LJP28        9.26948
P01_F09_LJP33 P01_F09_LJP33       18.02295
P01_F10_LJP35 P01_F10_LJP35       18.19052
P01_F11_LJP37 P01_F11_LJP37       35.36100
P01_G01_TS3     P01_G01_TS3       15.52937
P01_G03_TS7     P01_G03_TS7       13.38156
P01_G04_TS9     P01_G04_TS9       15.84190
P01_G05_TS12   P01_G05_TS12       14.87107
P01_G06_TS16   P01_G06_TS16       17.10798
P01_G07_TS20   P01_G07_TS20       13.69276
P01_G08_TS29   P01_G08_TS29       16.09990
P01_G09_TS34   P01_G09_TS34       16.56404
P01_G10_TS51   P01_G10_TS51       18.15728
P01_G11_TS52   P01_G11_TS52       20.89157
P01_G12_SR1     P01_G12_SR1       16.12251
P01_H01_SR11   P01_H01_SR11       14.15424
P01_H02_SR16   P01_H02_SR16       12.45462
P01_H03_SR17   P01_H03_SR17       14.93357
P01_H04_SR18   P01_H04_SR18       14.20478
P01_H05_SR25   P01_H05_SR25       19.78775
P01_H06_SR26   P01_H06_SR26       12.74853
P01_H07_SR29   P01_H07_SR29       16.00814
P01_H08_SR32   P01_H08_SR32       14.62902
P01_H09_SR36   P01_H09_SR36       12.92940
P01_H10_SR38   P01_H10_SR38       14.23271
P01_H12_SR40   P01_H12_SR40       15.01603

#Now missing data per population 
# Calculate missing data for each locus across all individuals within each population
missing_per_population <- apply(genind_obj@tab, 2, function(x) tapply(is.na(x), genind_obj@pop, mean))
View(missing_per_population)

# Convert to percentages
missing_per_population_percent <- missing_per_population * 100

# Calculate average missing data per population
average_missing_population <- colMeans(missing_per_population_percent, na.rm = TRUE)

# Create a data frame with the results
missing_population_df <- data.frame(
  Population = names(average_missing_population),
  PercentMissing = average_missing_population
)

# Calculate the average missing data per population
average_missing_data <- rowMeans(missing_per_population, na.rm = TRUE)

# Convert to a data frame for easier inspection
average_missing_data_df <- data.frame(
  Population = rownames(missing_per_population),
  AverageMissingData = average_missing_data
)

# View the results NB THIS IS PROPORTION 
print(average_missing_data_df)
Population AverageMissingData
P01_B03_K48     P01_B03_K48          0.1582080
P01_C08_PP21   P01_C08_PP21          0.1610843
P01_D12_CC39   P01_D12_CC39          0.1565762
P01_E02_M5       P01_E02_M5          0.1341215
P01_E04_S6       P01_E04_S6          0.1694240
P01_F11_LJP37 P01_F11_LJP37          0.1620659
P01_G11_TS52   P01_G11_TS52          0.1621374
P01_H12_SR40   P01_H12_SR40          0.1476844

# BayeScan () ####

#STEP 1. 
#In terminal. 
#Convert genepop file from pipeline to bayescan file using PGD spider (code in Notion)

#STEP 2. 
#In hpc. 
#Upload to hpc and run bayescan there (code in Notion)

#STEP 3. 
#In terminal. 
#Generate list of SNPs in the right order (code in Notion)

#STEP 4. 
#in R.
#Export the R _fst2.txt file to workspace and load below for analysis.
bayescan=read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/bayescan/_fst.txt")
View(bayescan)

#STEP 5:
#In R: 
#Then import the SNPs information.
snp=read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/bayescan/populations.snps.id.txt",header=FALSE) 
#75193 snps
View(snp)
View(bayescan)

#Merge the names of the outliers with the results from the bayescan dataframe.
snp_bayescan=cbind(snp, bayescan)

#Rename columns.
colnames(snp_bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")

#Write the results.
write.table(snp_bayescan, "BFS-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE)

View(snp_bayescan)

#Change the value of the Q_VALUE column: 0 == 0.0001.
class(snp_bayescan$Q_VALUE) ## [1] "numeric"

snp_bayescan$Q_VALUE <- as.numeric(snp_bayescan$Q_VALUE) 
snp_bayescan[snp_bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 
View(snp_bayescan)

#Round the values.
snp_bayescan$LOG_PO <- (round(snp_bayescan$LOG_PO, 4)) 
snp_bayescan$Q_VALUE <- (round(snp_bayescan$Q_VALUE, 4)) 
snp_bayescan$ALPHA <- (round(snp_bayescan$ALPHA, 4)) 
snp_bayescan$FST <- (round(snp_bayescan$FST, 6))

#Add a column for the type of selection grouping based on a Q-VALUE < 0.05. 
#You can also choose a Q-VALUE < 0.01 if you want to be more conservative.
snp_bayescan$SELECTION <- ifelse(snp_bayescan$ALPHA>=0&snp_bayescan$Q_VALUE<=0.05,"diversifying",
                                 ifelse(snp_bayescan$ALPHA>=0&snp_bayescan$Q_VALUE>0.05,"neutral","balancing")) 

snp_bayescan$SELECTION<- factor(snp_bayescan$SELECTION)
levels(snp_bayescan$SELECTION) 
View(snp_bayescan)

#Save the results of the SNPs potentially under positive (divergent) 
#and balancing selection (qvalue < 0.05).
positive <- snp_bayescan[snp_bayescan$SELECTION=="diversifying",] 
neutral <- snp_bayescan[snp_bayescan$SELECTION=="neutral",] 
balancing <- snp_bayescan[snp_bayescan$SELECTION=="balancing",]

#Check the number of SNPs belonging to each category.
xtabs(data=snp_bayescan, ~SELECTION) 
#SELECTION
#balancing diversifying      neutral 
#38302         2344        34547 

#Write the results of the SNPs potentially under selection (qvalue < 0.05).
write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 

#Transformation Log of the Q value in order to create the ggplot graph.
range(snp_bayescan$Q_VALUE)  ### 0.0001 0.8355
snp_bayescan$LOG10_Q <- -log10(snp_bayescan$Q_VALUE)

#Make the ggplot graph.
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
library(ggplot2)

#Create a title for the ggplot graph.
#x_title="Log(q-value)" 
#y_title="Fst" 

ggplot(snp_bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)", y="Fst")+
#  labs(x_title)+
#  labs(y_title)+
  theme_classic()

#balancing diversifying      neutral 
#38302         2344        34547 

#balancing and diversifying outliers 
bayescan_outliers_div_bal <- c(positive$SNP, balancing$SNP)
length(bayescan_outliers_div_bal) #40646

write.table(bayescan_outliers_div_bal, "bayescan_outliers_div_bal.txt", row.names=F, quote=F) 

#just diversifying outliers 
bayescan_outliers_div <- c(positive$SNP)
length(bayescan_outliers_div) #40646

write.table(bayescan_outliers_div, "bayescan_outliers_divs.txt", row.names=F, quote=F) 

# snmf () ####

#STEP 1
#In terminal
#grep -v "#" populations.snp.resubmission72.vcf | cut -f 3 > populations.snp.resubmission72.snp.id.txt
#and upload that to hpc for later use 

#STEP 2
#in R on laptop: change vcf to geno for use in LEA 

if (!requireNamespace("LEA", quietly = TRUE))
  BiocManager::install("LEA")
library(LEA)

vcf2geno("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.snps.vcf")
#sends to populations_outputs_92 but it keeps aborting when I try to place output in R_directory so gonna move on. Just move it manually. 
#and upload that to hpc

#STEP 3
#Create a snmf_full project and associated K plot in hpc to determine optimal K (code in Notion) 

#STEP 4
#Detect outliers: Use K from step 2 to plug into pvalue function to detect outliers in snmf_full project. 

#STEP 5 
#Download outlier list from hpc and upload to R. 
outliers_snmf <- read.delim("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/snmf/outliers_snmf.txt", na.strings="")
snmf_outliers <- outliers_snmf[ ,"LOCUS"]
length(snmf_outliers) #3760

### MAYBE AFTER LOOKING AT ALL OUTLIER? Step 4: removed outliers to create netural and outlier .geno files? AND repeat K analysis for 
#Do it in hpc qsubi to avoid downloading all of project_snmf_full data. 

# PCadapt () ####

#STEP 1
#in terminal 
#run PLINK: ped to bed
#populations.plink.ped is an output of STACKS populations program
#outputs files to R directory. I put them in new folder PCadapt_outputs but then moved them back.  

system ("/Users/Charlotte/Documents/Bioinformatics/Programs/plink_mac_developer.nosync/plink --file /Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.plink --make-bed --allow-extra-chr --out /Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/pcadapt/populations.bedfile")
#need to keep on same line 
#75193 variants, 92 people

#STEP 2 
#in R
if (!require("pcadapt", quietly = TRUE))
  install.packages("pcadapt")
library(pcadapt)

#read bed: 
bedfile <- system.file("extdata", "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/pcadapt/populations.bedfile.bed", package = "pcadapt")

mat <- bed2matrix("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/pcadapt/populations.bedfile.bed") 
dim(mat) #92 75193

read.pcadapt.output <- read.pcadapt(
  mat,
  type = c("lfmm"),
)

pc.adapt.output.findK <- pcadapt(
  read.pcadapt.output,
  K = 15,
  method = "mahalanobis",
  min.maf = 0.05,
  ploidy = 2,
  LD.clumping = NULL,
  pca.only = FALSE,
  tol = 1e-04
)

popmap <- read.delim("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/popmap_ref_92", header=FALSE, na.strings="")

poplist.int <-popmap [ , "V2"]

plot.manhattan_screeplot <- plot(
  pc.adapt.output.findK,
  option = "screeplot",
  i = 1,
  j = 2,
  pop = poplist.int,
  chr.info = NULL,
  snp.info = NULL,
  plt.pkg = "ggplot",
  K = NULL
) #K=3

plot.manhattan <- plot(
  pc.adapt.output.findK,
  option = "manhattan",
  i = 1,
  j = 2,
  pop = poplist.int,
  chr.info = NULL,
  snp.info = NULL,
  plt.pkg = "ggplot",
  K = 3
)

#In terminal (DONE): 
#grep -v "#" populations.snps.vcf | cut -f 3 > populations.snps.vcf.snp.id.txt
#same order as vcf input file 

if (!requireNamespace("qvalue", quietly = TRUE))
  BiocManager::install("qvalue")
library(qvalue)

qval <- qvalue(pc.adapt.output.findK$pvalues)$qvalues
alpha <- 0.05
outliers.pc <- which(qval < alpha)
length(outliers.pc) #7216

snp_list=read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/pcadapt/populations.snps.id.txt", header=FALSE)
snps_qvalues <- cbind(snp_list, qval)
View(snps_qvalues)

write.table(snps_qvalues, "pvalues_92.txt", sep="\t", quote=FALSE)
pvalues_92<-read.table("pvalues_92.txt")

snps_qvalues_no_na <- na.omit(snps_qvalues)
View(snps_qvalues_no_na)

snps_qvalues_0.5 <- snps_qvalues_no_na[snps_qvalues_no_na$qval <= 0.05, ]
pcadapt_outliers<-snps_qvalues_0.5[ ,1]
View(pcadapt_outliers)
length(pcadapt_outliers) #just checking. 7216.

write.table(pcadapt_outliers, "pcadapt_outliers.txt", sep="\t", quote=FALSE)

# Outlier comparison (QUESTION) ####
length(snmf_outliers) #3760
length(bayescan_outliers_div_bal) #40646 37488 37705
length(bayescan_outliers_div) #2344 
length(pcadapt_outliers) #7216

#With balancing loci 
intersect_bayescan_pcadapt_outliers_div_bal<-intersect(bayescan_outliers_div_bal,pcadapt_outliers)
intersect_bayescan_snmf_outliers_div_bal<-intersect(bayescan_outliers_div_bal,snmf_outliers)
intersect_pcadapt_snmf_outliers<-intersect(pcadapt_outliers,snmf_outliers)
intersect123_div_bal<-intersect(intersect_pcadapt_snmf_outliers,bayescan_outliers_div_bal)

length(intersect_bayescan_pcadapt_outliers_div_bal) #1932
length(intersect_bayescan_snmf_outliers_div_bal) #1612
length(intersect_pcadapt_snmf_outliers_div_bal) #890
length(intersect123_div_bal) #709

#Without balancing loci
intersect_bayescan_pcadapt_outliers_div<-intersect(bayescan_outliers_div,pcadapt_outliers)
intersect_bayescan_snmf_outliers_div<-intersect(bayescan_outliers_div,snmf_outliers)
intersect_pcadapt_snmf_outliers<-intersect(pcadapt_outliers,snmf_outliers)
intersect123_div<-intersect(intersect_pcadapt_snmf_outliers,bayescan_outliers_div)

length(intersect_bayescan_pcadapt_outliers_div) #966
length(intersect_bayescan_snmf_outliers_div) #1587
length(intersect_pcadapt_snmf_outliers) #890
length(intersect123_div) #709

#Venn diagram 
if (!requireNamespace("VennDiagram", quietly = TRUE))
  BiocManager::install("VennDiagram")
library(VennDiagram)

# Create Venn Diagram with customized text font and transparency
# Load the required libraries
library(VennDiagram)
library(grid)

#Venn diagram with balancing loci
# Clear the plotting window
grid.newpage()

# Define sets (your outliers)
venn.plot.div.bal <- draw.triple.venn(
  area1 = length(bayescan_outliers_div_bal),
  area2 = length(pcadapt_outliers),
  area3 = length(snmf_outliers),
  n12 = length(intersect_bayescan_pcadapt_outliers_div_bal),
  n23 = length(intersect_pcadapt_snmf_outliers),
  n13 = length(intersect_bayescan_snmf_outliers_div_bal),
  n123 = length(intersect123),
  category = c("Bayescan", "PCAdapt", "SNMF"),
  fill = c("skyblue", "pink1", "mediumorchid"),
  alpha = 0.5,  # Adjust transparency if needed
  label.col = "black",
  cex = 2,  # Adjust text size
  fontfamily = "Arial",  # Use Arial font for the text
  cat.cex = 2,  # Adjust category label size
  cat.fontfamily = "Arial"  # Use Arial font for category labels
)

# Venn diagram with balancing 
grid.draw(venn.plot.div.bal)

venn.plot.div.bal <- venn.diagram(
  x = list(
    SNMF = snmf_outliers,
    Bayescan = bayescan_outliers_div_bal,
    PCAdapt = pcadapt_outliers
  ),
  category.names = c("SNMF", "Bayescan", "PCAdapt"),
  filename = NULL,
  fill = c("red", "green", "blue"), # Customize colors if needed
  alpha = 1, # Adjust transparency (0 = fully transparent, 1 = fully opaque)
  cex = 1.5, # Font size
  fontfamily = "Arial", # Set font to Arial
  cat.cex = 1.2, # Font size for category labels
  cat.fontfamily = "Arial", # Set font to Arial for category labels
  margin = 0.1 # Adjust margins if needed
)

venn.plot.div.bal

# Venn diagram without balancing 
# Clear the plotting window
grid.newpage()

# Define sets (your outliers)
venn.plot.div <- draw.triple.venn(
  area1 = length(bayescan_outliers_div),
  area2 = length(pcadapt_outliers),
  area3 = length(snmf_outliers),
  n12 = length(intersect_bayescan_pcadapt_outliers_div),
  n23 = length(intersect_pcadapt_snmf_outliers),
  n13 = length(intersect_bayescan_snmf_outliers_div),
  n123 = length(intersect123),
  category = c("Bayescan", "PCAdapt", "SNMF"),
  fill = c("skyblue", "pink1", "mediumorchid"),
  alpha = 0.5,  # Adjust transparency if needed
  label.col = "black",
  cex = 2,  # Adjust text size
  fontfamily = "Arial",  # Use Arial font for the text
  cat.cex = 2,  # Adjust category label size
  cat.fontfamily = "Arial"  # Use Arial font for category labels
)

grid.draw(venn.plot.div)

venn.plot.div <- venn.diagram(
  x = list(
    SNMF = snmf_outliers,
    Bayescan = bayescan_outliers_div,
    PCAdapt = pcadapt_outliers
  ),
  category.names = c("SNMF", "Bayescan", "PCAdapt"),
  filename = NULL,
  fill = c("red", "green", "blue"), # Customize colors if needed
  alpha = 1, # Adjust transparency (0 = fully transparent, 1 = fully opaque)
  cex = 1.5, # Font size
  fontfamily = "Arial", # Set font to Arial
  cat.cex = 1.2, # Font size for category labels
  cat.fontfamily = "Arial", # Set font to Arial for category labels
  margin = 0.1 # Adjust margins if needed
)

venn.plot.div

#With balancing and diversifying loci
# Step 1: Combine outliers detected by at least two methods
shared_by_two_or_more_div_bal <- union(intersect_bayescan_pcadapt_outliers_div_bal, 
                               union(intersect_bayescan_snmf_outliers_div_bal, 
                                     intersect_pcadapt_snmf_outliers))


# Step 2: Check the length of the final list of outliers
length(shared_by_two_or_more_div_bal) #3016
outliers_div_bal<-shared_by_two_or_more_div_bal

#With diversifying loci
# Step 1: Combine outliers detected by at least two methods
shared_by_two_or_more_div <- union(intersect_bayescan_pcadapt_outliers_div, 
                                       union(intersect_bayescan_snmf_outliers_div, 
                                             intersect_pcadapt_snmf_outliers))


# Step 2: Check the length of the final list of outliers
length(shared_by_two_or_more_div) #2025
outliers_div<-shared_by_two_or_more_div

write.table(outliers_div_bal, "outliers.div.bal.txt", row.names=F, quote=F)
write.table(outliers_div, "outliers.div.txt", row.names=F, quote=F)

#Create outlier lists
write.table(outliers_div_bal, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci.div.bal.txt", row.names=F, quote=F)
write.table(outliers_div, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci.div.txt", row.names=F, quote=F)

#NB go out of R and manually remove :+, :- and change : to _ in sublime text using replace function. 

#load datasets
outliers.db.table<-read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci.div.bal.txt", header=TRUE)
outliers.d.table<-read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci.div.txt", header=TRUE)

#table to vector
outlier.db <- outliers.db.table[, 1]
outlier.d <- outliers.d.table[, 1]

# load Structure file
populations.structure<-read.structure("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.structure.stru", 
                                      n.ind = 92, n.loc = 75193, col.lab = 1, col.pop = 2, row.marknames = 1)
# Which other optional columns should be read (press 'return' when done)? ENTER
# Are genotypes coded by a single row (y/n)?   n
#if I say y, then Error in mat[, (ncol(mat) - p + 1):ncol(mat)] : only 0's may be mixed with negative subscripts

# create neutral lists
neutral.db <- setdiff(locNames(populations.structure), outlier.db)
length(neutral.db) #72177

neutral.d <- setdiff(locNames(populations.structure), outlier.d)
length(neutral.d) #73168

write.table(neutral.d, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/neutralloci.div.txt", row.names=F, quote=F)
write.table(neutral.db, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/neutralloci.div.bal.txt", row.names=F, quote=F)

# Private alleles (BUSY) ####
#FULL 
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

#poppr_genind called in missing data
poppr_genind <- import2genind("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.snps.genepop.gen")
levels(poppr_genind$pop) <- c("K", "PP", "CC", "M", "S", "LJP", "TS", "SR") #Based on order in popmap_ref_92
poppr_geneclone <- as.genclone(poppr_genind)

# Finding private alleles by population
private_alleles_result <- poppr::private_alleles(poppr_geneclone, level = "population")

# Convert the result to a data frame
private_alleles_df <- as.data.frame(private_alleles_result)
View(private_alleles_df)

#Don't run (print)
# Function to display private alleles by population
display_private_alleles <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    cat("\nPrivate alleles for population:", pop, "\n")
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print the private alleles and their count
    if (length(private_alleles_pop) > 0) {
      print(private_alleles_pop)
      cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
    } else {
      cat("No private alleles for population", pop, "\n")
    }
  }
}

#Don't run (print)
# Assuming private_alleles_result is your input data (e.g., allele frequencies or counts)
# Example of how to call the function
display_private_alleles(private_alleles_result)

#Run (counts)
# Function to display only the count of private alleles by population
display_private_allele_counts <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print only the count of private alleles
    cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
  }
}

#Run (counts)
# Example of how to call the function
display_private_allele_counts(private_alleles_result)

Number of private alleles for population K : 1977 
Number of private alleles for population PP : 2260 
Number of private alleles for population CC : 2352 
Number of private alleles for population M : 280 
Number of private alleles for population S : 239 
Number of private alleles for population LJP : 1286 
Number of private alleles for population TS : 626 
Number of private alleles for population SR : 853 

#Private alleles in SA and AUS
#poppr_genind called in missing data
poppr_genind <- import2genind("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.snps.genepop.gen")
poppr_genind_private<-poppr_genind
levels(poppr_genind_private$pop) <- c("SA", "SA", "SA", "SA", "SA", "AUS", "AUS", "AUS") #Based on order in popmap_ref_92
poppr_geneclone_private <- as.genclone(poppr_genind_private)

# Finding private alleles by population
private_alleles_result <- poppr::private_alleles(poppr_geneclone_private, level = "population")

# Convert the result to a data frame
private_alleles_df <- as.data.frame(private_alleles_result)
View(private_alleles_df)

# Function to display only the count of private alleles by population
display_private_allele_counts <- function(private_alleles_result) {
  
  # Convert the result to a data frame
  private_alleles_df <- as.data.frame(private_alleles_result)
  
  # Iterate through each population
  for (pop in rownames(private_alleles_df)) {
    
    # Select private alleles for the current population
    private_alleles_pop <- private_alleles_df[pop, ]
    
    # Filter out the non-private alleles (e.g., non-zero values)
    private_alleles_pop <- private_alleles_pop[private_alleles_pop != 0]
    
    # Print only the count of private alleles
    cat("Number of private alleles for population", pop, ":", length(private_alleles_pop), "\n")
  }
}

#Run (counts)
# Example of how to call the function
display_private_allele_counts(private_alleles_result)


#identifies alleles that are present in any of the AUS populations 
#but absent from all of the SA populations. 
#These alleles are considered unique to AUS in the context of the given data.







# Subset STRU ####
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#NB need to manually at ".stru" to your structure file from Stacks output and delete line 1 (header)
populations.structure<-read.structure("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosyn/populations.structure.stru") 
#How many genotypes are there? 92
#How many markers are there? 75193
#Which column contains labels for genotypes ('0' if absent)?   1
#Which column contains the population factor ('0' if absent)?   2
#Which other optional columns should be read (press 'return' when done)? ENTER 
#  Which row contains the marker names ('0' if absent)?   1
#Are genotypes coded by a single row (y/n)?   n
#if I say y, then 
#Error in mat[, (ncol(mat) - p + 1):ncol(mat)] : 
#only 0's may be mixed with negative subscripts
View(populations.structure)

# DIV BAL 
#Create outlier and netural lists
write.table(outliers_div_bal, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci_div_bal.txt", row.names=F, quote=F)

#NB Manually remove :+, :- and change : to _ in sublime text using replace function. 
outliers_db_table<-read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci_div_bal.txt", header=TRUE)

#table to vector
outlier_db <- outliers_db_table[, 1]
View(outliers_db) 

neutral_db <- setdiff(locNames(populations.structure), outlier_db)
length(neutral_db) #72177

write.table(neutral_db, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/neutralloci_div_bal.txt", row.names=F, quote=F)
neutral_db<-read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci_div_bal.txt", header=TRUE)

#Write neutral and outlier genind objects
populations.structure_db_neutral<-populations.structure[loc = neutral_db]
populations.structure_db_outlier<-populations.structure[loc = outlier_db]

View(populations.structure_db_neutral) #72177
View(populations.structure_db_outlier) #3016

# DIV
#Create outlier and netural lists
write.table(outliers_div, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci_div.txt", row.names=F, quote=F)

#NBNBNB Manually remove :+, :- and change : to _ in sublime text using replace function. 
outliers_d_table<-read.table("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/outlierloci_div.txt", header=TRUE)

#table to vector
outlier_d <- outliers_d_table[, 1]
length(outlier_d) #2025
View(outlier_d)

neutral_d <- setdiff(locNames(populations.structure), outlier_d)
length(neutral_d) #73168

write.table(neutral_d, "/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/R_directory/neutralloci_div.txt", row.names=F, quote=F)

#Write neutral and outlier genind objects
populations.structure_d_neutral<-populations.structure[loc = neutral_d]
populations.structure_d_outlier<-populations.structure[loc = outlier_d]

View(populations.structure_d_neutral) #73168
View(populations.structure_d_outlier) #2025

# PCA in adgenet () ####

if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#NB need to manually at ".stru" to your structure file from Stacks output and delete line 1 (header)
populations.structure<-read.structure("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.structure.stru", 
                                      n.ind = 92, n.loc = 75193, col.lab = 1, col.pop = 2, col.others = NULL, row.marknames = 1)
# Which row contains the marker names ('0' if absent)?   1
#Are genotypes coded by a single row (y/n)?   n
  #if I say y, then 
  #Error in mat[, (ncol(mat) - p + 1):ncol(mat)] : 
  #only 0's may be mixed with negative subscripts

View(populations.structure)

  # PCA DIV_BAL ####
#subset continents genind objects
sa_pops <- c("K", "PP", "CC", "S", "M")
aus_pops <- c("LJP", "TS", "SR")

sa_genind <- populations.structure[pop(populations.structure) %in% sa_pops, ]
aus_genind <- populations.structure[pop(populations.structure) %in% aus_pops, ]

sa_genind_db_neutral <- populations.structure_db_neutral[pop(populations.structure_db_neutral) %in% sa_pops, ]
aus_genind_db_neutral <- populations.structure_db_neutral[pop(populations.structure_db_neutral) %in% aus_pops, ]

sa_genind_db_outlier <- populations.structure_db_outlier[pop(populations.structure_db_outlier) %in% sa_pops, ]
aus_genind_db_outlier <- populations.structure_db_outlier[pop(populations.structure_db_outlier) %in% aus_pops, ]

#FULL
tab_db_full <- tab(populations.structure, freq=TRUE, NA.method="mean")
pca_db_full <- dudi.pca(tab_db_full, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_db_full_sa <- tab(sa_genind, freq=TRUE, NA.method="mean")
pca_db_full_sa <- dudi.pca(tab_db_full_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_db_full_aus <- tab(aus_genind, freq=TRUE, NA.method="mean")
pca_db_full_aus <- dudi.pca(tab_db_full_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#NEUTRAL
tab_db_neutral <- tab(populations.structure_db_neutral, freq=TRUE, NA.method="mean")
pca_db_neutral <- dudi.pca(tab_db_full_neutral, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_db_neutral_sa <- tab(sa_genind_db_neutral, freq=TRUE, NA.method="mean")
pca_db_neutra_sa <- dudi.pca(tab_db_neutral_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_db_neutral_aus <- tab(aus_genind_db_neutral, freq=TRUE, NA.method="mean")
pca_db_neutral_aus <- dudi.pca(tab_db_neutral_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#OUTLIER
tab_db_outlier <- tab(populations.structure_db_outlier, freq=TRUE, NA.method="mean")
pca_db_outlier <- dudi.pca(tab_db_outlier, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_db_outlier_sa <- tab(sa_genind_db_outlier, freq=TRUE, NA.method="mean")
pca_db_outlier_sa <- dudi.pca(tab_db_outlier_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_db_outlier_aus <- tab(aus_genind_db_outlier, freq=TRUE, NA.method="mean")
pca_db_outlier_aus <- dudi.pca(tab_db_outlier_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#axes 1 and 2
#full
s.class(pca_db_full$li, 
        fac=pop(populations.structure),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1)
        title="pca_db_full"

#full SA
s.class(pca_db_full_sa$li, 
        fac=pop(sa_genind),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#full AUS
s.class(pca_db_full_aus$li, 
        fac=pop(aus_genind),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 


#netural
s.class(pca_db_neutral$li, 
        fac=pop(populations.structure_db_neutral),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#neutral SA
s.class(pca_db_neutra_sa$li, 
        fac=pop(sa_genind_db_neutral),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#neutral AUS
s.class(pca_db_full_aus$li, 
        fac=pop(aus_genind_db_outlier),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier 
s.class(pca_db_outlier$li, 
        fac=pop(populations.structure_db_outlier),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier SA
s.class(pca_db_outlier_sa$li, 
        fac=pop(sa_genind_outlier),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier AUS
s.class(pca_db_outlier_aus$li, 
        fac=pop(aus_genind_outlier),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 



  # PCA DIV ####
#subset continents genind objects
sa_pops <- c("K", "PP", "CC", "S", "M")
aus_pops <- c("LJP", "TS", "SR")

sa_genind <- populations.structure[pop(populations.structure) %in% sa_pops, ]
aus_genind <- populations.structure[pop(populations.structure) %in% aus_pops, ]

sa_genind_d_neutral <- populations.structure_d_neutral[pop(populations.structure_d_neutral) %in% sa_pops, ]
aus_genind_d_neutral <- populations.structure_d_neutral[pop(populations.structure_d_neutral) %in% aus_pops, ]

sa_genind_d_outlier <- populations.structure_d_outlier[pop(populations.structure_d_outlier) %in% sa_pops, ]
aus_genind_d_outlier <- populations.structure_d_outlier[pop(populations.structure_d_outlier) %in% aus_pops, ]

#FULL
tab_d_full <- tab(populations.structure, freq=TRUE, NA.method="mean")
pca_d_full <- dudi.pca(tab_d_full, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_d_full_sa <- tab(sa_genind, freq=TRUE, NA.method="mean")
pca_d_full_sa <- dudi.pca(tab_d_full_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_d_full_aus <- tab(aus_genind, freq=TRUE, NA.method="mean")
pca_d_full_aus <- dudi.pca(tab_d_full_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#NEUTRAL
tab_d_neutral <- tab(populations.structure_d_neutral, freq=TRUE, NA.method="mean")
pca_d_neutral <- dudi.pca(tab_d_neutral, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_d_neutral_sa <- tab(sa_genind_d_neutral, freq=TRUE, NA.method="mean")
pca_d_neutral_sa <- dudi.pca(tab_d_neutral_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_d_neutral_aus <- tab(aus_genind_d_neutral, freq=TRUE, NA.method="mean")
pca_d_neutral_aus <- dudi.pca(tab_d_neutral_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#OUTLIER
tab_d_outlier <- tab(populations.structure_d_outlier, freq=TRUE, NA.method="mean")
pca_d_outlier <- dudi.pca(tab_d_outlier, center=TRUE, scale=FALSE, nf=4, scannf = FALSE)

tab_d_outlier_sa <- tab(sa_genind_d_outlier, freq=TRUE, NA.method="mean")
pca_d_outlier_sa <- dudi.pca(tab_d_outlier_sa, scannf = FALSE, nf = 2) #Adjust nf for the number of axes you want

tab_d_outlier_aus <- tab(aus_genind_d_outlier, freq=TRUE, NA.method="mean")
pca_d_outlier_aus <- dudi.pca(tab_d_outlier_aus, scannf = FALSE, nf = 2)  #Adjust nf for the number of axes you want

#axes 1 and 2
#full
s.class(pca_d_full$li, 
        fac=pop(populations.structure),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1)

#full SA
s.class(pca_d_full_sa$li, 
        fac=pop(sa_genind),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#full AUS
s.class(pca_d_full_aus$li, 
        fac=pop(aus_genind),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 


#netural
s.class(pca_d_neutral$li, 
        fac=pop(populations.structure_d_neutral),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#neutral SA
s.class(pca_d_neutral_sa$li, 
        fac=pop(sa_genind_d_neutral),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#neutral AUS
s.class(pca_d_full_aus$li, 
        fac=pop(aus_genind_d_outlier),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier 
s.class(pca_d_outlier$li, 
        fac=pop(populations.structure_d_outlier),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154","blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier SA
s.class(pca_d_outlier_sa$li, 
        fac=pop(sa_genind_outlier),
        col=transp((values=c("#addc30","#21918c","#CC7722","#3b528b","#440154")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 

#outlier AUS
s.class(pca_d_outlier_aus$li, 
        fac=pop(aus_genind_outlier),
        col=transp((values=c("blue","#5ec962","#fde725")),.8),
        axesel=FALSE, cstar=0, cpoint=3, clabel=1) 




# POPULATION STATS ####
#Loading Data and Packages 
#Already loaded? 
#adegenet for read.structure
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#NB need to manually at ".stru" to your structure file from Stacks output and delete line 1 (header)
populations.structure<-read.structure("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.structure.stru") 
#How many genotypes are there? 92
#How many markers are there? 75193
#Which column contains labels for genotypes ('0' if absent)?   1
#Which column contains the population factor ('0' if absent)?   2
#Which other optional columns should be read (press 'return' when done)? ENTER 
#  Which row contains the marker names ('0' if absent)?   1
#Are genotypes coded by a single row (y/n)?   n
#if I say y, then 
#Error in mat[, (ncol(mat) - p + 1):ncol(mat)] : 
#only 0's may be mixed with negative subscripts

#For investigating metrics of diversity, load (hierfstat) 
#https://cran.r-project.org/web/packages/hierfstat/vignettes/hierfstat.html
if (!require("hierfstat", quietly = TRUE))
  install.packages("hierfstat")
library(hierfstat)

#write a standard error function
se <- function(x) sd(x)/sqrt(length(x))

#Basic stats 
basic.stats.data <-basic.stats(populations.structure,diploid=TRUE,digits=4) #Takes long. Load from computer.

# Save the object to a file on your computer
save(basic.stats.data, file = "basic_stats_data.RData")

# Load the saved object from the file
load("basic_stats_data.RData")

#View the data
basic.stats.data

$overall
Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.0321 0.3695 0.3928 0.0233 0.3962 0.0266 0.0593 0.0672 0.9131 0.0422 

#Ho
  #FULL  ####
    #Ho ####
#mean observed heterozygosities

Ho(populations.structure) #takes long
K        PP        CC         M         S       LJP        TS        SR 
0.3914727 0.3931838 0.3900177 0.3683588 0.3962536 0.3478815 0.3438643 0.3345483 

Ho<-basic.stats.data[["Ho"]]
View(Ho)

#Prep data for se calulation 
Ho_PP<-Ho[,2,drop=FALSE]
Ho_CC<-Ho[,3,drop=FALSE]
Ho_K<-Ho[,1,drop=FALSE]
Ho_S<-Ho[,5,drop=FALSE]
Ho_M<-Ho[,4,drop=FALSE]
Ho_LJP<-Ho[,6,drop=FALSE]
Ho_TS<-Ho[,7,drop=FALSE]
Ho_SR<-Ho[,8,drop=FALSE]

#se calculation
se(Ho_PP) #0.0004655942
se(Ho_CC) #0.0004684425
se(Ho_K) #0.0004706186
se(Ho_S) #0.0005598928
se(Ho_M) #0.0005779488
se(Ho_LJP) #0.00045651
se(Ho_TS) #0.000466359
se(Ho_SR) #0.0004715281

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Ho_PP, main="Histogram of Ho_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_PP)
qqline(Ho_PP, col = "red")
#PP not normal

hist(Ho_CC, main="Histogram of Ho_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_CC)
qqline(Ho_CC, col = "red")
#CC not normal

hist(Ho_K, main="Histogram of Ho_K", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_K)
qqline(Ho_K, col = "red")
#K not normal

hist(Ho_S, main="Histogram of Ho_S", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_S)
qqline(Ho_S, col = "red")
#S not normal

hist(Ho_M, main="Histogram of Ho_M", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_M)
qqline(Ho_M, col = "red")
#M not normal

hist(Ho_LJP, main="Histogram of Ho_LJP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_LJP)
qqline(Ho_LJP, col = "red")
#LJP not normal

hist(Ho_TS, main="Histogram of Ho_TS", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_TS)
qqline(Ho_TS, col = "red")
#TS not normal

hist(Ho_SR, main="Histogram of Ho_SR", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_SR)
qqline(Ho_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Ho_df <- as.data.frame(Ho)

Ho_long <- Ho_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Ho") # Name of the new 'Ho' column
# View the reshaped data
View(Ho_long)

kruskal.test(Ho ~ factor(population), data = Ho_long)
#Kruskal-Wallis rank sum test
#data:  Ho by factor(population)
#Kruskal-Wallis chi-squared = 11816, df = 7, p-value < 2.2e-16

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test
dunn_results_Ho <- dunn.test(Ho_long$Ho, Ho_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

# Create a data frame from the results
dunn_results_Ho_df <- data.frame(
  Comparison = dunn_results_Ho$comparisons,  # Pairwise comparisons
  Z = dunn_results_Ho$Z,                     # Z-statistics
  P.adjusted = dunn_results_Ho$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Ho)

#All comparisons are significant except: PP-CC, PP-K, M-S
Comparison          Z    P.adjusted
1      CC - K   2.616448  4.975591e-03 
2    CC - LJP  35.392957 2.359934e-274
3     K - LJP  32.776508 1.113805e-235
4      CC - M  67.793734  0.000000e+00
5       K - M  65.177286  0.000000e+00
6     LJP - M  32.400778 2.204045e-230
7     CC - PP   1.201407  1.147966e-01
8      K - PP  -1.415041  8.456891e-02
9    LJP - PP -34.191549 3.228905e-256
10     M - PP -66.592327  0.000000e+00
11     CC - S  69.074825  0.000000e+00
12      K - S  66.458377  0.000000e+00
13    LJP - S  33.681869 9.943548e-249
14      M - S   1.281091  1.037875e-01
15     PP - S  67.873418  0.000000e+00
16    CC - SR  42.944612  0.000000e+00
17     K - SR  40.328164  0.000000e+00
18   LJP - SR   7.551656  2.615926e-14
19     M - SR -24.849122 1.947475e-136
20    PP - SR  41.743205  0.000000e+00
21     S - SR -26.130213 1.286713e-150
22    CC - TS  49.450234  0.000000e+00
23     K - TS  46.833786  0.000000e+00
24   LJP - TS  14.057277  4.423326e-45
25     M - TS -18.343500  2.480343e-75
26    PP - TS  48.248826  0.000000e+00
27     S - TS -19.624591  6.673829e-86
28    SR - TS   6.505622  4.513413e-11

Kruskal-Wallis rank sum test

data: x and group
Kruskal-Wallis chi-squared = 11816.0742, df = 7, p-value = 0


Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
          K |   2.616448
            |    0.0050*
            |
        LJP |   35.39295   32.77650
            |    0.0000*    0.0000*
            |
          M |   67.79373   65.17728   32.40077
            |    0.0000*    0.0000*    0.0000*
            |
         PP |   1.201407  -1.415040  -34.19154  -66.59232
            |     0.1148     0.0846    0.0000*    0.0000*
            |
          S |   69.07482   66.45837   33.68186   1.281091   67.87341
            |    0.0000*    0.0000*    0.0000*     0.1038    0.0000*
            |
         SR |   42.94461   40.32816   7.551655  -24.84912   41.74320  -26.13021
            |    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
            |
         TS |   49.45023   46.83378   14.05727  -18.34350   48.24882  -19.62459
            |    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
   Col Mean-|
   Row Mean |         SR
---------+-----------
         TS |   6.505621
            |    0.0000*


    #Hs ####
#mean gene diversities within population Hs
Hs(populations.structure) #Takes long
K         PP         CC          M          S        LJP         TS         SR 
0.03492746 0.03434253 0.03488855 0.03544213 0.03369330 0.02753722 0.02755842 0.02837595 

Hs<-basic.stats.data[["Hs"]]

#Prep data for se calulation 
Hs_PP<-Hs[,2,drop=FALSE]
Hs_CC<-Hs[,3,drop=FALSE]
Hs_K<-Hs[,1,drop=FALSE]
Hs_S<-Hs[,5,drop=FALSE]
Hs_M<-Hs[,4,drop=FALSE]
Hs_LJP<-Hs[,6,drop=FALSE]
Hs_TS<-Hs[,7,drop=FALSE]
Hs_SR<-Hs[,8,drop=FALSE]

#se calulation 
se(Hs_PP) #0.0007119324
se(Hs_CC) #0.0007186348
se(Hs_K) #0.0007260384
se(Hs_S) #0.001746817
se(Hs_M) #0.001717899
se(Hs_LJP) #0.000776976
se(Hs_TS) #0.0008557705
se(Hs_SR) #0.000828558

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Hs_PP, main="Histogram of Hs_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_PP)
qqline(Hs_PP, col = "red")
#PP not normal

hist(Hs_CC, main="Histogram of Hs_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_CC)
qqline(Hs_CC, col = "red")
#CC not normal

hist(Hs_K, main="Histogram of Hs_K", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_K)
qqline(Hs_K, col = "red")
#K not normal

hist(Hs_S, main="Histogram of Hs_S", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_S)
qqline(Hs_S, col = "red")
#S not normal

hist(Hs_M, main="Histogram of Hs_M", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_M)
qqline(Hs_M, col = "red")
#M not normal

hist(Hs_LJP, main="Histogram of Hs_LJP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_LJP)
qqline(Hs_LJP, col = "red")
#LJP not normal

hist(Hs_TS, main="Histogram of Hs_TS", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_TS)
qqline(Hs_TS, col = "red")
#TS not normal

hist(Hs_SR, main="Histogram of Hs_SR", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_SR)
qqline(Hs_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Hs_df <- as.data.frame(Hs)

Hs_long <- Hs_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Hs") # Name of the new 'Hs' column
# View the reshaped data
View(Hs_long)

kruskal.test(Hs ~ factor(population), data = Hs_long)
#	Kruskal-Wallis rank sum test
#data:  Hs by factor(population)
#Kruskal-Wallis chi-squared = 12134, df = 7, p-value < 2.2e-16

#PostHoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#Dunn's test
dunn_results_Hs <- dunn.test(Hs_long$Hs, Hs_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

Kruskal-Wallis rank sum test
data: x and group
Kruskal-Wallis chi-squared = 12134.2268, df = 7, p-value = 0

Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
         K |  -2.303996
           |    0.0110*
           |
       LJP |   31.90622   34.21021
          |    0.0000*    0.0000*
          |
          M |   74.97292   77.27692   43.06670
          |    0.0000*    0.0000*    0.0000*
          |
        PP |  -2.494391  -0.190395  -34.40061  -77.46732
          |    0.0068*     0.4245    0.0000*    0.0000*
          |
          S |   59.52958   61.83358   27.62336  -15.44334   62.02397
          |    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
          |
        SR |   42.13360   44.43760   10.22738  -32.83932   44.62800  -17.39597
          |    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
          |
        TS |   34.74814   37.05214   2.841928  -40.22477   37.24254  -24.78143
          |    0.0000*    0.0000*    0.0025*    0.0000*    0.0000*    0.0000*
  Col Mean-|
  Row Mean |         SR
---------+-----------
        TS |  -7.385458
          |    0.0000*
  
  alpha = 0.05
Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Hs_df <- data.frame(
  Comparison = dunn_results_Hs$comparisons,  # Pairwise comparisons
  Z = dunn_results_Hs$Z,                     # Z-statistics
  P.adjusted = dunn_results_Hs$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Hs_df)
Comparison           Z    P.adjusted
1      CC - K  -2.3039960  1.100445e-02
2    CC - LJP  31.9062208 1.702548e-223
3     K - LJP  34.2102168 1.491253e-256
4      CC - M  74.9729295  0.000000e+00
5       K - M  77.2769255  0.000000e+00
6     LJP - M  43.0667087  0.000000e+00
7     CC - PP  -2.4943915  6.793946e-03
8      K - PP  -0.1903955  4.244996e-01
9    LJP - PP -34.4006123 2.304291e-259
10     M - PP -77.4673210  0.000000e+00
11     CC - S  59.5295874  0.000000e+00
12      K - S  61.8335835  0.000000e+00
13    LJP - S  27.6233666 4.297193e-168
14      M - S -15.4433421  5.323247e-54
15     PP - S  62.0239790  0.000000e+00
16    CC - SR  42.1336087  0.000000e+00
17     K - SR  44.4376047  0.000000e+00
18   LJP - SR  10.2273879  9.099840e-25
19     M - SR -32.8393209 1.332538e-236
20    PP - SR  44.6280002  0.000000e+00
21     S - SR -17.3959788  5.900129e-68
22    CC - TS  34.7481498 1.478108e-264
23     K - TS  37.0521458 1.786103e-300
24   LJP - TS   2.8419290  2.511123e-03
25     M - TS -40.2247797  0.000000e+00
26    PP - TS  37.2425413 1.632364e-303
27     S - TS -24.7814376 9.949912e-136
28    SR - TS  -7.3854589  8.862496e-14

#All comparisons are significant except: PP-K

    #FIS ####
#Inbreeding Coeficcient 

#subset data to find FIS using basic.stats
#poppr for popsub
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

populations.structure_PP<-
  popsub(
    populations.structure,
    sublist = "PP",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_CC<-
  popsub(
    populations.structure,
    sublist = "CC",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_S<-
  popsub(
    populations.structure,
    sublist = "S",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_K<-
  popsub(
    populations.structure,
    sublist = "K",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_M<-
  popsub(
    populations.structure,
    sublist = "M",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_LJP<-
  popsub(
    populations.structure,
    sublist = "LJP",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_TS<-
  popsub(
    populations.structure,
    sublist = "TS",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_SR<-
  popsub(
    populations.structure,
    sublist = "SR",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE)

#Find FIS using basic.stats
basic.stats_PP<-basic.stats(populations.structure_PP,diploid=TRUE,digits=5)
basic.stats_PP
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03435 0.39319 0.39319 0.00000     NaN     NaN 0.00000     NaN 0.91265     NaN 

FIS_PP #0.91265

basic.stats_CC <-basic.stats(populations.structure_CC,diploid=TRUE,digits=5)
basic.stats_CC
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03489 0.39002 0.39002 0.00000     NaN     NaN 0.00000     NaN 0.91055     NaN 

FIS_CC #0.91055

basic.stats_S <-basic.stats(populations.structure_S,diploid=TRUE,digits=5)
basic.stats_S
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03369 0.39625 0.39625 0.00000     NaN     NaN 0.00000     NaN 0.91497     NaN 

FIS_S #0.91497

basic.stats_K <-basic.stats(populations.structure_K,diploid=TRUE,digits=5)
basic.stats_K
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03493 0.39148 0.39148 0.00000     NaN     NaN 0.00000     NaN 0.91078     NaN 

FIS_K #0.91078

basic.stats_M <-basic.stats(populations.structure_M,diploid=TRUE,digits=5)
basic.stats_M
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03544 0.36836 0.36836 0.00000     NaN     NaN 0.00000     NaN 0.90378     NaN 

FIS_M #0.90378

basic.stats_LJP <-basic.stats(populations.structure_LJP,diploid=TRUE,digits=5)
basic.stats_LJP
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02754 0.34788 0.34788 0.00000     NaN     NaN 0.00000     NaN 0.92085     NaN 

FIS_LJP #0.92085

basic.stats_TS <-basic.stats(populations.structure_TS,diploid=TRUE,digits=5)
basic.stats_TS
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02756 0.34386 0.34386 0.00000     NaN     NaN 0.00000     NaN 0.91986     NaN 

FIS_TS #0.91986

basic.stats_SR <-basic.stats(populations.structure_SR,diploid=TRUE,digits=5)
basic.stats_SR
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02838 0.33455 0.33455 0.00000     NaN     NaN 0.00000     NaN 0.91518     NaN

FIS_SR #0.91518

#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=populations.structure,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
View(populations.structure_boot_fis)
populations.structure_boot_fis
$fis.ci
ll     hl
K 1   0.9086 0.9129
PP 2  0.9104 0.9150
CC 3  0.9084 0.9124
M 4   0.9008 0.9064
S 5   0.9120 0.9177
LJP 6 0.9182 0.9231
TS 7  0.9168 0.9225
SR 8  0.9120 0.9179

#Subset datasets to prepare to test nomality 
#Can probably also use basic.stats_PP$FIS
FIS<-basic.stats.data[["Fis"]]
View(FIS)
View(basic.stats.data)

FIS_PP<-FIS[,2,drop=FALSE]
FIS_CC<-FIS[,3,drop=FALSE]
FIS_K<-FIS[,1,drop=FALSE]
FIS_S<-FIS[,5,drop=FALSE]
FIS_M<-FIS[,4,drop=FALSE]
FIS_LJP<-FIS[,6,drop=FALSE]
FIS_TS<-FIS[,7,drop=FALSE]
FIS_SR<-FIS[,8,drop=FALSE]

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(FIS_PP, main="Histogram of FIS_PP", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_PP)
qqline(FIS_PP, col = "red")
#PP not normal

hist(FIS_CC, main="Histogram of FIS_CC", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_CC)
qqline(FIS_CC, col = "red")
#CC not normal

hist(FIS_K, main="Histogram of FIS_K", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_K)
qqline(FIS_K, col = "red")
#K not normal

hist(FIS_S, main="Histogram of FIS_S", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_S)
qqline(FIS_S, col = "red")
#S not normal

hist(FIS_M, main="Histogram of FIS_M", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_M)
qqline(FIS_M, col = "red")
#M not normal

hist(FIS_LJP, main="Histogram of FIS_LJP", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_LJP)
qqline(FIS_LJP, col = "red")
#LJP not normal

hist(FIS_TS, main="Histogram of FIS_TS", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_TS)
qqline(FIS_TS, col = "red")
#TS not normal

hist(FIS_SR, main="Histogram of FIS_SR", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_SR)
qqline(FIS_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
FIS_df <- as.data.frame(FIS)

FIS_long <- FIS_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "FIS") # Name of the new 'FIS' column
# View the reshaped data
View(FIS_long)

kruskal.test(FIS ~ factor(population), data = FIS_long)
#	Kruskal-Wallis rank sum test
#data:  FIS by factor(population)
#Kruskal-Wallis chi-squared = 2669.9, df = 7, p-value < 2.2e-16

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test
dunn_results_FIS <- dunn.test(FIS_long$FIS, FIS_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

# Create a data frame from the results
dunn_results_FIS_df <- data.frame(
  Comparison = dunn_results_FIS$comparisons,  # Pairwise comparisons
  Z = dunn_results_FIS$Z,                     # Z-statistics
  P.adjusted = dunn_results_FIS$P.adjusted    # Adjusted p-values
)

Kruskal-Wallis rank sum test

data: x and group
Kruskal-Wallis chi-squared = 2669.8598, df = 7, p-value = 0


Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
          K |  -2.599589
            |    0.0054*
            |
        LJP |  -27.26165  -24.67677
            |    0.0000*    0.0000*
            |
          M |  -13.43958  -11.43235   7.789322
            |    0.0000*    0.0000*    0.0000*
            |
         PP |  -1.282677   1.323427   26.05244   12.47031
            |     0.0998     0.0963    0.0000*    0.0000*
            |
          S |  -19.75971  -17.70235   2.011939  -5.032265  -18.77635
            |    0.0000*    0.0000*    0.0238*    0.0000*    0.0000*
            |
        SR |  -29.30369  -26.77397  -2.599168  -9.746696  -28.12371  -4.089740
            |    0.0000*    0.0000*    0.0052*    0.0000*    0.0000*    0.0000*
            |
          TS |  -34.36567  -31.84765  -7.756127  -13.81973  -33.20220  -8.273210
            |    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  Col Mean-|
  Row Mean |         SR
---------+-----------
          TS |  -5.069313
             |    0.0000*
  
  alpha = 0.05
Reject Ho if p <= alpha/2

# Print the data frame
print(dunn_results_FIS)

#All comparisons are significant except: PP-CC, PP-K 

    #FST ####
#FST from Stamp 
if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

#MDependancies for dartR (not installing without them)
BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)


#FULL
populations.structure.gl<-gi2gl(populations.structure)
stamppFst_values <- stamppFst(populations.structure.gl, nboots = 100, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

# Save the object to a file on your computer
save(stamppFst_values, file = "stamppFst_values.RData")

# Load the saved object from the file
load("stamppFst_values.RData")

options(max.print = 100)
stamppFst_values
View(stamppFst_values)

FST_pvalues<-stamppFst_values[["Pvalues"]]
FST_pvalues

$Fsts
K           PP          CC           M         S         LJP          TS SR
K              NA           NA          NA          NA        NA          NA          NA NA
PP   0.0028593991           NA          NA          NA        NA          NA          NA NA
CC   0.0008211852  0.000500510          NA          NA        NA          NA          NA NA
M   -0.0140177640 -0.014560058 -0.01345659          NA        NA          NA          NA NA
S    0.0087417167  0.009834058  0.01045641 0.004807446        NA          NA          NA NA
LJP  0.1104935683  0.110943954  0.10973024 0.127627279 0.1429823          NA          NA NA
TS   0.1062377864  0.106273228  0.10670518 0.127996934 0.1410762 0.002987782          NA NA
SR   0.1094841204  0.109847457  0.10917434 0.135512272 0.1538608 0.001104700 0.002134911 NA

$Pvalues #I think the decimals are shown as the least necessary 
K   PP CC    M  S  LJP TS SR
K   NA   NA NA   NA NA   NA NA NA
PP   0   NA NA   NA NA   NA NA NA
CC   0 0.09 NA   NA NA   NA NA NA
M    1 1.00  1   NA NA   NA NA NA
S    0 0.00  0 0.14 NA   NA NA NA
LJP  0 0.00  0 0.00  0   NA NA NA
TS   0 0.00  0 0.00  0 0.00 NA NA
SR   0 0.00  0 0.00  0 0.01  0 NA

    #HDW ####

# Read the VCF file
if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

#Check full to see if same as Populations
#data <- read.vcf("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.snps.vcf")
populations.snps.vcf <- read.vcf("populations.snps.vcf")
populations.snps.vcfR <- read.vcfR("populations.snps.vcf")

# Convert VCF to genotype matrix
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#.gen file was unavailable so I used the vcf
genind_full <- vcfR2genind(populations.snps.vcfR)

# Perform HWE test on the genind object using the `pegas` package
if (!require("pegas", quietly = TRUE))
  install.packages("pegas")
library(pegas)

HWE_results <- hw.test(genind_full)
# Write the data frame to a CSV file
write.csv(hwe_results, file = "HWE_results_full.csv", row.names = FALSE)

#Now Neutral
#data <- read.vcf("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92.nosync/populations.snps.vcf")
populations.snps.netural.vcf <- read.vcf("populations.snps.neutral.vcf")

if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

populations.snps.netural.vcfR <- read.vcfR("populations.snps.neutral.vcf")

populations.snps.netural.vcfR <- read.vcfR("~/msc/00_info/populations.snps.neutral.vcf")

# Convert VCF to genotype matrix
if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#.gen file was unavailable so I used the vcf
genind_neutral <- vcfR2genind(populations.snps.netural.vcfR)

# Perform HWE test on the genind object using the `pegas` package
if (!require("pegas", quietly = TRUE))
  install.packages("pegas")
library(pegas)

HWE_results_netural <- hw.test(genind_neutral)
# Write the data frame to a CSV file

write.csv(hwe_results_neutral, file = "HWE_results_neutral.csv", row.names = FALSE)
 --> Run in hpc
#OUTLIER  ####
if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

outlier.db.vcf <- read.vcf("populations.snps.outlier.db.vcf")

outlier.stru <- vcf2structure(outlier.db.vcf, out = "outlier.db.stru")

if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

#NB need to manually at ".stru" to your structure file from Stacks output and delete line 1 (header)
populations.structure<-read.structure("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Bioinformatics/populations_outputs_92/populations.structure.stru") 



    #Ho ####
#mean observed heterozygosities

Ho(populations.structure) #takes long
K        PP        CC         M         S       LJP        TS        SR 
0.3914727 0.3931838 0.3900177 0.3683588 0.3962536 0.3478815 0.3438643 0.3345483 

Ho<-basic.stats.data[["Ho"]]
View(Ho)

#Prep data for se calulation 
Ho_PP<-Ho[,2,drop=FALSE]
Ho_CC<-Ho[,3,drop=FALSE]
Ho_K<-Ho[,1,drop=FALSE]
Ho_S<-Ho[,5,drop=FALSE]
Ho_M<-Ho[,4,drop=FALSE]
Ho_LJP<-Ho[,6,drop=FALSE]
Ho_TS<-Ho[,7,drop=FALSE]
Ho_SR<-Ho[,8,drop=FALSE]

#se calculation
se(Ho_PP) #0.0004655942
se(Ho_CC) #0.0004684425
se(Ho_K) #0.0004706186
se(Ho_S) #0.0005598928
se(Ho_M) #0.0005779488
se(Ho_LJP) #0.00045651
se(Ho_TS) #0.000466359
se(Ho_SR) #0.0004715281

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Ho_PP, main="Histogram of Ho_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_PP)
qqline(Ho_PP, col = "red")
#PP not normal

hist(Ho_CC, main="Histogram of Ho_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_CC)
qqline(Ho_CC, col = "red")
#CC not normal

hist(Ho_K, main="Histogram of Ho_K", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_K)
qqline(Ho_K, col = "red")
#K not normal

hist(Ho_S, main="Histogram of Ho_S", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_S)
qqline(Ho_S, col = "red")
#S not normal

hist(Ho_M, main="Histogram of Ho_M", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_M)
qqline(Ho_M, col = "red")
#M not normal

hist(Ho_LJP, main="Histogram of Ho_LJP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_LJP)
qqline(Ho_LJP, col = "red")
#LJP not normal

hist(Ho_TS, main="Histogram of Ho_TS", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_TS)
qqline(Ho_TS, col = "red")
#TS not normal

hist(Ho_SR, main="Histogram of Ho_SR", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_SR)
qqline(Ho_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Ho_df <- as.data.frame(Ho)

Ho_long <- Ho_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Ho") # Name of the new 'Ho' column
# View the reshaped data
View(Ho_long)

kruskal.test(Ho ~ factor(population), data = Ho_long)
#Kruskal-Wallis rank sum test
#data:  Ho by factor(population)
#Kruskal-Wallis chi-squared = 11816, df = 7, p-value < 2.2e-16

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test
dunn_results_Ho <- dunn.test(Ho_long$Ho, Ho_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

# Create a data frame from the results
dunn_results_Ho_df <- data.frame(
  Comparison = dunn_results_Ho$comparisons,  # Pairwise comparisons
  Z = dunn_results_Ho$Z,                     # Z-statistics
  P.adjusted = dunn_results_Ho$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Ho)

#All comparisons are significant except: PP-CC, PP-K, M-S
Comparison          Z    P.adjusted
1      CC - K   2.616448  4.975591e-03 
2    CC - LJP  35.392957 2.359934e-274
3     K - LJP  32.776508 1.113805e-235
4      CC - M  67.793734  0.000000e+00
5       K - M  65.177286  0.000000e+00
6     LJP - M  32.400778 2.204045e-230
7     CC - PP   1.201407  1.147966e-01
8      K - PP  -1.415041  8.456891e-02
9    LJP - PP -34.191549 3.228905e-256
10     M - PP -66.592327  0.000000e+00
11     CC - S  69.074825  0.000000e+00
12      K - S  66.458377  0.000000e+00
13    LJP - S  33.681869 9.943548e-249
14      M - S   1.281091  1.037875e-01
15     PP - S  67.873418  0.000000e+00
16    CC - SR  42.944612  0.000000e+00
17     K - SR  40.328164  0.000000e+00
18   LJP - SR   7.551656  2.615926e-14
19     M - SR -24.849122 1.947475e-136
20    PP - SR  41.743205  0.000000e+00
21     S - SR -26.130213 1.286713e-150
22    CC - TS  49.450234  0.000000e+00
23     K - TS  46.833786  0.000000e+00
24   LJP - TS  14.057277  4.423326e-45
25     M - TS -18.343500  2.480343e-75
26    PP - TS  48.248826  0.000000e+00
27     S - TS -19.624591  6.673829e-86
28    SR - TS   6.505622  4.513413e-11

Kruskal-Wallis rank sum test

data: x and group
Kruskal-Wallis chi-squared = 11816.0742, df = 7, p-value = 0


Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
  K |   2.616448
|    0.0050*
  |
  LJP |   35.39295   32.77650
|    0.0000*    0.0000*
  |
  M |   67.79373   65.17728   32.40077
|    0.0000*    0.0000*    0.0000*
  |
  PP |   1.201407  -1.415040  -34.19154  -66.59232
|     0.1148     0.0846    0.0000*    0.0000*
  |
  S |   69.07482   66.45837   33.68186   1.281091   67.87341
|    0.0000*    0.0000*    0.0000*     0.1038    0.0000*
  |
  SR |   42.94461   40.32816   7.551655  -24.84912   41.74320  -26.13021
|    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  |
  TS |   49.45023   46.83378   14.05727  -18.34350   48.24882  -19.62459
|    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  Col Mean-|
  Row Mean |         SR
---------+-----------
  TS |   6.505621
|    0.0000*
  
  
    #Hs ####
#mean gene diversities within population Hs
Hs(populations.structure) #Takes long
K         PP         CC          M          S        LJP         TS         SR 
0.03492746 0.03434253 0.03488855 0.03544213 0.03369330 0.02753722 0.02755842 0.02837595 

Hs<-basic.stats.data[["Hs"]]

#Prep data for se calulation 
Hs_PP<-Hs[,2,drop=FALSE]
Hs_CC<-Hs[,3,drop=FALSE]
Hs_K<-Hs[,1,drop=FALSE]
Hs_S<-Hs[,5,drop=FALSE]
Hs_M<-Hs[,4,drop=FALSE]
Hs_LJP<-Hs[,6,drop=FALSE]
Hs_TS<-Hs[,7,drop=FALSE]
Hs_SR<-Hs[,8,drop=FALSE]

#se calulation 
se(Hs_PP) #0.0007119324
se(Hs_CC) #0.0007186348
se(Hs_K) #0.0007260384
se(Hs_S) #0.001746817
se(Hs_M) #0.001717899
se(Hs_LJP) #0.000776976
se(Hs_TS) #0.0008557705
se(Hs_SR) #0.000828558

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Hs_PP, main="Histogram of Hs_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_PP)
qqline(Hs_PP, col = "red")
#PP not normal

hist(Hs_CC, main="Histogram of Hs_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_CC)
qqline(Hs_CC, col = "red")
#CC not normal

hist(Hs_K, main="Histogram of Hs_K", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_K)
qqline(Hs_K, col = "red")
#K not normal

hist(Hs_S, main="Histogram of Hs_S", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_S)
qqline(Hs_S, col = "red")
#S not normal

hist(Hs_M, main="Histogram of Hs_M", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_M)
qqline(Hs_M, col = "red")
#M not normal

hist(Hs_LJP, main="Histogram of Hs_LJP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_LJP)
qqline(Hs_LJP, col = "red")
#LJP not normal

hist(Hs_TS, main="Histogram of Hs_TS", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_TS)
qqline(Hs_TS, col = "red")
#TS not normal

hist(Hs_SR, main="Histogram of Hs_SR", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_SR)
qqline(Hs_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Hs_df <- as.data.frame(Hs)

Hs_long <- Hs_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Hs") # Name of the new 'Hs' column
# View the reshaped data
View(Hs_long)

kruskal.test(Hs ~ factor(population), data = Hs_long)
#	Kruskal-Wallis rank sum test
#data:  Hs by factor(population)
#Kruskal-Wallis chi-squared = 12134, df = 7, p-value < 2.2e-16

#PostHoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#Dunn's test
dunn_results_Hs <- dunn.test(Hs_long$Hs, Hs_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

Kruskal-Wallis rank sum test
data: x and group
Kruskal-Wallis chi-squared = 12134.2268, df = 7, p-value = 0

Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
  K |  -2.303996
|    0.0110*
  |
  LJP |   31.90622   34.21021
|    0.0000*    0.0000*
  |
  M |   74.97292   77.27692   43.06670
|    0.0000*    0.0000*    0.0000*
  |
  PP |  -2.494391  -0.190395  -34.40061  -77.46732
|    0.0068*     0.4245    0.0000*    0.0000*
  |
  S |   59.52958   61.83358   27.62336  -15.44334   62.02397
|    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  |
  SR |   42.13360   44.43760   10.22738  -32.83932   44.62800  -17.39597
|    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  |
  TS |   34.74814   37.05214   2.841928  -40.22477   37.24254  -24.78143
|    0.0000*    0.0000*    0.0025*    0.0000*    0.0000*    0.0000*
  Col Mean-|
  Row Mean |         SR
---------+-----------
  TS |  -7.385458
|    0.0000*
  
  alpha = 0.05
Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Hs_df <- data.frame(
  Comparison = dunn_results_Hs$comparisons,  # Pairwise comparisons
  Z = dunn_results_Hs$Z,                     # Z-statistics
  P.adjusted = dunn_results_Hs$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Hs_df)
Comparison           Z    P.adjusted
1      CC - K  -2.3039960  1.100445e-02
2    CC - LJP  31.9062208 1.702548e-223
3     K - LJP  34.2102168 1.491253e-256
4      CC - M  74.9729295  0.000000e+00
5       K - M  77.2769255  0.000000e+00
6     LJP - M  43.0667087  0.000000e+00
7     CC - PP  -2.4943915  6.793946e-03
8      K - PP  -0.1903955  4.244996e-01
9    LJP - PP -34.4006123 2.304291e-259
10     M - PP -77.4673210  0.000000e+00
11     CC - S  59.5295874  0.000000e+00
12      K - S  61.8335835  0.000000e+00
13    LJP - S  27.6233666 4.297193e-168
14      M - S -15.4433421  5.323247e-54
15     PP - S  62.0239790  0.000000e+00
16    CC - SR  42.1336087  0.000000e+00
17     K - SR  44.4376047  0.000000e+00
18   LJP - SR  10.2273879  9.099840e-25
19     M - SR -32.8393209 1.332538e-236
20    PP - SR  44.6280002  0.000000e+00
21     S - SR -17.3959788  5.900129e-68
22    CC - TS  34.7481498 1.478108e-264
23     K - TS  37.0521458 1.786103e-300
24   LJP - TS   2.8419290  2.511123e-03
25     M - TS -40.2247797  0.000000e+00
26    PP - TS  37.2425413 1.632364e-303
27     S - TS -24.7814376 9.949912e-136
28    SR - TS  -7.3854589  8.862496e-14

#All comparisons are significant except: PP-K

    #FIS ####
#Inbreeding Coeficcient 

#subset data to find FIS using basic.stats
#poppr for popsub
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

populations.structure_PP<-
  popsub(
    populations.structure,
    sublist = "PP",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_CC<-
  popsub(
    populations.structure,
    sublist = "CC",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_S<-
  popsub(
    populations.structure,
    sublist = "S",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_K<-
  popsub(
    populations.structure,
    sublist = "K",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_M<-
  popsub(
    populations.structure,
    sublist = "M",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_LJP<-
  popsub(
    populations.structure,
    sublist = "LJP",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_TS<-
  popsub(
    populations.structure,
    sublist = "TS",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

populations.structure_SR<-
  popsub(
    populations.structure,
    sublist = "SR",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE)

#Find FIS using basic.stats
basic.stats_PP<-basic.stats(populations.structure_PP,diploid=TRUE,digits=5)
basic.stats_PP
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03435 0.39319 0.39319 0.00000     NaN     NaN 0.00000     NaN 0.91265     NaN 

FIS_PP #0.91265

basic.stats_CC <-basic.stats(populations.structure_CC,diploid=TRUE,digits=5)
basic.stats_CC
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03489 0.39002 0.39002 0.00000     NaN     NaN 0.00000     NaN 0.91055     NaN 

FIS_CC #0.91055

basic.stats_S <-basic.stats(populations.structure_S,diploid=TRUE,digits=5)
basic.stats_S
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03369 0.39625 0.39625 0.00000     NaN     NaN 0.00000     NaN 0.91497     NaN 

FIS_S #0.91497

basic.stats_K <-basic.stats(populations.structure_K,diploid=TRUE,digits=5)
basic.stats_K
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03493 0.39148 0.39148 0.00000     NaN     NaN 0.00000     NaN 0.91078     NaN 

FIS_K #0.91078

basic.stats_M <-basic.stats(populations.structure_M,diploid=TRUE,digits=5)
basic.stats_M
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.03544 0.36836 0.36836 0.00000     NaN     NaN 0.00000     NaN 0.90378     NaN 

FIS_M #0.90378

basic.stats_LJP <-basic.stats(populations.structure_LJP,diploid=TRUE,digits=5)
basic.stats_LJP
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02754 0.34788 0.34788 0.00000     NaN     NaN 0.00000     NaN 0.92085     NaN 

FIS_LJP #0.92085

basic.stats_TS <-basic.stats(populations.structure_TS,diploid=TRUE,digits=5)
basic.stats_TS
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02756 0.34386 0.34386 0.00000     NaN     NaN 0.00000     NaN 0.91986     NaN 

FIS_TS #0.91986

basic.stats_SR <-basic.stats(populations.structure_SR,diploid=TRUE,digits=5)
basic.stats_SR
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.02838 0.33455 0.33455 0.00000     NaN     NaN 0.00000     NaN 0.91518     NaN

FIS_SR #0.91518

#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=populations.structure,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
View(populations.structure_boot_fis)
populations.structure_boot_fis
$fis.ci
ll     hl
K 1   0.9086 0.9129
PP 2  0.9104 0.9150
CC 3  0.9084 0.9124
M 4   0.9008 0.9064
S 5   0.9120 0.9177
LJP 6 0.9182 0.9231
TS 7  0.9168 0.9225
SR 8  0.9120 0.9179

#Subset datasets to prepare to test nomality 
#Can probably also use basic.stats_PP$FIS
FIS<-basic.stats.data[["Fis"]]
View(FIS)
View(basic.stats.data)

FIS_PP<-FIS[,2,drop=FALSE]
FIS_CC<-FIS[,3,drop=FALSE]
FIS_K<-FIS[,1,drop=FALSE]
FIS_S<-FIS[,5,drop=FALSE]
FIS_M<-FIS[,4,drop=FALSE]
FIS_LJP<-FIS[,6,drop=FALSE]
FIS_TS<-FIS[,7,drop=FALSE]
FIS_SR<-FIS[,8,drop=FALSE]

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(FIS_PP, main="Histogram of FIS_PP", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_PP)
qqline(FIS_PP, col = "red")
#PP not normal

hist(FIS_CC, main="Histogram of FIS_CC", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_CC)
qqline(FIS_CC, col = "red")
#CC not normal

hist(FIS_K, main="Histogram of FIS_K", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_K)
qqline(FIS_K, col = "red")
#K not normal

hist(FIS_S, main="Histogram of FIS_S", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_S)
qqline(FIS_S, col = "red")
#S not normal

hist(FIS_M, main="Histogram of FIS_M", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_M)
qqline(FIS_M, col = "red")
#M not normal

hist(FIS_LJP, main="Histogram of FIS_LJP", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_LJP)
qqline(FIS_LJP, col = "red")
#LJP not normal

hist(FIS_TS, main="Histogram of FIS_TS", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_TS)
qqline(FIS_TS, col = "red")
#TS not normal

hist(FIS_SR, main="Histogram of FIS_SR", xlab="Value", ylab="Frequency") #very left-skewed
qqnorm(FIS_SR)
qqline(FIS_SR, col = "red")
#SR not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
FIS_df <- as.data.frame(FIS)

FIS_long <- FIS_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "FIS") # Name of the new 'FIS' column
# View the reshaped data
View(FIS_long)

kruskal.test(FIS ~ factor(population), data = FIS_long)
#	Kruskal-Wallis rank sum test
#data:  FIS by factor(population)
#Kruskal-Wallis chi-squared = 2669.9, df = 7, p-value < 2.2e-16

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test
dunn_results_FIS <- dunn.test(FIS_long$FIS, FIS_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

# Create a data frame from the results
dunn_results_FIS_df <- data.frame(
  Comparison = dunn_results_FIS$comparisons,  # Pairwise comparisons
  Z = dunn_results_FIS$Z,                     # Z-statistics
  P.adjusted = dunn_results_FIS$P.adjusted    # Adjusted p-values
)

Kruskal-Wallis rank sum test

data: x and group
Kruskal-Wallis chi-squared = 2669.8598, df = 7, p-value = 0


Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |         CC          K        LJP          M         PP          S
---------+------------------------------------------------------------------
  K |  -2.599589
|    0.0054*
  |
  LJP |  -27.26165  -24.67677
|    0.0000*    0.0000*
  |
  M |  -13.43958  -11.43235   7.789322
|    0.0000*    0.0000*    0.0000*
  |
  PP |  -1.282677   1.323427   26.05244   12.47031
|     0.0998     0.0963    0.0000*    0.0000*
  |
  S |  -19.75971  -17.70235   2.011939  -5.032265  -18.77635
|    0.0000*    0.0000*    0.0238*    0.0000*    0.0000*
  |
  SR |  -29.30369  -26.77397  -2.599168  -9.746696  -28.12371  -4.089740
|    0.0000*    0.0000*    0.0052*    0.0000*    0.0000*    0.0000*
  |
  TS |  -34.36567  -31.84765  -7.756127  -13.81973  -33.20220  -8.273210
|    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*    0.0000*
  Col Mean-|
  Row Mean |         SR
---------+-----------
  TS |  -5.069313
|    0.0000*
  
  alpha = 0.05
Reject Ho if p <= alpha/2

# Print the data frame
print(dunn_results_FIS)

#All comparisons are significant except: PP-CC, PP-K 

    #FST ####
#FST from Stamp 
if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

#MDependancies for dartR (not installing without them)
BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)


#FULL
populations.structure.gl<-gi2gl(populations.structure)
stamppFst_values <- stamppFst(populations.structure.gl, nboots = 100, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

# Save the object to a file on your computer
save(stamppFst_values, file = "stamppFst_values.RData")

# Load the saved object from the file
load("stamppFst_values.RData")

options(max.print = 100)
stamppFst_values
View(stamppFst_values)

FST_pvalues<-stamppFst_values[["Pvalues"]]
FST_pvalues

$Fsts
K           PP          CC           M         S         LJP          TS SR
K              NA           NA          NA          NA        NA          NA          NA NA
PP   0.0028593991           NA          NA          NA        NA          NA          NA NA
CC   0.0008211852  0.000500510          NA          NA        NA          NA          NA NA
M   -0.0140177640 -0.014560058 -0.01345659          NA        NA          NA          NA NA
S    0.0087417167  0.009834058  0.01045641 0.004807446        NA          NA          NA NA
LJP  0.1104935683  0.110943954  0.10973024 0.127627279 0.1429823          NA          NA NA
TS   0.1062377864  0.106273228  0.10670518 0.127996934 0.1410762 0.002987782          NA NA
SR   0.1094841204  0.109847457  0.10917434 0.135512272 0.1538608 0.001104700 0.002134911 NA

$Pvalues #I think the decimals are shown as the least necessary 
K   PP CC    M  S  LJP TS SR
K   NA   NA NA   NA NA   NA NA NA
PP   0   NA NA   NA NA   NA NA NA
CC   0 0.09 NA   NA NA   NA NA NA
M    1 1.00  1   NA NA   NA NA NA
S    0 0.00  0 0.14 NA   NA NA NA
LJP  0 0.00  0 0.00  0   NA NA NA
TS   0 0.00  0 0.00  0 0.00 NA NA
SR   0 0.00  0 0.00  0 0.01  0 NA



#Other ####

