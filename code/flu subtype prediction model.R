# This R script is conducted using toydata in Singapore as example.

#install packages:
#typical installation time: 5-6 minutes 
install.packages(c("lubridate", "plyr", "dplyr", "tibble"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

# set working dictionary and library packages
setwd("./dominant subtype prediction/toydata")
library(Biostrings)  
library(lubridate)   
library(plyr)        
library(dplyr)       
library(tibble)  

###########################################################
# Step 1: prepare protein sequence
# The function 'Prepareseq' allows to input an aligned protein sequence data and output a csv format data with flu season

cat("  Step 1: Preparing protein sequences from FASTA files\n")

PrepareSeq <- function(input_fasta_seq){
  fasta_seq <- readAAStringSet(input_fasta_seq)
  name <- strsplit(fasta_seq@ranges@NAMES,split="|",fixed=T)
  seq <- strsplit(as.character(tolower(fasta_seq)),split="")
  aa_name <- c()
  complete_data <- matrix(0,ncol=(length(fasta_seq[[1]])+4),nrow=length(fasta_seq))
  for (i in 1:length(fasta_seq)){
    complete_data[i,1]<-name[[i]][1] #accession number#
    complete_data[i,2]<-name[[i]][3] #isolate name#
    complete_data[i,3]<-name[[i]][2] #date
    if(month(name[[i]][2])<10){
      complete_data[i,4]<-year(name[[i]][2])}else{
        complete_data[i,4]<-year(name[[i]][2])+1
      }
    for (j in 1:length(fasta_seq[[1]])){
      complete_data[i,j+4]<-seq[[i]][j]
    }
  }
  for (k in 1:length(fasta_seq[[1]])){
    aa_name <- c(aa_name,paste("x",k,sep="",collapse=""))
  }
  colnames(complete_data) <- c("accession number","name","date","season",aa_name)
  complete_data <- as.data.frame(complete_data)
  complete_data[complete_data=="-"] <- NA
  return(complete_data)
}

h1_seq <- PrepareSeq(input_fasta_seq = "SGH1_HA_sequence.fasta")
cat(sprintf("  H1N1: %d sequences, %d amino acid positions, seasons %d-%d\n",
    nrow(h1_seq), ncol(h1_seq) - 4,
    min(as.numeric(h1_seq$season), na.rm = TRUE),
    max(as.numeric(h1_seq$season), na.rm = TRUE)))

h3_seq <- PrepareSeq(input_fasta_seq = "SGH3_HA_sequence.fasta")
cat(sprintf("  H3N2: %d sequences, %d amino acid positions, seasons %d-%d\n",
    nrow(h3_seq), ncol(h3_seq) - 4,
    min(as.numeric(h3_seq$season), na.rm = TRUE),
    max(as.numeric(h3_seq$season), na.rm = TRUE)))

# write.csv(h1_seq,file="toydata_result/SGH1_HA_sequence.csv",row.names = F)
# write.csv(h3_seq,file="toydata_result/SGH3_HA_sequence.csv",row.names = F)

# Output: HA and NA sequences with sequence name, accession number, collection date and flu season

cat("\n  Step 1 completed.\n\n")


###########################################################
# step 2: calculate site wise amino acid prevalence
# The function 'SitePrev' allows to input csv format sequence data and output site wise prevalence through time
# Run time: 3 minutes

cat("  Step 2: Calculating site-wise amino acid prevalence\n")
cat("  (This may take 3 minutes per serotype...)\n")

SitePrev <- function(input_csv_seq, label = ""){
  seq <- input_csv_seq
  year <- sort(unique(seq$season))
  year <- as.vector(year)
  n_sites <- length(seq[1,]) - 4
  codata <- matrix(nrow=length(year),ncol=0)
  for (i in 1:n_sites){
    ppdata <- NULL
    for (j in year){
      data <- seq[which(seq$season == j),]
      datat <- table(data[,i+4])
      dataf <- prop.table(datat)
      dataf <- as.data.frame(dataf)
      dataf <- t(dataf)
      colnames <- NULL
      for (k in 1:length(dataf[1,])){
        colnames <- c(colnames,paste(i,dataf[1,k],seq=""))
      }
      dimnames(dataf)=list(NULL,colnames)
      dataf <- dataf[-1,]
      dataf <- as.data.frame(t(dataf))
      ppdata <- rbind.fill(ppdata,dataf)
    }
    codata <- cbind(codata,ppdata)
  }
  codata[is.na(codata)] <- 0
  codata <- as.data.frame((lapply(codata,as.numeric)))
  rownames(codata)<-year
  return(codata)
}

cat("  [1/2] H1N1 prevalence ...\n")
h1_prev <- SitePrev(input_csv_seq = h1_seq, label = "H1N1:")
cat(sprintf("  H1N1 prevalence done: %d sites x %d seasons\n\n", ncol(h1_prev), nrow(h1_prev)))

cat("  [2/2] H3N2 prevalence ...\n")
h3_prev <- SitePrev(input_csv_seq = h3_seq, label = "H3N2:")
cat(sprintf("  H3N2 prevalence done: %d sites x %d seasons\n\n", ncol(h3_prev), nrow(h3_prev)))

# write.csv(h1_prev,file="toydata_result/SGH1_HA_prevalence.csv")
# write.csv(h3_prev,file="toydata_result/SGH3_HA_prevalence.csv")

# Output: HA and NA yearly site prevalence table

cat("  Step 2 completed.\n\n")


###########################################################
# step 3: calculate g-measure to quantify yearly mutational activity
# For detailed method, please see reference:
# Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026
# The function 'Getgmeature' allows to input site wise prevalence and output g-measure and average transition time
# Run time: 2 minutes

cat("  Step 3: Calculating g-measure (mutational activity)\n")
cat("  (This may take 2 minutes per serotype...)\n")

Getgmeature <- function(input_prev,theta_sg,h_sg){
  myh <- input_prev
  gmeasure <- matrix(ncol=0,nrow=length(myh[,1]))
  colname <- c()
  transition_time <- c()
  for (theta in theta_sg){
    for (h in h_sg){
      mut <- matrix(0,nrow=length(myh[,1]),ncol=length(myh))
      trans_time <- c()
      for (i in 1:length(myh[1,])){
        start <- 1
        r<-1
        while(r <= length(myh[,1])){
          # detect effective mutation
          if (myh[r,i]>=theta && any(myh[start:r,i]==0)){
            c <- r # record fisrt effective mutation site#
            start <- r+1 # prepare to detect next effective mutation#
            for (s in 1:r){
              if (myh[s,i]==0){
                a <- s           }
            }# record last zero#
            mut[(a+1):c,i] <- 1
            trans_time <- c(trans_time,c-a)
            # extend h years
            if (h!=0){
              if (c+h > length(myh[,1])){
                fakeh <- length(myh[,1])-c 
                if (fakeh!=0){
                  for (j in 1:fakeh){
                    if (myh[(c+j),i] >= theta){
                      mut[(c+j),i] <- 1       }else
                      {break}    
                  }
                }  
              }else{
                for (j in 1:h){
                  if (myh[(c+j),i] >= theta){
                    mut[(c+j),i] <-1}else
                    {break}}}
            }}
          r=r+1               }
      }
      transition_time <- c(transition_time,mean(trans_time))
      product <- c()
      for (y in 1:length(myh)){
        for(x in 1:length(myh[,1])){
          product <- c(product, myh[x,y]*mut[x,y])
        }
      }
      myy <- matrix(product, ncol=length(myh),nrow=length(myh[,1]),byrow=F)
      gsum <- rowSums(myy)
      gmeasure <- cbind(gmeasure,gsum)
      colname <- c(colname,paste0("theta=",theta,",h=",h))
    }
  }
  rownames(gmeasure)<-rownames(myh)
  colnames(gmeasure)<-colname
  return(gmeasure)
}

h1_gmeasure <- as.data.frame(Getgmeature(input_prev = h1_prev, theta_sg = 0.8, h_sg = 2))
cat(sprintf("  H1N1 g-measure done: %d seasons\n", nrow(h1_gmeasure)))

h3_gmeasure <- as.data.frame(Getgmeature(input_prev = h3_prev, theta_sg = 0.9, h_sg = 2))
cat(sprintf("  H3N2 g-measure done: %d seasons\n", nrow(h3_gmeasure)))

cat("  Merging g-measure results ...\n")

gmeasure <- full_join(
  rownames_to_column(h1_gmeasure, "season"),
  rownames_to_column(h3_gmeasure, "season"),
  by = "season") %>% 
  mutate(season = as.numeric(season)) %>%  
  arrange(season) %>% 
  column_to_rownames("season")
colnames(gmeasure) <- c("h1gmeasure","h3gmeasure")

# write.csv(gmeasure,file="toydata_result/SG_gmeasure.csv")

# Output: H1N1 and H3N2 g-measure by season

cat("\n  Step 3 completed.\n\n")


###########################################################
# step 4: Regression fitting associating serotype positivity rates and mutational activity
# Model: log(H3_rate / H1_rate) = beta0 + beta1 * log(H3_g / H1_g)

cat("  Step 4: Regression fitting associating serotype positivity rates and mutational activity\n")

sero <- read.csv("SG_seropositivity_rate.csv")
gmeasure <- read.csv("toydata_result/SG_gmeasure.csv")
colnames(gmeasure)[1] <- "season"
data <- merge(sero, gmeasure, by = "season")

# Calculate log odds:
#   logh3h1_p = log(Pr_H3 / Pr_H1) = log odds of positive rates
#   logh3h1_g = log(g_H3 / g_H1)   = log odds of g-measure

data$logh3h1_p <- log(data$h3n2_rate / data$h1n1pdm09_rate)
data$logh3h1_g <- log(data$h3gmeasure / data$h1gmeasure)
excl_idx <- which(!is.finite(data$logh3h1_p) | !is.finite(data$logh3h1_g))
if (length(excl_idx) > 0) {
  cat("Excluded seasons:", paste(data$season[excl_idx], collapse = ", "),
      "(non-finite values due to zero rates or NA g-measure)\n")
}
data_clean <- data[is.finite(data$logh3h1_p) & is.finite(data$logh3h1_g), ]
data_clean <- data_clean[which(data_clean$season<=2019),]

# Determine observed dominant serotype (for classification evaluation)
data_clean$dominant_obs <- ifelse(data_clean$h3n2_rate > data_clean$h1n1pdm09_rate, "h3", "h1")

# Linear regression: log(H3/H1 positivity) ~ log(H3/H1 g-measure)
fit <- lm(logh3h1_p ~ logh3h1_g, data = data_clean)
cat("---------- Model Summary ----------\n")
print(summary(fit))
cat("\nR-squared:", round(summary(fit)$r.squared, 4), "\n")
cat("p-value (logh3h1_g):", format(summary(fit)$coefficients[2, 4], digits = 4), "\n")
cat("Equation: log(H3/H1)_p =", round(coef(fit)[1], 4), "+",
    round(coef(fit)[2], 4), "* log(H3/H1)_g\n\n")

# Predicted values from the model
data_clean$pred_logh3h1_p <- predict(fit, data_clean)
data_clean$pred_dominant <- ifelse(data_clean$pred_logh3h1_p > 0, "h3", "h1")

# In-sample discrimination
cat("---------- In-sample Prediction ----------\n")
cat("Predicted vs Observed:\n")
print(table(Predicted = data_clean$pred_dominant, Observed = data_clean$dominant_obs))



###########################################################
# step 5 (optional): Scatter plot with regression line and 95% CI

pdf(file = "toydata_result/SG_regression_fit.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 2.5, 1.5))
plot(data_clean$logh3h1_g, data_clean$logh3h1_p,
     xlab = expression(log(italic(g)[H3N2] / italic(g)[H1N1])),
     ylab = expression(log(Pr[H3N2] / Pr[H1N1])),
     main = "Singapore: log odds of positive rate vs log odds of g-measure",
     pch = 19, col = ifelse(data_clean$dominant_obs == "h3", "#E41A1C", "#377EB8"),
     cex = 1.3, xlim = range(data_clean$logh3h1_g) * c(0.9, 1.1))

# Add season labels
text(data_clean$logh3h1_g, data_clean$logh3h1_p,
     labels = data_clean$season, pos = 3, cex = 0.6, col = "grey40")

# Regression line and 95% CI
x_seq <- seq(min(data_clean$logh3h1_g) - 0.5, max(data_clean$logh3h1_g) + 0.5, length.out = 100)
pred <- predict(fit, newdata = data.frame(logh3h1_g = x_seq), interval = "confidence", level = 0.95)
lines(x_seq, pred[, "fit"], col = "darkred", lwd = 2)
polygon(c(x_seq, rev(x_seq)), c(pred[, "lwr"], rev(pred[, "upr"])),
        col = rgb(0.8, 0, 0, 0.12), border = NA)

# Reference lines at 0
abline(h = 0, lty = 3, col = "grey60")
abline(v = 0, lty = 3, col = "grey60")

# Legend
eq_text <- paste0("y = ", round(coef(fit)[1], 3), " + ", round(coef(fit)[2], 3), " x")
r2_text <- paste0("R^2 == ", round(summary(fit)$r.squared, 3))
legend("topleft",
       legend = c(eq_text, parse(text = r2_text)),
       bty = "n", cex = 0.85)
legend("bottomright",
       legend = c("H3N2 dominant", "H1N1 dominant"),
       col = c("#E41A1C", "#377EB8"), pch = 19, bty = "n", cex = 0.8)
dev.off()
cat("Scatter plot saved to toydata_result/SG_regression_fit.pdf\n\n")

