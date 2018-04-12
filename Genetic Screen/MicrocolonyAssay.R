#load required libraries
library(openxlsx)
library(ggplot2)



#list all raw microcolony data files the in folder

strains <- list.files("~/Tsai_et_al./MicrocolonyAssay/RawData_Round1")

#Read data

strains_eu <- lapply(files, function(x){ read.xlsx(x, sheet = 2)})
  
strains_aneu <- lapply(files, function(x){ read.xlsx(x, sheet = 1)})


#Compute growth rate for every microcolony in each well

growth_rate <- function(x)
{
  y <- log2(na.exclude(x))
  t <- 0:(length(y)-1)
  fit <- lm(y ~t)
  rsquare <- cor(t,y)
  res <- cbind(coef(fit)[2], rsquare)
  res
}


strain_rates <- function(z){
  gr <- t(apply(z,2, growth_rate))
  wells <- t(data.frame(strsplit(row.names(gr), split = "_")))
  row.names(gr) = row.names(wells) = NULL
  res2 <- data.frame(wells,gr)
  colnames(res2) <- c("Gene", "Well", "Microcolony", "Growth Rate", "Rsquared")
}


#Apply to all strains

Growth_euploids <- lapply(strains_eu, strains_rates)
Growth_euploids <- do.call("rbind", Growth_euploids)

Growth_aneuploids <- lapply(strains_aneu, strain_rates)
Growth_aneuploids <- do.call("rbind", Growth_aneuploids)

Results <- rbind(Growth_euploids, Growth_aneuploids)

#Save results to excel file

write.xlsx(Results, "Results_Microcolony_R1.xlsx")




