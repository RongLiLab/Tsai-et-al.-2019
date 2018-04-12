#####TECAN batch processing code
options(stringsAsFactors=FALSE) # Don't covert strings to factors

##Load packages
library(openxlsx)
library(ggplot2)
library(grofit)
library(reshape2)
library(growthcurve)

##List files in folder 
files <- list.files()

##Read data from all files 
raw.data <- lapply(files, function(x){ read.xlsx(x,  sheet = 1, colNames = FALSE) }) 
names(raw.data) <- files

##Save Plots for all files as pdfs & data as excel workbooks

TECAN_plot <- function(y)
  {
#name columns
colnames(y) <- paste0(rep(LETTERS[1:8], each = 12), rep(1:12,8))

#Arranged data
y$Time <- 0.5*seq(0,(length(y$A1)-1),1)
arrange.data <- melt(y, id.vars = "Time")
colnames(arrange.data) <- c("Time","Well","OD") 
arrange.data$Row <- gsub("[[:digit:]]", "", arrange.data$Well)
arrange.data$Column <- gsub("[^[:digit:]]", "", arrange.data$Well)

#arrange wells in correct order
arrange.data$Well <- factor(arrange.data$Well, levels = colnames(y)[-97])
arrange.data$Column <- factor(arrange.data$Column, levels = as.character(1:12))

#Result plot
g <- ggplot(arrange.data, aes(x = Time, y = OD)) + geom_point(color = "black", size = 0.5, shape = 1) + 
  stat_growthcurve(color = "red", model = "spline", show.legend = FALSE) + facet_grid(Row ~ Column) + theme_bw() +
  theme(panel.grid = element_blank())

#Split by well
well.data <- split(arrange.data, list(arrange.data$Well))

#Get Spline fit & growth curve parameters

#Parameters
get_spline <- function(w)
{
  x <- w$OD  
  t <- w$Time
  well <- as.character(w$Well[1])
  run <- gcFitSpline(t, x,gcID = "undefined", control = grofit.control())
  L <- summary(run)
  total <- matrix(c(well,unlist(summary(run), use.names = FALSE)), nrow = 1, ncol = 5, byrow = TRUE)
  total
}

#Result table

Res <- lapply(well.data, get_spline)
Res <- do.call("rbind", Res)
colnames(Res) <- c("Well","Growth Rate","Lag Phase","Maximum Growth","AUC")
Res <- data.frame(Res)

list(parameters = Res, plots = g)

}

All_Results <- lapply(raw.data,TECAN_plot)

#Save Plots & results
for(i in 1:length(All_Results))
{
  
  title_plot <- strsplit(files[i],".xlsx")[1]
 p <- All_Results[[i]]$plots + ggtitle(title_plot)
 p
 ggsave(paste0(title_plot,".pdf"), width = 12, height = 7, units = "in")
write.xlsx(All_Results[[i]]$parameters, paste0("Result_", files[i]))
}

