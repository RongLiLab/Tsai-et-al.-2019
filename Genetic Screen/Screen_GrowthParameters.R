#required libraries
library(minpack.lm)
library(pracma)
library(openxlsx)

#List files in Directory "Screen_rawdata"
files <- list.files("~/Tsai_et_al./Screen/PrimaryScreen/RawData")

#Function to fit a sigmoid curve to each excel file/plate
apply_sigmoid <- function(x)
{
  data.OD <- read.xlsx(x, sheet = 1)
  row.names(data.OD) <- data.OD[,1]
  data.OD <- data.OD[,-1]
  time <- as.numeric(colnames(data.OD))
  Res <- apply(data.OD,1,function(z){
    AUC <- trapz(time,as.numeric(z))
    fit <- tryCatch(nlsLM( (z/max(z)) ~  1/(1 + exp(-k * (time-t0))) , start = list(k= 0.2, t0=20), control = nls.lm.control(maxiter = 200))
                    ,error=function(fit) {list(asymptote=NaN, k=NaN, t0=NaN)})
    c(coef(fit), AUC)
  })
  Final <- t(Res)
  colnames(Final) <- c("Maximum Growth Rate", "Lag Phase", "AUC")
  write.xlsx(Final,paste0("Results_",x), col.names = TRUE, row.names = TRUE)
}

#Run over all plates
lapply(files, apply_sigmoid)