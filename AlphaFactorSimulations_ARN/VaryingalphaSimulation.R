#Load libraries
library(ggplot2)
library(openxlsx)
#Read the data produced by Alpha.py
data.sim <- read.csv("base.csv", header = TRUE)
data.sim$Run <- 0
data.sim$EqDiameter <- 8.6*data.sim$Volume^(0.3) #This was from the nls power law fit to th euploid
alpha.iter <- vector(mode = "list", length = 100)

#Vary alpha raandomly over 0.7 to 1 for each simulation. Repeat the simulation 100 times.
for(i in 1:100)
{
  alpha <- sample(seq(0.7,1,0.05),20, replace = TRUE)
  volume <- data.sim$Volume*alpha
  volume[1] <- 1
  volume[20] <- 2
  sd <- data.sim$SD*alpha
   ploidy <- data.sim$Ploidy
  alpha.iter[[i]] <- data.frame(Run = i, Ploidy = ploidy , Volume = volume, SD = sd)
}

#Combine the result of all the simulations
alpha.iter <- do.call("rbind", alpha.iter)
alpha.iter$EqDiameter <- 8.6*alpha.iter$Volume^(0.3)

#Plot data
alpha.plot.data <- rbind(data.sim,alpha.iter)
alpha.plot.data$EqDiameter <- 8.6*alpha.plot.data$Volume^(0.3)

#Plots from individual simulations
ggplot(alpha.plot.data, aes(x = Ploidy, y = EqDiameter)) + geom_line(color = "red") + 
 theme_bw() + facet_wrap(~Run,nrow = 10) + 
 theme(panel.grid = element_blank(), strip.background = element_rect(fill = NA))

#Average of all simulations compared to constant alpha = 0.7 and alpha = 1.0
ggplot() + geom_smooth(data = alpha.iter,aes(x = Ploidy, y = EqDiameter), method = "loess", span = 0.2, color = "red") + 
  geom_line(data = data.sim,aes(x = Ploidy, y = EqDiameter), color = "darkblue", linetype = "dashed") + 
  geom_line(data = data.sim,aes(x = Ploidy, y = 0.7^(1/3)*EqDiameter),color = "steelblue", linetype = "dashed" ) + 
  theme_classic() +
  scale_color_manual(values = c("red"="red","darkblue"="darkblue","steelblue"="steelblue"), labels = c("varying alpha", "alpha = 1.0", "alpha = 0.7"))

write.xlsx(alpha.plot.data,"VaryingAlpha_Simlation.xlsx")
                          
