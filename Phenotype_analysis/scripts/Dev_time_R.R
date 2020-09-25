# Set working directory
setwd("/Users/niclasbackstrom/Documents/Work/Projects/Leptidea/Host_plants/Results/")

# Read in the data
df <- read.delim("2018_dev_time.txt", 
               header = TRUE, 
               quote="\"", 
               stringsAsFactors= TRUE, 
               strip.white = TRUE)

# Print out `df`
print(df)

# Check structure of data
names(df)
str(df)
# Check summary of data
summary(df)

# Categorize data sets and make box plot (Sweden)
boxplot(df$Dev_time[df[,1] == "Swe" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Swe" & df[,3] == "Lotus"],
		col = rgb(0.8,0.1,0.3,0.6),
		names=c(expression(italic(Lotus_dorycnium)),expression(italic(Lotus_corniculatus))),
		ylab = "Development time (days)",
		main = "Sweden",
		cex.lab = 1.5,
		ylim=c(20,45))
		
# Categorize data sets and test for difference (Sweden)
wilcox.test(df$Dev_time[df[,1] == "Swe" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Swe" & df[,3] == "Lotus"])
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Dor"] and df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Lotus"]
# W = 5967.5, p-value = 5.819e-07
# alternative hypothesis: true location shift is not equal to 0
#################################################################################################################

# Categorize data sets and make box plot (Catalonia)
boxplot(df$Dev_time[df[,1] == "Cat" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Lotus"],
		col = rgb(0.1,0.1,0.7,0.5),
		names=c(expression(italic(Lotus_dorycnium)),expression(italic(Lotus_corniculatus))),
		ylab = "Development time (days)",
		main = "Catalonia",
		cex.lab = 1.5,
		ylim=c(20,45))
		
# Categorize data sets and test for difference (Catalonia)
wilcox.test(df$Dev_time[df[,1] == "Cat" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Lotus"])
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Dor"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Lotus"]
# W = 21425, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
#################################################################################################################		

# Categorize data sets and test for difference between populations on both host plants
wilcox.test(df$Dev_time[df[,1] == "Swe" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Dor"])
######################Results##########################Results##########################Results##################
#	Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Dor"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Dor"]
# W = 1541.5, p-value = 1.161e-08
# alternative hypothesis: true location shift is not equal to 0
#################################################################################################################	

wilcox.test(df$Dev_time[df[,1] == "Swe" & df[,3] == "Lotus"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Lotus"])

######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Lotus"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Lotus"]
# W = 7909.5, p-value = 1.273e-11
# alternative hypothesis: true location shift is not equal to 0
#################################################################################################################



# SUMMARIES OF COUNTS ETC


> head(df)
  Population  ID Host_plant Sex Dev_time  X
1        Cat  60      Lotus   F       30 NA
2        Cat  88      Lotus   F       31 NA
3        Cat 164        Dor   M       35 NA
4        Cat 165        Dor   M       35 NA
5        Cat 252        Dor   F       40 NA
6        Cat   3      Lotus   M       26 NA

> table(df$Population,df$Host_plant)
     
      Dor Lotus
  Cat 159   151
  Swe  44   183

> table(df$Population,df$Host_plant,df$Sex)
, ,  = F

     
      Dor Lotus
  Cat  79    64
  Swe  24    95

, ,  = M

     
      Dor Lotus
  Cat  80    87
  Swe  20    88

> tapply(df$Dev_time,list(df$Sex,df$Host_plant),mean,na.rm=T)
       Dor    Lotus
F 33.69903 29.44654
M 33.10000 28.16571

> tapply(df$Dev_time,list(df$Population,df$Host_plant),mean,na.rm=T)
         Dor    Lotus
Cat 34.13836 29.83444
Swe 30.75000 27.90164

> tapply(df$Dev_time,list(df$Population,df$Host_plant),sd,na.rm=T)
         Dor    Lotus
Cat 3.063625 2.423305
Swe 3.376561 2.557288
