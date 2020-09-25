# Set working directory
setwd("/Users/niclasbackstrom/Documents/Work/Projects/Leptidea/Research_projects/Host_plant/Results/Wing_size/Analysis/R")

# Read in the data
df <- read.delim("Input_table.txt", 
               header = TRUE, 
               quote="\"", 
               stringsAsFactors= TRUE, 
               strip.white = TRUE)

# Print out `df`
print(df)

# Check structure of data
names(df)
str(df)

######################Results##########################Results##########################Results##################
# 'data.frame':	537 obs. of  11 variables:
# $ Family      : Factor w/ 48 levels "C1","C12","C13",..: 19 23 24 24 24 24 28 1 28 19 ...
# $ ID          : int  3 4 5 6 7 8 9 10 11 12 ...
# $ Photo       : Factor w/ 297 levels "1","10","100",..: 1 112 218 229 240 251 262 273 284 2 ...
# $ Population  : Factor w/ 2 levels "Cat","Swe": 1 1 1 1 1 1 1 1 1 1 ...
# $ Wing_length : int  1273 1262 1268 1241 1219 1210 1214 1096 1187 1256 ...
# $ Wing_size   : int  35255 34386 34448 34094 32133 34167 32302 26934 33405 35381 ...
# $ Host_plant  : Factor w/ 2 levels "Dor","Lotus": 2 2 2 2 2 2 2 1 1 2 ...
# $ Egg_date    : Factor w/ 19 levels "02-May","04-May",..: 17 17 17 17 17 17 14 14 14 17 ...
# $ Adult_eclose: Factor w/ 33 levels "01-Jun","02-Jun",..: 24 24 24 24 24 24 24 24 23 26 ...
# $ Sex         : Factor w/ 2 levels "F","M": 2 2 2 2 2 1 2 2 1 2 ...
# $ Dev_Time    : int  26 26 26 26 26 26 30 29 30 27 ...
######################Results##########################Results##########################Results################## 

# Check summary of data
summary(df)

######################Results##########################Results##########################Results##################
#     Family          ID             Photo     Population  Wing_length     Wing_size     Host_plant     Egg_date    Adult_eclose Sex        Dev_Time    
# SK..8  : 34   Min.   :   3.0   broken : 23   Cat:310    Min.   : 881   Min.   :16596   Dor  :203   08-May :100   30-May : 37   F:262   Min.   :21.00  
# SK..11 : 30   1st Qu.: 137.0   missing:  4   Swe:227    1st Qu.:1129   1st Qu.:28798   Lotus:334   16-May : 72   28-May : 32   M:275   1st Qu.:28.00  
# C34    : 28   Median : 271.0   1      :  2              Median :1187   Median :32041               28-Apr : 69   16-Jun : 30           Median :30.00  
# SK..6  : 25   Mean   : 561.8   10     :  2              Mean   :1174   Mean   :31704               05-May : 57   10-Jun : 27           Mean   :30.53  
# SK..7  : 25   3rd Qu.:1093.0   100    :  2              3rd Qu.:1235   3rd Qu.:35317               20-May : 55   15-Jun : 25           3rd Qu.:33.00  
# C30    : 23   Max.   :1227.0   101    :  2              Max.   :1331   Max.   :42401               25-Apr : 37   31-May : 25           Max.   :45.00  
# (Other):372                    (Other):502              NA's   :28     NA's   :28                  (Other):147   (Other):361                          
######################Results##########################Results##########################Results##################

# Categorize data sets and make box plot (Sweden)
boxplot(df$Wing_length[df[,4] == "Swe" & df[,7] == "Dor"], 
		df$Wing_length[df[,4] == "Swe" & df[,7] == "Lotus"],
		col = rgb(0.8,0.1,0.3,0.6),
		names=c(expression(italic('Lotus dorycnium')),expression(italic('Lotus corniculatus'))),
		ylab = "Wing lenght (pixels)",
		main = "Sweden",
		cex.lab = 1.5,
		ylim=c(900,1350))
		
# Categorize data sets and test for difference (Sweden)
wilcox.test(df$Wing_length[df[,4] == "Swe" & df[,7] == "Dor"], 
		df$Wing_length[df[,4] == "Swe" & df[,7] == "Lotus"])
		
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Wing_length[df[, 4] == "Swe" & df[, 7] == "Dor"] and df$Wing_length[df[, 4] == "Swe" & df[, 7] == "Lotus"]
# W = 1403.5, p-value = 5.125e-10
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################

# Categorize data sets and make box plot (Catalonia)
boxplot(df$Wing_length[df[,4] == "Cat" & df[,7] == "Dor"], 
		df$Wing_length[df[,4] == "Cat" & df[,7] == "Lotus"],
		col = rgb(0.1,0.1,0.7,0.5),
		names=c(expression(italic('Lotus dorycnium')),expression(italic('Lotus corniculatus'))),
		ylab = "Wing lenght (pixels)",
		main = "Catalonia",
		cex.lab = 1.5,
		ylim=c(900,1350))
		
# Categorize data sets and test for difference (Catalonia)
wilcox.test(df$Wing_length[df[,4] == "Cat" & df[,7] == "Dor"], 
		df$Wing_length[df[,4] == "Cat" & df[,7] == "Lotus"])
		
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Dor"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Lotus"]
# W = 21425, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################		

# Categorize data sets and make box plots of wing areas (Sweden)
boxplot(df$Wing_size[df[,4] == "Swe" & df[,7] == "Dor"], 
		df$Wing_size[df[,4] == "Swe" & df[,7] == "Lotus"],
		col = rgb(0.8,0.1,0.3,0.6),
		names=c(expression(italic('Lotus dorycnium')),expression(italic('Lotus corniculatus'))),
		ylab = "Wing area (pixels)",
		main = "Sweden",
		cex.lab = 1.5,
		ylim=c(15000,45000))
		
# Categorize data sets and test for difference (Sweden)
wilcox.test(df$Wing_size[df[,4] == "Swe" & df[,7] == "Dor"], 
		df$Wing_size[df[,4] == "Swe" & df[,7] == "Lotus"])
		
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
# 
# data:  df$Wing_size[df[, 4] == "Swe" & df[, 7] == "Dor"] and df$Wing_size[df[, 4] == "Swe" & df[, 7] == "Lotus"]
# W = 1385, p-value = 3.715e-10
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################

# Categorize data sets and make box plots of wing areas (Catalonia)
boxplot(df$Wing_size[df[,4] == "Cat" & df[,7] == "Dor"], 
		df$Wing_size[df[,4] == "Cat" & df[,7] == "Lotus"],
		col = rgb(0.1,0.1,0.7,0.5),
		names=c(expression(italic('Lotus dorycnium')),expression(italic('Lotus corniculatus'))),
		ylab = "Wing area (pixels)",
		main = "Catalonia",
		cex.lab = 1.5,
		ylim=c(15000,45000))
		
# Categorize data sets and test for difference (Catalonia)
wilcox.test(df$Wing_size[df[,4] == "Cat" & df[,7] == "Dor"], 
		df$Wing_size[df[,4] == "Cat" & df[,7] == "Lotus"])
		
######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
# 
# data:  df$Wing_size[df[, 4] == "Cat" & df[, 7] == "Dor"] and df$Wing_size[df[, 4] == "Cat" & df[, 7] == "Lotus"]
# W = 2412.5, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################


# Categorize data sets and test for difference between populations on both host plants
wilcox.test(df$Dev_time[df[,1] == "Swe" & df[,3] == "Dor"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Dor"])
		
######################Results##########################Results##########################Results##################
#	Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Dor"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Dor"]
# W = 1541.5, p-value = 1.161e-08
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################	

wilcox.test(df$Dev_time[df[,1] == "Swe" & df[,3] == "Lotus"], 
		df$Dev_time[df[,1] == "Cat" & df[,3] == "Lotus"])

######################Results##########################Results##########################Results##################
# Wilcoxon rank sum test with continuity correction
#
# data:  df$Dev_time[df[, 1] == "Swe" & df[, 3] == "Lotus"] and df$Dev_time[df[, 1] == "Cat" & df[, 3] == "Lotus"]
# W = 7909.5, p-value = 1.273e-11
# alternative hypothesis: true location shift is not equal to 0
######################Results##########################Results##########################Results##################


# SUMMARIES OF COUNTS ETC

table(df$Population,df$Host_plant,df$Sex)

######################Results##########################Results##########################Results##################
# , ,  = F 
#       Dor Lotus
#   Cat  79    64
#   Swe  24    95
# , ,  = M
#       Dor Lotus
#   Cat  80    87
#   Swe  20    88
######################Results##########################Results##########################Results##################
 
table(df$Population[df$Wing_length!="NA"],df$Host_plant[df$Wing_length!="NA"],df$Sex[df$Wing_length!="NA"])

######################Results##########################Results##########################Results##################
# , ,  = F
#       Dor Lotus
#   Cat  76    60
#   Swe  23    91
# , ,  = M
#       Dor Lotus
#   Cat  77    79
#   Swe  19    84
######################Results##########################Results##########################Results##################

> head(df)

######################Results##########################Results##########################Results##################
#  Population  ID Host_plant Sex Dev_time  X
# 1        Cat  60      Lotus   F       30 NA
# 2        Cat  88      Lotus   F       31 NA
# 3        Cat 164        Dor   M       35 NA
# 4        Cat 165        Dor   M       35 NA
# 5        Cat 252        Dor   F       40 NA
# 6        Cat   3      Lotus   M       26 NA
######################Results##########################Results##########################Results##################

table(df$Population,df$Host_plant)

######################Results##########################Results##########################Results##################    
#       Dor Lotus
#   Cat 159   151
#   Swe  44   183
######################Results##########################Results##########################Results##################

tapply(df$Wing_length,list(df$Population,df$Host_plant),mean,na.rm=T)

######################Results##########################Results##########################Results##################    
#         Dor    Lotus
# Cat 1101.647 1214.554
# Swe 1115.095 1218.200
######################Results##########################Results##########################Results##################   

tapply(df$Wing_length,list(df$Population,df$Host_plant),sd,na.rm=T)

######################Results##########################Results##########################Results##################    
#         Dor    Lotus
# Cat 70.05867 54.41523
# Swe 93.99217 50.75114
######################Results##########################Results##########################Results##################  

tapply(df$Wing_size,list(df$Population,df$Host_plant),mean,na.rm=T)

######################Results##########################Results##########################Results##################    
#         Dor    Lotus
# Cat 27560.53 33383.11
# Swe 28738.17 34703.96
######################Results##########################Results##########################Results################## 

tapply(df$Wing_size,list(df$Population,df$Host_plant),sd,na.rm=T)

######################Results##########################Results##########################Results##################    
#          Dor    Lotus
# Cat 3697.864 3413.433
# Swe 5428.048 3167.473
######################Results##########################Results##########################Results################## 






tapply(df$Wing_length,list(df$Population,df$Sex,df$Host_plant),mean,na.rm=T)

######################Results##########################Results##########################Results##################
# , , Dor
#            F        M
# Cat 1105.526 1097.818
# Swe 1137.870 1087.526
# , , Lotus
#            F        M
# Cat 1225.900 1205.937
# Swe 1226.692 1209.000
######################Results##########################Results##########################Results##################

# Categorize data sets and check difference in median wing lengths
tapply(df$Wing_length,list(df$Population,df$Sex),median,na.rm=T)

######################Results##########################Results##########################Results##################
#       F      M
# Cat 1161 1158.5
# Swe 1229 1196.0
######################Results##########################Results##########################Results##################

# Categorize data sets and check difference in mean wing lengths
tapply(df$Wing_length,list(df$Population,df$Sex),mean,na.rm=T)

######################Results##########################Results##########################Results##################
#           F        M
# Cat 1158.632 1152.571
# Swe 1208.772 1186.592
######################Results##########################Results##########################Results##################