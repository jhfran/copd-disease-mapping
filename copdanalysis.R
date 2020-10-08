##### 
#Loading packages
library(tidyverse)
library(spdep)
library(sf)
library(CARBayes)
library(rgdal)
library(rgeos)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(knitr)
#####
#Reading in data 
COPDexpected <- read_csv("copdexpected.csv")
COPDobserved <- read_csv("copdobserved.csv")
COPDsf <- readOGR(dsn = '.',
                  layer = 'englandlocalauthority') #reading in shapefiles
#####
#Merging datasets 
COPD <- merge(COPDobserved, COPDexpected, by = 'Name')
#Renaming columns
COPD<- rename(COPD, "name" = "Name")
COPD <- rename(COPD, "O2001" = "Y2001", 
               "O2002" = "Y2002",
               "O2003" = "Y2003",
               "O2004" = "Y2004",
               "O2005" = "Y2005",
               "O2006" = "Y2006",
               "O2007" = "Y2007",
               "O2008" = "Y2008",
               "O2009" = "Y2009",
               "O2010" = "Y2010")
#####
##Question 1

#Summarising the number of hospital admissions 
COPDsummary<- summary(COPD)
Means <- c(mean(COPD$O2001), mean(COPD$O2002),
           mean(COPD$O2003), mean(COPD$O2004),
           mean(COPD$O2005), mean(COPD$O2006),
           mean(COPD$O2007), mean(COPD$O2008),
           mean(COPD$O2009), mean(COPD$O2010))
mdf <- data.frame(Year = 2001:2010, ObservedMeanHospitalAdmissions = Means)
kable (mdf, digits = 2, caption = "Table 1: Mean of Observed Hospital Admissions for COPD in England between 2001 - 2010 ")
#Plot Hospital admissions 
ggplot(data = mdf, aes(x = Year, y = ObservedMeanHospitalAdmissions))+
  geom_line(color = "blue") +
  ylab("Observed Mean Hospital Admissions for COPD") +
  xlab("Year") +
  theme_minimal() +
  ggtitle("Mean Observed Hospital Admissions for COPD", subtitle = "England, Years 2001 - 2010") +
  scale_x_continuous(breaks = c(2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010))
#####
##2A - Raw SMR 2001

#Calculating the raw SMRS for the Year 2001
COPD$SMR2001raw <- COPD$O2001/COPD$E2001
#Merging Raw SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Raw SMRs
raw2001plot <- ggplot(SMR,
                      aes(fill = SMR2001raw)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))
#####
##2B Smooth SMR 2001

# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2001 <- S.CARleroux(formula = O2001 ~ offset(log(E2001)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2001
COPD$FittedValue2001 <- smodel2001$fitted.values
COPD$SMR2001smooth <- COPD$FittedValue2001/COPD$E2001
#Summarising raw SMRs
summary(COPD$SMR2001raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2001smooth)
#Merging Smooth SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Smooth SMRs
smooth2001plot <- ggplot(SMR,
                         aes(fill = SMR2001smooth)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))
# Plot of raw vs smoothed SMRs
rvssm2001 <- ggplot(COPD, aes(SMR2001raw,SMR2001smooth)) +
  geom_point(colour = 'blue') +
  # Adding x = y line for comparison
  geom_abline(intercept = 0,
              slope = 1,
              colour = 'red') +
  labs(x = 'Raw SMR',
       y = 'Smooth SMR') +
  theme_minimal()
##### 
##Raw SMR 2003
#Calculating the raw SMRS for the Year 2003
COPD$SMR2003raw <- COPD$O2003/COPD$E2003
#Merging Raw SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Raw SMRs
raw2003plot <- ggplot(SMR,
                      aes(fill = SMR2003raw)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Raw SMRS for England in 2003")
#####
#2B Smooth SMR 2003
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2003 <- S.CARleroux(formula = O2003 ~ offset(log(E2003)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2003
COPD$FittedValue2003 <- smodel2003$fitted.values
COPD$SMR2003smooth <- COPD$FittedValue2003/COPD$E2003
#Summarising raw SMRs
summary(COPD$SMR2003raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2003smooth)
#Merging Smooth SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Smooth SMRs
smooth2003plot <- ggplot(SMR,
                         aes(fill = SMR2003smooth)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Smooth SMRS for England in 2003")
# Plot of raw vs smoothed SMRs
rvssm2003 <- ggplot(COPD, aes(SMR2003raw,SMR2003smooth)) +
  geom_point(colour = 'blue') +
  # Adding x = y line for comparison
  geom_abline(intercept = 0,
              slope = 1,
              colour = 'red') +
  labs(x = 'Raw SMR',
       y = 'Smooth SMR') +
  theme_minimal()+
  ggtitle("Raw vs Smooth SMRs 2003")

#####
##Raw SMR 2008
#Calculating the raw SMRS for the Year 2008
COPD$SMR2008raw <- COPD$O2008/COPD$E2008
#Merging Raw SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Raw SMRs
raw2008plot <- ggplot(SMR,
                      aes(fill = SMR2008raw)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Raw SMRS for England in 2008")
#Smooth SMR 2008
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2008 <- S.CARleroux(formula = O2008 ~ offset(log(E2008)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2008
COPD$FittedValue2008 <- smodel2008$fitted.values
COPD$SMR2008smooth <- COPD$FittedValue2008/COPD$E2008
#Summarising raw SMRs
summary(COPD$SMR2008raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2008smooth)
#Merging Smooth SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Smooth SMRs
smooth2008plot <- ggplot(SMR,
                         aes(fill = SMR2008smooth)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Smooth SMRS for England in 2008")
# Plot of raw vs smoothed SMRs
rvssm2008 <- ggplot(COPD, aes(SMR2008raw,SMR2008smooth)) +
  geom_point(colour = 'blue') +
  # Adding x = y line for comparison
  geom_abline(intercept = 0,
              slope = 1,
              colour = 'red') +
  labs(x = 'Raw SMR',
       y = 'Smooth SMR') +
  theme_minimal()+
  ggtitle("Raw vs Smooth SMRs 2008")
#####
##Raw SMR 2010
#Calculating the raw SMRS for the Year 2010
COPD$SMR2010raw <- COPD$O2010/COPD$E2010
#Merging Raw SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Raw SMRs
raw2010plot <- ggplot(SMR,
                      aes(fill = SMR2010raw)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Raw SMRS for England in 2010")
#Smooth SMR 2010
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2010 <- S.CARleroux(formula = O2010 ~ offset(log(E2010)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2010
COPD$FittedValue2010 <- smodel2010$fitted.values
COPD$SMR2010smooth <- COPD$FittedValue2010/COPD$E2010
#Summarising raw SMRs
summary(COPD$SMR2010raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2010smooth)
#Merging Smooth SMRs and the shapefile
SMR <- merge(COPDsf, COPD, by = 'name')
# Converting shapefile into something ggplot can work with
SMR <- st_as_sf(SMR)
# Creating map of Smooth SMRs
smooth2010plot <- ggplot(SMR,
                         aes(fill = SMR2010smooth)) +
  geom_sf(colour = NA) +
  theme_minimal() +
  labs(x = 'Longitude',
       y = 'Latitude',
       fill = 'SMR') +
  scale_fill_gradientn(colours = brewer.pal(9, 'PuBu'),
                       breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.5))+
  ggtitle("Smooth SMRS for England in 2010")
# Plot of raw vs smoothed SMRs
rvssm2010 <- ggplot(COPD, aes(SMR2010raw,SMR2010smooth)) +
  geom_point(colour = 'blue') +
  # Adding x = y line for comparison
  geom_abline(intercept = 0,
              slope = 1,
              colour = 'red') +
  labs(x = 'Raw SMR',
       y = 'Smooth SMR') +
  theme_minimal()+
  ggtitle("Raw vs Smooth SMRs 2010")

#####
##Raw SMR 2002
#Calculating the raw SMRS for the Year 2002
COPD$SMR2002raw <- COPD$O2002/COPD$E2002
#Smooth SMR 2002
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2002 <- S.CARleroux(formula = O2002 ~ offset(log(E2002)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2002
COPD$FittedValue2002 <- smodel2002$fitted.values
COPD$SMR2002smooth <- COPD$FittedValue2002/COPD$E2002
#Summarising raw SMRs
summary(COPD$SMR2002raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2002smooth)
#####
#Calculating SMRs for 2004
##Raw SMR 2004
#Calculating the raw SMRS for the Year 2004
COPD$SMR2004raw <- COPD$O2004/COPD$E2004
#Smooth SMR 2004
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2004 <- S.CARleroux(formula = O2004 ~ offset(log(E2004)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2004
COPD$FittedValue2004 <- smodel2004$fitted.values
COPD$SMR2004smooth <- COPD$FittedValue2004/COPD$E2004
#Summarising raw SMRs
summary(COPD$SMR2004raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2004smooth)
#####
#Calculating SMRs for 2005
##Raw SMR 2005
#Calculating the raw SMRS for the Year 2005
COPD$SMR2005raw <- COPD$O2005/COPD$E2005
#Smooth SMR 2005
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
#Calculating SMRs for 2005
smodel2005 <- S.CARleroux(formula = O2005 ~ offset(log(E2005)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2005
COPD$FittedValue2005 <- smodel2005$fitted.values
COPD$SMR2005smooth <- COPD$FittedValue2005/COPD$E2005
#Summarising raw SMRs
summary(COPD$SMR2005raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2005smooth)
#####
#Calculating SMRs for 2006
##Raw SMR 2006
#Calculating the raw SMRS for the Year 2006
COPD$SMR2006raw <- COPD$O2006/COPD$E2006
#Smooth SMR 2006
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2006 <- S.CARleroux(formula = O2006 ~ offset(log(E2006)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2006
COPD$FittedValue2006 <- smodel2006$fitted.values
COPD$SMR2006smooth <- COPD$FittedValue2006/COPD$E2006
#Summarising raw SMRs
summary(COPD$SMR2006raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2006smooth)
#####
#Calculating SMRs for 2007
##Raw SMR 2007
#Calculating the raw SMRS for the Year 2007
COPD$SMR2007raw <- COPD$O2007/COPD$E2007
#Smooth SMR 2007
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2007 <- S.CARleroux(formula = O2007 ~ offset(log(E2007)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2007
COPD$FittedValue2007 <- smodel2007$fitted.values
COPD$SMR2007smooth <- COPD$FittedValue2007/COPD$E2007
#Summarising raw SMRs
summary(COPD$SMR2007raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2007smooth)
#####
#Calculating SMRs for 2009
##Raw SMR 2009
#Calculating the raw SMRS for the Year 2009
COPD$SMR2009raw <- COPD$O2009/COPD$E2009
#Smooth SMR 2009
# Creates the neighbourhood
W.nb <- poly2nb(COPDsf, row.names = rownames(COPDsf))
# Creates a matrix for following function call
W.mat <- nb2mat(W.nb, style="B")
#Running a smooth model
smodel2009 <- S.CARleroux(formula = O2009 ~ offset(log(E2009)), # Model Formula
                          data = COPD, # Dataset name
                          family = "poisson", # Choosing Poisson Regression
                          W = W.mat, # Neighbourhood matrix
                          burnin = 20000, # Number of burn in samples (throw away)
                          n.sample = 100000, # Number of MCMC sample
                          thin = 10,
                          rho = 1)
# Creating a dataset with smoothed SMRs in 2009
COPD$FittedValue2009 <- smodel2009$fitted.values
COPD$SMR2009smooth <- COPD$FittedValue2009/COPD$E2009
#Summarising raw SMRs
summary(COPD$SMR2009raw) 
# Summarising smoothed SMRs
summary(COPD$SMR2009smooth)
#####
#Calculating the smooth mean over 10 years
COPD$SMRSmoothMean<-(COPD$SMR2001smooth+COPD$SMR2002smooth+COPD$SMR2003smooth+COPD$SMR2004smooth+COPD$SMR2005smooth+COPD$SMR2006smooth+COPD$SMR2007smooth+COPD$SMR2008smooth+COPD$SMR2009smooth+COPD$SMR2010smooth)/10
#Isalating mean smooth SMR for 10 years
COPDmeanSMR <- subset(COPD, select = c("name", "SMRSmoothMean"))
#Mean SMR for 10 years with SMR 1.5
above15<- grepl("^1.5", COPDmeanSMR$SMRSmoothMean)
above15<- COPDmeanSMR[above15, ]
#Mean SMR for 10 years with SMR 1.6
above16<- grepl("^1.6", COPDmeanSMR$SMRSmoothMean)
above16<- COPDmeanSMR[above16, ]
#Mean SMR for 10 years with SMR 1.7
above17<- grepl("^1.7", COPDmeanSMR$SMRSmoothMean)
above17<- COPDmeanSMR[above17, ]
#Mean SMR for 10 years with SMR 1.8
above18<- grepl("^1.8", COPDmeanSMR$SMRSmoothMean)
above18<- COPDmeanSMR[above18, ]
#none
#Mean SMR for 10 years with SMR 1.9
above19<- grepl("^1.9", COPDmeanSMR$SMRSmoothMean)
above19<- COPDmeanSMR[above19, ]
#none
#Mean SMR for 10 years with SMR 2
above2<- grepl("^2", COPDmeanSMR$SMRSmoothMean)
#none 
above2<- COPDmeanSMR[above2, ]
#Creating tables for boroughs with SMRs over 1.5
kable(above15, digits = 2)
kable(above16, digits = 2)
kable(above17, digits = 2)

#####
MeanSmoothSMR <- c(mean(COPD$SMR2001smooth), mean(COPD$SMR2002smooth),
           mean(COPD$SMR2003smooth), mean(COPD$SMR2004smooth),
           mean(COPD$SMR2005smooth), mean(COPD$SMR2006smooth),
           mean(COPD$SMR2007smooth), mean(COPD$SMR2008smooth),
           mean(COPD$SMR2009smooth), mean(COPD$SMR2010smooth))
sdf <- data.frame(Year = 2001:2010, MeanSmoothSMR = MeanSmoothSMR)

#Plot Mean Smooth SMR
ggplot(data = sdf, aes(x = Year, y = MeanSmoothSMR))+
  geom_line(color = "blue") +
  ylab("Mean Smooth SMR") +
  xlab("Year") +
  theme_minimal() +
  ggtitle("Mean Smooth SMR in England", subtitle = "Years 2001 - 2010") +
  scale_x_continuous(breaks = c(2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010))










