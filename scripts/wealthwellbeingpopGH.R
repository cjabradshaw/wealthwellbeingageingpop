## analysis to examine whether population trajectory affects wealth
## i.e., tests hypothesis that shrinking populations have compromised economies
## Corey Bradshaw, Flinders University
## July 2025

## load libraries
library(boot)
library(dismo)
library(gbm)
library(ggarrange)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(usdm)

# custom functions
modifyVecFunc <- function(obj_name, index, new_value) {
  tryCatch({
    if (!object_exists(obj_name)) {
      stop("object does not exist: ", obj_name)
    }
    
    obj <- get(obj_name)
    
    # check object type and handle accordingly
    if (is.vector(obj) && !is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[index] <- new_value
    } else if (is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[[index]] <- new_value
    } else if (is.data.frame(obj)) {
      if (!all(index <= dim(obj))) stop("index out of bounds")
      obj[index[1], index[2]] <- new_value
    } else {
      stop("unsupported object type")
    }
    
    # Save modified object back to its original name
    assign(obj_name, obj, envir = .GlobalEnv)
    
  }, error = function(e) {
    message("error: ", e$message)
    return(FALSE)
  })
}

object_exists <- function(obj_name) {
  exists(obj_name, envir = .GlobalEnv)
}

## import data
# population data
popdat <- read.csv("popXage.csv", header = TRUE, stringsAsFactors = FALSE)
head(popdat)

## pop total
popdat$Ntot <- round(1000 * apply(popdat[,4:dim(popdat)[2]], MARGIN=1, sum, na.rm=T), 0)
head(popdat)
dim(popdat)

## remove duplicates
popdat <- popdat[!duplicated(popdat[,c("cntry.code", "year")]),]
dim(popdat)

## first year by cntry check
cntry.vec <- attr(table(popdat$cntry.code), "names")
first.year <- rep(NA, length(cntry.vec)) # vector to hold first year for each country
for (c in 1:length(cntry.vec)) {
  cntry <- cntry.vec[c]
  first.year[c] <- min(popdat$year[popdat$cntry.code == cntry])
} # end c loop
first.year

## average population trend by country
r.mean <- r.med <- r.sd <- r.mad <- r.up <- r.lo <-
  r.mean2012.2021 <- r.med2012.2021 <- r.sd2012.2021 <- r.mad2012.2021 <- r.up2012.2021 <- r.lo2012.2021 <- rep(NA, length(cntry.vec))
for (c in 1:length(cntry.vec)) {
  cntry <- cntry.vec[c]
  popvec <- popdat[popdat$cntry.code == cntry, dim(popdat)[2]]
  rvec <- log(popvec[2:length(popvec)] / popvec[1:length(popvec)-1])
  popvec2012.2021 <- popdat[popdat$cntry.code == cntry & popdat$year >= 2012 & popdat$year <= 2021, dim(popdat)[2]]
  rvec2012.2021 <- log(popvec2012.2021[2:length(popvec2012.2021)] / popvec2012.2021[1:length(popvec2012.2021)-1])
  r.mean[c] <- mean(rvec, na.rm = TRUE)
  r.med[c] <- median(rvec, na.rm = TRUE)
  r.sd[c] <- sd(rvec, na.rm = TRUE)
  r.mad[c] <- mad(rvec, na.rm = TRUE)
  r.up[c] <- quantile(rvec, probs = 0.975, na.rm = TRUE)
  r.lo[c] <- quantile(rvec, probs = 0.025, na.rm = TRUE)
  r.mean2012.2021[c] <- mean(rvec2012.2021, na.rm = TRUE)
  r.med2012.2021[c] <- median(rvec2012.2021, na.rm = TRUE)
  r.sd2012.2021[c] <- sd(rvec2012.2021, na.rm = TRUE)
  r.mad2012.2021[c] <- mad(rvec2012.2021, na.rm = TRUE)
  r.up2012.2021[c] <- quantile(rvec2012.2021, probs = 0.975, na.rm = TRUE)
  r.lo2012.2021[c] <- quantile(rvec2012.2021, probs = 0.025, na.rm = TRUE)
} # end c loop
r.dat <- data.frame(cntry.code=cntry.vec, rMean=r.mean, rSD=r.sd, rMed=r.med, rMAD=r.mad, rup=r.up, rlo=r.lo,
                    rMean1221=r.mean2012.2021, rSD1221=r.sd2012.2021, rMed1221=r.med2012.2021,
                    rMAD1221=r.mad2012.2021, rup1221=r.up2012.2021, rlo1221=r.lo2012.2021)
head(r.dat)
dim(r.dat)

# remove duplicates
r.dat <- r.dat[!duplicated(r.dat$cntry.code),]
dim(r.dat)

hist(r.dat$rMean, main="", xlab="mean r 1950-2021", ylab="frequency")

# wealth data (domestic comprehensive wealth index)
# https://data360.worldbank.org/en/indicator/WB_CWON_NW_DOW?compBreak1=WB_CWON_PC
wealthdat <- read.csv("DCWI.csv", header = TRUE, stringsAsFactors = FALSE)
head(wealthdat)

## compare to PPP per-capita GDP
gdppcPPP <- read.csv("gdppcPPP.csv", header = TRUE, stringsAsFactors = FALSE)
head(gdppcPPP)
gdppcPPP2020 <- data.frame(cntry.code=gdppcPPP$cntry.cod, gdppcPPP2020=gdppcPPP$a2020)
head(gdppcPPP2020)

## take latest year only
cntry.vec2 <- attr(table(wealthdat$cntry.code), "names")
max.year <- rep(NA, length(cntry.vec2))
for (c in 1:length(cntry.vec2)) {
  cntry <- cntry.vec2[c]
  max.year[c] <- max(wealthdat$year[wealthdat$cntry.code == cntry])
} # end c loop
max.year

## latest year only (2020)
wealthdat2020 <- wealthdat[wealthdat$year == 2020,]
head(wealthdat2020)
dim(wealthdat2020)

## chained only
wealthdat2020chained <- wealthdat2020[wealthdat2020$UNIT_MEASURE == "USD_REAL_CHAINED_2019",]

## per capita only
wealthdat2020chainedpc <- wealthdat2020chained[wealthdat2020chained$COMP_BREAKDOWN_1_LABEL == "Aggregation: per capita",]
dim(wealthdat2020chainedpc)
head(wealthdat2020chainedpc)
wealthDCWI <- wealthdat2020chainedpc[,c("cntry.code", "DCWI")]
head(wealthDCWI)
dim(wealthDCWI)

## remove duplicates
wealthDCWI <- wealthDCWI[!duplicated(wealthDCWI$cntry.code),]
dim(wealthDCWI)

## compare to GDP pc PPP
head(wealthDCWI)
DCWIgdp <- merge(wealthDCWI, gdppcPPP2020, by = "cntry.code", all.x = TRUE)
head(DCWIgdp)

ggplot(DCWIgdp, aes(x = gdppcPPP2020, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "per-capita GDP PPP (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()
fitDCWIgdp <- lm(log10(DCWI) ~ log10(gdppcPPP2020), data = DCWIgdp)
summary(fitDCWIgdp)

## merge wealth and population data
wealthr <- merge(wealthDCWI, r.dat, by = "cntry.code", all.x = TRUE)
head(wealthr)
dim(wealthr)

## remove duplicates
wealthr <- wealthr[!duplicated(wealthr$cntry.code),]
dim(wealthr)

## plot x,y relationship, with x error bars in ggplot2
ggplot(wealthr, aes(x = rMean, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rMean-rSD, xmax = rMean+rSD), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "median population growth rate (1950-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthr, aes(x = rMed, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rMed-rMAD, xmax = rMed+rMAD), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "median population growth rate (1950-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthr, aes(x = rMean1221, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rlo1221, xmax = rup1221), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "mean population growth rate (1912-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

## poptot data
poptot <- popdat[,c("cntry.code", "year", "Ntot")]
poptot2020 <- poptot[poptot$year == 2020,]
head(poptot2020)
dim(poptot2020)

## remove duplicates
poptot2020 <- poptot2020[!duplicated(poptot2020$cntry.code),]
dim(poptot2020)

## merge wealth and r data
wealthrpop <- merge(wealthr, poptotD2020, by = "cntry.code", all.x = TRUE)
head(wealthrpop)

ggplot(wealthrpop, aes(x = Ntot, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "total population size (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthrpop, aes(x = D, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "total population density (/km2; 2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

## ggplot2 histograms
ggplot(wealthrpop, aes(x = rMean)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "mean population growth rate (1950-2021)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = rMean1221)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "mean population growth rate (1912-2021)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(DCWI))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 per-capita domestic comprehensive wealth index (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(Ntot))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 total population size (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(D))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 population density (2020)", y = "frequency") +
  theme_minimal()

## boosted regression tree
head(wealthrpop)
wealthrpop$lDCWI <- log10(wealthrpop$DCWI)
wealthrpop$lNtot <- log10(wealthrpop$Ntot)
wealthrpop$lD <- log10(wealthrpop$D)
head(wealthrpop)

# r from 1950-2021
head(wealthrpop)
wealthbrt5021 <- gbm.step(wealthrpop, gbm.x = attr(wealthrpop, "names")[c(3,22)],
                    gbm.y = attr(wealthrpop, "names")[20], family="gaussian", max.trees=100000,
                    tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthbrt5021)
barplot(summary(wealthbrt5021)$rel.inf, names.arg = summary(wealthbrt5021)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt5021.summ <- summary(wealthbrt5021)

wealthbrt5021.plot <- ggplot(wealthbrt5021.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt5021.plot.flip <- wealthbrt5021.plot + coord_flip()
wealthbrt5021.plot.flip

wealthbrt5021.CV.cor <- 100 * wealthbrt5021$cv.statistics$correlation.mean
wealthbrt5021.CV.cor.se <- 100 * wealthbrt5021$cv.statistics$correlation.se
print(c(wealthbrt5021.CV.cor, wealthbrt5021.CV.cor.se))

gbm.plot(wealthbrt5021, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))

# r from 2012-2021
head(wealthrpop)
wealthbrt1221 <- gbm.step(wealthrpop, gbm.x = attr(wealthrpop, "names")[c(9,22)],
                          gbm.y = attr(wealthrpop, "names")[20], family="gaussian", max.trees=100000,
                          tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                          tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthbrt1221)
barplot(summary(wealthbrt1221)$rel.inf, names.arg = summary(wealthbrt1221)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt1221.summ <- summary(wealthbrt1221)

wealthbrt1221.plot <- ggplot(wealthbrt1221.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt1221.plot.flip <- wealthbrt1221.plot + coord_flip()
wealthbrt1221.plot.flip

wealthbrt1221.CV.cor <- 100 * wealthbrt1221$cv.statistics$correlation.mean
wealthbrt1221.CV.cor.se <- 100 * wealthbrt1221$cv.statistics$correlation.se
print(c(wealthbrt1221.CV.cor, wealthbrt1221.CV.cor.se))

gbm.plot(wealthbrt1221, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))


##############################
## examine by age structure ##
##############################

## create dependency ratio >65/(16-65)
colnames(popdat)
pop1665.2020 <- apply(popdat[popdat$year == 2020,20:69], MARGIN=1, sum, na.rm=T)
pop66plus.2020 <- apply(popdat[popdat$year == 2020,70:(dim(popdat)[2]-1)], MARGIN=1, sum, na.rm=T)
depratio2020 <- pop66plus.2020 / pop1665.2020
depratio2020.dat <- data.frame(cntry.code = popdat$cntry.code[popdat$year == 2020], depratio = depratio2020)
head(depratio2020.dat)
dim(depratio2020.dat)

## remove duplicates
depratio2020.dat <- depratio2020.dat[!duplicated(depratio2020.dat$cntry.code),]
dim(depratio2020.dat)

## merge with wealth and r data
wealthrpopdepratio <- merge(wealthrpop, depratio2020.dat, by = "cntry.code", all.x = TRUE)
head(wealthrpopdepratio)
dim(wealthrpopdepratio)

## remove duplicates
wealthrpopdepratio <- wealthrpopdepratio[!duplicated(wealthrpopdepratio$cntry.code),]
dim(wealthrpopdepratio)

ggplot(wealthrpopdepratio, aes(x = depratio)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "dependency ratio (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = logit(depratio))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "logit dependency ratio (2020)", y = "frequency") +
  theme_minimal()

## create logit dependency ratio
wealthrpopdepratio$ldepratio <- logit(wealthrpopdepratio$depratio)

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = rMean)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_errorbar(aes(ymin = rMean - rSD, ymax = rMean + rSD), width = 0.01) +
  labs(x = "logit dependency ratio (2020)", y = "mean r (1950-2021") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = rMean1221)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_errorbar(aes(ymin = rMean1221 - rSD1221, ymax = rMean1221 + rSD1221), width = 0.01) +
  labs(x = "logit dependency ratio (2020)", y = "mean r (2012-2021") +
  theme_minimal()

## variance inflation
vif(wealthrpopdepratio[,c("lNtot", "rMean1221", "ldepratio")])
vif(wealthrpopdepratio[,c("lNtot", "rMean", "ldepratio")])

vif(wealthrpopdepratio[,c("lD", "rMean1221", "ldepratio")])
vif(wealthrpopdepratio[,c("lD", "rMean", "ldepratio")])

## correlation matrix
cor(wealthrpopdepratio[,c("lNtot", "lD", "rMean", "ldepratio")], use = "pairwise.complete.obs", method = "spearman")
cor(wealthrpopdepratio[,c("lNtot", "lD", "rMean1221", "ldepratio")], use = "pairwise.complete.obs", method = "spearman")


# r from 2012-2021
head(wealthrpopdepratio)
wealthbrt1221dr <- gbm.step(wealthrpopdepratio, gbm.x = attr(wealthrpopdepratio, "names")[c(9,22,24)],
                          gbm.y = attr(wealthrpopdepratio, "names")[20], family="gaussian", max.trees=100000,
                          tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                          tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthbrt1221dr)
barplot(summary(wealthbrt1221dr)$rel.inf, names.arg = summary(wealthbrt1221dr)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt1221dr.summ <- summary(wealthbrt1221dr)

wealthbrt1221dr.plot <- ggplot(wealthbrt1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt1221dr.plot.flip <- wealthbrt1221dr.plot + coord_flip()
wealthbrt1221dr.plot.flip

wealthbrt1221dr.CV.cor <- 100 * wealthbrt1221dr$cv.statistics$correlation.mean
wealthbrt1221dr.CV.cor.se <- 100 * wealthbrt1221dr$cv.statistics$correlation.se
print(c(wealthbrt1221dr.CV.cor, wealthbrt1221dr.CV.cor.se))

gbm.plot(wealthbrt1221dr, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221dr, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221dr, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratio)
wealthbrt5021dr <- gbm.step(wealthrpopdepratio, gbm.x = attr(wealthrpopdepratio, "names")[c(3,22,24)],
                            gbm.y = attr(wealthrpopdepratio, "names")[20], family="gaussian", max.trees=100000,
                            tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                            tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthbrt5021dr)
barplot(summary(wealthbrt5021dr)$rel.inf, names.arg = summary(wealthbrt5021dr)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt5021dr.summ <- summary(wealthbrt5021dr)

wealthbrt5021dr.plot <- ggplot(wealthbrt5021dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt5021dr.plot.flip <- wealthbrt5021dr.plot + coord_flip()
wealthbrt5021dr.plot.flip

wealthbrt5021dr.CV.cor <- 100 * wealthbrt5021dr$cv.statistics$correlation.mean
wealthbrt5021dr.CV.cor.se <- 100 * wealthbrt5021dr$cv.statistics$correlation.se
print(c(wealthbrt5021dr.CV.cor, wealthbrt5021dr.CV.cor.se))

gbm.plot(wealthbrt5021dr, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021dr, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021dr, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


## add regional categories for spatial resampling step
cont.cntry <- read.csv("continent.country2.csv", header = TRUE, stringsAsFactors = FALSE)
head(cont.cntry)
dim(cont.cntry)

# remove duplicates
cont.cntry <- cont.cntry[!duplicated(cont.cntry$cntry.code),]
dim(cont.cntry)

wealthrpopdepratio.reg <- merge(wealthrpopdepratio, cont.cntry[,c("cntry.code", "cont")], by = "cntry.code", all.x = TRUE)
head(wealthrpopdepratio.reg)
dim(wealthrpopdepratio.reg)

## remove duplicates
wealthrpopdepratio.reg <- wealthrpopdepratio.reg[!duplicated(wealthrpopdepratio.reg$cntry.code),]
dim(wealthrpopdepratio.reg)


table(wealthrpopdepratio.reg$cont)
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cntry.code == "MEX"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "NAM"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "CAR"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "OC"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cntry.code == "CHN"),]

## group regions for sample size increase
# group CAR with SA
wealthrpopdepratio.reg$cont2 <- wealthrpopdepratio.reg$cont
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "CAR"] <- "SACAR"
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "OC"] <- "ASIAOC"
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "ASIA"] <- "ASIAOC"
head(wealthrpopdepratio.reg)

## plot by region
ggplot(wealthrpopdepratio.reg, aes(x = depratio, y = DCWI, color = cont2)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_y_log10() +
  labs(x = "dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal() +
  theme(legend.position = "bottom")

theme1 = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"),
  panel.background = element_rect(fill = 'white'),
  panel.grid.major.y = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.grid.minor.y = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.grid.major.x = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.border = element_blank(),
  legend.title = element_text(size = 16, face="bold"),
  legend.text = element_text(size = 13),
  legend.key.size = unit(1.2, "cm"),
  legend.key.width = unit(0.8,"cm"))


# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=logit(depratio), y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size", color = "region") +
  theme1


# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=logit(depratio), y=DCWI, size = D, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population density") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population density (per km2")

# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=lNtot, y=DCWI, size = ldepratio, color=cont2)) +
  geom_point(alpha=0.6) +
  scale_y_log10() +
  scale_size(range = c(.1, 15), name="logit dependency ratio (2020)") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "log population size (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "logit dependency ratio (2020)")

ggplot(wealthrpopdepratio.reg, aes(x=rMean1221, y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.6) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "mean r (2012-2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size")


## resampled BRT loop
## 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100

table(wealthrpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.8 * min(table(wealthrpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthrpopdepratio.reg$cont2))
reg.vec <- unique(wealthrpopdepratio.reg$cont2)

head(wealthrpopdepratio.reg)
## r mean 1950-2021
traincols <- attr(wealthrpopdepratio.reg, "names")[c(3,18,20)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratio.reg.sub <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratio.reg.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lDCWI=dat.smp$lDCWI,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                            gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                            tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                            tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                     (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                     (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                        (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                            (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                          NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med <- median(CV.cor.update, na.rm=T)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo, CV.cor.med, CV.cor.up))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort

# plot
ri.plt <- ggplot(ri.sort) +
  geom_bar(aes(x=reorder(row.names(ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo,CV.cor.med,CV.cor.up), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort
top.ri.sort <- ri.sort[1:topNvars,]
topNvars.names <- rownames(top.ri.sort)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                         pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("log10 per-capita domestic comprehensive wealth index (2020)"))
}
ggarrange(plt1, plt2, plt3, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                        up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "Pred.csv", sep=""), row.names=F)
}


###################
## r mean 1912-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratio.reg)
traincols <- attr(wealthrpopdepratio.reg, "names")[c(9,18,20)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratio.reg.sub <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratio.reg.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lDCWI=dat.smp$lDCWI,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.1221, CV.cor.med.1221, CV.cor.up.1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.1221

# plot
ri.plt.1221 <- ggplot(ri.sort.1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.1221,CV.cor.med.1221,CV.cor.up.1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.1221
top.ri.sort.1221 <- ri.sort.1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("log10 per-capita domestic comprehensive wealth index (2020)"))
}
ggarrange(plt1.1221, plt2.1221, plt3.1221, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "Pred1221.csv", sep=""), row.names=F)
}


####################
## wellbeing data ##
####################
## Blanchflower, D.G., Bryson, A. Wellbeing Rankings. Soc Indic Res 171, 513–565 (2024)
## https://doi.org/10.1007/s11205-023-03262-y

wellbeing <- read.csv("wellbeingrank.csv", header=T, stringsAsFactors = F)
head(wellbeing)
head(wealthrpopdepratio.reg)

## re-rank overall wellbeing index
wellbeing$WB <- rank(wellbeing$overall, ties.method = "min")
head(wellbeing)
WB <- wellbeing[,c("cntry.code", "WB")]
head(WB)
range(wellbeing$WB)
wealthrpopdepratioWB <- merge(wealthrpopdepratio.reg, WB, by="cntry.code", all.x=T)
head(wealthrpopdepratioWB)

# bubble plots
ggplot(wealthrpopdepratioWB, aes(x=logit(depratio), y=WB, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "wellbeing rank",
       size = "population size", color = "region") +
  theme1

## merge with other data
wealthrpopdepratioWBgdp <- merge(wealthrpopdepratioWB, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdp)

ggplot(wealthrpopdepratioWBgdp, aes(x=gdppcPPP2020, y=WB, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "wellbeing rank", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

head(wealthrpopdepratioWBgdp)
fit.GDP_WB <- lm(WB ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdp)
summary(fit.GDP_WB)


# r from 2012-2021
head(wealthrpopdepratioWB)
dim(wealthrpopdepratioWB)

# remove WB NAs
wealthrpopdepratioWB.noNA <- na.omit(wealthrpopdepratioWB)
dim(wealthrpopdepratioWB.noNA)
str(wealthrpopdepratioWB.noNA)
wealthrpopdepratioWB.noNA$WB <- as.numeric(wealthrpopdepratioWB.noNA$WB)
hist(wealthrpopdepratioWB.noNA$WB, breaks=20, main="", xlab="wellbeing rank", ylab="frequency")

wealthrpopdepratioWB1221dr.brt <- gbm.step(wealthrpopdepratioWB.noNA, gbm.x = attr(wealthrpopdepratioWB.noNA, "names")[c(9,21,24)],
                            gbm.y = attr(wealthrpopdepratioWB.noNA, "names")[27], family="gaussian", max.trees=100000,
                            tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                            tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthrpopdepratioWB1221dr.brt)
barplot(summary(wealthrpopdepratioWB1221dr.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWB1221dr.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWB1221dr.summ <- summary(wealthrpopdepratioWB1221dr)

wealthrpopdepratioWB1221dr.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWB1221dr.plot.flip <- wealthrpopdepratioWB1221dr.plot + coord_flip()
wealthrpopdepratioWB1221dr.plot.flip

wealthrpopdepratioWB1221dr.CV.cor.CV.cor <- 100 * wealthrpopdepratioWB1221dr$cv.statistics$correlation.mean
wealthrpopdepratioWB1221dr.CV.cor.se <- 100 * wealthrpopdepratioWB1221dr$cv.statistics$correlation.se
print(c(wealthrpopdepratioWB1221dr.CV.cor.CV.cor, wealthrpopdepratioWB1221dr.CV.cor.se))

gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratioWB.noNA)
wealthrpopdepratioWB5021dr.brt <- gbm.step(wealthrpopdepratioWB.noNA, gbm.x = attr(wealthrpopdepratioWB.noNA, "names")[c(3,21,24)],
                                       gbm.y = attr(wealthrpopdepratioWB.noNA, "names")[27], family="gaussian", max.trees=100000,
                                       tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                                       tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthrpopdepratioWB5021dr.brt)
barplot(summary(wealthrpopdepratioWB5021dr.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWB5021dr.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWB5021dr.summ <- summary(wealthrpopdepratioWB5021dr.brt)

wealthrpopdepratioWB5021dr.plot <- ggplot(wealthrpopdepratioWB5021dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWB5021dr.plot.flip <- wealthrpopdepratioWB5021dr.plot + coord_flip()
wealthrpopdepratioWB5021dr.plot.flip

wealthrpopdepratioWB5021dr.CV.cor.CV.cor <- 100 * wealthrpopdepratioWB5021dr$cv.statistics$correlation.mean
wealthrpopdepratioWB5021dr.CV.cor.se <- 100 * wealthrpopdepratioWB5021dr$cv.statistics$correlation.se
print(c(wealthrpopdepratioWB5021dr.CV.cor.CV.cor, wealthrpopdepratioWB5021dr.CV.cor.se))

gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))



######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWB.noNA)
traincols <- attr(wealthrpopdepratioWB.noNA, "names")[c(9,21,24)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWB.noNA.sub <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWB.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, WB=dat.smp$WB,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.WB1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.WB1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.WB1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.WB1221, CV.cor.med.WB1221, CV.cor.up.WB1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.WB1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.WB1221

# plot
ri.plt.WB1221 <- ggplot(ri.sort.WB1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.WB1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.WB1221,CV.cor.med.WB1221,CV.cor.up.WB1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.WB1221
top.ri.sort.WB1221 <- ri.sort.WB1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.WB1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".WB1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.WB1221, plt2.WB1221, plt3.WB1221, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "WBPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWB.noNA)
traincols <- attr(wealthrpopdepratioWB.noNA, "names")[c(3,21,24)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWB.noNA.sub <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWB.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, WB=dat.smp$WB,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.WB5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.WB5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.WB5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.WB5021, CV.cor.med.WB5021, CV.cor.up.WB5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.WB5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.WB5021

# plot
ri.plt.WB5021 <- ggplot(ri.sort.WB5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.WB5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.WB5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.WB5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.WB5021,CV.cor.med.WB5021,CV.cor.up.WB5021), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.WB5021
top.ri.sort.WB5021 <- ri.sort.WB5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.WB5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".WB5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.WB5021, plt2.WB5021, plt3.WB5021, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "WBPred5021.csv", sep=""), row.names=F)
}



#############################
## Human Development Index ##
#############################
## https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
hdi <- read.csv("HDI.csv", header=T, stringsAsFactors=F)
head(hdi)

## merge with wealthrpopdepratioWB
wealthrpopdepratioWBhdi <- merge(wealthrpopdepratioWB, hdi, by.x="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdi)

# bubble plots
ggplot(wealthrpopdepratioWBhdi, aes(x=logit(depratio), y=HDI23, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "human development index (2023)",
       size = "population size", color = "region") +
  theme1

## relationship to GDP
gdpHDI <- merge(gdppcPPP2020, hdi, by="cntry.code", all.x=T)
head(gdpHDI)

## merge with other data
wealthrpopdepratioWBgdphdi <- merge(wealthrpopdepratioWB, logit(gdpHDI), by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdphdi)

ggplot(wealthrpopdepratioWBgdphdi, aes(x=gdppcPPP2020, y=logit(HDI23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "logit Human Development Index (2023)", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

fit.GDP_HDI <- lm(logit(HDI23) ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdphdi)
summary(fit.GDP_HDI)

# r from 2012-2021
head(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi)

# remove WB NAs
wealthrpopdepratioWBhdi.noNA <- na.omit(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi.noNA)
str(wealthrpopdepratioWBhdi.noNA)
wealthrpopdepratioWBhdi.noNA$lHDI23 <- logit(wealthrpopdepratioWBhdi.noNA$HDI23)
hist(wealthrpopdepratioWBhdi.noNA$lHDI23, breaks=20, main="", xlab="human development index", ylab="frequency")
head(wealthrpopdepratioWBhdi.noNA)
dim(wealthrpopdepratioWBhdi.noNA)

wealthrpopdepratioWBhdi.brt <- gbm.step(wealthrpopdepratioWBhdi.noNA, gbm.x = attr(wealthrpopdepratioWBhdi.noNA, "names")[c(9,21,24)],
                                       gbm.y = attr(wealthrpopdepratioWBhdi.noNA, "names")[29], family="gaussian", max.trees=100000,
                                       tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                                       tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthrpopdepratioWBhdi.brt)
barplot(summary(wealthrpopdepratioWBhdi.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBhdi.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBhdi.brt.summ <- summary(wealthrpopdepratioWBhdi.brt)

wealthrpopdepratioWBhdi.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBhdi.plot.flip <- wealthrpopdepratioWBhdi.plot + coord_flip()
wealthrpopdepratioWBhdi.plot.flip

wealthrpopdepratioWBhdi.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBhdi.CV.cor.se <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBhdi.CV.cor.CV.cor, wealthrpopdepratioWBhdi.CV.cor.se))

gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi)

wealthrpopdepratioWBhdi.brt <- gbm.step(wealthrpopdepratioWBhdi.noNA, gbm.x = attr(wealthrpopdepratioWBhdi.noNA, "names")[c(3,21,24)],
                                        gbm.y = attr(wealthrpopdepratioWBhdi.noNA, "names")[29], family="gaussian", max.trees=100000,
                                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                                        tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthrpopdepratioWBhdi.brt)
barplot(summary(wealthrpopdepratioWBhdi.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBhdi.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBhdi.brt.summ <- summary(wealthrpopdepratioWBhdi.brt)

wealthrpopdepratioWBhdi.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBhdi.plot.flip <- wealthrpopdepratioWBhdi.plot + coord_flip()
wealthrpopdepratioWBhdi.plot.flip

wealthrpopdepratioWBhdi.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBhdi.CV.cor.se <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBhdi.CV.cor.CV.cor, wealthrpopdepratioWBhdi.CV.cor.se))

gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWBhdi.noNA)
traincols <- attr(wealthrpopdepratioWBhdi.noNA, "names")[c(9,21,24)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBhdi.noNA.sub <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBhdi.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDI23=dat.smp$lHDI23,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDI1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDI1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDI1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDI1221, CV.cor.med.HDI1221, CV.cor.up.HDI1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDI1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDI1221

# plot
ri.plt.HDI1221 <- ggplot(ri.sort.HDI1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDI1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDI1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDI1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.HDI1221,CV.cor.med.HDI1221,CV.cor.up.HDI1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDI1221
top.ri.sort.HDI1221 <- ri.sort.HDI1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDI1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDI1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human development index"))
}
ggarrange(plt1.HDI1221, plt2.HDI1221, plt3.HDI1221, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWBhdi.noNA)
traincols <- attr(wealthrpopdepratioWBhdi.noNA, "names")[c(3,21,24)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBhdi.noNA.sub <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBhdi.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDI23=dat.smp$lHDI23,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDI5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDI5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDI5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDI5021, CV.cor.med.HDI5021, CV.cor.up.HDI5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDI5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDI5021

# plot
ri.plt.HDI5021 <- ggplot(ri.sort.HDI5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDI5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDI5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDI5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.HDI5021,CV.cor.med.HDI5021,CV.cor.up.HDI5021), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDI5021
top.ri.sort.HDI5021 <- ri.sort.HDI5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDI5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDI5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human development index"))
}
ggarrange(plt1.HDI5021, plt2.HDI5021, plt3.HDI5021, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPred5021.csv", sep=""), row.names=F)
}


######################################
## Human Development Index          ##
## adjusted for planetary pressures ##
######################################
## https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
hdipp <- read.csv("HDIPP.csv", header=T, stringsAsFactors=F)
head(hdipp)
hist(logit(hdipp$HDIPP23), main="")

## merge with wealthrpopdepratioWB
wealthrpopdepratioWBhdipp <- merge(wealthrpopdepratioWB, hdipp, by.x="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdipp)

# bubble plots
ggplot(wealthrpopdepratioWBhdipp, aes(x=logit(depratio), y=logit(HDIPP23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "logit planetary pressure-adjusted human development index (2023)",
       size = "population size", color = "region") +
  theme1

## relationship to GDP
gdpHDIPP <- merge(gdppcPPP2020, hdipp, by="cntry.code", all.x=T)
head(gdpHDIPP)

## merge with other data
wealthrpopdepratioWBgdphdipp <- merge(wealthrpopdepratioWB, gdpHDIPP, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdphdipp)

ggplot(wealthrpopdepratioWBgdphdipp, aes(x=gdppcPPP2020, y=logit(HDIPP23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "logit planetary pressure-adjusted Human Development Index (2023)", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

head(wealthrpopdepratioWBgdphdipp)
fit.GDP_HDIPP <- lm(logit(HDIPP23) ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdphdipp)
summary(fit.GDP_HDIPP)

# r from 1950-2021
head(wealthrpopdepratioWBgdphdipp)
dim(wealthrpopdepratioWBgdphdipp)

wealthrpopdepratioWBgdphdipp.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean", "ldepratio", "lNtot", "HDIPP23")])
wealthrpopdepratioWBgdphdipp.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp.noNA)

wealthrpopdepratioWBgdphdipp.brt <- gbm.step(wealthrpopdepratioWBgdphdipp.noNA, gbm.x = attr(wealthrpopdepratioWBgdphdipp.noNA, "names")[c(2,3,4)],
                                        gbm.y = attr(wealthrpopdepratioWBgdphdipp.noNA, "names")[6], family="gaussian", max.trees=100000,
                                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.75,
                                        tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(wealthrpopdepratioWBgdphdipp.brt)
barplot(summary(wealthrpopdepratioWBgdphdipp.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBgdphdipp.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBgdphdipp.brt.summ <- summary(wealthrpopdepratioWBgdphdipp.brt)

wealthrpopdepratioWBgdphdipp.plot <- ggplot(wealthrpopdepratioWBgdphdipp.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBgdphdipp.plot.flip <- wealthrpopdepratioWBgdphdipp.plot + coord_flip()
wealthrpopdepratioWBgdphdipp.plot.flip

wealthrpopdepratioWBgdphdipp.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBgdphdipp.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBgdphdipp.CV.cor.se <- 100 * wealthrpopdepratioWBgdphdipp.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBgdphdipp.CV.cor.CV.cor, wealthrpopdepratioWBgdphdipp.CV.cor.se))

gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

wealthrpopdepratioWBgdphdipp1221.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean1221", "rSD1221", "ldepratio", "lNtot", "HDIPP23","cont2")])
wealthrpopdepratioWBgdphdipp1221.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp1221.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp1221.noNA)

traincols <- attr(wealthrpopdepratioWBgdphdipp1221.noNA, "names")[c(2,5,4)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBgdphdipp1221.noNA.sub <- wealthrpopdepratioWBgdphdipp1221.noNA[wealthrpopdepratioWBgdphdipp1221.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBgdphdipp1221.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBgdphdipp1221.noNA[wealthrpopdepratioWBgdphdipp1221.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDIPP23=dat.smp$lHDIPP23,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDIPP1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDIPP1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDIPP1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDIPP1221, CV.cor.med.HDIPP1221, CV.cor.up.HDIPP1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDIPP1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDIPP1221

# plot
ri.plt.HDIPP1221 <- ggplot(ri.sort.HDIPP1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDIPP1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.HDIPP1221, CV.cor.med.HDIPP1221, CV.cor.up.HDIPP1221))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDIPP1221
top.ri.sort.HDIPP1221 <- ri.sort.HDIPP1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDIPP1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDIPP1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.HDIPP1221, plt2.HDIPP1221, plt3.HDIPP1221, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPPPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

wealthrpopdepratioWBgdphdipp5021.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean", "rSD", "ldepratio", "lNtot", "HDIPP23","cont2")])
wealthrpopdepratioWBgdphdipp5021.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp5021.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp5021.noNA)

traincols <- attr(wealthrpopdepratioWBgdphdipp5021.noNA, "names")[c(2,5,4)] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBgdphdipp5021.noNA.sub <- wealthrpopdepratioWBgdphdipp5021.noNA[wealthrpopdepratioWBgdphdipp5021.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBgdphdipp5021.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBgdphdipp5021.noNA[wealthrpopdepratioWBgdphdipp5021.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDIPP23=dat.smp$lHDIPP23,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                      gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                      tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.5,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step(dat.smp.rsmp, gbm.x = attr(dat.smp.rsmp, "names")[c(3,4,5)],
                        gbm.y = attr(dat.smp.rsmp, "names")[2], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0001, bag.fraction=0.5,
                        tree.complexity = 2, silent=T, step.size = 20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDIPP5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDIPP5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDIPP5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDIPP5021, CV.cor.med.HDIPP5021, CV.cor.up.HDIPP5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDIPP5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDIPP5021

# plot
ri.plt.HDIPP5021 <- ggplot(ri.sort.HDIPP5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDIPP5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDIPP5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDIPP5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.HDIPP5021, CV.cor.med.HDIPP5021, CV.cor.up.HDIPP5021))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDIPP5021
top.ri.sort.HDIPP5021 <- ri.sort.HDIPP5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDIPP5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDIPP5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.HDIPP5021, plt2.HDIPP5021, plt3.HDIPP5021, ncol=1, nrow=3)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPPPred5021.csv", sep=""), row.names=F)
}


###############
## Other plots
head(gdppcPPP2020)
wealthrpopdepratioWBhdiGDP <- merge(wealthrpopdepratioWBhdi, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdiGDP)

ggplot(wealthrpopdepratioWBhdiGDP, aes(x=gdppcPPP2020, y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "PPP-adjusted per-capita gross domestic product (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1
