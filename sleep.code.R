library(minfi)
library(IlluminaHumanMethylation450kmanifest)
path0 <- 'C:/Work/Sleep'
targets <- read.csv('C:/Work/Sleep/carskadon_infiniumhdmeth_sample_sheet.csv', 
                    stringsAsFactors = FALSE, skip = 7, sep='\t')
targets$Sample_Group <- substr(targets$Sample_Name, 8,9)
targets$sleep <- 'x'
targets$sleep[substr(targets$Sample_Name, 1,7) %in% c('1131194','1121078','1130195','1121169',
                                                      '1120462','1120581','1121065','1120591',
                                                      '1130376', '1120204')] <- 'y'
targets$Sample_Group1 <- targets$Sample_Group
targets$Sample_Group1[targets$Sample_Name == '1130962t1'] <- 'bad'
targets$Sample_Group1[targets$Sample_Name == '1120591t1'] <- 'bad'
targets$Sample_Group1[targets$Sample_Name == '1130195t1'] <- 'bad'
targets$sleep1 <- targets$sleep
targets$sleep1[targets$Sample_Name == '1130962t1'] <- 'bad'
targets$sleep1[targets$Sample_Name == '1120591t1'] <- 'bad'
targets$sleep1[targets$Sample_Name == '1130195t1'] <- 'bad'
good_sleep   <- targets$sleep1[targets$Sample_Group1!='bad']
good_smp_grp <- targets$Sample_Group1[targets$Sample_Group1!='bad']

all_idata <- file.path(path0, targets$Sentrix_ID, 
                       paste(targets$Sentrix_ID, targets$Sentrix_Position, 'Red.idat', sep='_'))
rgset <- read.450k(all_idata)#, verbose=TRUE)
#rgset1 <- read.450k.exp('/Users/q6600sl/Documents/Work/Sleep/', recursive = TRUE)
#not a good approach, missmatach may occur
densityPlot(rgset, sampGroups=targets$Sample_Group1, main='Density Plot')
densityBeanPlot(rgset, sampGroups=targets$Sample_Group1, sampNames=targets$Sample_Name)
controlStripPlot(rgset, sampNames=targets$Sample_Name, controls = "BISULFITE CONVERSION II")



#Normalization methods
mset.raw <- preprocessRaw(rgset)
mset.ilm <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls", reference = 2)
#reference: the index of the array used for control normalization (? array or control probe?)
mset.swr <- preprocessSWAN(rgset, mset.raw)
mset.swi <- preprocessSWAN(rgset, mset.ilm)



#Diganositc Plot for different normalization methods
par(mfrow=c(2,2))
plotBetasByType(mset.raw[,1], main = "Raw", cex.legend = 0.5)
plotBetasByType(mset.swr[,1], main = "SWAN on Raw", cex.legend = 0.5)
plotBetasByType(mset.ilm[,1], main = "Illumina Normalization", cex.legend = 0.5)
plotBetasByType(mset.swi[,1], main = "SWAN on Illumina", cex.legend = 0.5)



#Multi-demensional Plots of SWAN normalization
par(mfrow=c(1,2))
mdsPlot(mset.swi, numPositions = 1000, sampGroups = targets$Sample_Group1, main='By Time')
mdsPlot(mset.swi[,targets$Sample_Group1!='bad'], numPositions = 1000, 
        sampGroups = good_smp_grp, main='By Time, exclude bad samples')
mdsPlot(mset.swi, numPositions = 1000, sampGroups = targets$sleep1, main='By Sleep')
mdsPlot(mset.swi[,targets$Sample_Group1!='bad'], numPositions = 1000, 
        sampGroups = good_sleep, main='By Sleep, exclude bad samples')



#Simple test of Time effect on M score
msubset <- mset.swi[1:485512,targets$Sample_Group1!='bad'] #485512 means all the probes
M <- getM(msubset, type = "beta", betaThreshold = 0.001)
#by time
dmp <- dmpFinder(M, pheno=good_smp_grp, type="categorical")
head(dmp)
cpgs <- rownames(dmp)[1:10] 
par(mfrow=c(2,5)) 
plotCpg(msubset, cpg=cpgs, pheno=good_smp_grp)
#by sleep
dmps <- dmpFinder(M, pheno=good_sleep, type="categorical")
head(dmps)
cpgss <- rownames(dmps)[1:10] 
par(mfrow=c(2,5)) 
plotCpg(msubset, cpg=cpgss, pheno=good_sleep)



#Simple test of Time effect on Beta score
B <- getBeta(msubset, type = 'Illumina', offset=100)
#by time
dmpb <- dmpFinder(B, pheno=good_smp_grp, type="categorical")
head(dmpb)
cpgsb <- rownames(dmpb)[1:10] 
par(mfrow=c(2,5)) 
plotCpg(msubset, cpg=cpgsb, pheno=good_smp_grp)
#by sleep
dmpbs <- dmpFinder(B, pheno=good_sleep, type="categorical")
head(dmpbs)
cpgsbs <- rownames(dmpbs)[1:10] 
par(mfrow=c(2,5)) 
plotCpg(msubset, cpg=cpgsbs, pheno=good_sleep)