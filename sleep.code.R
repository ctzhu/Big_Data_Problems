library(minfi)
library(IlluminaHumanMethylation450kmanifest)
path0 <- '/Users/q6600sl/Documents/Work/Sleep/Infinium'
targets <- read.csv('/Users/q6600sl/Documents/Work/Sleep/carskadon_infiniumhdmeth_sample_sheet.csv', 
                    stringsAsFactors = FALSE, skip = 7, sep='\t')
targets$Sample_Group <- substr(targets$Sample_Name, 8,9)
targets$sleep <- 'x'
targets$sleep[substr(targets$Sample_Name, 1,7) %in% c('1131194','1121078','1130195','1121169',
                                                      '1120462','1120581','1121065','1120591',
                                                      '1130376','1120204')] <- 'y'
targets$mood  <- 'a'
targets$mood[substr(targets$Sample_Name, 1,7) %in% c('1120658','1130309','1130496','1131164',
                                                     '1120567','1130976','1120524','1130962','1130198',
                                                     '1120581','1121065','1120591','1130376','1120204')] <- 'b'
targets$Sample_Group1 <- targets$Sample_Group
targets$sleep1        <- targets$sleep
targets$mood1         <- targets$mood

for (v in c('1130962t1', '1120591t1', '1130195t1', '1120285t1', '1120581t2')){
  targets$Sample_Group1[targets$Sample_Name == v] <- 'bad'
  targets$sleep1[targets$Sample_Name == v]        <- 'bad'
  targets$mood1[targets$Sample_Name == v]         <- 'bad'
}

good_mood    <- targets$mood1[targets$Sample_Group1!='bad']
good_sleep   <- targets$sleep1[targets$Sample_Group1!='bad']
good_smp_grp <- targets$Sample_Group1[targets$Sample_Group1!='bad']
good_chip    <- targets$Sentrix_ID[targets$Sample_Group1!='bad']

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
mdsPlot(mset.swi, numPositions = 1000, sampGroups = targets$mood1, main='By Sleep')
mdsPlot(mset.swi[,targets$Sample_Group1!='bad'], numPositions = 1000, 
        sampGroups = good_mood, main='By mood, exclude bad samples')



#Simple test of Time effect on M score
MB_str  <- c('Using M score',' Using Beta score')
fct_str <- c('Time effect','Sleep effect','Mood effect')
factors <- list(good_smp_grp, good_sleep, good_mood)
list_MB <- list(list(0,0,0), list(0,0,0))
msubset <- mset.swi[1:485512,targets$Sample_Group1!='bad'] #485512 means all the probes
M <- getM(msubset, type = "beta", betaThreshold = 0.001)
B <- getBeta(msubset, type = 'Illumina', offset=100)
MandBeta<- list(M, B)
for (j in c(1,2)){
  #by time, by sleep, and by mood
  for (i in c(1:3)){
    dmp <- dmpFinder(MandBeta[[j]], pheno=factors[[i]], type="categorical")
    print(paste(MB_str[j], 'Test for', fct_str[i], 'Top 5', sep=' '), quote=FALSE)
    print('-----', quote=FALSE)
    print(head(dmp))
    print('-----', quote=FALSE)
    cpgs <- rownames(dmp)[1:10] 
    par(mfrow=c(2,5)) 
    plotCpg(msubset, cpg=cpgs, pheno=factors[[i]])
    title(main=paste(MB_str[j], 'Test for', fct_str[i], sep=' '), line=-1, outer=T)
    list_MB[[j]][[i]] <- dmp
  }  
}


#mset <- mapToGenome(mset)
#plot(as.matrix(getQC(mset)))
#plotQC(getQC(mset))
