library(minfi)
path0 <- '/Users/q6600sl/Documents/Work/Sleep/Infinium'
targets <- read.csv('/Users/q6600sl/Documents/Work/Sleep/carskadon_infiniumhdmeth_sample_sheet.csv', 
                     stringsAsFactors = FALSE, skip = 7, sep='\t')
targets$Sample_Group <- substr(targets$Sample_Name, 8,9)
all_idata <- file.path(path0, targets$Sentrix_ID, 
                      paste(targets$Sentrix_ID, targets$Sentrix_Position, 'Red.idat', sep='_'))
rgset <- read.450k(all_idata, verbose=TRUE)
#rgset1 <- read.450k.exp('/Users/q6600sl/Documents/Work/Sleep/', recursive = TRUE)
#not a good approach, missmatach may occur
densityPlot(rgset, sampGroups=targets$Sample_Group, main='Density Plot')
densityBeanPlot(rgset, sampGroups=targets$Sample_Group, sampNames=targets$Sample_Name)
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
par(mfrow=c(1,1))
mdsPlot(mset.swi, numPositions = 1000, sampGroups = targets$Sample_Group)

#Simple test of Time effect on M score
msubset <- mset.swi[1:485512,] #485512 means all the probes
M <- getM(msubset, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=targets$Sample_Group, type="categorical")
head(dmp)
cpgs <- rownames(dmp)[1:4] 
par(mfrow=c(2,2)) 
plotCpg(msubset, cpg=cpgs, pheno=targets$Sample_Group)

#Simple test of Time effect on Beta score
B <- getBeta(msubset, type = 'Illumina', offset=100)
dmpb <- dmpFinder(B, pheno=targets$Sample_Group, type="categorical")
head(dmpb)
cpgsb <- rownames(dmpb)[1:4] 
par(mfrow=c(2,2)) 
plotCpg(msubset, cpg=cpgsb, pheno=targets$Sample_Group)

#mset <- mapToGenome(mset)
#plot(as.matrix(getQC(mset)))
#plotQC(getQC(mset))
