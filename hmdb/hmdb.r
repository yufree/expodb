load('Z:/PM/UTMeta/CHEAR_PF08/UntargetedPipeline/RPN/AllBatches/xcms_CAMERA/2018_CHEAR_PF08_S048_AllBatch_RPN_fill.RData')
library(xcms)
fill@filepaths <- gsub('Y:','Z:',fill@filepaths)
eic <- getEIC(fill,rt = "corrected",groupidx = 1:nrow(fill@groups),step = 1)
saveRDS(eic,'PF08.RDS')

path <- 'Z:/PM/UTMeta/entact/hmdbnew.csv'
library(tidyverse)
hmdbnew <- readr::read_csv(path)

mz <- unique(hmdbnew$monisotopic_molecular_weight)
mz0 <- hmdbnew$monisotopic_molecular_weight
# mz <- mz0[!is.na(mz0)]

accession <- hmdbnew$accession[!is.na(mz0)]

formula <- hmdbnew$chemical_formula[!is.na(mz0)]

class <- hmdbnew$class[!is.na(mz0)]

supclass <- hmdbnew$super_class[!is.na(mz0)]

kegg <- hmdbnew$kegg[!is.na(mz0)]

mz <- mz[!is.na(mz0)]
mz2 <- mz0[!is.na(mz0)]

hmdbclass <- cbind.data.frame(mz=mz2,supclass)
hmdbclass0 <- cbind.data.frame(mz=mz2,supclass,class,formula,accession)

hmdbclass2 <- hmdbclass[complete.cases(hmdbclass),]
hmdbclass3 <- hmdbclass0[complete.cases(hmdbclass),]

splitclass <- split(hmdbclass2,hmdbclass2$supclass,drop = T)



getrda2 <- function(df0){
  m <- as.matrix(df0$mz)
  rownames(m) <- df0$mz
  dis <- dist(m, method = "manhattan")
  
  
  #dfx <- reshape2::melt(as.matrix(lower.tri(dis)), varnames = c("ms1", "ms2"))
  df <- data.frame(
    ms1 = df0$mz[which(lower.tri(dis), arr.ind = T)[,1]],
    ms2 = df0$mz[which(lower.tri(dis), arr.ind = T)[,2]],
    diff = as.numeric(dis),
    diff2 = round(as.numeric(dis),digits = 2)
  )
  
  freq <- table(df$diff2)[order(table(df$diff2), decreasing = T)]
  freq <- freq[freq>1]
  sda <- df[(df$diff2 %in% c(as.numeric(names(
    freq[freq > 1]
  )))),]
  pmd <- unique(sda$diff2)[order(unique(sda$diff2))]
  
  dfx <- NULL
  
  split <- split.data.frame(sda,sda$diff2)
  rtpmd <- function(bin, i) {
    mass <- unique(c(bin$ms1[bin$diff2==i],bin$ms2[bin$diff2==i]))
    index <- df0$mz %in% mass
    return(as.integer(index))
  }
  result <- mapply(rtpmd, split, as.numeric(names(split)))
  rownames(result) <- df0$mz
  return(result)
}
memory.limit(size=106000000)

class1 <- getrda2(splitclass[[1]])
class2 <- getrda2(splitclass[[2]])

class3 <- getrda2(splitclass[[3]])
class4 <- getrda2(splitclass[[4]])
class5 <- getrda2(splitclass[[5]])
class6 <- getrda2(splitclass[[6]])
class7 <- getrda2(splitclass[[7]])
class8 <- getrda2(splitclass[[8]])
class9 <- getrda2(splitclass[[9]])
class10 <- getrda2(splitclass[[10]])

class11 <- getrda2(splitclass[[11]])
class12 <- getrda2(splitclass[[12]])

class13 <- getrda2(splitclass[[13]])
class14 <- getrda2(splitclass[[14]])
class15 <- getrda2(splitclass[[15]])
class16 <- getrda2(splitclass[[16]])
class17 <- getrda2(splitclass[[17]])
class18 <- getrda2(splitclass[[18]])
class19 <- getrda2(splitclass[[19]])
class20 <- getrda2(splitclass[[20]])

class21 <- getrda2(splitclass[[21]])
class22 <- getrda2(splitclass[[22]])

class23 <- getrda2(splitclass[[23]])
class24 <- getrda2(splitclass[[24]])
class25 <- getrda2(splitclass[[25]])
class26 <- getrda2(splitclass[[26]])
class27 <- getrda2(splitclass[[27]])
class28 <- getrda2(splitclass[[28]])
# 1, 8, 10, 11, 14, 20, 26,
# Acetylides, Inorganic compounds, Lipids and lipid-like molecules, Miscellaneous inorganic compounds, Organic 1,3-dipolar compounds, Organic salts, Organophosphorus compounds  
class <- list(data = hmdbclass3,
              AlkaloidsAndDerivatives = class2, 
              Benzenoids = class3, 
              HomogeneousMetalCompounds = class4, 
              HomogeneousNonmetalCompounds = class5,
              HydrocarbonDerivatives = class6,
              Hydrocarbons = class7,
              LignansNeolignansAndRelatedCompounds = class9,
              MixedMetalNonmetalCompounds = class12,
              NucleosidesNucleotidesAndAnalogues = class13,
              OrganicAcidsAndDerivatives = class15,
              OrganicCompounds = class16,
              OrganicNitrogenCompounds = class17,
              OrganicOxygenCompounds = class18,
              OrganicPolymers = class19,
              OrganohalogenCompounds = class21,
              OrganoheterocyclicCompounds = class22,
              OrganometallicCompounds = class23,
              OrganonitrogenCompounds = class24,
              OrganooxygenCompounds = class25,
              OrganosulfurCompounds = class27,
              PhenylpropanoidsAndPolyketides = class28
              )
saveRDS(class,'hmdbclass.RDS',compress = 'xz')

mzr <- mz %% 1
m <- as.matrix(mz)
rownames(m) <- mz
dis <- parallelDist::parDist(m, method = "manhattan")


#dfx <- reshape2::melt(as.matrix(lower.tri(dis)), varnames = c("ms1", "ms2"))
df <- data.frame(
  ms1 = mz[which(lower.tri(dis), arr.ind = T)[,1]],
  ms2 = mz[which(lower.tri(dis), arr.ind = T)[,2]],
  diff = as.numeric(dis),
  diff2 = round(as.numeric(dis),digits = 2)
)

freq <- table(df$diff2)[order(table(df$diff2), decreasing = T)]

pmd <- unique(df$diff2)[order(unique(df$diff2))]
dfx <- NULL

split <- split.data.frame(df,df$diff2)

rtpmd <- function(bin, i) {
  mass <- unique(c(bin$ms1[bin$diff2==i],bin$ms2[bin$diff2==i]))
  index <- mz %in% mass
  return(as.integer(index))
}

result <- mapply(rtpmd, split, as.numeric(names(split)))

rownames(result) <- mz
saveRDS(result,'hmdb.RDS',compress = 'xz')

library(pmd)
t <- pmd::getrda(mz,top = 1000)
t2 <- pmd::getrda(mz,top = 100)

n <- as.numeric(colnames(t))
nr <- n %% 1

all <- apply(result,1,sum)
all2 <- apply(result,2,sum)
allp <- all/length(all2)
all2p <- all2/length(all)
hmdbp <- list(massp=allp,pmdp=all2p)
saveRDS(hmdbp,'hmdbp.RDS',compress = 'xz')
