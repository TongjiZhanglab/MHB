# step0. pre-defined parameters and functions
args <- commandArgs(T)
labelN <- args[1]
dataFN <- args[2]
outFN <- args[3]


# step1. read data and limit the maximum
data <- read.table(dataFN, sep = '\t', header = F, skip = 2)
firstRow <- read.table(dataFN, sep = '\t', header = F, nrows = 2)

# remove outliers
quan_k9 <- quantile(data[which(data[,1]>=0),1], 0.9999)
quan_methylDensity <- quantile(data[which(data[,2]>=0),2], 0.9999)
data[which(data[,1] > quan_k9),1] <- quan_k9
data[which(data[,2] > quan_methylDensity),2] <- quan_methylDensity


# step2. balance
# decide lambda according to mean
lambda <- mean(data[which(data[,2]>=0),2])/mean(data[which(data[,1]>=0),1])
# print(t(c(labelN, lambda)))
write.table(t(c(labelN, lambda)), file = outFN, sep = '\t', quote = F, row.names = F, col.names = F, append = T)