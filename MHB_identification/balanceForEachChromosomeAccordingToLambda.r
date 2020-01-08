# step0. pre-defined parameters and functions
args <- commandArgs(T)
lambda <- as.numeric(args[1])
dataFN <- args[2]
outFN <- args[3]


# step1. read data and limit the maximum
data <- read.table(dataFN, sep = '\t', header = F, skip = 2)
firstRow <- read.table(dataFN, sep = '\t', header = F, nrows = 2)


# step2. balance
data_balanced_k9 <- data[,1]*lambda

# set bins without signal to -1
data_balanced_k9[data_balanced_k9<0] <- -1

# merge into one data frame
data_balanced <- as.data.frame(cbind(as.character(round(data_balanced_k9,5)),as.character(data[,2])))
colnames(data_balanced) <- c('V1','V2')
data_balanced <- rbind(firstRow, data_balanced)

# write to file
write.table(data_balanced, file = outFN, sep = '\t', quote = F, row.names = F, col.names = F)