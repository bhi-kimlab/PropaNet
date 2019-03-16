# ----------------------------
# common template for python script
# ----------------------------

#.libPaths("~/lib/R")
library(limma)
args <- commandArgs(TRUE)
eset <- read.table(args[1], header=TRUE, row.names=1,check.names=FALSE)
genes <- rownames(eset)

group_info.raw <- read.table(args[2], header=FALSE, row.names=1)
group_info <- as.vector(group_info.raw$V2)
names(group_info) <- rownames(group_info.raw)

deg_info <- read.table(args[3], header=FALSE)

if (args[4] == 'absolute_binary'){
	isSigned <- FALSE
	isBinary <- TRUE
	isZstats <- FALSE
	isLogFC <- FALSE
} else if (args[4] == 'signed_binary'){
	isSigned <- TRUE
	isBinary <- TRUE
	isZstats <- FALSE
	isLogFC <- FALSE
} else if (args[4] == 'absolute_zstats'){
	isSigned <- FALSE
	isBinary <- FALSE
	isZstats <- TRUE
	isLogFC <- FALSE
} else if (args[4] == 'signed_zstats'){
	isSigned <- TRUE
	isBinary <- FALSE
	isZstats <- TRUE
	isLogFC <- FALSE
} else if (args[4] == 'absolute_logFC'){
	isSigned <- FALSE
	isBinary <- FALSE
	isZstats <- FALSE
	isLogFC <- TRUE
} else if (args[4] == 'signed_logFC'){
	isSigned <- TRUE
	isBinary <- FALSE
	isZstats <- FALSE
	isLogFC <- TRUE
}

for (i in 1:nrow(deg_info)){
	group1 <- as.character(deg_info$V1[i])
	group2 <- as.character(deg_info$V2[i])
	if (group1 == group2){
		result <- rep(0, nrow(eset))
		names(result) <- rownames(eset)
		result <- result[order(names(result))]
	} else {
		eset.sub <- eset[1:nrow(eset),group_info[names(eset)] == group1 | group_info[names(eset)] == group2]
		trts <- factor(c(rep("A",sum(group_info[names(eset)]  == group1)),rep("B",sum(group_info[names(eset)] == group2 ))))
		#design <- model.matrix(~factor(group_info[names(eset.sub)]))
		design <- model.matrix(~trts)
		colnames(design) <- c("A","B")
		fit <- lmFit(eset.sub, design)
		fit <- eBayes(fit)
		options(digits=2)
		table <- topTable(fit, coef=2, adjust="fdr", n=Inf)
		row.names(table) <- table$ID
		table <- table[order(row.names(table)),]
		#val.p <- table[1:nrow(table),4]
		val.p <- table$P.Value
		val.p[is.na(val.p)] <- 1.0
		val.adjp <- table$adj.P.Val
		if(is.null(val.adjp)){
			val.adjp <- rep(1.0, nrow(table))
		} else{
			val.adjp[is.na(val.adjp)] <- 1.0
		}
		#val.logFC <- table[1:nrow(table),1]
		val.logFC <- table$logFC
		#val.sign <- sign(table[1:nrow(table),1])
		val.sign <- sign(table$t)
		# ID  logFC AveExpr     t P.Value adj.P.Val    B
		# columns of table 1:logFC, 2:AvgExpr 3:t 4:P.value 5:adj.P.val 6:B
		if (isBinary == TRUE){
			result <- as.integer(val.adjp < 0.05)
		}
		else if (isZstats == TRUE){
			result <- pmin(abs(qnorm(val.p)),300)
		}
		else if (isLogFC == TRUE){
			result <- abs(val.logFC)
		}
		if (isSigned == TRUE){
			#result <- result * sign(table[1:nrow(table),1])
			result <- result * val.sign
		}
		names(result) <- rownames(table)
	}
	if (i==1){
		result.all <- cbind(result)
	} else{
		result.all <- cbind(result.all, result)
	}
	colnames(result.all)[ncol(result.all)] <- paste("",group2, sep="")
}

result.all[is.na(result.all)] <- 0

filename=args[5]
write.table(t(c("gene", colnames(result.all))), file=filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(result.all, file=filename, append=TRUE, sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
