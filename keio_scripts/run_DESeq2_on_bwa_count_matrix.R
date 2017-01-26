args = commandArgs(trailingOnly=TRUE)
library(DESeq2)

data<-read.delim(args[1], header=T, row.names=1)
samples=ncol(data)/3-1
for( i in 1:samples){
  #if( i == samples) next
  for( j in i:samples+1){
      if( i == j) next


      s1 <- (i-1)*3+1
      e1 <- (i-1)*3+3
      n1  <- strsplit(colnames(data)[s1], ".", fixed=T)[[1]][1]
      s2 <- (j-1)*3+1
      e2 <- (j-1)*3+3
      n2<- strsplit(colnames(data)[s2], ".", fixed=T)[[1]][1]

      out_f1 <- paste(paste("DESeq2", n1, n2, sep="-"), "out", sep=".")

      y<-data[,c(s1:e1, s2:e2)]

      print(paste( s1, e1,s2,e2,sep= "_"))
      print("Starting DESeq")
      #前処理(DESeqDataSetオブジェクトの作成)
      data.cl <- c(rep(1, 3), rep(2, 3))#G1群を1、G2群を2としたベクトルdata.clを作成
      colData <- data.frame(condition=as.factor(data.cl))#condition列にクラスラベル情報を格納したcolDataオブジェクトを作成
      d <- DESeqDataSetFromMatrix(countData=y, colData=colData, design=~condition)#DESeqDataSetオブジェクトdの作成

      #本番(DEG検出)
      #d <- estimateSizeFactors(d)
      #d <- estimateDispersions(d)
      #d <- nbinomLRT(d, full= ~condition, reduced= ~1)
      d <- DESeq(d)
      tmp <- results(d)
      p.value <- tmp$pvalue
      p.value[is.na(p.value)] <- 1
      q.value <- tmp$padj
      q.value[is.na(q.value)] <- 1
      ranking <- rank(p.value)

      tmp <- cbind( n1, n2, rownames(data), y, p.value, q.value, ranking)
      write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=F)
  }
}
