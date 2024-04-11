


tab=read.table("/Users/hoehnemi/Downloads/curhack/Control_siRNA_1.trimmomatic.histo")[,1:2]


  xmax=as.numeric(args[2])

  xmax=max(tab$V1)


tab$V2 <- as.numeric(tab$V2)

pdf("/Users/hoehnemi/Downloads/curhack/kmer_histo.pdf")
plot(tab$V1,tab$V2,xlim=c(0,xmax),type = "l",log="y",xlab="Kmer abundance",ylab="Number of distinct kmers",main="Kmer profile")
grid(lty=1)
lines(tab$V1,tab$V2)
dev.off()

