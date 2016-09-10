kap<-read.csv("KAP.out")
bus<-read.csv("BUS1.out")
nosex<-read.csv("SPS_out_noasex.txt")
small<-read.csv("SPS_out_small.txt")
pa<-read.csv("PA.out")

X=2.9
Y=6.3

col3=rgb(1,0,0)
col4=rgb(1,0,0)

pdf("f.pdf", width=X, height=Y)
q=list(nosex$fC[nosex$line.A==0 & nosex$min.1 > 2.5 ], small$fC[small$line.A==0 & small$min.1 > 2.5 ], kap$fC[kap$line.A==0 & kap$min.1 > 2.5], pa$fC[pa$line.A==0 & pa$min.1 > 2.5 ],  bus$fC[bus$line.A==0 & bus$min.1 > 2.5])
boxplot(q, ylab=expression(italic(f) ), xlab="population", border="white", names=c("", "", "", "", ""), outline=FALSE, ylim=c(-0.2, 0.4) )
stripchart(q, vertical=TRUE, pch='.', col=col3, method="jitter", add = TRUE )
boxplot(q, names=c("trim", "SPS", "KAP", "PA", "BUS"), add=TRUE, outline=FALSE, las=3 )
dev.off()


