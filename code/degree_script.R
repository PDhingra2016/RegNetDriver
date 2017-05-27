#input final network for degree centrality calculation

args <- commandArgs(TRUE)
library(igraph)
	netdata<-read.table(args[1],header=TRUE)
	net<-graph.data.frame(netdata,directed=TRUE)

	indegree<-degree(net,mode="in",loops = TRUE)
        fname_out=paste(args[1],"_indegree",sep="")
	write.table(file=fname_out,indegree,sep=" ",quote=FALSE)

	outdegree<-degree(net,mode="out",loops = TRUE)
        fname_out=paste(args[1],"_outdegree",sep="")
	write.table(file=fname_out,outdegree,sep=" ",quote=FALSE)

	totaldegree<-degree(net,mode="total",loops = TRUE)
        fname_out=paste(args[1],"_totaldegree",sep="")
	write.table(file=fname_out,totaldegree,sep=" ",quote=FALSE)


	TF_hub_num=0.25*sum(outdegree>0)
	TF_hub<-head(sort(outdegree,decreasing = TRUE),TF_hub_num)
 	write.table(file="TF_hubs",TF_hub,sep="\t",quote=FALSE)	

	pdf("plot_outdegree", width=6, height=6)
	d1<-read.table(paste(args[1],"_outdegree",sep=""))
	h<-hist(d1$x, breaks="Scott", plot=FALSE)
	 plot(h$mids, h$density, log="y", type='b',col="blue",main="Outdegree distribution")


	pdf("plot_indegree", width=6, height=6)
	d1<-read.table(paste(args[1],"_indegree",sep=""))
	hist(d1$x,breaks=50,xlab="in degree", ylab="frquency", main="In degree distribution",col="red")
	mx<-mean(d1$x)
	abline(v=mx,col="blue",lwd=3)

