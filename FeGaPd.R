#load XRD data
XRD <- read.table("data/FeGaPd_Exp_PP.txt",header=TRUE)
#load composition data
com <- read.table("data/FeGaPd_Exp_CMP.txt",header=TRUE)
x <- seq(24.40,56.70, by=0.02)

#display sample data
par(mfrow=c(1,1))
b <- as.numeric(XRD[1,])
n <-10
space <- 150
plot(x,b,type="n",ylim=c(0,n*space+5*space))
for (i in 1:n){
b <- as.numeric(XRD[i*27,])
points(x,b+(i-1)*space*1.2,type="l",col=i)
text(27,b[1]+(i-1)*space*1.2+80,paste("row",i*27))
}
dev.copy(png,'figures/display.png')
dev.off()

# set parameters
cluster_No = 5
method="euclidean" # choose distance metric

# plot function
plot_results = function (cluster_No, results, titles){
  triangle_X=1.00-com[ ,2]/2-com[ ,1]
  triangle_Y=0.866*com[ ,2]
  plot(triangle_X, triangle_Y, col=rainbow(cluster_No)[results], 
       pch=19, cex=1, xlim=c(-0.05,0.60), ylim=c(-0.05,0.55))

  segments(0, 0, 0.3, 0.6*0.866)
  segments(0, 0, 0.60, 0)
  segments(0.60, 0,0.3, 0.6*0.866)
  text(0, -0.03, colnames(com[1]), cex=1)
  text(0.6, -0.03, colnames(com[3]), cex=1)
  text(0.3, 0.44, colnames(com[2]), cex=1)
  title_name = paste(titles)
  text(0.3, 0.55, title_name, cex=0.9)
}

#calculate distance
XRD_euclidean_dist=dist(XRD, method)
library(stylo)
XRD_table <- t(XRD)%*%as.matrix(XRD)
XRD_cosine_dist = dist.cosine(XRD_table)

# k-medoidos
library(cluster)
results_kmedoidos = pam(XRD_euclidean_dist, cluster_No)$cluster
results_kmedoidos_cosine = pam(XRD_cosine_dist, cluster_No)$cluster
#hc
hc=hclust(XRD_euclidean_dist, method = "ward.D2")
results_hc <- cutree(hc, k = cluster_No)

# spectral clustering
library(kernlab)
results_spec <- specc(XRD, centers = cluster_No)


#plot
par(mfrow=c(2,2))
plot_results(cluster_No, results_kmedoidos, "euclidean with kmedoidos")
plot_results(cluster_No, results_hc,"hc")
plot_results(cluster_No, results_spec,"spectral clustering")
plot_results(cluster_No, results_kmedoidos_cosine,"cosine with kmedoidos")
dev.copy(png,'figures/clustering.png')
dev.off()


#compare results
library(fpc)
cluster.stats(results_kmedoidos, results_hc)



