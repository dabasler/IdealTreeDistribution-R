dev.off()
closeAllConnections()
rm(list = ls(all = TRUE))
R.Version ()$version.string
gc()


##############################################################
# load the predifined parameters
path<-"E:/canopy/ITD"
filename<-sprintf("%s/ITD_S3.csv",path)
source(sprintf("%s/ITD.r",path))
##############################################################

#An example dataset
n<-100
ds<-data.frame(dbh=rnorm(n,25,7),species=rep(c("Acer platanoides","Ulmus rubra"),n/2))
ds$dbh[ds$dbh<2]=2

a_plot<-(n/2)^2

# RUN ITD
ITD_par<-import_ITD_parameters(filename)
result<-ITD(ds,a_plot,ITD_par)
plot_treelineup(result)
z<-result[[1]]
ds<-result[[2]]


for (i in 1:nrow(result[[2]])){
  plot_tree(result,i)
  Sys.sleep(0.5)
}
