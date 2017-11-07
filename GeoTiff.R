library("raster")
library("rgdal")
new_elev <- as.matrix(read.table("B.txt"))
new_elev_rast <- raster(new_elev) 
new_elev_tiff <- writeRaster(new_elev_rast, filename="new_elevfel.tif", format="GTiff", overwrite=TRUE)
# DInf flow directions
system("mpiexec -n 8 DinfFlowdir -ang new_elevang.tif -slp new_elevslp.tif -fel new_elevfel.tif",show.output.on.console=F,invisible=F)
slp=raster("new_elevslp.tif")
slp[is.na(slp)] <- 0
slp_matrix<- as.matrix(slp,nrow= nrow(slp), ncol=ncol(slp),rownames = FALSE, colnames = FALSE)
flowdir=raster("new_elevang.tif")
flowdir[is.na(flowdir)] <- 8.0   ###The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...
flowdir <- as.matrix(flowdir,nrow= nrow(flowdir), ncol=ncol(flowdir), rownames=FALSE, colnames=FALSE)
write.table(flowdir, file="flowdir_new.txt", sep=",")
write.table(slp_matrix, file="slp_matrix.txt", sep=",")