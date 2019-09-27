#' rpmodel.grid
#'
#' rpmodel looping through pixels in parallel and writing on the go 
#' @param  same as rpmodel, tc, vpd, fapar, ppfd, elev, soilm, meanalpha as grid type objects
#' @import raster  
#' @keywords rpmodel
#' @export
#' @examples
#' rpmodel.grid()	
	
rpmodel.grid<-function(tc, vpd, co2, fapar, ppfd, elev = NA, 
		kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.0870, 0.0817), 0.0492), 
		beta = 146.0, soilm = 1.0, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.685, 
		c4 = FALSE, method_optci = "prentice14", method_jmaxlim = "wang17", 
		do_ftemp_kphio = TRUE, do_soilmstress = TRUE,inmem=FALSE,outdir=getwd()){
			
		###############################################################################################
		# 00. create array for results
		###############################################################################################
		rasterOptions(maxmemory=1e9, tmptime = 24, chunksize = 1e8,todisk = FALSE, overwrite=TRUE, tolerance = 0.5)
		y<-as.numeric(unique(format(getZ(tc),'%Y')))
		ny <- nlayers(tc)
		# gpp gC/m2
		gpp<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(gpp)<-extent(elev)
		gpp<-setZ(gpp,getZ(tc))
		# lue g C / mol photons
		lue<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(lue)<-extent(elev)
		gpp<-setZ(gpp,getZ(tc))
		# transpiration mm
		Tr<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(Tr)<-extent(elev)
		Tr<-setZ(Tr,getZ(tc))
		# wue g C / mol H2O
		wue<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(wue)<-extent(elev)
		wue<-setZ(wue,getZ(tc))
		gc()
		setwd(outdir)
		if(inmem){
			tc<-readAll(tc)
			vap<-readAll(vap)
			elev<-readAll(elev)
			fapar<-readAll(fa)
			ppfd<-readAll(ppfd)
			soilm<-readAll(soilm)
			meanalpha<-readAll(meanalpha)
			
		}
		# ************************************************************************
		# Name:     density_h2o
		# Inputs:   - double (tc), air temperature, degrees C
		#           - double (pa), atm pressure, Pa
		# Returns:  double, kg/m^3
		# Features: This function calculates the temperature and pressure
		#           dependent density of pure water
		# * Ref:    Chen, C.T., R.A. Fine, and F.J. Millero (1977), The equation
		#             of state of pure water determined from sound speeds, The
		#             Journal of Chemical Physics 66, 2142;
		#             doi:10.1063/1.434179
		# ************************************************************************
		density_h2o <- function(tc, pa) {
			# Calculate density of water at 1 atm, g/cm^3
			po <- 0.99983952 +
			(6.788260e-5)*tc +
			-(9.08659e-6)*tc*tc +
			(1.022130e-7)*tc*tc*tc +
			-(1.35439e-9)*tc*tc*tc*tc +
			(1.471150e-11)*tc*tc*tc*tc*tc +
			-(1.11663e-13)*tc*tc*tc*tc*tc*tc +
			(5.044070e-16)*tc*tc*tc*tc*tc*tc*tc +
			-(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc
			
			# Calculate the bulk modulus of water at 1 atm, atm
			ko <- 19652.17 +
			148.1830*tc +
			-2.29995*tc*tc +
			0.01281*tc*tc*tc +
			-(4.91564e-5)*tc*tc*tc*tc +
			(1.035530e-7)*tc*tc*tc*tc*tc
			
			# Calculate temperature-dependend coefficients
			ca <- 3.26138 +
			(5.223e-4)*tc +
			(1.324e-4)*tc*tc +
			-(7.655e-7)*tc*tc*tc +
			(8.584e-10)*tc*tc*tc*tc
			
			cb <- (7.2061e-5) +
			-(5.8948e-6)*tc +
			(8.69900e-8)*tc*tc +
			-(1.0100e-9)*tc*tc*tc +
			(4.3220e-12)*tc*tc*tc*tc
			
			# Convert pressure to bar (1 bar = 100000 Pa)
			pbar <- (1e-5)*pa
			
			pw <- (1e3)*po*(ko + ca*pbar + cb*pbar^2)/(ko + ca*pbar + cb*pbar^2 - pbar)
			return(pw)
		}
		###############################################################################################
		# 01. set the clusters for parallel computing
		###############################################################################################	
		cl <- getCluster()
		on.exit( returnCluster() )
		nodes <- length(cl)
		bs<- blockSize(tc, minblocks=nodes*5)
		parallel::clusterEvalQ(cl, library("rpmodel"))
		pmodel<-Vectorize(rpmodel,c("tc","vpd","co2","fapar","ppfd","meanalpha","soilm"))
		parallel::clusterExport(cl, varlist=c("tc","vpd","co2",'elev','fapar',"ppfd","kphio",'beta','soilm','meanalpha','apar_soilm','bpar_soilm','c4','method_optci','method_jmaxlim','do_ftemp_kphio','do_soilmstress','pmodel','ny','density_h2o','y','bs'),envir=environment()) 
		pb <- pbCreate(bs$n)
		pb <- txtProgressBar(min=1,max = bs$n, style = 3)
		###############################################################################################
		# 02. create the functions to send to the workers, split the data in chunks
		###############################################################################################	
		clFun <- function(i) {
			tcrow<-split(getValues(tc,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			vpdrow<-split(getValues(vpd,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			elevrow<-split(getValues(elev,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			faparrow<-split(getValues(fapar,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			ppfdrow<-split(getValues(ppfd,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			soilmrow<-split(getValues(soilm,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			alpharow<-split(getValues(meanalpha,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
			co2row<-lapply(1:(ncol(elev)*bs$nrows[i]),function(j) co2)
			# define function to get transpiration
			calc_transp<-function(mat,vpd,tc, elev){
				gs_mol<-as.numeric(mat[12,])*calc_patm(elev)
				gs_mol[gs_mol<0]<-NA
				T_mol<-1.6*gs_mol*(vpd)/calc_patm(elev)
				T_mm<-T_mol*(18/density_h2o(tc,calc_patm(elev)))
				wue<-as.numeric(mat[10,])/T_mol
				return(rbind(mat,T_mm,wue))
			}
			
			# apply vectorized rpmodel						
			result<-mapply(pmodel, 
				tc             = tcrow,           # temperature, deg C
				vpd            = vpdrow,         # Pa,
				co2            = co2row,          # ppm,
				elv            = elevrow,        # m.a.s.l.,
				kphio          = kphio,         # quantum yield efficiency,
				beta           = beta,         # unit cost ratio a/b,
				fapar          = faparrow,      # fraction  ,
				ppfd           = ppfdrow,      # mol/m2/d month?,
				method_optci   = method_optci,
				method_jmaxlim = method_jmaxlim,
				soilm = soilmrow,
				meanalpha=alpharow, SIMPLIFY = FALSE)
			# calc transpiration and wue	
			result<-mapply(calc_transp,result,vpdrow,tcrow,elevrow, SIMPLIFY = FALSE)
			gc()
			return(result)
		}
		###############################################################################################
		# 03. send tasks to the nodes
		###############################################################################################
		for (i in 1:nodes) {
			parallel:::sendCall(cl[[i]], clFun, i, tag=i)
		}
		###############################################################################################
		# 04. write to the disk on the go, or save to the ram
		###############################################################################################
		if(!inmem){
			gpp<-writeStart(gpp,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'gpp',".","nc"),
						format="CDF",overwrite=TRUE,varname='gpp', varunit='gC/m2',longname='gross primary production',
						xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
			Tr<-writeStart(Tr,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'Tr',".","nc"),
						format="CDF",overwrite=TRUE,varname='Tr', varunit='mm',longname='Transpiration',
						xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
			wue<-writeStart(wue,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'wue',".","nc"),
						format="CDF",overwrite=TRUE,varname='wue', varunit='gC/mol',longname='water use efficiency',
						xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
			lue<-writeStart(lue,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'lue',".","nc"),
						format="CDF",overwrite=TRUE,varname='wue', varunit='gC/mol',longname='light use efficiency',
						xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
			
			
		}else {
			matgpp <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
			matTr<- matrix(ncol=nlayers(tc), nrow=ncell(tc))
			matwue <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
			matlue <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
			endind<-cumsum(bs$nrows*tc@ncols)
			startind<-c(1,endind+1)    
		}
		
		# mbl1<-do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[10,])}))
		###############################################################################################
		# 05. receive results from the nodes
		###############################################################################################	
		for (i in 1:bs$n) {
			
			d <- parallel:::recvOneData(cl)
			# error?
			if (! d$value$success) {
				stop('error!! check the data...')
			}
			# which block is this?
			b <- d$value$tag
			# cat('received block: ',b,'\n'); flush.console();
			if (!inmem) {
				gpp <- writeValues(gpp,do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[10,])})), bs$row[b])
				Tr <- writeValues(Tr, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[16,])})), bs$row[b])
				wue <- writeValues(wue, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[17,])})), bs$row[b])
				lue <- writeValues(lue, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[9,])})), bs$row[b])
				
				
			} else {
				
				matgpp[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[10,])}))
				matTr[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[16,])}))
				matwue[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[17,])}))
				matlue[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[9,])}))
				
			}
			
			# need to send more data?
			ni <- nodes + i
			if (ni <= bs$n) {
				parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
			}
			gc()
			setTxtProgressBar(pb,i)
		}
		###############################################################################################
		# 06. close connection with the files, or assign valueas to the raster objects
		###############################################################################################
		
		if (!inmem) {
			gpp <- writeStop(gpp)
			Tr <- writeStop(Tr)
			wue <- writeStop(wue)
			lue <- writeStop(lue)
			
		} else {
			# gpp
			gpp<-setValues(gpp,matgpp)
			# transpiration
			Tr<-setValues(Tr,matTr)
			# water use efficiency
			wue<-setValues(wue,matwue)
			# light use efficiency
			lue<-setValues(lue,matlue)
						
		}
		close(pb)
		gc()
		return(list(gpp=gpp,lue=lue,Tr=Tr,wue=wue))
		}
