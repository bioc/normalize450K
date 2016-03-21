read450K <- function(idat_files){

    J = length(idat_files)
    ex = file.exists( paste0(idat_files,rep( c('_Grn.idat','_Red.idat'),each=J )))
    if(!all(ex)) stop('Some .idat files are missing')

    sample_names = strsplit(x=idat_files,split='/')
    sample_names = sapply(sample_names,tail,n=1L)

    M = U = matrix(NA_real_   ,nrow=nrow(hm450),ncol=J) #methylated (M) and unmethylated (U) signal intensities
    N = V = matrix(NA_integer_,nrow=nrow(hm450),ncol=J) #number of beads underlying methylated (N) and unmethylated (V) signal intensities
    ctrlGrn = ctrlRed = matrix(NA_real_,nrow=nrow(hm450.controls),ncol=J) # signal intensities of control probes

    ### indices of probes by probe type and color channel
    i1g = hm450[hm450$channel=='Grn' ,]
    i1r = hm450[hm450$channel=='Red' ,]
    i2  = hm450[hm450$channel=='Both',]

    ### read the intensities of the green channel
    for(j in 1:J){
        tmp = readIDAT(paste0(idat_files[j],'_Grn.idat'))$Quants
        means  = tmp[,'Mean']
        nbeads = tmp[,'NBeads']
        
        M[i1g$index,j] = means [i1g$Mi]
        N[i1g$index,j] = nbeads[i1g$Mi]
        
        U[i1g$index,j] = means [i1g$Ui]
        V[i1g$index,j] = nbeads[i1g$Ui]

        M[i2$index ,j] = means [ i2$Mi]
        N[i2$index ,j] = nbeads[ i2$Mi]

        ctrlGrn[,j]    = means[hm450.controls$i]
    }

    ### the red channel
    for(j in 1:J){
        tmp = readIDAT(paste0(idat_files[j],'_Red.idat'))$Quants
        means  = tmp[,'Mean']
        nbeads = tmp[,'NBeads']

        M[i1r$index,j] = means [i1r$Mi]
        N[i1r$index,j] = nbeads[i1r$Mi]

        U[i1r$index,j] = means [i1r$Ui]
        V[i1r$index,j] = nbeads[i1r$Ui]

        U[i2$index ,j] = means [ i2$Ui]
        V[i2$index ,j] = nbeads[ i2$Ui]

        ctrlRed[,j]    = means[hm450.controls$i]
    }

    intensities = list(M=M,U=U,N=N,V=V,ctrlGrn=ctrlGrn,ctrlRed=ctrlRed,sample_names=sample_names)
    return(intensities)
}


dont_normalize450K <- function(intensities){

    if(!all(names(intensities)==c('M','U','N','V','ctrlGrn','ctrlRed','sample_names'))) stop('Invalid argument')

    with(intensities,{
        M[N==0] = NA
        U[V==0] = NA

        M[M<1] = 1
        U[U<1] = 1
        
        meth = M/(M+U)
        
        rownames(meth) = hm450$probe_id
        colnames(meth) = sample_names
        meth = ExpressionSet(assayData=meth)
        return(meth)
    })
}

normalize450K <- function(intensities,tissue=''){

    if(!all(names(intensities)==c('M','U','N','V','ctrlGrn','ctrlRed','sample_names'))) stop('Invalid argument')

    with(intensities,{
        J = ncol(M)
        if(J<2) stop('More than one sample required to perform "between-array" normalization')

        i1g = hm450[hm450$channel=='Grn' ,]
        i1r = hm450[hm450$channel=='Red' ,]
        i2  = hm450[hm450$channel=='Both',]

        cat('[Correcting dye bias]\n')
        ### here I use the extension control probes
        A2C = ctrlRed[7,]/ctrlGrn[9,]
        i = i1g$index
        M[i,] = t(t(M[i,]) * A2C)
        U[i,] = t(t(U[i,]) * A2C)

        A2G = ctrlRed[7,]/ctrlGrn[10,]
        i = i2$index
        M[i,] = t(t(M[i,]) * A2G)

        A2T = ctrlRed[7,]/ctrlRed[8,]
        i = i1r[i1r$next_base=='T',]$index
        M[i,] = t(t(M[i,]) * A2T)
        U[i,] = t(t(U[i,]) * A2T)

        M[N==0] = NA
        U[V==0] = NA

        M[M<1] = 1
        U[U<1] = 1

        hk_1g = hk[hk$channel=='Grn' ,]
        hk_1r = hk[hk$channel=='Red' ,]
        hk_2  = hk[hk$channel=='Both',]

        cat('[Correcting intensity-dependent bias]\n')
        pb <- txtProgressBar(min=0,max=J,style=3)
        ### compute the reference values
        if(tissue=='Blood'){
            rmg = hk_1g$mref
            rug = hk_1g$uref
            rmr = hk_1r$mref
            rur = hk_1r$uref
        }else{
            rmg = log(rowMedians(M[hk_1g$index,],na.rm=TRUE))
            rug = log(rowMedians(U[hk_1g$index,],na.rm=TRUE))
            rmr = log(rowMedians(M[hk_1r$index,],na.rm=TRUE))
            rur = log(rowMedians(U[hk_1r$index,],na.rm=TRUE))
        }

        ### correct intensity-dependent bias
        for(j in 1:J){
            # green channel
            x = log(c(M[hk_1g$index,j],U[hk_1g$index,j]))
            y = x-c(rmg,rug)

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.2,degree=1,family='symmetric')

            x = M[i1g$index,j]; x = x * exp(-predict(f,log(x))); M[i1g$index,j] <- ifelse(is.na(x),M[i1g$index,j],x)
            x = U[i1g$index,j]; x = x * exp(-predict(f,log(x))); U[i1g$index,j] <- ifelse(is.na(x),U[i1g$index,j],x)
            x = M[ i2$index,j]; x = x * exp(-predict(f,log(x))); M[ i2$index,j] <- ifelse(is.na(x),M[ i2$index,j],x)

            # red channel
            x = log(c(M[hk_1r$index,j],U[hk_1r$index,j]))
            y = x-c(rmr,rur)

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.2,degree=1,family='symmetric')

            x = M[i1r$index,j]; x = x * exp(-predict(f,log(x))); M[i1r$index,j] <- ifelse(is.na(x),M[i1r$index,j],x)
            x = U[i1r$index,j]; x = x * exp(-predict(f,log(x))); U[i1r$index,j] <- ifelse(is.na(x),U[i1r$index,j],x)
            x = U[ i2$index,j]; x = x * exp(-predict(f,log(x))); U[ i2$index,j] <- ifelse(is.na(x),U[ i2$index,j],x)
            setTxtProgressBar(pb, j)
        }
        close(pb)

        meth = log(M/U)
        rm(M,U,rmg,rmr,rug,rur)

        ### correct methylation-dependent bias
        cat('[Correcting methylation-dependent bias]\n')
        pb <- txtProgressBar(min=0,max=J,style=3)

        if(tissue=='Blood'){
            rb = hk_2$mref
        }else{
            rb = rowMedians(meth[hk_2$index,],na.rm=TRUE)
        }

        for(j in 1:J){
            x = meth[hk_2$index,j]
            y = x - rb

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.2,family='symmetric',degree=1)
            x = meth[i2$index,j]
            x = x - predict(f,x)
            meth[i2$index,j] <- ifelse(is.na(x),meth[i2$index,j],x)
            setTxtProgressBar(pb,j)
        }

        close(pb)
        rm(rb)
   
        meth = exp(meth)
        meth = meth/(meth+1)

        rownames(meth) = hm450$probe_id
        colnames(meth) = sample_names
        meth = ExpressionSet(assayData=meth)

        return(meth)
        })
}
