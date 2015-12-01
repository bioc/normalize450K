read450K <- function(idat_files){

    J = length(idat_files)
    ex = file.exists( paste0(idat_files,rep( c('_Grn.idat','_Red.idat'),each=J )))
    if(!all(ex)) stop('Some .idat files are missing')

    sample_names = strsplit(x=idat_files,split='/')
    sample_names = sapply(sample_names,tail,n=1L)

    M = U = matrix(NA_real_,nrow=nrow(hm450),ncol=J) #methylated (M) and unmethylated (U) signal intensities
    ctrlGrn = ctrlRed = matrix(NA_real_,nrow=nrow(hm450.controls),ncol=J) # signal intensities of control probes

    ### indices of probes by probe type and color channel
    i1g = hm450[hm450$channel=='Grn' ,]
    i1r = hm450[hm450$channel=='Red' ,]
    i2  = hm450[hm450$channel=='Both',]

    ### read the intensities of the green channel
    for(j in 1:J){
        means = readIDAT(paste0(idat_files[j],'_Grn.idat'))$Quants[,'Mean']
        M[i1g$index,j] = means[i1g$Mi]
        U[i1g$index,j] = means[i1g$Ui]
        M[i2$index ,j] = means[i2$Mi ]
        ctrlGrn[,j]    = means[hm450.controls$i]
    }

    ### the red channel
    for(j in 1:J){
        means = readIDAT(paste0(idat_files[j],'_Red.idat'))$Quants[,'Mean']
        M[i1r$index,j] = means[i1r$Mi]
        U[i1r$index,j] = means[i1r$Ui]
        U[i2$index ,j] = means[i2$Ui ]
        ctrlRed[,j]    = means[hm450.controls$i]
    }

    intensities = list(M=M,U=U,ctrlGrn=ctrlGrn,ctrlRed=ctrlRed,sample_names=sample_names)
    return(intensities)
}

dont_normalize450K <- function(intensities){

    if(!all(names(intensities)==c('M','U','ctrlGrn','ctrlRed','sample_names'))) stop('Invalid argument')

    with(intensities,{
        M[M<1] = 1
        U[U<1] = 1
        meth = M/(M+U)
        rownames(meth) = hm450$probe_id
        colnames(meth) = sample_names
        meth = ExpressionSet(assayData=meth)
        return(meth)
    })
}

normalize450K <- function(intensities){

    if(!all(names(intensities)==c('M','U','ctrlGrn','ctrlRed','sample_names'))) stop('Invalid argument')

    with(intensities,{
        J = ncol(M)
        if(J<2) stop('More than one sample required to perform "between-array" normalization')

        i1g = hm450[hm450$channel=='Grn' ,]
        i1r = hm450[hm450$channel=='Red' ,]
        i2  = hm450[hm450$channel=='Both',]

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

        M[M<1] = 1
        U[U<1] = 1

        hk = hm450[hm450$probe_id%in%hk,]
        hk_1g = hk[hk$channel=='Grn' ,]$index
        hk_1r = hk[hk$channel=='Red' ,]$index
        hk_2  = hk[hk$channel=='Both',]$index
        rm(hk)

        ### compute the reference values
        rmg = log(rowMedians(M[hk_1g,]))
        rug = log(rowMedians(U[hk_1g,]))
        rmr = log(rowMedians(M[hk_1r,]))
        rur = log(rowMedians(U[hk_1r,]))

        ### correct intensity-dependent bias
        for(j in 1:J){
            # green channel
            x = log(c(M[hk_1g,j],U[hk_1g,j]))
            y = x-c(rmg,rug)

            f = loess(y~x,span=.2,degree=1,family='symmetric')

            x = M[i1g$index,j]; x = x * exp(-predict(f,log(x))); M[i1g$index,j] <- ifelse(is.na(x),M[i1g$index,j],x)
            x = U[i1g$index,j]; x = x * exp(-predict(f,log(x))); U[i1g$index,j] <- ifelse(is.na(x),U[i1g$index,j],x)
            x = M[ i2$index,j]; x = x * exp(-predict(f,log(x))); M[ i2$index,j] <- ifelse(is.na(x),M[ i2$index,j],x)

            # red channel
            x = log(c(M[hk_1r,j],U[hk_1r,j]))
            y = x-c(rmr,rur)

            f = loess(y~x,span=.2,degree=1,family='symmetric')

            x = M[i1r$index,j]; x = x * exp(-predict(f,log(x))); M[i1r$index,j] <- ifelse(is.na(x),M[i1r$index,j],x)
            x = U[i1r$index,j]; x = x * exp(-predict(f,log(x))); U[i1r$index,j] <- ifelse(is.na(x),U[i1r$index,j],x)
            x = U[ i2$index,j]; x = x * exp(-predict(f,log(x))); U[ i2$index,j] <- ifelse(is.na(x),U[ i2$index,j],x)
        }

        meth = log(M/U)
        rm(M,U,rmg,rmr,rug,rur)

        ### correct methylation-dependent bias
        rb = rowMedians(meth[hk_2,],na.rm=TRUE)

        for(j in 1:J){
            x = meth[hk_2,j]
            y = x - rb

            f = loess(y~x,span=.2,family='symmetric',degree=1)
            x = meth[i2$index,j]
            x = x - predict(f,x)
            meth[i2$index,j] <- ifelse(is.na(x),meth[i2$index,j],x)
        }

        rm(rb)

        meth = exp(meth)
        meth = meth/(meth+1)

        rownames(meth) = hm450$probe_id
        colnames(meth) = sample_names
        meth = ExpressionSet(assayData=meth)

        return(meth)
        })
}