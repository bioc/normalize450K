estimateLC <- function(eSet){
    
    J = ncol(eSet)

    ### model trained on Reinius dataset
    markers = match(rownames(coefs_reinius),hm450$probe_id)
    EST = sapply(1:J,function(j){
        tmp = exprs(eSet[markers,j])
        i = !is.na(tmp)
        solve.QP(
             t(coefs_reinius[i,]) %*% coefs_reinius[i,]
            ,t(coefs_reinius[i,]) %*% tmp[i]
            ,cbind(c(1,1,1,1,1,1),diag(6))
            ,c(1,0,0,0,0,0,0)
            ,meq=1
        )$sol
        })
    EST = t(EST)
    colnames(EST) = colnames(coefs_reinius)
    EST = data.frame(EST)


    ### model trained on whole blood samples (LOLIPOP study)
    markers = match(rownames(coefs_lolipop),hm450$probe_id)
    EST2 = sapply(1:J,function(j){
        tmp = exprs(eSet[markers,j])
        i = !is.na(tmp)
        solve.QP(
             t(coefs_lolipop[i,]) %*% coefs_lolipop[i,]
            ,t(coefs_lolipop[i,]) %*% tmp[i]
            ,cbind(c(1,1,1,1,1),diag(5))
            ,c(1,0,0,0,0,0)
            ,meq=1
        )$sol
        })
    EST2 = t(EST2)
    colnames(EST2) = paste0(colnames(coefs_lolipop),'2')
    EST2 = data.frame(EST2)


    pData(eSet) = cbind(pData(eSet),EST,EST2)
    
    return(eSet)
}