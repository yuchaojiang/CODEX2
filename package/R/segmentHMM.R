segmentHMM = function( Y_qc, Yhat, optK, K, sampname_qc, ref_qc, chr, mode){
  
  #Normalize reads by expected reads
  YhatOptK = Yhat[[which(K == optK)]]
  #Z=log(pmax(0.01,Y_qc)/YhatOptK)
  Z=(Y_qc-YhatOptK)/sqrt(YhatOptK)
  
  exomtarg=matrix(ncol=3,nrow=length(ref_qc))
  exomtarg[,1]=rep(chr,length(ref_qc))
  exomtarg[,2]=start(ref_qc)
  exomtarg[,3]=end(ref_qc)
  
  finalcalls = list( )
  length( finalcalls ) = length( sampname_qc )
  
  for (k in 1:length(sampname_qc)){
    message("Segmenting sample ", k, ": ", sampname_qc[k], ".")
    zval=Z[,k]
    
    #Define HMM--------------------------------------------------------------------
    
    # emission probability based on mixture of 3 normal distributions 
    # corresponding to deletion, normal, and duplication, respectively
    emission1=dnorm(zval,-3,1)
    emission1[emission1==0]=min(emission1[emission1>0])
    emission2=dnorm(zval,0,1)
    emission2[emission2==0]=min(emission2[emission2>0])
    emission3=dnorm(zval,3,1)
    emission3[emission3==0]=min(emission3[emission3>0])
    
    # transition probability based on physical chromosomal distance
    D=70000
    exomcenter=floor(0.5*(exomtarg[,2]+exomtarg[,3]))
    d=pmax(exomcenter[2:length(exomcenter)]-exomcenter[1:(length(exomcenter)-1)],1)
    f=exp(-d/D)
    p=10^-8
    q=1/6
    transition11=f*(1-q)+(1-f)*p
    transition12=f*q+(1-f)*(1-2*p)
    transition13=(1-f)*p
    transition21=rep(p,length(d))
    transition22=rep(1-2*p,length(d))
    transition23=rep(p,length(d))
    transition31=(1-f)*p
    transition32=f*q+(1-f)*(1-2*p)
    transition33=f*(1-q)+(1-f)*p
    n=length(zval)
    
    #Viterbi algorithm-----------------------------------------------------------------
    
    #populate matrices
    v=matrix(NA,nrow=n,ncol=3)
    pointer=matrix(NA,nrow=(n+1),ncol=3)
    v[1,1]=log(emission1[1]*1*p)
    v[1,2]=log(emission2[1]*1*(1-2*p))
    v[1,3]=log(emission3[1]*1*p)
    pointer[1,1]=0
    pointer[1,2]=0
    pointer[1,3]=0
    for (i in 2:n){
      max1=max((v[i-1,1]+log(transition11[i-1])),(v[i-1,2]+log(transition21[i-1])),(v[i-1,3]+log(transition31[i-1])))
      v[i,1]=log(emission1[i])+max1
      if ((v[i-1,1]+log(transition11[i-1]))==max1){
        pointer[i,1]=1
      } else if ((v[i-1,2]+log(transition21[i-1]))==max1){
        pointer[i,1]=2
      } else{
        pointer[i,1]=3
      }
      max2=max((v[i-1,1]+log(transition12[i-1])),(v[i-1,2]+log(transition22[i-1])),(v[i-1,3]+log(transition32[i-1])))
      v[i,2]=log(emission2[i])+max2
      if ((v[i-1,1]+log(transition12[i-1]))==max2){
        pointer[i,2]=1
      } else if ((v[i-1,2]+log(transition22[i-1]))==max2){
        pointer[i,2]=2
      } else {
        pointer[i,2]=3
      }
      max3=max((v[i-1,1]+log(transition13[i-1])),(v[i-1,2]+log(transition23[i-1])),(v[i-1,3]+log(transition33[i-1])))
      v[i,3]=log(emission3[i])+max3
      if ((v[i-1,1]+log(transition13[i-1]))==max3){
        pointer[i,3]=1
      } else if ((v[i-1,2]+log(transition23[i-1]))==max3){
        pointer[i,3]=2
      } else {
        pointer[i,3]=3
      }
    }
    
    # termination
    max4=max(v[n,1],v[n,2],v[n,3])
    if (v[n,1]==max4){
      pointer[n+1,1]=1
    } else if (v[n,2]==max4) {
      pointer[n+1,1]=2
    } else {
      pointer[n+1,1]=3
    }
    
    # traceback
    traceback=pointer[n+1,1]
    result=rep(NA,n)
    for (i in seq(n+1,2,-1)){
      if (traceback==1){
        result[i-1]="deletion"
        traceback=pointer[i-1,1]
      } else if (traceback==2) {
        result[i-1]="neutral"
        traceback=pointer[i-1,2]
      } else {
        result[i-1]='duplication'
        traceback=pointer[i-1,3]
      }
    }
    
    #Create Output-----------------------------------------------------------------
    delIndex = which(result=="deletion")
    dupIndex = which(result=="duplication")
    if( ( length( delIndex ) + length( dupIndex ) ) == 0 ) {
      finalcalls[[ k ]] = NULL
    } else {
      stDel = delIndex[ !delIndex %in% ( delIndex + 1 ) ]
      edDel = delIndex[ !delIndex %in% ( delIndex - 1 ) ]
      stDup = dupIndex[ !dupIndex %in% ( dupIndex + 1 ) ]
      edDup = dupIndex[ !dupIndex %in% ( dupIndex - 1 ) ]
      st_exon = c( stDel, stDup )
      ed_exon = c( edDel, edDup )
      cnv_type = c( rep( "del", length( stDel ) ),
                    rep( "dup", length( stDup ) ) )
      st_bp = ref_qc@start[ st_exon ]
      ed_bp = ref_qc@start[ ed_exon ] + ref_qc@width[ ed_exon ] - 1
      length_kb = ( ed_bp - st_bp +1 ) / 1000
      raw_cov = sapply( seq_along( st_exon ), function ( j ){
        sum( Y_qc[ st_exon[ j ]:ed_exon[ j ], k ] )
      })
      norm_cov = sapply( seq_along( st_exon ), function ( j ){
        sum( YhatOptK[ st_exon[ j ]:ed_exon[ j ], k ] )
      }) 
      copy_no = round( 2 * ( raw_cov / norm_cov ), 3 )
      if ( mode == "integer" ) copy_no = round( copy_no )
      finalcalls[[ k ]] = data.frame( sampname_qc[ k ], chr, cnv_type, st_bp, ed_bp, 
                                      length_kb, st_exon, ed_exon, raw_cov, norm_cov, 
                                      copy_no)
      colnames( finalcalls[[ k ]] ) = c( 'sample_name','chr','cnv','st_bp','ed_bp', 
                                         'length_kb', 'st_exon', 'ed_exon', 'raw_cov', 
                                         'norm_cov', 'copy_no' )
      finalcalls[[ k ]]  = finalcalls[[ k ]][ order(st_exon), ]
    }
    
    k=k+1
  }
  finalcall = do.call("rbind", finalcalls)
  if( mode == 'integer') finalcall = finalcall[finalcall$copy_no!=2,]
  rownames(finalcall) = paste("cnv", 1:nrow(finalcall), sep = "")
  finalcall
}
