########################################################################################

# Modified from the function etc/Sug_Inference8.R

# Rcpp::sourceCpp("Inference_C.cpp");
# 
# source("augmentation_method.R");
# source("replacement_method.R");

Inference_R = function(p1, p2, epi0, alpha){
  
  # Step1. Optimal rejection region FDR1 less than alpha.
  p = Data.Format.Change(p1, p2)
  
  # Data needs to be sorted by p2.
  RR_id = Rejection.Region.id( p, RR_ids = Reject.Region.ids(p, epi0), alpha)
  
  # Step2. Generate original DEG ids of the rejection region (RR).
  if( RR_id$DEG == 0 ){
    DEGs = NULL
  }else{
    # Rec-shape or L-shape.
    DEGs = deg.generate(p, RR_id)
  }
  
  return(list( RR_id = RR_id,
               DEGs = DEGs) )
}

########################################################################################

Data.Format.Change = function(p1, p2){
  
  m = length(p1);
  p = data.frame(p1, p2, org_id = 1:m);
  p = p %>% arrange(p1);
  p$id1 = 0:(m-1); # Generate id1
  p = p %>% arrange(p2); # Sorty by p2
  
  return(p)
}

Reject.Region.ids = function(p, epi0){
  
  # Run Algorithm.
  p = p %>% arrange(p2)
  RR_ids = data.frame(Inference_C(p$p1, p$p2, p$id1, p$org_id, epi0))
  # RR_ids2 = data.frame(Inference_C_int(p$p1, p$p2, p$id1, p$org_id, epi0))
  colnames(RR_ids) = c("DEG","FDR","i","j")
  RR_ids = RR_ids[RR_ids$i!=-2,]; # -2 means NA.
  RR_ids = rbind( c(DEG = 0,FDR = 0,i = -10,j = -10), RR_ids); # -10 means NULL.
  
  return(RR_ids)
}

Rejection.Region.id = function(p, RR_ids, alpha){
  
  # RR_id.
  p = p %>% arrange(org_id) # Sort by original id.
  Opt_id = rev(which(RR_ids$FDR<=alpha))[1]
  RR_id = RR_ids[Opt_id,]
  
  return(RR_id)
}

deg.generate = function(p, RR_id){
  
  p = p %>% arrange(org_id)
  
  if(RR_id$j==-1){
    # case1: Rectangular - Shape.
    DEGs = which((p$p1<=p$p1[RR_id$i])&(p$p2<=p$p2[RR_id$i]))
  }else{
    # case2: L - shape.
    DEGs = which(((p$p1<=p$p1[RR_id$i])&(p$p2<=p$p2[RR_id$i]))|(p$p2<=p$p2[RR_id$j]))
  }
  
  return(DEGs)
}

########################################################################################
