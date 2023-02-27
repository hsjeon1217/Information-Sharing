########################################################################################
# Motivated from replacement_method.R
# rep.alpha = function(p1, p2, alphas, rep_per = 0.1, M.mc = 500){
#   
#   m = length(p1)
#   
#   fdr_rel = lapply(1:M.mc, function(.){ 
#     
#     rep.fdr.relation( p1, p2, 
#                       alpha_seq = seq(0, 1, 0.001), rep_m = round(rep_per*m))
#     
#   }) %>% Reduce("+", .)/M.mc
#   
#   # head(fdr_rel,200); plot(fdr_rel, xlim=c(0,1), ylim=c(0,1))
#   
#   rep_alphas = sapply(alphas, function(alpha){ fdr_rel[max(which(fdr_rel$kFDR <= alpha)),]$hFDR })
#   
#   return(rep_alphas)
#   
# }

rep.alpha_0501 = function(p1, p2, alphas, rep_per = 0.1, M.mc = 500){
  
  m = length(p1)
  alpha_seq = seq(0, 1, 0.001)
  
  fdr_rel = lapply(1:M.mc, function(.){ 
    
    rep.fdr.relation( p1, p2, alpha_seq, rep_m = round(rep_per*m))
    
  }) %>% Reduce("+", .)/M.mc
  
  # head(fdr_rel,200); plot(fdr_rel, xlim=c(0,1), ylim=c(0,1))
  
  rep_alphas = sapply(alphas, function(alpha){ fdr_rel[max(which(fdr_rel$kFDR <= alpha)),]$hFDR })
  
  return(rep_alphas)
  
}


rep.fdr.relation = function(p1, p2, alpha_seq, rep_m){
  
  # Step0.
  m = length(p1);
  
  # Step1. Optimal rejection region FDR1 less than alpha.
  np1 = runif(rep_m, min = 0, max = 1);
  np2 = runif(rep_m, min = 0, max = 1);
  
  sampled_ids = sample(m, rep_m)
  tp1 = p1;  tp2 = p2;
  tp1[sampled_ids] = np1; tp2[sampled_ids] = np2;
  Tepi0 = orr.estimate.m0(tp1, tp2, B = 20)/m
  
  # Data needs to be sorted by p2.
  p = Data.Format.Change(p1 = tp1, p2 = tp2);
  RR_ids = Reject.Region.ids(p, epi0 = Tepi0)
  
  sapply(1:nrow(RR_ids), function(i){
    
    RR_id = RR_ids[i,];
    
    if(RR_id$DEG == 0){
      
      kFDR = 0
      
    }else if(RR_id$j==-1){
      
      NP = sum(as.logical(((np1<=tp1[RR_id$i])*(np2<=tp2[RR_id$i]))));
      kFDR = (NP*(Tepi0*m)/rep_m)/RR_id$DEG;
      
    }else{
      
      NP = sum(as.logical(((np1<=tp1[RR_id$i])*(np2<=tp2[RR_id$i])) + (np2<=tp2[RR_id$j])));
      kFDR = (NP*(Tepi0*m)/rep_m)/RR_id$DEG;
      
    }
    return(kFDR)
  }) -> kFDR
  
  TFDR = data.frame(kFDR = kFDR, hFDR = RR_ids[,2])
  
  FDR_mat = lapply(alpha_seq, function(alpha){
    
    TFDR[rev(which(TFDR$hFDR<=alpha))[1],] # VVV
    
  }) %>% do.call("rbind",.)
  
  return( FDR_mat )
}

########################################################################################
