###############################################################################################################
# Package preparation:
pkgs = c("Rcpp","parallel","doParallel","foreach","pipeR","tidyr","dplyr","stats")
all(sapply(pkgs , require, character.only = T))

# Check for missing packages and install them:
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Sources:
source("All_sources.R")
nominal.alpha = 0.01

###############################################################################################################

# Step 1. Input data:
# pA and pB are your pilot and main datasets, respectively.
data = readRDS("Ensembl_df_10858")

# Step 2. Inference
# Step 2-0. Orr estimator.
m = nrow(data); pi0_orr = orr.estimate.m0(data$pA, data$pB, B = 20)/m

# Step 2-1. alpha selection.
set.seed(1)
rep_alphas = rep.alpha_0501(p1 = data$pA, p2 = data$pB, alphas = nominal.alpha, rep_per = 0.05, M.mc = 10000)

# Step 2-2. DEG selection. # output = readRDS(paste0("extra/M_10000_alpha",0.05)); output$adj_alpha
# When M.mc = 10000, rep_alphas = 0.04255349.

DEG.proposed = Inference_R(p1 = data$pA, p2 = data$pB, epi0 = pi0_orr, alpha = rep_alphas)


###############################################################################################################

