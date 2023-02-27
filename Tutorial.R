###############################################################################################################
# Package preparation:
pkgs = c("Rcpp","parallel","doParallel","foreach","pipeR","tidyr","dplyr","stats")
all(sapply(pkgs , require, character.only = T))

# Check for missing packages and install them:
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Sources:
source("All_sources.R")

###############################################################################################################

# Step 1. Input data:
# pA and pB are your pilot and main datasets, respectively.
data = readRDS("Ensembl_df_10858")

# Step 2. Inference
# Step 2-0. Orr estimator.
m = nrow(data); pi0_orr = orr.estimate.m0(data$pA, data$pB, B = 20)/m

# Step 2-1. alpha selection.
set.seed(1)
rep_alphas = rep.alpha_0501(p1 = data$pA, p2 = data$pB, alphas = 0.05, rep_per = 0.05, M.mc = 100)

# Step 2-2. DEG selection. # output = readRDS(paste0("extra/M_10000_alpha",0.05)); output$adj_alpha
                         # When M.mc = 10000, rep_alphas = 0.04255349.

DEG.proposed = Inference_R(p1 = data$pA, p2 = data$pB, epi0 = pi0_orr, alpha = rep_alphas)$DEGs


###############################################################################################################

# c.f. Data Preparation:
# # Subset genes wrt Ensemble database.
# p_df = read.csv("extra/RST_qNurADG.csv")
# gene_df = read.csv("extra/2021-09-08 GeneINFO_LetsGo.csv")
# gene_df = gene_df[!is.na(gene_df$GeneType),]
# 
# data = gene_df %>% left_join(., p_df, by = "GeneID") %>%
#   dplyr::select(GeneID, pA, pB, pAB) %>%
#   na.omit(.) # Find gene-length for GBP5.
# 
# p_pilot = data$pA; p_main = data$pB
# 
# # Input:
# data = data.frame(pA = p_pilot, pB = p_main)
# # saveRDS(data, "Ensembl_df_10858")


