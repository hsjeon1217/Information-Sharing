###############################################################################################################
# Package preparation:
pkgs = c("ggplot2", "Rcpp","parallel","doParallel","foreach","pipeR","tidyr","dplyr","stats")
all(sapply(pkgs , require, character.only = T))

# Sources:
source("Q_value_sources.R")
source("R_data_generate_0826.R")
sourceCpp("Inference_C.cpp")
source("R_inference_0417.R") 
source("replacement_method_0501.R")

###############################################################################################################

# Step 0. Input data:
data = readRDS("Ensembl_df_10858")

# Step 1. Data exploration using a density plot
density.plot = ggplot(data, aes(x = pA, y = pB) )  +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(color = "black", size = 15, vjust = 2),
        legend.text = element_text(size = 13),
        axis.line = element_line(color = "darkblue", size = 1, linetype = "solid")) + 
  labs( x = TeX('$p_1$'), y = TeX('$p_2$')) + border() + coord_fixed(ratio = 1)+ labs(fill = "Density") 


# Step 2. Inference
# Step 2-0. Orr estimator.
m = nrow(data); pi0_orr = orr.estimate.m0(data$pA, data$pB, B = 20)/m

# Step 2-1. alpha selection.
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


