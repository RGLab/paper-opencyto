library(plyr)

# Constructs a data.frame containing the Gating Template that will be exported
# to an OpenCyto CSV file 

gating_template <- data.frame(
  alias = "singlet", pop = "singlet+", parent = "root", dims = "FSC-A,FSC-H",
  method = "singletGate", args = "prediction_level=0.999",
  stringsAsFactors = FALSE
)

gating_template <- rbind(gating_template,
                         c("viable", "viable-", "singlet", "ViViD", "flowClust.1d", "neg=1,pos=2"),
                         c("lymph", "lymph+", "viable", "FSC-A,SSC-A", "flowClust.2d", "K=2,target=c(5e4, 2.5e4)"),
                         c("cd3", "cd3+", "lymph", "cd3", "mindensity", ""),
                         c("cd4_pos", "cd4+", "cd3", "cd4", "flowClust.1d", ""),
                         c("cd4_neg", "cd4-", "cd3", "cd4", "flowClust.1d", ""),
                         c("cd8gate_pos", "cd8+", "cd3", "cd8", "flowClust.1d", ""),
                         c("cd8gate_neg", "cd8-", "cd4_neg", "cd8", "flowClust.1d", ""),
                         c("cd4", "cd4+cd8-", "cd3", "cd4,cd8", "refGate", "cd4_pos:cd8gate_pos"),
                         c("cd8", "cd4-cd8+", "cd3", "cd4,cd8", "refGate", "cd4_pos:cd8gate_neg"))


# Cytokine Gates
TNFa <- IFNg <- IL2 <- c(0.995, 0.999, 0.9999)

gt_TNFa <- lapply(TNFa, function(quant) {
  rbind(c(paste0("TNFa", quant), "TNFa+", "cd4", "TNFa", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)),
  c(paste0("TNFa", quant), "TNFa+", "cd8", "TNFa", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)))
})

gt_IFNg <- lapply(IFNg, function(quant) {
  rbind(c(paste0("IFNg", quant), "IFNg+", "cd4", "IFNg", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)),
  c(paste0("IFNg", quant), "IFNg+", "cd8", "IFNg", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)))
})

gt_IL2 <- lapply(IL2, function(quant) {
  rbind(c(paste0("IL2", quant), "IL2+", "cd4", "IL2", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)),
  c(paste0("IL2", quant), "IL2+", "cd8", "IL2", "flowClust.1d", paste0("neg=3,pos=0,quantile=", quant)))
})

gt_cytokines <- do.call(rbind, c(gt_TNFa, gt_IFNg, gt_IL2))
gt_cytokines <- as.data.frame(gt_cytokines, row.names = FALSE, stringsAsFactors = FALSE)
colnames(gt_cytokines) <- colnames(gating_template)

# Cytokine Polyfunctional Gates

cytokine_quantiles <- expand.grid(TNFa = TNFa, IFNg = IFNg, IL2 = IL2)

gt_polyfunction <- ddply(cytokine_quantiles, .(TNFa, IFNg, IL2), function(cyto_quantiles) {

  TNFa_quant <- paste0("TNFa", cyto_quantiles$TNFa)
  IFNg_quant <- paste0("IFNg", cyto_quantiles$IFNg)
  IL2_quant <- paste0("IL2", cyto_quantiles$IL2)

  # Polyfunctions
  cd4_double_pf <- c(paste0("cd4:", paste(IFNg_quant, IL2_quant, sep = "/")),
    paste("cd4", paste(IFNg_quant, IL2_quant, sep = "/"), "subsets"),
    "cd4", "", "polyFunctions",
    paste(paste0("cd4/", c(IFNg_quant, IL2_quant)), collapse = ":"))

  cd4_triple_pf <- c(paste0("cd4:", paste(TNFa_quant, IFNg_quant, IL2_quant, sep = "/")),
    paste("cd4", paste(TNFa_quant, IFNg_quant, IL2_quant, sep = "/"), "subsets"),
    "cd4", "", "polyFunctions",
    paste(paste0("cd4/", c(TNFa_quant, IFNg_quant, IL2_quant)), collapse = ":"))

  cd8_double_pf <- c(paste0("cd8:", paste(IFNg_quant, IL2_quant, sep = "/")),
    paste("cd8", paste(IFNg_quant, IL2_quant, sep = "/"), "subsets"),
    "cd8", "", "polyFunctions",
    paste(paste0("cd8/", c(IFNg_quant, IL2_quant)), collapse = ":"))

  cd8_triple_pf <- c(paste0("cd8:", paste(TNFa_quant, IFNg_quant, IL2_quant, sep = "/")),
    paste("cd8", paste(TNFa_quant, IFNg_quant, IL2_quant, sep = "/"), "subsets"),
    "cd8", "", "polyFunctions",
    paste(paste0("cd8/", c(TNFa_quant, IFNg_quant, IL2_quant)), collapse = ":"))

  gt_pf <- rbind(cd4_double_pf, cd4_triple_pf, cd8_double_pf, cd8_triple_pf)

  # Merges the two gating templates: 1) the quantiles. 2) the polyfunctional gates
  gt_pf <- as.data.frame(gt_pf, row.names = FALSE, stringsAsFactors = FALSE)

  colnames(gt_pf) <- colnames(gating_template)
  gt_pf
})

gt_polyfunction <- subset(gt_polyfunction, select = -c(TNFa, IFNg, IL2))

# Combines the upstream gating template, the cytokine gates, and the polyfunctional gates.
gating_template <- rbind(gating_template, gt_cytokines, gt_polyfunction)

# Writes the Gating Template to a CSV file
write.csv(gating_template, file = "HVTN065-GatingTemplate.csv", row.names = FALSE)


