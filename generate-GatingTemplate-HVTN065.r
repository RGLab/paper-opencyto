# Constructs a data.frame containing the Gating Template that will be exported
# to an OpenCyto CSV file 

gating_template <- data.frame(
  alias = "singlet", pop = "singlet+", parent = "root", dims = "FSC-A,FSC-H",
  method = "singletGate", args = "prediction_level=0.999",
  stringsAsFactors = FALSE
)

gating_template <- rbind(gating_template,
                         c("viable", "viable-", "singlet", "ViViD", "flowClust", "K=3"),
                         c("lymph", "lymph", "viable", "FSC-A,SSC-A", "flowClust", "K=2,target=c(5e4, 2.5e4)"),
                         c("cd3", "cd3+", "lymph", "cd3", "mindensity", ""),
                         c("cd4_pos", "cd4+", "cd3", "cd4", "flowClust", "K=2"),
                         c("cd4_neg", "cd4-", "cd3", "cd4", "flowClust", "K=2"),
                         c("cd8gate_pos", "cd8+", "cd3", "cd8", "flowClust", "K=2"),
                         c("cd8gate_neg", "cd8-", "cd4_neg", "cd8", "flowClust", "K=2"),
                         c("cd4", "cd4+cd8-", "cd3", "cd4,cd8", "refGate", "cd4_pos:cd8gate_pos"),
                         c("cd8", "cd4-cd8+", "cd3", "cd4,cd8", "refGate", "cd4_pos:cd8gate_neg"))

# Cytokine Gates
cytokine_tolerance <- seq_len(3)

gt_TNFa <- lapply(cytokine_tolerance, function(tol) {
  rbind(c(paste0("TNFa_tol", tol), "TNFa+", "cd4", "TNFa", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)),
  c(paste0("TNFa_tol", tol), "TNFa+", "cd8", "TNFa", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)))
})

gt_IFNg <- lapply(cytokine_tolerance, function(tol) {
  rbind(c(paste0("IFNg_tol", tol), "IFNg+", "cd4", "IFNg", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)),
  c(paste0("IFNg_tol", tol), "IFNg+", "cd8", "IFNg", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)))
})

gt_IL2 <- lapply(cytokine_tolerance, function(tol) {
  rbind(c(paste0("IL2_tol", tol), "IL2+", "cd4", "IL2", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)),
  c(paste0("IL2_tol", tol), "IL2+", "cd8", "IL2", "cytokine", paste0("split='PTID:VISITNO',tol=1e-", tol)))
})

gt_cytokines <- do.call(rbind, c(gt_TNFa, gt_IFNg, gt_IL2))
gt_cytokines <- as.data.frame(gt_cytokines, row.names = FALSE, stringsAsFactors = FALSE)
colnames(gt_cytokines) <- colnames(gating_template)

# Cytokine Polyfunctional Gates
gt_polyfunction <- lapply(cytokine_tolerance, function(tol) {

  TNFa_alias <- paste0("TNFa_tol", tol)
  IFNg_alias <- paste0("IFNg_tol", tol)
  IL2_alias <- paste0("IL2_tol", tol)

  # Polyfunctions
  cd4_double_pf <- c(paste0("cd4:", IL2_alias, "|", IFNg_alias),
                     paste0("cd4:", IL2_alias, "|", IFNg_alias),
                     "cd4", "", "boolGate",
                     paste(paste0("cd4/", c(IL2_alias, IFNg_alias)), collapse = "|"))

  cd8_double_pf <- c(paste0("cd8:", IL2_alias, "|", IFNg_alias),
                     paste0("cd8:", IL2_alias, "|", IFNg_alias),
                     "cd8", "", "boolGate",
                     paste(paste0("cd8/", c(IL2_alias, IFNg_alias)), collapse = "|"))

  cd4_triple_pf <- c(paste0("cd4:", paste(TNFa_alias, IFNg_alias, IL2_alias, sep = "/")),
    paste("cd4", paste(TNFa_alias, IFNg_alias, IL2_alias, sep = "/"), "subsets"),
    "cd4", "", "polyFunctions",
    paste(paste0("cd4/", c(TNFa_alias, IFNg_alias, IL2_alias)), collapse = ":"))

  cd8_triple_pf <- c(paste0("cd8:", paste(TNFa_alias, IFNg_alias, IL2_alias, sep = "/")),
    paste("cd8", paste(TNFa_alias, IFNg_alias, IL2_alias, sep = "/"), "subsets"),
    "cd8", "", "polyFunctions",
    paste(paste0("cd8/", c(TNFa_alias, IFNg_alias, IL2_alias)), collapse = ":"))

  gt_pf <- rbind(cd4_double_pf, cd8_double_pf, cd4_triple_pf, cd8_triple_pf)
  gt_pf <- as.data.frame(gt_pf, row.names = FALSE, stringsAsFactors = FALSE)
  colnames(gt_pf) <- colnames(gating_template)

  gt_pf
})

gt_polyfunction <- do.call(rbind, gt_polyfunction)

# Combines the upstream gating template, the cytokine gates, and the polyfunctional gates.
gating_template <- rbind(gating_template, gt_cytokines, gt_polyfunction)

# Writes the Gating Template to a CSV file
write.csv(gating_template, file = "gt-HVTN065.csv", row.names = FALSE)


