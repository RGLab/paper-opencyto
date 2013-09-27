library(flowWorkspace)

# Path where the HVTN 080 GatingSet will be stored
gs_path <- "/loc/no-backup/ramey/HVTN/080/gating-set"

# Constructs a GatingSetList containing the manual gates
manual_path <- '/shared/silo_researcher/Gottardo_R/mike_working/HVTN/080/gslist'
gs_manual_paths <- list.dirs(manual_path, recursive = FALSE)
gslist_manual <- lapply(gs_manual_paths, function(gs_path) {
  load_gs(gs_path)
})
gslist_manual <- GatingSetList(gslist_manual)

# Creates a GatingSet containing all HVTN 080 data without any gates
fs_080 <- getData(gslist_manual, "root", max = 400)
gs_HVTN080 <- GatingSet(fs_080)

# Archives the HVTN 080 GatingSet
save_gs(gs_HVTN080, path = gs_path)
