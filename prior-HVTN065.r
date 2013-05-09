library(ProjectTemplate)
load.project()

# Loads the necessary pipeline code
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

# Customizes the path where the GS is unarchived.
# By default, /tmp is used and often can run out of space, due to the large nc files.
archive_path <- '/loc/no-backup/ramey'


gs_negctrl <- unarchive(file.path(archive_path, "HVTN065-negctrl-42patients-cytokines995.tar"), archive_path)
gs_ENV <- unarchive(file.path(archive_path, "HVTN065-ENV-42patients-cytokines995.tar"), archive_path)
gs_GAG <- unarchive(file.path(archive_path, "HVTN065-GAG-42patients-cytokines995.tar"), archive_path)

prior_list <- list()
prior_list$negctrl <- list()
prior_list$ENV <- list()
prior_list$GAG <- list()

gating(gating_template, gs_negctrl, batch = TRUE)
gating(gating_template, gs_ENV, batch = TRUE)
gating(gating_template, gs_GAG, batch = TRUE)

stimulation <- "GAG"

channel_IL2 <- markers2channels(getData(wf)[[1]], "IL2")
channel_IFNg <- markers2channels(getData(wf)[[1]], "IFNg")
channel_TNFa <- markers2channels(getData(wf)[[1]], "TNFa")

parent <- "CD4"
prior_IL2 <- prior_flowClust1d(flow_set = getData(wf, parent),
                               channel = channel_IL2, K = 2)
prior_IFNg <- prior_flowClust1d(flow_set = getData(wf, parent),
                                channel = channel_IFNg, K = 2)
prior_TNFa <- prior_flowClust1d(flow_set = getData(wf, parent),
                                channel = channel_TNFa, K = 2)

prior_list[[stimulation]]$CD4 <<- list(IL2 = prior_IL2, IFNg = prior_IFNg, TNFa = prior_TNFa)

parent <- "CD8"
prior_IL2 <- prior_flowClust1d(flow_set = getData(wf, parent),
                               channel = channel_IL2, K = 2)
prior_IFNg <- prior_flowClust1d(flow_set = getData(wf, parent),
                                channel = channel_IFNg, K = 2)
prior_TNFa <- prior_flowClust1d(flow_set = getData(wf, parent),
                                channel = channel_TNFa, K = 2)

prior_list[[stimulation]]$CD8 <<- list(IL2 = prior_IL2, IFNg = prior_IFNg, TNFa = prior_TNFa)
