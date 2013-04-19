library(ProjectTemplate)
load.project()

num_cores <- 12

# The path of the root directory containing the FCS files
fcs_path <- '/shared/silo_researcher/Gottardo_R/gfinak_working/HVTN065/FACSData'

# Subdirectories containing FCS files
# BUG: The 'list.dirs' function returns only the full.names.
# We want only the folder name -- not the entire path.
fcs_subdirs <- list.dirs(fcs_path, full.names = FALSE, recursive = FALSE)
fcs_subdirs <- sapply(strsplit(fcs_subdirs, "/"), tail, n = 1)

# We traverse each of the subdirectories above and read the headers of each FCS
# file.  In particular, we care about the headers: 'EXPERIMENT NAME' and 'Sample
# Order' We extract these keywords and construct a data.frame summarizing both
# headers for all FCS files.
fcs_headers <- mclapply(fcs_subdirs, function(current_subdir) {
  subdir_fullpath <- file.path(fcs_path, current_subdir)
  fcs_subdir <- dir(subdir_fullpath, "fcs")

  # If no FCS files are in the directory, we return NULL.
  fcs_headers <- NULL
  if (length(fcs_subdir) > 0) {
    fcs_headers <- read.FCSheader(files = fcs_subdir, path = subdir_fullpath,
                                  keyword = c("EXPERIMENT NAME", "Sample Order","Stim"))
    fcs_headers <- data.frame(do.call(rbind, fcs_headers))
    fcs_headers$FCS_file <- file.path(subdir_fullpath, rownames(fcs_headers))
    fcs_headers
  }
  fcs_headers
}, mc.cores = num_cores)

# Remove the NULL entries.
fcs_headers <- fcs_headers[!sapply(fcs_headers, is.null)]

# Some of the folders have no "Sample Order" keyword. In this case, a column of
# NA values is stored. In this case, we remove them.
fcs_headers <- do.call(rbind, lapply(fcs_headers, na.omit))


#Filter out columns of analysis.plan that we don't need
#We are only interested  in mapping fcs files to ptids.
ap<-analysis.plan[,c("ASSAYID","PTID","SAMP_ORD","ANTIGEN","VISITNO")]
ap<-na.omit(unique(ap))
ap<-subset(ap,VISITNO%in%c("2","12"))

#data.table
library(data.table)
ap<-data.table(ap)
fcs_headers<-data.table(fcs_headers)
setkey(ap,SAMP_ORD,ASSAYID,ANTIGEN)
ap[,PTID:=factor(PTID)]
ap[,SAMP_ORD:=factor(SAMP_ORD)]
setnames(fcs_headers,c("ASSAYID","SAMP_ORD","ANTIGEN","FCS_file"))
setkey(fcs_headers,SAMP_ORD,ASSAYID,ANTIGEN)
merged<-merge(ap,fcs_headers,by=c("SAMP_ORD","ASSAYID","ANTIGEN"))
setkey(merged,PTID)
merged.sub<-merged[ANTIGEN%in%c("ENV-1-PTEG","GAG-1-PTEG","sebctrl","negctrl","POL-1-PTEG")]
merged.sub<-merged.sub[,.SD[length(FCS_file)==1],c("PTID,VISITNO,ANTIGEN")]
merged.sub<-merged.sub[,.SD[length(unique(VISITNO))==2],c("PTID,ANTIGEN")]
analysis_plan<-merged.sub
analysis_plan<-as.data.frame(lapply(analysis_plan,factor))
# Next, we determine which 'sebctrl' FCS files must be read. We do this by
# comparing the FCS headers from above with the corresponding information in the
# 'analysis.plan', which is based on a query Greg ran in LabKey.
# colnames(fcs_headers) <- c("ASSAYID", "SAMP_ORD", "FCS_file")
# analysis_plan<-merge(fcs_headers,ap,by=c("SAMP_ORD","ASSAYID"))
# analysis_plan<-as.data.frame(lapply(analysis_plan,factor))
# #analysis_plan <- plyr:::join(fcs_headers, ap, type="inner")
# analysis_plan <- subset(analysis_plan, ANTIGEN == "sebctrl")
# analysis_plan <- subset(analysis_plan, !is.na(PTID))
# analysis_plan <- subset(analysis_plan, PTID %in% unique(pData_HVTN065$PTID))
# analysis_plan <- subset(analysis_plan, VISITNO %in% c(2, 12))

# TODO: It may be smarter/faster to read headers in from subdirs having names
# in ASSAYID only.

# TODO: For some reason I have too many FCS files left over...1488 of them.
#   Is the join not working as expected?
#   Do I need to subset further?

# 24 files per patient. The 16 corresponds to the rows for the 2^3 = 8 cytokine
# combos and the 2 markers CD4 and CD8.
#table(analysis_plan$PTID) / 16

# Given that there are 2 visit numbers (2 and 12), this means there are 12 = 24
# / 2 files per patient visit. Is this right?

# To see if this is correct, let's examine the first patient.
#foo <- subset(analysis_plan, PTID == unique(analysis_plan$PTID)[1])

# There are 2 assays
#sum(table(foo$ASSAYID) > 0)




