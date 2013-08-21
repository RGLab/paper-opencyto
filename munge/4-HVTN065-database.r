library(ProjectTemplate)
load.project()

db_creation_script <- "./munge/database-HVTN065.sql"
sql_db <- "./cache/db-HVTN065.sqlite"

if (!file.exists(sql_db)) {
  system(paste("cat", db_creation_script, "| sqlite3", sql_db))
}

# Opens connection to SQLite database
db_conn <- dbConnect(SQLite(), dbname = sql_db)

# Populates Treatment table
group1_treatment <- extract_treatment_groups(treatment.HVTN065, "Group 1")
group2_treatment <- extract_treatment_groups(treatment.HVTN065, "Group 2")
group3_treatment <- extract_treatment_groups(treatment.HVTN065, "Group 3")
group4_treatment <- extract_treatment_groups(treatment.HVTN065, "Group 4")
treatment_data <- rbind(group1_treatment, group2_treatment, group3_treatment, group4_treatment)

sql_insert <- sql_insert_from_df(df = treatment_data, table_name = "Treatment")
treatment_data$Treatment_Primary_Key <- bulk_insert(sql_insert, treatment_data, db_conn)

treatment.HVTN065 <- merge(treatment.HVTN065, treatment_data,
                           by.x = c("group", "product1_rx", "product1_dose", "product2_rx", "product2_dose"),
                           by.y = c("Group_Name", "Product1_RX", "Product1_Dose", "Product2_RX", "Product2_Dose"))

# Populates Patient table
patient_data <- data.frame(PTID = gsub("-", "", treatment.HVTN065$Ptid),
                           FK_Treatment_Group = treatment.HVTN065$Treatment_Primary_Key)
sql_insert <- sql_insert_from_df(df = patient_data, table_name = "Patient")
treatment.HVTN065$Patient_Primary_Key <- bulk_insert(sql_insert, patient_data, db_conn)

# Populates Patient_FCS_Files table
patient_data$Patient_Primary_Key <- treatment.HVTN065$Patient_Primary_Key[match(patient_data$PTID, gsub("-", "", treatment.HVTN065$Ptid))]
patient_FCS_data <- merge(patient_data, pData_HVTN065)
patient_FCS_data <- subset(patient_FCS_data, select = c(Patient_Primary_Key, Stim, VISITNO, name))
colnames(patient_FCS_data) <- c("FK_Patient", "Stimulation_Group", "Visit_Number", "FCS_Filename")
sql_insert <- sql_insert_from_df(df = patient_FCS_data, table_name = "Patient_FCS_Files")
patient_FCS_data$Primary_Key <- bulk_insert(sql_insert, patient_FCS_data, db_conn)


# Populates Gating_Statistics table
m_popstats <- melt(popstats_HVTN065, varnames = c("Population", "FCS_Filename"))
colnames(m_popstats)[3] <- "Proportion"

m_counts <- melt(counts_HVTN065, varnames = c("Population", "FCS_Filename"))
colnames(m_counts)[3] <- "Count"

gating_statistics <- merge(m_popstats, m_counts, sort = FALSE)
gating_statistics$Marker <- get_marker(gating_statistics$Population)
gating_statistics$Parent <- get_parent(gating_statistics$Population)

match_FCS_Files <- match(gating_statistics$FCS_Filename, patient_FCS_data$FCS_Filename)
gating_statistics$FK_Patient_FCS_Files <- patient_FCS_data$Primary_Key[match_FCS_Files]

# TODO: Determine parent's count...this will have to wait. It will take a bit to get right.
gating_statistics$Parent_Count <- 0

gating_statistics <- subset(gating_statistics, select = c(FK_Patient_FCS_Files, Marker, Count, Proportion, Parent, Parent_Count))

sql_insert <- sql_insert_from_df(df = gating_statistics, table_name = "Gating_Statistics")
gating_statistics$Primary_Key <- bulk_insert(sql_insert, gating_statistics, db_conn)


# Populates ELISPOT table
elispot_data <- subset(elispot065, select = c(ptid, assayid, antigen, visitno,
                                              delta, response, filt_flag, reason,
                                              drawdt, testdt))
colnames(elispot_data) <- c("PTID", "AssayID", "Antigen", "Visit_Number", "Delta",
                            "Response", "Filter_Flag", "Filter_Reason",
                            "Draw_Date", "Test_Date")
elispot_data$PTID <- gsub("-", "", elispot_data$PTID)
elispot_data$Filter_Flag <- !is.na(elispot_data$Filter_Flag)

match_Patient_PK <- match(elispot_data$PTID, patient_data$PTID)
elispot_Patient_FK <- patient_data$Patient_Primary_Key[match_Patient_PK]
elispot_data <- data.frame(FK_Patient = elispot_Patient_FK,
                           subset(elispot_data, select = -PTID))
sql_insert <- sql_insert_from_df(df = elispot_data, table_name = "ELISPOT")
elispot_data$Primary_Key <- bulk_insert(sql_insert, elispot_data, db_conn)

# Populates Antibody_E065 table
antibody_E065_data <- subset(neutralizing.antibody.e065,
                             select = c(ptid, assayid, assaytyp, celltype,
                                        isolate, visitno, titer2, response,
                                        filt_flag, reason, drawdt, testdt))

colnames(antibody_E065_data) <- c("PTID", "AssayID", "Assay_Type", "Cell_Type",
                                  "Isolate", "Visit_Number", "Titer", "Response",
                                  "Filter_Flag", "Filter_Reason", "Draw_Date",
                                  "Test_Date")
antibody_E065_data$PTID <- gsub("-", "", antibody_E065_data$PTID)
antibody_E065_data$Filter_Flag <- !is.na(antibody_E065_data$Filter_Flag)

match_Patient_PK <- match(antibody_E065_data$PTID, patient_data$PTID)
antibody_Patient_FK <- patient_data$Patient_Primary_Key[match_Patient_PK]
antibody_E065_data <- data.frame(FK_Patient = antibody_Patient_FK,
                           subset(antibody_E065_data, select = -PTID))

sql_insert <- sql_insert_from_df(df = antibody_E065_data, table_name = "Antibody_E065")
antibody_E065_data$Primary_Key <- bulk_insert(sql_insert, antibody_E065_data, db_conn)

# Populates Antibody_S065 table
antibody_S065_data <- subset(neutralizing.antibody.s065,
                             select = c(ptid, assayid, assaytyp, celltype,
                                        isolate, visitno, titer2, response,
                                        filt_flag, reason, drawdt, testdt))

colnames(antibody_S065_data) <- c("PTID", "AssayID", "Assay_Type", "Cell_Type",
                                  "Isolate", "Visit_Number", "Titer", "Response",
                                  "Filter_Flag", "Filter_Reason", "Draw_Date",
                                  "Test_Date")
antibody_S065_data$PTID <- gsub("-", "", antibody_S065_data$PTID)
antibody_S065_data$Filter_Flag <- !is.na(antibody_S065_data$Filter_Flag)

match_Patient_PK <- match(antibody_S065_data$PTID, patient_data$PTID)
antibody_Patient_FK <- patient_data$Patient_Primary_Key[match_Patient_PK]
antibody_S065_data <- data.frame(FK_Patient = antibody_Patient_FK,
                                 subset(antibody_S065_data, select = -PTID))

sql_insert <- sql_insert_from_df(df = antibody_S065_data, table_name = "Antibody_S065")
antibody_S065_data$Primary_Key <- bulk_insert(sql_insert, antibody_S065_data, db_conn)

# Closes connection to SQLite database
dbDisconnect(db_conn)
