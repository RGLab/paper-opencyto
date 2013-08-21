#' Function to generate SQL Insert statement from data.frame
#' 
#' An SQL Insert statement is generated from the column names of a data.frame.
#' This string is meant to be compatible with SQLite databases interfaced via
#' \code{\link{RSQLite}}.
#' 
#' @param df data.frame
#' @param table_name string containing the SQL table's name
#' @return a string containing the SQL Insert statement generated
sql_insert_from_df <- function(df, table_name) {
  paste("insert into",
        table_name,
        paste0("(", paste(colnames(df), collapse = ", "), ")"),
        "values",
        paste0("(", paste(rep("?", ncol(df)), collapse = ", "), ")"))
}

#' Bulk insert of data.frame into a database
#'
#' This function inserts multiple rows into the database from a data.frame.
#' The code is borrowed and modified from an example in the
#' \code{\link{RSQLite}} documentation.
#'
#' @param sql_insert the SQL insert statement
#' @param df a data.frame containing the records to insert
#' @param db_conn a database connection (aka the SQLite connection)
#' @param primary_keys logical value indicating the primary keys should be
#' returned for the records inserted
#' @return vector containing primary keys by default. Otherwise, \code{NULL}
bulk_insert <- function(sql_insert, df, db_conn, primary_keys = TRUE) {
  dbBeginTransaction(db_conn)
  dbGetPreparedQuery(db_conn, sql_insert, bind.data = df)
  dbCommit(db_conn)
  
  pk_values <- NULL
  if (primary_keys) {
    last_id <- dbGetQuery(db_conn, "SELECT last_insert_rowid()")[[1]]
    pk_values <- seq.int(to = last_id, length = nrow(df))
  }
  pk_values
}

#' From HVTN065 treatment data.frame, extracts and formats relevant data
#' 
#' @param treatment_data data.frame of HVTN065 treatment data
#' @param group_name string containing the group name
#' @return data.frame containing the relevant data
extract_treatment_groups <- function(treatment_data, group_name = "Group 1") {
  treatment_data <- subset(treatment_data, group == group_name,
                           select = c(Control, product1_rx, product1_dose, product2_rx, product2_dose))
  treatment_data <- treatment_data[!duplicated(treatment_data), ]
  treatment_groups <- ifelse(treatment_data$Control == 0, "Treatment", "Placebo")
  treatment_data <- subset(treatment_data, select = -Control)
  colnames(treatment_data) <- c("Product1_RX", "Product1_Dose", "Product2_RX", "Product2_Dose")
  data.frame(Treatment_Group = treatment_groups, Group_Name = group_name, treatment_data)
}

#' For a full cellular population path, we extract the population of interest
#' 
#' @param population_path string containing the full population path
#' @return string containing the population node
get_marker <- Vectorize(function(population_path) {
  population_path <- as.character(population_path)
  pop_split <- strsplit(population_path, "/", fixed = TRUE)[[1]]
  tail(pop_split, n = 1)
})

#' For a full cellular population path, we extract the parent node
#' 
#' If the root node is passed, then \code{NA} is returned.
#' 
#' @param population_path string containing the full population path
#' @return string containing the parent node
get_parent <- Vectorize(function(population_path) {
  population_path <- as.character(population_path)
  pop_split <- strsplit(population_path, "/", fixed = TRUE)[[1]]
  parent_node <- pop_split[length(pop_split) - 1]
  
  if (length(parent_node) == 0) {
    parent_node <- NA
  } else if (parent_node == "") {
    parent_node <- "root"
  }
  
  parent_node
})
