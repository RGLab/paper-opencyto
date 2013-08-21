CREATE TABLE Antibody_E065
(
PK_Antibody_E065 INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Patient INTEGER REFERENCES Patient (PK_Patient),
AssayID TEXT,
Assay_Type TEXT,
Cell_Type TEXT,
Isolate TEXT,
Visit_Number TEXT,
Titer NUMERIC,
Response TEXT,
Filter_Flag TEXT,
Filter_Reason TEXT,
Draw_Date TEXT,
Test_Date TEXT
);

CREATE TABLE Antibody_S065
(
PK_Antibody_S065 INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Patient INTEGER REFERENCES Patient (PK_Patient),
AssayID TEXT,
Assay_Type TEXT,
Cell_Type TEXT,
Isolate TEXT,
Visit_Number TEXT,
Titer NUMERIC,
Response TEXT,
Filter_Flag TEXT,
Filter_Reason TEXT,
Draw_Date TEXT,
Test_Date TEXT
);

CREATE TABLE ELISPOT
(
PK_ELISPOT INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Patient INTEGER REFERENCES Patient (PK_Patient),
AssayID TEXT,
Antigen TEXT,
Visit_Number TEXT,
Delta NUMERIC,
Response TEXT,
Filter_Flag TEXT,
Filter_Reason TEXT,
Draw_Date TEXT,
Test_Date TEXT
);

CREATE TABLE Patient_FCS_Files
(
PK_Patient_FCS_Files INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Patient INTEGER REFERENCES Patient (PK_Patient),
Stimulation_Group TEXT,
Visit_Number TEXT,
FCS_Filename TEXT
);

CREATE TABLE Gating_Statistics
(
PK_Gating_Statistics INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Patient_FCS_Files INTEGER REFERENCES Patient_FCS_Files (PK_Patient_FCS_Files),
Marker TEXT,
Count INTEGER,
Proportion NUMERIC,
Parent TEXT,
Parent_Count INTEGER
);

CREATE TABLE Treatment
(
PK_Treatment INTEGER PRIMARY KEY AUTOINCREMENT,
Treatment_Group TEXT,
Group_Name TEXT,
Product1_Dose TEXT,
Product1_RX TEXT,
Product2_Dose TEXT,
Product2_RX TEXT
);

CREATE TABLE Patient
(
PK_Patient INTEGER PRIMARY KEY AUTOINCREMENT,
FK_Treatment_Group INTEGER REFERENCES Treatment (PK_Treatment),
PTID TEXT
);

