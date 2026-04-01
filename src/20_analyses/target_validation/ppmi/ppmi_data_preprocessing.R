## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/ppmi/ppmi_data_preprocessing.R")

# ==============================================================================
# PPMI Project 9000 Plasma Proteomics Data Cleaning Script
# T.Nanasi 2026.02
# ==============================================================================
# Description:
#   This script processes PPMI Project 9000 plasma proteomics data and metadata
#   from the September 2025 data release. It performs nomenclature updates,
#   extracts measurements, and compiles comprehensive donor and sample-level
#   phenotypic information.
#
# Dataset (p9k) structure:
#   - dataProt: Protein expression matrix (NPX values)
#   - dataMiss: Missingness indicators
#   - dataQC: Quality control flags
#   - pheno.samp: Sample-level phenotype data
#   - pheno.don: Donor-level phenotype data
#   - nom: Updated nomenclature table
# ==============================================================================

output_filename <- glue::glue(
  "{get_config()$checkpoint_dir}/name={name}/version={version}/dataset.rds",
  name = "PPMI_p9000_Plasma_MS",
  version = 29985
)

# ==============================================================================
# CONFIGURATION: File Paths and Working Directory
# ==============================================================================
# IMPORTANT: Modify these paths for your environment before running this script.
# All file paths below are relative to the base_directory.

#####
# Data files obtained from PPMI via https://www.ppmi-info.org/access-data-specimens/download-data
#####
path_ppmi_curated <- "_DATA/PPMI_curated_data_cuts/Proteomic_Working_Group_Curated_Dataset(multi-project)/ppmi_project_9000_plasma.xlsx"
path_data_dict <- "_DATA/PPMI2509/Data___Databases/Data_Dictionary_-_Harmonized_15Sep2025.csv"
path_participant_status <- "_DATA/PPMI2509/Pheno/Subject_Characteristics/Patient_Status/Participant_Status_15Sep2025.csv"
path_demographics <- "_DATA/PPMI2509/Pheno/Subject_Characteristics/Subject_Demographics/Demographics_15Sep2025.csv"
path_socioeconomics <- "_DATA/PPMI2509/Pheno/Subject_Characteristics/Subject_Demographics/Socio-Economics_15Sep2025.csv"
path_pd_diagnosis <- "_DATA/PPMI2509/Pheno/Medical_History/Medical/PD_Diagnosis_History_15Sep2025.csv"
path_smoking_online <- "_DATA/PPMI2509/Pheno/Online/Medical/Smoking_History__Online__18Sep2025.csv"
path_smoking_found <- "_DATA/PPMI2509/Pheno/FOUND/Follow_Up_persons_w_Neurologic_Disease/FOUND_RFQ_Smoking_History_19Sep2025.csv"
path_alcohol_found <- "_DATA/PPMI2509/Pheno/FOUND/Follow_Up_persons_w_Neurologic_Disease/FOUND_RFQ_Alcohol_19Sep2025.csv"
path_medical_conditions <- "_DATA/PPMI2509/Pheno/Medical_History/Medical/Medical_Conditions_Log_15Sep2025.csv"
path_falls_freezing <- "_DATA/PPMI2509/Pheno/Medical_History/Medical/Determination_of_Freezing_and_Falls_15Sep2025.csv"
path_vital_signs <- "_DATA/PPMI2509/Pheno/Medical_History/Medical/Vital_Signs_15Sep2025.csv"
path_age_at_visit <- "_DATA/PPMI2509/Pheno/Subject_Characteristics/Subject_Demographics/Age_at_visit_15Sep2025.csv"
# ==============================================================================


# ==============================================================================
# LOAD RAW DATA: PPMI Project 9000 Proteomic Data
# ==============================================================================
# Read Excel file containing Results, Lookup table, and Notes sheets
{
  Xraw <- list()
  sheet_names <- openxlsx::getSheetNames(path_ppmi_curated)

  for (j in 1:length(sheet_names)) {
    Xraw[[j]] <- openxlsx::read.xlsx(path_ppmi_curated, sheet_names[j])
  }
  names(Xraw) <- sheet_names
}
# ==============================================================================


# ==============================================================================
# NOMENCLATURE CLEANUP AND UPDATE
# ==============================================================================
# Build OLINKID-level nomenclature table from LUT + master UniProt mapping
#
# Process:
#   1. Extract lookup table (LUT) from raw data
#   2. Fetch UniProt annotations from file
#   3. Apply manual token patches to UNIPROT IDs
#   4. Map UNIPROT to SYMBOL, ENTREZID, ENSEMBL using annotations
#   5. Handle non-bijective mappings by collapsing with '|'
#
# Output: nom (1 row per OLINKID)
#   - LUT fields kept at OLINKID level (e.g. PANEL, LABEL)
#   - Original LUT-derived fields preserved as *.orig
#   - UNIPROT updated with manual token patches
#   - SYMBOL / ENTREZID / ENSEMBL pulled from annotations_uniprot
# ==============================================================================
{
  # ----------------------------------------------------------------------------
  # Helper Functions
  # ----------------------------------------------------------------------------
  {
    drop_bad <- function(x) {
      x <- as.character(x)
      x <- x[!is.na(x)]
      x <- trimws(x)
      x <- x[x != ""]
      x <- x[x != "NA"]
      x
    }

    collapse_unique <- function(x) {
      x <- drop_bad(x)
      if (length(x) == 0) {
        return(NA_character_)
      }
      paste(sort(unique(x)), collapse = "|")
    }

    split_lut_field <- function(x) {
      x <- as.character(x)
      parts <- unlist(strsplit(x, ";|_", perl = TRUE))
      parts <- drop_bad(parts)
      collapse_unique(parts)
    }

    # this is a manual patch based on expert review
    patch_uniprot_token <- function(x) {
      x <- as.character(x)
      x <- ifelse(x == "NTproBNP", "P16860",
        ifelse(x == "O43521-2", "O43521", x)
      )
      x
    }

    assert_bijective <- function(df, key_col, val_col) {
      key <- as.character(df[[key_col]])
      val <- as.character(df[[val_col]])

      nuniq <- tapply(
        val,
        key,
        function(z) length(unique(z[!is.na(z)]))
      )

      if (any(nuniq > 1, na.rm = TRUE)) {
        bad <- names(nuniq)[which(nuniq > 1)]
        stop(val_col, " is not unique per ", key_col, ". Examples: ", paste(head(bad, 10), collapse = ", "))
      }
      TRUE
    }
  }

  # ----------------------------------------------------------------------------
  # Load Inputs
  # ----------------------------------------------------------------------------

  lookup_table <- Xraw$Lookup
  head(lookup_table)
  Hmisc::describe(lookup_table)

  annotations_uniprot <- glue::glue(
    "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
    name = "uniprot_links_updating",
    version = 29235
  ) |>
    arrow::open_dataset() |>
    dplyr::collect()

  # ----------------------------------------------------------------------------
  # Sanity Checks on Lookup Table
  # ----------------------------------------------------------------------------

  if (!("OLINKID" %in% colnames(lookup_table))) {
    stop("Lookup table is missing OLINKID column.")
  }
  if (max(table(lookup_table$OLINKID)) != 1) {
    stop("Lookup table OLINKID is not unique.")
  }

  lut_carry_cols <- c("PANEL", "LABEL")
  missing_cols <- setdiff(lut_carry_cols, colnames(lookup_table))
  if (length(missing_cols) > 0) {
    stop("Lookup table is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  for (cc in lut_carry_cols) {
    assert_bijective(lookup_table, "OLINKID", cc)
  }

  # ----------------------------------------------------------------------------
  # Build nom.orig (1 row per OLINKID)
  # ----------------------------------------------------------------------------

  nom.orig <- data.frame(
    OLINKID = lookup_table$OLINKID,
    PANEL = lookup_table$PANEL,
    LABEL = lookup_table$LABEL,
    SYMBOL = NA_character_,
    UNIPROT = NA_character_,
    stringsAsFactors = FALSE
  )

  for (k in seq_len(nrow(nom.orig))) {
    nom.orig$SYMBOL[k] <- split_lut_field(lookup_table$GENE_SYMBOL[k])
    nom.orig$UNIPROT[k] <- split_lut_field(lookup_table$UNIPROT[k])
  }

  # ----------------------------------------------------------------------------
  # Create nom with *.orig columns
  # ----------------------------------------------------------------------------

  nom <- data.frame(
    OLINKID = nom.orig$OLINKID,
    PANEL = nom.orig$PANEL,
    LABEL = nom.orig$LABEL,
    SYMBOL.orig = nom.orig$SYMBOL,
    UNIPROT.orig = nom.orig$UNIPROT,
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------------------------
  # Update UNIPROT: Apply manual token patches, then collapse back to '|'
  # ----------------------------------------------------------------------------

  uniprot_long <- tidyr::separate_rows(
    data.frame(
      OLINKID = nom$OLINKID,
      UNIPROT = nom$UNIPROT.orig,
      stringsAsFactors = FALSE
    ),
    UNIPROT,
    sep = "\\|"
  )

  uniprot_long$UNIPROT <- trimws(as.character(uniprot_long$UNIPROT))
  uniprot_long$UNIPROT <- patch_uniprot_token(uniprot_long$UNIPROT)
  uniprot_long$UNIPROT <- drop_bad(uniprot_long$UNIPROT)

  uniprot_aggregated <- dplyr::summarise(
    dplyr::group_by(uniprot_long, OLINKID),
    UNIPROT = collapse_unique(UNIPROT),
    .groups = "drop"
  )

  nom <- dplyr::left_join(nom, uniprot_aggregated, by = "OLINKID")

  # ----------------------------------------------------------------------------
  # Build Annotation Lookup: Combine uniprot_id + uniprot_id_alt
  # ----------------------------------------------------------------------------

  annotations_primary <- data.frame(
    uniprot_key = annotations_uniprot$uniprot_id,
    gene_symbol = annotations_uniprot$gene_symbol,
    entrez_gene_id = annotations_uniprot$entrez_gene_id,
    ensembl_gene_id = annotations_uniprot$ensembl_gene_id,
    stringsAsFactors = FALSE
  )

  annotations_alternate <- data.frame(
    uniprot_key = annotations_uniprot$uniprot_id_alt,
    gene_symbol = annotations_uniprot$gene_symbol,
    entrez_gene_id = annotations_uniprot$entrez_gene_id,
    ensembl_gene_id = annotations_uniprot$ensembl_gene_id,
    stringsAsFactors = FALSE
  )

  annotations_combined <- rbind(annotations_primary, annotations_alternate)
  annotations_combined$uniprot_key <- drop_bad(annotations_combined$uniprot_key)

  annotations_aggregated <- dplyr::summarise(
    dplyr::group_by(annotations_combined, uniprot_key),
    gene_symbol = collapse_unique(gene_symbol),
    entrez_gene_id = collapse_unique(entrez_gene_id),
    ensembl_gene_id = collapse_unique(ensembl_gene_id),
    .groups = "drop"
  )

  if (max(table(annotations_aggregated$uniprot_key)) > 1) {
    stop("annotations_aggregated is not unique by uniprot_key (unexpected).")
  }

  # ----------------------------------------------------------------------------
  # Map UNIPROT Tokens to Gene Annotations
  # ----------------------------------------------------------------------------

  nomenclature_uniprot_long <- tidyr::separate_rows(
    data.frame(
      OLINKID = nom$OLINKID,
      UNIPROT = nom$UNIPROT,
      stringsAsFactors = FALSE
    ),
    UNIPROT,
    sep = "\\|"
  )

  nomenclature_uniprot_long$UNIPROT <- drop_bad(trimws(as.character(nomenclature_uniprot_long$UNIPROT)))

  nomenclature_mapped <- dplyr::left_join(
    nomenclature_uniprot_long,
    annotations_aggregated,
    by = c("UNIPROT" = "uniprot_key"),
    relationship = "many-to-one"
  )

  if (nrow(nomenclature_mapped) != nrow(nomenclature_uniprot_long)) {
    stop("Join changed row count (expected many-to-one).")
  }

  nomenclature_final <- dplyr::summarise(
    dplyr::group_by(nomenclature_mapped, OLINKID),
    SYMBOL = collapse_unique(gene_symbol),
    ENTREZID = collapse_unique(entrez_gene_id),
    ENSEMBL = collapse_unique(ensembl_gene_id),
    .groups = "drop"
  )

  nom <- dplyr::left_join(nom, nomenclature_final, by = "OLINKID")

  # ----------------------------------------------------------------------------
  # Coverage and Tracking Statistics
  # ----------------------------------------------------------------------------

  unique_uniprot_tokens <- unique(nomenclature_uniprot_long$UNIPROT)
  unique_uniprot_tokens <- drop_bad(unique_uniprot_tokens)

  tracked_uniprots <- intersect(unique_uniprot_tokens, annotations_aggregated$uniprot_key)

  print(
    paste(
      "Original probes:", nrow(nom.orig),
      "--> UniProt rows:", length(unique_uniprot_tokens),
      "--> Unique UniProts:", length(unique(unique_uniprot_tokens)),
      "--> Tracked by master:", length(tracked_uniprots)
    )
  )

  if (length(tracked_uniprots) != length(unique(unique_uniprot_tokens))) {
    warning("Some UNIPROT tokens are not tracked by master mapping.")
  }

  # ----------------------------------------------------------------------------
  # Stable Column Order for nom
  # ----------------------------------------------------------------------------

  nom <- nom[, c(
    "OLINKID", "PANEL", "LABEL",
    "SYMBOL.orig", "UNIPROT.orig", "UNIPROT",
    "SYMBOL", "ENTREZID", "ENSEMBL"
  )]
}
# ==============================================================================


# ==============================================================================
# EXTRACT MEASUREMENTS: Protein Data, Missingness, and QC Flags
# ==============================================================================
# Parse Results sheet into three sub-tables:
#   - dataProt: NPX (Normalized Protein eXpression) values
#   - dataMiss: Missingness indicators
#   - dataQC: Quality control flags
#
# Process:
#   1. Extract sample metadata (PATNO, EVENT_ID)
#   2. Separate columns by suffix: _NPX, _MISS, _QC
#   3. Synchronize column order with nomenclature table (nom)
#   4. Convert to numeric and validate data integrity
# ==============================================================================
{
  results_raw <- Xraw$Results

  # Extract sample metadata (first 2 columns: PATNO, EVENT_ID)
  dataOther <- results_raw[, c(1:2)]
  dataOther$SampleId <- paste(dataOther$PATNO, dataOther$EVENT_ID, sep = "_")
  max(table(dataOther$SampleId)) == 1 # Verify samples are unique

  # Extract measurement columns (columns 3 to end)
  results_measurements <- results_raw[, 3:ncol(results_raw)]
  rownames(results_measurements) <- paste(dataOther$PATNO, dataOther$EVENT_ID, sep = "_")

  # Identify column indices by suffix
  idx_npx <- grep("_NPX", colnames(results_measurements))
  idx_miss <- grep("_MISS", colnames(results_measurements))
  idx_qc <- grep("_QC", colnames(results_measurements))

  # Verify column ordering is as expected (NPX, QC, MISS)
  identical(c(idx_npx, idx_qc, idx_miss), c(1:ncol(results_measurements))) == TRUE

  # Extract sub-tables (NOT overwriting, creating new variables)
  dataProt <- results_measurements[, idx_npx]
  dataMiss <- results_measurements[, idx_miss]
  dataQC <- results_measurements[, idx_qc]

  # Verify column labels match across the three sub-tables
  labels_prot <- strsplit(colnames(dataProt), "_")
  for (k in 1:length(labels_prot)) {
    labels_prot[[k]] <- paste(labels_prot[[1]][1:(length(labels_prot[[1]]) - 1)], collapse = "_")
  }
  labels_prot <- unlist(labels_prot)

  labels_miss <- strsplit(colnames(dataMiss), "_")
  for (k in 1:length(labels_miss)) {
    labels_miss[[k]] <- paste(labels_miss[[1]][1:(length(labels_miss[[1]]) - 1)], collapse = "_")
  }
  labels_miss <- unlist(labels_miss)

  labels_qc <- strsplit(colnames(dataQC), "_")
  for (k in 1:length(labels_qc)) {
    labels_qc[[k]] <- paste(labels_qc[[1]][1:(length(labels_qc[[1]]) - 1)], collapse = "_")
  }
  labels_qc <- unlist(labels_qc)

  identical(labels_prot, labels_miss)
  identical(labels_prot, labels_qc)

  # Verify all nomenclature labels are present
  identical(
    sort(colnames(dataProt)),
    sort(paste(nom$LABEL, "NPX", sep = "_"))
  ) == TRUE

  # NOTE: Column order mismatch detected at row 85
  # Synchronize column order with nom table

  label_index_prot <- data.frame(
    LABEL = colnames(dataProt),
    idx = 1:ncol(dataProt)
  )
  label_index_prot <- dplyr::left_join(
    data.frame(LABEL = paste(nom$LABEL, "NPX", sep = "_")),
    label_index_prot,
    by = "LABEL"
  )
  dataProt <- dataProt[, label_index_prot$idx]

  label_index_miss <- data.frame(
    LABEL = colnames(dataMiss),
    idx = 1:ncol(dataMiss)
  )
  label_index_miss <- dplyr::left_join(
    data.frame(LABEL = paste(nom$LABEL, "MISS", sep = "_")),
    label_index_miss,
    by = "LABEL"
  )
  dataMiss <- dataMiss[, label_index_miss$idx]

  label_index_qc <- data.frame(
    LABEL = colnames(dataQC),
    idx = 1:ncol(dataQC)
  )
  label_index_qc <- dplyr::left_join(
    data.frame(LABEL = paste(nom$LABEL, "QC", sep = "_")),
    label_index_qc,
    by = "LABEL"
  )
  dataQC <- dataQC[, label_index_qc$idx]

  # Verify synchronization succeeded
  identical(colnames(dataProt), paste(nom$LABEL, "NPX", sep = "_"))
  identical(colnames(dataMiss), paste(nom$LABEL, "MISS", sep = "_"))
  identical(colnames(dataQC), paste(nom$LABEL, "QC", sep = "_"))

  # Convert to numeric and validate
  for (k in 1:ncol(dataProt)) {
    dataProt[, k] <- as.numeric(dataProt[, k])
  }
  identical(which(!is.finite(as.matrix(dataProt))), which(is.na(dataProt))) == TRUE

  for (k in 1:ncol(dataMiss)) {
    dataMiss[, k] <- as.numeric(dataMiss[, k])
  }
  identical(which(!is.finite(as.matrix(dataMiss))), which(is.na(dataMiss))) == TRUE

  table(unlist(c(dataQC)))
  sum(is.na(unlist(c(dataQC))))
  identical(which(is.na(dataQC)), which(is.na(dataMiss))) == TRUE
}
# ==============================================================================


# ==============================================================================
# COMPILE METADATA: Donor and Sample Phenotypes
# ==============================================================================
# Build comprehensive phenotype tables:
#   - pheno.don: Donor-level (one row per patient)
#   - pheno.samp: Sample-level (one row per visit/timepoint)
#
# Data sources:
#   - Enrollment status and genetic mutations
#   - Demographics (sex, race, handedness, birthday)
#   - Socioeconomic (education)
#   - Clinical (PD diagnosis, symptoms, falls)
#   - Lifestyle (smoking, alcohol)
#   - Medical conditions (hypertension, diabetes, etc.)
#   - Vital signs (BMI, blood pressure, heart rate)
# ==============================================================================
{
  # Read data dictionary
  variable_dictionary <- read.csv(path_data_dict)

  # ----------------------------------------------------------------------------
  # Initialize Phenotype Tables: Enrollment Status and Cohort Assignment
  # ----------------------------------------------------------------------------
  {
    participant_status_data <- read.csv(path_participant_status)

    # Filter to relevant columns
    participant_status_filtered <- dplyr::select(
      participant_status_data,
      c("PATNO", "ENROLL_STATUS", "ENRLGBA", "ENRLLRRK2", "ENRLPINK1", "ENRLPRKN", "ENRLSNCA")
    )
    participant_status_filtered$PATNO <- as.character(participant_status_filtered$PATNO)

    # Create enrollment and mutation dataframe
    enrollment_mutations <- data.frame(
      PATNO = participant_status_filtered$PATNO,
      ENROLL_STATUS = participant_status_filtered$ENROLL_STATUS,
      mutation_GBA = participant_status_filtered$ENRLGBA,
      mutation_LRRK2 = participant_status_filtered$ENRLLRRK2,
      mutation_PINK1 = participant_status_filtered$ENRLPINK1,
      mutation_PARKIN = participant_status_filtered$ENRLPRKN,
      mutation_SNCA = participant_status_filtered$ENRLSNCA
    )

    # Initialize sample-level phenotype table
    pheno.samp <- dataOther
    pheno.samp$cohort <- "removed"

    # Initialize donor-level phenotype table
    pheno.don <- data.frame(
      PATNO = names(table(dataOther$PATNO)),
      n_samples = as.numeric(table(dataOther$PATNO))
    )
    pheno.don$cohort <- "removed"

    # Merge enrollment/mutation data
    pheno.don <- dplyr::left_join(pheno.don, enrollment_mutations, by = "PATNO")

    # Assign cohort labels (PD, control, prodromal)
    {
      cohort_data <- read.csv(path_participant_status)
      mean(dataOther$PATNO %in% cohort_data$PATNO) == 1 # Verify all samples present
      table(cohort_data$COHORT, cohort_data$COHORT_DEFINITION)
      max(table(cohort_data$PATNO)) == 1 # Verify unique patients

      donors_pd <- cohort_data$PATNO[which(cohort_data$COHORT_DEFINITION == "Parkinson's Disease")]
      donors_hc <- cohort_data$PATNO[which(cohort_data$COHORT_DEFINITION == "Healthy Control")]
      donors_prodromal <- cohort_data$PATNO[which(cohort_data$COHORT_DEFINITION == "Prodromal")]

      pheno.don$cohort[which(pheno.don$PATNO %in% donors_pd)] <- "PD"
      pheno.don$cohort[which(pheno.don$PATNO %in% donors_hc)] <- "control"
      pheno.don$cohort[which(pheno.don$PATNO %in% donors_prodromal)] <- "prodromal"

      table(pheno.don$cohort, pheno.don$ENROLL_STATUS)

      # Copy cohort labels to sample table
      pheno.samp$cohort[which(pheno.samp$PATNO %in% donors_pd)] <- "PD"
      pheno.samp$cohort[which(pheno.samp$PATNO %in% donors_hc)] <- "control"
      pheno.samp$cohort[which(pheno.samp$PATNO %in% donors_prodromal)] <- "prodromal"

      table(pheno.samp$cohort)
      table(pheno.don$cohort)
    }
  }

  # ----------------------------------------------------------------------------
  # Donor-Specific Features
  # ----------------------------------------------------------------------------
  {
    # --------------------------------------------------------------------------
    # Demographics: Birthday, Sex, Childbearing Potential, Handedness, Race
    # --------------------------------------------------------------------------
    {
      demographics_data <- read.csv(path_demographics)
      mean(dataOther$PATNO %in% demographics_data$PATNO) == 1
      max(table(demographics_data$PATNO)) == 1 # 1 entry per patient
      demographics_data$PATNO <- as.character(demographics_data$PATNO)

      demographics_merged <- dplyr::left_join(pheno.don, demographics_data, by = "PATNO")
      identical(pheno.don$PATNO, demographics_merged$PATNO)

      # Sex: 0 = female, 1 = male
      pheno.don$sex <- NA
      pheno.don$sex[which(pheno.don$PATNO %in% demographics_data$PATNO[which(demographics_data$SEX == 0)])] <- "F"
      pheno.don$sex[which(pheno.don$PATNO %in% demographics_data$PATNO[which(demographics_data$SEX == 1)])] <- "M"

      # Childbearing potential
      pheno.don$childbearing_potential <- demographics_merged$CHLDBEAR

      # Handedness
      pheno.don$handedness <- demographics_merged$HANDED

      # Birthday (convert to months since 1900, then to year)
      pheno.don$t_birth <- (as.numeric(sapply(strsplit(demographics_merged$BIRTHDT, "/"), "[[", 2)) - 1900) * 12 +
        as.numeric(sapply(strsplit(demographics_merged$BIRTHDT, "/"), "[[", 1))
      pheno.don$birth <- pheno.don$t_birth / 12 + 1900

      # Race categorization
      {
        race_white_cols <- c("RAWHITE", "ASHKJEW", "BASQUE", "AFRICBERB", "HISPLAT")
        race_nonwhite_cols <- c("RABLACK", "RAASIAN", "RAHAWOPI", "RAINDALS")

        patients_white <- demographics_merged$PATNO[rowSums(demographics_merged[, colnames(demographics_merged) %in% race_white_cols], na.rm = T) > 0]
        patients_nonwhite <- demographics_merged$PATNO[rowSums(demographics_merged[, colnames(demographics_merged) %in% race_nonwhite_cols], na.rm = T) > 0]

        pheno.don$race <- "unknown"
        pheno.don$race[which(pheno.don$PATNO %in% patients_white)] <- "caucasian"
        pheno.don$race[which(pheno.don$PATNO %in% patients_nonwhite)] <- "non_caucasian"
        pheno.don$race[which(pheno.don$PATNO %in% intersect(patients_white, patients_nonwhite))] <- "mixed"
      }
    }

    # --------------------------------------------------------------------------
    # Education: Years of Education
    # --------------------------------------------------------------------------
    {
      socioeconomics_data <- read.csv(path_socioeconomics)
      mean(dataOther$PATNO %in% socioeconomics_data$PATNO) == 1
      socioeconomics_data$PATNO <- as.character(socioeconomics_data$PATNO)

      pheno.don$education_years <- NA
      for (k in 1:nrow(pheno.don)) {
        patient_records <- which(socioeconomics_data$PATNO == pheno.don$PATNO[k])
        if (length(patient_records) > 0) {
          pheno.don$education_years[k] <- max(socioeconomics_data$EDUCYRS[patient_records], na.rm = T)
        }
      }
    }

    # --------------------------------------------------------------------------
    # PD Diagnosis and First Symptom Date
    # --------------------------------------------------------------------------
    {
      pd_diagnosis_data <- read.csv(path_pd_diagnosis)
      mean(pheno.don$PATNO[which(pheno.don$cohort == "PD")] %in% pd_diagnosis_data$PATNO) == 1
      max(table(pd_diagnosis_data$PATNO)) == 1
      pd_diagnosis_data$PATNO <- as.character(pd_diagnosis_data$PATNO)

      pd_diagnosis_merged <- dplyr::left_join(pheno.don, pd_diagnosis_data, by = "PATNO")
      identical(pheno.don$PATNO, pd_diagnosis_merged$PATNO)

      # PD diagnosis date
      pheno.don$t_pddiag <- NA
      idx_pddiag <- which(!is.na(pd_diagnosis_merged$PDDXDT))
      pheno.don$t_pddiag[idx_pddiag] <- (as.numeric(sapply(strsplit(pd_diagnosis_merged$PDDXDT[idx_pddiag], "/"), "[[", 2)) - 1900) * 12 +
        as.numeric(sapply(strsplit(pd_diagnosis_merged$PDDXDT[idx_pddiag], "/"), "[[", 1))
      pheno.don$pddiag <- pheno.don$t_pddiag / 12 + 1900

      # First symptom date
      pheno.don$t_pdsymptom <- NA
      idx_pdsymptom <- which(!is.na(pd_diagnosis_merged$SXDT))
      pheno.don$t_pdsymptom[idx_pdsymptom] <- (as.numeric(sapply(strsplit(pd_diagnosis_merged$SXDT[idx_pdsymptom], "/"), "[[", 2)) - 1900) * 12 +
        as.numeric(sapply(strsplit(pd_diagnosis_merged$SXDT[idx_pdsymptom], "/"), "[[", 1))
      pheno.don$pdsymptom <- pheno.don$t_pdsymptom / 12 + 1900

      table(is.na(pheno.don$pdsymptom), is.na(pheno.don$pddiag))
      barplot(sort(pheno.don$t_pdsymptom - pheno.don$t_pddiag))
      mean(pheno.don$t_pdsymptom - pheno.don$t_pddiag, na.rm = T) / 12 # Symptoms predate diagnosis by ~1.7 years

      # Derived variables
      pheno.don$age_at_pd_diag <- (pheno.don$t_pddiag - pheno.don$t_birth) / 12
      pheno.don$early_onset_pd <- pheno.don$age_at_pd_diag < 50 # Early onset defined as age < 50
    }

    # --------------------------------------------------------------------------
    # Smoking History
    # --------------------------------------------------------------------------
    {
      # Online smoking data
      smoking_online_data <- read.csv(path_smoking_online)

      smoking_online_ever <- rep(NA, nrow(pheno.don))
      smoking_online_now <- smoking_online_ever

      variables_smoking_ever <- c(
        "FIVE_PACKS_IN_LIFE_OL", "SMOKE_REG_IN_LIFE_OL",
        "SMOKE_PIPESCIGARS_REG_IN_LIFE_OL", "SMOKE_SMKLESS_REG_IN_LIFE_OL",
        "SMOKINGHX_ECIG_REG_IN_LIFE_OL"
      )
      variables_smoking_now <- c(
        "SMOKE_CIG_CURR_OL", "SMOKE_PIPESCIGARS_CURR_OL",
        "SMOKINGHX_SMKLESS_CURR_OL", "SMOKINGHX_ECIG_CURR_OL"
      )

      for (k in 1:nrow(pheno.don)) {
        patient_records <- which(smoking_online_data$PATNO == pheno.don$PATNO[k])
        if (length(patient_records) > 0) {
          patient_data <- smoking_online_data[patient_records, ]

          values_ever <- unlist(patient_data[, which(colnames(patient_data) %in% variables_smoking_ever)])
          values_ever <- values_ever[!is.na(values_ever)]
          if (length(values_ever) > 0 & mean(values_ever != 1) == 1) {
            smoking_online_ever[k] <- 0
          }
          if (length(values_ever) > 0 & sum(values_ever == 1) > 0) {
            smoking_online_ever[k] <- 1
          }

          values_now <- unlist(patient_data[, which(colnames(patient_data) %in% variables_smoking_now)])
          values_now <- values_now[!is.na(values_now)]
          if (length(values_now) > 0 & mean(values_now != 1) == 1) {
            smoking_online_now[k] <- 0
          }
          if (length(values_now) > 0 & sum(values_now == 1) > 0) {
            smoking_online_now[k] <- 1
          }
        }
      }

      # FOUND smoking data
      smoking_found_data <- read.csv(path_smoking_found)
      smoking_found_ever_data <- dplyr::select(smoking_found_data, c("patno", "smq1", "smq2", "smq7"))
      smoking_found_now_data <- dplyr::select(smoking_found_data, c("patno", "smq4", "smq9"))

      colnames(smoking_found_ever_data)[1] <- "PATNO"
      smoking_found_ever_data$PATNO <- as.character(smoking_found_ever_data$PATNO)
      colnames(smoking_found_now_data)[1] <- "PATNO"
      smoking_found_now_data$PATNO <- as.character(smoking_found_now_data$PATNO)

      smoking_found_ever_merged <- dplyr::left_join(
        data.frame(PATNO = pheno.don$PATNO),
        smoking_found_ever_data,
        by = "PATNO"
      )
      smoking_found_now_merged <- dplyr::left_join(
        data.frame(PATNO = pheno.don$PATNO),
        smoking_found_now_data,
        by = "PATNO"
      )

      smoking_found_ever_matrix <- as.matrix(smoking_found_ever_merged[, 2:ncol(smoking_found_ever_merged)])
      smoking_found_now_matrix <- as.matrix(smoking_found_now_merged[, 2:ncol(smoking_found_now_merged)])
      smoking_found_ever_matrix[smoking_found_ever_matrix > 1] <- NA
      smoking_found_now_matrix[smoking_found_now_matrix > 1] <- NA

      smoking_found_ever <- rep(NA, nrow(pheno.don))
      smoking_found_now <- smoking_found_ever
      smoking_found_ever[which(rowMeans(is.na(smoking_found_ever_matrix)) != 1 & rowSums(smoking_found_ever_matrix == 1, na.rm = T) > 0)] <- 1
      smoking_found_ever[which(rowMeans(is.na(smoking_found_ever_matrix)) != 1 & rowSums(smoking_found_ever_matrix == 1, na.rm = T) == 0)] <- 0
      smoking_found_now[which(rowMeans(is.na(smoking_found_now_matrix)) != 1 & rowSums(smoking_found_now_matrix == 1, na.rm = T) > 0)] <- 1
      smoking_found_now[which(rowMeans(is.na(smoking_found_now_matrix)) != 1 & rowSums(smoking_found_now_matrix == 1, na.rm = T) == 0)] <- 0

      # Merge smoking data from both sources
      pheno.don$smoking <- NA
      pheno.don$smoking[which(smoking_online_ever == 0 | smoking_found_ever == 0)] <- "never"
      pheno.don$smoking[which(smoking_online_ever == 1 | smoking_found_ever == 1)] <- "past"
      pheno.don$smoking[which(smoking_online_now == 1 | smoking_found_now == 1)] <- "now"
    }

    # --------------------------------------------------------------------------
    # Alcohol Consumption
    # --------------------------------------------------------------------------
    {
      alcohol_found_data <- read.csv(path_alcohol_found)
      alcohol_found_selected <- dplyr::select(alcohol_found_data, c("patno", "alq1", "alq2", "alq4"))
      colnames(alcohol_found_selected)[1] <- "PATNO"
      alcohol_found_selected$PATNO <- as.character(alcohol_found_selected$PATNO)

      alcohol_found_merged <- dplyr::left_join(
        data.frame(PATNO = pheno.don$PATNO),
        alcohol_found_selected,
        by = "PATNO"
      )
      alcohol_found_matrix <- as.matrix(alcohol_found_merged[, 2:ncol(alcohol_found_merged)])
      alcohol_found_matrix[alcohol_found_matrix > 1] <- NA

      # Check for conflicts
      if (sum(alcohol_found_matrix[, 3] == 1 & alcohol_found_matrix[, 1] == 0, na.rm = T) != 0) {
        stop("Alcohol conflict detected")
      }
      if (sum(alcohol_found_matrix[, 3] == 1 & alcohol_found_matrix[, 2] == 0, na.rm = T) != 0) {
        stop("Alcohol conflict detected")
      }

      alcohol_status <- rep(NA, nrow(pheno.don))
      alcohol_status[which(pmax(alcohol_found_matrix[, 1], alcohol_found_matrix[, 2], na.rm = T) == 0)] <- "never"
      alcohol_status[which(pmax(alcohol_found_matrix[, 1], alcohol_found_matrix[, 2], na.rm = T) > 0)] <- "past"
      alcohol_status[which(alcohol_found_matrix[, 3] > 0)] <- "now"

      pheno.don$alcohol <- alcohol_status
    }

    # --------------------------------------------------------------------------
    # Medical Conditions: Hypertension, Diabetes, Metabolic, Psych, Cognitive
    # --------------------------------------------------------------------------
    {
      medical_conditions_data <- read.csv(path_medical_conditions)
      medical_conditions_data$PATNO <- as.character(medical_conditions_data$PATNO)
      medical_conditions_filtered <- medical_conditions_data[which(medical_conditions_data$PATNO %in% pheno.don$PATNO), ]

      all_medical_terms <- data.frame(
        term = names(table(tolower(stringr::str_trim(medical_conditions_filtered$MHTERM)))),
        count = as.numeric(table(tolower(stringr::str_trim(medical_conditions_filtered$MHTERM))))
      )
      all_medical_terms <- all_medical_terms[order(all_medical_terms$term), ]

      # Define medical condition search terms
      terms_cognitive_deficit <- c("alzheimer", "cognitive", "dementia")
      terms_psychiatric <- c("anxiety", "depressed", "depression", "depressive", "panic", "bipolar")
      terms_hypertension <- c("hypertension", "hypertonia", "blood pressure")
      terms_diabetes <- c("diabetes", "diabetic", "diabetis")
      terms_metabolic <- c(
        "cholesterol", "hypercholesterinemia", "hypercholesterolaemia", "hypercholesterolemia",
        "hypercholesterolemia/hyperlipidemia", "hypercholesterolenia", "hyperlipemia",
        "hyperlipidemia", "hypertriglyceridemia"
      )

      # Build medical condition indicators
      medical_condition_keys <- list()
      medical_condition_keys[[1]] <- terms_hypertension
      medical_condition_keys[[2]] <- terms_diabetes
      medical_condition_keys[[3]] <- terms_metabolic
      medical_condition_keys[[4]] <- terms_psychiatric
      medical_condition_keys[[5]] <- terms_cognitive_deficit

      medical_conditions_matrix <- as.data.frame(matrix(NA, nrow = nrow(pheno.don), ncol = 5))
      colnames(medical_conditions_matrix) <- c("hypertension", "diabetes", "lipid_cholesterol", "depression_anxiety", "cognitive_deficit")
      medical_conditions_matrix[which(pheno.don$PATNO %in% medical_conditions_filtered$PATNO), ] <- 0

      medical_terms_lower <- tolower(stringr::str_trim(medical_conditions_filtered$MHTERM))

      for (k in 1:nrow(pheno.don)) {
        patient_records <- which(medical_conditions_filtered$PATNO == pheno.don$PATNO[k])
        if (length(patient_records) > 0) {
          for (w in 1:length(medical_condition_keys)) {
            medical_conditions_matrix[k, w] <- as.numeric(
              length(grep(paste(medical_condition_keys[[w]], collapse = "|"), medical_terms_lower[patient_records])) > 0
            )
          }
        }
      }

      pheno.don <- cbind(pheno.don, medical_conditions_matrix)

      table(pheno.don$hypertension, pheno.don$cohort)
      table(pheno.don$diabetes, pheno.don$cohort)
      table(pheno.don$lipid_cholesterol, pheno.don$cohort)
      table(pheno.don$depression_anxiety, pheno.don$cohort)
      table(pheno.don$cognitive_deficit, pheno.don$cohort)
    }

    # --------------------------------------------------------------------------
    # Falls and Freezing Events
    # --------------------------------------------------------------------------
    {
      # Manual curation of fall terms from medical conditions log
      # Exclude fall risk, sports injuries, and ladder/horse accidents
      {
        falls_medical_conditions <- read.csv(path_medical_conditions)
        falls_medical_conditions$PATNO <- as.character(falls_medical_conditions$PATNO)

        falls_terms_all <- data.frame(
          term = unique(falls_medical_conditions$MHTERM[grep("fall", tolower(falls_medical_conditions$MHTERM))])
        )

        falls_terms_relevant <- c(
          "following a fall", "from fall", "due to fall", "from a fall", "after fall",
          "after a fall", "due to falling", "accidental fall", "after falling",
          "falling due to PD", "Fall resulting in stitches", "falls (PD related)",
          "falls (PD meds affecting balance)", "Fall with possible concussion",
          "Fall (Tripped and fell)", "Severe Back Pain (due to fall)", "Lesion post fall",
          "resulted from fall", "Frequent Falls", "Intermittent falls", "A fall(accidental?)",
          "Fall-related", "Bad fall", "Recurring Falls", "Resulting from fall",
          "Fall - Cracked 3 ribs", "Fall with cervical fracture", "d/t Fall",
          "same fall", "Fall - bruise"
        )

        falls_terms_matched <- falls_terms_all[grep(paste(tolower(falls_terms_relevant), collapse = "|"), tolower(falls_terms_all$term)), ]

        falls_terms_excluded <- c(
          "LWK fractures due to fall from ladder",
          "An accidental fall when playing pickleball, resulting in a hairline fracture of the nose, soft tissue injury of the knees, and damage to breast implants(s/p breast cancer sugery) that required removal ",
          "An accidental fall on ice causing a right ankle fracture, treated non-surgically. ",
          "Epilepsy due to fall on head",
          "Left wrist fracture caused by an accidental fall during dancing, s/p in a cast.  ",
          "Head Trauma (Concussion after fall from roof)"
        )

        falls_terms_ok <- setdiff(falls_terms_matched, falls_terms_excluded)

        falls_terms_ok <- c(falls_terms_ok, c(
          "Frequency of falls (PD related) ",
          "Frequency of falls (PD meds affecting balance)",
          "Fall (Tripped and fell)",
          "A fall(accidental?) resulting in left wrist fracture, in casting",
          "Fall",
          "Subarachnoid hemorrhage due to a fall ",
          "R leg pain - s/p fall"
        ))

        # Add back ice and dance falls
        falls_terms_ok <- c(falls_terms_ok, c(
          "Left wrist contusion/fracture from a slip and fall on ice; treated with immobilization followed by physical therapy",
          "Parkinson's disease Symptoms began in 2018, pt noticed that my hand writing became difficult December 2019/January 2020, left hand occasionally shaking, hand writing is difficult. In 2021, left hand lost some flexibility and lost some strength.  Left hand shaking increase frequency.  However, pt reports that she can do everything by herself. In November 2021,  pt reported noticing right hand had lost some strength.   Sleeping pattern has changed in last two years. Pt feels it is difficult to fall asleep during night.  Now pt takes melatonin 5 mg per day at 10:30 pm.  Sleep improved to some degree but some days pt still cannot fall in sleep",
          "An accidental fall on ice causing a right ankle fracture, treated non-surgically. ",
          "Left wrist fracture caused by an accidental fall during dancing, s/p in a cast.  "
        ))

        falls_medical_filtered <- falls_medical_conditions[which(falls_medical_conditions$MHTERM %in% falls_terms_ok), ]
        falls_medical_filtered$t_fall_medcond <- NA # Months since 1900

        for (k in 1:nrow(falls_medical_filtered)) {
          date_parts <- unlist(strsplit(falls_medical_filtered$MHDIAGDT[k], "/"))
          if (length(date_parts) > 1) {
            falls_medical_filtered$t_fall_medcond[k] <- (as.numeric(date_parts[2]) - 1900) * 12 + as.numeric(date_parts[1])
          }
        }

        first_fall_medical_conditions <- rep(NA, nrow(pheno.don))
        for (k in 1:nrow(pheno.don)) {
          patient_records <- which(falls_medical_filtered$PATNO == pheno.don$PATNO[k])
          if (length(patient_records) > 0) {
            first_fall_medical_conditions[k] <- min(falls_medical_filtered$t_fall_medcond[patient_records], na.rm = T)
          }
        }
      }

      # Processed falls and freezing data
      {
        falls_freezing_data <- read.csv(path_falls_freezing)
        falls_freezing_data$PATNO <- as.character(falls_freezing_data$PATNO)

        matrixStats::colMaxs(as.matrix(falls_freezing_data[6:(ncol(falls_freezing_data) - 2)]), na.rm = T)

        # Index-based selection: columns 7 to (ncol-2) contain fall/freeze indicators
        # NOTE: Kept as index-based for stability, as column names may vary
        falls_freezing_sum <- rowSums(as.matrix(falls_freezing_data[, 7:(ncol(falls_freezing_data) - 2)]), na.rm = T)
        falls_freezing_sum[which(rowMeans(is.na(falls_freezing_data[, 7:(ncol(falls_freezing_data) - 2)])) == 1)] <- NA

        falls_freezing_data$sum.freezing.falls <- falls_freezing_sum

        falls_freezing_filtered <- falls_freezing_data[which(falls_freezing_data$PATNO %in% pheno.don$PATNO), ]
        falls_freezing_date_check <- falls_freezing_filtered[which(falls_freezing_filtered$INFODT != falls_freezing_filtered$ORIG_ENTRY), ]
        falls_freezing_date_check <- dplyr::select(falls_freezing_date_check, c("INFODT", "ORIG_ENTRY"))
        falls_freezing_date_check$t_INFODT <- as.numeric(sapply(sapply(falls_freezing_date_check$INFODT, strsplit, "/"), "[[", 1)) +
          as.numeric(sapply(sapply(falls_freezing_date_check$INFODT, strsplit, "/"), "[[", 2)) * 12
        falls_freezing_date_check$t_ORIG_ENTRY <- as.numeric(sapply(sapply(falls_freezing_date_check$ORIG_ENTRY, strsplit, "/"), "[[", 1)) +
          as.numeric(sapply(sapply(falls_freezing_date_check$ORIG_ENTRY, strsplit, "/"), "[[", 2)) * 12
        range(falls_freezing_date_check$t_ORIG_ENTRY - falls_freezing_date_check$t_INFODT) # INFODT is earlier

        falls_freezing_filtered$t_entry <- as.numeric(sapply(sapply(falls_freezing_filtered$INFODT, strsplit, "/"), "[[", 1)) +
          (as.numeric(sapply(sapply(falls_freezing_filtered$INFODT, strsplit, "/"), "[[", 2)) - 1900) * 12

        first_freezing_fall <- rep(NA, nrow(pheno.don))
        for (k in 1:nrow(pheno.don)) {
          patient_records <- which(falls_freezing_filtered$PATNO == pheno.don$PATNO[k])
          if (length(patient_records) > 0) {
            first_freezing_fall[k] <- min(falls_freezing_filtered$t_entry[patient_records], na.rm = T)
          }
        }
      }

      # Combine medical conditions falls and processed falls
      pheno.don$t_first_fall_or_freezing <- pmin(first_fall_medical_conditions, first_freezing_fall, na.rm = T)
      pheno.don$first_fall_or_freezing <- pheno.don$t_first_fall_or_freezing / 12 + 1900
    }
  }

  # ----------------------------------------------------------------------------
  # Sample-Specific Features
  # ----------------------------------------------------------------------------
  {
    # Merge static donor-level times into sample table
    pheno.samp <- dplyr::left_join(
      pheno.samp,
      dplyr::select(pheno.don, c("PATNO", "t_birth", "t_pddiag", "t_pdsymptom", "t_first_fall_or_freezing")),
      by = "PATNO"
    )

    # --------------------------------------------------------------------------
    # Vital Signs: BMI, Blood Pressure, Body Temperature
    # --------------------------------------------------------------------------
    {
      vital_signs_data <- read.csv(path_vital_signs)
      vital_signs_data$SampleId <- paste(vital_signs_data$PATNO, vital_signs_data$EVENT_ID, sep = "_")
      mean(dataOther$SampleId %in% vital_signs_data$SampleId) # Missing only one prodromal sample
      dataOther[which(!dataOther$SampleId %in% vital_signs_data$SampleId), ]

      vital_signs_merged <- dplyr::left_join(pheno.samp, vital_signs_data, by = "SampleId")
      identical(vital_signs_merged$SampleId, pheno.samp$SampleId) == T

      # Visit date (INFODT)
      pheno.samp$t_visit <- NA
      idx_visit <- which(!is.na(vital_signs_merged$INFODT))
      pheno.samp$t_visit[idx_visit] <- (as.numeric(sapply(strsplit(vital_signs_merged$INFODT[idx_visit], "/"), "[[", 2)) - 1900) * 12 +
        as.numeric(sapply(strsplit(vital_signs_merged$INFODT[idx_visit], "/"), "[[", 1))
      pheno.samp$visit_date <- pheno.samp$t_visit / 12 + 1900

      # Derived temporal variables
      pheno.samp$t_age <- pheno.samp$t_visit - pheno.samp$t_birth
      pheno.samp$age_at_visit <- pheno.samp$t_age / 12
      pheno.samp$months_since_PD_diagnosis <- pheno.samp$t_visit - pheno.samp$t_pddiag
      pheno.samp$months_since_PD_symptomps <- pheno.samp$t_visit - pheno.samp$t_pdsymptom
      pheno.samp$months_since_first_fall_or_freezing <- pheno.samp$t_visit - pheno.samp$t_first_fall_or_freezing

      # Basic vitals
      pheno.samp$body_temperature <- vital_signs_merged$TEMPC
      pheno.samp$blood_pressure_systolic <- vital_signs_merged$SYSSUP
      pheno.samp$blood_pressure_diastolic <- vital_signs_merged$DIASUP
      pheno.samp$pulse_pressure <- pheno.samp$blood_pressure_systolic - pheno.samp$blood_pressure_diastolic
      pheno.samp$heart_rate <- vital_signs_merged$HRSUP
      pheno.samp$standing_blood_pressure <- vital_signs_merged$SYSSTND
      pheno.samp$standing_heart_rate <- vital_signs_merged$HRSTND
      pheno.samp$height <- vital_signs_merged$HTCM
      pheno.samp$weight <- vital_signs_merged$WGTKG
      pheno.samp$BMI <- pheno.samp$weight / (pheno.samp$height / 100)^2
    }

    # --------------------------------------------------------------------------
    # Age at Visit: High-Precision Version from Source File
    # --------------------------------------------------------------------------
    {
      age_at_visit_data <- read.csv(path_age_at_visit)
      age_at_visit_data$SampleId <- paste(age_at_visit_data$PATNO, age_at_visit_data$EVENT_ID, sep = "_")
      mean(dataOther$SampleId %in% age_at_visit_data$SampleId) # Missing one prodromal sample
      dataOther[which(!dataOther$SampleId %in% age_at_visit_data$SampleId), ]

      age_at_visit_merged <- dplyr::left_join(
        pheno.samp,
        dplyr::select(age_at_visit_data, c("SampleId", "AGE_AT_VISIT")),
        by = "SampleId"
      )
      identical(pheno.samp$SampleId, age_at_visit_merged$SampleId) == T

      plot(age_at_visit_merged$age_at_visit, age_at_visit_merged$AGE_AT_VISIT)
      hist((age_at_visit_merged$age_at_visit - age_at_visit_merged$AGE_AT_VISIT) * 12) # Jitter < 1.5 months
      max(abs(age_at_visit_merged$age_at_visit - age_at_visit_merged$AGE_AT_VISIT), na.rm = T) * 12

      identical(is.na(age_at_visit_merged$age_at_visit), is.na(age_at_visit_merged$AGE_AT_VISIT)) == T

      # Save calculated version, use source file version as primary
      pheno.samp$age_at_visit_calculated <- pheno.samp$age_at_visit
      pheno.samp$age_at_visit <- age_at_visit_merged$AGE_AT_VISIT
    }
  }
}
# ==============================================================================


# ==============================================================================
# ASSEMBLE FINAL DATABASE
# ==============================================================================
# Combine all processed data into a single list object:
#   - dataProt: Protein expression (NPX)
#   - pheno.samp: Sample-level phenotypes
#   - pheno.don: Donor-level phenotypes
#   - nom: Updated nomenclature
#   - nom.original: Original nomenclature
#   - timestamp: Processing timestamp
#   - dataset.name: Dataset identifier
# ==============================================================================
{
  p9k <- list()
  p9k$dataProt <- dataProt
  p9k$pheno.samp <- pheno.samp
  p9k$pheno.don <- pheno.don
  p9k$nom <- nom
  p9k$nom.original <- nom.orig
  p9k$timestamp <- Sys.time()
  p9k$dataset.name <- "PPMI_project_9000_from_2025_September"
}
# ==============================================================================

p9k |>
  readr::write_rds(output_filename, compress = "gz")
