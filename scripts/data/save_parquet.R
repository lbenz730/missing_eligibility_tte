### Script to convert large .sas7bdat files to parquet for easier loading
library(tidyverse)
library(glue)
library(haven)
library(arrow)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
if(!dir.exists(glue('{ehr_dir}/parquet_files/'))) {
  dir.create(glue('{ehr_dir}/parquet_files/'))
}
if(!dir.exists(glue('{ehr_dir}/parquet_files/cases'))) {
  dir.create(glue('{ehr_dir}/parquet_files/cases'))
}

if(!dir.exists(glue('{ehr_dir}/parquet_files/controls'))) {
  dir.create(glue('{ehr_dir}/parquet_files/controls'))
}

### Cases
blood_pressure_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/blood_pressure_cases.sas7bdat'))
write_parquet(blood_pressure_cases, glue('{ehr_dir}/parquet_files/cases/blood_pressure_cases.parquet'))

hba1c_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/hba1c_glucose_labs_cases.sas7bdat'))
write_parquet(hba1c_cases, glue('{ehr_dir}/parquet_files/cases/hba1c_glucose_labs_cases.parquet'))

diabetes_icd9_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_diabetes_dx_cases.sas7bdat'))
write_parquet(diabetes_icd9_cases, glue('{ehr_dir}/parquet_files/cases/raw_diabetes_dx_cases.parquet'))

diabetes_rx_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_diabetes_rx.sas7bdat'))
write_parquet(diabetes_rx_cases, glue('{ehr_dir}/parquet_files/cases/all_diabetes_rx.parquet'))

kidney_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/kidney_lab_cases_final.sas7bdat'))
write_parquet(kidney_cases, glue('{ehr_dir}/parquet_files/cases/kidney_lab_cases_final.parquet'))

neuropathy_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_neuropathy_dx_cases.sas7bdat'))
write_parquet(neuropathy_cases, glue('{ehr_dir}/parquet_files/cases/raw_neuropathy_dx_cases.parquet'))

retinopathy_dx_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/retinopathy_dx.sas7bdat'), encoding = 'latin1')
write_parquet(retinopathy_dx_cases, glue('{ehr_dir}/parquet_files/cases/retinopathy_dx.parquet'))
retinopathy_px_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/retinopathy_px.sas7bdat'), encoding = 'latin1')
write_parquet(retinopathy_px_cases, glue('{ehr_dir}/parquet_files/cases/retinopathy_px.parquet'))

smoking_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/smoking_vdw.sas7bdat'))
write_parquet(smoking_cases, glue('{ehr_dir}/parquet_files/cases/smoking_vdw.parquet'))

### Controls
blood_pressure_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/blood_pressure_controls.sas7bdat'))
write_parquet(blood_pressure_controls, glue('{ehr_dir}/parquet_files/controls/blood_pressure_controls.parquet'))

hba1c_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/hba1c_glucose_labs_controls.sas7bdat'))
write_parquet(hba1c_controls, glue('{ehr_dir}/parquet_files/controls/hba1c_glucose_labs_controls.parquet'))

diabetes_icd9_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_diabetes_dx_controls.sas7bdat'))
write_parquet(diabetes_icd9_controls, glue('{ehr_dir}/parquet_files/controls/raw_diabetes_dx_controls.parquet'))

diabetes_rx_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/rx_diabetes_controls.sas7bdat'))
write_parquet(diabetes_rx_controls, glue('{ehr_dir}/parquet_files/controls/rx_diabetes_controls.parquet'))

kidney_controls <- 
  bind_rows(
    read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/kidney_lab_control_final_pt1.sas7bdat')),
    read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/kidney_lab_control_final_pt2.sas7bdat'))
  ) 
write_parquet(kidney_controls, glue('{ehr_dir}/parquet_files/controls/kidney_lab_control_final.parquet'))

control_pregnancy <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_pregnancy_dx_controls.sas7bdat'))
write_parquet(control_pregnancy, glue('{ehr_dir}/parquet_files/controls/raw_pregnancy_dx_controls.parquet'))

neuropathy_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_neuropathy_dx_controls.sas7bdat'))
write_parquet(neuropathy_controls, glue('{ehr_dir}/parquet_files/controls/raw_neuropathy_dx_controls.parquet'))

retinopathy_dx_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/retinopathy_dx_controls.sas7bdat'), encoding = 'latin1')
write_parquet(retinopathy_dx_controls, glue('{ehr_dir}/parquet_files/controls/retinopathy_dx_controls.parquet'))
retinopathy_px_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/retinopathy_px_controls.sas7bdat'), encoding = 'latin1')
write_parquet(retinopathy_px_controls, glue('{ehr_dir}/parquet_files/controls/retinopathy_px_controls.parquet'))

smoking_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/smoking_vdw_controls.sas7bdat'))
write_parquet(smoking_controls, glue('{ehr_dir}/parquet_files/controls/smoking_vdw_controls.parquet'))