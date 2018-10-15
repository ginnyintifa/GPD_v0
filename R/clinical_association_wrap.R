

#' Build univariate cox regression model with confounders adjusted (age, gender, race)
#' 
#' This function generates one file as the result, the coefficients and significance level estimated for the model.
#' 
#' @param piu_filename The name of the file containing PIU of interest. 
#' @param patient_info_df The dataframe containing patient baseline information.
#' @param patient_outcome_df The dataframe containing survival information.
#' @param mutation_type Somatic or germline mutation, default is "somatic".
#' @param piu_of_interest Domain or PTM sites as the PIU of interest, default is "domain".
#' @param row_sum_min The minimum number of total occurrence for each PIU.
#' @param output_dir The directory you would like to have your output files in.
#' @import dplyr magrittr data.table survival qvalue
#' @export
#' @details 
#' @examples 


univariate_cox_model_for_piu = function(piu_filename,
                                patient_info_df,
                                patient_outcome_df,
                                mutation_type = "somatic",
                                piu_of_interest = "domain",
                                row_sum_min = 1,
                                output_dir)
{
  
  
  #system(paste0("mkdir ", output_dir))
  
  
  
  piu_unite = piu_counts_clinical_unite (piu_count_filename = piu_filename,
                                                 patient_info_df = patient_info_df,
                                                 patient_outcome_df = patient_outcome_df,
                                                 piu_of_interest = piu_of_interest,
                                                 row_sum_min = row_sum_min)
  
  fit_survival_model(interest_variable_info = piu_unite[[3]],
                     unite_data = piu_unite[[1]],
                     output_dir = output_dir,
                     output_name = paste0(mutation_type, "_",piu_of_interest,"_univariate_cox.tsv"))
  
 
  cat("Univariate model on piu fitted!", "\n")
  

}




#' Build univariate cox regression model for gene level counts with confounders adjusted (age, gender, race)
#' 
#' This function generates one file as the result, the coefficients and significance level estimated for the model.
#' 
#' @param gene_level_count_filename The name of the file having gene level count of interest. 
#' @param quiry_gene_id The set of gene ids of interest. Because we try to avoid analysing all the genes.
#' @param patient_info_df The dataframe containing patient baseline information.
#' @param patient_outcome_df The dataframe containing survival information.
#' @param row_sum_min The minimum number of total occurrence for each PIU.
#' @param output_dir The directory you would like to have your output files in.
#' @param output_name The name of the output file.
#' @import dplyr magrittr data.table survival qvalue
#' @export
#' @details 
#' @examples 

univariate_cox_model_for_gene_level_summary_count = function(gene_level_count_filename,
                                                             quiry_gene_id,
                                                             patient_info_df,
                                                             patient_outcome_df,
                                                             row_sum_min = 1,
                                                             output_dir,
                                                             output_name)
  
{
  #system(paste0("mkdir ", output_dir))
  
  
  gene_level_unite = gene_level_counts_clinical_unite (gene_level_count_filename = gene_level_count_filename, 
                                                       quiry_gene_id = quiry_gene_id,
                                                       patient_info_df = patient_info_df,
                                                       patient_outcome_df = patient_outcome_df,
                                                       row_sum_min = row_sum_min)
  
  fit_survival_model (interest_variable_info = gene_level_unite[[2]],
                      unite_data = gene_level_unite[[1]],
                      output_dir = output_dir,
                      output_name = output_name)
  
  
  
  cat("Univariate model on gene level counts fitted!", "\n")
  
}



#' Build univariate cox regression model for gene level counts with confounders adjusted (age, gender, race)
#' 
#' This function generates one file as the result, the coefficients and significance level estimated for the model.
#' 
#' @param somatic_piu_filename The name of the file having somatic PIU counts.
#' @param germline_piu_filename The name of the file having germline PIU counts.
#' @param somatic_bpiu_filename The name of the file having somatic bPIU counts.
#' @param germline_bpiu_filename The name of the file having germline bPIU counts.
#' @param germline_npc_filename The name of the file having germline npc counts.
#' @param patient_info_df The dataframe containing patient baseline information.
#' @param patient_outcome_df The dataframe containing survival information.
#' @param piu_of_interest Domain or PTM sites as the PIU of interest, default is "domain".
#' @param row_sum_min The minimum number of total occurrence for each PIU.
#' @param output_dir The directory you would like to have your output files in.
#' @param output_name The name of the output file.
#' @import dplyr magrittr data.table survival qvalue
#' @export
#' @details 
#' @examples 
#' multivariate_cox_model_with_interaction(somatic_piu_filename = "/data/ginny/tcga_pancan/STAD_somatic/STAD_summarise_mutation/piu_mapping_count.tsv",
#'                                        germline_piu_filename = "/data/ginny/tcga_pancan/germline_raw_process/STAD_summarise_mutation/piu_mapping_count.tsv",
#'                                        somatic_bpiu_filename = "/data/ginny/tcga_pancan/STAD_somatic/STAD_summarise_mutation/bpiu_summarising_count.tsv",
#'                                        germline_bpiu_filename = "/data/ginny/tcga_pancan/germline_raw_process/STAD_summarise_mutation/bpiu_summarising_count.tsv",
#'                                        germline_npc_filename = "/data/ginny/tcga_pancan/germline_raw_process/STAD_summarise_mutation/npc_summarising_count.tsv",
#'                                        patient_info_df = patient_info,
#'                                        patient_outcome_df = patient_outcome,
#'                                        piu_of_interest = "domain",
#'                                        row_sum_min = 1,
#'                                        output_dir = "/data/ginny/tcga_pancan/STAD_somatic/cox_model/",
#'                                        output_name = "mutlivariate_interaction_model.tsv")



multivariate_cox_model_with_interaction = function(somatic_piu_filename,
                                                   germline_piu_filename,
                                                   somatic_bpiu_filename,
                                                   germline_bpiu_filename,
                                                   germline_npc_filename,
                                                   patient_info_df,
                                                   patient_outcome_df,
                                                   piu_of_interest,
                                                   row_sum_min,
                                                   output_dir,
                                                   output_name)
{
  # 
  # somatic_piu_filename = "/data/ginny/tcga_pancan/THCA_somatic/THCA_summarise_mutation/piu_mapping_count.tsv"
  # germline_piu_filename = "/data/ginny/tcga_pancan/germline_raw_process/THCA_summarise_mutation/piu_mapping_count.tsv"
  # somatic_bpiu_filename = "/data/ginny/tcga_pancan/THCA_somatic/THCA_summarise_mutation/bpiu_summarising_count.tsv"
  # germline_bpiu_filename = "/data/ginny/tcga_pancan/germline_raw_process/THCA_summarise_mutation/bpiu_summarising_count.tsv"
  # germline_npc_filename = "/data/ginny/tcga_pancan/germline_raw_process/THCA_summarise_mutation/npc_summarising_count.tsv"
  # patient_info_df = patient_info
  # patient_outcome_df = patient_outcome
  # piu_of_interest = "domain"
  # row_sum_min = 1
  # output_dir = "/data/ginny/tcga_pancan/THCA_somatic/cox_model/"
  # output_name = "mutlivariate_interaction_model.tsv"
  # 
  # 
  
  
  
   somatic_piu_unite = piu_counts_clinical_unite (piu_count_filename = somatic_piu_filename,
                                                 patient_info_df = patient_info_df,
                                                 patient_outcome_df = patient_outcome_df,
                                                 piu_of_interest = piu_of_interest,
                                                 row_sum_min = row_sum_min)
  
germline_npc_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_npc_filename,
                                                             quiry_gene_id = unique(somatic_piu_unite[[4]]),
                                                             patient_info_df = patient_info,
                                                             patient_outcome_df = patient_outcome,
                                                             row_sum_min = row_sum_min)
  
  
  
  summarise_piu_towards_gene(piu_filename = germline_piu_filename,
                             piu_of_interest = "domain",
                             output_dir = output_dir,
                             output_name = "germline_domain_mapping_gene_level_count.tsv")
  
  germline_domain_gene_level_filename = paste0(output_dir, "germline_domain_mapping_gene_level_count.tsv")
  
  summarise_piu_towards_gene(piu_filename = germline_piu_filename,
                             piu_of_interest = "ptm",
                             output_dir = output_dir,
                             output_name = "germline_ptm_mapping_gene_level_count.tsv")
  
  germline_ptm_gene_level_filename = paste0(output_dir, "germline_ptm_mapping_gene_level_count.tsv")
  
  
 germline_domain_gene_level_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_domain_gene_level_filename,
                                                                           quiry_gene_id = somatic_piu_unite[[4]],
                                                                           patient_info_df = patient_info,
                                                                           patient_outcome_df = patient_outcome,
                                                                           row_sum_min = row_sum_min)
 
 
 germline_ptm_gene_level_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_ptm_gene_level_filename,
                                                                     quiry_gene_id = somatic_piu_unite[[4]],
                                                                     patient_info_df = patient_info,
                                                                     patient_outcome_df = patient_outcome,
                                                                     row_sum_min = row_sum_min)
 
 
 germline_bpiu_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_bpiu_filename,
                                                       quiry_gene_id = somatic_piu_unite[[4]],
                                                       patient_info_df = patient_info_df,
                                                       patient_outcome_df = patient_outcome_df,
                                                       row_sum_min = row_sum_min)
 
 somatic_bpiu_unite = gene_level_counts_clinical_unite(gene_level_count_filename = somatic_bpiu_filename,
                                                       quiry_gene_id = somatic_piu_unite[[4]],
                                                       patient_info_df = patient_info_df,
                                                       patient_outcome_df = patient_outcome_df,
                                                       row_sum_min = row_sum_min)
 
 fit_survival_model_with_bpiu_germline(interest_variable_info = somatic_piu_unite[[3]],
                                       unite_data = somatic_piu_unite[[1]],
                                       piu_gene_df = somatic_piu_unite[[2]],
                                       germline_domain_gene_level_count = germline_domain_gene_level_unite[[1]],
                                       germline_ptm_gene_level_count = germline_ptm_gene_level_unite[[1]],
                                       germline_npc_count = germline_npc_unite[[1]],
                                       germline_bpiu_count = germline_bpiu_unite[[1]],
                                       somatic_bpiu_count = somatic_bpiu_unite[[1]],
                                       output_dir = output_dir,
                                       output_name = "somatic_domain_with_bpiu_germline_piu_npc.tsv")
 
  
 
 cat("Multivariate model on with interaction fitted!", "\n")
 
  
}
  
  
  
  


  
  
  
  



multivariate_cox_model_with_interaction_no_gender = function(somatic_piu_filename,
                                                   germline_piu_filename,
                                                   somatic_bpiu_filename,
                                                   germline_bpiu_filename,
                                                   germline_npc_filename,
                                                   patient_info_df,
                                                   patient_outcome_df,
                                                   piu_of_interest,
                                                   row_sum_min,
                                                   output_dir,
                                                   output_name)
{
  
  somatic_piu_unite = piu_counts_clinical_unite (piu_count_filename = somatic_piu_filename,
                                                 patient_info_df = patient_info_df,
                                                 patient_outcome_df = patient_outcome_df,
                                                 piu_of_interest = piu_of_interest,
                                                 row_sum_min = row_sum_min)
  
  germline_npc_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_npc_filename,
                                                        quiry_gene_id = unique(somatic_piu_unite[[4]]),
                                                        patient_info_df = patient_info,
                                                        patient_outcome_df = patient_outcome,
                                                        row_sum_min = row_sum_min)
  
  
  
  summarise_piu_towards_gene(piu_filename = germline_piu_filename,
                             piu_of_interest = "domain",
                             output_dir = output_dir,
                             output_name = "germline_domain_mapping_gene_level_count.tsv")
  
  germline_domain_gene_level_filename = paste0(output_dir, "germline_domain_mapping_gene_level_count.tsv")
  
  summarise_piu_towards_gene(piu_filename = germline_piu_filename,
                             piu_of_interest = "ptm",
                             output_dir = output_dir,
                             output_name = "germline_ptm_mapping_gene_level_count.tsv")
  
  germline_ptm_gene_level_filename = paste0(output_dir, "germline_ptm_mapping_gene_level_count.tsv")
  
  
  germline_domain_gene_level_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_domain_gene_level_filename,
                                                                      quiry_gene_id = somatic_piu_unite[[4]],
                                                                      patient_info_df = patient_info,
                                                                      patient_outcome_df = patient_outcome,
                                                                      row_sum_min = row_sum_min)
  
  
  germline_ptm_gene_level_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_ptm_gene_level_filename,
                                                                   quiry_gene_id = somatic_piu_unite[[4]],
                                                                   patient_info_df = patient_info,
                                                                   patient_outcome_df = patient_outcome,
                                                                   row_sum_min = row_sum_min)
  
  
  germline_bpiu_unite = gene_level_counts_clinical_unite(gene_level_count_filename = germline_bpiu_filename,
                                                         quiry_gene_id = somatic_piu_unite[[4]],
                                                         patient_info_df = patient_info_df,
                                                         patient_outcome_df = patient_outcome_df,
                                                         row_sum_min = row_sum_min)
  
  somatic_bpiu_unite = gene_level_counts_clinical_unite(gene_level_count_filename = somatic_bpiu_filename,
                                                        quiry_gene_id = somatic_piu_unite[[4]],
                                                        patient_info_df = patient_info_df,
                                                        patient_outcome_df = patient_outcome_df,
                                                        row_sum_min = row_sum_min)
  
  fit_survival_model_with_bpiu_germline_no_gender(interest_variable_info = somatic_piu_unite[[3]],
                                        unite_data = somatic_piu_unite[[1]],
                                        piu_gene_df = somatic_piu_unite[[2]],
                                        germline_domain_gene_level_count = germline_domain_gene_level_unite[[1]],
                                        germline_ptm_gene_level_count = germline_ptm_gene_level_unite[[1]],
                                        germline_npc_count = germline_npc_unite[[1]],
                                        germline_bpiu_count = germline_bpiu_unite[[1]],
                                        somatic_bpiu_count = somatic_bpiu_unite[[1]],
                                        output_dir = output_dir,
                                        output_name = "somatic_domain_with_bpiu_germline_piu_npc.tsv")
  
  
  
  cat("Multivariate model on with interaction fitted!", "\n")
  
  
}










  
  
  
