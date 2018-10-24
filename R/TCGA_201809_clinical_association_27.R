### functions for association 


# get the piu counts and CDR clinical info united before modelling ----------------


piu_counts_cdr_clinical_unite = function(piu_count_filename,
                                     cdr_clinical,
                                     piu_of_interest,
                                     row_sum_min,
                                     output_dir,
                                     output_name)
  
{
  
  # piu_count_filename = "/data/ginny/tcga_pancan/stad_somatic/stad_summarise_mutation/piu_mapping_count.tsv"
  # cdr_clinical = stad_cdr
  # piu_of_interest = "domain"
  # row_sum_min = 1
  # output_dir = 
  
  piu_count_df = fread(piu_count_filename,
                       stringsAsFactors = F)
  
  piu_count_sel = piu_count_df%>%
    dplyr::filter(unit_label == piu_of_interest)%>%
    dplyr::mutate(piu_info = paste(uniprot_accession, start_position,end_position, unit_name, gene_name, gene_id, sep = "_"))%>%
    dplyr::mutate(gene_info = paste(gene_name, gene_id,sep = "_"))%>%
    dplyr::filter(row_sum > row_sum_min) %>%
    dplyr::select(uniprot_accession, start_position, end_position, unit_name, gene_name, gene_id, unit_label,gene_info, piu_info, row_sum,
                  everything())
  
  #piu_info = piu_count_sel$piu_info
  
  
  which_barcode = grep("TCGA", colnames(piu_count_sel))
  
  piu_matrix = t(as.matrix(piu_count_sel[,which_barcode]))
  
  piu_count = data.frame(barcode = colnames(piu_count_sel)[which_barcode],
                         piu_matrix,
                         stringsAsFactors = F)
  
  colnames(piu_count) = c("barcode", piu_count_sel$piu_info)
  rownames(piu_count) = NULL
  
  piu_gene_df = data.frame(piu_info = piu_count_sel$piu_info,
                           gene_info = piu_count_sel$gene_info,
                           stringsAsFactors = F)
  
  
  
  
  
  piu_clinical_unite_data = piu_count %>%
    dplyr::left_join(cdr_clinical, by = c("barcode" = "bcr_patient_barcode")) %>%
    dplyr::select(barcode,colnames(cdr_clinical)[3:34],everything())
  
   
  # piu_clinical_unite_data = piu_count %>%
  #   dplyr::left_join(patient_outcome, by = c("barcode" = "bcr_patient_barcode")) %>%
  #   na.omit() %>%
  #   dplyr::left_join(patient_info, by  = c("barcode" = "bcr_patient_barcode"))%>%
  #   na.omit() %>% 
  #   dplyr::arrange(desc(surv_time))%>%
  #   dplyr::arrange(vital_status) 
  
  
  return_list = vector(mode = "list", length = 4)
  
  return_list[[1]] = piu_clinical_unite_data
  return_list[[2]] = piu_gene_df
  return_list[[3]] = piu_count_sel$piu_info
  return_list[[4]] = piu_count_sel$gene_id
  
  
  
  write.table(return_list[[1]], paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
  
  return(return_list)
  
}





# get the locus level counts and CDR clinical info united before modelling ----------------


locus_counts_cdr_clinical_unite = function(locus_count_filename,
                                         cdr_clinical,
                                         row_sum_min,
                                         output_dir,
                                         output_name)
  
{
  
  # locus_count_filename = "/data/ginny/tcga_pancan/STAD_somatic/STAD_summarise_mutation/stad_mc3_count_matrix.tsv"
  # cdr_clinical = stad_cdr
  # row_sum_min = 1
  # output_dir =  "/data/ginny/tcga_pancan/STAD_somatic/cox_model/"
  # output_name = "stad_locus_level_cdr_clinical_unite.tsv"
  # 
  locus_count_df = fread(locus_count_filename,
                       stringsAsFactors = F)
  
  locus_count_sel = locus_count_df%>%
    dplyr::filter(row_sum > row_sum_min)
    
  #piu_info = piu_count_sel$piu_info
  
  
  which_barcode = grep("TCGA", colnames(locus_count_sel))
  
  locus_matrix = t(as.matrix(locus_count_sel[,which_barcode]))
  
  locus_count = data.frame(barcode = colnames(locus_count_sel)[which_barcode],
                         locus_matrix,
                         stringsAsFactors = F)
  
  colnames(locus_count) = c("barcode", locus_count_sel$locus_info)
  rownames(locus_count) = NULL
  
  
  locus_clinical_unite_data = locus_count %>%
    dplyr::left_join(cdr_clinical, by = c("barcode" = "bcr_patient_barcode")) %>%
    dplyr::select(barcode,colnames(cdr_clinical)[3:34],everything())
  
  
  
  return_list = vector(mode = "list", length = 2)
  
  return_list[[1]] = locus_clinical_unite_data
  return_list[[2]] = locus_count_sel$locus_info
 
  
  
  write.table(return_list[[1]], paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
  
  return(return_list)
  
}




# get the germline npc counts and clinical info united before modelling ----------------
## used in generating results for point3
gene_level_counts_cdr_clinical_unite = function(gene_level_count_filename,
                                            quiry_gene_id,
                                            cdr_clinical,
                                            row_sum_min)
  
{
  
  # gene_level_count_filename = "/data/ginny/tcga_pancan/raw_process/stad_summarise_mutation/gene_level_summarising_count.tsv"
  # quiry_gene_id = unique(stad_somatic_domain_unite[[4]])
  # cdr_clinical = stad_cdr,
  # row_sum_min = 1
  
  
  gene_level_count_df = fread(gene_level_count_filename,
                              stringsAsFactors = F)
  
  
  gene_level_count_sel = gene_level_count_df %>%
    dplyr::filter(gene_id %in% quiry_gene_id)%>%
    dplyr::mutate(gene_info = paste(gene_name,gene_id,sep = "_"))%>%
    dplyr::filter(row_sum >= row_sum_min) %>%
    dplyr::select(gene_name, gene_id, gene_info, row_sum,
                  everything())
  
  #stad_mind_ncu_info = gene_level_count_sel$gene_info
  which_barcode = grep("TCGA", colnames(gene_level_count_sel))
  
  gene_level_matrix = t(as.matrix(gene_level_count_sel[,which_barcode]))
  
  gene_level_count = data.frame(barcode = colnames(gene_level_count_sel)[which_barcode], gene_level_matrix,
                                stringsAsFactors = F)
  
  colnames(gene_level_count) = c("barcode", gene_level_count_sel$gene_info)
  rownames(gene_level_count) = NULL
  
  
  gene_level_clinical_unite_data = gene_level_count %>%
    dplyr::left_join(cdr_clinical, by = c("barcode" = "bcr_patient_barcode")) 
  
  
  return_list = vector(mode = "list", length = 2)
  
  
  
  return_list[[1]] = gene_level_clinical_unite_data
  return_list[[2]] = gene_level_count_sel$gene_info
  
  
  
  return(return_list)
  
}










# add germline pius of the same gene together -----------------------------

summarise_piu_towards_gene = function(piu_filename,
                                               piu_of_interest,
                                               output_dir,
                                               output_name)
{
  
  #piu_filename = "/data/ginny/tcga_pancan/germline_raw_process/stad_summarise_mutation/piu_mapping_count.tsv"
  
  piu_mapping = fread(piu_filename,stringsAsFactors = F)
  sample_col = grep("TCGA", colnames(piu_mapping), value = T)
  sel_col = c("gene_name", "gene_id", sample_col,"row_sum")
  
  gene_level_piu = piu_mapping %>%
    dplyr::filter(unit_label == piu_of_interest)%>%
    dplyr::select(one_of(sel_col))%>%
    dplyr::group_by(gene_name, gene_id)%>%
    dplyr::summarise_all(funs(sum)) %>%
    dplyr::arrange(desc(row_sum))
  
  write.table(gene_level_piu, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}



# univariate survival analysis  -------------------------------------------
# I will produce result for all four types of clinical end points OS, DSS, DFI, PFI
# so four files I will generate
# think about subtype/tumor stage later

# 
# stad_somatic_domain_unite = piu_counts_cdr_clinical_unite(piu_count_filename = "/data/ginny/tcga_pancan/STAD_somatic/STAD_summarise_mutation/piu_mapping_count.tsv",
#   cdr_clinical = stad_cdr,
#   piu_of_interest = "domain",
#   row_sum_min = 1,
#   output_dir = "/data/ginny/tcga_pancan/STAD_somatic/cox_model/",
#   output_name = "stad_somatic_domain_cdr_clinical_unite.tsv")
# 




cdr_tidy_up_for_model = function(interest_variable_info, unite_data, race_group_min, 
                                 output_dir,
                                 output_name)
{
  
  # interest_variable_info = locus_unite[[2]]
  # unite_data = locus_unite[[1]]
  # race_group_min = 5
  # output_dir = output_dir
  # output_name = paste0(mutation_type, "_locus_level_survival_info.tsv")
  
  patient_count = unite_data %>%
    dplyr::select(barcode, age_at_initial_pathologic_diagnosis,gender,one_of(interest_variable_info)) %>%
    dplyr::rename(age = age_at_initial_pathologic_diagnosis)
  
  unite_data$race[grep("\\[", unite_data$race)] = "OTHER"
  
  os_data = unite_data %>%
    dplyr::select(barcode, race, OS, OS.time) %>%
    na.omit()
  os_race_freq = as.data.frame(table(os_data$race))
  os_race_minor = os_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  os_data$race[which(os_data$race %in% os_race_minor$Var1)] = "OTHER"
  colnames(os_data)[which(colnames(os_data)=="race")] = "os_race"
  
  
  dss_data = unite_data %>%
    dplyr::select(barcode, race, DSS, DSS.time) %>%
    na.omit()
  dss_race_freq = as.data.frame(table(dss_data$race))
  dss_race_minor = dss_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  dss_data$race[which(dss_data$race %in% dss_race_minor$Var1)] = "OTHER"
  colnames(dss_data)[which(colnames(dss_data)=="race")] = "dss_race"
  
  
  dfi_data = unite_data %>%
    dplyr::select(barcode, race, DFI, DFI.time) %>%
    na.omit()
  
  dfi_race_freq = as.data.frame(table(dfi_data$race))
  dfi_race_minor = dfi_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  dfi_data$race[which(dfi_data$race %in% dfi_race_minor$Var1)] = "OTHER"
  colnames(dfi_data)[which(colnames(dfi_data)=="race")] = "dfi_race"
  
  
  pfi_data = unite_data %>%
    dplyr::select(barcode, race, PFI, PFI.time) %>%
    na.omit()
  pfi_race_freq = as.data.frame(table(pfi_data$race))
  pfi_race_minor = pfi_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  pfi_data$race[which(pfi_data$race %in% pfi_race_minor$Var1)] = "OTHER"
  colnames(pfi_data)[which(colnames(pfi_data)=="race")] = "pfi_race"
  
  
  
  
  patient_count_survival = patient_count %>%
    dplyr::left_join(os_data, by = "barcode") %>%
    dplyr::left_join(dss_data, by = "barcode") %>%
    dplyr::left_join(dfi_data, by =  "barcode") %>%
    dplyr::left_join(pfi_data, by = "barcode") %>%
    dplyr::select(barcode, age, gender,
                  os_race, OS, OS.time,
                  dss_race, DSS, DSS.time,
                  dfi_race, DFI, DFI.time,
                  pfi_race, PFI, PFI.time,
                  one_of(interest_variable_info))
  
  write.table(patient_count_survival, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  return(patient_count_survival)

}
  


### the four end points may not be all meaningful 
### some does not have enough data 
### so in this case I need a choosing procedure 


fit_survival_model = function(surv_info_data,
                              interest_variable_info,
                              min_surv_time,
                              min_surv_people,
                              output_dir,
                              output_name)
  
{
  
  
  # surv_info_data = piu_info
  # interest_variable_info = piu_unite[[3]]
  # min_surv_time = 30
  # min_surv_people = 5
  # output_dir = output_dir
  # output_name = paste0(mutation_type,"_", piu_of_interest, "_cdr_univariate_piu.tsv")
   
  endpoint_flag = data.frame(OS = T, DSS = T, DFI = T, PFI = T, stringsAsFactors = F)
  
  os_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender, os_race, OS, OS.time) %>%
    na.omit() %>%
    dplyr::arrange(OS, OS.time) %>%
    dplyr::filter(OS.time >= min_surv_time)
  
  os_surv_table = as.data.frame(table(os_surv_data$OS))
  if(min(os_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$OS = F
  }else{
    os_age = as.numeric(os_surv_data$age)
    os_gender = relevel(as.factor(os_surv_data$gender), ref = unique(os_surv_data$gender)[1])
    os_race = relevel(as.factor(os_surv_data$os_race), ref = unique(os_surv_data$os_race)[1])
    os_surv_object = Surv(time = os_surv_data$OS.time, event =  os_surv_data$OS)
  }
  ###
  dss_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender, dss_race, DSS, DSS.time) %>%
    na.omit() %>%
    dplyr::arrange(DSS, DSS.time)%>%
    dplyr::filter(DSS.time >= min_surv_time)
  dss_surv_table = as.data.frame(table(dss_surv_data$DSS))
  
  if(min(dss_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$DSS = F
  }else{
  dss_age = as.numeric(dss_surv_data$age)
  dss_gender = relevel(as.factor(dss_surv_data$gender), ref = unique(dss_surv_data$gender)[1])
  dss_race = relevel(as.factor(dss_surv_data$dss_race), ref = unique(dss_surv_data$dss_race)[1])
  dss_surv_object = Surv(time = dss_surv_data$DSS.time, event =  dss_surv_data$DSS)
  }
  ###
  dfi_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender, dfi_race, DFI, DFI.time) %>%
    na.omit() %>%
    dplyr::arrange(DFI, DFI.time)%>%
    dplyr::filter(DFI.time >= min_surv_time)
  dfi_surv_table = as.data.frame(table(dfi_surv_data$DFI))
  if(min(dfi_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$DFI = F
  }else{
  dfi_age = as.numeric(dfi_surv_data$age)
  dfi_gender = relevel(as.factor(dfi_surv_data$gender), ref = unique(dfi_surv_data$gender)[1] )
  dfi_race = relevel(as.factor(dfi_surv_data$dfi_race), ref = unique(dfi_surv_data$dfi_race)[1])
  dfi_surv_object = Surv(time = dfi_surv_data$DFI.time, event =  dfi_surv_data$DFI)
  }
  ###
  
  pfi_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender, pfi_race, PFI, PFI.time) %>%
    na.omit() %>%
    dplyr::arrange(PFI, PFI.time)%>%
    dplyr::filter(PFI.time >= min_surv_time)
  pfi_surv_table = as.data.frame(table(pfi_surv_data$PFI))
  if(min(pfi_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$PFI = F
  }else{
  pfi_age = as.numeric(pfi_surv_data$age)
  pfi_gender = relevel(as.factor(pfi_surv_data$gender), ref = unique(pfi_surv_data$gender)[1])
  pfi_race = relevel(as.factor(pfi_surv_data$pfi_race), ref = unique(pfi_surv_data$pfi_race)[1])
  pfi_surv_object = Surv(time = pfi_surv_data$PFI.time, event =  pfi_surv_data$PFI)
  }
  
  
  #####
  
 surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(i)
      # for(i in 600:length(interest_variable_info))
  {
  #  i = 552

    this_count_df = surv_info_data %>%
      dplyr::select(barcode, one_of(interest_variable_info[i]))
    
    ###
    os_count = os_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_os = sum(os_count[,7]!=0)
    total_patients_os = nrow(os_count)
    count_coeff_os = NA
    count_exp_coeff_os = NA
    count_pval_os = NA
    
    if(endpoint_flag$OS == T)
    {
      
      
      if(length(unique(os_gender))>1)
      {
        if(length(unique(os_race))>1)
        {
          os_model = coxph(os_surv_object ~  os_age  + os_gender + os_race+ 
                             os_count[,7])   
        }else{
          os_model = coxph(os_surv_object ~  os_age  + os_gender +
                             os_count[,7])   
        }

      }else{
        if(length(unique(os_race))>1)
        {
          os_model = coxph(os_surv_object ~  os_age  + os_race+ 
                             os_count[,7])   
        }else{
          os_model = coxph(os_surv_object ~  os_age +
                             os_count[,7])   
        }
        
      }
      
     # r_os = 1 + length(unique(os_gender))-1 +length(unique(os_race))-1 +1
      os = summary(os_model)
      os_coef = os$coefficients
      r_os =  length(unique(os_gender)) +length(unique(os_race))
      count_coeff_os = os_coef[r_os,1]
      count_exp_coeff_os = os_coef[r_os,2]
      count_pval_os = os_coef[r_os,5]
      
    }
    
  
    ###
    
    dss_count = dss_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_dss = sum(dss_count[,7]!=0)
    total_patients_dss = nrow(dss_count)
    count_coeff_dss = NA
    count_exp_coeff_dss = NA
    count_pval_dss = NA
    
    if(endpoint_flag$DSS == T)
    {
      if(length(unique(dss_gender))>1)
      {
        if(length(unique(dss_race))>1)
        {
          
          dss_model = coxph(dss_surv_object ~  dss_age  + dss_gender + dss_race+ 
                              dss_count[,7])   
          
        }else{
          
          dss_model = coxph(dss_surv_object ~  dss_age  + dss_gender +
                              dss_count[,7])   
          
        }
      }else{
        if(length(unique(dss_race))>1)
        {
          
          dss_model = coxph(dss_surv_object ~  dss_age  + dss_race+ 
                              dss_count[,7])   
          
        }else{
          
          dss_model = coxph(dss_surv_object ~  dss_age  +
                              dss_count[,7])   
          
        }
      }
      
      dss = summary(dss_model)
      dss_coef = dss$coefficients
      r_dss = length(unique(dss_gender))+length(unique(dss_race))
      count_coeff_dss = dss_coef[r_dss,1]
      count_exp_coeff_dss = dss_coef[r_dss,2]
      count_pval_dss = dss_coef[r_dss,5]
    }
    
    ###
    
    dfi_count = dfi_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_dfi = sum(dfi_count[,7]!=0)
    total_patients_dfi = nrow(dfi_count)
    count_coeff_dfi = NA
    count_exp_coeff_dfi = NA
    count_pval_dfi = NA
    
    if(endpoint_flag$DFI == T)
    {
      if(length(unique(dfi_gender))>1)
      {
        if(length(unique(dfi_race))>1)
        {
          
          dfi_model = coxph(dfi_surv_object ~  dfi_age  + dfi_gender + dfi_race+ 
                              dfi_count[,7])   
          
        }else{
          
          dfi_model = coxph(dfi_surv_object ~  dfi_age  + dfi_gender +
                              dfi_count[,7])   
          
        }
        
      }else{
        if(length(unique(dfi_race))>1)
        {
          
          dfi_model = coxph(dfi_surv_object ~  dfi_age  +  dfi_race+ 
                              dfi_count[,7])   
          
        }else{
          
          dfi_model = coxph(dfi_surv_object ~  dfi_age  + 
                              dfi_count[,7])   
          
        }
        
      }
     
        
      dfi = summary(dfi_model)
      dfi_coef = dfi$coefficients
      r_dfi = length(unique(dfi_gender))+length(unique(dfi_race))
      
      count_coeff_dfi = dfi_coef[r_dfi,1]
      count_exp_coeff_dfi = dfi_coef[r_dfi,2]
      count_pval_dfi = dfi_coef[r_dfi,5]
    }
      
    ###
    
    
    pfi_count = pfi_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_pfi = sum(pfi_count[,7]!=0)
    total_patients_pfi = nrow(pfi_count)
    count_coeff_pfi = NA
    count_exp_coeff_pfi = NA
    count_pval_pfi = NA
    
    if(endpoint_flag$PFI == T)
    {
      
      if(length(unique(pfi_gender))>1)
      {
        if(length(unique(pfi_race))>1)
        {
          
          pfi_model = coxph(pfi_surv_object ~  pfi_age  + pfi_gender + pfi_race+ 
                              pfi_count[,7])   
          
        }else{
          
          pfi_model = coxph(pfi_surv_object ~  pfi_age  + pfi_gender +
                              pfi_count[,7])   
          
        }
      }else{
        if(length(unique(pfi_race))>1)
        {
          pfi_model = coxph(pfi_surv_object ~  pfi_age  + pfi_race+ 
                              pfi_count[,7])   
        }else{
          pfi_model = coxph(pfi_surv_object ~  pfi_age  +
                              pfi_count[,7])   
        }
      }
         
      pfi = summary(pfi_model)
      pfi_coef = pfi$coefficients
      r_pfi = length(unique(pfi_gender))+length(unique(pfi_race))
      
      count_coeff_pfi = pfi_coef[r_pfi,1]
      count_exp_coeff_pfi = pfi_coef[r_pfi,2]
      count_pval_pfi = pfi_coef[r_pfi,5]
  
    }
    
    
    
    this_surv_result = data.frame(count_info = interest_variable_info[i],
                                  num_patients_os,total_patients_os,
                                  count_coeff_os, count_exp_coeff_os,count_pval_os,
                                  num_patients_dss, total_patients_dss,
                                  count_coeff_dss , count_exp_coeff_dss,count_pval_dss ,
                                  num_patients_dfi,total_patients_dfi,
                                  count_coeff_dfi, count_exp_coeff_dfi,count_pval_dfi,
                                  num_patients_pfi,total_patients_pfi,
                                  count_coeff_pfi, count_exp_coeff_pfi,count_pval_pfi,
                                  count_qval_os = NA,count_qval_dss = NA,count_qval_dfi = NA,count_qval_pfi = NA,
                                  stringsAsFactors = F)
    
    if(i%%100 == 0 )
      cat(i, "\n")
    
    return(this_surv_result)
    
  }))
  
  
  # q_val = p.adjust(surv_result_df$count_pval,method = "BH")
  ## prevent truncated p value distribution 
  if(length(interest_variable_info)>5)
  {
    if(endpoint_flag$OS == T)
    {
      if(length(surv_result_df$count_pval_os)>300 & min(surv_result_df$count_pval_os,na.rm = T)<0.05 & max(surv_result_df$count_pval_os,na.rm = T)>0.95)
      {
        os_qval = qvalue(surv_result_df$count_pval_os)
      }else{
        os_qval = qvalue(surv_result_df$count_pval_os, pi0 = 1)
        
      }
      surv_result_df$count_qval_os = os_qval$qvalues
      
    }
    
    if(endpoint_flag$DSS == T)
    {  
      if(length(surv_result_df$count_pval_dss)>300 &min(surv_result_df$count_pval_dss,na.rm = T)<0.05 & max(surv_result_df$count_pval_dss,na.rm = T)>0.95)
      {
        dss_qval = qvalue(surv_result_df$count_pval_dss)
      }else{
        dss_qval = qvalue(surv_result_df$count_pval_dss, pi0 = 1)
      }
      surv_result_df$count_qval_dss = dss_qval$qvalues
    }
    
    
    if(endpoint_flag$DFI == T)
    {
      if(length(surv_result_df$count_pval_dfi)>300&min(surv_result_df$count_pval_dfi,na.rm = T)<0.05 & max(surv_result_df$count_pval_dfi,na.rm = T)>0.95)
      {
        dfi_qval = qvalue(surv_result_df$count_pval_dfi)
      }else{
        dfi_qval = qvalue(surv_result_df$count_pval_dfi, pi0 = 1)
      }
      surv_result_df$count_qval_dfi = dfi_qval$qvalues
    }
    if(endpoint_flag$PFI == T)
    {
      if(length(surv_result_df$count_pval_pfi)>300&min(surv_result_df$count_pval_pfi,na.rm = T)<0.05 & max(surv_result_df$count_pval_pfi,na.rm = T)>0.95)
      {
        pfi_qval = qvalue(surv_result_df$count_pval_pfi)
      }else{
        pfi_qval = qvalue(surv_result_df$count_pval_pfi, pi0 = 1)
        
      }
      surv_result_df$count_qval_pfi = pfi_qval$qvalues
    }
    
    
  }else{
    
    cat("Small unit size, no q-val estimation.","\n")
  }
  
  
  
  surv_df = surv_result_df %>%
    dplyr::arrange(desc(num_patients_os)) %>%
    replace(is.na(.),"")
  
  
  write.table(surv_df, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}







fit_survival_model_no_gender = function(surv_info_data,
                              interest_variable_info,
                              min_surv_time,
                              min_surv_people,
                              output_dir,
                              output_name)
  
{
  
  # surv_info_data = locus_info
  # interest_variable_info = locus_unite[[2]]
  # min_surv_time = min_surv_days
  # min_surv_people = 5
  # output_dir =  output_dir
  # output_name = paste0(mutation_type, "_cdr_univariate_locus.tsv")
  # 
  
  ### I should set a variable as a flag to indicate the inclusion of the 4 endpoints
  
  endpoint_flag = data.frame(OS = T, DSS = T, DFI = T, PFI = T, stringsAsFactors = F)
  
  os_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender, os_race, OS, OS.time) %>%
    na.omit() %>%
    dplyr::arrange(OS, OS.time) %>%
    dplyr::filter(OS.time >= min_surv_time)
  
  os_surv_table = as.data.frame(table(os_surv_data$OS))
  if(min(os_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$OS = F
  }else{
    os_age = as.numeric(os_surv_data$age)
    os_gender = relevel(as.factor(os_surv_data$gender), ref = unique(os_surv_data$gender)[1])
    os_race = relevel(as.factor(os_surv_data$os_race), ref = unique(os_surv_data$os_race)[1])
    os_surv_object = Surv(time = os_surv_data$OS.time, event =  os_surv_data$OS)
  }
  ###
  dss_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender,dss_race, DSS, DSS.time) %>%
    na.omit() %>%
    dplyr::arrange(DSS, DSS.time)%>%
    dplyr::filter(DSS.time >= min_surv_time)
  dss_surv_table = as.data.frame(table(dss_surv_data$DSS))
  
  if(min(dss_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$DSS = F
  }else{
    dss_age = as.numeric(dss_surv_data$age)
    dss_gender = relevel(as.factor(dss_surv_data$gender), ref = unique(dss_surv_data$gender)[1])
    dss_race = relevel(as.factor(dss_surv_data$dss_race), ref = unique(dss_surv_data$dss_race)[1])
    dss_surv_object = Surv(time = dss_surv_data$DSS.time, event =  dss_surv_data$DSS)
  }
  ###
  dfi_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender,dfi_race, DFI, DFI.time) %>%
    na.omit() %>%
    dplyr::arrange(DFI, DFI.time)%>%
    dplyr::filter(DFI.time >= min_surv_time)
  dfi_surv_table = as.data.frame(table(dfi_surv_data$DFI))
  if(min(dfi_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$DFI = F
  }else{
    dfi_age = as.numeric(dfi_surv_data$age)
    dfi_gender = relevel(as.factor(dfi_surv_data$gender), ref = unique(dfi_surv_data$gender)[1] )
    dfi_race = relevel(as.factor(dfi_surv_data$dfi_race), ref = unique(dfi_surv_data$dfi_race)[1])
    dfi_surv_object = Surv(time = dfi_surv_data$DFI.time, event =  dfi_surv_data$DFI)
  }
  ###
  
  pfi_surv_data = surv_info_data %>%
    dplyr::select(barcode, age, gender,pfi_race, PFI, PFI.time) %>%
    na.omit() %>%
    dplyr::arrange(PFI, PFI.time)%>%
    dplyr::filter(PFI.time >= min_surv_time)
  pfi_surv_table = as.data.frame(table(pfi_surv_data$PFI))
  if(min(pfi_surv_table$Freq)<min_surv_people)
  {
    endpoint_flag$PFI = F
  }else{
    pfi_age = as.numeric(pfi_surv_data$age)
    pfi_gender = relevel(as.factor(pfi_surv_data$gender), ref = unique(pfi_surv_data$gender)[1])
    pfi_race = relevel(as.factor(pfi_surv_data$pfi_race), ref = unique(pfi_surv_data$pfi_race)[1])
    pfi_surv_object = Surv(time = pfi_surv_data$PFI.time, event =  pfi_surv_data$PFI)
  }
  
  
  #####
  
  surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(i)
  {
    #i = 1
    
    this_count_df = surv_info_data %>%
      dplyr::select(barcode, one_of(interest_variable_info[i]))
    
    ###
    os_count = os_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_os = sum(os_count[,7]!=0)
    total_patients_os = nrow(os_count)
    count_coeff_os = NA
    count_exp_coeff_os = NA
    count_pval_os = NA
    
    if(endpoint_flag$OS == T)
    {
      ### a potential bug, for both gender/no gender functions, race may not be there >=2 types.
      ### I should prepare for this 
      
      if(length(unique(os_race))>1)
      {
        
        os_model = coxph(os_surv_object ~  os_age  + os_race+ 
                           os_count[,7])   
        
      }else{
        
        os_model = coxph(os_surv_object ~  os_age  + 
                           os_count[,7])   
        
      }
      
      os = summary(os_model)
      os_coef = os$coefficients
      r_os = 1+length(unique(os_race))
      
      count_coeff_os = os_coef[r_os,1]
      count_exp_coeff_os = os_coef[r_os,2]
      count_pval_os = os_coef[r_os,5]
      
    }
    
    
    ###
    
    dss_count = dss_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_dss = sum(dss_count[,7]!=0)
    total_patients_dss = nrow(dss_count)
    count_coeff_dss = NA
    count_exp_coeff_dss = NA
    count_pval_dss = NA
    
    if(endpoint_flag$DSS == T)
    {
      
      if(length(unique(dss_race))>1)
      {
        
        dss_model = coxph(dss_surv_object ~  dss_age  + dss_race+ 
                           dss_count[,7])   
        
      }else{
        
        dss_model = coxph(dss_surv_object ~  dss_age  + 
                           dss_count[,7])   
        
      }
      
      dss = summary(dss_model)
      dss_coef = dss$coefficients
      r_dss = 1+length(unique(dss_race))
      count_coeff_dss = dss_coef[r_dss,1]
      count_exp_coeff_dss = dss_coef[r_dss,2]
      count_pval_dss = dss_coef[r_dss,5]
    }
    
    ###
    
    dfi_count = dfi_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_dfi = sum(dfi_count[,7]!=0)
    total_patients_dfi = nrow(dfi_count)
    count_coeff_dfi = NA
    count_exp_coeff_dfi = NA
    count_pval_dfi = NA
    
    if(endpoint_flag$DFI == T)
    {
      
      
      if(length(unique(dfi_race))>1)
      {
        
        dfi_model = coxph(dfi_surv_object ~  dfi_age  + dfi_race+ 
                            dfi_count[,7])   
        
      }else{
        
        dfi_model = coxph(dfi_surv_object ~  dfi_age  + 
                            dss_count[,7])   
        
      }
      
      dfi = summary(dfi_model)
      dfi_coef = dfi$coefficients
      r_dfi = 1+length(unique(dfi_race))
      
      count_coeff_dfi = dfi_coef[r_dfi,1]
      count_exp_coeff_dfi = dfi_coef[r_dfi,2]
      count_pval_dfi = dfi_coef[r_dfi,5]
    }
    
    ###
    
    
    pfi_count = pfi_surv_data %>%
      dplyr::left_join(this_count_df, by = "barcode")
    
    num_patients_pfi = sum(pfi_count[,7]!=0)
    total_patients_pfi = nrow(pfi_count)
    count_coeff_pfi = NA
    count_exp_coeff_pfi = NA
    count_pval_pfi = NA
    
    if(endpoint_flag$PFI == T)
    {
      if(length(unique(pfi_race))>1)
      {
        
        pfi_model = coxph(pfi_surv_object ~  pfi_age  + pfi_race+ 
                            pfi_count[,7])   
        
      }else{
        
        pfi_model = coxph(pfi_surv_object ~  pfi_age  + 
                            pfi_count[,7])   
        
      }
      
      pfi = summary(pfi_model)
      pfi_coef = pfi$coefficients
      r_pfi = 1+length(unique(pfi_race))
      
      count_coeff_pfi = pfi_coef[r_pfi,1]
      count_exp_coeff_pfi = pfi_coef[r_pfi,2]
      count_pval_pfi = pfi_coef[r_pfi,5]
      
    }
    
    
    
    this_surv_result = data.frame(count_info = interest_variable_info[i],
                                  num_patients_os,total_patients_os,
                                  count_coeff_os, count_exp_coeff_os,count_pval_os,
                                  num_patients_dss, total_patients_dss,
                                  count_coeff_dss , count_exp_coeff_dss,count_pval_dss ,
                                  num_patients_dfi,total_patients_dfi,
                                  count_coeff_dfi, count_exp_coeff_dfi,count_pval_dfi,
                                  num_patients_pfi,total_patients_pfi,
                                  count_coeff_pfi, count_exp_coeff_pfi,count_pval_pfi,
                                  count_qval_os = NA,count_qval_dss = NA,count_qval_dfi = NA,count_qval_pfi = NA,
                                  stringsAsFactors = F)
    
    if(i%%100 == 0 )
      cat(i, "\n")
    
    return(this_surv_result)
    
  }))
  
  
  # q_val = p.adjust(surv_result_df$count_pval,method = "BH")
  ## prevent truncated p value distribution 
  
  if(length(interest_variable_info)>5)
  {
    if(endpoint_flag$OS == T)
    {
      if(length(surv_result_df$count_pval_os)>300 & min(surv_result_df$count_pval_os,na.rm = T)<0.05 & max(surv_result_df$count_pval_os,na.rm = T)>0.95)
      {
        os_qval = qvalue(surv_result_df$count_pval_os)
      }else{
        os_qval = qvalue(surv_result_df$count_pval_os, pi0 = 1)
        
      }
      surv_result_df$count_qval_os = os_qval$qvalues
      
    }
    if(endpoint_flag$DSS == T)
    {  
      if(length(surv_result_df$count_pval_dss)>300&min(surv_result_df$count_pval_dss,na.rm = T)<0.05 & max(surv_result_df$count_pval_dss,na.rm = T)>0.95)
      {
        dss_qval = qvalue(surv_result_df$count_pval_dss)
      }else{
        dss_qval = qvalue(surv_result_df$count_pval_dss, pi0 = 1)
      }
      surv_result_df$count_qval_dss = dss_qval$qvalues
    }
    if(endpoint_flag$DFI == T)
    {
      if(length(surv_result_df$count_pval_dfi)>300&min(surv_result_df$count_pval_dfi,na.rm = T)<0.05 & max(surv_result_df$count_pval_dfi,na.rm = T)>0.95)
      {
        dfi_qval = qvalue(surv_result_df$count_pval_dfi)
      }else{
        dfi_qval = qvalue(surv_result_df$count_pval_dfi, pi0 = 1)
      }
      surv_result_df$count_qval_dfi = dfi_qval$qvalues
    }
    if(endpoint_flag$PFI == T)
    {
      if(length(surv_result_df$count_pval_pfi)>300&min(surv_result_df$count_pval_pfi,na.rm = T)<0.05 & max(surv_result_df$count_pval_pfi,na.rm = T)>0.95)
      {
        pfi_qval = qvalue(surv_result_df$count_pval_pfi)
      }else{
        pfi_qval = qvalue(surv_result_df$count_pval_pfi, pi0 = 1)
        
      }
      surv_result_df$count_qval_pfi = pfi_qval$qvalues
    }
    
    
  }else{
    cat("Small unit size, no q-val estimation.", "\n")
  }
  
  
  surv_df = surv_result_df %>%
    dplyr::arrange(desc(num_patients_os)) %>%
    replace(is.na(.),"")
  
  
  write.table(surv_df, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}










# multivariate model  ------------------------------------------------------


fit_survival_model_with_germline = function(interest_variable_info,
                                            unite_data,
                                            piu_gene_df,
                                            germline_domain_gene_level_count,
                                            germline_ptm_gene_level_count,
                                            germline_npc_count,
                                            output_dir,
                                            output_name)
{
  
  # interest_variable_info = stad_somatic_domain_unite[[3]]
  # unite_data = stad_somatic_domain_unite[[1]]
  # piu_gene_df = stad_somatic_domain_unite[[2]]
  # germline_domain_gene_level_count = stad_germline_domain_gene_level_unite[[1]]
  # germline_ptm_gene_level_count = stad_germline_ptm_gene_level_unite[[1]]
  # germline_npc_count = stad_germline_npc_unite[[1]]
  # output_dir = "/data/ginny/tcga_pancan/stad_somatic/stad_clinical_association_27/"
  # output_name = "with_interatcion_stad_somatic_uni_domain_total_survival.tsv"
  # 
  
  survival_data = unite_data %>%
    dplyr::filter(vital_status %in% c("Alive","Dead")) %>%
    dplyr::mutate(surv_status = dplyr::recode(vital_status,
                                              "Alive" = 0, 
                                              "Dead" = 1 ))
  
  surv_patient_df = data.frame(barcode = survival_data$barcode, stringsAsFactors = F)
  
  
  get_age = as.numeric(survival_data$age)
  get_gender = relevel(as.factor(survival_data$gender), ref = "MALE")
  get_race = relevel(as.factor(survival_data$crace), ref = "WHITE")
  
  surv_object = Surv(time = survival_data$surv_time, event =  survival_data$surv_status)
  
  
  survival_count = survival_data %>%
    dplyr::select(one_of(interest_variable_info))
  
  germline_ptm_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_ptm_gene_level_count, by = "barcode") %>%
    replace(is.na(.), 0)

  
  germline_domain_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_domain_gene_level_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  
  
  germline_npc_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_npc_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  
  int_total_surv_result_df = rbindlist(lapply(1:ncol(survival_count), function(i)
    
  {
   # i = 2
    
    this_piu_info = colnames(survival_count)[i]
    this_gene_info = piu_gene_df %>%
      dplyr::filter(piu_info == this_piu_info) %>%
      dplyr::select(gene_info) %>%
      unique()
    
    if(this_gene_info$gene_info %in% colnames(germline_ptm_count_psort))
    {

      this_gene_germline_ptm = germline_ptm_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_ptm = this_gene_germline_ptm[,1]
    }else{
      this_germline_ptm = rep(0, nrow(germline_ptm_count_psort))
    }

    
    if(this_gene_info$gene_info %in% colnames(germline_domain_count_psort))
    {
      
      this_gene_germline_domain = germline_domain_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_domain = this_gene_germline_domain[,1]
    }else{
      this_germline_domain = rep(0, nrow(germline_domain_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_npc_count_psort))
    {
      
      this_gene_germline_npc = germline_npc_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_npc = this_gene_germline_npc[,1]
    }else{
      this_germline_npc = rep(0, nrow(germline_npc_count_psort))
    }
    
    count_variable = survival_count[,i]
    
    
    this_germline_piu = this_germline_ptm + this_germline_domain
    
   
    
    surv_model = coxph(surv_object ~  get_age + get_gender + get_race+ 
                         count_variable +
                         this_germline_piu +
                         this_germline_npc +
                         count_variable*this_germline_piu +
                         count_variable*this_germline_npc)
    
   
    no_interaction_surv_model = coxph(surv_object ~  get_age + get_gender + get_race+ 
                                        count_variable +
                                        this_germline_piu +
                                        this_germline_npc)
    
    ss = summary(surv_model)
    ss_coef = ss$coefficients
    
     niss = summary(no_interaction_surv_model)
    niss_coef = niss$coefficients
    
    this_surv_result = data.frame(count_info = colnames(survival_count)[i],
                                  num_patient = sum(count_variable!=0), 
                                  num_patient_gpiu = sum(this_germline_piu!=0),
                                  num_patient_gnpc = sum(this_germline_npc!=0),
                                  count_coeff = ss_coef[6,1], count_exp_coeff = ss_coef[6,2],count_pval = ss_coef[6,5],
                                  gpiu_coeff = ss_coef[7,1], gpiu_exp_coeff = ss_coef[7,2],gpiu_pval = ss_coef[7,5],
                                  gnpc_coeff = ss_coef[8,1], gnpc_exp_coeff = ss_coef[8,2],gnpc_pval = ss_coef[8,5],
                                  sgpiu_coeff = ss_coef[9,1], sgpiu_exp_coeff = ss_coef[9,2],sgpiu_pval = ss_coef[9,5],
                                  sgnpc_coeff = ss_coef[10,1], sgnpc_exp_coeff = ss_coef[10,2],sgnpc_pval = ss_coef[10,5],
                                 NIcount_coeff = niss_coef[6,1], NIcount_exp_coeff = niss_coef[6,2],NIcount_pval = niss_coef[6,5],
                                  NIgpiu_coeff = niss_coef[7,1], NIgpiu_exp_coeff = niss_coef[7,2],NIgpiu_pval = niss_coef[7,5],
                                  NIgnpc_coeff = niss_coef[8,1], NIgnpc_exp_coeff = niss_coef[8,2],NIgnpc_pval = niss_coef[8,5],
                                  stringsAsFactors = F)
    
    if(i%%100 == 0)
      cat(i, "\n")
    
    return(this_surv_result)
    
    
    
  }))
  
  
  
  
  
  
  q_val = qvalue(int_total_surv_result_df$count_pval)
 # bq_val = qvalue(int_total_surv_result_df$bcount_pval)
  
  
  niq_val = qvalue(int_total_surv_result_df$NIcount_pval)
  #bniq_val = qvalue(int_total_surv_result_df$bNIcount_pval)
  
  
  surv_df = int_total_surv_result_df %>%
    dplyr::mutate(count_qval = q_val$qvalues,
                 # bcount_qval =bq_val$qvalues,
                  NIcount_qval = niq_val$qvalues)%>%
                  #bNIcount_qval = bniq_val$qvalues) %>%
    dplyr::arrange(desc(num_patient)) %>%
    replace(is.na(.), "")
  
  write.table(surv_df,paste0(output_dir, output_name),
              quote = F,row.names = F, sep = "\t")
  
  
}



# with germline and bPIU --------------------------------------------------



fit_survival_model_with_bpiu_germline = function(interest_variable_info,
                                            unite_data,
                                            piu_gene_df,
                                            germline_domain_gene_level_count,
                                            germline_ptm_gene_level_count,
                                            germline_npc_count,
                                            somatic_bpiu_count,
                                            germline_bpiu_count,
                                            output_dir,
                                            output_name)
{
  
  
  # 
  # interest_variable_info = somatic_piu_unite[[3]]
  # unite_data = somatic_piu_unite[[1]]
  # piu_gene_df = somatic_piu_unite[[2]]
  # germline_domain_gene_level_count = germline_domain_gene_level_unite[[1]]
  # germline_ptm_gene_level_count = germline_ptm_gene_level_unite[[1]]
  # germline_npc_count = germline_npc_unite[[1]]
  # germline_bpiu_count = germline_bpiu_unite[[1]]
  # somatic_bpiu_count = somatic_bpiu_unite[[1]]
  # output_dir = output_dir
  # output_name = "somatic_domain_with_bpiu_germline_piu_npc.tsv"
  # 
   survival_data = unite_data %>%
    dplyr::filter(vital_status %in% c("Alive","Dead")) %>%
    dplyr::mutate(surv_status = dplyr::recode(vital_status,
                                              "Alive" = 0, 
                                              "Dead" = 1 ))
  
  surv_patient_df = data.frame(barcode = survival_data$barcode, stringsAsFactors = F)
  
  
  get_age = as.numeric(survival_data$age)
  get_gender = relevel(as.factor(survival_data$gender), ref = "MALE")
  get_race = relevel(as.factor(survival_data$crace), ref = "WHITE")
  
  surv_object = Surv(time = survival_data$surv_time, event =  survival_data$surv_status)
  
  
  survival_count = survival_data %>%
    dplyr::select(one_of(interest_variable_info))
  
  germline_ptm_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_ptm_gene_level_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  germline_domain_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_domain_gene_level_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  
  
  germline_npc_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_npc_count, by = "barcode")%>%
    replace(is.na(.), 0)
  
  
  germline_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_bpiu_count, by = "barcode")%>%
    replace(is.na(.), 0)
  
  
  somatic_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(somatic_bpiu_count, by = "barcode")%>%
    replace(is.na(.), 0)
  
  
  int_total_surv_result_df = rbindlist(lapply(1:ncol(survival_count), function(i)
    
  {
   # i = 1
    
    this_piu_info = colnames(survival_count)[i]
    this_gene_info = piu_gene_df %>%
      dplyr::filter(piu_info == this_piu_info) %>%
      dplyr::select(gene_info) %>%
      unique()
    
    if(this_gene_info$gene_info %in% colnames(germline_ptm_count_psort))
    {
      
      this_gene_germline_ptm = germline_ptm_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_ptm = this_gene_germline_ptm[,1]
    }else{
      this_germline_ptm = rep(0, nrow(germline_ptm_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_domain_count_psort))
    {
      
      this_gene_germline_domain = germline_domain_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_domain = this_gene_germline_domain[,1]
    }else{
      this_germline_domain = rep(0, nrow(germline_domain_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_npc_count_psort))
    {
      
      this_gene_germline_npc = germline_npc_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_npc = this_gene_germline_npc[,1]
    }else{
      this_germline_npc = rep(0, nrow(germline_npc_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_bpiu_count_psort))
    {
      
      this_gene_germline_bpiu = germline_bpiu_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_bpiu = this_gene_germline_bpiu[,1]
    }else{
      this_germline_bpiu = rep(0, nrow(germline_bpiu_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(somatic_bpiu_count_psort))
    {
      
      this_gene_somatic_bpiu = somatic_bpiu_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_somatic_bpiu = this_gene_somatic_bpiu[,1]
    }else{
      this_somatic_bpiu = rep(0, nrow(somatic_bpiu_count_psort))
    }
    
  
    count_variable = survival_count[,i]
    
    
    this_germline_piu = this_germline_ptm + this_germline_domain
    
    
    
    surv_model = coxph(surv_object ~  get_age + get_gender + get_race+ 
                         count_variable +
                         this_germline_piu +
                         this_germline_npc +
                         this_germline_bpiu +
                         this_somatic_bpiu +
                         count_variable*this_germline_piu +
                         count_variable*this_germline_npc +
                         count_variable*this_germline_bpiu +
                         count_variable*this_somatic_bpiu 
                         )
    
    
    no_interaction_surv_model = coxph(surv_object ~  get_age + get_gender + get_race+ 
                                        count_variable +
                                        this_germline_piu +
                                        this_germline_npc +
                                        this_germline_bpiu +
                                        this_somatic_bpiu)
    
    ss = summary(surv_model)
    ss_coef = ss$coefficients
    
    niss = summary(no_interaction_surv_model)
    niss_coef = niss$coefficients
    
    
    rs =  1+ 1 +length(unique(survival_data$crace)) -1 +1
    
    
    
    this_surv_result = data.frame(count_info = colnames(survival_count)[i],
                                  num_patient = sum(count_variable!=0), 
                                  num_patient_gpiu = sum(this_germline_piu!=0),
                                  num_patient_gnpc = sum(this_germline_npc!=0),
                                  count_coeff = ss_coef[rs,1], count_exp_coeff = ss_coef[rs,2],count_pval = ss_coef[rs,5],
                                  gpiu_coeff = ss_coef[rs+1,1], gpiu_exp_coeff = ss_coef[rs+1,2],gpiu_pval = ss_coef[rs+1,5],
                                  gnpc_coeff = ss_coef[rs+2,1], gnpc_exp_coeff = ss_coef[rs+2,2],gnpc_pval = ss_coef[rs+2,5],
                                  gbpiu_coeff = ss_coef[rs+3,1], gbpiu_exp_coeff = ss_coef[rs+3,2],gbpiu_pval = ss_coef[rs+3,5],
                                  sbpiu_coeff = ss_coef[rs+4,1], sbpiu_exp_coeff = ss_coef[rs+4,2],sbpiu_pval = ss_coef[rs+4,5],
                                  sgpiu_coeff = ss_coef[rs+5,1], sgpiu_exp_coeff = ss_coef[rs+5,2],sgpiu_pval = ss_coef[rs+5,5],
                                  sgnpc_coeff = ss_coef[rs+6,1], sgnpc_exp_coeff = ss_coef[rs+6,2],sgnpc_pval = ss_coef[rs+6,5],
                                  sgbpiu_coeff = ss_coef[rs+7,1], sgbpiu_exp_coeff = ss_coef[rs+7,2],sgbpiu_pval = ss_coef[rs+7,5],
                                  ssbpiu_coeff = ss_coef[rs+8,1], ssbpiu_exp_coeff = ss_coef[rs+8,2],ssbpiu_pval = ss_coef[rs+8,5],
                                  NIcount_coeff = niss_coef[rs,1], NIcount_exp_coeff = niss_coef[rs,2],NIcount_pval = niss_coef[rs,5],
                                  NIgpiu_coeff = niss_coef[rs+1,1], NIgpiu_exp_coeff = niss_coef[rs+1,2],NIgpiu_pval = niss_coef[rs+1,5],
                                  NIgnpc_coeff = niss_coef[rs+2,1], NIgnpc_exp_coeff = niss_coef[rs+2,2],NIgnpc_pval = niss_coef[rs+2,5],
                                  NIgbpiu_coeff = niss_coef[rs+3,1], NIgbpiu_exp_coeff = niss_coef[rs+3,2],NIgbpiu_pval = niss_coef[rs+3,5],
                                  NIsbpiu_coeff = niss_coef[rs+4,1], NIsbpiu_exp_coeff = niss_coef[rs+4,2],NIsbpiu_pval = niss_coef[rs+4,5],
                                  stringsAsFactors = F)
    
    
    
    
    if(i%%100 == 0)
      cat(i, "\n")
    
    return(this_surv_result)
    
    
    
  }))
  
  
  
  
  
  
  q_val = qvalue(int_total_surv_result_df$count_pval)
  # bq_val = qvalue(int_total_surv_result_df$bcount_pval)
  
  
  niq_val = qvalue(int_total_surv_result_df$NIcount_pval)
  #bniq_val = qvalue(int_total_surv_result_df$bNIcount_pval)
  
  
  surv_df = int_total_surv_result_df %>%
    dplyr::mutate(count_qval = q_val$qvalues,
                  # bcount_qval =bq_val$qvalues,
                  NIcount_qval = niq_val$qvalues)%>%
    #bNIcount_qval = bniq_val$qvalues) %>%
    dplyr::arrange(desc(num_patient)) %>%
    replace(is.na(.), "")
  
  write.table(surv_df,paste0(output_dir, output_name),
              quote = F,row.names = F, sep = "\t")
  
  
}



###


fit_survival_model_with_bpiu_germline_no_gender = function(interest_variable_info,
                                                 unite_data,
                                                 piu_gene_df,
                                                 germline_domain_gene_level_count,
                                                 germline_ptm_gene_level_count,
                                                 germline_npc_count,
                                                 somatic_bpiu_count,
                                                 germline_bpiu_count,
                                                 output_dir,
                                                 output_name)
{
  
  # interest_variable_info = stad_somatic_domain_unite[[3]]
  # unite_data = stad_somatic_domain_unite[[1]]
  # piu_gene_df = stad_somatic_domain_unite[[2]]
  # germline_domain_gene_level_count = stad_germline_domain_gene_level_unite[[1]]
  # germline_ptm_gene_level_count = stad_germline_ptm_gene_level_unite[[1]]
  # germline_npc_count = stad_germline_npc_unite[[1]]
  # germline_bpiu_count = stad_germline_bpiu_unite[[1]]
  # somatic_bpiu_count = stad_somatic_bpiu_unite[[1]]
  # output_dir = "/data/ginny/tcga_pancan/stad_somatic/stad_clinical_association_27/"
  # output_name = "somatic_domain_with_bpiu_germline_piu_npc.tsv"
  # 
  # 
  survival_data = unite_data %>%
    dplyr::filter(vital_status %in% c("Alive","Dead")) %>%
    dplyr::mutate(surv_status = dplyr::recode(vital_status,
                                              "Alive" = 0, 
                                              "Dead" = 1 ))
  
  surv_patient_df = data.frame(barcode = survival_data$barcode, stringsAsFactors = F)
  
  
  get_age = as.numeric(survival_data$age)
 # get_gender = relevel(as.factor(survival_data$gender), ref = "MALE")
  get_race = relevel(as.factor(survival_data$crace), ref = "WHITE")
  
  surv_object = Surv(time = survival_data$surv_time, event =  survival_data$surv_status)
  
  
  survival_count = survival_data %>%
    dplyr::select(one_of(interest_variable_info))
  
  germline_ptm_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_ptm_gene_level_count, by = "barcode") %>%
    replace(is.na(.),0)
  
  
  germline_domain_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_domain_gene_level_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  
  
  germline_npc_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_npc_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  germline_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_bpiu_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  somatic_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(somatic_bpiu_count, by = "barcode") %>%
    replace(is.na(.), 0)
  
  
  int_total_surv_result_df = rbindlist(lapply(1:ncol(survival_count), function(i)
    
  {
    # i = 1
    
    this_piu_info = colnames(survival_count)[i]
    this_gene_info = piu_gene_df %>%
      dplyr::filter(piu_info == this_piu_info) %>%
      dplyr::select(gene_info) %>%
      unique()
    
    if(this_gene_info$gene_info %in% colnames(germline_ptm_count_psort))
    {
      
      this_gene_germline_ptm = germline_ptm_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_ptm = this_gene_germline_ptm[,1]
    }else{
      this_germline_ptm = rep(0, nrow(germline_ptm_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_domain_count_psort))
    {
      
      this_gene_germline_domain = germline_domain_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_domain = this_gene_germline_domain[,1]
    }else{
      this_germline_domain = rep(0, nrow(germline_domain_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_npc_count_psort))
    {
      
      this_gene_germline_npc = germline_npc_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_npc = this_gene_germline_npc[,1]
    }else{
      this_germline_npc = rep(0, nrow(germline_npc_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(germline_bpiu_count_psort))
    {
      
      this_gene_germline_bpiu = germline_bpiu_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_germline_bpiu = this_gene_germline_bpiu[,1]
    }else{
      this_germline_bpiu = rep(0, nrow(germline_bpiu_count_psort))
    }
    
    
    if(this_gene_info$gene_info %in% colnames(somatic_bpiu_count_psort))
    {
      
      this_gene_somatic_bpiu = somatic_bpiu_count_psort %>%
        dplyr::select(this_gene_info$gene_info)
      this_somatic_bpiu = this_gene_somatic_bpiu[,1]
    }else{
      this_somatic_bpiu = rep(0, nrow(somatic_bpiu_count_psort))
    }
    
    
    count_variable = survival_count[,i]
    
    
    this_germline_piu = this_germline_ptm + this_germline_domain
    
    
    
    surv_model = coxph(surv_object ~  get_age +  get_race+ 
                         count_variable +
                         this_germline_piu +
                         this_germline_npc +
                         this_germline_bpiu +
                         this_somatic_bpiu +
                         count_variable*this_germline_piu +
                         count_variable*this_germline_npc +
                         count_variable*this_germline_bpiu +
                         count_variable*this_somatic_bpiu 
    )
    
    
    no_interaction_surv_model = coxph(surv_object ~  get_age  + get_race+ 
                                        count_variable +
                                        this_germline_piu +
                                        this_germline_npc +
                                        this_germline_bpiu +
                                        this_somatic_bpiu)
    
    ss = summary(surv_model)
    ss_coef = ss$coefficients
    
    niss = summary(no_interaction_surv_model)
    niss_coef = niss$coefficients
    
    
    ### count the number of rows before the counting variable 
    # , age, race(level)
    rs = 1+ length(unique(survival_data$crace)) -1 +1
    
    
    
    this_surv_result = data.frame(count_info = colnames(survival_count)[i],
                                  num_patient = sum(count_variable!=0), 
                                  num_patient_gpiu = sum(this_germline_piu!=0),
                                  num_patient_gnpc = sum(this_germline_npc!=0),
                                  count_coeff = ss_coef[rs,1], count_exp_coeff = ss_coef[rs,2],count_pval = ss_coef[rs,5],
                                  gpiu_coeff = ss_coef[rs+1,1], gpiu_exp_coeff = ss_coef[rs+1,2],gpiu_pval = ss_coef[rs+1,5],
                                  gnpc_coeff = ss_coef[rs+2,1], gnpc_exp_coeff = ss_coef[rs+2,2],gnpc_pval = ss_coef[rs+2,5],
                                  gbpiu_coeff = ss_coef[rs+3,1], gbpiu_exp_coeff = ss_coef[rs+3,2],gbpiu_pval = ss_coef[rs+3,5],
                                  sbpiu_coeff = ss_coef[rs+4,1], sbpiu_exp_coeff = ss_coef[rs+4,2],sbpiu_pval = ss_coef[rs+4,5],
                                  sgpiu_coeff = ss_coef[rs+5,1], sgpiu_exp_coeff = ss_coef[rs+5,2],sgpiu_pval = ss_coef[rs+5,5],
                                  sgnpc_coeff = ss_coef[rs+6,1], sgnpc_exp_coeff = ss_coef[rs+6,2],sgnpc_pval = ss_coef[rs+6,5],
                                  sgbpiu_coeff = ss_coef[rs+7,1], sgbpiu_exp_coeff = ss_coef[rs+7,2],sgbpiu_pval = ss_coef[rs+7,5],
                                  ssbpiu_coeff = ss_coef[rs+8,1], ssbpiu_exp_coeff = ss_coef[rs+8,2],ssbpiu_pval = ss_coef[rs+8,5],
                                  NIcount_coeff = niss_coef[rs,1], NIcount_exp_coeff = niss_coef[rs,2],NIcount_pval = niss_coef[rs,5],
                                  NIgpiu_coeff = niss_coef[rs+1,1], NIgpiu_exp_coeff = niss_coef[rs+1,2],NIgpiu_pval = niss_coef[rs+1,5],
                                  NIgnpc_coeff = niss_coef[rs+2,1], NIgnpc_exp_coeff = niss_coef[rs+2,2],NIgnpc_pval = niss_coef[rs+2,5],
                                  NIgbpiu_coeff = niss_coef[rs+3,1], NIgbpiu_exp_coeff = niss_coef[rs+3,2],NIgbpiu_pval = niss_coef[rs+3,5],
                                  NIsbpiu_coeff = niss_coef[rs+4,1], NIsbpiu_exp_coeff = niss_coef[rs+4,2],NIsbpiu_pval = niss_coef[rs+4,5],
                                  stringsAsFactors = F)
    
    if(i%%100 == 0)
      cat(i, "\n")
    
    return(this_surv_result)
    
    
    
  }))
  
  
  
  
  
  
  q_val = qvalue(int_total_surv_result_df$count_pval)
  # bq_val = qvalue(int_total_surv_result_df$bcount_pval)
  
  
  niq_val = qvalue(int_total_surv_result_df$NIcount_pval)
  #bniq_val = qvalue(int_total_surv_result_df$bNIcount_pval)
  
  
  surv_df = int_total_surv_result_df %>%
    dplyr::mutate(count_qval = q_val$qvalues,
                  # bcount_qval =bq_val$qvalues,
                  NIcount_qval = niq_val$qvalues)%>%
    #bNIcount_qval = bniq_val$qvalues) %>%
    dplyr::arrange(desc(num_patient)) %>%
    replace(is.na(.), "")
  
  write.table(surv_df,paste0(output_dir, output_name),
              quote = F,row.names = F, sep = "\t")
  
  
}



