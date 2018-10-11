### functions for association 


# get the piu counts and clinical info united before modelling ----------------


piu_counts_clinical_unite = function(piu_count_filename,
                                     patient_info_df,
                                     patient_outcome_df,
                                     piu_of_interest,
                                     row_sum_min)
  
{
  
  # piu_count_filename = "/data/ginny/tcga_pancan/stad_somatic/stad_summarise_mutation/piu_mapping_count.tsv"
  # patient_info_df = patient_info
  # patient_outcome_df = patient_outcome
  # piu_of_interest = "domain"
  # row_sum_min = 1
  # 
  # piu_count_filename = "/data/ginny/tcga_pancan/germline_raw_process/stad_summarise_mutation/piu_mapping_count.tsv"
  # patient_info_df = patient_info
  # patient_outcome_df = patient_outcome
  # piu_of_interest = "domain"
  # row_sum_min = 1
  # 
  piu_count_df = fread(piu_count_filename,
                       stringsAsFactors = F)
  
  piu_count_sel = piu_count_df%>%
    dplyr::filter(unit_label == piu_of_interest)%>%
    dplyr::mutate(piu_info = paste(uniprot_accession, start_position,end_position, unit_name, gene_name, gene_id, sep = "_"))%>%
    dplyr::mutate(gene_info = paste(gene_name, gene_id,sep = "_"))%>%
    dplyr::filter(row_sum >= row_sum_min) %>%
    dplyr::select(uniprot_accession, start_position, end_position, unit_name, gene_name, gene_id, unit_label,gene_info, piu_info, row_sum,
                  everything())
  
  #piu_info = piu_count_sel$piu_info
  
  
  which_barcode = grep("TCGA", colnames(piu_count_sel))
  
  piu_matrix = t(as.matrix(piu_count_sel[,which_barcode]))
  
  piu_count = data.frame(barcode = colnames(piu_count_sel)[12:454],
                         piu_matrix,
                         stringsAsFactors = F)
  
  colnames(piu_count) = c("barcode", piu_count_sel$piu_info)
  rownames(piu_count) = NULL
  
  piu_gene_df = data.frame(piu_info = piu_count_sel$piu_info,
                           gene_info = piu_count_sel$gene_info,
                           stringsAsFactors = F)
  
  
  
  
  piu_clinical_unite_data = piu_count %>%
    dplyr::left_join(patient_outcome, by = c("barcode" = "bcr_patient_barcode")) %>%
    na.omit() %>%
    dplyr::left_join(patient_info, by  = c("barcode" = "bcr_patient_barcode"))%>%
    na.omit() %>% 
    dplyr::arrange(desc(surv_time))%>%
    dplyr::arrange(vital_status) 
  
  
  return_list = vector(mode = "list", length = 4)
  
  return_list[[1]] = piu_clinical_unite_data
  return_list[[2]] = piu_gene_df
  return_list[[3]] = piu_count_sel$piu_info
  return_list[[4]] = piu_count_sel$gene_id
  
  
  
  
  
  return(return_list)
  
}


# get the germline npc counts and clinical info united before modelling ----------------




gene_level_counts_clinical_unite = function(gene_level_count_filename,
                                              quiry_gene_id,
                                              patient_info_df,
                                              patient_outcome_df,
                                              row_sum_min)
  
{
  
  # gene_level_count_filename = "/data/ginny/tcga_pancan/raw_process/stad_summarise_mutation/gene_level_summarising_count.tsv"
  # quiry_gene_id = unique(stad_somatic_domain_unite[[4]])
  # patient_info_df = patient_info
  # patient_outcome_df = patient_outcome
  # row_sum_min = 1
  # 
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
    dplyr::left_join(patient_outcome, by = c("barcode" = "bcr_patient_barcode")) %>%
    na.omit() %>%
    dplyr::left_join(patient_info, by  = c("barcode" = "bcr_patient_barcode"))%>%
    na.omit() %>% 
    dplyr::arrange(desc(surv_time))%>%
    dplyr::arrange(vital_status) 
  
  
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




fit_survival_model = function(interest_variable_info,
                              unite_data,
                              output_dir,
                              output_name)
  
{
  
  # interest_variable_info = stad_germline_npc_unite[[2]]
  # unite_data = stad_germline_npc_unite[[1]]
  # output_dir = "/data/ginny/tcga_pancan/stad_somatic/stad_clinical_association_27/"
  # output_name = "univariate_germline_npc.tsv"
  # 

  survival_data = unite_data %>%
    dplyr::filter(vital_status %in% c("Alive","Dead")) %>%
    dplyr::mutate(surv_status = dplyr::recode(vital_status,
                                              "Alive" = 0, 
                                              "Dead" = 1 ))
  
  
  get_age = as.numeric(survival_data$age)
  get_gender = relevel(as.factor(survival_data$gender), ref = "MALE")
  get_race = relevel(as.factor(survival_data$crace), ref = "WHITE")
  
  
  survival_count = survival_data %>%
    dplyr::select(one_of(interest_variable_info))
  
  
  surv_object = Surv(time = survival_data$surv_time, event =  survival_data$surv_status)
  
  
  #for(i in 1:ncol(survival_count))
  #{
  
  surv_result_df = rbindlist(lapply(1:ncol(survival_count), function(i)
  {
    
    
    #i = 4345
    count_variable = survival_count[,i]
    
    #count_variable[which(count_variable>0)]=1
    
    surv_model = coxph(surv_object ~  get_age + get_gender + get_race+ 
                         count_variable)
    
    ss = summary(surv_model)
    ss_coef = ss$coefficients
    
    this_surv_result = data.frame(count_info = colnames(survival_count)[i],
                                  num_patients = sum(count_variable!=0),
                                  count_coeff = ss_coef[6,1], count_exp_coeff = ss_coef[6,2],
                                  count_pval = ss_coef[6,5],
                                  stringsAsFactors = F)
    
    if(i%%100 == 0 )
      cat(i, "\n")
    
    return(this_surv_result)
    
  }))
  
  
  # q_val = p.adjust(surv_result_df$count_pval,method = "BH")
  
  
  q_val = qvalue(surv_result_df$count_pval)
  
  surv_df = surv_result_df %>%
    dplyr::mutate(count_qval = q_val$qvalues) %>%
    na.omit()%>%
    dplyr::arrange(desc(num_patients)) 
  
  
  
  #}
  
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
    dplyr::left_join(germline_ptm_gene_level_count, by = "barcode")

  
  germline_domain_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_domain_gene_level_count, by = "barcode")
  
  
  
  
  germline_npc_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_npc_count, by = "barcode")
  
  
  
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
  get_gender = relevel(as.factor(survival_data$gender), ref = "MALE")
  get_race = relevel(as.factor(survival_data$crace), ref = "WHITE")
  
  surv_object = Surv(time = survival_data$surv_time, event =  survival_data$surv_status)
  
  
  survival_count = survival_data %>%
    dplyr::select(one_of(interest_variable_info))
  
  germline_ptm_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_ptm_gene_level_count, by = "barcode")
  
  
  germline_domain_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_domain_gene_level_count, by = "barcode")
  
  
  
  
  germline_npc_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_npc_count, by = "barcode")
  
  
  germline_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(germline_bpiu_count, by = "barcode")
  
  
  somatic_bpiu_count_psort = surv_patient_df %>%
    dplyr::left_join(somatic_bpiu_count, by = "barcode")
  
  
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
    
    this_surv_result = data.frame(count_info = colnames(survival_count)[i],
                                  num_patient = sum(count_variable!=0), 
                                  num_patient_gpiu = sum(this_germline_piu!=0),
                                  num_patient_gnpc = sum(this_germline_npc!=0),
                                  count_coeff = ss_coef[6,1], count_exp_coeff = ss_coef[6,2],count_pval = ss_coef[6,5],
                                  gpiu_coeff = ss_coef[7,1], gpiu_exp_coeff = ss_coef[7,2],gpiu_pval = ss_coef[7,5],
                                  gnpc_coeff = ss_coef[8,1], gnpc_exp_coeff = ss_coef[8,2],gnpc_pval = ss_coef[8,5],
                                  gbpiu_coeff = ss_coef[9,1], gbpiu_exp_coeff = ss_coef[9,2],gbpiu_pval = ss_coef[9,5],
                                  sbpiu_coeff = ss_coef[10,1], sbpiu_exp_coeff = ss_coef[10,2],sbpiu_pval = ss_coef[10,5],
                                  sgpiu_coeff = ss_coef[11,1], sgpiu_exp_coeff = ss_coef[11,2],sgpiu_pval = ss_coef[11,5],
                                  sgnpc_coeff = ss_coef[12,1], sgnpc_exp_coeff = ss_coef[12,2],sgnpc_pval = ss_coef[12,5],
                                  sgbpiu_coeff = ss_coef[13,1], sgbpiu_exp_coeff = ss_coef[13,2],sgbpiu_pval = ss_coef[13,5],
                                  ssbpiu_coeff = ss_coef[14,1], ssbpiu_exp_coeff = ss_coef[14,2],ssbpiu_pval = ss_coef[14,5],
                                  NIcount_coeff = niss_coef[6,1], NIcount_exp_coeff = niss_coef[6,2],NIcount_pval = niss_coef[6,5],
                                  NIgpiu_coeff = niss_coef[7,1], NIgpiu_exp_coeff = niss_coef[7,2],NIgpiu_pval = niss_coef[7,5],
                                  NIgnpc_coeff = niss_coef[8,1], NIgnpc_exp_coeff = niss_coef[8,2],NIgnpc_pval = niss_coef[8,5],
                                  NIgbpiu_coeff = niss_coef[9,1], NIgbpiu_exp_coeff = niss_coef[9,2],NIgbpiu_pval = niss_coef[9,5],
                                  NIsbpiu_coeff = niss_coef[10,1], NIsbpiu_exp_coeff = niss_coef[10,2],NIsbpiu_pval = niss_coef[10,5],
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




