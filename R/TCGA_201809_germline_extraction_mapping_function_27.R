# 
# library(data.table)
# library(dplyr)
# library(magrittr)





parse_row_col = function(s_number, total_row, total_col)
{
  
  if(s_number > total_row*total_col)
  {
    stop("serial number not match!")
    
  }else{
    row_col = rep(0,2)
    
    r = s_number%%total_row
    
    if(r==0)
    {
      row_col[1] = total_row
    }else
    {
      row_col[1] = r
    }
    
    c = floor(s_number/total_row) +1
    
    row_col[2] = c
    
    
    return(row_col)
    
  }
  
}

# I need to get more info out of this 

parse_annotation = function(whole_INFO)
{
  ann_df = rbindlist(lapply(1:length(whole_INFO), function(x) {
    ta = unlist(strsplit(whole_INFO[x], split = "ANN="))[2]
    get_split = unlist(strsplit(ta, split = "|", fixed = T))
    this_df = data.frame(allele = get_split[1],
                         annotation = get_split[2],
                         annotation_impact = get_split[3],
                         gene_name = get_split[4],
                         gene_id = get_split[5],
                         feature_type = get_split[6],
                         feature_id = get_split[7],
                         transcript_biotype = get_split[8],
                         rank = get_split[9],
                         hgvs.c = get_split[10],
                         hgvs.p = get_split[11],
                         cdna.pos_length = get_split[12],
                         cds.pos_length = get_split[13],
                         aa.pos_length = get_split[14],
                         distance = get_split[15],
                         #    error_warnings_info = get_split[16],
                         stringsAsFactors = F)
    
    return(this_df)
    
  }))
  
  return(ann_df)
  
}



get_cancer_type_col_id = function(sample_cn_id, cancer_barcode)
{
  
  which_col =  sample_cn_id %in% cancer_barcode
  
  return(which_col)
  
}



# get whole mutation data for each cancer type ----------------------------------

snpeff_record_whole_mutation_in_cancer = function(all_names, chr,
                                                  cancer_type,
                                                  cancer_barcode, sample_cn_id,
                                                  output_dir)
  
{
  to_grep = paste0("canon_", chr, "_")
  
  this_chr_names = grep(to_grep, all_names, value = T)
  
  tmp = gsub("snpeff_annotation_canon","mut_store",this_chr_names)
  
  this_chr_mut_file_names = gsub(".vcf", ".Rds", tmp)
  
  
  whole_chr_type =  rbindlist(lapply(1:length(this_chr_names), function(i){
    
    # i = 1
    
    this_mut_dir = paste0("/data/ginny/tcga_pancan/germline_raw_process/mut_files/", 
                          this_chr_mut_file_names[i])
    this_mut = readRDS(this_mut_dir)
    
    ts = this_mut[[2]]
    tr = this_mut[[3]][1]
    tc = this_mut[[3]][2]
    
    rc_in_mut = unlist(lapply(1:length(ts), function(x)
    {
      row_col =   parse_row_col(ts[x], tr, tc)
      
      return(row_col)
      
    }))
    
    rc_row = rc_in_mut[c(T,F)]
    rc_col = rc_in_mut[c(F,T)]
    
    rc_df = data.frame(mut_row = rc_row, sample_col = rc_col, stringsAsFactors = F)
    
    
    this_piece_dir = paste0("/data/ginny/tcga_pancan/germline_raw_process/snpeff_annotation_output_files/", this_chr_names[i])
    this_piece = fread(this_piece_dir, stringsAsFactors = F)
    colnames(this_piece)[1] = "CHROM"
    
    info_part = parse_annotation(this_piece$INFO) 
    whole_snpeff = cbind(this_piece, info_part)
    
    sel_col = seq(1:tc)[get_cancer_type_col_id(sample_cn_id, cancer_barcode)]
    
    whole_type_snpeff = whole_snpeff %>%
      dplyr::select(CHROM,POS,ID,QUAL, FORMAT, annotation, gene_name, gene_id, hgvs.c, hgvs.p, aa.pos_length) %>%
      unique()
    
    
    sel_row = whole_type_snpeff$ID
    rc_in_mut_sel = rc_df%>%
      dplyr::filter(mut_row %in% sel_row, sample_col %in% sel_col)
    
    
    sel_whole_type_snpeff = whole_type_snpeff %>%
      dplyr::filter(ID %in% rc_in_mut_sel$mut_row)
    
    
    
    rc_agg = rc_in_mut_sel %>%
      dplyr::group_by(mut_row) %>%
      dplyr::mutate(agg_sample_id = paste0(sample_col, collapse = "_"))%>%
      dplyr::mutate(mut_freq = n())%>%
      dplyr::select(mut_row, agg_sample_id, mut_freq) %>%
      unique()
    
    
    info_agg = sel_whole_type_snpeff %>%
      dplyr::left_join(rc_agg, by = c("ID" = "mut_row")) %>%
      na.omit()
    
    
    
    cat(i, " nrow ", nrow(info_agg), " ~~~~~~~~~~~~~~", "\n" )
    
    return(info_agg)
    
  }))
  
  
  out_name = paste0("whole_types_chr_",chr, "_", cancer_type, "_snpeff_mutations.tsv")
  write.table(whole_chr_type,paste0(output_dir,out_name),
              sep = "\t", quote = F, row.names = F)
}



# combine all chromosome together -----------------------------------------



combine_chr_data = function(dir_to_combine, quality_filter, dir_for_output, output_name)
{
  
  to_combine_file_names  = list.files(dir_to_combine)
  
  combine_df = rbindlist(lapply(1:length(to_combine_file_names), function(x)
    #combine_df = rbindlist(lapply(1:3, function(x)
  {
    
    this_tsv = fread(paste0(dir_to_combine, to_combine_file_names[x]),
                     stringsAsFactors = F)
    
    return(this_tsv)
    
    
  }))
  
  write.table(combine_df, paste0(dir_for_output, output_name), 
              quote = F, row.names = F, sep = "\t")
  
  high_qual_combine_df = combine_df %>%
    dplyr::filter(QUAL >= quality_filter)
  
  #hq_name = paste0("higher_",quality_filter,"_",output_name)
  write.table(high_qual_combine_df, paste0(dir_for_output, output_name), 
              quote = F, row.names = F, sep = "\t")
  
  
}



# divide into pc and npc --------------------------------------------------



divide_germline_to_pc_npc = function(combine_data_name,
                                     output_dir,
                                     pc_output_name,
                                     npc_output_name)
{
  
  
  combine_data = fread(combine_data_name, stringsAsFactors = F)
  
  ### first, remove synonymous mutation 
  
  
  pc_data = combine_data %>%
    dplyr::filter(!grepl("synonymous", annotation)) %>%
    dplyr::filter(nchar(hgvs.p)>0) 
  
  
  npc_data = combine_data %>%
    dplyr::filter(!grepl("synonymous", annotation)) %>%
    dplyr::filter(nchar(hgvs.p)==0) 
  
  
  write.table(pc_data, paste0(output_dir, pc_output_name), 
              quote = F, row.names = F, sep = "\t")
  
  
  write.table(npc_data, paste0(output_dir, npc_output_name), 
              quote = F, row.names = F, sep = "\t")
  
  
  
}



# annotate position for pc ------------------------------------------------


annotate_snpeff_combine_pc_position_info= function(pc_data_name, 
                                                      output_dir, output_name)
{
  
  pc_data = fread(pc_data_name, stringsAsFactors = F)
  
  ### ok , for snpeff there is no cofusion, each mutation only has one annotation  cool
  
  p_pos =  gsub("[[:alpha:]]","", pc_data$hgvs.p)
  fp_pos = gsub("\\.","",p_pos)
  start_end_pos = rbindlist(lapply(1:length(fp_pos), function(t){
    
    get_it = unlist(strsplit(fp_pos[t], split = "_"))
    
    start_end_df = data.frame(pos_start = 0, pos_end = 0)
    
    get_it[1] = gsub("\\D", "", get_it[1])
    
    
    start_end_df$pos_start = as.integer(get_it[1])
    
    
    if(length(get_it)!=2)
    {
      start_end_df$pos_end =  as.integer(get_it[1])
    }else{
      get_it[2] = gsub("\\D", "", get_it[2])
      start_end_df$pos_end =  as.integer(get_it[2])
    }
    
    return(start_end_df)
  }
  ))
  ### I need to get the protein ID mappable 
  
  sel_pc = pc_data %>%
    dplyr::mutate(prot_start_pos = start_end_pos$pos_start, 
                  prot_end_pos = start_end_pos$pos_end)%>%
    dplyr::select(CHROM, POS, QUAL,ID, FORMAT, annotation, gene_name, gene_id, 
                  hgvs.p, prot_start_pos, prot_end_pos, agg_sample_id, mut_freq)%>%
    unique()
  

  write.table(sel_pc, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}


# PIU mapping and bPIU summarising for pc ---------------------------------




snpeff_combine_map_uni_piu = function(ptm_domain_filename,
                                      pc_data_name,
                                      barcode_seq_df,
                                      output_dir,
                                      piu_output_filename,
                                      bpiu_output_filename)
  
  
 {
#   ptm_domain_filename = "/data/ginny/tcga_pancan/important_files/ptm_domain_combine_df.tsv"
#   pc_data_name = "/data/ginny/tcga_pancan/germline_raw_process/stad_germline_variant_annotation/higher_30_snpeff_variant_anno_pc_pos.tsv"
#   barcode_seq_df = stad_barcode_seq_df
#   output_dir  = "/data/ginny/tcga_pancan/germline_raw_process/stad_summarise_mutation/"
#   piu_output_filename = "piu_mapping_count.tsv"
#   bpiu_output_filename = "bpiu_summarising_count.tsv"
  
  
  
  
  
  pc_mut = fread(pc_data_name, stringsAsFactors = F)
  piu = fread(ptm_domain_filename, stringsAsFactors = F)%>%
    na.omit()
  
  unique_prot = intersect(piu$gene_id, pc_mut$gene_id)
  unique_sample = barcode_seq_df$cancer_code
  
  watch_piu  = piu %>%
    dplyr::filter(gene_id %in% unique_prot)
  
  if(length(unique_prot) == 0)
  {
    cat("no single protien mapped", "\n")
  }else{
    
   
piu_df = rbindlist(lapply(1:length(unique_prot), function(x)
     #piu_df = rbindlist(lapply(1:100, function(x)
        
    {

      # x=  1

      this_prot = unique_prot[x]
      this_prot_piu = piu%>%
        dplyr::filter(gene_id == this_prot)



      count_piu_matrix = matrix(0, nrow = nrow(this_prot_piu), ncol = length(unique_sample))
      colnames(count_piu_matrix) = unique_sample

      local = pc_mut %>%
        dplyr::filter(gene_id == this_prot)

      if(nrow(local)>0)
      {
        for(i in 1:nrow(local))
        {
          # i = 1
          get_sample = data.frame(this_seq = unlist(strsplit(local$agg_sample_id[i], split = "_")), stringsAsFactors = F)
          get_sample_barcode = get_sample %>%
            dplyr::left_join(barcode_seq_df, by = c("this_seq" = "cancer_seq")) %>%
            na.omit()%>%
            dplyr::select(cancer_code)

          these_patient_col = which(unique_sample %in% get_sample_barcode$cancer_code)


          piu_row_flag <- (local$prot_start_pos[i]>=this_prot_piu$start_position & local$prot_start_pos[i]<= this_prot_piu$end_position |local$prot_end_pos[i]>= this_prot_piu$start_position & local$prot_end_pos[i]<= this_prot_piu$end_position)


          count_piu_matrix[piu_row_flag, these_patient_col] =
            count_piu_matrix[piu_row_flag, these_patient_col] +1

          # if(sum(piu_row_flag)==0)
          # {
          #   count_bpiu_matrix[x, these_patient_col] = count_bpiu_matrix[x, these_patient_col] +1
          # }


        }
      }

      row_sum = rowSums(count_piu_matrix)


      this_prot_local_df = cbind(this_prot_piu, count_piu_matrix, row_sum)

    if(x%%100 == 0)
        cat(x, "\n")

      return(this_prot_local_df)

    }))
      
    with_info_piu_df = piu_df %>%
      dplyr::filter(row_sum>0) %>%
      dplyr::arrange(desc(row_sum))
      write.table(with_info_piu_df, paste0(output_dir, piu_output_filename), quote = F, row.names = F, sep = "\t")
 
      
       #####################################################
      
      bpiu_df = rbindlist(lapply(1:length(unique_prot), function(x)
        
      {
        
        # x=  1
        
        this_prot = unique_prot[x]
        this_prot_piu = piu%>%
          dplyr::filter(gene_id == this_prot)
        
        this_bpiu_info = this_prot_piu %>%
          dplyr::select(gene_id,gene_name)%>%
          unique()%>%
          dplyr::top_n(n = 1, wt = gene_name)
        
        
        count_bpiu_matrix = matrix(0, nrow = 1, ncol = length(unique_sample))
        colnames(count_bpiu_matrix) = unique_sample
        
        local = pc_mut %>%
          dplyr::filter(gene_id == this_prot)
        
        if(nrow(local)>0)
        {
          for(i in 1:nrow(local))
          {
            # i = 1
            get_sample = data.frame(this_seq = unlist(strsplit(local$agg_sample_id[i], split = "_")), stringsAsFactors = F)
            get_sample_barcode = get_sample %>%
              dplyr::left_join(barcode_seq_df, by = c("this_seq" = "cancer_seq")) %>%
              na.omit()%>%
              dplyr::select(cancer_code)
            
            these_patient_col = which(unique_sample %in% get_sample_barcode$cancer_code)
            
            
            piu_row_flag <- (local$prot_start_pos[i]>=this_prot_piu$start_position & local$prot_start_pos[i]<= this_prot_piu$end_position |local$prot_end_pos[i]>= this_prot_piu$start_position & local$prot_end_pos[i]<= this_prot_piu$end_position)
            
          
            if(sum(piu_row_flag)==0)
            {
              count_bpiu_matrix[1, these_patient_col] = count_bpiu_matrix[1, these_patient_col] +1
            }
            
            
          }
        }
        
        row_sum = rowSums(count_bpiu_matrix)
        
        
        this_prot_bpiu_df = cbind(this_bpiu_info, count_bpiu_matrix, row_sum)
        
       if(x%%100 == 0)
        cat(x, "\n")
        
        return(this_prot_bpiu_df)
        
      }))
      
      with_info_bpiu_df = bpiu_df %>%
        dplyr::filter(row_sum>0) %>%
        dplyr::arrange(desc(row_sum))
      write.table(with_info_bpiu_df, paste0(output_dir, bpiu_output_filename), quote = F, row.names = F, sep = "\t")
      
      
      
  }
  
}

###next step is to get ncu

# npc summarising ---------------------------------------------------------



snpeff_combine_map_npc = function(npc_data_name,
                                  barcode_seq_df,
                                  output_dir,
                                  output_filename)
  
  
{ 
  
  # npc_data_name = "/data/ginny/tcga_pancan/germline_raw_process/stad_germline_variant_annotation/higher_30_snpeff_variant_anno_npc.tsv"
  # barcode_seq_df = stad_barcode_seq_df
  # output_dir  = "/data/ginny/tcga_pancan/germline_raw_process/stad_summarise_mutation/"
  # output_filename = "npc_summarising_count.tsv"
  # 
  
  non_impact = fread(npc_data_name, stringsAsFactors = F)
  gene = non_impact %>%
    dplyr::select(gene_name, gene_id)%>%
    unique()
  
  
  unique_sample = barcode_seq_df$cancer_code
  
  if(nrow(gene) == 0)
  {
    cat("no single protein mapped", "\n")
  }else{
    
    gene_df = rbindlist(lapply(1:nrow(gene), function(x)
    {
      
      this_gene_id = gene$gene_id[x]
      
      
      count_gene = gene[x,]
      
      count_gene_matrix = matrix(0, nrow = 1, ncol = length(unique_sample))
      colnames(count_gene_matrix) = unique_sample
      
      this_non_impact = non_impact %>%
        dplyr::filter(gene_id == this_gene_id)
      
      if(nrow(this_non_impact)>0)
      {
        for(i in 1:nrow(this_non_impact))
        {
          # i = 1
          get_sample = data.frame(this_seq = unlist(strsplit(this_non_impact$agg_sample_id[i], split = "_")), stringsAsFactors = F)
          get_sample_barcode = get_sample %>%
            dplyr::left_join(barcode_seq_df, by = c("this_seq" = "cancer_seq")) %>%
            na.omit()%>%
            dplyr::select(cancer_code)
          
          these_patient_col = which(colnames(count_gene_matrix) %in% get_sample_barcode$cancer_code)
          
          count_gene_matrix[1, these_patient_col] = 
            count_gene_matrix[1, these_patient_col] +1
          
        }
      }
      
      row_sum = rowSums(count_gene_matrix)
      
      
      this_gene_df = cbind(count_gene, count_gene_matrix, row_sum)
      
      if(x%%100 == 0)
        cat(x, "\n")
      
      return(this_gene_df)
      
    }))
    
    
    with_info_gene_df = gene_df %>%
      dplyr::filter(row_sum>0) %>%
      dplyr::arrange(desc(row_sum))
    
    
    write.table(with_info_gene_df, paste0(output_dir, output_filename), quote = F, row.names = F, sep = "\t")
    
    
  }
  
}



