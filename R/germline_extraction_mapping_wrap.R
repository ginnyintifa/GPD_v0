# 
# library(data.table)
# library(dplyr)
# library(magrittr)






#' Extract germline mutations annotated by the snpEff tool for a cancer type of your interest.
#' 
#' This function generates two files as the final products. First, mutations which result in protein-consequence with
#' effect on proteins annotated; Second, mutations without protein-consequence.
#' 
#' @param raw_snpeff_output_dir Directory where raw snpeff outputs are stored。 
#' @param cancer_type TCGA cancer type identifier, e.g. STAD, BRCA. 
#' @param cancer_barcode TCGA barcodes for this cancer type cohort.
#' @param sample_cn_id A barcode codebook provided by the packge.
#' @param quality_filter The minimiu value of QUAL for selection of mutations.
#' @param original_dir The parent directory you would like to have all you output files in.
#' @param output_dir The directory you would like to have your raw extraction files in.
#' @import dplyr magrittr data.table
#' @export
#' @details 
#' @examples 
#'germline_extraction_annotation_pos (raw_snpeff_output_dir = "/data/ginny/tcga_pancan/germline_raw_process/snpeff_annotation_output_files",
#'                                    cancer_type = "STAD,    
#'                                    cancer_barcode = stad_barcode,
#'                                    sample_cn_id = sample_cn_id,
#'                                    quality_filter = 30,
#'                                    original_dir = "/data/ginny/tcga_pancan/germline_raw_process/" ,
#'                                    output_dir = "/data/ginny/tcga_pancan/germline_raw_process/stad_snpeff_type_muts/")




germline_extraction_annotation_pos = function(raw_snpeff_output_dir,
                               cancer_type,
                               cancer_barcode,
                               sample_cn_id,
                               quality_filter = 30,
                               original_dir,
                               output_dir)
{
  
  
  
  
  
  snpeff_output_file_names = list.files(raw_snpeff_output_dir)
  
  
  system(paste0("mkdir ", output_dir))
  
  
  chr_name = c(seq(1:22),"X","Y")
  
  for(t in 1:length(chr_name))
  {
    
    snpeff_record_whole_mutation_in_cancer (all_names = snpeff_output_file_names,
                                            chr = chr_name[t] ,
                                            cancer_type = cancer_type,
                                            cancer_barcode = cancer_barcode,
                                            sample_cn_id = sample_cn_id,
                                            output_dir = output_dir
    )
    
    
    cat("finish chr ", chr_name[t], "\n")
  }
  
    
    combine_chr_data_output_dir = paste0(original_dir,cancer_type,"_snpeff_annotation/")

    system(paste0("mkdir ", combine_chr_data_output_dir))
    
    
    
    combine_chr_data(dir_to_combine = output_dir,
                     quality_filter = quality_filter,
                     dir_for_output = combine_chr_data_output_dir,
                     output_name = "snpeff_variant_anno_combine.tsv")
    
    
    
    divide_germline_to_pc_npc(combine_data_name = paste0(combine_chr_data_output_dir,"snpeff_variant_anno_combine.tsv"),
                              output_dir = combine_chr_data_output_dir,
                              pc_output_name = "snpeff_variant_anno_pc.tsv",
                              npc_output_name = "snpeff_variant_anno_npc.tsv")
    
    
    
    
    annotate_snpeff_combine_pc_position_info(pc_data_name = paste0(combine_chr_data_output_dir, "snpeff_variant_anno_pc.tsv"),
                                             output_dir = combine_chr_data_output_dir,
                                             output_name = "snpeff_variant_anno_pc_pos.tsv")
    
    

 cat("Germline extraction and annotation finished!","\n")
 
}









