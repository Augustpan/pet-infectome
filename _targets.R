library(targets)

# libraries used in functions
library(tidyverse)
library(readxl)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(cowplot)
library(nlme)
library(vegan)
library(ggvenn)
library(DESeq2)
library(Cairo)
library(coin)

# load functions
source("R/functions.R")
source("R/functions_production.R")

# create an output folder
if (!dir.exists("output"))
  dir.create("data")

# libraries used in functions
tar_option_set(packages = c("tidyverse", "readxl", "circlize", "RColorBrewer", 
                            "ComplexHeatmap", "cowplot", "nlme", "vegan",
                            "ggvenn", "DESeq2", "Cairo", "coin"))

# targets list
list(
  # define paths to data files
  tar_target(file_lib_meta_raw, "raw_data/lib_meta.csv", format = "file"),
  tar_target(file_fastq_manifest, "raw_data/fastq_manifest.txt", format = "file"),
  tar_target(file_species_info, "raw_data/species_identification.txt", format = "file"),
  tar_target(file_raw_vir_profile_amy, "raw_data/Virus_profiling_amy.xlsx", format = "file"),
  tar_target(file_raw_bac_profile_amy, "raw_data/Bacteria_mataphlan_profiling_amy.xlsx", format = "file"),
  #tar_target(file_raw_bac_profile_genus_pan, "raw_data/bacteria_genus_profile_pan.csv", format = "file"),
  #tar_target(file_raw_fun_profile_genus_pan, "raw_data/fungi_genus_profile_pan.csv", format = "file"),
  tar_target(file_fun_missing, "raw_data/fun_missing.txt", format = "file"),
  tar_target(file_bac_missing, "raw_data/new_busco_results_bac.txt", format = "file"),
  tar_target(file_raw_bac_mapping_pan, "raw_data/gtdb_norRNA_all_stat.parsed.txt", format = "file"),
  tar_target(file_raw_fun_mapping_pan, "raw_data/fungi_norRNA_all_stat.parsed.txt", format = "file"),
  tar_target(file_lane_info, "raw_data/lane_info.txt", format = "file"),
  tar_target(file_norRNA_reads, "raw_data/norRNA_reads.txt", format = "file"),
  tar_target(file_pathogen_manifest, "raw_data/KEGG_pathogen_lineages.txt", format = "file"),
  tar_target(file_virus_metadata, "raw_data/virus_metadata.txt", format = "file"),
  tar_target(file_tax_map, "raw_data/taxonomy_map.txt.gz", format = "file"),
  
  # lefse results
  tar_target(file_lefse_res_anal_cat, "output/lefse_input_anal_cat.res", format = "file"),
  tar_target(file_lefse_res_anal_dog, "output/lefse_input_anal_dog.res", format = "file"),
  tar_target(file_lefse_res_throat_cat, "output/lefse_input_throat_cat.res", format = "file"),
  tar_target(file_lefse_res_throat_dog, "output/lefse_input_throat_dog.res", format = "file"),
  
  # load data
  tar_target(lib_meta_raw, read_csv(file_lib_meta_raw)),
  tar_target(fastq_manifest, readLines(file_fastq_manifest)),
  tar_target(species_info, read_tsv(file_species_info)),
  tar_target(raw_vir_profile_amy, read_excel(file_raw_vir_profile_amy, skip=1)),
  tar_target(raw_bac_profile_amy, read_excel(file_raw_bac_profile_amy)),
  #tar_target(raw_bac_profile_genus_pan, read_csv(file_raw_bac_profile_genus_pan)),
  #tar_target(raw_fun_profile_genus_pan, read_csv(file_raw_fun_profile_genus_pan)),
  tar_target(missing_markers, load_busco_results(file_bac_missing, file_fun_missing)),
  tar_target(lane_info, read_tsv(file_lane_info)),
  tar_target(norRNA_reads, read_tsv(file_norRNA_reads)),
  tar_target(raw_map_table_gtdb, read_and_filter_raw_mapping(file_raw_bac_mapping_pan)),
  tar_target(raw_map_table_fungi, read_and_filter_raw_mapping(file_raw_fun_mapping_pan)),
  tar_target(pathogen_manifest, read_tsv("raw_data/KEGG_pathogen_lineages.txt",
                                         col_names = c("taxid", "kingdom", "phylum", "class", 
                                                       "order", "family", "genus", "species"))),
  tar_target(virus_metadata, read_tsv(file_virus_metadata)),
  tar_target(tax_map, read_tsv(file_tax_map)),
  
  # processing raw data
  tar_target(lib_meta, gen_clean_lib_meta(lib_meta_raw, species_info)),
  tar_target(out_lib_stat, write_lib_stat(lib_meta), format = "file"),
  tar_target(tbl_validation, verify_libraries(lib_meta_raw, fastq_manifest, species_info)),
  tar_target(raw_bac_profile_genus_pan, calc_rpm(raw_map_table_gtdb, lane_info, norRNA_reads, RPM_threshold=100)),
  tar_target(raw_fun_profile_genus_pan, calc_rpm(raw_map_table_fungi, lane_info, norRNA_reads, RPM_threshold=100)),
  tar_target(rpm_table_vir_amy, gen_rpm_table_vir_amy(lib_meta, raw_vir_profile_amy)),
  tar_target(rpm_table_bac_amy, gen_rpm_table_bac_amy(lib_meta, raw_bac_profile_amy)),
  tar_target(rpm_table_bac_genus_pan, gen_rpm_table_pan(lib_meta, raw_bac_profile_genus_pan, missing_markers)),
  tar_target(rpm_table_fun_genus_pan, gen_rpm_table_pan(lib_meta, raw_fun_profile_genus_pan, missing_markers)),
  
  tar_target(data_pack_pan, list(lib_meta = lib_meta,
                                 rpm_table_vir = rpm_table_vir_amy,
                                 rpm_table_bac = rpm_table_bac_genus_pan,
                                 rpm_table_fun = rpm_table_fun_genus_pan)),
  tar_target(lefse_results, load_lefse_results(c("anal_cat" = file_lefse_res_anal_cat,
                                                 "anal_dog" = file_lefse_res_anal_dog,
                                                 "throat_cat" = file_lefse_res_throat_cat,
                                                 "throat_dog" = file_lefse_res_throat_dog))),
  
  # analyzing and plotting
  tar_target(out_vir_heatmap_amy, draw_vir_heatmap(lib_meta, rpm_table_vir_amy), format = "file"),
  tar_target(out_bac_heatmap_amy, draw_bac_heatmap(lib_meta, rpm_table_bac_amy, "genus_amy"), format = "file"),
  tar_target(out_bac_heatmap_pan, draw_bac_heatmap(lib_meta, rpm_table_bac_genus_pan, "genus_pan"), format = "file"),
  tar_target(out_fun_heatmap_pan, draw_bac_heatmap(lib_meta, rpm_table_fun_genus_pan, "fun_genus_pan"), format = "file"),
  
  tar_target(out_alpha_div_plots_amy, draw_alpha_div_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_alpha_div_stats_amy, stat_alpha_div_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_beta_div_plots_amy, draw_beta_div_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_beta_div_stats_amy, stat_beta_div_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_venn_diagrams_amy, draw_venn_diagrams_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_deseq_analyses_amy, deseq_analysis_amy(lib_meta, rpm_table_vir_amy, rpm_table_bac_amy), format = "file"),
  tar_target(out_infectome_heatmap, draw_infectome_heatmap(data_pack_pan, pathogen_manifest), format = "file"),
  tar_target(out_beta_div_pcoa, draw_beta_div_pcoa(data_pack_pan), format = "file"),
  tar_target(out_alpha_div, calc_and_comp_alpha_div(data_pack_pan), format = "file"),
  tar_target(out_alpha_div_plots, draw_alpha_div(data_pack_pan), format = "file"),
  tar_target(out_beta_div_tests, test_beta_div(data_pack_pan), format = "file"),
  tar_target(out_taxa_compositions, draw_taxa_composition(data_pack_pan, virus_metadata, tax_map), format = "file"),
  tar_target(out_lefse_inputs, prepare_lefse_data(data_pack_pan, virus_metadata, tax_map), format = "file"),
  #tar_target(out_deseq_analyses_pan, deseq_analysis_pan(lib_meta, rpm_table_bac_genus_pan, rpm_table_fun_genus_pan), format = "file"),

  # deseq analysis
  tar_target(rpm_table_all, merge_rpm_tables(rpm_table_vir_amy, rpm_table_bac_genus_pan, rpm_table_fun_genus_pan)),
  tar_target(out_diff_all_cat_anal, diff_expr(lib_meta, rpm_table_all, "Cat", "Anal swab", "all"), format = "file"),
  tar_target(out_diff_all_cat_throat, diff_expr(lib_meta, rpm_table_all, "Cat", "Throat swab", "all"), format = "file"),
  tar_target(out_diff_all_dog_anal, diff_expr(lib_meta, rpm_table_all, "Dog", "Anal swab", "all"), format = "file"),
  tar_target(out_diff_all_dog_throat, diff_expr(lib_meta, rpm_table_all, "Dog", "Throat swab", "all"), format = "file"),
  
  tar_target(out_diff_lefse_cat_anal, draw_diff_heatmap_lefse(rpm_table_all, lib_meta, lefse_results, "Cat", "Anal swab"), format = "file"),
  tar_target(out_diff_lefse_cat_throat, draw_diff_heatmap_lefse(rpm_table_all, lib_meta, lefse_results, "Cat", "Throat swab"), format = "file"),
  tar_target(out_diff_lefse_dog_anal, draw_diff_heatmap_lefse(rpm_table_all, lib_meta, lefse_results, "Dog", "Anal swab"), format = "file"),
  tar_target(out_diff_lefse_dog_throat, draw_diff_heatmap_lefse(rpm_table_all, lib_meta, lefse_results, "Dog", "Throat swab"), format = "file"),
  
  tar_target(out_diff_lefse_deseq_cat_anal, draw_diff_heatmap_deseq2_and_lefse(rpm_table_all, lib_meta, lefse_results, "Cat", "Anal swab"), format = "file"),
  tar_target(out_diff_lefse_deseq_cat_throat, draw_diff_heatmap_deseq2_and_lefse(rpm_table_all, lib_meta, lefse_results, "Cat", "Throat swab"), format = "file"),
  tar_target(out_diff_lefse_deseq_dog_anal, draw_diff_heatmap_deseq2_and_lefse(rpm_table_all, lib_meta, lefse_results, "Dog", "Anal swab"), format = "file"),
  tar_target(out_diff_lefse_deseq_dog_throat, draw_diff_heatmap_deseq2_and_lefse(rpm_table_all, lib_meta, lefse_results, "Dog", "Throat swab"), format = "file"),
  
  # manual DAA
  tar_target(
    out_daa_manual, 
    diffential_abundance_test(data_pack_pan), 
    format = "file",
    packages = c("tidyverse", "ggpubr", "cowplot")
  ),
  
  # 3methods differential abundance analysis heatmap
  tar_target(
    out_heatmap_daa_3methods_cat_anal, 
    draw_diff_heatmap_deseq2_and_lefse_and_manual_daa(rpm_table_all, lib_meta, lefse_results, "Cat", "Anal swab"), 
    format = "file",
    packages = c("tidyverse", "ComplexHeatmap", "DESeq2")
  ),
  tar_target(
    out_heatmap_daa_3methods_cat_throat, 
    draw_diff_heatmap_deseq2_and_lefse_and_manual_daa(rpm_table_all, lib_meta, lefse_results, "Cat", "Throat swab"), 
    format = "file",
    packages = c("tidyverse", "ComplexHeatmap", "DESeq2")
  ),
  tar_target(
    out_heatmap_daa_3methods_dog_anal, 
    draw_diff_heatmap_deseq2_and_lefse_and_manual_daa(rpm_table_all, lib_meta, lefse_results, "Dog", "Anal swab"), 
    format = "file",
    packages = c("tidyverse", "ComplexHeatmap", "DESeq2")
  ),
  tar_target(
    out_heatmap_daa_3methods_dog_throat, 
    draw_diff_heatmap_deseq2_and_lefse_and_manual_daa(rpm_table_all, lib_meta, lefse_results, "Dog", "Throat swab"), 
    format = "file",
    packages = c("tidyverse", "ComplexHeatmap", "DESeq2")
  ),
  
  # writing infectome profiles to files
  tar_target(out_vir_profile, write_microb_profile(lib_meta, rpm_table_vir_amy, "vir_species_rpm1.csv"), format = "file"),
  tar_target(out_bac_profile_genus_amy, write_microb_profile(lib_meta, rpm_table_bac_amy, "bac_genus_metaphlan.csv"), format = "file"),
  tar_target(out_bac_profile_genus_pan, write_microb_profile(lib_meta, rpm_table_bac_genus_pan, "bac_genus_gtdb_rpm100.csv"), format = "file"),
  tar_target(out_fun_profile_genus_pan, write_microb_profile(lib_meta, rpm_table_fun_genus_pan, "fun_genus_gtdb_rpm100.csv"), format = "file")
)
