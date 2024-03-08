verify_libraries = function(lib_meta_raw, fastq_manifest, species_info) {
  bind_rows(
    tibble(
      lib_id = setdiff(fastq_manifest, lib_meta_raw$`Libraries-sever`),
      reason = "Libs not analysed."
    ),
    tibble(
      lib_id = setdiff(lib_meta_raw$`Libraries-sever`, species_info$lib_id),
      reason = "Libs not identified."
    )
  )
}

gen_clean_lib_meta = function(lib_meta_raw, species_info) {
  
  libs_to_drop = c(
    # contaminated
    "CHEBT527", "CXXT191", "DGDA026", 
    "DHEBA439-0221", "DHEBA491-0310", "DNBT423",
    "DGDA545",
    
    # duplicated
    "DXXA299-0221", "DXXA294-0221", "DNBT419"
  )
  
  lib_meta_raw %>%
    left_join(species_info, by = c("Libraries-sever"="lib_id")) %>%
    mutate(Species = str_to_sentence(coi_species)) %>%
    filter(Species != "Mix") %>%
    filter(!Libraries %in% libs_to_drop)
}

read_and_filter_raw_mapping = function(filename) {
  read_tsv(filename) %>% 
    filter(mapped_len > 300) %>%
    group_by(lib_id, genus) %>%
    summarise(numreads = sum(mapped_reads)) %>%
    ungroup() %>%
    dplyr::rename(virus_name=genus)
}

calc_rpm = function(read_mapping_table, 
                    lane_id_table,
                    read_count_table, 
                    index_hopping_threshold=0.001, 
                    RPM_threshold=1) {
  
  read_table_wide = read_mapping_table %>%
    select(lib_id, virus_name, numreads) %>%
    pivot_wider(names_from="virus_name", values_from="numreads", values_fill=0)
  
  lib_id_arr = read_table_wide$lib_id
  virus_name_arr = colnames(read_table_wide)[2:ncol(read_table_wide)]
  
  # Merge lane id info into the transposed table
  reads_mapping_t_lane = inner_join(lane_id_table, read_table_wide, by="lib_id")
  # check if the two tables are correctly joined
  flag = setequal(reads_mapping_t_lane$lib_id, lib_id_arr) 
  if (!flag) { warning("join failed.") }
  
  reads_sum_by_lane = reads_mapping_t_lane %>%
    group_by(lane_id) %>%
    summarise(across(all_of(virus_name_arr), sum))
  
  # set threshold to reduce false positive due to index hopping
  threshold_matrix = reads_mapping_t_lane %>%
    select(lib_id, lane_id) %>%
    left_join(reads_sum_by_lane, by="lane_id") %>%
    select(all_of(virus_name_arr)) %>%
    mutate_all(function(x) {x*index_hopping_threshold}) # threshold value set at 0.1%
  
  raw_reads_matrix = select(reads_mapping_t_lane, all_of(virus_name_arr))
  mask_matrix = raw_reads_matrix > threshold_matrix
  
  filtered_count_lane = sum((raw_reads_matrix <= threshold_matrix) & (raw_reads_matrix != 0))
  print(str_glue("{filtered_count_lane} possible false positives deleted (index hopping filter)"))
  
  filtered_reads_matrix = mask_matrix * raw_reads_matrix
  
  reads_mapping_t_lane_filtered = bind_cols(
    select(reads_mapping_t_lane, lib_id, lane_id), 
    filtered_reads_matrix)
  
  # load norRNA reads count
  norRNA_reads = read_count_table %>% select(lib_id, norRNA_reads)
  lane_filtered_with_reads = inner_join(norRNA_reads, reads_mapping_t_lane_filtered, by="lib_id")
  # check if the two tables are correctly joined
  flag = assertthat::are_equal(nrow(lane_filtered_with_reads), nrow(reads_mapping_t_lane_filtered)) 
  if (!flag) { warning("join failed.") }
  
  # calc RPM and apply RPM > 1 filter
  virus_table_filtered = lane_filtered_with_reads %>%
    mutate(across(all_of(virus_name_arr), function(x) {x/norRNA_reads * 1e6})) %>%
    mutate(across(all_of(virus_name_arr), function(x) {ifelse(x>1, x, 0)})) %>%
    select(-norRNA_reads, -lane_id) %>%
    select(where(~ifelse(is.numeric(.x), sum(.x)>0, T)))
  
  RPM_mat = lane_filtered_with_reads %>% 
    mutate(across(all_of(virus_name_arr), function(x) {x/norRNA_reads * 1e6})) %>% 
    select(all_of(virus_name_arr))
  filtered_count_RPM = sum(RPM_mat >0 & RPM_mat <= 1)
  print(str_glue("{filtered_count_RPM} possible false positives deleted (RPM filter)"))
  
  return(virus_table_filtered)
}

gen_rpm_table_vir_amy = function(lib_meta_env, raw_vir_profile) {
  # select Healthy and Diseased subset
  raw_vir_profile = raw_vir_profile %>%
    filter(`Health condition` %in% c("Healthy", "Diseased"))
  # do index-hopping correction: virus
  ih_threshold_matrix = raw_vir_profile %>%
    group_by(lane) %>%
    summarise(across(`Alphacoronavirus 1`:`Feline Foamy Virus`, ~sum(.x)*0.001)) %>%
    left_join(select(raw_vir_profile, lane), ., by="lane")
  
  # confirm that rows matches with original table
  if (any(ih_threshold_matrix$lane != raw_vir_profile$lane)) stop()
  
  lib_meta = select(raw_vir_profile, `Correction`:`Libraries-sever`)
  rpm_table = select(raw_vir_profile, `Alphacoronavirus 1`:`Feline Foamy Virus`)
  
  mask_matrix_ih = rpm_table > select(ih_threshold_matrix, -lane)
  mask_matrix_rpm1 = rpm_table > 10
  
  rpm_table_filtered = rpm_table * mask_matrix_ih * mask_matrix_rpm1
  rpm_table_final = bind_cols(Libraries=lib_meta$Libraries, rpm_table_filtered)
  
  as_tibble(rpm_table_final) %>%
    filter(Libraries %in% lib_meta_env$Libraries)
}

gen_rpm_table_bac_amy = function(lib_meta_env, raw_bac_profile) {
  # column name correction, change to the same names as in virus profile 
  raw_bac_profile = raw_bac_profile %>%
    dplyr::rename(Correction=Notes, `Libraries-sever`=`Libraries--sever`) %>%
    filter(`Health condition` %in% c("Healthy", "Diseased"))
  
  # Correct library name DXXT277 -> DXXT227
  row_id = which(raw_bac_profile$Libraries=="DXXT277")
  raw_bac_profile[row_id, "Libraries"] = "DXXT227"
  
  # fix missing values (treat as zero)
  bac_colnames = raw_bac_profile %>%
    select(`g__Aerococcus`:`g__Xanthomonadaceae_unclassified`) %>%
    colnames()
  narep_list = setNames(as.list(rep(0, length(bac_colnames))), bac_colnames)
  raw_bac_profile = raw_bac_profile %>%
    replace_na(narep_list)
  
  ih_threshold_matrix = raw_bac_profile %>%
    group_by(lane) %>%
    summarise(across(`g__Aerococcus`:`g__Xanthomonadaceae_unclassified`, ~sum(.x)*0.001)) %>%
    left_join(select(raw_bac_profile, lane), ., by="lane")
  
  # confirm that rows matches with original table
  if (any(ih_threshold_matrix$lane != raw_bac_profile$lane)) stop()
  
  lib_meta = select(raw_bac_profile, `Correction`:`Libraries-sever`)
  rpm_table = select(raw_bac_profile, `g__Aerococcus`:`g__Xanthomonadaceae_unclassified`)
  
  mask_matrix_ih = rpm_table > select(ih_threshold_matrix, -lane)
  mask_matrix_rpm1 = rpm_table > 100
  
  rpm_table_filtered = rpm_table * mask_matrix_ih * mask_matrix_rpm1
  rpm_table_final = bind_cols(Libraries=lib_meta$Libraries, rpm_table_filtered)
  
  as_tibble(rpm_table_final) %>%
    filter(Libraries %in% lib_meta_env$Libraries)
}

draw_vir_heatmap = function(lib_meta, vir_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  #### Virus ####
  # log1p transformation
  logRPM_virus = as.data.frame(vir_rpm_table)
  rownames(logRPM_virus) = logRPM_virus$Libraries
  logRPM_virus$Libraries = NULL
  logRPM_virus = log10(logRPM_virus+1)
  
  # re-ordering rows (samples or libraries)
  logRPM_virus = logRPM_virus[lib_ordered$Libraries,]
  
  # transpose the table (viruses are rows)
  logRPM_virus = t(logRPM_virus)
  
  # plot settings
  # must manually set factor levels to match the row ordering in lib_ordered
  col_split = paste0(lib_ordered$Species, lib_ordered$`Health condition`) %>%
    factor(levels=c("CatHealthy", "CatDiseased", "DogHealthy", "DogDiseased"))
  
  # RPM color mapper
  RPM_col_fun = colorRamp2(c(0, ceiling(max(logRPM_virus))), 
                           c("#F4F2ED", "#DD1332"))
  
  # location strip
  locations = unique(lib_ordered$Location)
  color_loc = brewer.pal(length(locations), "Set2")
  legend_loc = Legend(labels = locations, 
                      title="Location", 
                      legend_gp=gpar(fill=color_loc))
  
  # sample strip
  sample_type = unique(lib_ordered$Type)
  color_type = brewer.pal(length(sample_type), "Set1")
  legend_type = Legend(labels = sample_type, 
                       title="Sample type", 
                       legend_gp=gpar(fill=color_type))
  
  # virus family strip
  viral_fam = c("RNA", "DNA", "Retro")
  color_vfam = c(brewer.pal(3, "Set3"))
  legend_type = Legend(labels = viral_fam, 
                       title="Viral group", 
                       legend_gp=gpar(fill=color_vfam))
  
  col_anno = columnAnnotation(
    location = anno_simple(lib_ordered$Location, 
                           col=setNames(color_loc, locations),
                           height=unit(0.3, "cm")),
    sample_type = anno_simple(lib_ordered$Type, 
                              col=setNames(color_type, sample_type),
                              height=unit(0.3, "cm")),
    gap = unit(1, "mm")
  )
  
  row_anno = rowAnnotation(
    virus_group = anno_simple(c(
      rep("RNA", 17),
      rep("DNA", 6),
      rep("Retro", 1)
    ), 
    col=setNames(color_vfam, viral_fam),
    width=unit(0.3, "cm"))
  )
  
  output_filename = "output/virus_heatmap.pdf"
  # draw heatmap: virus
  pdf(file=output_filename, width=14, height=5)
  ht = Heatmap(logRPM_virus,
               cluster_columns = F, 
               cluster_rows = F,
               col = RPM_col_fun,
               column_split = col_split,
               top_annotation = col_anno,
               left_annotation = row_anno,
               heatmap_legend_param = list(title = "logRPM+1"),
               column_labels = rep("", ncol(logRPM_virus)))
  draw(ht, annotation_legend_list = list(legend_loc, legend_type))
  dev.off()
  
  output_filename
}

draw_bac_heatmap = function(lib_meta, bac_rpm_table, output_suffix) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  logRPM_bac = bac_rpm_table %>%
    select(!(starts_with("g__GGB") | contains("unclassified"))) %>%
    as.data.frame()
  rownames(logRPM_bac) = logRPM_bac$Libraries
  logRPM_bac$Libraries = NULL
  
  # ranking bacteria genera according to relative abundance
  cc = sort(colSums(logRPM_bac)) %>% cumsum()
  cc = cc / max(cc)
  
  ggplot(aes(x=1:length(cc), y=cc), data=NULL) + 
    geom_bar(stat="identity") +
    geom_hline(yintercept=0.01, linetype=2) +
    xlab("Bacteria geneus rank") +
    ylab("Cumulative relative abundance") +
    theme_bw()
  bac_rank = str_glue("output/bacteria_genera_rank_{output_suffix}.pdf")
  ggsave(bac_rank, width=4, height=3)
  
  # top-ranked bacteria genera which account for 95% of the total bacterial reads
  bac_rank_text = str_glue("output/bacteria_genera_rank99_{output_suffix}.txt")
  rank99 = names(cc)[cc>0]
  write_lines(rank99, bac_rank_text)
  
  logRPM_bac = select(logRPM_bac, all_of(rev(rank99)))
  # log1p transformation
  logRPM_bac = log10(logRPM_bac+1)
  
  # re-ordering rows (samples or libraries)
  logRPM_bac = logRPM_bac[lib_ordered$Libraries,]
  
  rownames(logRPM_bac) = lib_ordered$Libraries
  logRPM_bac = apply(logRPM_bac, 2, function(x) {x[is.na(x)] <- 0; x})
  
  # transpose the table (viruses are rows)
  logRPM_bac = t(logRPM_bac)
  
  # plot settings
  # must manually set factor levels to match the row ordering in lib_ordered
  col_split = paste0(lib_ordered$Species, lib_ordered$`Health condition`) %>%
    factor(levels=c("CatHealthy", "CatDiseased", "DogHealthy", "DogDiseased"))
  
  # RPM color mapper
  RPM_col_fun = colorRamp2(c(0, ceiling(max(logRPM_bac))), 
                           c("#F4F2ED", "#DD1332"))
  
  # location strip
  locations = unique(lib_ordered$Location)
  color_loc = brewer.pal(length(locations), "Set2")
  legend_loc = Legend(labels = locations, 
                      title="Location", 
                      legend_gp=gpar(fill=color_loc))
  
  # sample strip
  sample_type = unique(lib_ordered$Type)
  color_type = brewer.pal(length(sample_type), "Set1")
  legend_type = Legend(labels = sample_type, 
                       title="Sample type", 
                       legend_gp=gpar(fill=color_type))
  
  col_anno = columnAnnotation(
    location = anno_simple(lib_ordered$Location, 
                           col=setNames(color_loc, locations),
                           height=unit(0.3, "cm")),
    sample_type = anno_simple(lib_ordered$Type, 
                              col=setNames(color_type, sample_type),
                              height=unit(0.3, "cm")),
    gap = unit(1, "mm")
  )
  
  # draw heatmap: bacteria
  bac_heatmap = str_glue("output/bacteria_heatmap_{output_suffix}.pdf")
  
  pdf(file=bac_heatmap, width=14, height=0.1785714 * nrow(logRPM_bac))
  ht = Heatmap(logRPM_bac,
               cluster_columns = F, 
               cluster_rows = F,
               col = RPM_col_fun,
               column_split = col_split,
               top_annotation = col_anno,
               heatmap_legend_param = list(title = "logRPM+1"),
               column_labels = rep("", ncol(logRPM_bac)),
               row_labels = str_replace(rownames(logRPM_bac), "g__", ""))
  draw(ht, annotation_legend_list = list(legend_loc, legend_type))
  dev.off()
  
  c(bac_rank, bac_rank_text, bac_heatmap)
}

draw_alpha_div_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  # calculate viral and bacterial load
  vload = tibble(
    Libraries = vir_rpm_table$Libraries, 
    viral_load = rowSums(select(vir_rpm_table, -Libraries)),
    num_viral_species = rowSums(select(vir_rpm_table, -Libraries)>0))
  bload = tibble(
    Libraries = bac_rpm_table$Libraries, 
    bacterial_load = rowSums(select(bac_rpm_table, -Libraries)),
    num_bac_genera = rowSums(select(bac_rpm_table, -Libraries)>0))
  
  load_tb = full_join(vload, bload, by="Libraries")
  
  # re-ordering
  load_ordered = left_join(lib_ordered, load_tb, by="Libraries") %>%
    mutate(Group=paste(Species, `Health condition`, Type))
  
  # draw bar plots
  g1 = ggplot(aes(x=1:nrow(load_ordered), 
                  y=log10(viral_load+1), 
                  fill=Group), 
              data=load_ordered) + 
    geom_bar(stat="identity") +
    ylab("Viral load\n(logRPM+1)") + 
    xlab("Sample") +
    theme(legend.position = "none")
  
  g2 = ggplot(aes(x=1:nrow(load_ordered), 
                  y=log10(bacterial_load+1), 
                  fill=Group), 
              data=load_ordered) + 
    geom_bar(stat="identity") +
    ylab("Bacterial load\n(logRPM+1)") + 
    xlab("Sample") +
    theme(legend.position = "none")
  
  g3 = ggplot(aes(x=1:nrow(load_ordered), 
                  y=num_viral_species, 
                  fill=Group), 
              data=load_ordered) + 
    geom_bar(stat="identity") +
    ylab("Num. of\nvirus species") + 
    xlab("Sample") +
    theme(legend.position = "none")
  
  g4 = ggplot(aes(x=1:nrow(load_ordered), 
                  y=num_bac_genera, 
                  fill=Group), 
              data=load_ordered) + 
    geom_bar(stat="identity") +
    ylab("Num. of\nbac. genera") + 
    xlab("Sample") +
    theme(legend.position = "none")
  
  plot_grid(g1,g2,g3,g4, nrow=4, align="h")
  ggsave("output/alpha_diversity.pdf", width=6, height=7)
  
  "output/alpha_diversity.pdf"
}

stat_alpha_div_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  # calculate viral and bacterial load
  vload = tibble(
    Libraries = vir_rpm_table$Libraries, 
    viral_load = rowSums(select(vir_rpm_table, -Libraries)),
    num_viral_species = rowSums(select(vir_rpm_table, -Libraries)>0))
  bload = tibble(
    Libraries = bac_rpm_table$Libraries, 
    bacterial_load = rowSums(select(bac_rpm_table, -Libraries)),
    num_bac_genera = rowSums(select(bac_rpm_table, -Libraries)>0))
  
  load_tb = full_join(vload, bload, by="Libraries")
  
  # re-ordering
  load_ordered = left_join(lib_ordered, load_tb, by="Libraries") %>%
    mutate(Group=paste(Species, `Health condition`, Type))
  
  model_viral_load = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp)
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    f_fe = gls(log10(viral_load+1) ~ HealthCondition, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    f_me = lme(log10(viral_load+1) ~ HealthCondition, 
               random=~1|Location, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    # check the model comparison result manually before going ahead!
    model_comp = anova(f_fe, f_me)
    
    aov_fe = anova(f_fe)
    aov_me = anova(f_me)
    
    return(list(model_comparison=model_comp, fixed_effect=aov_fe, mixed_effect=aov_me))
  }
  model_viral_richness = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp)
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    f_fe = gls(num_viral_species ~ HealthCondition, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    f_me = lme(num_viral_species ~ HealthCondition, 
               random=~1|Location, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    # check the model comparison result manually before going ahead!
    model_comp = anova(f_fe, f_me)
    
    aov_fe = anova(f_fe)
    aov_me = anova(f_me)
    
    return(list(model_comparison=model_comp, fixed_effect=aov_fe, mixed_effect=aov_me))
  }
  model_bacterial_load = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp)
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    f_fe = gls(log10(bacterial_load+1) ~ HealthCondition, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    f_me = lme(log10(bacterial_load+1) ~ HealthCondition, 
               random=~1|Location, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    # check the model comparison result manually before going ahead!
    model_comp = anova(f_fe, f_me)
    
    aov_fe = anova(f_fe)
    aov_me = anova(f_me)
    
    return(list(model_comparison=model_comp, fixed_effect=aov_fe, mixed_effect=aov_me))
  }
  model_bacterial_richness = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp)
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    f_fe = gls(num_bac_genera ~ HealthCondition, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    f_me = lme(num_bac_genera ~ HealthCondition, 
               random=~1|Location, 
               weights=varIdent(~1|HealthCondition), 
               data=tmp_data, 
               method="REML")
    
    # check the model comparison result manually before going ahead!
    model_comp = anova(f_fe, f_me)
    
    aov_fe = anova(f_fe)
    aov_me = anova(f_me)
    
    return(list(model_comparison=model_comp, fixed_effect=aov_fe, mixed_effect=aov_me))
  }
  
  candidates = expand.grid(species = c("Dog", "Cat"), sample = c("Anal swab", "Throat swab"))
  
  output_filenames = NULL
  for (i in 1:nrow(candidates)) {
    species = candidates$species[i]
    sample = candidates$sample[i]
    
    fname = str_glue("output/HC_effect_VL_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    model_viral_load(species, sample) %>% 
      capture.output(file=fname)
    
    fname = str_glue("output/HC_effect_VR_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    model_viral_richness(species, sample) %>% 
      capture.output(file=fname)
    
    fname = str_glue("output/HC_effect_BL_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    model_bacterial_load(species, sample) %>% 
      capture.output(file=fname)
    
    fname = str_glue("output/HC_effect_BR_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    model_bacterial_richness(species, sample) %>% 
      capture.output(file=fname)
  }
  
  output_filenames
}

draw_beta_div_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  # virus
  vir_mat = select(vir_rpm_table, -Libraries) %>% as.data.frame()
  rownames(vir_mat) = vir_rpm_table$Libraries
  vir_mat = decostand(vir_mat[rowSums(vir_mat)>0,], "hellinger")
  
  xx = left_join(vir_rpm_table, lib_meta) %>%
    group_by(Species, `Health condition`, Type, Location) %>%
    summarise(across(`Alphacoronavirus 1`:`Feline Foamy Virus`, mean))
  
  mm = xx[,1:4]
  xx[,5:ncol(xx)]
  
  ca = vegan::capscale(log10(xx[,5:ncol(xx)]+1)~1, distance="jaccard", binary=T)
  
  mm$PCoA1 = ca$CA$u[,1]
  mm$PCoA2 = ca$CA$u[,2]
  #qplot(x=PCoA1, y=PCoA2, color=`Health condition`, data=mm) + 
  #  facet_wrap(~Type+Species)
  
  
  xx = left_join(bac_rpm_table, lib_meta) %>%
    group_by(Species, `Health condition`, Type, Location) %>%
    summarise(across(`g__Aerococcus`:`g__Xanthomonadaceae_unclassified`, mean))
  
  mm = xx[,1:4]
  ca = vegan::capscale(log10(xx[,5:ncol(xx)]+1)~1, distance="jaccard")
  
  mm$PCoA1 = ca$CA$u[,1]
  mm$PCoA2 = ca$CA$u[,2]
  #qplot(x=PCoA1, y=PCoA2, color=`Health condition`, data=mm) + 
  #  facet_wrap(~Type+Species)
  
  meta_df = as.data.frame(lib_ordered)
  rownames(meta_df) = lib_ordered$Libraries
  meta_df = meta_df[rownames(vir_mat), ]
  vir_cca = rda(vir_mat)
  meta_df$PC1 = vir_cca$CA$u[,1]
  meta_df$PC2 = vir_cca$CA$u[,2]
  g1 = ggplot(data=meta_df) +
    geom_point(aes(x=PC1, 
                   y=PC2, 
                   color=paste(Species, `Health condition`), 
                   shape=Type), 
               size = 3) +
    ggtitle("Virus") +
    theme(legend.position = "noen")
  
  # bacteria
  bac_mat = select(bac_rpm_table, -Libraries) %>% as.data.frame()
  rownames(bac_mat) = bac_rpm_table$Libraries
  #bac_mat = decostand(bac_mat[rowSums(bac_mat)>0,], "hellinger")
  bac_mat = bac_mat[rowSums(bac_mat)>0,]
  meta_df = as.data.frame(lib_ordered)
  rownames(meta_df) = lib_ordered$Libraries
  meta_df = meta_df[rownames(bac_mat), ]
  bac_cca = capscale(log10(bac_mat+1)~1, distance="bray")
  meta_df$PC1 = bac_cca$CA$u[,1]
  meta_df$PC2 = bac_cca$CA$u[,2]
  g2 = ggplot(data=meta_df) +
    geom_point(aes(x=PC1, 
                   y=PC2, 
                   color=paste(Species, `Health condition`), 
                   shape=Type), 
               size = 3) +
    facet_wrap(~Type) +
    stat_ellipse(aes(x=PC1, 
                     y=PC2, 
                     color=paste(Species, `Health condition`))) + 
    ggtitle("Bacteria")
  
  plot_grid(g1, g2, ncol=2, align="h", rel_widths = c(1,1.55))
  ggsave("output/beta_diversity.pdf", width=12, height=5)
  
  "output/beta_diversity.pdf"
}

stat_beta_div_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  # calculate viral and bacterial load
  vload = tibble(
    Libraries = vir_rpm_table$Libraries, 
    viral_load = rowSums(select(vir_rpm_table, -Libraries)),
    num_viral_species = rowSums(select(vir_rpm_table, -Libraries)>0))
  bload = tibble(
    Libraries = bac_rpm_table$Libraries, 
    bacterial_load = rowSums(select(bac_rpm_table, -Libraries)),
    num_bac_genera = rowSums(select(bac_rpm_table, -Libraries)>0))
  
  load_tb = full_join(vload, bload, by="Libraries")
  
  # re-ordering
  load_ordered = left_join(lib_ordered, load_tb, by="Libraries") %>%
    mutate(Group=paste(Species, `Health condition`, Type))
  
  adonis_virus = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp) %>% as.data.frame()
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    rownames(tmp_data) = tmp_data$Libraries
    
    vir_mat = vir_rpm_table %>% as.data.frame()
    rownames(vir_mat) = vir_rpm_table$Libraries
    vir_mat$Libraries = NULL
    
    vir_mat = vir_mat[tmp_data$Libraries,]
    vir_mat = vir_mat[rowSums(vir_mat)>0,]
    tmp_data = tmp_data[rownames(vir_mat),]
    
    adonis2(vir_mat ~ HealthCondition, data=tmp_data)
  }
  adonis_bacteria = function(sp, tp) {
    tmp_data = filter(load_ordered, Species == sp, Type==tp) %>% as.data.frame()
    tmp_data$HealthCondition = tmp_data$`Health condition`
    
    rownames(tmp_data) = tmp_data$Libraries
    
    vir_mat = bac_rpm_table %>% as.data.frame()
    rownames(vir_mat) = vir_rpm_table$Libraries
    vir_mat$Libraries = NULL
    
    vir_mat = vir_mat[tmp_data$Libraries,]
    vir_mat = vir_mat[rowSums(vir_mat)>0,]
    tmp_data = tmp_data[rownames(vir_mat),]
    
    adonis2(vir_mat ~ HealthCondition, data=tmp_data)
  }
  
  candidates = expand.grid(species = c("Dog", "Cat"), sample = c("Anal swab", "Throat swab"))
  
  output_filenames = NULL
  for (i in 1:nrow(candidates)) {
    species = candidates$species[i]
    sample = candidates$sample[i]
    
    fname = str_glue("output/HC_effect_VA_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    adonis_virus(species, sample) %>% 
      capture.output(file=fname)
  
    fname = str_glue("output/HC_effect_BA_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    adonis_bacteria(species, sample) %>% 
      capture.output(file=fname)
  }
  
  output_filenames
}

draw_venn_diagrams_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  # virus
  vir_names = colnames(vir_rpm_table)[2:ncol(vir_rpm_table)]
  vir_venn_data = left_join(lib_ordered, vir_rpm_table, by="Libraries") %>%
    group_by(Species, `Health condition`, Type) %>%
    summarise(across(all_of(vir_names), ~sum(.x)>0))
  vir_venn_data_notp = left_join(lib_ordered, vir_rpm_table, by="Libraries") %>%
    group_by(Species, `Health condition`) %>%
    summarise(across(all_of(vir_names), ~sum(.x)>0))
  
  dc = vir_venn_data_notp[1,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  hc = vir_venn_data_notp[2,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  
  diseased_cat_only_vir = "output/diseased_cat_only_vir.txt"
  vir_names[setdiff(dc, hc)] %>% write_lines(diseased_cat_only_vir)
  
  dd = vir_venn_data_notp[3,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  hd = vir_venn_data_notp[4,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  
  diseased_dog_only_vir = "output/diseased_dog_only_vir.txt"
  vir_names[setdiff(dd, hd)] %>% write_lines(diseased_dog_only_vir)
  
  draw_venn_sptp = function(sp, tp) {
    tmp = filter(vir_venn_data, Species==sp, Type==tp) %>% 
      select(-Species, -Type) %>% 
      t()
    colnames(tmp) = tmp[2,]
    tmp = as.data.frame(tmp[3:nrow(tmp),])
    tmp[,1] = ifelse(tmp[,1]=="TRUE", T, F)
    tmp[,2] = ifelse(tmp[,2]=="TRUE", T, F)
    
    ggvenn(tmp, c("Diseased", "Healthy")) + ggtitle(str_glue(sp, tp))
  }
  draw_venn_sp = function(sp) {
    tmp = filter(vir_venn_data_notp, Species==sp) %>% 
      select(-Species) %>% 
      t()
    colnames(tmp) = tmp[2,]
    tmp = as.data.frame(tmp[3:nrow(tmp),])
    tmp[,1] = ifelse(tmp[,1]=="TRUE", T, F)
    tmp[,2] = ifelse(tmp[,2]=="TRUE", T, F)
    
    ggvenn(tmp, c("Diseased", "Healthy")) + ggtitle(sp)
  }
  
  g1 = draw_venn_sptp("Cat", "Anal swab")
  g2 = draw_venn_sptp("Cat", "Throat swab")
  g3 = draw_venn_sptp("Dog", "Anal swab")
  g4 = draw_venn_sptp("Dog", "Throat swab")
  
  g5 = draw_venn_sp("Cat")
  g6 = draw_venn_sp("Dog")
  
  plot_grid(g1,g2,g3,g4, nrow=2, ncol=2)
  venn_vir_species_type = "output/venn_virus_species_type.pdf"
  ggsave(venn_vir_species_type, width=8, height=7)
  
  plot_grid(g5,g6, nrow=1, ncol=2)
  venn_vir_species = "output/venn_virus_species.pdf"
  ggsave(venn_vir_species, width=6, height=7)
  
  # bacteria
  vir_names = colnames(bac_rpm_table)[2:ncol(bac_rpm_table)]
  vir_venn_data = left_join(lib_ordered, bac_rpm_table, by="Libraries") %>%
    group_by(Species, `Health condition`, Type) %>%
    summarise(across(all_of(vir_names), ~sum(.x)>0))
  vir_venn_data_notp = left_join(lib_ordered, bac_rpm_table, by="Libraries") %>%
    group_by(Species, `Health condition`) %>%
    summarise(across(all_of(vir_names), ~sum(.x)>0))
  
  dc = vir_venn_data_notp[1,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  hc = vir_venn_data_notp[2,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  
  diseased_cat_only_bac = "output/diseased_cat_only_bac.txt"
  vir_names[setdiff(dc, hc)] %>% write_lines(diseased_cat_only_bac)
  
  dd = vir_venn_data_notp[3,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  hd = vir_venn_data_notp[4,3:ncol(vir_venn_data_notp)] %>% t() %>% which()
  
  diseased_dog_only_bac = "output/diseased_dog_only_bac.txt"
  vir_names[setdiff(dd, hd)] %>% write_lines(diseased_dog_only_bac)
  
  draw_venn_sptp = function(sp, tp) {
    tmp = filter(vir_venn_data, Species==sp, Type==tp) %>% 
      select(-Species, -Type) %>% 
      t()
    colnames(tmp) = tmp[2,]
    tmp = as.data.frame(tmp[3:nrow(tmp),])
    tmp[,1] = ifelse(tmp[,1]=="TRUE", T, F)
    tmp[,2] = ifelse(tmp[,2]=="TRUE", T, F)
    
    ggvenn(tmp, c("Diseased", "Healthy")) + ggtitle(str_glue(sp, tp))
  }
  draw_venn_sp = function(sp) {
    tmp = filter(vir_venn_data_notp, Species==sp) %>% 
      select(-Species) %>% 
      t()
    colnames(tmp) = tmp[2,]
    tmp = as.data.frame(tmp[3:nrow(tmp),])
    tmp[,1] = ifelse(tmp[,1]=="TRUE", T, F)
    tmp[,2] = ifelse(tmp[,2]=="TRUE", T, F)
    
    ggvenn(tmp, c("Diseased", "Healthy")) + ggtitle(sp)
  }
  
  g1 = draw_venn_sptp("Cat", "Anal swab")
  g2 = draw_venn_sptp("Cat", "Throat swab")
  g3 = draw_venn_sptp("Dog", "Anal swab")
  g4 = draw_venn_sptp("Dog", "Throat swab")
  
  g5 = draw_venn_sp("Cat")
  g6 = draw_venn_sp("Dog")
  
  
  plot_grid(g1,g2,g3,g4, nrow=2, ncol=2)
  venn_bac_species_type = "output/venn_bac_species_type.pdf"
  ggsave(venn_bac_species_type, width=8, height=7)
  
  plot_grid(g5,g6, nrow=1, ncol=2)
  venn_bac_species = "output/venn_bac_species.pdf"
  ggsave(venn_bac_species, width=6, height=7)
  
  
  output_filenames = c(
    diseased_cat_only_vir,
    diseased_dog_only_vir,
    diseased_cat_only_bac,
    diseased_dog_only_bac,
    venn_vir_species_type,
    venn_vir_species,
    venn_bac_species_type,
    venn_bac_species
  )
  
  output_filenames
}

deseq_analysis_amy = function(lib_meta, vir_rpm_table, bac_rpm_table) {
  # ordering samples (libraries)
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  deseq_analysis = function(rpm_table, species, sample) {
    vir_rpm_table = rpm_table
    
    vir_mat = as.data.frame(vir_rpm_table)
    rownames(vir_mat) = vir_mat$Libraries
    vir_mat$Libraries = NULL
    vir_mat = vir_mat[rowSums(vir_mat)>0,]
    vir_mat = t(vir_mat)
    
    meta_df = as.data.frame(lib_ordered)
    rownames(meta_df) = meta_df$Libraries
    
    meta_df = meta_df[colnames(vir_mat),]
    meta_df$HealthCondition = meta_df$`Health condition`
    
    meta_df = meta_df[meta_df$Species==species & meta_df$Type==sample,]
    vir_mat = vir_mat[,rownames(meta_df)]
    dds = DESeqDataSetFromMatrix(countData = vir_mat+1,
                                 colData = meta_df,
                                 design = ~ HealthCondition)
    
    dds <- DESeq(dds)
    res <- results(dds)
    res = res[order(res$pvalue),]
    res
  }
  
  candidates = expand.grid(species = c("Dog", "Cat"), sample = c("Anal swab", "Throat swab"))
  
  output_filenames = NULL
  for (i in 1:nrow(candidates)) {
    species = candidates$species[i]
    sample = candidates$sample[i]
    
    fname = str_glue("output/DESeq_vir_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    deseq_analysis(vir_rpm_table, species, sample) %>% 
      capture.output(file=fname)
    
    fname = str_glue("output/DESeq_bac_{species}{sample}.txt")
    output_filenames = c(output_filenames, fname)
    deseq_analysis(bac_rpm_table, species, sample) %>% 
      capture.output(file=fname)
  }
  
  output_filenames
}

deseq_analysis_pan = function(lib_meta, bac_rpm_table, fun_rpm_table) {
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  deseq_analysis = function(rpm_table, species, sample) {
    vir_rpm_table = rpm_table
    
    vir_mat = as.data.frame(vir_rpm_table)
    rownames(vir_mat) = vir_mat$Libraries
    vir_mat$Libraries = NULL
    vir_mat = vir_mat[rowSums(vir_mat)>0,]
    vir_mat = t(vir_mat)
    
    meta_df = as.data.frame(lib_ordered)
    rownames(meta_df) = meta_df$Libraries
    
    meta_df = meta_df[colnames(vir_mat),]
    meta_df$HealthCondition = meta_df$`Health condition`
    
    meta_df = meta_df[meta_df$Species==species & meta_df$Type==sample,]
    vir_mat = vir_mat[,rownames(meta_df)]
    dds = DESeqDataSetFromMatrix(countData = round(vir_mat+1),
                                 colData = meta_df,
                                 design = ~ HealthCondition)
    
    dds <- DESeq(dds)
    res <- results(dds)
    res = res[order(res$pvalue),] %>%
      as.data.frame() %>%
      mutate(taxa = rownames(.), .before = 1) %>%
      as_tibble()
  }
  
  candidates = expand.grid(species = c("Dog", "Cat"), sample = c("Anal swab", "Throat swab"))
  
  output_filenames = NULL
  for (i in 1:nrow(candidates)) {
    species = candidates$species[i]
    sample = candidates$sample[i]
    
    fname = str_glue("output/DESeq_pan_bac_{species}{sample}.csv")
    output_filenames = c(output_filenames, fname)
    deseq_analysis(bac_rpm_table, species, sample) %>% 
      write_csv(file=fname)
    
    fname = str_glue("output/DESeq_pan_fun_{species}{sample}.csv")
    output_filenames = c(output_filenames, fname)
    deseq_analysis(fun_rpm_table, species, sample) %>% 
      write_csv(file=fname)
  }
  
  output_filenames
}

gen_rpm_table_pan = function(lib_meta, raw_rpm_table, missing_markers) {
  
  missing_markers = missing_markers %>%
    left_join(select(lib_meta, lib_id = `Libraries-sever`, Libraries))
  # RPM > 100 filter
  rpm_table = raw_rpm_table %>%
    inner_join(select(lib_meta, lib_id = `Libraries-sever`, Libraries), ., by = "lib_id") %>%
    select(-lib_id) %>% 
    pivot_longer(cols = 2:ncol(.)) %>%
    filter(value > 100) %>%
    filter(name != "NA") %>%
    # remove Escherichia (potential contamination)
    filter(name != "Escherichia") %>% 
    # remove unclassified
    filter(!str_detect(name, "^unclassified")) %>%
    semi_join(missing_markers, by = c("name" ="genus", "Libraries")) %>%
    pivot_wider(values_fill = 0)
  
  if ("Mycoplasma" %in% colnames(rpm_table)) {
    rpm_table = rpm_table %>%
      mutate(
        Mycoplasma = Mycoplasma + Metamycoplasma + Mesomycoplasma + Mycoplasmopsis
      )
    rpm_table$Mesomycoplasma = NULL
    rpm_table$Metamycoplasma = NULL
    rpm_table$Mycoplasmopsis = NULL
  }
  rpm_table
}

write_microb_profile = function(lib_meta, rpm_table, filename) {
  outpath = str_glue("output/{filename}")
  lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    inner_join(rpm_table, by = "Libraries") %>%
    write_csv(outpath)
  outpath
}

write_lib_stat = function(lib_meta) {
  x = mutate(lib_meta, tag = paste(Species, `Health condition`, Type, sep="_"))
  table(x$Location, x$tag) %>%
    as.data.frame() %>%
    pivot_wider(names_from = "Var2", values_from = "Freq", values_fill = 0) %>%
    write_csv("output/lib_stat.csv")
  "output/lib_stat.csv"
}

diff_expr = function(lib_meta, rpm_table, species, sample, sp) {
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)
  
  deseq_analysis = function(rpm_table, species, sample) {
    vir_rpm_table = rpm_table
    
    vir_mat = as.data.frame(vir_rpm_table)
    rownames(vir_mat) = vir_mat$Libraries
    vir_mat$Libraries = NULL
    vir_mat = vir_mat[rowSums(vir_mat)>0,]
    vir_mat = t(vir_mat)
    
    meta_df = as.data.frame(lib_ordered)
    rownames(meta_df) = meta_df$Libraries
    
    meta_df = meta_df[colnames(vir_mat),]
    meta_df$HealthCondition = meta_df$`Health condition`
    
    meta_df = meta_df[meta_df$Species==species & meta_df$Type==sample,]
    vir_mat = vir_mat[,rownames(meta_df)]
    dds = DESeqDataSetFromMatrix(countData = round(vir_mat+1),
                                 colData = meta_df,
                                 design = ~ HealthCondition)
    
    dds <- DESeq(dds)
    res <- results(dds)
    res = res[order(res$pvalue),] %>%
      as.data.frame() %>%
      mutate(taxa = rownames(.), .before = 1) %>%
      as_tibble()
    res
  }

  res = deseq_analysis(rpm_table, species, sample)
  
  res_f = res %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(log2FoldChange)
  
  x = rpm_table %>%
    select(Libraries, all_of(res_f$taxa)) %>%
    pivot_longer(cols=2:ncol(.)) %>%
    filter(value > 0) %>%
    pivot_wider(values_fill=0) %>%
    as.data.frame()
  rownames(x) = x$Libraries
  
  x$Libraries = NULL
  x = as.matrix(x)[,res_f$taxa]
  
  
  df_lib_meta = lib_meta %>%
    filter(Species == species, Type == sample) %>%
    filter(Libraries %in% rownames(x)) %>%
    as.data.frame() %>%
    arrange(`Health condition`, Location)
  rownames(df_lib_meta) = df_lib_meta$Libraries
  
  x = log10(x[rownames(df_lib_meta), ]+1) %>% t()
  
  drop_singletons = which(rowSums(x>0)>3)
  res_f = res_f[drop_singletons, ]
  x = x[drop_singletons, ]
  
  col_fun = colorRamp2(c(0, 8), c("#F4F2ED", "#DD1332"))
  
  row_anno = rowAnnotation(log2FC= anno_barplot(res_f$log2FoldChange))
  
  deseq_filename = str_glue("output/diff_stat_{sp}_{species}_{sample}.csv")
  write_csv(res, deseq_filename)
  
  heatmap_filename = str_glue("output/diff_heatmap_{sp}_{species}_{sample}.pdf")
  
  pdf(file = heatmap_filename,       
           width = 6.5,#0.1578947 * ncol(x),
           height = 0.2 * nrow(x))
  ht = Heatmap(x,
          name = "logRPM+1",
          row_order = 1:nrow(x),
          column_order = 1:ncol(x),
          row_split = res_f$log2FoldChange < 0,
          #row_labels = rep("", nrow(x)),
          column_split = df_lib_meta$`Health condition`,
          column_labels = rep("", ncol(x)),
          left_annotation = row_anno,
          col = col_fun
  )
  draw(ht)
  dev.off()
  
  c(deseq_filename, heatmap_filename)
}

merge_rpm_tables = function(vir_rpm_table, bac_rpm_table, fun_rpm_table) {
  vir_rpm_table %>%
    full_join(bac_rpm_table) %>%
    full_join(fun_rpm_table) %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    filter(!is.na(value) & value > 0) %>%
    pivot_wider(values_fill = 0)
}

load_busco_results = function(bac_missing, fun_missing) {
  tbl_fun_missing = read_tsv(fun_missing, col_names = c("path", "missing")) %>%
    mutate(lib_id = str_extract(path, "^busco_fun/(.+?)_.+/.+?", group=1),
           genus = str_extract(path, "^busco_fun/.+?_(.+?)/.+?", group=1)) %>%
    select(-path) %>%
    filter(missing < max(missing))
  
  tbl_bac_missing = read_tsv(bac_missing, col_names = c("path", "missing")) %>%
    mutate(lib_id = str_extract(path, "^(.+?)_.+$", group=1),
           genus = str_extract(path, "^.+?_(.+)$", group=1)) %>%
    select(-path) %>%
    filter(missing < max(missing))
  
  bind_rows(tbl_bac_missing, tbl_fun_missing)
}
