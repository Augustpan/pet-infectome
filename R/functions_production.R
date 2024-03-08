draw_infectome_heatmap = function(data_pack_pan, pathogenic_taxa) {
  
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  prevalence = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 1000) %>%
    left_join(lib_meta_short) %>%
    group_by(genus) %>%
    summarise(n = n(), n_dogs=sum(Species=="Dog"), n_cats=sum(Species=="Cat"), median_rpm = median(rpm)) %>%
    arrange(desc(n), desc(median_rpm))
  
  # transforming bacteria rpm table
  rpm_bac_pathogens = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 1000) %>%
    filter(genus %in% filter(prevalence, n > 40)$genus) %>%
    ##semi_join(pathogenic_taxa, by = "genus") %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac_pathogens) = rpm_bac_pathogens$Libraries
  rpm_bac_pathogens$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun_pathogens = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    filter(genus %in% c("Malassezia", "Candida")) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun_pathogens) = rpm_fun_pathogens$Libraries
  rpm_fun_pathogens$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir_pathogens = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir_pathogens) = rpm_vir_pathogens$Libraries
  rpm_vir_pathogens$Libraries = NULL
  
  # column ordering
  rpm_bac_pathogens = rpm_bac_pathogens[,order(colSums(rpm_bac_pathogens > 0),  decreasing = T)]
  rpm_fun_pathogens = rpm_fun_pathogens[,order(colSums(rpm_fun_pathogens > 0),  decreasing = T)]
  rpm_vir_pathogens = rpm_vir_pathogens[,order(colSums(rpm_vir_pathogens > 0),  decreasing = T)]
  
  # force row ordering
  rpm_bac_pathogens = rpm_bac_pathogens[lib_meta_short$Libraries,]
  rpm_fun_pathogens = rpm_fun_pathogens[lib_meta_short$Libraries,]
  rpm_vir_pathogens = rpm_vir_pathogens[lib_meta_short$Libraries,]
  
  # merge all tables
  all_rpm_tables = cbind(
    rpm_vir_pathogens,
    rpm_bac_pathogens, 
    rpm_fun_pathogens 
  )
  all_rpm_tables = log10(all_rpm_tables + 1) %>% t()
  
  # heatmap annotations
  
  location_names = unique(lib_meta_short$Location)
  location_colors = brewer.pal(length(location_names), "Set1")
  location_legends = Legend(labels = location_names, 
                            title = "Location", 
                            legend_gp = gpar(fill = location_colors))
  
  health_names = unique(lib_meta_short$`Health condition`)
  health_colors = c("#96C37D", "#FA7F6F")
  health_legends = Legend(labels = health_names, 
                          title = "Health", 
                          legend_gp = gpar(fill = health_colors))
  
  species_names = unique(lib_meta_short$Species)
  species_colors = c("#9AC9DB", "#FFBE7A")
  species_legends = Legend(labels = species_names, 
                           title = "Species", 
                           legend_gp = gpar(fill = species_colors))
  
  sample_annotations = columnAnnotation(
    #  location = anno_simple(lib_meta_short$Location, 
    #                         col = setNames(location_colors, location_names),
    #                         height = unit(0.3, "cm")),
    species = anno_simple(as.vector(lib_meta_short$Species), 
                          col = setNames(species_colors, species_names),
                          height = unit(0.3, "cm")),
    health = anno_simple(as.vector(lib_meta_short$`Health condition`), 
                         col = setNames(health_colors, health_names),
                         height = unit(0.3, "cm")),
    
    gap = unit(1, "mm")
  )
  
  sample_split = paste0(str_replace(lib_meta_short$Type, " swab", ""), "_", lib_meta_short$Species)
  taxa_split = c(
    rep("viruses", ncol(rpm_vir_pathogens)),
    rep("bacteria", ncol(rpm_bac_pathogens)),
    rep("fungi", ncol(rpm_fun_pathogens))
  ) %>%
    factor(levels = c("viruses", "bacteria", "fungi"))
  
  color_func = colorRamp2(c(0, max(all_rpm_tables)), c("#F4F2ED", "#DD1332"))
  
  # drawing heatmap
  output_filename = "output/infectome_heatmap.pdf"
  cairo_pdf(file=output_filename, width=10, height=7.5)
  ht = all_rpm_tables %>%
    Heatmap(., 
            name = "log(RPM+1)", 
            row_order = 1:nrow(.), 
            column_order = 1:ncol(.),
            row_split = taxa_split,
            column_split = sample_split,
            col = color_func,
            column_labels = rep("", ncol(.)),
            top_annotation = sample_annotations,
            row_names_gp = gpar(fontsize = 10, fontface = "italic")
            #row_gap = unit(0, "mm"), 
            #column_gap = unit(0, "mm"), 
            #border = TRUE
    )
  draw(ht, annotation_legend_list = list(species_legends, health_legends))
  dev.off()
  
  output_filename
}

draw_beta_div_pcoa = function(data_pack_pan) {
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac, 
    rpm_fun 
  )

  beta_div_plot = function(rpm_table, name) {
    rpm_table = rpm_table[rowSums(rpm_table) > 0, ]
    
    meta_all = filter(lib_meta_short, Libraries %in% names(rowSums(rpm_table)>0))
    meta_anal = filter(meta_all, Type == "Anal swab")
    meta_throat = filter(meta_all, Type == "Throat swab")
    meta_anal_cat = filter(meta_all, Type == "Anal swab", Species == "Cat")
    meta_throat_cat = filter(meta_all, Type == "Throat swab", Species == "Cat")
    meta_anal_dog = filter(meta_all, Type == "Anal swab", Species == "Dog")
    meta_throat_dog = filter(meta_all, Type == "Throat swab", Species == "Dog")
    
    all_libs = meta_all$Libraries
    anal_libs = meta_anal$Libraries
    throat_libs = meta_throat$Libraries
    anal_cat_libs = meta_anal_cat$Libraries
    throat_cat_libs = meta_throat_cat$Libraries
    anal_dog_libs = meta_anal_dog$Libraries
    throat_dog_libs = meta_throat_dog$Libraries
    
    pcoa_all = vegan::capscale(rpm_table~1, distance = "bray")
    pcoa_anal = vegan::capscale(rpm_table[anal_libs,]~1, distance = "bray")
    pcoa_throat = vegan::capscale(rpm_table[throat_libs,]~1, distance = "bray")
    pcoa_anal_cat = vegan::capscale(rpm_table[anal_cat_libs,]~1, distance = "bray")
    pcoa_throat_cat = vegan::capscale(rpm_table[throat_cat_libs,]~1, distance = "bray")
    pcoa_anal_dog = vegan::capscale(rpm_table[anal_dog_libs,]~1, distance = "bray")
    pcoa_throat_dog = vegan::capscale(rpm_table[throat_dog_libs,]~1, distance = "bray")
    
    axis1_explained = ((pcoa_all$CA$eig / sum(pcoa_all$CA$eig))[1] * 100) %>% round(digits=1) 
    axis2_explained = ((pcoa_all$CA$eig / sum(pcoa_all$CA$eig))[2] * 100) %>% round(digits=1) 
    graph_all = ggplot(meta_all) + 
      geom_point(aes(x=pcoa_all$CA$u[,1], y=pcoa_all$CA$u[,2], color=Type),
                 size = 2) +
      xlab(str_glue("PCoA1 ({axis1_explained}%)")) +
      ylab(str_glue("PCoA2 ({axis2_explained}%)"))
    
    axis1_explained = ((pcoa_anal$CA$eig / sum(pcoa_anal$CA$eig))[1] * 100) %>% round(digits=1) 
    axis2_explained = ((pcoa_anal$CA$eig / sum(pcoa_anal$CA$eig))[2] * 100) %>% round(digits=1) 
    graph_anal = ggplot(meta_anal) + 
      geom_point(aes(x=pcoa_anal$CA$u[,1], y=pcoa_anal$CA$u[,2], 
                     color=Species, shape=`Health condition`),
                 size = 2.5, stroke = 1) +
      stat_ellipse(aes(x=pcoa_anal$CA$u[,1], y=pcoa_anal$CA$u[,2], 
                       color=Species)) +
      scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 1)) +
      scale_color_manual(values = c("Cat" = "#66C2A5", "Dog" = "#FC8D62")) +
      xlab(str_glue("PCoA1 ({axis1_explained}%)")) +
      ylab(str_glue("PCoA2 ({axis2_explained}%)"))
    
    axis1_explained = ((pcoa_throat$CA$eig / sum(pcoa_throat$CA$eig))[1] * 100) %>% round(digits=1) 
    axis2_explained = ((pcoa_throat$CA$eig / sum(pcoa_throat$CA$eig))[2] * 100) %>% round(digits=1) 
    graph_throat = ggplot(meta_throat) + 
      geom_point(aes(x=pcoa_throat$CA$u[,1], y=pcoa_throat$CA$u[,2], 
                     color=Species, shape=`Health condition`),
                 size = 2.5, stroke = 1) +
      stat_ellipse(aes(x=pcoa_throat$CA$u[,1], y=pcoa_throat$CA$u[,2], 
                       color=Species)) +
      scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 1)) +
      scale_color_manual(values = c("Cat" = "#66C2A5", "Dog" = "#FC8D62")) +
      xlab(str_glue("PCoA1 ({axis1_explained}%)")) +
      ylab(str_glue("PCoA2 ({axis2_explained}%)"))
    
    filename_graph_all = str_glue("output/beta_div_{name}_all.pdf")
    filename_graph_anal = str_glue("output/beta_div_{name}_anal.pdf")
    filename_graph_throat = str_glue("output/beta_div_{name}_throat.pdf")
    
    ggsave(filename_graph_all, plot = graph_all, 
           width = 120, height = 90, units = "mm")
    ggsave(filename_graph_anal, plot = graph_anal, 
           width = 120, height = 90, units = "mm")
    ggsave(filename_graph_throat, plot = graph_throat, 
           width = 120, height = 90, units = "mm")
    
    c(filename_graph_all,
      filename_graph_anal,
      filename_graph_throat)
  }
  beta_div_plot_tsne = function(rpm_table, name) {
    rpm_table = rpm_table[rowSums(rpm_table) > 0, ]
    
    meta_all = filter(lib_meta_short, Libraries %in% names(rowSums(rpm_table)>0))
    meta_anal = filter(meta_all, Type == "Anal swab")
    meta_throat = filter(meta_all, Type == "Throat swab")
    
    all_libs = meta_all$Libraries
    anal_libs = meta_anal$Libraries
    throat_libs = meta_throat$Libraries
    
    set.seed(872467)
    
    tsne_all = Rtsne::Rtsne(log1p(rpm_table), check_duplicates = F, max_iter = 5000, perplexity=30)
    tsne_anal = Rtsne::Rtsne(log1p(rpm_table[anal_libs,]), check_duplicates = F, max_iter = 5000, perplexity=10)
    tsne_throat = Rtsne::Rtsne(log1p(rpm_table[throat_libs,]), check_duplicates = F, max_iter = 5000, perplexity=10)
    
    graph_all = ggplot(meta_all) + 
      geom_point(aes(x=tsne_all$Y[,1], y=tsne_all$Y[,2], color=Type),
                 size = 2) +
      xlab(str_glue("t-SNE Axis1")) +
      ylab(str_glue("t-SNE Axis2"))
    
    graph_anal = ggplot(meta_anal) + 
      geom_point(aes(x=tsne_anal$Y[,1], y=tsne_anal$Y[,2], 
                     color=Species, shape=`Health condition`),
                 size = 2.5, stroke = 1) +
      scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 1)) +
      scale_color_manual(values = c("Cat" = "#66C2A5", "Dog" = "#FC8D62")) +
      xlab(str_glue("t-SNE Axis1")) +
      ylab(str_glue("t-SNE Axis2"))
    
    graph_throat = ggplot(meta_throat) + 
      geom_point(aes(x=tsne_throat$Y[,1], y=tsne_throat$Y[,2], 
                     color=Species, shape=`Health condition`),
                 size = 2.5, stroke = 1) +
      scale_shape_manual(values = c("Healthy" = 16, "Diseased" = 1)) +
      scale_color_manual(values = c("Cat" = "#66C2A5", "Dog" = "#FC8D62")) +
      xlab(str_glue("t-SNE Axis1")) +
      ylab(str_glue("t-SNE Axis2"))
    
    filename_graph_all = str_glue("output/beta_div_{name}_all_tsne.pdf")
    filename_graph_anal = str_glue("output/beta_div_{name}_anal_tsne.pdf")
    filename_graph_throat = str_glue("output/beta_div_{name}_throat_tsne.pdf")
    
    ggsave(filename_graph_all, plot = graph_all, 
           width = 120, height = 90, units = "mm")
    ggsave(filename_graph_anal, plot = graph_anal, 
           width = 120, height = 90, units = "mm")
    ggsave(filename_graph_throat, plot = graph_throat, 
           width = 120, height = 90, units = "mm")
    
    c(filename_graph_all,
      filename_graph_anal,
      filename_graph_throat)
  }
  
  c(
    beta_div_plot(rpm_bac, "bac"),
    beta_div_plot_tsne(rpm_vir, "vir")
  )
} 

test_beta_div = function(data_pack_pan) {
  data_pack_pan = tar_read(data_pack_pan)
  
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac, 
    rpm_fun 
  )
  
  rpm_table = rpm_bac
  beta_div_test = function(rpm_table, name) {
    rpm_table = rpm_table[rowSums(rpm_table) > 0, ]
    
    meta_all = filter(lib_meta_short, Libraries %in% names(rowSums(rpm_table)>0))
    meta_anal = filter(meta_all, Type == "Anal swab")
    meta_throat = filter(meta_all, Type == "Throat swab")
    meta_anal_cat = filter(meta_all, Type == "Anal swab", Species == "Cat")
    meta_throat_cat = filter(meta_all, Type == "Throat swab", Species == "Cat")
    meta_anal_dog = filter(meta_all, Type == "Anal swab", Species == "Dog")
    meta_throat_dog = filter(meta_all, Type == "Throat swab", Species == "Dog")
    
    all_libs = meta_all$Libraries
    anal_libs = meta_anal$Libraries
    throat_libs = meta_throat$Libraries
    anal_cat_libs = meta_anal_cat$Libraries
    throat_cat_libs = meta_throat_cat$Libraries
    anal_dog_libs = meta_anal_dog$Libraries
    throat_dog_libs = meta_throat_dog$Libraries
    
    dbrda_all = vegan::dbrda(rpm_table~Type, data=meta_all, distance = "bray")
    dbrda_anal = vegan::dbrda(rpm_table[anal_libs,]~Species*`Health condition`, data=meta_anal, distance = "bray")
    dbrda_throat = vegan::dbrda(rpm_table[throat_libs,]~Species*`Health condition`, data=meta_throat, distance = "bray")
    dbrda_anal_cat = vegan::dbrda(rpm_table[anal_cat_libs,]~`Health condition`, data=meta_anal_cat, distance = "bray")
    dbrda_throat_cat = vegan::dbrda(rpm_table[throat_cat_libs,]~`Health condition`, data=meta_throat_cat, distance = "bray")
    dbrda_anal_dog = vegan::dbrda(rpm_table[anal_dog_libs,]~`Health condition`, data=meta_anal_dog, distance = "bray")
    dbrda_throat_dog = vegan::dbrda(rpm_table[throat_dog_libs,]~`Health condition`, data=meta_throat_dog, distance = "bray")
    
    output_filename = str_glue("output/beta_div_tests_{name}.txt")
    capture.output(list(
      all = anova(dbrda_all),
      anal = anova(dbrda_anal, by = "term"),
      throat = anova(dbrda_throat, by = "term"),
      anal_cat = anova(dbrda_anal_cat),
      anal_dog = anova(dbrda_anal_dog),
      throat_cat = anova(dbrda_throat_cat),
      throat_dog = anova(dbrda_throat_dog)
    ), file = output_filename)
    
    output_filename
  }
  
  c(
    beta_div_test(rpm_bac, "bac"),
    beta_div_test(rpm_vir, "vir")
  )
}

calc_and_comp_alpha_div = function(data_pack_pan) {

  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac, 
    rpm_fun 
  )
  
  alpha_diversity = cbind(
    shannon_bac = vegan::diversity(rpm_bac, index = "shannon"),
    shannon_fun = vegan::diversity(rpm_fun, index = "shannon"),
    shannon_vir = vegan::diversity(rpm_vir, index = "shannon"),
    shannon_all = vegan::diversity(rpm_all, index = "shannon"),
    richness_bac = vegan::specnumber(rpm_bac),
    richness_fun = vegan::specnumber(rpm_fun),
    richness_vir = vegan::specnumber(rpm_vir),
    richness_all = vegan::specnumber(rpm_all),
    reads_bac = rowSums(rpm_bac),
    reads_fun = rowSums(rpm_fun),
    reads_vir = rowSums(rpm_vir),
    reads_all = rowSums(rpm_all)) %>%
    as.data.frame() %>%
    mutate(., Libraries = rownames(.), .before = 1) %>%
    mutate(Libraries = factor(Libraries, levels = lib_meta_short$Libraries)) %>%
    as_tibble() %>%
    left_join(lib_meta_short, by = "Libraries")
  
  alpha_diversity_long = alpha_diversity %>%
    pivot_longer(cols = shannon_bac:reads_all, names_to = "index", values_to = "value")
  
  compare_alpha_diversity = function(index = "richness_bac") {
    alpha_diversity = mutate(alpha_diversity, index_to_test = alpha_diversity[,index][[1]])
    test_results = NULL
    ## SAMPLT TYPE COMPARISON: ANAL SWAB VS THROAT SWAB
    # anal swab vs throat swab - all
    test = coin::wilcox_test(formula = index_to_test ~ Type, data = filter(alpha_diversity))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab vs throat swab - all", statistic=W, pvalue=p))
    
    # anal swab vs throat swab - healthy
    test = coin::wilcox_test(formula = index_to_test ~ Type, data = filter(alpha_diversity, `Health condition` == "Healthy"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab vs throat swab - healthy", statistic=W, pvalue=p))
    
    ## SPECIES COMPARISON: CATS VS DOGS
    # anal swab - healthy cats vs healthy dogs
    test = coin::wilcox_test(formula = index_to_test ~ Species, 
                             data = filter(alpha_diversity, Type == "Anal swab", `Health condition` == "Healthy"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab - healthy cats vs healthy dogs", statistic=W, pvalue=p))
    
    # anal swab - all cats vs all dogs
    test = coin::wilcox_test(formula = index_to_test ~ Species, 
                             data = filter(alpha_diversity, Type == "Anal swab"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab - all cats vs all dogs", statistic=W, pvalue=p))
    
    
    # throat swab - healthy cats vs healthy dogs
    test = coin::wilcox_test(formula = index_to_test ~ Species, 
                             data = filter(alpha_diversity, Type == "Throat swab", `Health condition` == "Healthy"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="throat swab - healthy cats vs healthy dogs", statistic=W, pvalue=p))
    
    # throat swab - all cats vs all dogs
    test = coin::wilcox_test(formula = index_to_test ~ Species, 
                             data = filter(alpha_diversity, Type == "Throat swab"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="throat swab - all cats vs all dogs", statistic=W, pvalue=p))
    
    
    ## HEALTH CONDITIONS COMPARISON: HEALTHY VS DISEASED
    test = coin::wilcox_test(formula = index_to_test ~ `Health condition`, 
                             data = filter(alpha_diversity, Type == "Anal swab", Species == "Cat"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab - healthy cats vs diseased cats", statistic=W, pvalue=p))
    
    test = coin::wilcox_test(formula = index_to_test ~ `Health condition`, 
                             data = filter(alpha_diversity, Type == "Throat swab", Species == "Cat"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="throat swab - healthy cats vs diseased cats", statistic=W, pvalue=p))
    
    
    test = coin::wilcox_test(formula = index_to_test ~ `Health condition`, 
                             data = filter(alpha_diversity, Type == "Anal swab", Species == "Dog"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="anal swab - healthy dogs vs diseased dogs", statistic=W, pvalue=p))
    
    test = coin::wilcox_test(formula = index_to_test ~ `Health condition`, 
                             data = filter(alpha_diversity, Type == "Throat swab", Species == "Dog"))
    W = coin::statistic(test)
    p = coin::pvalue(test)
    test_results = rbind(test_results, c(comparison="throat swab - healthy dogs vs diseased dogs", statistic=W, pvalue=p))
    
    test_results = as_tibble(test_results) %>% 
      mutate(statistic = as.double(statistic), pvalue = as.double(pvalue)) %>%
      mutate(is_significant = pvalue < 0.05) %>%
      mutate(index = index)
    
    test_results
  }
  
  output_filename_tests = "output/alpha_div_comparisons.csv"
  bind_rows(
    compare_alpha_diversity("richness_bac"),
    compare_alpha_diversity("richness_vir")) %>%
    write_csv(output_filename_tests)
  
  output_filename_alpha_div = "output/alpha_div.csv"
  write_csv(alpha_diversity, output_filename_alpha_div)
  
  c(output_filename_tests, output_filename_alpha_div)
}

draw_alpha_div = function(data_pack_pan) {
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac, 
    rpm_fun 
  )
  
  alpha_diversity = cbind(
    shannon_bac = vegan::diversity(rpm_bac, index = "shannon"),
    shannon_fun = vegan::diversity(rpm_fun, index = "shannon"),
    shannon_vir = vegan::diversity(rpm_vir, index = "shannon"),
    shannon_all = vegan::diversity(rpm_all, index = "shannon"),
    richness_bac = vegan::specnumber(rpm_bac),
    richness_fun = vegan::specnumber(rpm_fun),
    richness_vir = vegan::specnumber(rpm_vir),
    richness_all = vegan::specnumber(rpm_all),
    reads_bac = rowSums(rpm_bac),
    reads_fun = rowSums(rpm_fun),
    reads_vir = rowSums(rpm_vir),
    reads_all = rowSums(rpm_all)) %>%
    as.data.frame() %>%
    mutate(., Libraries = rownames(.), .before = 1) %>%
    mutate(Libraries = factor(Libraries, levels = lib_meta_short$Libraries)) %>%
    as_tibble() %>%
    left_join(lib_meta_short, by = "Libraries")
  
  alpha_diversity_long = alpha_diversity %>%
    pivot_longer(cols = shannon_bac:reads_all, names_to = "index", values_to = "value")
  
  
  draw_alpha_div_index = function(index) {
    alpha_diversity = mutate(alpha_diversity, index_to_draw = alpha_diversity[,index][[1]])
    
    type_comparison = alpha_diversity %>%
      ggplot() + 
      geom_boxplot(aes(x = Type, y = index_to_draw, fill = Type, color = Type)) + 
      scale_color_manual(values = c("Anal swab" = "#eb736c", "Throat swab" = "#00bfc3")) +
      scale_fill_manual(values = c("Anal swab" = "#ffa174", "Throat swab" = "#33d7db")) +
      ylab(index) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
    
    species_comparison = alpha_diversity %>% 
      filter(`Health condition` == "Healthy") %>%
      ggplot() + 
      geom_boxplot(aes(x = Type, color = Species, fill=Species, y = index_to_draw), alpha = 0.5) + 
      scale_color_manual(values = c("Cat" = "#66C2A5", "Dog" = "#FC8D62")) +
      scale_fill_manual(values = c("Cat" = "#77e8c6", "Dog" = "#ffac68")) +
      theme(legend.position = "none", 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    
    health_comparison = alpha_diversity %>% 
      ggplot() + 
      geom_boxplot(aes(x = paste0(Species, "_", Type), color = `Health condition`, fill = `Health condition`, y = index_to_draw), alpha = 0.5) + 
      scale_color_manual(values = c("Healthy" = "#00bfc3", "Diseased" = "#f97670")) +
      scale_fill_manual(values =c("Healthy" = "#00bfc3", "Diseased" = "#f97670")) +
      ylab(index) +
      theme(legend.position = "none", 
            axis.title.x = element_blank())
    
    cowplot::plot_grid(
      cowplot::plot_grid(type_comparison, species_comparison, nrow = 1, rel_widths = c(1, 1.65)),
      health_comparison, 
      nrow = 2)
    output_filename = str_glue("output/alpha_div_{index}.pdf")
    ggsave(output_filename, width = 95, height = 88, units = "mm")
    
  }
  
  c(
    draw_alpha_div_index("richness_bac"),
    draw_alpha_div_index("richness_vir")
  )
}

draw_taxa_composition = function(data_pack_pan, virus_metadata, tax_map) {
  data_pack_pan = tar_read(data_pack_pan)
  
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac,
    rpm_fun
  )
  
  taxa_rank = rpm_all %>%
    mutate(Libraries = rownames(.), .before = 1) %>%
    as_tibble() %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    filter(value > 0) %>%
    left_join(distinct(bind_rows(tax_map,
                                 select(virus_metadata, genus = virus_name, family)), 
                       name=genus, taxa=family), by = "name") %>%
    group_by(Libraries, taxa) %>%
    summarise(rpm = sum(value)) %>%
    ungroup() %>%
    left_join(lib_meta_short, by = "Libraries") %>%
    group_by(Species, `Health condition`, Type, taxa) %>%
    summarise(n = n(), median_rpm = median(rpm)) %>%
    ungroup() %>%
    mutate(is_virus = taxa %in% virus_metadata$family) %>%
    group_by(Species, `Health condition`, Type, is_virus) %>%
    slice_max(tibble(n, median_rpm), 
              n = 5, 
              with_ties = F) %>%
    ungroup() %>%
    mutate(mark = T)
  
  
  rpm_all %>%
    mutate(Libraries = rownames(.), .before = 1) %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    left_join(distinct(tax_map, genus, taxa=family), by = "genus") %>%
    group_by(Libraries, taxa) %>%
    summarise(rpm = sum(rpm)) %>%
    ungroup() %>%
    mutate(Libraries = factor(Libraries, levels = lib_meta_short$Libraries)) %>%
    left_join(lib_meta_short) %>%
    left_join(select(taxa_rank, -n, -median_rpm)) %>%
    replace_na(list(mark = F)) %>%
    mutate(taxa_new = ifelse(mark, taxa, "Others")) %>%
    ggplot() + 
    geom_bar(aes(x=Libraries, y=rpm,fill=taxa_new), stat="identity", position = "fill") +
    facet_wrap(~Type+Species+`Health condition`, nrow=2, scales = "free_x") +
    #scale_fill_manual(values = c("Others" = "gray")) + 
    guides(fill=guide_legend(ncol=2)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(), panel.background = element_blank())
  output_filename = "output/taxa_composition.pdf"
  ggsave(output_filename, width = 220, height = 80, units= "mm")
  output_filename
}

prepare_lefse_data = function(data_pack_pan, virus_metadata, tax_map) {
  
  tax_info = distinct(
      rbind(select(tax_map, name=genus, family, kingdom),
            select(virus_metadata, name=virus_name, family) %>% 
              mutate(kingdom="Viruses"))
    )
  
  lib_meta_short = data_pack_pan$lib_meta %>%
    select(Libraries, Species, `Health condition`, Location, Type) %>%
    mutate(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")),
           Species = factor(Species, levels = c("Cat", "Dog")),
           Type = factor(Type, levels = c("Throat swab", "Anal swab"))) %>%
    arrange(Type, Species, `Health condition`, Location, Libraries)
  lib_meta_ids = select(lib_meta_short, Libraries)
  
  # transforming bacteria rpm table
  rpm_bac = data_pack_pan$rpm_table_bac %>%
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_bac) = rpm_bac$Libraries
  rpm_bac$Libraries = NULL
  
  # transforming fungi rpm table
  rpm_fun = data_pack_pan$rpm_table_fun %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "genus", values_to = "rpm") %>%
    filter(rpm > 0) %>%
    pivot_wider(names_from = "genus", values_from = "rpm", values_fill = 0) %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_fun) = rpm_fun$Libraries
  rpm_fun$Libraries = NULL
  
  # transforming virus rpm table
  rpm_vir = data_pack_pan$rpm_table_vir  %>%
    left_join(lib_meta_ids, ., by = "Libraries") %>%
    mutate(across(2:ncol(.), ~replace_na(.x, 0))) %>%
    as.data.frame()
  rownames(rpm_vir) = rpm_vir$Libraries
  rpm_vir$Libraries = NULL
  
  # force row ordering
  rpm_bac = rpm_bac[lib_meta_short$Libraries,]
  rpm_fun = rpm_fun[lib_meta_short$Libraries,]
  rpm_vir = rpm_vir[lib_meta_short$Libraries,]
  
  rpm_all = cbind(
    rpm_vir,
    rpm_bac, 
    rpm_fun 
  ) %>%
    mutate(Libraries = rownames(.), .before = 1) %>%
    as_tibble() %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    filter(value > 0) %>%
    left_join(tax_info, by = "name") %>%
    mutate(new_name = str_replace_all(str_glue("{kingdom}|{family}|{name}"), " ", "_")) %>%
    select(Libraries, col_name=new_name, value) %>%
    pivot_wider(names_from = Libraries, values_from= "value", values_fill = 0)
  
  meta_info = bind_rows(
    lib_meta_short %>%
      select(Libraries, Species) %>%
      pivot_wider(names_from = "Libraries", values_from = "Species") %>%
      mutate(col_name = "species",.before = 1),
    lib_meta_short %>%
      select(Libraries, `Health condition`) %>%
      pivot_wider(names_from = "Libraries", values_from = "Health condition") %>%
      mutate(col_name = "health",.before = 1),
    lib_meta_short %>%
      select(Libraries, Type) %>%
      pivot_wider(names_from = "Libraries", values_from = "Type") %>%
      mutate(col_name = "type",.before = 1),
    lib_meta_short %>%
      select(Libraries, Location) %>%
      pivot_wider(names_from = "Libraries", values_from = "Location") %>%
      mutate(col_name = "location",.before = 1)
  )
  
  meta_all = lib_meta_short
  meta_anal_cat = filter(meta_all, Type == "Anal swab", Species == "Cat")
  meta_throat_cat = filter(meta_all, Type == "Throat swab", Species == "Cat")
  meta_anal_dog = filter(meta_all, Type == "Anal swab", Species == "Dog")
  meta_throat_dog = filter(meta_all, Type == "Throat swab", Species == "Dog")
  
  all_libs = meta_all$Libraries
  anal_cat_libs = meta_anal_cat$Libraries
  throat_cat_libs = meta_throat_cat$Libraries
  anal_dog_libs = meta_anal_dog$Libraries
  throat_dog_libs = meta_throat_dog$Libraries
  
  data = rbind(meta_info, rpm_all) %>%
    dplyr::rename(lib_id = col_name)
  
  ds = select(data, lib_id, all_of(anal_cat_libs))
  f = rowSums(mutate_all(ds[5:nrow(ds),2:ncol(ds)], as.double)) > 0
  f = c(F, T, F, F, f)
  ds = ds[f,]
  write_tsv(ds, "output/lefse_input_anal_cat.txt")
  
  ds = select(data, lib_id, all_of(anal_dog_libs))
  f = rowSums(mutate_all(ds[5:nrow(ds),2:ncol(ds)], as.double)) > 0
  f = c(F, T, F, F, f)
  ds = ds[f,]
  write_tsv(ds, "output/lefse_input_anal_dog.txt")
  
  ds = select(data, lib_id, all_of(throat_cat_libs))
  f = rowSums(mutate_all(ds[5:nrow(ds),2:ncol(ds)], as.double)) > 0
  f = c(F, T, F, F, f)
  ds = ds[f,]
  write_tsv(ds, "output/lefse_input_throat_cat.txt")
  
  ds = select(data, lib_id, all_of(throat_dog_libs))
  f = rowSums(mutate_all(ds[5:nrow(ds),2:ncol(ds)], as.double)) > 0
  f = c(F, T, F, F, f)
  ds = ds[f,]
  write_tsv(ds, "output/lefse_input_throat_dog.txt")
  
  c(
    "output/lefse_input_anal_cat.txt",
    "output/lefse_input_anal_dog.txt",
    "output/lefse_input_throat_cat.txt",
    "output/lefse_input_throat_dog.txt"
  )
}

load_lefse_results = function(filenames) {
  lefse_results = NULL
  for (type in names(filenames)) {
    tmp = read_tsv(filenames[type], 
                   col_names = c("taxa", "raw_lda_score", "enrich", "lda_score", "pvalue"),
                   col_types = "cdcdd") %>%
      filter(pvalue < 0.05) %>%
      filter(!is.na(enrich)) %>%
      select(-raw_lda_score) %>%
      mutate(level = str_count(taxa, "\\.")+1) %>%
      mutate(kingdom = str_split(taxa, "\\.", simplify = T)[,1]) %>%
      mutate(family = str_split(taxa, "\\.", simplify = T)[,2]) %>%
      mutate(name = str_split(taxa, "\\.", simplify = T)[,3]) %>%
      mutate(name = str_replace_all(name, "_", " ")) %>%
      mutate(split = type)
    lefse_results = rbind(lefse_results, tmp)
  }
  lefse_results
}

draw_diff_heatmap_lefse = function(rpm_table, lib_meta, lefse_results, sel_species, sample) {
  if (sel_species == "Cat" & sample == "Anal swab") {
    s_type = "anal_cat"
  } else if (sel_species == "Cat" & sample == "Throat swab") {
    s_type = "throat_cat"
  } else if (sel_species == "Dog" & sample == "Anal swab") {
    s_type = "anal_dog"
  } else if (sel_species == "Dog" & sample == "Throat swab") {
    s_type = "throat_dog"
  } else {
    NULL
  }
  
  lib_ordered = arrange(lib_meta, Species, desc(`Health condition`), desc(Type), Location)

  res_f = lefse_results %>%
    filter(level == 3) %>%
    mutate(taxa = name) %>%
    filter(split == s_type) %>%
    mutate(lda_score = lda_score * ifelse(enrich == "Diseased", -1, 1)) %>%
    arrange(desc(lda_score))
  
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
    filter(Species == sel_species, Type == sample) %>%
    filter(Libraries %in% rownames(x)) %>%
    as.data.frame() %>%
    arrange(`Health condition`, Location)
  rownames(df_lib_meta) = df_lib_meta$Libraries
  
  x = log10(x[rownames(df_lib_meta), ]+1) %>% t()
  
  res_f = res_f[which(rowSums(x)>0), ]
  x = x[which(rowSums(x)>0), ]
  
  col_fun = colorRamp2(c(0, 8), c("#F4F2ED", "#DD1332"))
  
  row_anno = rowAnnotation(`LDA score` = anno_barplot(res_f$lda_score))
  
  heatmap_filename = str_glue("output/lefse_heatmap_{sel_species}_{sample}.pdf")
  
  pdf(file = heatmap_filename,       
      width = 6.5,#0.1578947 * ncol(x),
      height = max(0.2 * nrow(x),2))
  ht = Heatmap(x,
               name = "logRPM+1",
               row_order = 1:nrow(x),
               column_order = 1:ncol(x),
               row_split = res_f$enrich,
               #row_labels = rep("", nrow(x)),
               column_split = df_lib_meta$`Health condition`,
               column_labels = rep("", ncol(x)),
               left_annotation = row_anno,
               col = col_fun
  )
  draw(ht)
  dev.off()
  heatmap_filename
}

# merge DESeq2 and LefSe results
draw_diff_heatmap_deseq2_and_lefse = function(rpm_table, lib_meta, lefse_results, sel_species, sample) {
  if (sel_species == "Cat" & sample == "Anal swab") {
    s_type = "anal_cat"
  } else if (sel_species == "Cat" & sample == "Throat swab") {
    s_type = "throat_cat"
  } else if (sel_species == "Dog" & sample == "Anal swab") {
    s_type = "anal_dog"
  } else if (sel_species == "Dog" & sample == "Throat swab") {
    s_type = "throat_dog"
  } else {
    NULL
  }
  species = sel_species
  
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
  
  res_deseq = deseq_analysis(rpm_table, species, sample)
  
  res_f_deseq = res_deseq %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(log2FoldChange) %>%
    mutate(enrich = ifelse(log2FoldChange > 0, "Healthy", "Diseased"))
  
  res_f = lefse_results %>%
    filter(level == 3) %>%
    mutate(taxa = name) %>%
    filter(split == s_type) %>%
    mutate(lda_score = lda_score * ifelse(enrich == "Diseased", -1, 1)) %>%
    arrange(desc(lda_score)) %>%
    left_join(select(res_deseq, taxa, log2FoldChange))
  
  res_f = bind_rows(res_f_deseq, res_f) %>%
    select(taxa, log2FoldChange, enrich) %>%
    distinct() %>%
    arrange(log2FoldChange) %>%
    mutate(flag = ifelse(taxa %in% res_f_deseq$taxa, 
                         "deseq2",
                         "")) %>%
    mutate(flag = ifelse(taxa %in% res_f$taxa, 
                         paste0(flag, "lefse"),
                         flag))
  
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
    filter(Species == sel_species, Type == sample) %>%
    filter(Libraries %in% rownames(x)) %>%
    as.data.frame() %>%
    arrange(`Health condition`, Location)
  rownames(df_lib_meta) = df_lib_meta$Libraries
  
  x = log10(x[rownames(df_lib_meta), ]+1) %>% t()
  
  res_f = res_f[which(rowSums(x)>0), ]
  x = x[which(rowSums(x)>0), ]
  
  col_fun = colorRamp2(c(0, 8), c("#F4F2ED", "#DD1332"))
  
  row_anno = rowAnnotation(`log2FC` = anno_barplot(res_f$log2FoldChange))
  
  heatmap_filename = str_glue("output/lefse_add_deseq_heatmap_{sel_species}_{sample}.pdf")
  
  pdf(file = heatmap_filename,       
      width = 6.5,#0.1578947 * ncol(x),
      height = max(0.2 * nrow(x),2))
  ht = Heatmap(x,
               name = "logRPM+1",
               row_order = 1:nrow(x),
               column_order = 1:ncol(x),
               row_split = res_f$enrich,
               #row_labels = rep("", nrow(x)),
               row_labels = res_f$flag,
               column_split = df_lib_meta$`Health condition`,
               column_labels = rep("", ncol(x)),
               left_annotation = row_anno,
               col = col_fun
  )
  draw(ht)
  dev.off()
  heatmap_filename
}

diffential_abundance_test = function(data_pack) {
  # input:
  rpm_table = data_pack$rpm_table_vir %>%
    left_join(data_pack$rpm_table_bac, by = "Libraries") %>%
    left_join(data_pack$rpm_table_fun, by = "Libraries")
  lib_meta = data_pack$lib_meta
  
  do_test = function(`Health condition`, value){
    data = tibble(`Health condition` = `Health condition`, value = log10(value + 1))
    
    x = filter(data, `Health condition` == "Healthy")$value
    y = filter(data, `Health condition` == "Diseased")$value
    
    test_result = try(wilcox.test(x, y), silent = T)
    p.value = try(test_result$p.value, silent = T)
    if (inherits(p.value, "try-error")) {
      p.value = NA 
    }
    
    return(p.value)
  }
  
  test_results = rpm_table %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    #filter(value > 0) %>%
    left_join(select(lib_meta, Libraries, Species, Type, `Health condition`), by = join_by(Libraries)) %>%
    group_by(Species, Type, name) %>%
    summarise(wilcox_p = do_test(`Health condition`, value), 
              n_healthy = sum(`Health condition` == "Healthy" & value > 0),
              n_diseased = sum(`Health condition` == "Diseased" & value > 0),
              n_total = n_healthy + n_diseased) %>%
    ungroup() %>%
    arrange(wilcox_p)
  
  draw_boxplot = function(`Health condition`, value) {
    data = tibble(`Health condition` = factor(`Health condition`, levels = c("Healthy", "Diseased")), 
                  value = log10(value+1))
    ggplot(data) +
      geom_boxplot(aes(x = `Health condition`, y = value)) +
      geom_jitter(aes(x = `Health condition`, y = value), height = 0, width = 0.1, size = 1.5) + 
      ggpubr::stat_compare_means(aes(x = `Health condition`, y = value), label.y = max(data$value) * 1.1) + 
      ylim(c(0, max(data$value) * 1.2)) +
      ylab("logRPM+1")
  }
  
  boxplots = rpm_table %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    #filter(value > 0) %>%
    left_join(select(lib_meta, Libraries, Species, Type, `Health condition`), by = join_by(Libraries)) %>%
    #semi_join(filter(test_results, wilcox_p < 0.05), by = join_by(name, Species, Type)) %>%
    group_by(Species, Type, name) %>%
    summarise(plot = list(draw_boxplot(`Health condition`, value) + 
                            ggtitle(str_glue("{name}\n{Species} - {Type}")) + 
                            theme(plot.title = element_text(size = 10)))) %>%
    ungroup() %>%
    left_join(test_results, by = join_by(name, Species, Type)) %>%
    arrange(Species, Type, name)
  
  n_signif = nrow(filter(boxplots, !is.na(wilcox_p) & wilcox_p < 0.05))
  n_non_signif = nrow(filter(boxplots, !is.na(wilcox_p) & wilcox_p >= 0.05))
  
  # save figures
  dir.create("output/daa", showWarnings = FALSE)
  g1 = cowplot::plot_grid(plotlist = filter(boxplots, !is.na(wilcox_p) & wilcox_p < 0.05)$plot, ncol = 6)
  s = 1.75
  output_plot_signif = "output/daa/boxplot-signif.pdf"
  ggsave(filename = output_plot_signif, 
         plot = g1, 
         width = 180 * s, 
         height = 270 * s, 
         units = "mm",
         limitsize = FALSE)
  
  batch_size = 40
  s = 1.35
  output_plot_non_signif = boxplots %>%
    filter(!is.na(wilcox_p) & wilcox_p >= 0.05) %>%
    mutate(batch = cut(1:n_non_signif, breaks = 0:(ceiling(n_non_signif / batch_size)) * batch_size)) %>%
    group_by(batch) %>%
    summarise(graph = list(cowplot::plot_grid(plotlist = plot, ncol = 4))) %>%
    mutate(batch_id = 1:nrow(.)) %>%
    rowwise() %>%
    mutate(filename = str_glue("output/daa/boxplot-non-signif-batch{batch_id}.pdf"),
           ggsave(filename = filename, 
                  plot = graph, 
                  width = 180 * s, 
                  height = ceiling(batch_size / 4) * 55 * s, 
                  units = "mm", 
                  limitsize = FALSE))
  
  # save table
  output_table = "output/daa/test-results.csv"
  boxplots %>%
    select(-plot) %>%
    write_csv(output_table)
  
  c(
    output_plot_non_signif$filename,
    output_plot_signif,
    output_table
  )
}

draw_diff_heatmap_deseq2_and_lefse_and_manual_daa = function(rpm_table, lib_meta, lefse_results, sel_species, sample) {
  if (sel_species == "Cat" & sample == "Anal swab") {
    s_type = "anal_cat"
  } else if (sel_species == "Cat" & sample == "Throat swab") {
    s_type = "throat_cat"
  } else if (sel_species == "Dog" & sample == "Anal swab") {
    s_type = "anal_dog"
  } else if (sel_species == "Dog" & sample == "Throat swab") {
    s_type = "throat_dog"
  } else {
    NULL
  }
  species = sel_species
  
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
    
    dds = DESeq(dds)
    res = results(dds)
    res = res[order(res$pvalue),] %>%
      as.data.frame() %>%
      mutate(taxa = rownames(.), .before = 1) %>%
      as_tibble()
    res
  }
  
  res_deseq = deseq_analysis(rpm_table, species, sample)
  
  res_f_deseq = res_deseq %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(log2FoldChange) %>%
    mutate(enrich = ifelse(log2FoldChange > 0, "Healthy", "Diseased"))
  
  res_f = lefse_results %>%
    filter(level == 3) %>%
    mutate(taxa = name) %>%
    filter(split == s_type) %>%
    mutate(lda_score = lda_score * ifelse(enrich == "Diseased", -1, 1)) %>%
    arrange(desc(lda_score)) %>%
    left_join(select(res_deseq, taxa, log2FoldChange))
  
  daa = read_csv("output/daa/test-results.csv") %>%
    filter(Type == sample, Species == sel_species) %>%
    filter(wilcox_p < 0.05) %>%
    select(taxa = name) %>%
    left_join(res_deseq)
  
  res_f = bind_rows(res_f_deseq, res_f, daa) %>%
    distinct(taxa, .keep_all = T) %>%
    select(taxa, log2FoldChange, enrich) %>%
    distinct() %>%
    arrange(log2FoldChange) %>%
    mutate(flag = ifelse(taxa %in% res_f_deseq$taxa, 
                         "deseq2/",
                         "")) %>%
    mutate(flag = ifelse(taxa %in% res_f$taxa, 
                         paste0(flag, "lefse/"),
                         flag)) %>%
    mutate(flag = ifelse(taxa %in% daa$taxa, 
                         paste0(flag, "daa/"),
                         flag)) %>%
    mutate(enrich = ifelse(log2FoldChange < 0, "Diseased", "Healthy"))
  
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
    filter(Species == sel_species, Type == sample) %>%
    filter(Libraries %in% rownames(x)) %>%
    as.data.frame() %>%
    arrange(`Health condition`, Location)
  rownames(df_lib_meta) = df_lib_meta$Libraries
  
  x = log10(x[rownames(df_lib_meta), ]+1) %>% t()
  
  res_f = res_f[which(rowSums(x)>0), ]
  x = x[which(rowSums(x)>0), ]
  res_f = res_f[rowSums(x > 0) >= 5, ]
  x = x[rowSums(x > 0) >= 5, ]
  
  col_fun = colorRamp2(c(0, 8), c("#F4F2ED", "#DD1332"))
  
  row_anno = rowAnnotation(`log2FC` = anno_barplot(res_f$log2FoldChange))

  draw__ = function(x, heatmap_filename, row_labels) {
    pdf(file = heatmap_filename,       
        width = 6.5,#0.1578947 * ncol(x),
        height = max(0.2 * nrow(x)+0.6, 2))
    ht = Heatmap(x,
                 name = "logRPM+1",
                 row_order = 1:nrow(x),
                 column_order = 1:ncol(x),
                 row_split = res_f$enrich,
                 row_labels = row_labels,
                 #row_labels = rep("", nrow(x)),
                 #row_labels = res_f$flag,
                 column_split = df_lib_meta$`Health condition`,
                 column_labels = rep("", ncol(x)),
                 left_annotation = row_anno,
                 col = col_fun
    )
    draw(ht)
    dev.off()
    return(heatmap_filename)
  }
  dir.create("output/daa/heatmap-tag", showWarnings = FALSE)
  dir.create("output/daa/heatmap-name", showWarnings = FALSE)
  dir.create("output/daa/heatmap-plot", showWarnings = FALSE)
  c(
    draw__(x, str_glue("output/daa/heatmap-plot/3methods_heatmap_{sel_species}_{sample}.pdf"),
         rep("", nrow(x))),
    draw__(x, str_glue("output/daa/heatmap-tag/3methods_heatmap_{sel_species}_{sample}.pdf"),
         res_f$flag),
    draw__(x, str_glue("output/daa/heatmap-name/3methods_heatmap_{sel_species}_{sample}.pdf"),
         rownames(x))
  )
}
