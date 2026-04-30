
#--------------
# ESTIMATE IBS
#---------------
estimate_IBS <- function(vcf){
   
   ibs.matrix <- MIPanalyzer::get_IBS_distance(x = vcf,
                                 ignore_het = FALSE,
                                 report_progress = FALSE)
   
   return(ibs.matrix)
}

#------------------
# ESTIMATE MLE IBD
#------------------
estimate_MLE_IBD <- function(vcf){
   
   ibd.matrix <- MIPanalyzer::inbreeding_mle(x = vcf,
                                             f = seq(0.01, 0.99, 0.01),
                                             ignore_het = FALSE,
                                             report_progress = FALSE)
   
   return(ibd.matrix)
}


#------------------
# ESTIMATE HMM IBD
#------------------
estimate_HMM_IBD <- function(vcf){
   
   #----------------
   # hmmIBD setup
   #----------------
   
   #----------
   # Part 1
   #----------
   ad <- vcfR::extract.gt(vcf, element = "AD")
   dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = T)
   altad <- vcfR::masplit(ad, record = 2, sort = F)
   wsaf <- altad/dp
   colnames(wsaf) <- colnames(vcf@gt)[2:ncol(vcf@gt)]
   
   # Impute missing values
   locus_impute <- apply(wsaf, 2, median, na.rm = TRUE)
   locus_impute <- outer(rep(1, nrow(wsaf)), locus_impute)
   wsaf[is.na(wsaf)] <- locus_impute[is.na(wsaf)]
   
   #--------------
   # Part 2
   #--------------
   gtmat <- dplyr::bind_cols("chrom" = vcfR::getCHROM(vcf), "pos" = vcfR::getPOS(vcf), wsaf)
   gtmat <- gtmat %>%
      tidyr::pivot_longer(., cols = !c("chrom", "pos")) %>%
      dplyr::mutate(chrom = as.numeric(factor(chrom)),
                    gthmm = round(value, 0)) %>%
      dplyr::select(-c("value")) %>%
      tidyr::pivot_wider(data = ., names_from = "name", values_from = "gthmm")
   
   # Save genotype matrix for hmmIBD
   readr::write_tsv(x = gtmat, file = "gt_matrix_for_hmmIBD.txt", col_names = T)
   
   #----------
   # Part 3 
   #----------
   altafvec <- rowMeans(wsaf, na.rm = T)
   afmat <- dplyr::bind_cols("chrom" = vcfR::getCHROM(vcf),
                             "pos" = vcfR::getPOS(vcf),
                             1-altafvec, altafvec) %>%
      dplyr::mutate(chrom = as.numeric(factor(chrom)))
   
   # Save allele frequency file for hmmIBD
   readr::write_tsv(afmat, file = "af_matrix_for_hmmIBD.txt", col_names = F)
   
   #----------------
   # Running hmmIBD
   #----------------
   tf <- tempfile(pattern = "output")
   out <- hmmibdr::hmm_ibd(input_file = "gt_matrix_for_hmmIBD.txt",
                           allele_freqs =  "af_matrix_for_hmmIBD.txt",
                           max_fit_iterations = 20,
                           rec_rate = 1e-2, # note the small recombo rate relative to what would be expected in malaria
                           output_file = tf)
   
   return(out)
}




format_matrix <- function(ibd.data, metadata, gcdistance){
   
   # --------------
   # data munging
   # --------------
   mtdt_x <- metadata %>% 
      dplyr::select(-VisitDate, -latitude, -longitude) %>% 
      dplyr::rename(p1 = SampleID)
   
   mtdt_y <- metadata %>% 
      dplyr::select(-VisitDate, -latitude, -longitude) %>% 
      dplyr::rename(p2 = SampleID)
   
   # -------------------------
   # Add metadata information
   # -------------------------
   ibd_metadata <- dplyr::left_join(ibd.data, mtdt_x, by = "p1") %>%
      dplyr::left_join(., mtdt_y, by = "p2") %>%
      rename_with(., ~gsub(".x", "_p1", .x, fixed = T)) %>%
      rename_with(., ~gsub(".y", "_p2", .x, fixed = T))
   
   
   ibd_metadata <- ibd_metadata  %>%
      dplyr::left_join(rbind(gcdistance, transform(gcdistance, 
                                                   VillageCode_p1 = VillageCode_p2, 
                                                   VillageCode_p2 = VillageCode_p1))) %>% 
      mutate(dist_km = if_else(is.na(dist_km), 0, dist_km)) %>% 
      mutate(Compoud = case_when((VillageCode_p1 == VillageCode_p2 & CompoundCode_p1 == CompoundCode_p2) ~ 'Within Compound', 
                                 (VillageCode_p1 == VillageCode_p2 & CompoundCode_p1 != CompoundCode_p2) ~ 'Between Compouds'),
             Village = if_else(VillageCode_p1 == VillageCode_p2, 'Within Village', 'Between Village'),
             region = if_else(region_p1 == region_p2, 'Within Region', 'Between Region'),
             village_pair = if_else(dist_km <= 10, 'Village Pair', 'Non-pair Village'))
   
   return(ibd_metadata)
}


recap <- function(ibd.data){
   
   household <- ibd.data %>%
      filter((VillageCode_p1 == VillageCode_p2) & (HHCode_p1 == HHCode_p2)) %>% 
      summarise(n = n(), avg = mean(ibd), med = median(ibd),
                sd = sd(ibd), se = sd/sqrt(n), ymin = avg - se, 
                ymax = avg + se, U95CI = avg + 1.96 * se, 
                L95CI = avg - 1.96 * se) %>% 
      mutate(site = "household") %>% relocate(site, .before = n)
   
   compound <- ibd.data %>%
      filter((VillageCode_p1 == VillageCode_p2) & (CompoundCode_p1 == CompoundCode_p2)) %>% 
      summarise(n = n(), avg = mean(ibd), med = median(ibd),
                sd = sd(ibd), se = sd/sqrt(n), ymin = avg - se, 
                ymax = avg + se, U95CI = avg + 1.96 * se, 
                L95CI = avg - 1.96 * se) %>% 
      mutate(site = "compound") %>% relocate(site, .before = n)
   
   village <- ibd.data %>%
      filter(VillageCode_p1 == VillageCode_p2) %>% 
      summarise(n = n(), avg = mean(ibd), med = median(ibd),
                sd = sd(ibd), se = sd/sqrt(n), ymin = avg - se, 
                ymax = avg + se, U95CI = avg + 1.96 * se, 
                L95CI = avg - 1.96 * se) %>% 
      mutate(site = "village") %>% relocate(site, .before = n)
   
   region <- ibd.data %>%
      filter(region_p1 == region_p2) %>% 
      summarise(n = n(), avg = mean(ibd), med = median(ibd),
                sd = sd(ibd), se = sd/sqrt(n), 
                ymin = avg - se, ymax = avg + se, 
                U95CI = avg + 1.96 * se, L95CI = avg - 1.96 * se) %>% 
      mutate(site = "region") %>% relocate(site, .before = n)
   
   results <- rbind.data.frame(household, compound, village, region)
   return(results)
}

#---------------------------------

get_stats_for_index <- function(data){
   
   results <- data %>% 
   summarise(n = n(), avg = mean(ibd), med = median(ibd), sd = sd(ibd), se = sd/sqrt(n), 
          ymin = avg - se, ymax = avg + se, U95CI = avg + 1.96 * se, L95CI = avg - 1.96 * se,
          .groups = "drop")
   
   return(results)
}

#---------------------------------

my_theme <- theme(legend.position = "none",
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 10, color = 'black', face = "bold"),
            axis.title = element_text(size = 12, color = 'black', face = "bold"),
            legend.text = element_text(size = 10, color = 'black', face = "bold"),
            legend.title = element_text(size = 12, color = 'black', face = "bold", hjust = .5))

my_theme2 <- theme(legend.position = "top",
   axis.line = element_line(color = 'black'),
                  axis.text = element_text(size = 10, color = 'black', face = "bold"),
                  axis.title = element_text(size = 12, color = 'black', face = "bold"),
                  legend.text = element_text(size = 10, color = 'black', face = "bold"),
                  legend.title = element_text(size = 12, color = 'black', face = "bold", hjust = .5))

compare_methods_v1 <- function(df){
   
   ggplot(df, aes(x = method, y = ibd, fill = method)) +
      stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = .5) + # mean_se mean_cl_boot median_hilow
      stat_summary(fun.data = mean_se, geom = "bar", show.legend = FALSE) +
      viridis::scale_fill_viridis(discrete = T, option = "E") +
      xlab("") + ylab("Average IBD values") +
      guides(fill = guide_legend(title = "Methods")) +
      theme_bw() +
      my_theme
}

compare_methods_v2 <- function(df){
   
      ggplot(df) +
      geom_errorbar(aes(x=method, ymin=L95CI, ymax=U95CI), linewidth = 0.4,
                    width = 0.2, position = position_dodge(0.9)) +
      geom_bar(aes(x = method, y = avg, fill = method),
               stat="identity", position="dodge", width = .7) +
      viridis::scale_fill_viridis(discrete = T, option = "E") +
      xlab("") + ylab("Average IBD values") +
      guides(fill = guide_legend(title = "Methods")) +
      theme_bw() +
      my_theme
}


compare_methods_v3 <- function(df){
   
   ggplot(df, aes(x = method, y = avg, fill = method)) +
      geom_errorbar(aes(x=method, ymin=L95CI, ymax=U95CI), linewidth = 0.4,
                    width = 0.2, position = position_dodge(0.9)) +
      geom_col() +
      viridis::scale_fill_viridis(discrete = T, option = "E") +
      xlab("") + ylab("Average IBD values") +
      guides(fill = guide_legend(title = "Methods")) +
      theme_bw() +
      my_theme
}

#---------------------------------

barplot.fill <- function(df){
   
   ggplot(df) +
      geom_bar(aes(x = {method}, y = {avg}, fill = {site}), stat="identity", position="fill", width = .6) +
      viridis::scale_fill_viridis(discrete = T, option = "E") +
      xlab("") + ylab("Mean IBD values") +
      guides(fill = guide_legend(title = "Relatedness level")) +
      theme_bw() +
      my_theme2
}


barplot.dodge <- function(df){
   
   ggplot(df) +
      geom_bar(aes(x = method, y = avg, fill = site), stat="identity", position="dodge", width = .6) +
      viridis::scale_fill_viridis(discrete = T, option = "E") +
      xlab("") + ylab("Mean IBD values") +
      guides(fill = guide_legend(title = "Relatedness level")) +
      theme_bw() +
      my_theme2
}


recap2 <- function(ibd.data){
   
   village <- ibd.data %>%
      filter((VillageCode_p1 == VillageCode_p2)) %>% 
      group_by(VillageCode_p1) %>% 
      summarise(n = n(), avg = mean(ibd), 
                med = median(ibd), sd = sd(ibd), se = sd/sqrt(n),
                ymin = avg - se, ymax = avg + se,
                U95CI = avg + 1.96 * se, L95CI = avg - 1.96 * se,
                .groups = 'drop') %>% 
      rename(site = VillageCode_p1) %>% mutate_all(~replace_na(., 0))
   
   
   regions <- ibd.data %>%
      filter(region_p1 == region_p2) %>%
      group_by(region_p1) %>%
      summarise(n = n(), avg = mean(ibd),
                med = median(ibd), sd = sd(ibd), se = sd/sqrt(n),
                ymin = avg - se, ymax = avg + se,
                U95CI = avg + 1.96 * se, L95CI = avg - 1.96 * se) %>%
      rename(site = region_p1) %>% mutate_all(~replace_na(., 0))
   
   results <- rbind.data.frame(village, regions)
   return(results)
}
