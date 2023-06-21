library(tidyverse)
## This is a loop that constructs plots which can then be called in individual chunks in the RMarkdown.
## Using Bray-Curtis or Jaccard distance for Bacteria and Fungi, respectfully. 
curr.test.type = "pa"

## Also compiling pairwise adonis and betadispersion results

PAdiffs <- list()

count <- 1


for(ds_name in names(ds_list)){
      
      Orddf <- ds_list[[ds_name]]
      rownames(Orddf) <- NULL
      
      tax_toadd <- Orddf[,c("#OTU ID", "taxonomy")]
      ## select just abundance columns & set ESV to rownames
      Orddf <- Orddf %>%
            column_to_rownames("#OTU ID")
      
      Orddf <- Orddf[, colnames(Orddf) %in% metadat$SampleID]
      
      ## Calculate relative abundance for Bacterial, bray-curtis distance, P/A for Fungi jaccard distance
      # if(ds_name == "Bac"){
      #       Orddf <- t(apply(Orddf, 2, function(x){x/sum(x)}))
      #       dist.metric = "bray"
      # }
      # if(ds_name == "Fun"){
      Orddf <- t(apply(Orddf, 2, function(x){ifelse(x > 0, 1, 0)}))
      dist.metric ="jaccard"
      #}
      
      # Add back taxonomy
      Orddf <- left_join(t(Orddf) %>% as.data.frame() %>% rownames_to_column("#OTU ID"), tax_toadd, by = "#OTU ID") %>%
            pivot_longer(cols = any_of(metadat$SampleID), names_to = "SampleID", values_to = "abund")
      
      
      # Test 1 compare 1 week presence/absece from Immediate extraction -----
      
      test.treatments <- unique(metadat$Treatment)[!unique(metadat$Treatment) %in% c("Extraction control", "PCR control")]
      
      test1.data <- Orddf %>% left_join(metadat, by = "SampleID") %>%
            mutate(keep = (Treatment == "Freeze" & Conservation_time_week == 3) | Conservation_time_week == 1 | Conservation_time_week == 0) %>%
            filter(keep & Treatment %in% test.treatments) %>% 
            dplyr::select(-SampleID, -CFL_ID, -Extraction_date, -Depth, -nrow, -keep, -Conservation_time_week) %>%
            pivot_wider(names_from = "Treatment", values_from = "abund") %>%
            pivot_longer(any_of(test.treatments[!test.treatments == "Immediate extraction"]), names_to = "Treatment", values_to = "abund") %>% 
            filter(!is.na(abund))
      
      
      ## Different models for each distinct subset
      
      ## Models by sampletype
      test1.overall <- test1.data %>%
            split(list(.$Sample_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x %>% mutate(Treatment = "Treatment") %>% 
                      pivot_wider(names_from = "Treatment", values_from = "abund") %>% 
                      pivot_longer(all_of(c("Immediate extraction", "Treatment")), names_to = "Treatment", values_to = "abund")
            ) %>%
            map(~ glm(abund~Treatment, data = .x, family = ifelse(ds_name == "Fun", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test1.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Treatment) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      ## Models by sampletype and subtype
      test1.subtype <- test1.data %>%
            split(list(.$Sample_type, .$Soil.tree_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x %>% mutate(Treatment = "Treatment") %>% 
                      pivot_wider(names_from = "Treatment", values_from = "abund") %>% 
                      pivot_longer(all_of(c("Immediate extraction", "Treatment")), names_to = "Treatment", values_to = "abund")
            ) %>%
            map(~ glm(abund~Treatment, data = .x, family = ifelse(curr.test.type == "pa", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   Soil.tree_type = str_extract(Treatment, str_c(unique(metadat$Soil.tree_type), collapse = "|")),
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test1.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Soil.tree_type, Treatment) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      
      #Test 2 compare presence/absence to freezing ----
      
      test.treatments <- unique(metadat$Treatment)[!unique(metadat$Treatment) %in% c("Extraction control", "PCR control", "Immediate extraction")]
      
      test2.data <- Orddf %>% left_join(metadat, by = "SampleID") %>%
            mutate(keep = (Treatment == "Freeze" & Conservation_time_week == 3) | Conservation_time_week == 1 | Conservation_time_week == 0) %>%
            filter(keep & Treatment %in% test.treatments) %>% 
            dplyr::select(-SampleID, -CFL_ID, -Extraction_date, -Depth, -nrow, -keep, -Conservation_time_week) %>%
            pivot_wider(names_from = "Treatment", values_from = "abund") %>%
            pivot_longer(any_of(test.treatments[!test.treatments == "Freeze"]), names_to = "Treatment", values_to = "abund") %>% 
            filter(!is.na(abund))
      
      ## Models by sampletype
      test2.overall <- test2.data %>%
            split(list(.$Sample_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x %>% mutate(Treatment = "Treatment") %>% 
                      pivot_wider(names_from = "Treatment", values_from = "abund") %>% 
                      pivot_longer(all_of(c("Freeze", "Treatment")), names_to = "Treatment", values_to = "abund")
            ) %>%
            map(~ glm(abund~Treatment, data = .x, family = ifelse(curr.test.type == "pa", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test2.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Treatment) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      ## Models by sampletype and subtype
      test2.subtype <- test2.data %>%
            split(list(.$Sample_type, .$Soil.tree_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x %>% mutate(Treatment = "Treatment") %>% 
                      pivot_wider(names_from = "Treatment", values_from = "abund") %>% 
                      pivot_longer(all_of(c("Freeze", "Treatment")), names_to = "Treatment", values_to = "abund")
            ) %>%
            map(~ glm(abund~Treatment, data = .x, family = ifelse(curr.test.type == "pa", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   Soil.tree_type = str_extract(Treatment, str_c(unique(metadat$Soil.tree_type), collapse = "|")),
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test2.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Soil.tree_type, Treatment) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      
      
      
      ### NEED TO RUN TESTS 3 & 6 now.
      # Test 3: Compare presence/absence of all time points for each treatment ----
      test.treatments <- unique(metadat$Treatment)[!unique(metadat$Treatment) %in% c("Extraction control", "PCR control", "Immediate extraction", "Freeze")]
      
      test3.data <- Orddf %>% left_join(metadat, by = "SampleID") %>%
            filter(Treatment %in% test.treatments) %>% 
            dplyr::select(-SampleID, -CFL_ID, -Depth, -nrow, -Extraction_date)  %>%
            filter(!is.na(abund)) %>%
            mutate(Conservation_time_week = as.factor(Conservation_time_week))
      
      ## Models by sampletype
      test3.overall <- test3.data %>%
            split(list(.$Sample_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x ) %>%
            map(~ glm(abund~Conservation_time_week, data = .x, family = ifelse(curr.test.type == "pa", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test3.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   Conservation_time_week = str_extract(parameter, "[[:digit:]]"),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Treatment, Conservation_time_week) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      ## Models by sampletype and subtype
      test3.subtype <- test3.data %>%
            split(list(.$Sample_type, .$Soil.tree_type, .$Treatment, .$`#OTU ID`), drop = TRUE) %>%
            discard(function(x){nrow(x) == 0|sum(x$abund) == 0}) %>%
            map(~ .x %>% mutate(Treatment = "Treatment") %>% 
                      pivot_wider(names_from = "Treatment", values_from = "abund") %>% 
                      pivot_longer(all_of(c("Freeze", "Treatment")), names_to = "Treatment", values_to = "abund")
            ) %>%
            map(~ glm(abund~Treatment, data = .x, family = ifelse(curr.test.type == "pa", "binomial", 'gaussian'))) %>%
            map(summary) %>%
            map('coefficients') %>%
            map(as.data.frame) %>%
            map(~ .x %>% rownames_to_column('parameter')) %>%
            bind_rows(.id="Treatment") %>% 
            filter(!parameter == "(Intercept)") %>%
            rename_with(~ gsub("Pr......", "Pr", .x)) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>%
            mutate(Pr.BH = p.adjust(Pr, method = "BH")) %>% 
            mutate(reg.sig = ifelse(Pr < 0.001, "***", ifelse(Pr<0.01, "**", ifelse(Pr < 0.05, "*", ""))),
                   BH.sig = ifelse(Pr.BH < 0.001, "***", ifelse(Pr.BH<0.01, "**", ifelse(Pr.BH < 0.05, "*", "")))) %>%
            mutate(Sample_type = str_extract(Treatment, str_c(unique(metadat$Sample_type), collapse = "|")), 
                   Soil.tree_type = str_extract(Treatment, str_c(unique(metadat$Soil.tree_type), collapse = "|")),
                   `#OTU ID` = str_extract(Treatment, str_c(unique(test3.data$`#OTU ID`), collapse = "|")), 
                   Treatment = str_extract(Treatment, str_c(unique(metadat$Treatment), collapse = "|")),
                   Conservation_time_week = str_extract(parameter, "[[:digit:]]"),
                   added = Estimate > 0 & Pr < 0.05,
                   removed = Estimate <=0 & Pr <0.05, 
                   added.BH = Estimate > 0 & Pr.BH < 0.05,
                   removed.BH = Estimate <=0 & Pr.BH <0.05) %>%
            left_join(tax_toadd, by = "#OTU ID") %>%
            mutate(genus = gsub("g__", "", str_extract(taxonomy, "g__[[:alnum:]-]+"))) %>%
            group_by(Sample_type, Soil.tree_type, Treatment, Conservation_time_week) %>%
            summarize(changes = paste("+", sum(added), ", -", sum(removed)), 
                      genuschange = paste("+", str_c(unique(genus[added]), collapse = ", "), "; -", str_c(unique(genus[removed]), collapse = ", ")), 
                      changes.BH = paste("+", sum(added.BH), ", -", sum(removed.BH)), 
                      genuschange.BH = paste("+", str_c(unique(genus[added.BH]), collapse = ", "), "; -", str_c(unique(genus[removed.BH]), collapse = ", ")))
      
      ### Compile the info into list 
      
      PAdiffs[[count]] <- list(test1.overall = test1.overall, test1.subtype = test1.subtype,
                               test2.overall = test2.overall, test2.subtype = test2.subtype,
                               test3.overall = test3.overall, test3.subtype = test3.subtype,)
      names(PAdiffs)[count] <- ds_name
      
      count = count + 1
}



rm(ds_name, Orddf, count, test1.data, test2.data, test3.data, test1.overall, test2.overall, test3.overall, test1.subtype, test2.subtype, test3.subtype)

