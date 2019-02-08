#Script to plot and analyze ephys results on intrinsic membrane properties
#Author: Kirill Zaslavsky

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
                                                             #if not using R-studio, set working directory manually
                                                             #using setwd()

library(data.table)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(dplyr)
library(scales)
library(ggthemes)
library(purrr)
source("SparCon_fn.R")
dir.create("Intrinsic Ephys")

ephys_db <- fread("intrinsic_ephys_db.txt")



#Define categorical variables
#Define factor levels
levels_SHANK2_gtype <- c("WT", "HET", "NULL")

levels_comparison_var <- c( "SHANK2 R841X","SHANK2 DEL" ,"SHANK2 KO")

levels_comparison_label <- c("CTRL","SHANK2 R841X",
                             
                             "SHANK2 DEL","SHANK2 KO"
)


#prepare lists for pipe
dat_list <- list(ephys_db)
cb_list <- list(list("SHANK2_gtype",levels_SHANK2_gtype),
                list("comparison_label",levels_comparison_label),
                list("comparison_var", levels_comparison_var))

#pipe
ephys_db_cat <- dat_list %>% 
  map(sparc_extract) %>%
  map(sparc_categorize, fac = cb_list)

ephys_db <- ephys_db_cat[[1]]

#####
# Make graphs


####list of variables
var_list <- as.list(names(dplyr::select(ephys_db, AP_height:Resting_membrane_potential)))
title_list <- unlist(lapply((str_replace_all(names(dplyr::select(ephys_db, AP_height:Resting_membrane_potential)), "_", " ")), 
                            simpleCap))
ylab_list <- c("mV", "mV", "mV", "ms", "ms", "ms", expression("G"*Omega), "mV")
filenames <- file.path("Intrinsic Ephys",  paste0(title_list, ".pdf"))

##looping function to generate the graphs
pl_loop <- function(i){  
 
pl <- ggplot(data = ephys_db, aes_string(x = "comparison_label", y = var_list[[i]]))

out_pl <- pl %>% exp.plfn.geomPOINT(size = 3.5, alpha = 0.6, jitterwidth = 0.3, meanwidth = 0.5, sewidth = 0.2,
                                    scale = 1.25)+
  labs(x="Group", y = ylab_list[i])+
  scale_y_continuous(labels = fmt.lengthen())+
  scale_fill_manual("Genotype", breaks=c("WT","HET","NULL"), values = color2[c(3,1,2)])+
  ggtitle(title_list[i])+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(.~comparison_var, scale = "free")


  ggsave(filenames[i], out_pl, width = 6, height = 5)
}

for (i in seq_along(var_list))
{
  pl_loop(i)
}

##do t-tests for each paired comparison
#get groups for comparison, var_list from above, create results table
comp_groups <- as.character(unique(ephys_db[, comparison_var]))
results_table <- data.frame(CTRL_grp_name = character(0),
                            MUT_grp_name = character(0),
                            Variable = character(0),
                            CTRL_estimate = numeric(0),
                            MUT_estimate = numeric(0),
                            CI_lower = numeric(0),
                            CI_upper = numeric(0),
                            df = numeric(0),
                            t = numeric(0),
                            p_val = numeric(0)
                            )

for (i in seq_along(comp_groups))
{
  comp_group_var <- comp_groups[i]
  for (j in seq_along(var_list))
  {
    ctrl <- ephys_db[comparison_var == comp_group_var & ind_type == "CTRL", .SD, .SDcols = var_list[[j]] ]
    mut <- ephys_db[comparison_var == comp_group_var & ind_type == "MUT", .SD, .SDcols = var_list[[j]] ]
    test_result <- t.test(ctrl,mut)
    row <- data.frame(CTRL_grp_name = (unique(ephys_db[comparison_var == comp_group_var & ind_type == "CTRL", (comparison_label)])),
                      MUT_grp_name = (unique(ephys_db[comparison_var == comp_group_var & ind_type == "MUT", (comparison_label)])),
                      Variable = var_list[[j]],
                      CTRL_estimate = test_result[["estimate"]][[1]],
                      MUT_estimate = test_result[["estimate"]][[2]],
                      CI_lower = test_result[["conf.int"]][[1]],
                      CI_upper = test_result[["conf.int"]][[2]],
                      df = test_result[["parameter"]],
                      t = test_result[["statistic"]],
                      p_val = test_result[["p.value"]]
            )
    results_table <- rbind(results_table, row)
  }
}

write.table(results_table, file = "Intrinsic Ephys/Intrinsic_ephys_t_tests.txt", row.names = FALSE, sep = "\t")

