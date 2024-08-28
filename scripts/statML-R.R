library(tidyverse)
library(ggpubr)
library(readxl)
library(dunn.test)
library(ggplot2)
library(ggsignif)
library(readxl)
library(car)
library(patchwork)




get_stats <- function(df, value_column, grouping_column) {
  
  plot_df <- df %>% 
    group_by(!!sym(grouping_column)) %>%
    summarise( 
      n = n(),
      mean = mean(!!sym(value_column)),
      sd = sd(!!sym(value_column))
    ) %>%
    mutate(SEM = sd / sqrt(n)) %>%
    mutate(MOE = SEM * qt((1 - 0.05) / 2 + 0.5, n - 1))
  
  return(plot_df)
  
}



test_multi_groups <- function(df, value_column, grouping_column, parametric = TRUE, paired = FALSE, adjustment.method = 'bonferroni') {
  
  setClass(
    "statistic",
    representation(
      test = "character",
      leven_var_test = "list",
      posthoc_test = "character",
      test_data = "list",
      posthoc_data = 'list'
      
    )
  )
  
  if (parametric == FALSE) {
    
    
    levene_result <- leveneTest(df[[value_column]] ~ df[[grouping_column]])
    
    if (levene_result$`Pr(>F)`[1] < 0.05) {
      
      info <- "Levene's test (p < 0.05): variance not equal."
      
    } else {
      
      info <- "Levene's test (p > 0.05): variance equal."
      
    }
    
    formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
    
    kruskal_result <- kruskal.test(formula, data = df)
    
    kruskal_list <- list(
      statistic = kruskal_result$statistic,
      df = kruskal_result$parameter,
      p.value = kruskal_result$p.value
    )
    
    
    
   
    pairwise_results <- pairwise.wilcox.test(df[[value_column]], df[[grouping_column]], p.adjust.method = 'none', paired = FALSE)
      
    
    
    p_val = c()
    pair1 = c()
    pair2 = c()
    
    for (i in 1:nrow(pairwise_results$p.value)) {
      for (j in 1:ncol(pairwise_results$p.value)) {
        if (!is.na(pairwise_results$p.value[i, j])) {
          pair1 <- c(pair1, rownames(pairwise_results$p.value)[i])
          pair2 <- c(pair2,  colnames(pairwise_results$p.value)[j])
          p_val <- c(p_val, pairwise_results$p.value[i, j])
        }
      }
    }
    
    posthc_results <- list(
      p.adjusted = p_val,
      pair1 = pair1,
      pair2 = pair2,
      adjustment = tolower(adjustment.method)
    )
    
    
    if (paired == TRUE) {
      posthoc_test = 'Wilcoxon Signed-Rank'
    } else if (paired == FALSE) {
      posthoc_test = 'Mann-Whitney U'
    }
    
    statistic <- new("statistic", 
                     test ='Kruskal-Wallis', 
                     leven_var_test = list('levene_results' = levene_result, 'response' = info),
                     posthoc_test = posthoc_test,
                     test_data = as.list(kruskal_list),
                     posthoc_data = posthc_results)
  
  } else if (parametric == TRUE) {
    
    
    levene_result <- leveneTest(df[[value_column]] ~ df[[grouping_column]])
    
    if (levene_result$`Pr(>F)`[1] <= 0.05) {
      
      info <- "Levene's test (p < 0.05): variance not equal."
      
      formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
      
      welch_anova <- oneway.test(formula, data = df, var.equal = FALSE)
      
      
      aov_results <- list(
        statistic =  welch_anova$statistic,
        df = welch_anova$parameter,
        p.value = welch_anova$p.value
        
      )
      
      
      pairwise_results <- pairwise.t.test(df[[value_column]], df[[grouping_column]], p.adjust.method = 'none', paired = FALSE, pool.sd = FALSE)
      
      p_val = c()
      pair1 = c()
      pair2 = c()
      
      for (i in 1:nrow(pairwise_results$p.value)) {
        for (j in 1:ncol(pairwise_results$p.value)) {
          if (!is.na(pairwise_results$p.value[i, j])) {
            pair1 <- c(pair1, rownames(pairwise_results$p.value)[i])
            pair2 <- c(pair2,  colnames(pairwise_results$p.value)[j])
            p_val <- c(p_val, pairwise_results$p.value[i, j])
          }
        }
      }
      
      
      posthc_results <- list(
        p.adjusted = p_val,
        pair1 = pair1,
        pair2 = pair2,
        adjustment = tolower(adjustment.method)
      )
      
      
      
      statistic <- new("statistic", 
                       test = "Welch's ANOVA",
                       leven_var_test = list('levene_results' = levene_result, 'response' = info),
                       posthoc_test = "Welch's T-test",
                       test_data = as.list(aov_results),
                       posthoc_data = posthc_results)
    
    } else {
      
      info <- "Levene's test (p > 0.05): variance equal."
      formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
      
      aov <- aov(formula, data = df)
      
      aov_results = summary(aov)
      aov_results <- list(
        statistic = aov_results[[1]][["F value"]][1],
        df = aov_results[[1]][["Df"]][1],
        p.value = aov_results[[1]][["Pr(>F)"]][1]
      )
      
      
      pairwise_results <- pairwise.t.test(df[[value_column]], df[[grouping_column]], p.adjust.method = 'none', paired = FALSE, pool.sd = FALSE)
      
      p_val = c()
      pair1 = c()
      pair2 = c()
      

            for (i in 1:nrow(pairwise_results$p.value)) {
        for (j in 1:ncol(pairwise_results$p.value)) {
          if (!is.na(pairwise_results$p.value[i, j])) {
            pair1 <- c(pair1, rownames(pairwise_results$p.value)[i])
            pair2 <- c(pair2,  colnames(pairwise_results$p.value)[j])
            p_val <- c(p_val, pairwise_results$p.value[i, j])
          }
        }
      }
      
      
      posthc_results <- list(
        p.adjusted = p_val,
        pair1 = pair1,
        pair2 = pair2,
        adjustment = tolower(adjustment.method)
      )
      
      
      
      
      statistic <- new("statistic", 
                       test = 'ANOVA',
                       leven_var_test = list('levene_results' = levene_result, 'response' = info),
                       posthoc_test = 'T-test',
                       test_data = as.list(aov_results),
                       posthoc_data = posthc_results)
      
      
    }
    

    
    
  }
  
  return(statistic)
  
  
}
  




test_two_groups <- function(df, value_column, grouping_column, parametric = TRUE, paired = FALSE) {
  
  setClass(
    "statistic",
    representation(
      test = "character",
      p.val = "ANY",
      statistic = "ANY",
      paired = 'ANY'
      
    )
  )
  
  if (parametric == FALSE) {
    
   
    

    wcox <- wilcox.test(df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[1]], df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[2]], alternative = 'two.side', paired)
    
   
    
    if (paired == TRUE) {
      test_type = 'Wilcoxon Signed-Rank'
    } else if (paired == FALSE) {
      test_type = 'Mann-Whitney U'
    }
    
    statistic <- new("statistic", 
                     test = test_type,
                     p.val = wcox$p.value,
                     statistic = wcox$statistic,
                     paired = paired
                     )
    
    
    
  } else if (parametric == TRUE) {
    
    
    tt <- t.test(df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[1]], df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[2]], alternative = 'two.side', paired, var.equal = TRUE)
    
    
    
    
    statistic <- new("statistic", 
                     test = "T-test",
                     p.val = tt$p.value,
                     statistic = tt$statistic,
                     paired = paired
    )
    
    
  }
  
  return(statistic)
  
}




two_groups_analysis <- function(value_column, grouping_column, data, bar_queue = NaN, x_label = '', x_angle = 30, y_label = '', size = 10, parametric = TRUE, paired = FALSE, brew_colors = 'Dark2') {
  
  
  if (length(unique(data[[grouping_column]])) == 2) {
    
    plot_df = get_stats(data, value_column, grouping_column)
    
    
    if (!unique(is.na(bar_queue)) & !is.list(bar_queue) & length(bar_queue) == length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = bar_queue)
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = bar_queue)
      
    } else if (!unique(is.na(bar_queue)) & !is.list(bar_queue) & length(bar_queue) != length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])
      
      
      print('Warning! The `bar_queue` length in not equal with number of groups!')
      
    } else {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])
      
    }
    
    
    results = test_two_groups(data, value_column, grouping_column, parametric, paired)
    
    
    list_of_comparison <- list(c(levels(plot_df[[grouping_column]])))
    p_value <- results@p.val
    
    if (p_value < 0.001) {
      sig <- '***'
    } else if (p_value < 0.01) {
      sig <- '**'
    } else if (p_value < 0.05) {
      sig <- '*'
    } else {
      sig <- 'ns'
      
    }
    
    
  
    
    
    bar_plot = ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
      geom_bar(stat = "identity", show.legend = FALSE, width = 0.4, color = "black")+
      geom_errorbar(aes(ymin = mean-SEM, ymax = mean+SEM), width = 0.2)
    
    
    box_plot <- ggplot(data, mapping =  aes(y = !!sym(value_column), x = !!sym(grouping_column), fill =  !!sym(grouping_column))) +
      geom_boxplot() + 
      theme_minimal()
    
    max_y <- max(plot_df$mean + plot_df$SEM)
    min_y <- min(plot_df$mean - plot_df$SEM)
    
    
    y_pos <- c(max_y*1.05)
    
    
    bar_plot = bar_plot + geom_signif(comparisons = list_of_comparison,
                                      annotations  = sig,
                                      y_position = y_pos,
                                      map_signif_level = FALSE, textsize = 3) +
      coord_cartesian(ylim = c(0, max_y*1.1 )) + 
      annotate("text", x = -Inf, y = Inf, label = paste(results@test, 'p.val =', results@p.val), hjust = -0.05, vjust = 1.25, size = 3)
    
    
    
    max_y <- max(data$size + plot_df$SEM)
    min_y <- min(data$size - plot_df$SEM)
    
    
    y_pos <- c(max_y*1.05)
    
    box_plot = box_plot +  geom_signif(comparisons = list_of_comparison,
                                       annotations  = sig,
                                       y_position = y_pos,
                                       map_signif_level = FALSE, textsize = 3) +
      coord_cartesian(ylim = c(min_y, max_y*1.1)) + 
      annotate("text", x = -Inf, y = Inf, label = paste(results@test, 'p.val =', results@p.val), hjust = -0.05, vjust = 1.25, size = 3)
    
    
    bar_plot = bar_plot + ylab(y_label)
    bar_plot = bar_plot + xlab(x_label)
    
    box_plot = box_plot + ylab(x_label)
    box_plot = box_plot + xlab(y_label)
    
    
    bar_plot = bar_plot +  theme_classic() +
      theme(axis.title.y = element_text(size = size),   
            axis.title.x = element_text(size = size)) + 
      scale_x_discrete(guide = guide_axis(angle = x_angle)) +
      scale_fill_brewer(palette=brew_colors)
    
    box_plot = box_plot +  theme_classic() +
      theme(axis.title.y = element_text(size = size),   
            axis.title.x = element_text(size = size)) + 
      scale_x_discrete(guide = guide_axis(angle = x_angle)) +
      scale_fill_brewer(palette=brew_colors) +
      theme(legend.position="none") 
    
    
    
    
    
    setClass(
      "two_groups_analysis",
      representation(
        bar_plot = "ANY",
        box_plot = 'ANY',
        statistic_tests = "ANY",
        statistic_data = "list"
        
      )
    )
    
    results <- new("two_groups_analysis",
                   bar_plot = bar_plot,
                   box_plot = box_plot,
                   statistic_tests = results,
                   statistic_data = plot_df)
    
    return(results)
  
  } else if (length(unique(data[[grouping_column]])) > 2) {
    
    stop("The number of groups in the analysis is greater than 2.\n   For more than two groups use the multi_groups_analysis() function")
    
  } else {
    
    stop("The number of groups in the analysis is wrong. Check grouping_column")
    
    
  }
  
}





################################################################################

multi_groups_analysis <- function(value_column, grouping_column, data, bar_queue = NaN, x_label = '', x_angle = 30, y_label = '', size = 10, bar_size = 0.5, parametric = FALSE, paired = FALSE, include_ns = FALSE, bars = 'sem', bars_size = 1, adjustment.method = 'bonferroni', y_break = NaN, brew_colors = 'Dark2') {


  if (length(unique(data[[grouping_column]])) > 2) {

    plot_df = get_stats(data, value_column, grouping_column)

    if (!unique(is.na(bar_queue)) & !is.list(bar_queue) & length(bar_queue) == length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = bar_queue)
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = bar_queue)

    } else if (!unique(is.na(bar_queue)) & !is.list(bar_queue) & length(bar_queue) != length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])


      print('Warning! The `bar_queue` length in not equal with number of groups!')

    } else {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])

    }



    results = test_multi_groups(data, value_column, grouping_column, parametric, paired, adjustment.method)


    results_tmp <- as.data.frame(results@posthoc_data)

    order_sequence <- levels(plot_df[[grouping_column]])


    results_tmp <- results_tmp[order(match(results_tmp$pair1, order_sequence), match(results_tmp$pair2, (order_sequence))), ]




    list_of_comparison <- list()
    p_values <- c()
    sigs <- c()

    for (i in 1:length(results_tmp$pair1)) {


      list_of_comparison[[i]] <- c(results_tmp$pair1[i], results_tmp$pair2[i])
      p_value <- results_tmp$p.adjusted[i]

      if (p_value < 0.001) {
        sig <- '***'
      } else if (p_value < 0.01) {
        sig <- '**'
      } else if (p_value < 0.05) {
        sig <- '*'
      } else {
        sig <- 'ns'

      }

      p_values <- c(p_values, p_value)
      sigs <- c(sigs, sig)

    }


    if (include_ns == FALSE) {

      to_rm <- c()

      for (s in 1:length(sigs)) {
        if (sigs[s] == 'ns') {

          to_rm <- c(to_rm, s)
        }

      }


      if (length(to_rm) > 0) {
        sigs <- sigs[-to_rm]
        p_values <- p_values[-to_rm]
        list_of_comparison <- list_of_comparison[-to_rm]
      }

    }




    ################################################################################

    if (bars == 'sd') {
      
      bar_plot = ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = bar_size - 0.1, size = bars_size)
      
    } else {
      
      bar_plot = ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-SEM, ymax = mean+SEM), width = bar_size - 0.1, size = bars_size)
      
    }
    
    



    if (!is.na(y_break)) {

      bar_plot <- bar_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }


    MinMeanSEMMax <- function(x) {
      v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
      names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
      v
    }
    
    if (bars == 'sd') {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_boxplot(width = bar_size - 0.1, size = bars_size) +
        theme_minimal()
      
    } else {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_boxplot(width = bar_size - 0.1, size = bars_size) +
        stat_summary(fun.data=MinMeanSEMMax, geom="boxplot") +
        theme_minimal()
      
      
    }
    


   

    if (!is.na(y_break)) {

      box_plot <- box_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }



    if (length(list_of_comparison) > 0) {
      
      
      if (bars == 'sd') {
        
        max_y <- max(plot_df$mean + plot_df$sd)
        min_y <- min(plot_df$mean - plot_df$sd)
        
      } else {
        
        max_y <- max(plot_df$mean + plot_df$SEM)
        min_y <- min(plot_df$mean - plot_df$SEM)
        
        
      }
      



      y_pos <- c()
      fc = 0
      for (o in 1:length(list_of_comparison)) {
        if (o == 1) {
          y_pos <- c(y_pos, 0)
  
        } else {
          fc = fc + 10
          y_pos <- c(y_pos, 0 + fc)
  
        }
  
      }

      

    signif_plot <- ggplot(plot_df, aes(x = !!sym(grouping_column), y = 0)) +
      geom_blank() +
      geom_signif(comparisons = list_of_comparison,
                  annotations  = sigs,
                  y_position = y_pos,
                  map_signif_level = FALSE, textsize = 4) +
      coord_cartesian(ylim = c(0, round(max(y_pos)*1.5))) +
      annotate("text", x = -Inf, y = Inf, label = results@leven_var_test$response, hjust = -0.05, vjust = 1.25, size = 3) +
      annotate("text", x = -Inf, y = Inf, label = paste(results@test, 'p.val =', results@test_data$p.value), hjust = -0.05, vjust = 2.75, size = 3) +
      annotate("text", x = -Inf, y = Inf, label = paste('Post-hoc:',results@posthoc_test, ' p.adj:', results@posthoc_data$adjustment), hjust = -0.06, vjust = 4.25, size = 3) +
      theme_void()




    } else {
      
      
      if (bars == 'sd') {
        
        max_y <- max(plot_df$mean + plot_df$sd)
        min_y <- min(plot_df$mean - plot_df$sd)
        
      } else {
        
        max_y <- max(plot_df$mean + plot_df$SEM)
        min_y <- min(plot_df$mean - plot_df$SEM)
        
        
      }




      bar_plot = bar_plot +
        coord_cartesian(ylim = c(0, round(max(max_y) * 1.08))) +
        annotate("text", x = -Inf, y = Inf, label = results@leven_var_test$response, hjust = -0.05, vjust = 1.25, size = 3) +
        annotate("text", x = -Inf, y = Inf, label = paste(results@test, 'p.val =', results@test_data$p.value), hjust = -0.05, vjust = 2.75, size = 3) +
        annotate("text", x = -Inf, y = Inf, label = paste('Post-hoc:',results@posthoc_test, ' p.adj:', results@posthoc_data$adjustment), hjust = -0.06, vjust = 4.25, size = 3)



      max_y <- max(data[[value_column]])
      min_y <- min(data[[value_column]])





      box_plot = box_plot +
        coord_cartesian(ylim = c(min_y, round(max(max_y)* 1.08))) +
        annotate("text", x = -Inf, y = Inf, label = results@leven_var_test$response, hjust = -0.05, vjust = 1.25, size = 3) +
        annotate("text", x = -Inf, y = Inf, label = paste(results@test, 'p.val =', results@test_data$p.value), hjust = -0.05, vjust = 2.75, size = 3) +
        annotate("text", x = -Inf, y = Inf, label = paste('Post-hoc:',results@posthoc_test, ' p.adj:', results@posthoc_data$adjustment), hjust = -0.06, vjust = 4.25, size = 3)




    }

    bar_plot = bar_plot + ylab(y_label)
    bar_plot = bar_plot + xlab(x_label)

    box_plot = box_plot + ylab(y_label)
    box_plot = box_plot + xlab(x_label)


    bar_plot = bar_plot +  theme_classic() +
             theme(axis.title.y = element_text(size = size),
                   axis.title.x = element_text(size = size)) +
             scale_x_discrete(guide = guide_axis(angle = x_angle)) +
             scale_fill_brewer(palette=brew_colors)

    box_plot = box_plot +  theme_classic() +
      theme(axis.title.y = element_text(size = size),
            axis.title.x = element_text(size = size)) +
      scale_x_discrete(guide = guide_axis(angle = x_angle)) +
      scale_fill_brewer(palette=brew_colors) +
      theme(legend.position="none")

    if (length(list_of_comparison) > 0) {
      box_plot <- signif_plot + box_plot  + plot_layout(ncol = 1)
      bar_plot <- signif_plot + bar_plot  + plot_layout(ncol = 1)
    }


    setClass(
      "multi_groups_analysis",
      representation(
        bar_plot = "ANY",
        box_plot = 'ANY',
        statistic_tests = "ANY",
        statistic_data = "list"

      )
    )

    results <- new("multi_groups_analysis",
                       bar_plot = bar_plot,
                       box_plot = box_plot,
                       statistic_tests = results,
                       statistic_data = plot_df)

    return(results)

  } else if (length(unique(data[[grouping_column]])) == 2) {

    stop("The number of groups in the analysis is equal to 2.\n   For two groups use the two_groups_analysis() function")

  } else {

    stop("The number of groups in the analysis is wrong. Check grouping_column")


  }

}







