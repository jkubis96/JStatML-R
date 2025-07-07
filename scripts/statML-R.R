local({
  
  packages <- c("tidyverse", "ggpubr", "readxl", "ggplot2", "dplyr", "purrr", "rlang",
                "ggsignif", "car", "patchwork", "ARTool", "stats", "RColorBrewer")
  
  
  installed <- packages %in% installed.packages()
  if (any(!installed)) {
    install.packages(packages[!installed])
  }
  
  lapply(packages, library, character.only = TRUE)
})




#' Compute summary statistics by group
#'
#' The \code{get_stats} function calculates summary statistics for a specified numeric column,
#' grouped by a given categorical column. It returns the number of observations, mean,
#' standard deviation, standard error of the mean (SEM), and margin of error (MOE) for each group.
#'
#' @param df A data frame containing the input data.
#' @param value_column A string specifying the name of the numeric column for which statistics are computed.
#' @param grouping_column A string specifying the name of the column used to group the data.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{grouping_column} — the column of grouping variables,
#'   \item \code{n} — number of observations in the group,
#'   \item \code{mean} — mean value,
#'   \item \code{sd} — standard deviation,
#'   \item \code{SEM} — standard error of the mean,
#'   \item \code{MOE} — margin of error (95% confidence interval).
#' }
#'
#' @import dplyr
#' @importFrom stats qt
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B"), each = 10),
#'   value = rnorm(20)
#' )
#' get_stats(df, "value", "group")
get_stats <- function(df, 
                      value_column, 
                      grouping_column) {
  
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


#' Calculate average fold changes between group means
#'
#' The \code{avg_FC} function computes pairwise fold changes between group means
#' from a summary data frame. It assumes the input data contains a grouping variable
#' and corresponding means (e.g., from \code{get_stats()}), and returns a data frame
#' with fold changes and their log2 transformations between all group combinations.
#'
#' @param data A data frame where the first column is the group identifier
#' and one of the columns must be named \code{mean}, representing group means.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{group1} — name of the first group,
#'   \item \code{group2} — name of the second group,
#'   \item \code{fold_change} — ratio of the mean of \code{group1} to the mean of \code{group2},
#'   \item \code{avg_logFC} — log2-transformed fold change.
#' }
#'
#' @details This is a **supporting function** used internally by \code{\link{two_groups_analysis}}
#' and \code{\link{multi_groups_analysis}}. It computes all pairwise fold changes and their reciprocals
#' (i.e., both group1/group2 and group2/group1), then applies a log2 transformation to produce \code{avg_logFC}.
#'
#' @importFrom utils combn
#' @export
#'
#' @examples
#' df_summary <- data.frame(
#'   group = c("A", "B", "C"),
#'   mean = c(5, 10, 20)
#' )
#' avg_FC(df_summary)
avg_FC <- function(data) {
  
  colnames(data)[1] <- 'group'
  
  results <- data.frame(
    group1 = character(),
    group2 = character(),
    avg_fold_change = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  combinations <- combn(data$group, 2)

    for (i in 1:ncol(combinations)) {
    group1 <- combinations[1, i]
    group2 <- combinations[2, i]
    
    mean1 <- data$mean[data$group == group1]
    mean2 <- data$mean[data$group == group2]
    
    fc1 <- mean1 / mean2
    fc2 <- mean2 / mean1
    
    # Store both comparisons
    results <- rbind(results, data.frame(group1 = group1, group2 = group2, fold_change = fc1))
    results <- rbind(results, data.frame(group1 = group2, group2 = group1, fold_change = fc2))
    
    
    }
  
  results$avg_logFC <- log2(results$fold_change)
  
  return(results)
}


#' Perform multi-group statistical testing
#'
#' The \code{test_multi_groups} function performs a statistical comparison of multiple groups
#' using either parametric or non-parametric methods, depending on the data and user preference.
#' It includes Levene's test for variance homogeneity, followed by one of: ANOVA, Welch's ANOVA,
#' or Kruskal-Wallis test. Post-hoc pairwise comparisons are also performed without p-value adjustment.
#'
#' @param df A data frame containing the input data.
#' @param value_column A string specifying the name of the numeric column with values to test.
#' @param grouping_column A string specifying the name of the column with grouping information.
#' @param parametric Logical. If \code{TRUE}, parametric tests (ANOVA or Welch’s ANOVA) are used.
#' If \code{FALSE}, the Kruskal-Wallis test is used.
#' @param paired Logical. Indicates whether post-hoc comparisons should be paired (relevant for non-parametric tests).
#' @param adjustment.method A string indicating the method for p-value adjustment (currently stored in the output
#' but not applied during post-hoc tests; default is \code{"bonferroni"}).
#'
#' @return An S4 object of class \code{statistic} with the following slots:
#' \itemize{
#'   \item \code{test} — name of the main statistical test used (\code{"ANOVA"}, \code{"Welch's ANOVA"}, or \code{"Kruskal-Wallis"}),
#'   \item \code{leven_var_test} — a list with results from Levene’s test and interpretation,
#'   \item \code{posthoc_test} — name of the post-hoc test performed,
#'   \item \code{test_data} — a list containing the main test statistic, degrees of freedom, and p-value,
#'   \item \code{posthoc_data} — a list containing raw p-values and group comparisons from post-hoc tests.
#' }
#'
#' @details
#' - If \code{parametric = FALSE}, Kruskal-Wallis and pairwise Wilcoxon tests are used.  
#' - If \code{parametric = TRUE}, Levene’s test is used to choose between:
#'   \itemize{
#'     \item \strong{Equal variances} → classical ANOVA and pairwise t-tests,
#'     \item \strong{Unequal variances} → Welch’s ANOVA and Welch t-tests.
#'   }
#' - P-values in post-hoc tests are returned unadjusted but can be manually corrected later using \code{adjustment.method}.
#'
#' This function is primarily used internally by \code{\link{multi_groups_analysis}} to encapsulate statistical logic.
#'
#' @import stats
#' @importFrom car leveneTest
#' @importFrom methods new setClass
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 10),
#'   value = c(rnorm(10, mean = 5), rnorm(10, mean = 6), rnorm(10, mean = 7))
#' )
#' test_multi_groups(df, value_column = "value", grouping_column = "group", parametric = TRUE)
test_multi_groups <- function(df, 
                              value_column, 
                              grouping_column, 
                              parametric = TRUE, 
                              paired = FALSE, 
                              adjustment.method = 'bonferroni') {
  
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
      
      info <- "Levene's p < 0.05: variance not equal"
      
    } else {
      
      info <- "Levene's p > 0.05: variance equal"
      
    }
    
    formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
    
    kruskal_result <- kruskal.test(formula, data = df)
    
    kruskal_list <- list(
      statistic = kruskal_result$statistic,
      df = kruskal_result$parameter,
      p.value = kruskal_result$p.value
    )
    
    
    
   
    pairwise_results <- pairwise.wilcox.test(df[[value_column]], df[[grouping_column]], p.adjust.method = adjustment.method, paired = FALSE)
      
    
    
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
      
      info <- "Levene's p < 0.05: variance not equal"
      
      formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
      
      welch_anova <- oneway.test(formula, data = df, var.equal = FALSE)
      
      
      aov_results <- list(
        statistic =  welch_anova$statistic,
        df = welch_anova$parameter,
        p.value = welch_anova$p.value
        
      )
      
      
      pairwise_results <- pairwise.t.test(df[[value_column]], df[[grouping_column]], p.adjust.method = adjustment.method, paired = FALSE, pool.sd = FALSE)
      
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
                       posthoc_test = "Welch's t-test",
                       test_data = as.list(aov_results),
                       posthoc_data = posthc_results)
    
    } else {
      
      info <- "Levene's p > 0.05: variance equal"
      formula <- as.formula(paste0(sym(value_column), " ~ ", sym(grouping_column)))
      
      aov <- aov(formula, data = df)
      
      aov_results = summary(aov)
      aov_results <- list(
        statistic = aov_results[[1]][["F value"]][1],
        df = aov_results[[1]][["Df"]][1],
        p.value = aov_results[[1]][["Pr(>F)"]][1]
      )
      
      
      pairwise_results <- pairwise.t.test(df[[value_column]], df[[grouping_column]], p.adjust.method = adjustment.method, paired = FALSE, pool.sd = FALSE)
      
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
                       posthoc_test = 't-test',
                       test_data = as.list(aov_results),
                       posthoc_data = posthc_results)
      
      
    }
    

    
    
  }
  
  return(statistic)
  
  
}
  



#' Perform statistical testing between two groups
#'
#' The \code{test_two_groups} function performs a statistical comparison between two groups
#' using either a parametric (t-test) or non-parametric (Wilcoxon or Mann-Whitney) test.
#' It returns a custom S4 object with the test result details.
#'
#' @param df A data frame containing the data for analysis.
#' @param value_column A string giving the name of the numeric column with values to compare.
#' @param grouping_column A string specifying the name of the column with group labels (must contain exactly two groups).
#' @param parametric Logical. If \code{TRUE}, a t-test is used; otherwise, a Wilcoxon or Mann-Whitney test is used.
#' @param paired Logical. Whether a paired test should be used. Applies to both parametric and non-parametric tests.
#'
#' @return An S4 object of class \code{statistic}, containing:
#' \itemize{
#'   \item \code{test} — the name of the statistical test performed,
#'   \item \code{p.val} — the p-value from the test,
#'   \item \code{statistic} — the test statistic value,
#'   \item \code{paired} — whether the test was paired.
#' }
#'
#' @details
#' - If \code{parametric = TRUE}, a two-sample or paired t-test is performed.
#' - If \code{parametric = FALSE}, the function uses either the Wilcoxon Signed-Rank Test (for paired data)
#' or the Mann-Whitney U Test (for independent groups).
#'
#' This is a **supporting function** used internally by \code{\link{two_groups_analysis}}.
#'
#' @import stats
#' @importFrom methods new setClass
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B"), each = 10),
#'   value = c(rnorm(10, mean = 5), rnorm(10, mean = 6))
#' )
#' test_two_groups(df, value_column = "value", grouping_column = "group", parametric = TRUE)
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
    
   
    

    wcox <- wilcox.test(df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[1]], df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[2]], alternative = 'two.sided', paired)
    
   
    
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
    
    
    tt <- t.test(df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[1]], df[[value_column]][df[[grouping_column]] %in% unique(df[[grouping_column]])[2]], alternative = 'two.sided', paired, var.equal = TRUE)
    
    
    
    
    statistic <- new("statistic", 
                     test = "t-test",
                     p.val = tt$p.value,
                     statistic = tt$statistic,
                     paired = paired
    )
    
    
  }
  
  return(statistic)
  
}



#' Two-group statistical analysis and plotting
#'
#' This function performs statistical analysis and generates summary plots (violin, bar, and box plots)
#' for comparisons between two groups. It supports both parametric (t-test) and non-parametric (Wilcoxon/Mann-Whitney) tests,
#' and can handle paired data.
#'
#' @param value_column A string. Name of the column containing the numeric values to be analyzed.
#' @param grouping_column A string. Name of the column indicating group membership (must contain exactly two unique values).
#' @param data A data frame containing the data.
#' @param bar_queue Optional character vector. Specifies the desired order of the groups on the x-axis. Must match group names.
#' @param x_label A string for the x-axis label. Default is an empty string.
#' @param x_angle Numeric. Rotation angle of x-axis labels (in degrees). Default is 30.
#' @param y_label A string for the y-axis label. Default is an empty string.
#' @param size Numeric. Font size for axis labels. Default is 10.
#' @param bar_size Numeric. Width of bars in bar and box plots. Default is 0.5.
#' @param parametric Logical. If \code{TRUE}, a t-test is performed; if \code{FALSE}, a Wilcoxon/Mann-Whitney test is used.
#' @param paired Logical. Whether to perform a paired test. Applies to both parametric and non-parametric tests.
#' @param bars Character. Type of error bars: either \code{"sem"} or \code{"sd"} for standard error or standard deviation. Default is \code{"sem"}.
#' @param bars_size Numeric. Size (thickness) of error bars. Default is 1.
#' @param stat_plot_ratio Numeric between 0 and 1. Proportion of plot height reserved for displaying significance annotations. Default is 0.15.
#' @param y_break Optional numeric. Interval for y-axis tick marks. Default is \code{NaN} (automatic).
#' @param brew_colors Character. RColorBrewer palette name. Default is \code{"Dark2"}.
#'
#' @return An S4 object of class \code{two_groups_analysis}, containing:
#' \itemize{
#'   \item \code{violin_plot} — Violin plot with overlaid statistics.
#'   \item \code{bar_plot} — Bar plot with error bars.
#'   \item \code{box_plot} — Box plot with optional statistical annotations.
#'   \item \code{statistic_tests} — An S4 object from \code{test_two_groups()}, with test details.
#'   \item \code{statistic_data} — A summary data frame with group-wise statistics.
#'   \item \code{statistic_txt_resum} — Character summary of the statistical test result.
#'   \item \code{avg_FC_results} — Output of \code{avg_FC()}, representing average fold change results.
#' }
#'
#' @details
#' This function is designed for exploratory visualization and analysis of two-group comparisons.
#' If more than two groups are present in the grouping column, an error is thrown and the user is advised to use \code{\link{multi_groups_analysis}}.
#'
#'
#' @import ggplot2
#' @import patchwork
#' @import ggsignif
#' @import stats
#' @importFrom methods new setClass
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = rep(c("Control", "Treatment"), each = 10),
#'   value = c(rnorm(10, 5), rnorm(10, 7))
#' )
#' result <- two_groups_analysis(value_column = "value", grouping_column = "group", data = df)
#' result@bar_plot  # View bar plot
#' result@statistic_txt_resum  # Summary of the statistical test
two_groups_analysis <- function(value_column, grouping_column, data, bar_queue = NaN, x_label = '', x_angle = 30, y_label = '', size = 10, bar_size = 0.5, parametric = FALSE, paired = FALSE, bars = 'sem', bars_size = 1, stat_plot_ratio = 0.15, y_break = NaN, brew_colors = 'Dark2') {
  
  
  if (length(unique(data[[grouping_column]])) == 2) {
    
    
    plot_df = get_stats(data, value_column, grouping_column)
    
    if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & length(bar_queue) == length(plot_df[[grouping_column]]) & identical(sort(bar_queue), sort(plot_df[[grouping_column]]))) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = bar_queue)
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = bar_queue)
      
    } else if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & length(bar_queue) != length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])
      
      
      print('Warning! The `bar_queue` length in not equal with number of groups!')
      
    } else if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & !identical(sort(bar_queue), sort(plot_df[[grouping_column]]))) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])
      
      
      print('Warning! The `bar_queue` vaqlue is not included in groups!')
      
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
    
    sigs <- c(sig)
    
    
    
    
    
    ################################################################################
    
    
    MinMeanSEMMax <- function(x) {
      v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
      names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
      v
    }
    
    MinMeanSDMax <- function(x) {
      v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
      names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
      return(v)
    }
    
    
    if (bars == 'sd') {
      
      violin_plot <- ggplot(data, aes(x = !!sym(grouping_column), y = !!sym(value_column), fill = !!sym(grouping_column))) +
        geom_violin(trim = FALSE, show.legend = FALSE, color = "black") +  # Plot the distribution
        geom_point() +
        stat_summary(fun.data=MinMeanSDMax, geom="boxplot",width = bar_size *0.3, color = "black", size = 0.5) +
        theme_minimal()
      
    } else {
      
      violin_plot <- ggplot(data, aes(x = !!sym(grouping_column), y = !!sym(value_column), fill = !!sym(grouping_column))) +
        geom_violin(trim = FALSE, show.legend = FALSE, color = "black") +  # Plot the distribution
        geom_point() +
        stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", width = bar_size *0.3, color = "black", size = 0.5) +
        theme_minimal()
    }
    
    
    
    if (!is.na(y_break)) {
      
      violin_plot <- violin_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }
    
    
    
    
    if (bars == 'sd') {
      
      bar_plot = ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = bar_size - 0.1, size = 0.5)
      
    } else {
      
      bar_plot <- ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-SEM, ymax = mean+SEM), width = bar_size - 0.1, size = 0.5)
      
    }
    
    
    
    
    
    if (!is.na(y_break)) {
      
      bar_plot <- bar_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }
    
    
    
    if (bars == 'sd') {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_point() +
        stat_summary(fun.data=MinMeanSDMax, geom="boxplot",width = bar_size, color = "black", size = 0.5) +
        theme_minimal()
      
    } else {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_point() +
        stat_summary(fun.data=MinMeanSEMMax, geom="boxplot",width = bar_size, color = "black", size = 0.5) +
        theme_minimal()
      
      
    }
    
    
    
    
    
    if (!is.na(y_break)) {
      
      box_plot <- box_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }
    
    
    
    
    
  
      
      
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
      coord_cartesian(ylim = c(0, round(max(y_pos)))) +
      annotate("text", x = -Inf, y = Inf, label =  paste(' ', results@test, 'p =', results@p.val), 
               hjust = 0, vjust = 1.25, size = 2.8) +
      theme_void()
      
      
      
      
 
    
    
    bar_plot = bar_plot + ylab(y_label)
    bar_plot = bar_plot + xlab(x_label)
    
    box_plot = box_plot + ylab(y_label)
    box_plot = box_plot + xlab(x_label)
    
    violin_plot = violin_plot + ylab(y_label)
    violin_plot = violin_plot + xlab(x_label)
    
    
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
    
    
    violin_plot = violin_plot +  theme_classic() +
      theme(axis.title.y = element_text(size = size),
            axis.title.x = element_text(size = size)) +
      scale_x_discrete(guide = guide_axis(angle = x_angle)) +
      scale_fill_brewer(palette=brew_colors) +
      theme(legend.position="none")  
    
    
    
    if (length(list_of_comparison) > 0) {
      box_plot <- signif_plot + box_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      bar_plot <- signif_plot + bar_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      violin_plot <- signif_plot + violin_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      
    }
    
    
    results_text <- paste0('Group test: ', results@test,"\n")
    results_text <- paste0(results_text,'p-val: ',results@p.val,"\n")
    results_text <- paste0(results_text,'statistic: ', results@statistic,"\n")
    



    
    setClass(
      "two_groups_analysis",
      representation(
        violin_plot = "ANY",
        bar_plot = "ANY",
        box_plot = 'ANY',
        statistic_tests = "ANY",
        statistic_data = "list",
        statistic_txt_resum = 'ANY',
        avg_FC_results = 'ANY'
        
      )
    )
    
    results <- new("two_groups_analysis",
                   violin_plot = violin_plot,
                   bar_plot = bar_plot,
                   box_plot = box_plot,
                   statistic_tests = results,
                   statistic_data = plot_df,
                   statistic_txt_resum = results_text,
                   avg_FC_results = avg_FC(plot_df))
    
    return(results)
  
  } else if (length(unique(data[[grouping_column]])) > 2) {
    
    stop("The number of groups in the analysis is greater than 2.\n   For more than two groups use the multi_groups_analysis() function")
    
  } else {
    
    stop("The number of groups in the analysis is wrong. Check grouping_column")
    
    
  }
  
}





#' Multiple Group Statistical Analysis and Visualization
#'
#' Performs statistical analysis for more than two groups, including parametric or non-parametric tests, post-hoc comparisons, and generates bar, box, and violin plots with error bars and significance annotations.
#'
#' @param value_column Character. Name of the column in `data` representing numeric values to be compared.
#' @param grouping_column Character. Name of the column in `data` representing grouping factor.
#' @param data Data frame. The dataset containing `value_column` and `grouping_column`.
#' @param bar_queue Optional character vector. Custom order of groups for plotting (must match group names).
#' @param x_label Character. Label for the x-axis.
#' @param x_angle Numeric. Angle for x-axis text labels (default = 30).
#' @param y_label Character. Label for the y-axis.
#' @param size Numeric. Font size for axis titles (default = 10).
#' @param bar_size Numeric. Width of the bars and boxes (default = 0.5).
#' @param parametric Logical. If TRUE, uses parametric tests (ANOVA + t-test); if FALSE, uses Kruskal-Wallis + Wilcoxon / Mann-Whitney test.
#' @param paired Logical. If TRUE, assumes paired samples (only for parametric).
#' @param include_ns Logical. If FALSE, excludes non-significant comparisons from plot.
#' @param bars Character. Error bars type: `"sem"` (default) or `"sd"`.
#' @param bars_size Numeric. Thickness of error bars (not currently used).
#' @param adjustment.method Character. Method for p-value adjustment in post-hoc test (e.g., `"bonferroni"`, `"holm"`).
#' @param stat_plot_ratio Numeric. Relative height of significance plot above main plot (default = 0.2).
#' @param stat_hight Numeric. Controls spacing between significance brackets (default = 10).
#' @param y_break Numeric. Break interval for y-axis. If NA, default ggplot scale is used.
#' @param brew_colors Character. RColorBrewer palette name (default = `"Dark2"`).
#'
#' @return An S4 object of class `multi_groups_analysis` with the following slots:
#' \describe{
#'   \item{violin_plot}{ggplot object — violin plot with error bars.}
#'   \item{bar_plot}{ggplot object — bar plot with error bars.}
#'   \item{box_plot}{ggplot object — box plot with error bars.}
#'   \item{statistic_tests}{S4 object — result of group tests (ANOVA/KW + post-hoc).}
#'   \item{statistic_data}{data.frame — summary statistics (mean, SEM, SD) for each group.}
#'   \item{statistic_txt_resum}{Character — textual summary of statistical results.}
#'   \item{avg_FC_results}{data.frame — fold-change comparisons between groups.}
#' }
#'
#' @details
#' This function is designed for exploratory data analysis involving multiple groups. It automatically handles appropriate statistical tests based on the `parametric` and `paired` flags, formats plots, and annotates p-values with asterisks.
#'
#' If only 2 groups are detected, an error is returned recommending the use of `two_groups_analysis()`.
#'
#' @import ggplot2
#' @import dplyr
#' @import stats
#' @import ggsignif
#' @import patchwork
#' @importFrom RColorBrewer brewer.pal
#' 
#' @examples
#' result <- multi_groups_analysis(
#'   value_column = "Score",
#'   grouping_column = "Group",
#'   data = my_data,
#'   parametric = TRUE,
#'   bar_queue = c("Control", "Treatment1", "Treatment2"),
#'   x_label = "Groups",
#'   y_label = "Response",
#'   include_ns = FALSE
#' )
#' result@bar_plot
#' result@statistic_txt_resum
#'
#' @export
multi_groups_analysis <- function(value_column, 
                                  grouping_column, 
                                  data, 
                                  bar_queue = NaN, 
                                  x_label = '',
                                  x_angle = 30, 
                                  y_label = '', 
                                  size = 10, 
                                  bar_size = 0.5, 
                                  parametric = FALSE, 
                                  paired = FALSE, 
                                  include_ns = FALSE, 
                                  bars = 'sem', 
                                  bars_size = 1, 
                                  adjustment.method = 'bonferroni', 
                                  stat_plot_ratio = 0.2, 
                                  stat_hight = 10,
                                  y_break = NaN, 
                                  brew_colors = 'Dark2') {


  if (length(unique(data[[grouping_column]])) > 2) {

    plot_df = get_stats(data, value_column, grouping_column)

    if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & length(bar_queue) == length(plot_df[[grouping_column]]) & identical(sort(bar_queue), sort(plot_df[[grouping_column]]))) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = bar_queue)
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = bar_queue)

    } else if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & length(bar_queue) != length(plot_df[[grouping_column]])) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])


      print('Warning! The `bar_queue` length in not equal with number of groups!')

    } else if (!TRUE %in% unique(is.na(bar_queue)) & is.vector(bar_queue) & !identical(sort(bar_queue), sort(plot_df[[grouping_column]]))) {
      plot_df[[grouping_column]] <- factor(plot_df[[grouping_column]], levels = plot_df[[grouping_column]])
      data[[grouping_column]] <- factor(data[[grouping_column]], levels = plot_df[[grouping_column]])
      
      
      print('Warning! The `bar_queue` vaqlue is not included in groups!')
      
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
    
    
    MinMeanSEMMax <- function(x) {
      v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
      names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
      v
    }
    
    MinMeanSDMax <- function(x) {
      v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
      names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
      return(v)
    }
    
    
    if (bars == 'sd') {
      
      violin_plot <- ggplot(data, aes(x = !!sym(grouping_column), y = !!sym(value_column), fill = !!sym(grouping_column))) +
        geom_violin(trim = FALSE, show.legend = FALSE, color = "black") +  # Plot the distribution
        geom_point() +
        stat_summary(fun.data=MinMeanSDMax, geom="boxplot",width = bar_size *0.3, color = "black", size = 0.5) +
        theme_minimal()
      
    } else {
      
      violin_plot <- ggplot(data, aes(x = !!sym(grouping_column), y = !!sym(value_column), fill = !!sym(grouping_column))) +
        geom_violin(trim = FALSE, show.legend = FALSE, color = "black") +  # Plot the distribution
        geom_point() +
        stat_summary(fun.data=MinMeanSEMMax, geom="boxplot", width = bar_size *0.3, color = "black", size = 0.5) +
        theme_minimal()
    }
    
    
    
    if (!is.na(y_break)) {
      
      violin_plot <- violin_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }
    
    
    
    
    if (bars == 'sd') {
      
      bar_plot = ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = bar_size - 0.1, size = 0.5)
      
    } else {
      
      bar_plot <- ggplot(plot_df, aes(x = !!sym(grouping_column), y = mean, fill = !!sym(grouping_column)))+
        geom_bar(stat = "identity", show.legend = FALSE, width = bar_size, color = "black") +
        geom_errorbar(aes(ymin = mean-SEM, ymax = mean+SEM), width = bar_size - 0.1, size = 0.5)
      
    }
    
    



    if (!is.na(y_break)) {

      bar_plot <- bar_plot + scale_y_continuous(breaks = seq(0, max(data[[value_column]], na.rm = TRUE), by = y_break))
    }


    
    if (bars == 'sd') {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_point() +
        stat_summary(fun.data=MinMeanSDMax, geom="boxplot",width = bar_size, color = "black", size = 0.5) +
        theme_minimal()
      
    } else {
      
      box_plot <- ggplot(data, aes(y = !!sym(value_column), x = !!sym(grouping_column), fill = !!sym(grouping_column))) +
        geom_point() +
        stat_summary(fun.data=MinMeanSEMMax, geom="boxplot",width = bar_size, color = "black", size = 0.5) +
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
                  y_position = y_pos*(1+stat_hight),
                  map_signif_level = FALSE, textsize = 4) +
      coord_cartesian(ylim = c(0, round(max(y_pos) + fc + 5))) +
      annotate("text", x = -Inf, y = Inf, label = paste(' ', results@leven_var_test$response,
                                                        ' | ' ,results@test, 'p =', 
                                                        results@test_data$p.value, ' | ', 
                                                        ' post-hoc:',results@posthoc_test, ' | ', 
                                                        ' p.adj:', results@posthoc_data$adjustment), 
                                                        hjust = 0, vjust = 1.25, size = 2.8) +
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
        coord_cartesian(ylim = c(0, max(max_y) * 1.08)) +
        annotate("text", x = -Inf, y = Inf, label = paste(' ',results@leven_var_test$response,
                                                          ' | ' ,results@test, 'p =', 
                                                          results@test_data$p.value), 
                 hjust = 0, vjust = 1.25, size = 2.8) 



      max_y <- max(data[[value_column]])
      min_y <- min(data[[value_column]])





      box_plot = box_plot +
        coord_cartesian(ylim = c(min_y, max(max_y)* 1.08)) +
        annotate("text", x = -Inf, y = Inf, label = paste(' ', results@leven_var_test$response,
                                                          ' | ' ,results@test, 'p =', 
                                                          results@test_data$p.value), 
                 hjust = 0, vjust = 1.25, size = 2.8) 
      
      
      violin_plot = violin_plot +
        # coord_cartesian(ylim = c(min_y, max(max_y)* 1.08)) +
        annotate("text", x = -Inf, y = Inf, label = paste(' ', results@leven_var_test$response,
                                                          ' | ' ,results@test, 'p =', 
                                                          results@test_data$p.value), 
                 hjust = 0, vjust = 1.20, size = 2.8) 
      



    }

    
    bar_plot = bar_plot + ylab(y_label)
    bar_plot = bar_plot + xlab(x_label)

    box_plot = box_plot + ylab(y_label)
    box_plot = box_plot + xlab(x_label)
    
    violin_plot = violin_plot + ylab(y_label)
    violin_plot = violin_plot + xlab(x_label)


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
    
    
    violin_plot = violin_plot +  theme_classic() +
      theme(axis.title.y = element_text(size = size),
            axis.title.x = element_text(size = size)) +
      scale_x_discrete(guide = guide_axis(angle = x_angle)) +
      scale_fill_brewer(palette=brew_colors) +
      theme(legend.position="none")  
    
    

    if (length(list_of_comparison) > 0) {
      box_plot <- signif_plot + box_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      bar_plot <- signif_plot + bar_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      violin_plot <- signif_plot + violin_plot  + plot_layout(ncol = 1, heights = c(10*stat_plot_ratio, 10*(1-stat_plot_ratio)))
      
    }

    
    results_text <- paste0('Group test: ', results@test,"\n")
    results_text <- paste0(results_text,'Post-hoc test: ', results@posthoc_test,"\n")
    results_text <- paste0(results_text,'Post-hoc p-val adjustment: ', results@posthoc_data$adjustment,"\n")
    results_text <- paste0(results_text, 'Pair1   |   Pair2   |   p-val   ',"\n")
    
    
    for (p in 1:length(results@posthoc_data$p.adjusted)) {
      results_text <- paste0(results_text, ' ', results@posthoc_data$pair1[p], ' ', results@posthoc_data$pair2[p], ' ', results@posthoc_data$p.adjusted[p] ,"\n")
      
    }

    setClass(
      "multi_groups_analysis",
      representation(
        violin_plot = 'ANY',
        bar_plot = "ANY",
        box_plot = 'ANY',
        statistic_tests = "ANY",
        statistic_data = "list",
        statistic_txt_resum = 'ANY',
        avg_FC_results = 'ANY'

      )
    )

    results <- new("multi_groups_analysis",
                       violin_plot = violin_plot,
                       bar_plot = bar_plot,
                       box_plot = box_plot,
                       statistic_tests = results,
                       statistic_data = plot_df,
                       statistic_txt_resum = results_text,
                       avg_FC_results = avg_FC(plot_df)
                   
                      
    )
    
    gc()

    return(results)

  } else if (length(unique(data[[grouping_column]])) == 2) {

    stop("The number of groups in the analysis is equal to 2.\n   For two groups use the two_groups_analysis() function")

  } else {

    stop("The number of groups in the analysis is wrong. Check grouping_column")


  }

}




#' Multi-variable Groups Analysis with Statistical Tests and Plotting
#'
#' This function performs a statistical analysis (ANOVA or non-parametric alternatives) 
#' on a given dataset grouped by specified variables, runs post-hoc tests for each interval,
#' adjusts p-values, and generates a ggplot2 plot showing means with error bars and significance labels.
#' It also calculates average fold changes between groups.
#'
#' @param data A data.frame containing the dataset.
#' @param stat_col A string specifying the name of the column with the response variable (numeric).
#' @param interval_col A string specifying the name of the column representing intervals or time points (factor or character).
#' @param group_col A string specifying the name of the grouping variable column (factor or character).
#' @param parametric Logical; if TRUE, parametric tests (ANOVA, t-test) are used, otherwise non-parametric tests (ART-ANOVA, Wilcoxon/Kruskal-Wallis) are applied. Default is TRUE.
#' @param paired Logical; whether to perform paired tests (currently not implemented in detail). Default is FALSE.
#' @param adj Character or NA; method for p-value adjustment in post-hoc tests ("bf" for Bonferroni, "bh" for Benjamini-Hochberg, or NA for none). Default is NA.
#' @param error Character; type of error bars to show on the plot: 'sem' (standard error of mean) or 'sd' (standard deviation). Default is 'sem'.
#' @param tx_pos Numeric; position adjustment factor for placing significance text above error bars. Default is 0.02.
#'
#' @return An S4 object of class \code{multi_var_groups_analysis} containing:
#' \item{plot}{A ggplot2 object of the plotted means with error bars and significance annotations.}
#' \item{statistic_group}{Summary statistics of the main group test (ANOVA or ART-ANOVA).}
#' \item{test_name}{Name of the main statistical test used (character).}
#' \item{statistic_post_hoc}{Data frame of post-hoc test results for each interval.}
#' \item{stats}{Summary statistics of the data by group and interval.}
#' \item{avg_FC_results}{Data frame with average fold change results between groups for each interval.}
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import ggpubr
#' @import purrr
#' @import stats
#' @import ARTool
#' @import rlang
#' 
#' @examples
#' # Example usage
#' results <- multi_var_groups_analysis(data = my_data,
#'                                     stat_col = "measurement",
#'                                     interval_col = "timepoint",
#'                                     group_col = "treatment",
#'                                     parametric = TRUE,
#'                                     adj = "bh",
#'                                     error = "sem")
#' print(results@plot)
#'
#' @export
multi_var_groups_analysis <- function(data,
                                      stat_col, 
                                      interval_col, 
                                      group_col,
                                      parametric = TRUE,
                                      paired = FALSE,
                                      adj = NA,
                                      error = 'sem',
                                      tx_pos = 0.02
                                      
                                      
) {
  
  
  setClass(
    "multi_var_groups_analysis",
    representation(
      plot = 'ANY',
      statistic_group = 'ANY',
      test_name = 'character',
      statistic_post_hoc = 'ANY',
      stats = 'ANY',
      avg_FC_results = 'ANY'
      
    )
  )
  
  data[[interval_col]] <- factor(data[[interval_col]], levels = unique(data[[interval_col]]))
  data[[group_col]] <- as.factor(data[[group_col]])
  
  # anova ~ groups
  formula_text <- paste0("`", stat_col, "` ~ `", interval_col, "` * `", group_col, "`")
  formula <- as.formula(formula_text)  
  
  if (parametric) {
    
    aov_result <- aov(formula, data = data)
    summary_groups <- summary(aov_result)
    p_value_1 <- summary_groups[[1]]$`Pr(>F)`[1]
    p_nam_1 <- trimws(rownames(summary_groups[[1]])[1])
    p_value_2 <- summary_groups[[1]]$`Pr(>F)`[2]
    p_nam_2 <- trimws(rownames(summary_groups[[1]])[2])
    p_value_3 <- summary_groups[[1]]$`Pr(>F)`[3]
    p_nam_3 <- trimws(rownames(summary_groups[[1]])[3])
    
    test_name = 'ANOVA'
    
    
  } else if (parametric == FALSE) {
    
    safe_anova <- tryCatch({
      summary_art <- art(formula, data = data)
      summary_groups <- anova(summary_art)
      p_value_1 <- summary_groups$`Pr(>F)`[1]
      p_nam_1 <- summary_groups[[1]][1]
      p_value_2 <- summary_groups$`Pr(>F)`[2]
      p_nam_2 <- summary_groups[[1]][2]
      p_value_3 <- summary_groups$`Pr(>F)`[3]
      p_nam_3 <- summary_groups[[1]][3]
      
      statistic <- summary_groups$F[1]
      list(p_nam_1 = p_nam_1, 
           p_value_1 = p_value_1,
           p_nam_2 = p_nam_2,
           p_value_2 = p_value_2,
           p_nam_3 = p_nam_3,
           p_value_3 = p_value_3,
           statistic = statistic)
    }, error = function(e) {
      list(p_nam_1 = interval_col, 
           p_value_1 = 1,
           p_nam_2 = group_col,
           p_value_2 = 1,
           p_nam_3 = paste0(interval_col, ':', group_col),
           p_value_3 = 1,
           statistic = NaN)
    })
    
    summary_groups <- safe_anova
    
    p_value_1 <- summary_groups$p_value_1
    p_nam_1 <- summary_groups$p_nam_1
    p_value_2 <- summary_groups$p_value_2
    p_nam_2 <- summary_groups$p_nam_2
    p_value_3 <- summary_groups$p_value_3
    p_nam_3 <- summary_groups$p_nam_3
   
    # Aligned Rank Transform
    test_name = 'ART-ANOVA'
    
    
  }
  
  
    p_value_1 <- case_when(
      is.na(p_value_1)    ~ "ns",
      p_value_1 < 0.001   ~ "***",
      p_value_1 < 0.01    ~ "**",
      p_value_1 < 0.05    ~ "*",
      TRUE                ~ "ns"
    )
    
    
    p_value_2 <- case_when(
      is.na(p_value_2)    ~ "ns",
      p_value_2 < 0.001   ~ "***",
      p_value_2 < 0.01    ~ "**",
      p_value_2 < 0.05    ~ "*",
      TRUE                ~ "ns"
    )
    
    
    p_value_3 <- case_when(
      is.na(p_value_3)    ~ "ns",
      p_value_3 < 0.001   ~ "***",
      p_value_3 < 0.01    ~ "**",
      p_value_3 < 0.05    ~ "*",
      TRUE                ~ "ns"
    )
    
  
  
  
  # intervals
  
    
    run_test_for_interval <- function(data_slice, interval_value, parametric, adj, paired) {
      n_groups <- data_slice %>% pull(!!sym(group_col)) %>% unique() %>% length()
      
      test_type <- NA
      test_stat <- NA
      p_value <- NA
      
      group_vals <- unique(data_slice[[group_col]])
      
      # parametric test
      if (parametric) {
        if (n_groups == 2) {
          test_result <- tryCatch({
            if (paired) {
              
              # paired data
              x <- data_slice %>% filter(!!sym(group_col) == group_vals[1]) %>% pull(!!sym(stat_col))
              y <- data_slice %>% filter(!!sym(group_col) == group_vals[2]) %>% pull(!!sym(stat_col))
              
              if (length(x) != length(y)) stop("Lengths of x and y must be equal for a paired test.")
              
              t.test(x, y, alternative = "two.sided", paired = TRUE)
              
            } else {
              
              x <- data_slice %>% filter(!!sym(group_col) == group_vals[1]) %>% pull(!!sym(stat_col))
              y <- data_slice %>% filter(!!sym(group_col) == group_vals[2]) %>% pull(!!sym(stat_col))
              
              t.test(x, y, alternative = "two.sided", paired = FALSE, var.equal = TRUE)
            }
          }, error = function(e) return(NULL))
          
          test_type <- if (paired) "Paired t-test" else "t-test"
          if (!is.null(test_result)) {
            test_stat <- test_result$statistic
            p_value <- test_result$p.value
          } else {
            test_stat <- NaN
            p_value <- 1
          }
          
        } else if (n_groups > 2) {
          formula <- reformulate(group_col, response = stat_col)
          aov_result <- tryCatch(
            aov(formula, data = data_slice),
            error = function(e) return(NULL)
          )
          test_type <- "ANOVA"
          if (!is.null(aov_result)) {
            aov_summary <- summary(aov_result)
            test_stat <- aov_summary[[1]]$`F value`[1]
            p_value <- aov_summary[[1]]$`Pr(>F)`[1]
          } else {
            test_stat <- NaN
            p_value <- 1
          }
        }
        
      } else {
        if (n_groups == 2) {
          test_result <- tryCatch({
            if (paired) {
              x <- data_slice %>% filter(!!sym(group_col) == group_vals[1]) %>% pull(!!sym(stat_col))
              y <- data_slice %>% filter(!!sym(group_col) == group_vals[2]) %>% pull(!!sym(stat_col))
              
              if (length(x) != length(y)) stop("Lengths of x and y must be equal for a paired test.")
              
              wilcox.test(x, y, alternative = "two.sided", paired = TRUE)
            } else {
              
              x <- data_slice %>% filter(!!sym(group_col) == group_vals[1]) %>% pull(!!sym(stat_col))
              y <- data_slice %>% filter(!!sym(group_col) == group_vals[2]) %>% pull(!!sym(stat_col))
              
              wilcox.test(x, y, alternative = "two.sided", paired = FALSE)
            }
          }, error = function(e) return(NULL))
          
          test_type <- if (paired) "Wilcoxon signed-rank" else "Mann-Whitney U"
          if (!is.null(test_result)) {
            test_stat <- test_result$statistic
            p_value <- test_result$p.value
          } else {
            test_stat <- NaN
            p_value <- 1
          }
          
        } else if (n_groups > 2) {
          formula <- reformulate(group_col, response = stat_col)
          kruskal_result <- tryCatch(
            kruskal.test(formula, data = data_slice),
            error = function(e) return(NULL)
          )
          test_type <- "Kruskal-Wallis"
          if (!is.null(kruskal_result)) {
            test_stat <- kruskal_result$statistic
            p_value <- kruskal_result$p.value
          } else {
            test_stat <- NaN
            p_value <- 1
          }
        }
      }
      
      return(tibble(
        interval = interval_value,
        test_type = test_type,
        test_stat = test_stat,
        p_value_raw = p_value,
        p_adjust_method = ifelse(is.na(adj), "none", adj)
      ))
    }
  
  
  
  
  result_df <- data %>%
    group_split(!!sym(interval_col)) %>%
    map_dfr(~ run_test_for_interval(.x, unique(.x[[interval_col]]), parametric = parametric, adj = adj, paired = paired))
  
  result_df$p_value_raw[is.nan(result_df$p_value_raw)] = 1
  
  if (!all(is.na(result_df$p_value_raw))) {
    result_df <- result_df %>%
      mutate(
        p_value = case_when(
          tolower(p_adjust_method) == "bf" ~ p.adjust(p_value_raw, method = "bonferroni"),
          tolower(p_adjust_method) == "bh" ~ p.adjust(p_value_raw, method = "BH"),
          TRUE                   ~ p_value_raw
        ),
        signif_label = case_when(
          is.na(p_value)        ~ "ns",
          p_value < 0.001       ~ "***",
          p_value < 0.01        ~ "**",
          p_value < 0.05        ~ "*",
          TRUE                  ~ "ns"
        )
      )
  }
  
  
  #######################################################################
  
  summary_data <- data %>%
    group_by(!!sym(group_col), !!sym(interval_col)) %>%
    summarise(
      mean = mean(!!sym(stat_col), na.rm = TRUE),
      sd   = sd(!!sym(stat_col), na.rm = TRUE),
      sem  = sd / sqrt(n()),
      .groups = 'drop'
    )
  
  
  result_df2 <- result_df %>%
    rename(!!interval_col := interval) %>%  
    select(!!sym(interval_col), signif_label)
  
  summary_data_signif <- summary_data %>%
    left_join(result_df2, by = interval_col)
  
  summary_data_signif$signif_label[summary_data_signif$signif_label %in% 'ns'] = ''
  
  
  
  if (tolower(error) == 'sd') {
    
    max_points <- max(summary_data_signif$mean + summary_data_signif$sd)*(1+(tx_pos/2))
    
    plot <- ggplot(summary_data_signif, aes(x = !!sym(interval_col), y = mean, color = !!sym(group_col), group = !!sym(group_col))) +
      geom_point(position = position_dodge(width = 0.3), size = 3) +
      geom_line(position = position_dodge(width = 0.3), linewidth = 1) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                    width = 0.2, position = position_dodge(width = 0.3)) +
      theme_minimal() +
      labs(
        x = interval_col,
        y = stat_col
      ) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      scale_color_brewer(palette = "Set1") +
      
      geom_text(
        data = summary_data_signif,
        aes(x = !!sym(interval_col), y = max_points, label = signif_label),
        color = "black",
        size = 6,
        vjust = 0,
        inherit.aes = FALSE
      ) +
      
    
        annotate("text",
               x = -Inf,  
               y = max_points * (1+(tx_pos)), 
               label = paste0(' ', test_name, ': ', 
                               p_nam_1, ' p = ', p_value_1,
                              ';  ',p_nam_2, ' p = ', p_value_2,
                              ';  ',p_nam_3, ' p = ', p_value_3,'  |  ', 
                              'post-hoc: ', result_df$test_type[1], '  |  ', 
                              'p.adj: ', result_df$p_adjust_method[1]), 
               hjust = 0,
               vjust = 0,  
               size = 2.8)
    
    
    
  } else {
    
    max_points <- max(summary_data_signif$mean + summary_data_signif$sem)*(1+(tx_pos/2))
    
    plot <- ggplot(summary_data_signif, aes(x = !!sym(interval_col), y = mean, color = !!sym(group_col), group = !!sym(group_col))) +
      geom_point(position = position_dodge(width = 0.3), size = 3) +
      geom_line(position = position_dodge(width = 0.3), linewidth = 1) +
      geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                    width = 0.2, position = position_dodge(width = 0.3)) +
      theme_minimal() +
      labs(
        x = interval_col,
        y = stat_col
      ) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      scale_color_brewer(palette = "Set1") +
      
      geom_text(
        data = summary_data_signif,
        aes(x = !!sym(interval_col), y = max_points, label = signif_label),
        color = "black",
        size = 6,
        vjust = 0,
        inherit.aes = FALSE
      ) +
      
      annotate("text",
               x = -Inf,  
               y = max_points * (1+(tx_pos)), 
               label = paste0('  ', test_name, ': ', 
                               p_nam_1, ' p = ', p_value_1,
                              ';  ',p_nam_2, ' p = ', p_value_2,
                              ';  ',p_nam_3, ' p = ', p_value_3,'  |  ', 
                              'post-hoc: ', result_df$test_type[1], '  |  ', 
                              'p.adj: ', result_df$p_adjust_method[1]), 
               hjust = 0,
               vjust = 0,  
               size = 2.8)
    
    
    
    
    
  }
  
  
  
  
  
  ################################################################################
  
  results_FC <- data.frame(
    group1 = character(),
    group2 = character(),
    avg_fold_change = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  combinations <- combn(unique(summary_data[[group_col]]), 2)
  
  results <- data.frame()
  for (i in 1:ncol(combinations)) {
    group1 <- combinations[1, i]
    group2 <- combinations[2, i]
    
    
    s1 <- summary_data[summary_data[[group_col]] == group1,]
    s2 <- summary_data[summary_data[[group_col]] == group2,]
    
    fc1 <- s1$mean / s2$mean
    fc2 <- s2$mean / s1$mean
    
    # Store both comparisons
    results <- rbind(results, data.frame(group1 = group1, group2 = group2, fold_change = fc1, interval = s1[[interval_col]]))
    results <- rbind(results, data.frame(group1 = group2, group2 = group1, fold_change = fc2, interval = s2[[interval_col]]))
    
    
  }
  
  results$avg_logFC <- log2(results$fold_change)
  

  
  
  results <- new("multi_var_groups_analysis",
                 plot = plot,
                 statistic_group = summary_groups,
                 test_name = test_name,
                 statistic_post_hoc = result_df,
                 stats = summary_data,
                 avg_FC_results = results)
  
  
  return(results)
}


