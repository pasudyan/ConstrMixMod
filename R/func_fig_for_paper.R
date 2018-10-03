# ========================================================================================= #
# Collection of functions to for Visualization of Gibbs Sampling Results                    #
# ========================================================================================= #

library(cowplot)
library(plyr)
library(magrittr)
library(ggplot2)
library(reshape2)

# =============== Function to create Posterior Predictive Plot ============== #
#' Plot predictive density for Gaussian
#'
#' This function generates the posterior predictive density for the Gaussian case
#' @param test_pred predictive probability at the grid
#' @param thr maximum number of rejected proposals for each observation
#' @param mu_r mode of Gaussian
#' @param labs whether labels should be included in the figure
#' @return posterior predictive density plot
#' @export
#'
fig_pred_prob_gauss <- function(test_pred, thr, mu_r = 0, labs = T){
  # ========================================================================================= #
  # Function to create Posterior Predictive Plot for 1-dimensional data                       #
  # Input: Posterior predictive results, threshold, data type                                 #
  # Output: Plot of Posterior predictive results                                              #
  # ========================================================================================= #

  thr_l <- length(thr)

  # Grid points for predictions
  delta   <- 0.005
  pts     <- seq(0, 1, delta)

  df_pts     <- data.frame(pts = pts)
  df_pts_nm  <- c("pts")

  norm_c  <- function(fx, deltax){
    sum(deltax*fx)
  }

  lower <- 0
  upper <- 1

  for (ss in 1:thr_l){
    if (mu_r == 0){
      mean_pred <- test_pred[[ss]][[1]]$pred_pt
    }

    if (mu_r == 0.5){
      mean_pred <- test_pred[[ss]]$pred_pt
    }

    mp_norm   <- mean_pred/norm_c(mean_pred, delta)
    df_pts    <- cbind(df_pts, mp_norm)
    df_pts_nm <- c(df_pts_nm, paste0(thr[ss], sep=""))
  }

  sd_r <- 0.1
  db_sim_df    <- dnorm(pts, mu_r, sd_r)

  db_sim_norm <- db_sim_df
  db_sim_norm <- db_sim_df/norm_c(db_sim_df, delta)

  df_sim_norm <- cbind(df_pts, db_sim_norm)
  df_sim_norm <- data.frame(df_sim_norm)
  colnames(df_sim_norm) <- c(df_pts_nm, "Real Density")

  if (mu_r == 0){
    cut <- which(df_sim_norm[,1] == 0.80)
    df_sim_norm <- df_sim_norm[1:cut,]
  }

  df_sim_melt <- cbind(pts=rep(df_sim_norm[,1]), melt(df_sim_norm[,-1], variable.name = "thresh",
                                                      value.name = "value"))

  if (labs == T){

    pp <- ggplot(df_sim_melt, aes(x = pts, y = value, color = thresh)) + geom_line() +
      # ggtitle("Posterior Predictive Density") +
      xlab("Values") + ylab("Posterior Predictive Density") + labs(colour = "Threshold") +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            text = element_text(size = 12))

  } else {
    pp <- ggplot(df_sim_melt, aes(x = pts, y = value, color = thresh)) + geom_line() +
      # ggtitle("Posterior Predictive Density") +
      xlab("Values") + ylab("Posterior Predictive Density") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line = element_line(colour = "black"),
            legend.position = "none",
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            text = element_text(size = 12))
  }

  return(pp)
}

#' Plot predictive density for Beta distribution
#'
#' This function generates the posterior predictive density for the Beta distribution
#' @param test_pred predictive probability at the grid
#' @param thr maximum number of rejected proposals for each observation
#' @param labs whether labels should be included in the figure
#' @return posterior predictive density plot
#' @export
#'
fig_pred_prob_beta <- function(test_pred, thr, labs = T){
  # ========================================================================================= #
  # Function to create Posterior Predictive Plot for 1-dimensional data                       #
  # Input: Posterior predictive results, threshold, data type                                 #
  # Output: Plot of Posterior predictive results                                              #
  # ========================================================================================= #

  thr_l <- length(thr)

  # Grid points for predictions
  delta   <- 0.005
  pts     <- seq(0, 1, delta)

  df_pts     <- data.frame(pts = pts)
  df_pts_nm  <- c("pts")

  norm_c  <- function(fx, deltax){
    sum(deltax*fx)
  }

  lower <- 0
  upper <- 1

  for (ss in 1:thr_l){
    mean_pred <- test_pred[[ss]][[1]]$pred_pt
    mp_norm   <- mean_pred/norm_c(mean_pred, delta)
    df_pts    <- cbind(df_pts, mp_norm)
    df_pts_nm <- c(df_pts_nm, paste0(thr[ss], sep=""))
  }

  p1_beta <- 0.1
  p2_beta <- 0.1
  db_sim_df    <- dbeta(pts, p1_beta, p2_beta)
  db_sim_norm  <- db_sim_df

  # Fix infinite case
  db_sim_norm[db_sim_norm == Inf] = 20
  df_pts[df_pts > 30] = 20

  df_sim_norm <- cbind(df_pts, db_sim_norm)
  df_sim_norm <- data.frame(df_sim_norm)
  colnames(df_sim_norm) <- c(df_pts_nm, "Real Density")

  df_sim_melt <- cbind(pts=rep(df_sim_norm[,1]), melt(df_sim_norm[,-1], variable.name = "thresh",
                                                      value.name = "value"))

  if (labs == T){

    pp <- ggplot(df_sim_melt, aes(x = pts, y = value, color = thresh)) + geom_line() +
      # ggtitle("Posterior Predictive Density") +
      xlab("Values") + ylab("Posterior Predictive Density") + labs(colour = "Threshold") +
      scale_y_continuous(limits = c(0, 20)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            text = element_text(size = 14))

  } else {
    pp <- ggplot(df_sim_melt, aes(x = pts, y = value, color = thresh)) + geom_line() +
      # ggtitle("Posterior Predictive Density") +
      xlab("Values") + ylab("Posterior Predictive Density") +
      scale_y_continuous(limits = c(0, 20)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none",
            # legend.text = element_text(colour = "white"),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            text = element_text(size = 14))
  }

  return(pp)
}

# ===================== Function to create Performance Plot ===================== #
#' Algorithm performance in Beta distribution
#'
#' This function generates a plot for the performance on a Beta distribution
#' @param time_all_result list of times for all result
#' @param ppred_all_result list of prediction for all result
#' @param thr vector of maximum number of rejections
#' @param labs whether a legend should be included in the figures
#' @return plot of the performance on a Beta distribution
#' @export
#'
func_perf_plot_beta <- function(time_all_result, ppred_all_result, thr, labs){
  # ========================================================================================= #
  # Function to create performance plot of all runs for Beta experiments                      #
  # Input: list of time results, list of posterior predictive results, threshold              #
  # Output: Plot of performance                                                               #
  # ========================================================================================= #

  thr_l  <- length(thr)
  num_cv <- length(time_all_result)

  time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
  time_all_mat_motg <- time_all_mat_motg/50
  time_all_avg <- apply(time_all_mat_motg, 2, median)
  time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  perf_all_mat <- unlist(lapply(1:num_cv, function(x){
    temp <- ppred_all_result[[x]]
  }))
  perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
  perf_all_avg <- apply(perf_all_mat_motg, 2, median)
  perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
               "perf_avg","5%", "25%", "75%", "95%")

  df_time  <- data.frame(thr = factor(thr), time_all_avg)
  df_time  <- cbind(df_time, t(time_all_qtl))
  df_time_perf  <- cbind(df_time, perf_all_avg)
  df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))

  colnames(df_time_perf) <- df_name

  limit_x <- c(floor(min(time_all_qtl[1,])), ceiling(max(time_all_qtl[4,])))

  a <- df_time_perf

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24, 25)) +
    guides(shape=FALSE, fill=FALSE) +
    # geom_errorbar(aes(ymin = `5%`, ymax = `95%`), stat = "identity",
    #               position = "identity", width = 0.01, size = 0.3) +
    # geom_errorbarh(aes(xmin = `5%_t`, xmax = `95%_t`), stat = "identity",
    #                position = "identity", size = 0.3, height = 3) +
    geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                  position = "identity", width = 0.1, size = 0.3) +
    geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                   position = "identity", size = 0.3, height = 4) +
    scale_x_log10(breaks = c(2, 5, 50, 200, 1000)) +
    xlab("Time in seconds per 100 samples") +
    ylab("Test Log Likelihood") +
    # ggtitle("Performance of Algorithm") +
    labs(color = "Threshold")

  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             text = element_text(size = 12),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             text = element_text(size = 12),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             legend.position = "none")

  return(p_time)
}

#' Algorithm performance in Gaussian distribution
#'
#' This function generates a plot for the performance on a Gaussian distribution
#' @param time_all_result list of times for all result
#' @param ppred_all_result list of prediction for all result
#' @param thr vector of maximum number of rejections
#' @param labs whether a legend should be included in the figures
#' @return plot of the performance on a Beta distribution
#' @export
#'
func_perf_plot_gauss <- function(time_all_result, ppred_all_result, thr, labs){
  # ========================================================================================= #
  # Function to create performance plot of all runs for Gaussian experiments                  #
  # Input: list of time results, list of posterior predictive results, threshold              #
  # Output: Plot of performance                                                               #
  # ========================================================================================= #

  thr_l  <- length(thr)
  num_cv <- length(time_all_result)

  time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
  time_all_mat_motg <- time_all_mat_motg/50
  time_all_avg <- apply(time_all_mat_motg, 2, median)
  time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  perf_all_mat <- unlist(lapply(1:num_cv, function(x){
    temp <- ppred_all_result[[x]]
  }))
  perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
  perf_all_avg <- apply(perf_all_mat_motg, 2, median)
  perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
               "perf_avg","5%", "25%", "75%", "95%")

  df_time  <- data.frame(thr = factor(thr), time_all_avg)
  df_time  <- cbind(df_time, t(time_all_qtl))
  df_time_perf  <- cbind(df_time, perf_all_avg)
  df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))

  colnames(df_time_perf) <- df_name

  limit_x <- c(floor(min(time_all_qtl[1,])), ceiling(max(time_all_qtl[4,])))

  a <- df_time_perf

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24, 25)) +
    guides(shape=FALSE, fill=FALSE)

  if (mix == "motg"){
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                     position = "identity", width = 7, size = 0.3) +
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", height = 0.5, size = 0.3) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      labs(color = "")

  } else {
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                  position = "identity", width = 0.05, size = 0.3) +
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", height = 0.5, size = 0.3) +
      # geom_errorbar(aes(ymin = `5%`, ymax = `95%`), stat = "identity",
      #               position = "identity", width = 0.5, size = 0.3) +
      # geom_errorbarh(aes(xmin = `5%_t`, xmax = `95%_t`), stat = "identity",
      #                position = "identity", size = 0.3, height = 0.5) +
      # scale_x_log10() +
      # scale_x_log10(breaks = c(3, 8, 100)) +
      # scale_x_log10(breaks = c(2, 2.5, 3)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "")
  }


  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             text = element_text(size = 12),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12))

  return(p_time)
}

#' Algorithm performance in bivariate Gaussian distribution
#'
#' This function generates a plot for the performance on a bivariate Gaussian distribution
#' @param time_all_result list of times for all result
#' @param ppred_all_result list of prediction for all result
#' @param thr vector of maximum number of rejections
#' @param labs whether a legend should be included in the figures
#' @return plot of the performance on a Beta distribution
#' @export
#'
func_perf_plot_bvgauss <- function(time_all_result, ppred_all_result, thr, labs){
  # ========================================================================================= #
  # Function to create performance plot of all runs for Bivariate Gaussian experiments        #
  # Input: list of time results, list of posterior predictive results, threshold              #
  # Output: Plot of performance                                                               #
  # ========================================================================================= #

  thr_l  <- length(thr)
  num_cv <- length(time_all_result)

  time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
  time_all_mat_motg <- time_all_mat_motg/50
  time_all_avg <- apply(time_all_mat_motg, 2, median)
  time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  perf_all_mat <- unlist(lapply(1:num_cv, function(x){
    temp <- ppred_all_result[[x]]
  }))
  perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
  perf_all_avg <- apply(perf_all_mat_motg, 2, median)
  perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
               "perf_avg","5%", "25%", "75%", "95%")

  df_time  <- data.frame(thr = factor(thr), time_all_avg)
  df_time  <- cbind(df_time, t(time_all_qtl))
  df_time_perf  <- cbind(df_time, perf_all_avg)
  df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))

  colnames(df_time_perf) <- df_name

  limit_x <- c(floor(min(time_all_qtl[1,])), ceiling(max(time_all_qtl[4,])))

  a <- df_time_perf

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24, 25)) +
    guides(shape=FALSE, fill=FALSE)

    if (mix == "tmogt"){
      # TMoG bigSD and smallSD
      # p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
      #               position = "identity", width = 0.5, size = 0.3) +
      #   geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
      #                  position = "identity", size = 0.3, height = 1) +

      # Crime
        p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                         position = "identity", width = 10, size = 0.3) +
          geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                         position = "identity", size = 0.3, height = 1) +

        # scale_x_log10() +
        # scale_x_log10(breaks = c(3, 8, 100)) +
        # scale_x_log10(breaks = c(2, 2.5, 3)) +
        xlab("Time in seconds per 100 samples") +
        ylab("Test Log Likelihood") +
        # ggtitle("Performance of Algorithm") +
        labs(color = "Threshold")



    } else {
      # MoTG bigSD
      p_time <- p_time +

      #   geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
      #               position = "identity", width = 12, size = 0.3) +
      # geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
      #                  position = "identity", size = 0.3, height = 10) +

        # MoTG smallSD
        geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                      position = "identity", width = 20, size = 0.3) +
        geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                       position = "identity", size = 0.3, height = 10) +

        # scale_x_log10() +
        # scale_x_log10(breaks = c(3, 8, 100)) +
        # scale_x_log10(breaks = c(2, 2.5, 3)) +
        xlab("Time in seconds per 100 samples") +
        ylab("Test Log Likelihood") +
        # ggtitle("Performance of Algorithm") +
        labs(color = "Threshold")
    }


  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12),
                             legend.position = "none")

  return(p_time)
}

#' Algorithm performance in application to flowcytometry
#'
#' This function generates a plot for the performance in application to flowcytometry data
#' @param time_all_result list of times for all result
#' @param ppred_all_result list of prediction for all result
#' @param thr vector of maximum number of rejections
#' @param labs whether a legend should be included in the figures
#' @return plot of the performance on a Beta distribution
#' @export
#'
func_perf_plot_flowcyt <- function(time_all_result, ppred_all_result, thr, labs){
  # ========================================================================================= #
  # Function to create performance plot of all runs for Flowcytometry experiments             #
  # Input: list of time results, list of posterior predictive results, threshold              #
  # Output: Plot of performance                                                               #
  # ========================================================================================= #

  thr_l  <- length(thr)
  num_cv <- length(time_all_result)

  time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
  time_all_mat_motg <- time_all_mat_motg/50
  time_all_avg <- apply(time_all_mat_motg, 2, median)
  time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  perf_all_mat <- unlist(lapply(1:num_cv, function(x){
    temp <- ppred_all_result[[x]]
  }))
  perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
  perf_all_avg <- apply(perf_all_mat_motg, 2, median)
  perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
    quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

  df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
               "perf_avg","5%", "25%", "75%", "95%")

  df_time  <- data.frame(thr = factor(thr), time_all_avg)
  df_time  <- cbind(df_time, t(time_all_qtl))
  df_time_perf  <- cbind(df_time, perf_all_avg)
  df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))

  colnames(df_time_perf) <- df_name

  limit_x <- c(floor(min(time_all_qtl[1,])), ceiling(max(time_all_qtl[4,])))

  a <- df_time_perf[-5,]

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24)) +
    guides(shape=FALSE, fill=FALSE)
    # geom_errorbar(aes(ymin = `5%`, ymax = `95%`), stat = "identity",
    #               position = "identity", width = 0.01, size = 0.3) +

  if (mix == "tmogt"){
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                         position = "identity", width = 3, size = 0.3) + #w =100
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 3) + #w = 20
      # scale_x_log10() +
      # scale_x_log10(breaks = c(50, 100, 200, 1000)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")
  } else {
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                     position = "identity", width = 400, size = 0.3) + #w =100
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 20) + #w = 20
      # scale_x_log10() +
      # scale_x_log10(breaks = c(50, 100, 200, 1000)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")
  }


  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12),
                             legend.position = "none")

  return(p_time)
}


# Substitute time in sub test to value in all test
perf_plot_edit_time_sub_flowcyto <- function(time_all_result_1, ppred_all_result_1, time_all_result_3, ppred_all_result_3, thr_1, thr_3, labs = T){

  thr_l  <- length(thr_1)
  thr_ll <- length(thr_3)
  num_cv <- length(time_all_result_1)

  summarize_time_pred <- function(time_all_result, ppred_all_result, thr_l, thr){

    time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
    time_all_mat_motg <- time_all_mat_motg/50
    time_all_avg <- apply(time_all_mat_motg, 2, median)
    time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
      quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

    perf_all_mat <- unlist(lapply(1:num_cv, function(x){
      temp <- ppred_all_result[[x]]
    }))
    perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
    perf_all_avg <- apply(perf_all_mat_motg, 2, median)
    perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
      quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

    df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
                 "perf_avg","5%", "25%", "75%", "95%")

    df_time  <- data.frame(thr = factor(thr), time_all_avg)
    df_time  <- cbind(df_time, t(time_all_qtl))
    df_time_perf  <- cbind(df_time, perf_all_avg)
    df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))
    colnames(df_time_perf) <- df_name

    return(df_time_perf)
  }

  df_time_perf_1 <- summarize_time_pred(time_all_result_1, ppred_all_result_1, thr_l, thr_1)
  df_time_perf_3 <- summarize_time_pred(time_all_result_3, ppred_all_result_3, thr_ll, thr_3)

  df_time_perf_3[, 2:6] <- df_time_perf_1[-5, 2:6]

  a <- df_time_perf_3[-5,]

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24)) +
    guides(shape=FALSE, fill=FALSE)
  # geom_errorbar(aes(ymin = `5%`, ymax = `95%`), stat = "identity",
  #               position = "identity", width = 0.01, size = 0.3) +

  if (mix == "tmogt"){
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                     position = "identity", width = 5, size = 0.3) + #w =100
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 2) + #w = 20
      # scale_x_log10() +
      # scale_x_log10(breaks = c(50, 100, 200, 1000)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")
  } else {
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                     position = "identity", width = 100, size = 0.3) + #w =100
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 20) + #w = 20
      # scale_x_log10() +
      # scale_x_log10(breaks = c(50, 100, 200, 1000)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")
  }


  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12),
                             legend.position = "none")

  return(p_time)

}

# Substitute time in sub test to value in all test
perf_plot_edit_time_sub_crime <- function(time_all_result_1, ppred_all_result_1, time_all_result_3, ppred_all_result_3, thr, labs = T){

  thr_l  <- length(thr)
  num_cv <- length(time_all_result_1)

  summarize_time_pred <- function(time_all_result, ppred_all_result, thr_l, thr){

    time_all_mat_motg <- t(matrix(unlist(time_all_result), nrow = thr_l, ncol = num_cv))
    time_all_mat_motg <- time_all_mat_motg/50
    time_all_avg <- apply(time_all_mat_motg, 2, median)
    time_all_qtl <- apply(time_all_mat_motg, 2, function(x)
      quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

    perf_all_mat <- unlist(lapply(1:num_cv, function(x){
      temp <- ppred_all_result[[x]]
    }))
    perf_all_mat_motg <- t(matrix(perf_all_mat, nrow = thr_l, ncol = num_cv))
    perf_all_avg <- apply(perf_all_mat_motg, 2, median)
    perf_all_qtl <- apply(perf_all_mat_motg, 2, function(x)
      quantile(x, probs = c(0.05, 0.25, 0.75, 0.95)))

    df_name <- c("thr", "time_avg", "5%_t", "25%_t", "75%_t", "95%_t",
                 "perf_avg","5%", "25%", "75%", "95%")

    df_time  <- data.frame(thr = factor(thr), time_all_avg)
    df_time  <- cbind(df_time, t(time_all_qtl))
    df_time_perf  <- cbind(df_time, perf_all_avg)
    df_time_perf  <- cbind(df_time_perf, t(perf_all_qtl))
    colnames(df_time_perf) <- df_name

    return(df_time_perf)
  }

  df_time_perf_1 <- summarize_time_pred(time_all_result_1, ppred_all_result_1, thr_l, thr)
  df_time_perf_3 <- summarize_time_pred(time_all_result_3, ppred_all_result_3, thr_l, thr)

  df_time_perf_3[, 2:6] <- df_time_perf_1[, 2:6]

  a <- df_time_perf_3

  p_time <- ggplot(a, aes(x = time_avg, y = perf_avg, color = thr)) +
    geom_point(data = a, aes(x = time_avg, y = perf_avg, shape = thr, fill = thr), size = 2, shape = c(21, 22, 23, 24, 25)) +
    guides(shape=FALSE, fill=FALSE)

  if (mix == "tmogt"){
    # TMoG bigSD and smallSD
    # p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
    #               position = "identity", width = 0.5, size = 0.3) +
    #   geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
    #                  position = "identity", size = 0.3, height = 1) +

    # Crime
    p_time <- p_time + geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                                     position = "identity", width = 10, size = 0.3) +
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 1) +

      # scale_x_log10() +
      # scale_x_log10(breaks = c(3, 8, 100)) +
      # scale_x_log10(breaks = c(2, 2.5, 3)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")



  } else {
    # MoTG bigSD
    p_time <- p_time +

      #   geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
      #               position = "identity", width = 12, size = 0.3) +
      # geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
      #                  position = "identity", size = 0.3, height = 10) +

      # MoTG smallSD
      geom_errorbar(aes(ymin = `25%`, ymax = `75%`), stat = "identity",
                    position = "identity", width = 0.2, size = 0.3) +
      geom_errorbarh(aes(xmin = `25%_t`, xmax = `75%_t`), stat = "identity",
                     position = "identity", size = 0.3, height = 1) +

      # scale_x_log10() +
      # scale_x_log10(breaks = c(3, 8, 100)) +
      # scale_x_log10(breaks = c(2, 2.5, 3)) +
      xlab("Time in seconds per 100 samples") +
      ylab("Test Log Likelihood") +
      # ggtitle("Performance of Algorithm") +
      labs(color = "Threshold")
  }


  if (labs == T)
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12))
  else
    p_time <- p_time + theme(plot.title = element_text(hjust = 0.5, size = 12),
                             axis.line = element_line(colour = "black"),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 12),
                             text = element_text(size = 12),
                             legend.position = "none")

}

# ===================== Function to create Error Plot ===================== #
#' Algorithm for Error in the posterior
#'
#' This function generates an error plot of the posterior distribution
#' @param test_pred predictive probability for the test set
#' @return plot of error for the predictive density on the test set
#' @export
#'
fig_error_post <- function(test_pred){
  # ========================================================================================= #
  # Function to create error plot for posterior predictive results                            #
  # Input: list of posterior predictive results                                               #
  # Output: Plot of error bars                                                                #
  # ========================================================================================= #

  tp_l  <- length(test_pred)

  thr    <- c(0, 0.1, 0.5, 1, 5, 50, Inf)
  thr_l  <- length(thr)

  tpred_sq_err <- lapply(1:tp_l, function(x)
    test_pred[[x]]$sq_err_summary)
  tpred_sq_err <- matrix(unlist(lapply(tpred_sq_err,
                                       function(x) apply(x, 1, mean))), nrow = 5)

  tpred_ab_err <- lapply(1:tp_l, function(x)
    test_pred[[x]]$abs_err_summary)
  tpred_ab_err <- matrix(unlist(lapply(tpred_ab_err,
                                       function(x) apply(x, 1, mean))), nrow = 5)

  ind    <- 1:thr_l

  par(mfrow=c(1,2))
  p <- plot(ind, tpred_sq_err[3,],
            ylim = range(c(-0.01 , max(tpred_sq_err[5,]))),
            pch=19, xaxt = 'n', frame.plot = T, xlab="Threshold",
            ylab="Median, 10%, 25%, 75%, and 90%",
            main="Median Absolute Error") +
    axis(1, ind, thr) +
    arrows(ind, tpred_sq_err[1, ], ind, tpred_sq_err[5,], length=0.05, angle=90, code=3) +
    arrows(ind, tpred_sq_err[2, ], ind, tpred_sq_err[4,], length=0.05, angle=90, code=3)

  p <- p + plot(ind, tpred_ab_err[3,],
                ylim = range(c(-0.01 , max(tpred_ab_err[5,]))),
                pch=19, xaxt = 'n', frame.plot = T, xlab="Threshold",
                ylab="Median, 5%, 25%, 75%,and 95%",
                main="Median Squared Error") +
    axis(1, ind, thr) +
    arrows(ind, tpred_ab_err[1,], ind, tpred_ab_err[5,], length=0.05, angle=90, code=3) +
    arrows(ind, tpred_ab_err[2, ], ind, tpred_ab_err[4,], length=0.05, angle=90, code=3)

  return(p)
}


# ===================== Plot Number of Rejections ===================== #
#' Number of rejection plot
#'
#' This function generates a plot of the number of rejections
#' @param count_bins number of rejections in each bins
#' @param mix one of "tmogt" and "motg"
#' @param thr_l lenght of the threshold vector
#' @return histogram for the number of rejections in each threshold
#' @export
#'
fig_numrejects <- function(count_bins, mix, thr_l){
  # ========================================================================================= #
  # Function to create histogram of the number of rejections                                  #
  # Input: matrix of count in bins, length of threshold and mixture type                      #
  # Output: Plot of performance plot                                                          #
  # ========================================================================================= #

  numIter <- length(count_bins[[1]])
  df_rej  <- data.frame(matrix(unlist(count_bins), nrow = numIter, ncol = thr_l))

  df_name <- c()
  for (tt in 1:thr_l){
    df_name  <- c(df_name, paste("Threshold_", thr[tt], sep=""))
  }

  colnames(df_rej) <- df_name

  barfill <- "#4271AE"
  barline <- "#1F3552"

  for (uu in 1:(thr_l)){
    pq <- ggplot(df_rej, aes_string(x = df_name[uu])) +
      geom_histogram(binwidth = 10, colour = barline, fill = barfill) +
      ggtitle(paste("Threshold", thr[uu])) +
      theme(plot.title = element_text(hjust = 0.5, size = 9),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            text = element_text(size = 8)) +
      xlab("Number of rejections")
    fig <- paste("p_", uu, sep = "")
    assign(fig, pq)
  }

  plot1 <- plot_grid(p_2, p_3, p_4, p_5)

  return(plot1)
}





