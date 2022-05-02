#' Print and summary function for objects of class "mod2rm"
#'
#' Prints a summary of a list object of class "mod2rm"
#' @param object An object of class "mod2rm"
#' @param ... Additional parameters
#' @return Prints summary of the object, then returns NULL
#' @examples 
#' df = data.frame(out1 = c(2,4,5,6,6), out2 = c(7,4,5,1,1), w = c(9,4,4,2,3))
#' res = mod2rm(df, out1, out2, w)
#' res = mod2rm(df, out1, out2, w, c(2,5))
#' summary.mod2rm(res)
#' @export summary.mod2rm
#' @export
#' @importFrom stats confint lm na.omit pt qt sd t.test vcov
#' @importFrom utils capture.output

summary.mod2rm <- function(object, ...) 
{
  if(class(object) != "mod2rm") stop('Error: Please provide object of class "mod2rm"')
  
  # create output
  
  cat("\n\n\n\n") 
  cat("**********Moderation Analysis for Two-Instance Repeated Measures Designs:**************\n\n")
  cat("Based on MEMORE by Amanda Montoya (akmontoya.com)\n\n")
  cat("v0.02 - For a single continuous or binary moderator\n\n")
  cat("2022 - Matthias Forstmann (matthias.forstmann@uzh.ch)\n\n")
  
  cat("**************************************************************************************\n\n")
  cat("    Model: one single moderator\n\n")
  cat(paste("Variables: Y = ",object$name_y1,", ", object$name_y2,"    W = ", object$name_mod1,"\n\n"))
  cat(paste("Computed Variables: Ydiff = ",object$name_y1," - ", object$name_y2,"\n\n"))
  cat(paste("Sample Size: ",object$sample_size,"\n\n"))
  
  # print object for comparison of Y1 and Y2
  
  cat("**************************************************************************************\n\n")
  cat("    Testing for difference between Y1 and Y2 \n\n")
  
  
  out_t <- capture.output(print(object$y1y2_diff))
  out_t <- gsub("y1", object$name_y1 , out_t)
  out_t <- gsub("y2", object$name_y2 , out_t)
  out_t <- out_t[-c(1:2)]
  for(i in out_t) cat(paste(i,"\n"))
  
  
  # print object for test for moderation
  
  cat("**************************************************************************************\n\n")
  cat("    Test for Moderation\n\n")
  cat(paste("Outcome: Ydiff = ",  object$name_y1, " - ", object$name_y2,"\n\n"))
  
  out1 <- capture.output(print(object$res_mod))
  out1 <- gsub("mod", object$name_mod1 , out1)
  out1 <- out1[-c(1:9)]
  for(i in out1) cat(paste(i,"\n"))
  
  cat("\n\n**************************************************************************************\n\n")
  cat("    Conditional Effect of Moderator(s) on Y in each Condition\n\n")
  
  # print results of simple regression for condition 1
  
  cat(paste("Condition 1 Outcome: ",object$name_y1,"\n\n"))
  out2a <- capture.output(print(object$simp_eff1))
  out2a <- gsub("mod", object$name_mod1 , out2a)
  out2a <- out2a[-c(1:9)]
  for(i in out2a) cat(paste(i,"\n"))
  cat("-----\n")
  
  # print results of simple regression for condition 2
  
  cat(paste("\nCondition 2 Outcome: ",object$name_y2,"\n\n"))
  out2b <- capture.output(print(object$simp_eff2))
  out2b <- gsub("mod", object$name_mod1, out2b)
  out2b <- out2b[-c(1:9)]
  for(i in out2b) cat(paste(i,"\n"))
  
  # print results for conditional effects
  
  cat("\n\n**************************************************************************************\n\n")
  cat("   Conditional Effect of repeated measures factor ('X') on Y at values of moderator(s)\n\n")
  
  cat(paste(sprintf("%10s%10s%10s%10s%10s%10s%10s",object$name_mod1,"Effect", "SE","t-value", "p-value", "LLCI", "ULCI")),"\n")
  cat(paste(sprintf("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",object$cond_eff$mod_val_low, object$cond_eff$mod_grad_low, object$cond_eff$mod_grad_low_se,object$cond_eff$tval_low, object$cond_eff$pval_low,object$cond_eff$conf_low_lower, object$cond_eff$conf_low_upper)),"\n")
  if(!object$mod_binary) cat(paste(sprintf("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",object$cond_eff$mod_val_mean, object$cond_eff$mod_grad_mean, object$cond_eff$mod_grad_mean_se,object$cond_eff$tval_mean, object$cond_eff$pval_mean,object$cond_eff$conf_mean_lower, object$cond_eff$conf_mean_upper)),"\n")
  cat(paste(sprintf("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",object$cond_eff$mod_val_high, object$cond_eff$mod_grad_high, object$cond_eff$mod_grad_high_se,object$cond_eff$tval_high, object$cond_eff$pval_high,object$cond_eff$conf_high_lower, object$cond_eff$conf_high_upper)),"\n\n")
  cat(paste("Degrees of freedom for all conditional effects:\n\t",as.numeric(object$sample_size)-2,"\n\n"))
  if(!object$mod_binary) cat("Values for quantitative moderators are the mean and plus/minus one SD from the mean.\n\n")
  if(object$mod_binary) cat("Values for binary moderators are conditional effects at both values.\n\n")
  
  #return invisible
  
  invisible(NULL)
  
}



