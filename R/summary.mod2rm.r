#' Print and summary function for objects of class "mod2rm"
#'
#' Prints a summary of a list object of class "mod2rm"
#' @param object An object of class "mod2rm"
#' @param ... Additional parameters
#' @return Prints summary of the object, then returns NULL
#' @examples 
#' df = data.frame(out1 = c(2,4,5,6,6), out2 = c(7,4,5,1,1), w1 = c(9,4,4,2,3), w2 = c(7,4,1,1,2))
#' res = mod2rm(df, out1, out2, w1)
#' res = mod2rm(df, out1, out2, w1, MOD1val = c(2,5))
#' res = mod2rm(df, out1, out2, w1, w2, method = 2, standardize = TRUE)
#' summary.mod2rm(res)
#' summary.mod2rm(res)
#' @export summary.mod2rm
#' @export
#' @importFrom stats confint lm na.omit pt qt sd t.test vcov
#' @importFrom utils capture.output
#' @importFrom methods is

summary.mod2rm <- function(object, ...) 
{
  if(!is(object,"mod2rm")) stop('Error: Please provide object of class "mod2rm"')
  
  # get important variables
  name_y1 = object$var_names["y1"]
  name_y2 = object$var_names["y2"]
  name_df = object$var_names["df"]
  num_mods = object$info["num_mods"]
  sample_size = object$info["sample_size"]
  method = object$info["method"]
  name_method = ifelse(method == 1, "additive","multiplicative")
  if(num_mods == 1) name_method = "single"
  
  # create output1
  
  cat("\n\n\n\n") 
  cat("**********Moderation Analysis for Two-Instance Repeated Measures Designs:**************\n\n")
  cat("v0.10 - For up to three simultaneous moderators (dichotomous or continuous).\n\n")
  cat("        Supports both additive (method = 1) and multiplicative (method = 2) moderation.\n\n")
  cat("2022  - Matthias Forstmann (matthias.forstmann@uzh.ch)\n\n")
  
  cat("**************************************************************************************\n\n")
  cat(paste("Model: ",num_mods," moderator(s)"," (",name_method," moderation)\n\n",sep = ""))
  cat(paste("Outcome Variables: Y = ",name_y1,", ", name_y2,"\n\n"))
  if (num_mods == 1) cat(paste("Moderator Variables: W1 =", object$var_names["mod1"],"\n\n"))
  if (num_mods == 2) cat(paste("Moderator Variables: W1 =", object$var_names["mod1"], "W2 =", object$var_names["mod2"],"\n\n"))
  if (num_mods == 3) cat(paste("Moderator Variables: W1 =", object$var_names["mod1"], "W2 =", object$var_names["mod2"], "W3 =", object$var_names["mod3"],"\n\n"))
  cat(paste("Computed Variables: Ydiff = ",name_y1," - ", name_y2,"\n\n"))
  cat(paste("Sample Size: ",sample_size,"\n\n"))
  
  # print object for comparison of Y1 and Y2
  
  cat("**************************************************************************************\n\n")
  cat("    Testing for difference between Y1 and Y2 \n\n")
  
  object$res_y1y2_diff$data.name =  paste(name_y1,"and", name_y2)
  out_t <- capture.output(print(object$res_y1y2_diff))
  out_t <- out_t[-c(1:2)]
  for(i in out_t) cat(paste(i,"\n"))

  # print object for test for moderation

  cat("**************************************************************************************\n\n")
  cat("    Test for Moderation of Repeated-Measures Factor\n\n")
  cat(paste("Outcome: Ydiff = ",  name_y1, " - ", name_y2,"\n\n"))

  if(num_mods == 1) 
    names(object$res_mod$coefficients)[2] = object$var_names["mod1"]
  if(num_mods == 2) {
    names(object$res_mod$coefficients)[2] = object$var_names["mod1"]
    names(object$res_mod$coefficients)[3] = object$var_names["mod2"]
    if(method == 2)
      names(object$res_mod$coefficients)[4] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
  }
  if(num_mods == 3) {
    names(object$res_mod$coefficients)[2] = object$var_names["mod1"]
    names(object$res_mod$coefficients)[3] = object$var_names["mod2"]
    names(object$res_mod$coefficients)[4] = object$var_names["mod3"]
    if(method == 2)
      names(object$res_mod$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
      names(object$res_mod$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
      names(object$res_mod$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
      names(object$res_mod$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
  }
  
  out1_ci = cbind(summary(object$res_mod)$coefficients, confint(object$res_mod))
  print(out1_ci)
  out1 <- capture.output(print(summary(object$res_mod)))
  out1 <- out1[-c(1:(length(out1)-5))]
  for(i in out1) cat(paste(i,"\n"))
  cat("-----\n")
  
  

  cat("\n\n**************************************************************************************\n\n")
  cat("    Conditional Effect of Moderator(s) on Y in each Condition\n\n")

  # print results of simple regression for condition 1

  cat(paste("Condition 1 Outcome: ",name_y1,"\n\n"))
  if(num_mods == 1) 
    names(object$res_simple_y1$coefficients)[2] = object$var_names["mod1"]
  if(num_mods == 2) {
    names(object$res_simple_y1$coefficients)[2] = object$var_names["mod1"]
    names(object$res_simple_y1$coefficients)[3] = object$var_names["mod2"]
    if(method == 2)
      names(object$res_simple_y1$coefficients)[4] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
  }
  if(num_mods == 3) {
    names(object$res_simple_y1$coefficients)[2] = object$var_names["mod1"]
    names(object$res_simple_y1$coefficients)[3] = object$var_names["mod2"]
    names(object$res_simple_y1$coefficients)[4] = object$var_names["mod3"]
    if(method == 2)
      names(object$res_simple_y1$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
    names(object$res_simple_y1$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
    names(object$res_simple_y1$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
    names(object$res_simple_y1$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
  }
  
  out2a_ci = cbind(summary(object$res_simple_y1)$coefficients, confint(object$res_simple_y1))
  print(out2a_ci)
  out2a <- capture.output(print(summary(object$res_simple_y1)))
  out2a <- out2a[-c(1:(length(out2a)-5))]
  for(i in out2a) cat(paste(i,"\n"))
  cat("-----\n")

  # print results of simple regression for condition 2

  cat(paste("Condition 2 Outcome: ",name_y2,"\n\n"))
  if(num_mods == 1) 
    names(object$res_simple_y2$coefficients)[2] = object$var_names["mod1"]
  if(num_mods == 2) {
    names(object$res_simple_y2$coefficients)[2] = object$var_names["mod1"]
    names(object$res_simple_y2$coefficients)[3] = object$var_names["mod2"]
    if(method == 2)
      names(object$res_simple_y2$coefficients)[4] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
  }
  if(num_mods == 3) {
    names(object$res_simple_y2$coefficients)[2] = object$var_names["mod1"]
    names(object$res_simple_y2$coefficients)[3] = object$var_names["mod2"]
    names(object$res_simple_y2$coefficients)[4] = object$var_names["mod3"]
    if(method == 2)
      names(object$res_simple_y2$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
    names(object$res_simple_y2$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
    names(object$res_simple_y2$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
    names(object$res_simple_y2$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
  }
  
  out2b_ci = cbind(summary(object$res_simple_y2)$coefficients, confint(object$res_simple_y2))
  print(out2b_ci)
  out2b <- capture.output(print(summary(object$res_simple_y2)))
  out2b <- out2b[-c(1:(length(out2b)-5))]
  for(i in out2b) cat(paste(i,"\n"))
  cat("-----\n")

  # print results for conditional effects

  cat("\n\n**************************************************************************************\n\n")
  cat("   Conditional Effects of Condition on Y at Values of the Moderator(s)\n\n")
  
  colnames(object$res_cond_eff)[1] = paste(object$var_names["mod1"])
  if(num_mods > 1) 
    colnames(object$res_cond_eff)[2] = paste(object$var_names["mod2"])
  if(num_mods > 2) 
    colnames(object$res_cond_eff)[3] = paste(object$var_names["mod3"])
  
  #round predictor values if standardized (for display purposes), only if standardized
  if(object$info["standardize"] == TRUE) object$res_cond_eff[,1] = round(object$res_cond_eff[,1],7)
  if(num_mods > 1) 
    if(object$info["standardize"] == TRUE) object$res_cond_eff[,2] = round(object$res_cond_eff[,2],7)
  if(num_mods > 2) 
    if(object$info["standardize"] == TRUE) object$res_cond_eff[,3] = round(object$res_cond_eff[,3],7)
  

  print(object$res_cond_eff)
  
  cat("\n\nValues for quantitative moderators are the mean and plus/minus one SD from the mean (unless manually defined).\n\n")
  cat("Values for binary moderators are conditional effects at both values.\n\n")

  invisible(NULL)
  
}



