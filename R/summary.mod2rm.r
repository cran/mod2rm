#' Print and summary function for objects of class "mod2rm"
#'
#' Prints a summary of a list object of class "mod2rm", and (if requested) plot the results of the Johnson-Neyman procedure.
#' @param object An object of class "mod2rm"
#' @param ... Additional parameters. "plotjn = TRUE" produces a ggplot for the Johnson-Neyman procedure, "plotstyle" can be set to "simple" or "points" (including data points in the plot)
#' @return Prints summary of the object, then returns NULL or (when requested) a ggplot2 object for the Johnson-Neyman plot
#' @details This function produces a summary for the results of an object of the type mod2rm, and can further be used to plot a graph showing the results of the JN procedure, if it is included in the mod2rm object.
#' #' Results include number and name(s) of moderator(s), sample size, results of a paired t-test between both dependent variables, the results of the moderation analysis, conditional effects of the moderator on each of the dependent variables, conditional effects at values of the moderator, and the results of the Johnson-Neyman procedure (including critical values, proportion of the sample above/below these values, and conditional effects around the significance regions. 
#' @examples 
#' 
#' # Generate a dataset with a Johnson-Neyman (non-)significance region within the response range:
#' 
#' repeat{
#'   df = data.frame(out1 = runif(n = 100, min = 1, max = 9), 
#'                   out2 = runif(n = 100, min = 1, max = 9), 
#'                   w1 = runif(n = 100, min = 1, max = 9),  
#'                   w2 = runif(n = 100, min = 1, max = 9),
#'                   w3 = runif(n = 100, min = 1, max = 9))
#'   res = mod2rm(df, out1, out2, w1, jn = TRUE)
#'   if(res$res_jn_area["num_jn"] == 2 & res$res_jn_area["center_significant"] == FALSE)
#'     break
#' }
#' 
#' # Show summary including plot
#' summary.mod2rm(res, plotjn = TRUE, plotstyle = "simple")
#' 
#' # Multiple regression (3 moderators, additive)
#' res1 = mod2rm(df, out1, out2, w1, w2, w3, method = 1)
#' summary.mod2rm(res1)
#' 
#' # Multiple regression (2 moderators, multiplicative, manually defined conditional effects)
#' res2 = mod2rm(df, out1, out2, w1, w2, MOD1val = c(2,3,4), MOD2val = c(4,5), method = 2)
#' summary.mod2rm(res2)
#' 
#' 
#' @export summary.mod2rm
#' @export
#' @importFrom stats confint lm na.omit pt qt sd t.test vcov
#' @importFrom utils capture.output
#' @importFrom methods is
#' @importFrom ggplot2 .data ggplot aes element_text labs xlab ylab geom_hline geom_vline stat_smooth geom_point scale_x_continuous scale_y_continuous scale_color_manual theme_minimal theme 
#' @importFrom scales pretty_breaks
#' 
summary.mod2rm <- function(object, ...) 
{
  if(!is(object,"mod2rm")) stop('Error: Please provide object of class "mod2rm"')
  
  args <- list(...)
  
  # test if and how JN plot should be plotted
  
  if(!is.null(args[["plotjn"]]))
    plotjn = ifelse(args$plotjn == TRUE, TRUE, FALSE)
  if(is.null(args[["plotjn"]])) plotjn = FALSE
  
  if(!is.null(args[["plotstyle"]])) {
    plotstyle = args$plotstyle
    if(plotstyle != "simple" & plotstyle != "points")
      stop("plotstyle musst be 'simple' or 'points'")
  }
  if(is.null(args[["plotstyle"]])) plotstyle = "simple"
  
  
  # get important variables
  name_y1 = object$var_names["y1"]
  name_y2 = object$var_names["y2"]
  name_df = object$var_names["df"]
  mod1_values = object$res_mod$model$mod1    # individual moderator values, needed for plot
  num_mods = object$info["num_mods"]
  sample_size = object$info["sample_size"]
  method = object$info["method"]
  name_method = ifelse(method == 1, "additive","multiplicative")
  if(num_mods == 1) name_method = "single"
  
  # create output1
  
  cat("\n\n\n\n") 
  cat("********* Moderation Analysis for Two-Instance Repeated Measures Designs: ************\n\n")
  cat("v0.2.1 - For up to three simultaneous moderators (dichotomous or continuous).\n\n")
  cat("         Supports both additive (method = 1) and multiplicative (method = 2) moderation,\n\n")
  cat("         conditional effects, simple slopes, and Johnson-Neyman procedure (jn = 1).\n\n")
  cat("2022   - Matthias Forstmann (matthias.forstmann@uzh.ch)\n\n")
  
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
    if(method == 2) {
      names(object$res_mod$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
      names(object$res_mod$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
      names(object$res_mod$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
      names(object$res_mod$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
    }
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
    if(method == 2) {
      names(object$res_simple_y1$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
      names(object$res_simple_y1$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
      names(object$res_simple_y1$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
      names(object$res_simple_y1$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
    }
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
    if(method == 2) {
      names(object$res_simple_y2$coefficients)[5] = paste(object$var_names["mod1"],":",object$var_names["mod2"])
      names(object$res_simple_y2$coefficients)[6] = paste(object$var_names["mod1"],":",object$var_names["mod3"])
      names(object$res_simple_y2$coefficients)[7] = paste(object$var_names["mod2"],":",object$var_names["mod3"])
      names(object$res_simple_y2$coefficients)[8] = paste(object$var_names["mod1"],":",object$var_names["mod2"],":",object$var_names["mod3"])
    }
    
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


  print(object$res_cond_eff, row.names = F)

  cat("\n\nValues for quantitative moderators are the mean and plus/minus 1 SD from the mean (unless manually defined).\n\n")
  if(object$info["num_binary_mods"] > 0)
    cat("Values for binary moderators are conditional effects at both values.\n\n")

  # Johnson-Neyman Procedure
  
  
  if(object$info["jn"] == TRUE) {
    cat("\n\n***************************** Johnson-Neyman-Procedure *******************************\n\n")
    cat(paste("Number of JN significance points (p = .05) identified: ",object$res_jn_area["num_jn"],"\n\n"))
    if (object$res_jn_area["num_jn"] == 0)
      cat("No JN significance region has been identified, or both values are outside of the moderator's value range.\n\n")
    if (object$res_jn_area["num_jn"] > 0) {
      cat("Moderator value(s) defining JN significance region(s):\n\n")
      
      out_jn = data.frame(value = object$res_jn_area$jn_values[,1],"percent-below" = 100 - object$res_jn_area$jn_values[,2], "percent-above" = object$res_jn_area$jn_values[,2])
      out_jn = subset(out_jn, object$res_jn_area$jn_values[,2] != 100 & object$res_jn_area$jn_values[,2] != 0) # remove out of bounds      
      colnames(out_jn)[1] = paste(object$var_names["mod1"])
      
      print(out_jn, row.names = F)
      
      cat("\nPercentages indicate the proportion of the sample with higher/lower scores.")
      
      cat("\n\n**************************************************************************************\n\n")
      cat("   Conditional Effects around the JN significance region(s)\n\n")
      
      colnames(object$res_jn_cond_eff)[1] = paste(object$var_names["mod1"])
      if(object$info["standardize"] == TRUE) object$res_jn_cond_eff[,1] = round(object$res_jn_cond_eff[,1],7)
      
      print(object$res_jn_cond_eff, row.names = F)
      
      cat("\n\nValues for quantitative moderators represent the sample response rage in 20 intervals.\n\n")
      
      }
    
  }

  # plot JN with ggplot2 
  
  if(plotjn == TRUE) {

      # define which colors to use, depending on whether the JN region is significant or non-significant
      
      col_outer = ifelse(object$res_jn_area$center_significant == FALSE,"#00BFC4","#F8766D")
      col_inner = ifelse(object$res_jn_area$center_significant == TRUE,"#00BFC4","#F8766D")
      
      p <- ggplot(object$res_mod, aes(x = .data$mod1, y = .data$diff), environment = environment()) + 
        labs(title = "Johnson-Neyman Plot\n") +
        xlab(paste("\n", object$var_names["mod1"])) +
        ylab(paste("diff: ",object$var_names["y1"], " - ", object$var_names["y2"], "\n")) +
        geom_hline(yintercept=0) +
        geom_vline(aes(xintercept=object$res_jn_area$jn_values[1,1], color = "non-significant\narea"), linetype="dashed", size=0.8) +
        geom_vline(aes(xintercept=object$res_jn_area$jn_values[2,1], color = "non-significant\narea"), linetype="dashed", size=0.8) +
        
   
        {if(object$res_jn_area$jn_values[1,1] > min(mod1_values))
          stat_smooth(formula = "y ~ x", method = "lm", color = col_outer, fill = col_outer, alpha = 0.2, xseq = seq(min(mod1_values), object$res_jn_area$jn_values[1,1], length = 2000))
        }+
        
        {if(object$res_jn_area$jn_values[2,1] < max(mod1_values))
          stat_smooth(formula = "y ~ x", method = "lm", color = col_outer, fill = col_outer, alpha = 0.2, xseq = seq(object$res_jn_area$jn_values[2,1], max(mod1_values), length = 2000)) 
        }+
        
        stat_smooth(formula = "y ~ x", method = "lm", color = col_inner, fill = col_inner, alpha = 0.2, xseq = seq(object$res_jn_area$jn_values[1,1],object$res_jn_area$jn_values[2,1], length = 2000)) +
        
        geom_vline(aes(xintercept=max(mod1_values), color = "data range"), linetype="dashed", size=0.8) +
        geom_vline(aes(xintercept=min(mod1_values), color = "data range"), linetype="dashed", size=0.8) +
        
        {if(plotstyle == "points")
          geom_point(color = "blue", alpha = .5)  
        }+
        
        scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 15)) +
        
        scale_color_manual(name = "", values = c("data range" = "green", "non-significant\narea" = "#F8766D")) +
        
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) 
    
  }
  
  #return ggplot to environment if it was drawn, otherwise just return NULL
  
  if(plotjn == TRUE)
      return(p)
  
  if(plotjn == FALSE)
    invisible(NULL)

}






