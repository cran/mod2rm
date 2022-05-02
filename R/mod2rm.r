#' Moderation Analysis for Two-Instance Repeated Measures Designs
#'
#' Moderation analysis for two-instance repeated measures designs, including analyses of simple slopes and conditional effects at values of the moderator.\cr\cr
#' Currently supports a single continuous or binary moderator. For continuous moderators, two values for conditional effects tests can be specified manually.\cr\cr
#' Method and design based on: Montoya, A. K. (2018). Moderation analysis in two-instance repeated measures designs: Probing methods and multiple moderator models. \emph{Behavior Research Methods, 51(1)}, 61-82.\cr\cr
#' @param data A data frame
#' @param Y1 Name of the first outcome variable
#' @param Y2 Name of the second outcome variable
#' @param MOD1 Name of the moderator variable
#' @param MOD1val A vector containing two values of the moderator at which to test for conditional effects (optional).
#' @param method Method for dealing with two or more moderators (1 = addition, 2 = multiplication) (currently not supported)
#' @return 
#'   \item{total}{A list of class "mod2rm" containing:}
#'   \item{num_mods}{Numeric variable indicating the number of moderators (currently only 1 moderator is supported)}
#'   \item{mod_binary}{Boolean value indicating whether the moderator is binary or not}
#'   \item{name_y1}{String variable providing the name of first dependent variable (y1)}
#'   \item{name_y2}{String variable providing the name of second dependent variable (y2)}
#'   \item{name_mod}{String variable providing the name of moderator}
#'   \item{res_mod}{A list including the results of a simple regression, regressing the difference between y1 and y2 on the moderator}
#'   \item{simp_eff1}{A list including the results of a simple regression, regressing the y1 on the moderator}
#'   \item{simp_eff2}{A list including the results of a simple regression, regressing the y2 on the moderator}
#'   \item{cond_eff}{A list including the results of an analysis of conditional effects at different levels of the moderator}
#'   \item{y1y2_diff}{A list including the results of a repeated measures t-test for y1 and y2}
#'   \item{sample_size}{Numeric variable indicating the number of complete cases that the analysis is based on}
#' 
#' @examples 
#' df = data.frame(out1 = c(2,4,5,6,6), out2 = c(7,4,5,1,1), w = c(9,4,4,2,3))
#' res = mod2rm(df, out1, out2, w)
#' res = mod2rm(df, out1, out2, w, c(2,5))
#' summary.mod2rm(res)
#' @export
#' @importFrom stats confint lm na.omit pt qt sd t.test vcov
#' @importFrom utils capture.output

mod2rm <- function (data, Y1, Y2, MOD1, MOD1val = NULL, method = 1) 
{
  
  # test arguments
  
  if( method != 1)  stop('only method 1 (addititive relationship between moderators) is currently supported')
  if(length(MOD1val) != 2 & !missing(MOD1val)) stop('vector of self-set values for moderators must have a length of 2')
  
  num_mods = 1 # total number of moderators. not needed right now
  
  
  # testing if moderator(s) is/are binary. Give error if specific values are assigned
  
  mod_binary = ifelse((length(unique(na.omit(eval(substitute(MOD1), data)))) < 3),TRUE, FALSE)  # detect if moderator 1 is dichotomous or not
  
  if(mod_binary == TRUE & !missing(MOD1val)) stop('it is not possible to define values for binary moderators')
  
  
  
  # create dataframe and critical variables for use within this function
  
  df = na.omit(data.frame(y1 = eval(substitute(Y1), data), y2 = eval(substitute(Y2), data), mod = eval(substitute(MOD1), data)))
  y1 = df$y1
  y2 = df$y2
  mod = df$mod
  
  diff = y1 - y2
  
  #run regression, predicting difference score from moderator (i.e., testing the interaction)
  
  erg = lm(diff ~ mod)
  erg_sum = summary(erg)
  erg_sum$coefficients <- cbind(erg_sum$coefficients, confint(erg))  # add CIs to output
  erg_sum$coefficients <- round(erg_sum$coefficients, digits = 4) # round everything to 4 digits
  dimnames(erg_sum$coefficients)[[2]][5:6] = c("LLCI","ULCI")  # rename confidence intervals
  
  # look at simple effects of moderator on condition
  
  erg2a = lm(y1 ~ mod)
  erg2a_sum = summary(erg2a)
  erg2a_sum$coefficients <- cbind(erg2a_sum$coefficients, confint(erg2a))  # add CIs to output
  erg2a_sum$coefficients <- round(erg2a_sum$coefficients, digits = 4) # round everything to 4 digits
  dimnames(erg2a_sum$coefficients)[[2]][5:6] = c("LLCI","ULCI")  # rename confidence intervals
  
  erg2b = lm(y2 ~ mod)
  erg2b_sum = summary(erg2b)
  erg2b_sum$coefficients <- cbind(erg2b_sum$coefficients, confint(erg2b))  # add CIs to output
  erg2b_sum$coefficients <- round(erg2b_sum$coefficients, digits = 4) # round everything to 4 digits
  dimnames(erg2b_sum$coefficients)[[2]][5:6] = c("LLCI","ULCI")  # rename confidence intervals
  
  
  # critical: effects at different values of the moderator (-+ 1 SD in this case)
  
  low_val_of_moderator = mean(mod) - sd(mod)      # choose -1 SD of moderator
  high_val_of_moderator = mean(mod) + sd(mod)     # choose +1 SD of moderator
  
  if(mod_binary == TRUE) {                            # in case of binary moderator, select both scores
    low_val_of_moderator = sort(unique(mod))[1]
    high_val_of_moderator =sort (unique(mod))[2]   
  }
  
  if(length(MOD1val) == 2 & !missing(MOD1val)) {  # in case user specified two testing points
    low_val_of_moderator = min(MOD1val)
    high_val_of_moderator = max(MOD1val)
    mod_binary = TRUE
  }
  
  # manually calculate conditional effects at low and high values of the moderator
  
  unst_coeff_interc = as.numeric(erg$coefficients[1])      # get unstandardized coefficient intercept
  unst_coeff_moderator = as.numeric(erg$coefficients[2])   # get unstandardized coefficient moderator
  
  gradient_low = unst_coeff_interc + (low_val_of_moderator * unst_coeff_moderator)   # calculate gradient at -1 SD
  gradient_high = unst_coeff_interc + (high_val_of_moderator * unst_coeff_moderator)   # calculate gradient at -1 SD
  
  var_intercp = vcov(erg)[1,1]                      # get variance for intercept
  var_coeff_moderator = vcov(erg)[2,2]              # get variance for moderator
  covar_coeffi_mod_and_intercept = vcov(erg)[1,2]   # get covariance for both
  
  gradient_low_se = sqrt(var_intercp + low_val_of_moderator^2 * var_coeff_moderator + 2*low_val_of_moderator*covar_coeffi_mod_and_intercept)
  gradient_high_se = sqrt(var_intercp + high_val_of_moderator^2 * var_coeff_moderator + 2*high_val_of_moderator*covar_coeffi_mod_and_intercept)
  
  tvalue_low = gradient_low / gradient_low_se        #calculate t-value of slopes
  tvalue_high = gradient_high / gradient_high_se
  
  
  pvalue_low = 2*pt(q=abs(tvalue_low), df=nrow(df) - 2, lower.tail=FALSE)     # look up p value from that t value and the degrees of freedom (n minus 2) 
  pvalue_high = 2*pt(q=abs(tvalue_high), df=nrow(df) - 2, lower.tail=FALSE)
  
  # calculate 95% confidence intervals of the conditional effects
  
  margin_low <- qt(p = 0.05 / 2, df = (nrow(df)-2), lower.tail=F) * gradient_low_se
  margin_high <- qt(p = 0.05 / 2, df = (nrow(df)-2), lower.tail=F) * gradient_high_se
  
  # look up p value from that t value and the degrees of freedom (n minus 2, for one moderator)
  
  conf_low_lower = gradient_low - margin_low
  conf_low_upper = gradient_low + margin_low
  conf_high_lower = gradient_high - margin_high
  conf_high_upper = gradient_high + margin_high 
  
  # do the same thing for mean values, if moderator is not binary!
  if (mod_binary == FALSE) {
    mean_val_of_moderator = mean(mod)   
    gradient_mean = unst_coeff_interc + (mean_val_of_moderator * unst_coeff_moderator)   
    gradient_mean_se = sqrt(var_intercp + mean_val_of_moderator^2 * var_coeff_moderator + 2*mean_val_of_moderator*covar_coeffi_mod_and_intercept)
    tvalue_mean = gradient_mean / gradient_mean_se
    pvalue_mean = 2*pt(q=abs(tvalue_mean), df=nrow(df) - 2, lower.tail=FALSE)
    margin_mean <- qt(p = 0.05 / 2, df = (nrow(df)-2), lower.tail=F) * gradient_mean_se
    conf_mean_lower = gradient_mean - margin_mean
    conf_mean_upper = gradient_mean + margin_mean 
  }
  
  # save manual conditional effects into list
 
  if (mod_binary == TRUE) {
    cond = list(mod_val_low = low_val_of_moderator, mod_grad_low = gradient_low, mod_grad_low_se = gradient_low_se,tval_low = tvalue_low, pval_low = pvalue_low,conf_low_lower = conf_low_lower, conf_low_upper = conf_low_upper, 
                mod_val_high = high_val_of_moderator, mod_grad_high = gradient_high, mod_grad_high_se = gradient_high_se,tval_high = tvalue_high, pval_high = pvalue_high,conf_high_lower = conf_high_lower, conf_high_upper = conf_high_upper)
  } 
  if (mod_binary == FALSE) {
    cond = list(mod_val_low = low_val_of_moderator, mod_grad_low = gradient_low, mod_grad_low_se = gradient_low_se,tval_low = tvalue_low, pval_low = pvalue_low,conf_low_lower = conf_low_lower, conf_low_upper = conf_low_upper, 
                mod_val_mean = mean_val_of_moderator, mod_grad_mean = gradient_mean, mod_grad_mean_se = gradient_mean_se,tval_mean = tvalue_mean, pval_mean = pvalue_mean,conf_mean_lower = conf_mean_lower, conf_mean_upper = conf_mean_upper, 
                mod_val_high = high_val_of_moderator, mod_grad_high = gradient_high, mod_grad_high_se = gradient_high_se,tval_high = tvalue_high, pval_high = pvalue_high,conf_high_lower = conf_high_lower, conf_high_upper = conf_high_upper)
    
  }
 
  # t-test difference Y1 and Y2
  
  diff_ttest = t.test(y1, y2, paired = T)
  
  # save everything in new object
  
  results = list(num_mods = num_mods, mod_binary = mod_binary, name_y1 = paste(deparse(substitute(Y1))), name_y2 = paste(deparse(substitute(Y2))), name_mod1 = paste(deparse(substitute(MOD1))), res_mod = erg_sum, simp_eff1 = erg2a_sum, simp_eff2 = erg2b_sum, cond_eff = cond, y1y2_diff = diff_ttest, sample_size = nrow(df))
  
  class(results) <- c("mod2rm")
  
  return(results)
  
}


