#' Moderation Analysis for Two-Instance Repeated Measures Designs
#'
#' Multiple moderation analysis for two-instance repeated measures designs, including analyses of simple slopes and conditional effects at values of the moderator(s).\cr\cr
#' Currently supports both single- and multi-moderator models, with up to three simultaneous moderators (continuous and/or binary). Multi-moderator models support both additive (method = 1) and multiplicative (method = 2) moderation.\cr\cr
#' Also supports the Johnson-Neyman procedure for determining regions of significance in single moderator models (jn = T). Plots of the JN region can be obtained from the summary function (plotjn = T).\cr\cr
#' Moderator values at which to test for conditional effects are determined automatically (at -1, 0, and +1 SD of the mean if the moderator is countinuous, and at both values of the moderator if it is binary), but any number of test values can also be set manually for each moderator.\cr\cr
#' Method and output based on: Montoya, A. K. (2018). Moderation analysis in two-instance repeated measures designs: Probing methods and multiple moderator models. \emph{Behavior Research Methods, 51(1)}, 61-82.\cr\cr
#' @param data A data frame
#' @param Y1 Name of the first outcome variable
#' @param Y2 Name of the second outcome variable
#' @param MOD1 Name of moderator1 variable
#' @param MOD2 Name of moderator2 variable (optional)
#' @param MOD3 Name of moderator3 variable (optional)
#' @param MOD1val A vector containing values of moderator1 at which to test for conditional effects (even when variables have been standardized!)(optional)
#' @param MOD2val A vector containing values of moderator2 at which to test for conditional effects (even when variables have been standardized!)(optional)
#' @param MOD3val A vector containing values of moderator3 at which to test for conditional effects (even when variables have been standardized!)(optional)
#' @param method Method for dealing with two or more moderators (1 = additive, 2 = multiplicative) (default: additive)
#' @param standardize boolean variable indicating whether all predictor variables (moderators) should be standardized prior to the analyses (default: FALSE)
#' @param jn boolean variable indicating whether the Johnson-Neyman procedure should be calculated (only available for single moderator models)
#' @return 
#'   \item{total}{A list of class "mod2rm" containing:}
#'   \item{info}{A named number vector containing values for the number of moderators in the model (num_mods), the number of binary moderators (num_binary_mods), the sample size (sample_size), the method of moderation (method; 1 = additive, 2 = multiplicative), and whether the Johnson-Neyman procedure was run (jn)}
#'   \item{var_names}{A named character vector containing the name of the original dataframe (dataframe), the two outcome variables (y1,y2), and up to three moderators (mod1,mod2,mod3)}
#'   \item{res_mod}{A list including the results of a simple regression, regressing the difference between y1 and y2 on the moderator}
#'   \item{res_simple_y1}{A list including the results of a simple regression, regressing the y1 on the moderator}
#'   \item{res_simple_y2}{A list including the results of a simple regression, regressing the y2 on the moderator}
#'   \item{res_cond_eff}{A list including the results of an analysis of conditional effects at different levels of the moderator(s)}
#'   \item{res_y1y2_diff}{A list including the results of a repeated measures t-test for y1 and y2}
#'   \item{res_jn_area}{A list containing information on the Johnson-Neyman procedure, incluing the number of significance points identified within the data range (num_jn), the moderator values of these points, as well as the proportion of the sample scoring higher than these values (jn_values), and information on whether the JN region is significant or non-significent (center_significant; used for plotting.)}
#'   \item{res_jn_cond_eff}{A list containing additional conditonal effects at levels of the moderator around the JN region. Values span the entir data range in 20 steps.}
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
#' @export
#' @importFrom stats confint lm na.omit pt qt sd t.test vcov
#' @importFrom utils capture.output



mod2rm <- function (data, Y1, Y2, MOD1, MOD2 = NULL, MOD3 = NULL, MOD1val = NULL, MOD2val = NULL, MOD3val = NULL, method = 1, standardize = FALSE, jn = FALSE) 
{
  
  # test arguments: check if manually set values are of length 2, check if moderators are specified in order (i.e., no skipping mod2)
  
   if((length(MOD1val) < 2 & !missing(MOD1val)) | (length(MOD2val) < 2 & !missing(MOD2val)) | (length(MOD3val) < 2 & !missing(MOD3val))) 
     stop('Vectors of manually set values for conditional effects must have a length of at least 2')
  if(!missing(MOD3) & missing(MOD2)) 
    stop('You cannot specify moderator 3 (MOD3) without also specifiying moderator 2 (MOD2) ')
  
  # determine the total number of moderators.  
  
  num_mods = sum(c(!missing(MOD1), !missing(MOD2), !missing(MOD3)))
  
  # make sure JN is only performed on one moderator
  if (num_mods > 1 & jn == TRUE) 
    stop('Johnson-Neyman procedure is only supported for one moderator.')
  
  # testing if (each) moderator is binary. Give error if specific values are assigned to any binary moderator
  
  m1_bin = ifelse((length(unique(na.omit(eval(substitute(MOD1), data)))) < 3),TRUE, FALSE)  # detect if moderator 1 is dichotomous or not
  m2_bin = ifelse((length(unique(na.omit(eval(substitute(MOD2), data)))) < 3),TRUE, FALSE)  
  m3_bin = ifelse((length(unique(na.omit(eval(substitute(MOD3), data)))) < 3),TRUE, FALSE)  
   
  if((m1_bin == TRUE & !missing(MOD1val)) | (m2_bin == TRUE & !missing(MOD2val)) | (m3_bin == TRUE & !missing(MOD3val))) 
    stop('it is not possible to manually set conditional effects values for binary (or missing) moderators')
  
  
  # create dataframe (only including non-missing cases) and critical variables for use within the function, depending on number of moderators
  
  if(num_mods == 1) {
    df = na.omit(data.frame(y1 = eval(substitute(Y1), data), y2 = eval(substitute(Y2), data), mod1 = eval(substitute(MOD1), data)))
    }
  if(num_mods == 2) {
    df = na.omit(data.frame(y1 = eval(substitute(Y1), data), y2 = eval(substitute(Y2), data), mod1 = eval(substitute(MOD1), data), mod2 = eval(substitute(MOD2), data)))
  }
  if(num_mods == 3) {
    df = na.omit(data.frame(y1 = eval(substitute(Y1), data), y2 = eval(substitute(Y2), data), mod1 = eval(substitute(MOD1), data), mod2 = eval(substitute(MOD2), data), mod3 = eval(substitute(MOD3), data)))
  }

  if(standardize == FALSE){
    y1 = as.numeric(df$y1)
    y2 = as.numeric(df$y2)
    mod1 = as.numeric(df$mod1)
    if(num_mods > 1) 
      mod2 = as.numeric(df$mod2)
    if(num_mods > 2) {
      mod3 = as.numeric(df$mod3)
    }
  }
  if(standardize == TRUE){
    mod1 = scale(as.numeric(df$mod1))
    if(num_mods > 1) 
      mod2 = scale(as.numeric(df$mod2))
    if(num_mods > 2) {
      mod3 = scale(as.numeric(df$mod3))
    }
    
  }
  
  diff = y1 - y2
  
  #run regression, predicting difference score from moderator(s) (i.e., testing the interaction)

  if (num_mods == 1) erg = lm(diff ~ mod1)
  if (num_mods == 2) {
    if(method == 1) erg = lm(diff ~ mod1 + mod2)
    if(method == 2) erg = lm(diff ~ mod1*mod2)
  }
  if (num_mods == 3) {
    if(method == 1) erg = lm(diff ~ mod1 + mod2 + mod3)
    if(method == 2) erg = lm(diff ~ mod1*mod2*mod3)
  }


  # look at simple effects of moderator(s) on condition

  if (num_mods == 1) erg2a = lm(y1 ~ mod1)
  if (num_mods == 2) {
    if(method == 1) erg2a = lm(y1 ~ mod1 + mod2)
    if(method == 2) erg2a = lm(y1 ~ mod1*mod2)
  }
  if (num_mods == 3) {
    if(method == 1) erg2a = lm(y1 ~ mod1 + mod2 + mod3)
    if(method == 2) erg2a = lm(y1 ~ mod1*mod2*mod3)
  }


  if (num_mods == 1) erg2b = lm(y2 ~ mod1)
  if (num_mods == 2) {
    if(method == 1) erg2b = lm(y2 ~ mod1 + mod2)
    if(method == 2) erg2b = lm(y2 ~ mod1*mod2)
  }
  if (num_mods == 3) {
    if(method == 1) erg2b = lm(y2 ~ mod1 + mod2 + mod3)
    if(method == 2) erg2b = lm(y2 ~ mod1*mod2*mod3)
  }


  ### Conditional effects:
    
 
  # if testpoints were manually specified
  
  if(!missing(MOD1val)) { 
    m1_val = sort(unique(MOD1val))
  }
  
  if(!missing(MOD2val) & num_mods > 1) {  
    m2_val = sort(unique(MOD2val))
  }
  
  if(!missing(MOD3val) & num_mods > 2) {  
    m3_val = sort(unique(MOD3val))
  }
  
  #if nothing was manually specified, use mean +- 1 SD
  
  # Set values of moderator (first, -+ 1 SD)
  
  if(missing(MOD1val)) { 
    m1_val = c(mean(mod1) - sd(mod1), mean(mod1), mean(mod1) + sd(mod1))
  }
  
  if (num_mods > 1 & missing(MOD2val)) {
    m2_val = c(mean(mod2) - sd(mod2), mean(mod2), mean(mod2) + sd(mod2))
  }
  
  if (num_mods > 2 & missing(MOD3val)) {
    m3_val = c(mean(mod3) - sd(mod3), mean(mod3), mean(mod3) + sd(mod3))
  }
  
  if(m1_bin == TRUE) {                  # in case of binary moderator, select both scores that exist
    m1_val = c(sort(unique(mod1))[1],  sort(unique(mod1))[2])
  }
  
  if(num_mods > 1 & m2_bin == TRUE) {                            
    m2_val = c(sort(unique(mod2))[1],  sort(unique(mod2))[2])
  }
  
  if(num_mods > 2 & m3_bin == TRUE) {                            
    m3_val = c(sort(unique(mod3))[1],  sort(unique(mod3))[2])
  }
  
  
  # calculate number of binary moderators (num_bin), thus far only needed for the info output
  num_bin = sum(m1_bin)
  if(num_mods == 2) num_bin = sum(c(m1_bin, m2_bin))
  if(num_mods == 3) num_bin = sum(c(m1_bin, m2_bin, m3_bin))
  
  # create list of conditional effects
  
  if (num_mods == 1) val_list = do.call(expand.grid, list(m1_val)) # a, b, c = vectors of values
  if (num_mods == 2) val_list = do.call(expand.grid, list(m1_val, m2_val)) # a, b, c = vectors of values
  if (num_mods == 3) val_list = do.call(expand.grid, list(m1_val, m2_val, m3_val)) # a, b, c = vectors of values
  
  val_list = cbind(val_list, "Effect" = NA, "SE" = NA, "t-value" = NA, "p-value" = NA, "2.5 %" = NA, "97.5 %" = NA)
  
  for(i in 1:nrow(val_list)) {
    val1 = val_list[i,1]
    if (num_mods > 1) val2 = val_list[i,2]   # add up to three columns for test-values (for up to 3 moderators) to the conditional effects list
    if (num_mods > 2) val3 = val_list[i,3]

    ## ADDITIVE MODERATION
    
    if (method == 1) {

      ic_b = as.numeric(erg$coefficients[1])   # get unstandardized coefficient intercept
      m1_b = as.numeric(erg$coefficients[2])   # get unstandardized coefficient moderators
      
      if(num_mods == 1)           # calculate gradient at value(s)
        g = ic_b + (val1 * m1_b)   
      if(num_mods == 2) {
        m2_b = as.numeric(erg$coefficients[3]) 
        g = ic_b + (val1 * m1_b) + (val2 * m2_b)
      }
      if(num_mods == 3) {
        m2_b = as.numeric(erg$coefficients[3]) 
        m3_b = as.numeric(erg$coefficients[4]) 
        g = ic_b + (val1 * m1_b) + (val2 * m2_b) + (val3 * m3_b)
      }
      
      ic_var = vcov(erg)[1,1]              # get variances
      m1_var = vcov(erg)[2,2]              
      if (num_mods > 1) m2_var = vcov(erg)[3,3]             
      if (num_mods > 2) m3_var = vcov(erg)[4,4]
      
      ic_m1_covar = vcov(erg)[1,2]         # get covariances
      if (num_mods > 1) {
        ic_m2_covar = vcov(erg)[1,3]   
        m1_m2_covar = vcov(erg)[2,3]   
      }            
      if (num_mods > 2) {
        ic_m3_covar = vcov(erg)[1,4]   
        m2_m3_covar = vcov(erg)[3,4]   
        m1_m3_covar = vcov(erg)[2,4]   
      }

      #calculate standard error of gradient
      
      if(num_mods == 1)
        g_se = sqrt(ic_var + val1^2 * m1_var + 2*val1*ic_m1_covar)  # calculate standard error of gradient

      if(num_mods == 2)
        g_se = sqrt(ic_var +  
                    2*val1*ic_m1_covar +             # twice the covar of the intercept with the respective moderator, weighted by the moderator value
                    2*val2*ic_m2_covar + 
                    2*val1*val2*m1_m2_covar +   # twice the covar between the moderators, weighted by BOTH values
                    val1^2 * m1_var +                      # variances of both moderator1 and moderator2, each multiplied by the squared moderator values
                    val2^2 * m2_var  
               )
      
      if(num_mods == 3)
        g_se = sqrt(ic_var +  
                      2*val1*ic_m1_covar +             # twice the covar of the intercept with the respective moderator, weighted by the moderator value
                      2*val2*ic_m2_covar + 
                      2*val3*ic_m3_covar + 
                      
                      2*val1*val2*m1_m2_covar +   # twice the covar between the moderators, weighted by BOTH values
                      2*val1*val3*m1_m3_covar +   # twice the covar between the moderators, weighted by BOTH values
                      2*val2*val3*m2_m3_covar +   # twice the covar between the moderators, weighted by BOTH values
                      
                      val1^2 * m1_var +                      # variances of both moderator1 and moderator2, each multiplied by the squared moderator values
                      val2^2 * m2_var +
                      val3^2 * m3_var 
                )
    }
    
    ## MULTIPLICATIVE MODERATION
    
    if (method == 2) {
      
      ic_b = as.numeric(erg$coefficients[1])   # get unstandardized coefficient intercept
      m1_b = as.numeric(erg$coefficients[2])   # get unstandardized coefficient moderator
    
      if(num_mods == 1)           # calculate gradient at value(s)
        g = ic_b + (val1 * m1_b)   
      if(num_mods == 2) {
        m2_b = as.numeric(erg$coefficients[3]) 
        m1m2_b = as.numeric(erg$coefficients[4]) 
        val1x2 = val1 * val2                                        #set test value of interaction(s) (val1 * val2)
        g = ic_b + (val1 * m1_b) + (val2 * m2_b) + (val1x2 * m1m2_b)
      }
      
      if(num_mods == 3) {
        
        val1x2 = val1 * val2
        val1x3 = val1 * val3
        val2x3 = val2 * val3
        val1x2x3 = val1*val2*val3
        
        m2_b = as.numeric(erg$coefficients[3])   # get unstandardized coefficients of moderators
        m3_b = as.numeric(erg$coefficients[4])   
        m1m2_b = as.numeric(erg$coefficients[5])
        m1m3_b = as.numeric(erg$coefficients[6])
        m2m3_b = as.numeric(erg$coefficients[7])
        m1m2m3_b = as.numeric(erg$coefficients[8])
        
        g = ic_b + (val1 * m1_b) + (val2 * m2_b) + (val3 * m3_b) + (val1x2 * m1m2_b) + (val1x3 * m1m3_b) + (val2x3 * m2m3_b) +  (val1x2x3 * m1m2m3_b)   # calculate gradient val1
        
      }

      ic_var = vcov(erg)[1,1]              # get variances
      m1_var = vcov(erg)[2,2]              
      if (num_mods == 2) {
        m2_var = vcov(erg)[3,3]
        m1m2_var = vcov(erg)[4,4]
      }
      if (num_mods == 3) {
        m2_var = vcov(erg)[3,3]
        m3_var = vcov(erg)[4,4]
        m1m2_var = vcov(erg)[5,5]             
        m1m3_var = vcov(erg)[6,6]              
        m2m3_var = vcov(erg)[7,7]             
        m1m2m3_var = vcov(erg)[8,8]             
      }
      
      
      ic_m1_covar = vcov(erg)[1,2]         # get covariances
      if (num_mods == 2) {
        ic_m1_covar  = vcov(erg)[1,2]   
        ic_m2_covar = vcov(erg)[1,3]   
        ic_m1m2_covar = vcov(erg)[1,4]   
        m1_m2_covar = vcov(erg)[2,3]    
        m1_m1m2_covar = vcov(erg)[2,4] 
        m2_m1m2_covar = vcov(erg)[3,4] 
      }            
      
      if (num_mods == 3) {
        ic_m1_covar  = vcov(erg)[1,2]   
        ic_m2_covar = vcov(erg)[1,3]   
        ic_m3_covar = vcov(erg)[1,4]   
        ic_m1m2_covar  = vcov(erg)[1,5]   
        ic_m1m3_covar = vcov(erg)[1,6]   
        ic_m2m3_covar = vcov(erg)[1,7]   
        ic_m1m2m3_covar = vcov(erg)[1,8]   
        
        m1_m2_covar = vcov(erg)[2,3]   
        m1_m3_covar = vcov(erg)[2,4]   
        m1_m1m2_covar  = vcov(erg)[2,5]   
        m1_m1m3_covar = vcov(erg)[2,6]   
        m1_m2m3_covar = vcov(erg)[2,7]   
        m1_m1m2m3_covar = vcov(erg)[2,8]   
        
        m2_m3_covar = vcov(erg)[3,4]   
        m2_m1m2_covar  = vcov(erg)[3,5]   
        m2_m1m3_covar = vcov(erg)[3,6]   
        m2_m2m3_covar = vcov(erg)[3,7]   
        m2_m1m2m3_covar = vcov(erg)[3,8]   
        
        m3_m1m2_covar  = vcov(erg)[4,5]   
        m3_m1m3_covar = vcov(erg)[4,6]   
        m3_m2m3_covar = vcov(erg)[4,7]   
        m3_m1m2m3_covar = vcov(erg)[4,8]   

        m1m2_m1m3_covar = vcov(erg)[5,6]   
        m1m2_m2m3_covar = vcov(erg)[5,7]   
        m1m2_m1m2m3_covar = vcov(erg)[5,8]   
        
        m1m3_m2m3_covar = vcov(erg)[6,7]   
        m1m3_m1m2m3_covar = vcov(erg)[6,8]   
        
        m2m3_m1m2m3_covar = vcov(erg)[7,8]   
        
      }

      #calculate standard error of gradient
      
      if(num_mods == 1)
        g_se = sqrt(ic_var + val1^2 * m1_var + 2*val1*ic_m1_covar)  # calculate standard error of gradient
      
      if(num_mods == 2)
        g_se = sqrt(ic_var +  
                      2*val1*ic_m1_covar +             # twice the covar of the intercept with the respective moderator, weighted by the moderator value
                      2*val2*ic_m2_covar + 
                      2*val1x2*ic_m1m2_covar + 
                      
                      2*val1*val2*m1_m2_covar +   # twice the covar between the moderators, weighted by BOTH values
                      2*val1*val1x2*m1_m1m2_covar +   # twice the covar between the moderators, weighted by BOTH values
                      2*val2*val1x2*m2_m1m2_covar +   # twice the covar between the moderators, weighted by BOTH values
                      
                      val1^2 * m1_var +                      # variances of both moderator1 and moderator2, each multiplied by the squared moderator values
                      val2^2 * m2_var +
                      val1x2^2 * m1m2_var 
                    )
      
      if(num_mods == 3)
        g_se = sqrt(ic_var +  
                      2*val1*ic_m1_covar +             
                      2*val2*ic_m2_covar +             
                      2*val3*ic_m3_covar +             
                      
                      2*val1x2*ic_m1m2_covar + 
                      2*val1x3*ic_m1m3_covar + 
                      2*val2x3*ic_m2m3_covar + 
                      2*val1x2x3*ic_m1m2m3_covar + 
                      
                      2*val1*val2*m1_m2_covar + 
                      2*val1*val3*m1_m3_covar + 
                      2*val2*val3*m2_m3_covar + 
                      
                      2*val1*val1x2*m1_m1m2_covar +
                      2*val1*val1x3*m1_m1m3_covar +
                      2*val1*val2x3*m1_m2m3_covar +
                      2*val1*val1x2x3*m1_m1m2m3_covar +
                      
                      2*val2*val1x2*m2_m1m2_covar +
                      2*val2*val1x3*m2_m1m3_covar +
                      2*val2*val2x3*m2_m2m3_covar +
                      2*val2*val1x2x3*m2_m1m2m3_covar +
                      
                      2*val3*val1x2*m3_m1m2_covar +
                      2*val3*val1x3*m3_m1m3_covar +
                      2*val3*val2x3*m3_m2m3_covar +
                      2*val3*val1x2x3*m3_m1m2m3_covar +
                      
                      
                      2*val1x2*val1x3*m1m2_m1m3_covar +
                      2*val1x2*val2x3*m1m2_m2m3_covar +
                      2*val1x2*val1x2x3*m1m2_m1m2m3_covar +
                      
                      2*val1x3*val2x3*m1m3_m2m3_covar +
                      2*val1x3*val1x2x3*m1m3_m1m2m3_covar +
                      
                      2*val2x3*val1x2x3*m2m3_m1m2m3_covar +
  
                      val1^2 * m1_var +                     
                      val2^2 * m2_var +
                      val3^2 * m3_var +
                      val1x2^2 * m1m2_var +
                      val1x3^2 * m1m3_var +
                      val2x3^2 * m2m3_var +
                      val1x2x3^2 * m1m2m3_var
              )
    }
    
    t = g / g_se                         #calculate t-value of slopes

    # look up p value from that t value and the degrees of freedom (n minus 2, for one moderator)

    p = 2*pt(q=abs(t), df=nrow(df) - length(erg$coefficients), lower.tail=FALSE)     # look up p value from that t value and the degrees of freedom (n minus 2)

    # calculate 95% confidence intervals of the conditional effect

    margin <- qt(p = 0.05 / 2, df = nrow(df) - length(erg$coefficients), lower.tail=F) * g_se

    ci_lower = g - margin
    ci_upper = g + margin
    
    # name the moderator variables first
    colnames(val_list)[1] = "mod1"
    if(num_mods > 1) colnames(val_list)[2] = "mod2"
    if(num_mods > 2) colnames(val_list)[3] = "mod3"
    
    val_list[i,"Effect"] = g
    val_list[i,"SE"] = g_se
    val_list[i,"t-value"] = t
    val_list[i,"p-value"] = p
    val_list[i,"2.5 %"] = ci_lower
    val_list[i,"97.5 %"] = ci_upper
  
  }

  ## Johnson-Neyman Procedure for one moderator
  
  if(jn == TRUE) {
    
    # define some variables
    
    # calculate critical t-value
    
    critt = qt(p=.05/2, df=nrow(df) - length(erg$coefficients), lower.tail=FALSE)
    
    p3 = (ic_b^2) - (critt^2)*ic_var 
    p2 = 2*ic_b*m1_b - 2*(critt^2)*ic_m1_covar 
    p1 = (m1_b^2) - (critt^2 * m1_var) 
    jn_values = data.frame(values = NA, percent = NA)

    if ((p2^2 - 4*p3*p1) >= 0) {     # (only calculate JN values if applicable)

      
      # Moderator values defining the JN significance region(s) are:
      
      jn_val1 = (-1*p2 + sqrt(p2^2 - 4*p3*p1)) /  (2*p1) 
      jn_val2 = (-1*p2 - sqrt(p2^2 - 4*p3*p1)) /  (2*p1) 
      
      # calculate percentage of sample above values 
      
      perc_jn_val1 = round(nrow(df[mod1 > jn_val1,]) / nrow(df) * 100, digits = 2)
      perc_jn_val2 = round(nrow(df[mod1 > jn_val2,]) / nrow(df) * 100, digits = 2)  
      
      
      # determine numer of JN significance points
      
      num_jn = 2
      
      if (jn_val1 < min(mod1) | jn_val1 > max(mod1)) {
        num_jn = num_jn - 1
      }
      
      if (jn_val2 < min(mod1) | jn_val2 > max(mod1)) {
        num_jn = num_jn -1
      }
      
      #determine, whether the area between the two JN intervals is significant or not (important for plot)
      center_significant = ifelse(jn_val1 < jn_val2, TRUE, FALSE)
      
      # save everything to return vectors
      jn_values[1,1] = min(c(jn_val1, jn_val2))
      jn_values[2,1] = max(c(jn_val1, jn_val2))
      jn_values[1,2] = max(c(perc_jn_val1, perc_jn_val2))
      jn_values[2,2] = min(c(perc_jn_val1, perc_jn_val2))
      
    } 
    
    if ((p2^2 - 4*p3*p1) < 0) {
      num_jn = 0
      jn_values = NA
      center_significant = NA
    }
   

    # calculate conditional effects for JN procedure
    
    mod1_range = max(mod1) - min(mod1) # define ranger of moderator values
    jn_val = seq(min(mod1), max(mod1), by = mod1_range / 20) # ... and create 20 values spanning from lowest to highest moderator value
    
    jn_val_list = do.call(expand.grid, list(jn_val)) 
    jn_val_list = cbind(jn_val_list, "Effect" = NA, "SE" = NA, "t-value" = NA, "p-value" = NA, "2.5 %" = NA, "97.5 %" = NA)
    
    for(i in 1:nrow(jn_val_list)) {
      jn_val1 = jn_val_list[i,1]
      
      ic_b = as.numeric(erg$coefficients[1])   # get unstandardized coefficient intercept
      m1_b = as.numeric(erg$coefficients[2])   # get unstandardized coefficient moderators
      
      g = ic_b + (jn_val1 * m1_b)   
      
      ic_var = vcov(erg)[1,1]              # get variances
      m1_var = vcov(erg)[2,2]              
      
      ic_m1_covar = vcov(erg)[1,2]         # get covariances
      
      g_se = sqrt(ic_var + jn_val1^2 * m1_var + 2*jn_val1*ic_m1_covar)  # calculate standard error of gradient
      
      t = g / g_se                         #calculate t-value of slopes
      
      # look up p value from that t value and the degrees of freedom (n minus 2, for one moderator)
      
      p = 2*pt(q=abs(t), df=nrow(df) - length(erg$coefficients), lower.tail=FALSE)    
      
      # calculate 95% confidence intervals of the conditional effect
      
      margin <- qt(p = 0.05 / 2, df = nrow(df) - length(erg$coefficients), lower.tail=F) * g_se
      
      ci_lower = g - margin
      ci_upper = g + margin
      
      # name the moderator variables first
      colnames(jn_val_list)[1] = "mod1"
      if(num_mods > 1) colnames(jn_val_list)[2] = "mod2"
      if(num_mods > 2) colnames(jn_val_list)[3] = "mod3"
      
      jn_val_list[i,"Effect"] = g
      jn_val_list[i,"SE"] = g_se
      jn_val_list[i,"t-value"] = t
      jn_val_list[i,"p-value"] = p
      jn_val_list[i,"2.5 %"] = ci_lower
      jn_val_list[i,"97.5 %"] = ci_upper
    }
    
  }
  if(jn == FALSE) {
    num_jn = NA
    jn_values = NA
    jn_val_list = NA
    center_significant = NA
  }
  
  # t-test difference Y1 and Y2

  diff_ttest = t.test(y1, y2, paired = T)

  # save everything in new object
  
  names = c(dataframe = paste(deparse(substitute(data))), y1 = paste(deparse(substitute(Y1))), y2 = paste(deparse(substitute(Y2))), mod1 = paste(deparse(substitute(MOD1))), mod2 = paste(deparse(substitute(MOD2))), mod3 = paste(deparse(substitute(MOD3))))
  numbers = c(num_mods = num_mods, num_binary_mods = num_bin, sample_size = nrow(df), method = method, standardize = standardize, jn = jn)
  
  results = list(info = numbers, 
                 var_names = names, 
                 res_mod = erg, 
                 res_simple_y1 = erg2a, 
                 res_simple_y2 = erg2b, 
                 res_cond_eff = val_list, 
                 res_y1y2_diff = diff_ttest,
                 res_jn_area = list(num_jn = num_jn, jn_values = jn_values, center_significant = center_significant),
                 res_jn_cond_eff = jn_val_list
    )
                 
  class(results) <- c("mod2rm") # set class of return object

  return(results)
  
}

