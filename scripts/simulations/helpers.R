library(tidyverse) 

### Inverse Logit Function
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

### Logit Function
logit <- function(x) {
  return(log(x/(1-x)))
}

### Quick way to compute logistic regression coefficients when fitting saturated model with one coefficient
glm_quick <- function(Y, A, W) {
  beta0 <- logit(weighted.mean(Y[A == 0], w = W[A == 0])) 
  beta1 <- logit(weighted.mean(Y[A == 1], w = W[A == 1])) - beta0
  return(list('beta0' = beta0, 'beta1' = beta1))
}

### Function to compute the fitted value of a model given a list of coefficients
compute_model <- function(df, beta) {
  ### Initialize fitted values to be 0
  preds <- rep(0, nrow(df))
  
  ### Iterate over coefficients
  for(i in 1:length(beta)) {
    coeff <- names(beta)[i]
    
    if(coeff == '(Intercept)') { ### Intercept
      preds <- preds + beta[[coeff]]
    } else if(grepl('\\[', coeff) & grepl(':', coeff)) { ### Level of categorical variable w/ interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeffs[1]))
      variable <- gsub('\\[.*$', '', coeffs[1])
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + df[[ coeffs[2] ]][ix] * beta[[ coeff ]]
    } else if(grepl('\\[', coeff)) { ### Level of categorical variable
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeff))
      variable <- gsub('\\[.*$', '', coeff)
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + beta[[coeff]]
    } else if(grepl(':', coeff)) { ### Interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      preds <- preds + df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * beta[[coeff]]
    } else if(grepl('\\^2\\)$', coeff)) { ### quadratic 
      variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
      preds <- preds + df[[ variable ]]^2 * beta[[coeff]]
    } else if(grepl('^I\\(exp\\(', coeff)) { ### exponential
      variable <- gsub('\\)\\)$', '', gsub('^I\\(exp\\(', '', coeff))
      preds <- preds + exp(df[[ variable ]]) * beta[[coeff]]
    } else { ### Linear Term 
      preds <- preds + df[[coeff]] * beta[[coeff]]
    }
  }
  
  return(preds)
}

### Function to compute the matrix of fitted group probabilities for a 
### multinomial given a list of coefficients
compute_multinomial_model <- function(df, beta) {
  ### Initialize to logit of cumumlative probability of being in each group
  probs <- matrix(beta$intercepts, nrow = nrow(df), ncol = length(beta$intercepts), byrow = T)
  
  ### Add in rest of covariate effects
  for(i in 1:length(beta)) {
    coeff <- names(beta)[i]
    
    ### Note that positive beta --> subtraction 
    ### See: https://www.bookdown.org/rwnahhas/RMPH/blr-ordinal.html
    if(!(coeff %in% c('intercepts', '(Intercept)'))) {
      if(grepl(':', coeff)) { ### Interaction
        coeffs <- unlist(strsplit(coeff, ':'))
        probs <- probs - df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * beta[[coeff]]
      } else if(grepl('\\[', coeff) & grepl(':', coeff)) { ### Level of categorical variable w/ interaction
        coeffs <- unlist(strsplit(coeff, ':'))
        level <- gsub('^.*\\[', '',  gsub('\\]', '', coeffs[1]))
        variable <- gsub('\\[.*$', '', coeffs[1])
        ix <- df[[variable]] == level
        probs[ix,] <- probs[ix] + df[[ coeffs[2] ]][ix] * beta[[ coeff ]]
      } else if(grepl('\\[', coeff)) { ### Level of categorical variable
        level <- gsub('^.*\\[', '',  gsub('\\]', '', coeff))
        variable <- gsub('\\[.*$', '', coeff)
        ix <- df[[variable]] == level
        probs[ix,] <- probs[ix,] - beta[[coeff]]
      } else if(grepl('\\^2\\)$', coeff)) { ### quadratic 
        variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
        probs <- probs - df[[ variable ]]^2 * beta[[coeff]]
      } else if(grepl('^I\\(exp\\(', coeff)) { ### exponential
        variable <- gsub('\\)\\)$', '', gsub('^I\\(exp\\(', '', coeff))
        probs <- probs - exp(df[[ variable ]]) * beta[[coeff]]
      } else { ### Linear Term 
        probs <- probs - df[[coeff]] * beta[[coeff]]
      }
    }
  }
  
  return(probs)
}

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 12),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))


### Clean names so we know what variables to keep in our simulations
clean_names <- function(params, df_names) {
  x <- unique(gsub('^.*\\.', '', names(unlist(params))))
  x <- gsub('\\[.*', '', unlist(strsplit(x, ':')))
  x <- c(x, unique(unlist(map(params$truth_model, as.character))))
  x <- x[x %in% df_names]
  return(x)
}

### Function to replace null values w/ useful string
replace_null <- function(x, replace_string = NA) {
  x[is.null(x)] <- replace_string
  return(x)
}
