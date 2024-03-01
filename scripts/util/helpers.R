library(tidyverse)

### Function to strip model object of everything except ability to do 
### prediction (which is useful for reducing its size on large datasets)
### https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/
strip_model <- function(model) {
  model$y <- NULL
  model$model <- NULL
  
  model$residuals <- NULL
  model$fitted.values <- NULL
  model$effects <- NULL
  model$qr$qr <- NULL  
  model$linear.predictors <- NULL
  model$weights <- NULL
  model$prior.weights <- NULL
  model$data <- NULL
  
  
  model$family$variance <- NULL
  model$family$dev.resids <- NULL
  model$family$aic <- NULL
  model$family$validmu <- NULL
  model$family$simulate <- NULL
  attr(model$terms,".Environment") <- NULL
  attr(model$formula,".Environment") <- NULL
  
  return(model)
}

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 12),
                  strip.text.y = element_text(size = 8),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

winsorize <- function(x, q) {
  ### Get upper and lower quantiles of the data
  q_data <- quantile(x[!is.na(x)], q)
  lower <- q_data[1]
  upper <- q_data[2]
  
  ### Replace extreme values by quantiles
  x[which(x < lower)] <- lower
  x[which(x > upper)] <- upper
  
  return(x)
  
}

outlier_weight <- function(wt, date, n_sd, weight_tol) {
  ### Don't even bother smoothing on too few obs
  if(length(wt) <= 5) {
    return(rep(F, length(wt)))
  } else if(n_distinct(date) <= 5) {
    return(rep(F, length(wt)))
  }
  
  
  ### finer-grained smoothing if we have more
  s <- ifelse(length(wt) > 10, 0.5, 0.75) 
  
  date <- as.numeric(date)
  smoother <- loess(wt ~ date, span = s)
  wt_smooth <- tryCatch(predict(smoother, date, se = T))
  
  if(class(wt_smooth) != 'try-error') {
    upper_bound <- wt_smooth$fit + pmax(n_sd * wt_smooth$se.fit, weight_tol, na.rm = T)
    lower_bound <- wt_smooth$fit - pmax(n_sd * wt_smooth$se.fit, weight_tol, na.rm = T)
    
    ix <- wt < lower_bound | wt > upper_bound
    
  } else {
    return(rep(F, length(wt)))
  }
  
  return(ix)
}

### Method to save out survplots
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}