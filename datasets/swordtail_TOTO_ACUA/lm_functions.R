library(data.table)

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 


# 2. ----- r squared Functions ------

# From linear model
r.squared_lm <- function(data, z, formula1, formula2){
  # assumes data has coefficients for only a single scale
  fitData <- data[!is.na(z)]
  
  # fit models with intercept and either one or two predictors
  fit1 <- lm(formula1, fitData)
  r.squared_1 <- summary(fit1)$adj.r.squared
  anova1 <- anova(fit1)
  p1 <- anova1$`Pr(>F)`[1]
  
  fit2 <- lm(formula2, fitData)
  r.squared_2 <- summary(fit2)$adj.r.squared
  anova2 <- anova(fit2)
  p2 <- anova2$`Pr(>F)`[1:2]
  
  return(list(model = c("x","xy"),
              r.squared_adj = c(r.squared_1, r.squared_2),
              p.val_rec = c(p1, p2[1]),
              p.val_cds = c(NA, p2[2])))
}


# ===== bootstrap confidence intervals =====

r.sqrd_boot_single_level <- function(data, indices, formula1, formula2){
  d <- data[indices,] # allows boot to select sample
  fit1 <- lm(formula1, data=d)
  fit2 <- lm(formula2, data=d)
  return(c(summary(fit1)$adj.r.squared, summary(fit2)$adj.r.squared))
}

r.sqrd_boot_all <- function(data){
  # runs bootstrap sampling on single level and outputs
  # to be applied to data table grouped by level for output per level 
  
  d <- data[!is.na(ancestry_coeff)]
  
  result <- boot(d, statistic=r.sqrd_boot_single_level, 
                 formula1=ancestry_coeff~rec_coeff,
                 formula2=ancestry_coeff~rec_coeff+cds_coeff,
                 R=1000,
                 parallel="multicore",
                 ncpus=n.cores)
  
  x_CI <- boot.ci(result, type="perc", index=1) # recomb
  xy_CI <- boot.ci(result, type="perc", index=2) # recomb + gene density
  
  return(data.table(model = c("x","xy"),
                    CI_l = c(x_CI$bca[4],xy_CI$bca[4]),
                    CI_u = c(x_CI$bca[5],xy_CI$bca[5]))
  )
  
}

