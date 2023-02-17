setwd("~/workspace/selection-against-introgression/theory_and_simulations/var_among_individuals/")

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]}   

wv <- rbindlist(
  lapply(list.files(path = "results/neutral_sims/equilibrium/", pattern = "*wavelet_results.RData", full.names=T), function(x){
    return(loadFrom(x, "true_ancestry_allWV"))
  })
)


ind_gnom_avg <- rbindlist(
  lapply(list.files(path = "results/neutral_sims/equilibrium/", pattern = "*wavelet_results.RData", full.names=T), function(x){
    return(loadFrom(x, "ind_gnom_avg"))
  })
)



# -----look at among individual variance
v_among <- ind_gnom_avg[, .(v = var(gnom_avg)), by = .(rep, alpha, alpha_empirical, gen)]

# mean and standard error of among-individual variance
v_among_mn <- v_among[, .(v_mn = mean(v),
                          v_se = sd(v)/sqrt(nrow(.SD)),
                          alpha_empirical = mean(alpha_empirical)), 
                      by = .(gen, alpha)]


#ggplot(v_among, aes(x = gen, y = v, group = rep)) + geom_line() + facet_wrap(~alpha)
ggplot(v_among, aes(x = log10(gen), y = v, group = rep)) + geom_line() + facet_wrap(~alpha)


# ===== variance along sequence =====


v_along <- wv[signal == 'single_hap', .(v = sum(variance)), 
              by = .(rep, alpha, alpha_empirical, gen)]



ggplot(v_along,
       aes(x =gen)) +
  geom_line(aes(group = rep, y = v)) + facet_wrap(~alpha) 

# mean and se
v_along_mn <- v_along[, .(v_mn = mean(v), 
                          v_among_pred1 = alpha - alpha^2 - mean(v),
                          v_among_pred2 = mean(alpha_empirical)- mean(alpha_empirical)^2 - mean(v),
                          v_among_pred_se1 = sd(alpha-alpha^2 - v)/sqrt(nrow(.SD)), 
                          v_among_pred_se2 = sd(mean(alpha_empirical)-mean(alpha_empirical)^2 - v)/sqrt(nrow(.SD)),
                          alpha_empirical = mean(alpha_empirical)), by = .(gen, alpha)]




# ----- plot

# preliminary
ggplot(v_along_mn, aes(x = log10(gen+1))) + 
  facet_wrap(~alpha, scales = 'free_y') + 
  geom_point(aes(y = v_among_pred1)) +
  geom_errorbar(aes(ymin = v_among_pred2-1.96*v_among_pred_se2, ymax = v_among_pred2+1.96*v_among_pred_se2)) +
  geom_line(aes(y = v_among_pred2)) +
  geom_point(data = v_among_mn, aes(x = log10(gen+1), y = v_mn), color = 'red') + 
  geom_line(data = v_among_mn, aes(x = log10(gen+1), y = v_mn), color = 'red') + 
  geom_errorbar(data = v_among_mn, aes(ymin=v_mn-1.96*v_se, ymax=v_mn+1.96*v_se), color = 'red') + 
  #theme_classic() + 
  theme(aspect.ratio = 1, legend.title = element_text(size = 20)) +
  labs(y = "Variance among individuals", color = "Legend") + 
  scale_color_manual(values = c('true'="red", 'predicted' = 'black'))


plot_data <- rbind(v_among_mn[, .(gen, alpha, v_among = v_mn, v_among_se = v_se, type = "true")], 
      v_along_mn[, .(gen, alpha, v_among = v_among_pred2, v_among_se = v_among_pred_se2, type = "predicted")])

ggplot(plot_data, aes(x = log10(gen+1), y = v_among, group = type, color = type)) + 
  geom_point() + geom_line() + 
  facet_wrap(~alpha, scales = 'free_y') + 
  geom_errorbar(aes(ymin=v_among-1.96*v_among_se, ymax = v_among+1.96*v_among_se)) + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  labs(y = "Variance among individuals", 
       color = "") + 
  scale_color_manual(values = c("darkred", "darkgrey"))


