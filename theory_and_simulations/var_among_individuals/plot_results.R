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
                          v_se = sd(v)/sqrt(nrow(.SD))), 
                      by = .(gen, alpha)]


#ggplot(v_among, aes(x = gen, y = v, group = rep)) + geom_line() + facet_wrap(~alpha)
ggplot(v_among, aes(x = log10(gen), y = v, group = rep)) + geom_line() + facet_wrap(~alpha)


# ===== variance along sequence =====
v_along <- wv[signal == 'single_hap', .(v = sum(variance)), 
              by = .(rep, alpha, alpha_empirical, gen)]

# ggplot(v_along, 
#        aes(x =gen, y = v)) + 
#   geom_line(aes(group = rep)) + facet_wrap(~alpha)

# mean and se
v_along_mn <- v_along[, .(v_mn = mean(v), 
            v_se = sd(alpha-alpha^2 - v)/sqrt(nrow(.SD)), 
            alpha_empirical = mean(alpha_empirical)), by = .(gen, alpha)]

# ----- plot
ggplot(v_along_mn, aes(x = log10(gen+1), y =  alpha-alpha^2 - v_mn)) + 
  facet_wrap(~alpha) + geom_point() +
  geom_errorbar(aes(ymin = (alpha-alpha^2-v_mn)-1.96*v_se, ymax = (alpha-alpha^2-v_mn)+1.96*v_se)) +
  geom_line() +
  geom_point(data = v_among_mn, aes(x = log10(gen+1), y = v_mn), color = 'red') + 
  geom_line(data = v_among_mn, aes(x = log10(gen+1), y = v_mn), color = 'red') + 
  geom_errorbar(data = v_among_mn, aes(ymin=v_mn-1.96*v_se, ymax=v_mn+1.96*v_se), color = 'red')



