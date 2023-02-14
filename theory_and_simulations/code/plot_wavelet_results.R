library(data.table)
library(ggplot2)
library(cubature)
library(magrittr)
source('~/workspace/gnomwav/R/theory.R')

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# ===== Load Results =====

# ---- wavelet variance of true ancestry 
eqAncestryWV <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/', 
                    pattern="*wavelet_results.RData", full.names = T),
         function(x){
           wv <- loadFrom(x, "true_ancestry_allWV")
           wv
         }
  )
)


bnAncestryWV <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/bottleneck/', 
                    pattern="*wavelet_results.RData", full.names = T),
         function(x){
           wv <- loadFrom(x, "true_ancestry_allWV")
           wv
         }
  )
)

eqAncestryWV[, model := "equilibrium"]
bnAncestryWV[, model := 'bottleneck']
allAncestryWV <- rbind(eqAncestryWV, bnAncestryWV)


# ----- wavelet variance of snp statistic
eqSNPWV <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/', 
                    pattern="*wavelet_results.RData",
                    full.names = T),
         function(x){
           wv <- loadFrom(x, "allSNPWV")
           wv
           }
         )
)

bnSNPWV <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/bottleneck/', 
                    pattern="*wavelet_results.RData",
                    full.names=T),
         function(x){
           wv <- loadFrom(x, "allSNPWV")
           wv
         }
  )
)

eqSNPWV[, model := "equilibrium"]
bnSNPWV[, model := 'bottleneck']
allSNPWV <- rbind(eqSNPWV, bnSNPWV)


# ----- wavelet covariance of true ancestry
eq_trueWC <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/',
                    pattern="*wavelet_results.RData",
                    full.names = T),
         function(x){
           wc <- loadFrom(x, "true_ancestry_wc")
           wc
         }
  )
)


bn_trueWC <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/bottleneck/',
                    pattern="*wavelet_results.RData",
                    full.names = T),
         function(x){
           wc <- loadFrom(x, "true_ancestry_wc")
           wc
         }
  )
)

# wavelet covariance of SNP stat,
eq_snpWC <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/',
                    pattern="*wavelet_results.RData",
                    full.names = T),
         function(x){
           wc <- loadFrom(x, "all_snp_wc")
           wc
         }
  )
)

bn_snpWC <- rbindlist( 
  lapply(list.files(path='~/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/bottleneck/',
                    pattern="*wavelet_results.RData",
                    full.names = T),
         function(x){
           wc <- loadFrom(x, "all_snp_wc")
           wc
         }
  )
)

eq_trueWC[, model := "equilibrium"][, signal := 'true_ancestry'][]
bn_trueWC[, model := 'bottleneck'][, signal := 'true_ancestry'][]

all_trueWC <- rbind(eq_trueWC, bn_trueWC)

eq_snpWC[, model := 'equilibrium'][, signal := 'snp_stat']
bn_snpWC[, model := 'bottleneck'][, signal := 'snp_stat']

all_snpWC <- rbind(eq_snpWC, bn_snpWC)


# ===== Theoretical Expectations for Wavelet Variance ======
eq_thry <- wavelet_variance_equilbrium(n.pop = 20000, n.sample = c(1,20,20000), unit.dist = 2^-14, level = 1:14, gen = c(1,10,100,1000),alpha = 0.5)
eq_thry[, propvar := variance/sum(variance), by = .(n.sample, gen)]

bn_thry <- wavelet_variance_general(n.pop = c(20,20000), epochs = c(10,990), n.sample = c(1,20,20000), unit.dist = 2^-14, level = 1:14, gen = c(1,10,100,1000), alpha = 0.5)
bn_thry[, propvar := variance/sum(variance), by = .(n.sample, gen)]

# ====== True Ancestry Wavelet Variance ===== 

# ----- equilibrium: single haplotype

allAncestryWV[gen == 100 & model == 'equilibrium' & signal == 'single_hap' & ancestry_measure == 'direct', ] %>%
  ggplot(., aes(x = level, y = variance)) + 
  facet_wrap(~gen) + 
  geom_point() + geom_line(aes(group = rep)) + 
  geom_point(data=eq_thry[gen==100&n.sample==1], aes(x = level, y = variance), color = 'red') + 
  geom_line(data=eq_thry[gen==100&n.sample==1], aes(x = level, y = variance), color = 'red',size=1.2) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') +
  theme_classic()

# ----- equilibrium: sample mean

allAncestryWV[gen != 0 & model == 'equilibrium' & signal == 'smpl_mean' ] %>%
  ggplot(., aes(x = level, y = variance)) + 
  facet_wrap(~gen) + 
  geom_point() + geom_line(aes(group = rep)) + 
  geom_point(data=eq_thry[n.sample==20], aes(x = level, y = variance), color = 'red') + 
  geom_line(data=eq_thry[n.sample==20], aes(x = level, y = variance), color = 'red')+ 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') 

# ----- bottleneck: single haplotype
allAncestryWV[gen != 0 & model == 'bottleneck' & signal == 'single_hap' & ancestry_measure == 'direct', ] %>%
  ggplot(., aes(x = level, y = variance)) + 
  facet_wrap(~gen) + 
  geom_point() + geom_line(aes(group = rep)) + 
  geom_point(data=bn_thry[n.sample==1], aes(x = level, y = variance), color = 'red') + 
  geom_line(data=bn_thry[n.sample==1], aes(x = level, y = variance), color = 'red') +
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') 

allAncestryWV[gen != 0 & model == 'bottleneck' & signal == 'smpl_mean' ] %>%
  ggplot(., aes(x = level, y = variance)) + 
  facet_wrap(~gen) + 
  geom_point() + geom_line(aes(group = rep)) + 
  geom_point(data=bn_thry[n.sample==20], aes(x = level, y = variance), color = 'red') + 
  geom_line(data=bn_thry[n.sample==20], aes(x = level, y = variance), color = 'red') + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') 


# ===== SNP Stat single haplotype wavelet variance =====

# raw variance before correction
ggplot(allSNPWV[gen == 100 & ascertainment == 'gen0' & signal == 'single_hap' & pop=='p2'&model == 'equilibrium'], aes(x = level, y = variance)) + 
  geom_point() + geom_line(aes(group = rep)) +
  facet_wrap(~gen) + 
  theme(aspect.ratio = 1) +
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') + 
  theme_classic()+
  geom_point(data=eq_thry[gen==100&n.sample==1], aes(x = level, y = variance), color = 'red') + 
  geom_line(data=eq_thry[gen==100&n.sample==1], aes(x = level, y = variance), color = 'red',size=1.2) 
  

# ----- correction
h_wv_pops <- dcast(allSNPWV[signal=='single_hap'], rep+gen+level+ascertainment+ model ~ pop, value.var='variance')
h_wv_pops[, correction := p2 - sum(p0, p1)/2, by = .(rep, gen, level, ascertainment, model)]
h_wv_pops[, correction_propvar := correction/sum(correction), by = .(rep, gen, ascertainment, model)]

# average propvar 
h_mean_wv_pops <- h_wv_pops[, .(mean_correction_propvar = mean(correction_propvar), mean_correction = mean(correction), se_correction_propvar = sd(correction_propvar)/sqrt(10)), by= .(gen, level, ascertainment, model)]

# or take propvar of average
#h_mean_wv_pops2 <- h_wv_pops[, .(mean_correction = mean(correction)), by = .(gen, level, ascertainment, model)]
#h_mean_wv_pops2[, correction_propvar := mean_correction/sum(mean_correction), by = .(gen, ascertainment, model)]

# equilibrium
#ggplot(h_mean_wv_pops[model == 'equilibrium' & ascertainment == 'gen0' & gen %in% c(10,100,1000)], 
#       aes(x = level, y = correction_propvar)) +  
ggplot(h_mean_wv_pops[model == 'equilibrium' & ascertainment == 'gen0' & gen %in% c(10,100,1000)], 
       aes(x = level, y = mean_correction_propvar)) +  
  geom_line(data=eq_thry[n.sample==1 & gen %in% c(10,100,1000)], aes(x = level, y = propvar), color = 'red', size = 1.2) + 
  
  geom_point() + 
  facet_wrap(~gen) + 
  geom_errorbar(aes(ymin=mean_correction_propvar - 1.96*se_correction_propvar, ymax=mean_correction_propvar+1.96*se_correction_propvar)) +
  #geom_line(group=1) + 
  theme_classic() + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') + 
  #geom_point(data=eq_thry[n.sample==1 & gen %in% c(10,100,1000)], aes(x = level, y = propvar), color = 'red') +
  #geom_line(data=eq_thry[n.sample==1 & gen %in% c(10,100,1000)], aes(x = level, y = propvar), color = 'red', size = 1.2) + 
  theme(aspect.ratio=1, 
        axis.text.x = element_text(angle = 90))

# compare magnitudes of covariance and single haplotype variance
ggplot(all_snpWC[model == 'equilibrium' & pop == 'p0' & level != 's14'], aes(x = level, y = mean_cov)) + 
  facet_wrap(~gen) + geom_point() + 
  geom_point(data = allSNPWV[model == 'equilibrium' & pop == 'p0' & level != 's14' & signal == 'single_hap'], aes(x = level, y = variance, color = 'red'))

# bottleneck

allAncestryWV[, .(variance=mean(variance)), by = .(ancestry_measure, signal, model, level, gen)]
ggplot(h_mean_wv_pops[model == 'bottleneck' & ascertainment == 'gen0' & gen != 0], 
       aes(x = level, y = mean_correction_propvar)) +  
  geom_point() + 
  facet_wrap(~gen) + 
  geom_line(group = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') + 
  geom_point(data=bn_thry[n.sample==1], aes(x = level, y = propvar), color = 'red') +
  geom_line(data=bn_thry[n.sample==1], aes(x = level, y = propvar), color = 'red') + 
  theme(aspect.ratio=1, 
        axis.text.x = element_text(angle = 90))

  

# ===== SNP Stat sample mean wavelet variance =====

# raw variance before correction
# get mean over replicates
meanSNPWV <- allSNPWV[level != 's14', .(variance = mean(variance)), by = .(gen, pop, level, ascertainment, signal, model)]

ggplot(meanSNPWV[gen != 0  & model == 'equilibrium' & signal == 'smpl_mean'], aes(x = level, y = variance, color = pop)) + 
  geom_point() + 
  facet_wrap(~gen) + 
  theme(aspect.ratio = 1) +
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') 

ggplot(meanSNPWV[gen != 0  & model == 'equilibrium' & signal == 'pop_mean'], aes(x = level, y = variance, color = ascertainment)) + 
  geom_point() + 
  facet_wrap(~gen, scales = 'free_y') + 
  theme(aspect.ratio = 1) +
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') 

# check whether ascertainment has differential effects across scales
meanSNPWV[, propvar := variance/sum(variance), by = .(gen, pop, ascertainment, signal, model)]

ggplot(meanSNPWV[gen != 0  & model == 'equilibrium' & signal == 'pop_mean'], aes(x = level, y = propvar, color = ascertainment)) + 
  geom_point() + 
  facet_wrap(~gen, scales = 'free_y') + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') +
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 90)) 


# ===== True Ancestry Wavelet Covariance =====
# For true ancestry, the only covariance among hybrid haplotypes is that generated by drift and selection
# post-admixture
ggplot(all_trueWC[gen != 0], aes(x = level, y = mean_cov, color = model)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, model))) + 
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=90)) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Covariance true ancestry') 

# note there is no fine scale covariance of true ancestry 
  
# ===== SNP Stat Wavelet Covariance =====
# for alpha 0.5, the initial wavelet covariance of hybrid haplotypes due to sampling is average of mean parental within-pop cov and
# parental between-pop cov. additional covariance is generated by drift and selection
ggplot(eq_snpWC[gen != 0], aes(x = level, y = mean_cov, color = pop)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, pop)))

ggplot(bn_snpWC[gen != 0], aes(x = level, y = mean_cov, color = pop)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, pop)))     
# can clearly see bottleneck generates additional covariance

# ----- subtract covariance due to sampling -----
cov_initial <- all_snpWC[gen == 0]
cov_initial <- dcast(cov_initial, rep + level + model ~ pop, value.var = 'mean_cov')
cov_initial[, p2 := NULL]
cov_initial <- cov_initial[, .(cov_initial = 0.25*p0 + 0.25*p1 + 0.5*p0_p1), by = .(rep, level, model)]

cov_contemp <- all_snpWC[pop == 'p2'][, .(rep, gen, level, model, mean_cov)]
setnames(cov_contemp, 'mean_cov', 'cov_contemp')
cov_correction_data <- merge(cov_initial, cov_contemp, by = c("rep", "level", "model"))
cov_correction_data[, post_admix_cov := cov_contemp - cov_initial]

# covariance generated after admixture
ggplot(cov_correction_data[gen !=0 & model == 'equilibrium'], aes(x = level, y = post_admix_cov, color = model)) + 
  facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, model))) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'delta_cov') + 
  theme(axis.text.x = element_text(angle=90)) 

# we've now subtracted covariance due to initial sampling, but there is remaining noise
# not due to admixture covariance. This comes from parental single haplotype variance
# resulting from parental haplotypes coalescing in the hybrid population. obtain probability of coalescing as weight
cov_correction_data[model == 'equilibrium', prob_coal := 1-exp(-gen/20000)]

# merge covariances and variances
h_wv_pops[, prnt_wv := (p0+p1)/2, by = .(rep, gen, level, ascertainment, model)]
cov_correction_data2 <- merge(cov_correction_data, h_wv_pops[ascertainment == 'gen0'], by = c("rep", "level", "model", "gen"))

cov_correction_data2[, cov_correction := cov_contemp - (1-prob_coal)*cov_initial - prob_coal*prnt_wv ]

# covariance correction
ggplot(cov_correction_data2[gen !=0 & model == 'equilibrium'], aes(x = level, y = cov_correction, color = model)) + 
#ggplot(cov_correction_data2[gen !=0 & model == 'equilibrium'], aes(x = level, y = cov_contemp-cov_initial, color = model)) + 
  
  facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, model))) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'delta_cov') + 
  theme(axis.text.x = element_text(angle=90))  #+
  #geom_line(aes(y = post_admix_cov, group = rep), color = 'black')


# check against theory for sample mean
smpl_mean_correction_data <- merge(allSNPWV[signal == 'smpl_mean' &  ascertainment == 'gen0' & pop == 'p2', .(smpl_mean_wv = variance, rep, gen, level, model)], 
      cov_correction_data2, by = c("rep", "gen", "level", "model"))

smpl_mean_correction_data[, smpl_mean_correction := smpl_mean_wv - (1/20)*prnt_wv - (19/20)*((1-prob_coal)*cov_initial + prob_coal*prnt_wv)]
smpl_mean_correction_data[, smpl_mean_correction2 := smpl_mean_wv - (1/20)*prnt_wv - (19/20)*(cov_initial)]

smpl_mean_correction_data[, smpl_mean_correction_propvar := smpl_mean_correction/sum(smpl_mean_correction), by = .(rep, gen,  model)]
#smpl_mean_correction_data[, smpl_mean_correction_propvar2 := smpl_mean_correction2/sum(smpl_mean_correction2), by = .(rep, gen,  model)] # for examining effect of only subtracting off sampling covariance


mean_smpl_mean_correction_data <- smpl_mean_correction_data[, .(mean_correction_propvar = mean(smpl_mean_correction_propvar), 
                                                                se_correction_propvar = sd(smpl_mean_correction_propvar)/sqrt(10)) , by = .(level, gen, model)]


ggplot(mean_smpl_mean_correction_data[ model == 'equilibrium' & gen %in% c(10,100,1000)], aes(x = level, y = mean_correction_propvar)) + 
  geom_line(data = eq_thry[n.sample == 20 & gen != 1], aes(x = level, y = propvar), color = 'red', size=1.2) + 
  
  facet_wrap(~gen, scales = 'free_y') + geom_point() + #geom_line(aes(group = rep)) +
  geom_errorbar(aes(ymin=mean_correction_propvar - 1.96*se_correction_propvar, ymax = mean_correction_propvar + 1.96*se_correction_propvar)) +
  #geom_line(data = eq_thry[n.sample == 20 & gen != 1], aes(x = level, y = propvar), color = 'red', size=1.2) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=90), aspect.ratio = 1)
 
# problem: fine scale covariance being generated. 
# is there fine scale covariance generated in the parental pops over time?




prnt_cov_initial <- all_snpWC[pop %in% c("p0", "p1", "p0_p1") & gen == 0 & model == 'equilibrium']
setnames(prnt_cov_initial, 'mean_cov', 'cov_initial')
prnt_cov_initial[, gen := NULL]
 
prnt_cov_contemp <- all_snpWC[pop %in% c("p0", "p1", "p0_p1") & gen != 0 & model == 'equilibrium'] 
setnames(prnt_cov_contemp, 'mean_cov', 'cov_contemp')

prnt_delta_cov <-  merge(prnt_cov_initial, prnt_cov_contemp, by = c("rep", "level", "pop", "model", "signal"))
prnt_delta_cov[, delta_cov := cov_contemp - cov_initial]
prnt_delta_cov <- prnt_delta_cov[, .(prnt_delta_cov = mean(delta_cov)), by = .(level, pop, gen, rep, model)]

# additional covariance is generated between pairs of haplotypes
ggplot(prnt_delta_cov, aes(x = level, y = prnt_delta_cov)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, pop), color = pop)) +
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'delta_cov') + 
  theme(axis.text.x = element_text(angle=90))

prnt_within_delta_cov <- prnt_delta_cov[pop %in% c("p0", "p1"), .(prnt_within_delta_cov = mean(prnt_delta_cov)), by = .(level, gen, rep, model)]
prnt_between_delta_cov <- prnt_delta_cov[pop == 'p0_p1', .(prnt_between_delta_cov = mean(prnt_delta_cov)), by = .(level, gen, rep, model)]

prnt_delta_cov_v2 <- merge(prnt_within_delta_cov, prnt_between_delta_cov)

cov_correction_data_v2 <- merge(cov_correction_data[model == 'equilibrium' & gen !=0], 
      prnt_delta_cov_v2, by = c("rep", "level", "gen", "model"))

cov_correction_data_v2[, correction := post_admix_cov - prnt_within_delta_cov - prnt_between_delta_cov, by = .(rep, level, gen, model)]

cov_correction_data_v2[, .(mean_correction = mean(correction)), 
                       by = .(level, gen, model)] %>%
  ggplot(., aes(x= level, y = mean_correction)) + 
  facet_wrap(~gen) + 
  geom_line(aes(group = 1))

dcast(all_snpWC[model == 'equilibrium'], rep + level + pop ~ gen, value.var = 'mean_cov')
setnames(all_snpWC_gen_wide, c('0','1','10','100','1000'), paste0('gen', c(0,1,10,100,1000)))

all_snpWC_gen_wide[, ]

  

# ----- similar approach, but subtracting off parental covs in contemporary gens -----
all_snpWC_wide <- dcast(all_snpWC, rep + level + gen + model ~ pop, value.var = 'mean_cov')
all_snpWC_wide[, correction := p2 - 0.25*p0 - 0.25*p1 - 0.5*p0_p1, by = .(rep, level, gen, model)]
ggplot(all_snpWC_wide[model == 'equilibrium'], aes(x=level, y = correction, color = model)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, model)))


all_snpWC[model == 'equilibrium', .(mean_cov = mean(mean_cov)), by = .(gen, level, pop, model)] %>%
  ggplot(., aes(x = level, y = mean_cov)) + 
  geom_line(aes(group = interaction(pop, gen, model), linetype=as.factor(gen), color = pop))

all_snpWC_gen_wide <- dcast(all_snpWC[model == 'equilibrium'], rep + level + pop ~ gen, value.var = 'mean_cov')
setnames(all_snpWC_gen_wide, c('0','1','10','100','1000'), paste0('gen', c(0,1,10,100,1000)))
ggplot(all_snpWC_gen_wide, aes(x = level, y = gen1000-gen0, color = pop)) + geom_line(aes(group = interaction(pop,rep)))

p0p1_cov <- dcast(all_snpWC[pop=='p0_p1' & model == 'equilibrium'], rep + level ~ gen, value.var = 'mean_cov')
setnames(p0p1_cov, c('1000', '0'), c('g1000','g0'))
ggplot(p0p1_cov, aes(x = level, y = g1000-g0)) + geom_line(aes(group = rep))

# ---- try using this in the correction for sample mean

smpl_mean_correction_data <- merge(merge(allSNPWV[pop == 'p2' & ascertainment == 'contemporary' & signal == 'smpl_mean' & model == 'equilibrium',
               .(rep, gen, level, smpl_mean_wv = variance)],
      h_wv_pops[ascertainment == 'contemporary' & model == 'equilibrium', 
                .(rep, gen, level, h_prnt_wv = (p0+p1)/2)]), 
      cov_initial)

smpl_mean_correction_data[, correction := smpl_mean_wv - (1/n.sample)*h_prnt_wv - ((n.sample-1)/n.sample)*cov_initial]

ggplot(smpl_mean_correction_data, aes(x = level, y= correction)) + facet_wrap(~gen) + geom_line(aes(group = rep))

cov_contemp <- allWC[gen != 0 & pop == 'p2'][, pop := NULL]
setnames(cov_contemp, 'mean_cov', 'cov_contemp')

cov_correction_data <- merge(cov_initial, cov_contemp, by = c("rep", "level", "model"))

cov_correction_data[, delta_cov := cov_contemp - cov_initial, by = .(rep, level, model, gen)]

ggplot(cov_correction_data[ model == 'bottleneck'], aes(x = level, y = delta_cov)) + facet_wrap(~gen) + 
  geom_line(aes(group = interaction(rep, model))) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'delta_cov') + 
  theme(axis.text.x = element_text(angle=90))
  

meanWC <- allWC[, .(mean_cov = mean(mean_cov)), by = .(gen, level, pop, model)]

meanWC_wide <- dcast(meanWC, gen + level + model ~ pop, value.var = 'mean_cov')
meanWC_wide[, prnt_cov := sum(p0, p1)/2, by = .(gen, level, model)]
meanWC_wide[, p2sum := 0.5*prnt_cov + 0.5*p0_p1]
ggplot(meanWC_wide[model == 'equilibrium'], aes(x = level)) + facet_wrap(~gen) + 
  geom_line(aes(y = p2, group = 1), color = 'black') + 
  geom_line(aes(y = prnt_cov, group = 1, color = 'red')) + 
  geom_line(aes(y = p0_p1, group = 1), color = 'blue')+
  geom_line(aes(y = p2sum, group = 1), color = 'green')
  


# sample mean correction 
prnt_hap_wv <- meanSNPWV[signal == 'single_hap' & pop != 'p2' & ascertainment == 'contemporary', .(prnt_hap_wv = mean(variance)), by = .(gen, level,  model)]
#prnt_hap_wv[, gen := NULL]
hyb_hap_wv <- meanSNPWV[signal == 'single_hap' & pop == 'p2' & ascertainment == 'contemporary', .(hyb_hap_wv = variance), by= .(gen, level,  model)]
hyb_smpl_mean_wv <- meanSNPWV[signal == 'smpl_mean' & pop == 'p2' & ascertainment == 'contemporary', .(hyb_smpl_mean_wv = variance), by = .(gen, level,  model)]

prnt_same_pop_cov <- meanWC[ pop %in% c("p0", "p1"), .(prnt_same_pop_cov = mean(mean_cov)), by = .(gen, level, model)]
#prnt_same_pop_cov[, gen := NULL]
prnt_diff_pop_cov <- meanWC[pop == 'p0_p1', .(prnt_diff_pop_cov = mean_cov), by = .(gen, level, model)]
#prnt_diff_pop_cov[, gen := NULL]

all_correction_data <- merge(merge(merge(prnt_same_pop_cov, prnt_diff_pop_cov), prnt_hap_wv), merge(hyb_hap_wv, hyb_smpl_mean_wv))
#all_correction_data <- merge(merge(merge(merge(prnt_hap_wv, hyb_hap_wv), hyb_smpl_mean_wv), prnt_same_pop_cov),prnt_diff_pop_cov)

all_correction_data[, correction := hyb_smpl_mean_wv - (1/20)*prnt_hap_wv - (19/20)*(0.5*prnt_same_pop_cov + 0.5*prnt_diff_pop_cov)]
all_correction_data[, correction_propvar := correction/sum(abs(correction)), by = .(gen, model)]

ggplot(all_correction_data[gen != 0 & model == 'equilibrium'], aes(x = level, y = correction_propvar)) + 
  facet_wrap(~gen, scales = 'free_y') + geom_point()+ 
  geom_line(data=eq_thry[n.sample==20], aes(x = level, y = propvar), color = 'red')


ggplot(all_correction_data[model == 'bottleneck'], aes(x = level)) + facet_wrap(~gen) + 
  geom_line(aes(y = (1/20)*prnt_hap_wv, group = 1), color = 'red') + 
  geom_line(aes(y = (19/20)*0.5*prnt_same_pop_cov, group = 1), color = 'darkgreen') + 
  geom_line(aes(y = (19/20)*0.5*prnt_diff_pop_cov, group = 1), color = 'purple')


smpl_wv_pops <- dcast(allSNPWV[signal=='smpl_mean'], rep+gen+level+ascertainment+ model ~ pop, value.var='variance')
smpl_wv_pops[, ]
setnames(smpl_wv_pops, c("p0", "p1", "p2"), c("p0_smpl_mean_wv", "p1_smpl_mean_wv", "p2_smpl_mean_wv"))

setnames(meanWC_wide, "p0_p1", "p0_p1_cov")
correction_data <- merge(smpl_wv_pops, meanWC_wide[, .(gen, level, model, p0_p1_cov, prnt_cov)])

correction_data[, correction := p2_smpl_mean_wv - (p0_p1_cov + prnt_cov)/2, by = .(gen, level, model, rep, ascertainment)]
ggplot(correction_data, aes(x = level, y = correction)) + facet_wrap(~gen) + geom_point()

correction_data[, mean(.SD), .SDcols = .(p2_smpl_mean_wv, p0_p1_cov, prnt_cov, correction), 
                by = .(gen, level, model, ascertainment)]

# smpl_wv_pops[, correction := p2 - sum(p0, p1)/2, by = .(rep, gen, level, ascertainment, model)]
# smpl_wv_pops_means <- smpl_wv_pops[, .(mean_correction = mean(correction)), by = .(gen, level, ascertainment, model)]
# smpl_wv_pops_means[, correction_propvar := mean_correction/sum(mean_correction), by = .(gen, ascertainment, model)]
# 
# ggplot(smpl_wv_pops[ascertainment == 'gen0' & model == 'equilibrium'], aes(x = level)) + geom_point(aes(y = p2-p0)) + facet_wrap(~gen)
# 
# # potential correction - variance -  (1/n) * single hap parental variance - ((n-1)/n) * prnt cov
# # examine covariances - looks roughly same magnitude
# ggplot(allWC, aes(x = level, y = mean_cov)) + geom_point() + facet_wrap(~gen)
# 
# # get mean covariance among parental chromosomes
# prnt_cov <- allWC[pop != 'p2', .(prnt_cov = mean(mean_cov)), by = .(gen, level, model)]
# ggplot(prnt_cov, aes(x = level, y = prnt_cov)) + geom_point() + facet_wrap(~gen)
# 
# smpl_mean_wv <- allSNPWV[pop=='p2' & signal == 'smpl_mean' & level != 's14' & gen != 0, .(smpl_mean_variance = mean(variance)), by = .(gen, level, model)]
# ggplot(smpl_mean_wv, aes(x = level, y = smpl_mean_variance, color = model)) + geom_point() + facet_wrap(~gen)
# prnt_single_hap_wv <- allSNPWV[pop!='p2' & signal == 'single_hap', .(single_hap_variance = mean(variance)), by = .(gen, level, model)]
# 
# smpl_mean_data <- merge(merge(smpl_mean_wv, prnt_cov), prnt_single_hap_wv)
# smpl_mean_data[, correction := smpl_mean_variance - (19/20)*prnt_cov - (1/20)*single_hap_variance ]
# 
# ggplot(smpl_mean_data, aes(x = level, y= correction, color = model)) + facet_wrap(~gen) + geom_point()
# 
# # equilibrium
# ggplot(smpl_wv_pops_means[model == 'equilibrium' & ascertainment == 'gen0' & gen != 0], 
#        aes(x = level, y = correction_propvar)) +  
#   geom_point() + 
#   facet_wrap(~gen) + 
#   geom_line(group = 1) + 
#   theme_classic() + 
#   scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') + 
#   geom_point(data=eq_thry[n.sample==1], aes(x = level, y = propvar), color = 'red') +
#   geom_line(data=eq_thry[n.sample==1], aes(x = level, y = propvar), color = 'red') + 
#   theme(aspect.ratio=1, 
#         axis.text.x = element_text(angle = 90))




# ===== Sample Mean Wavelet Variance =====

mean_wf_smpl_mean <- allWV[, mean(variance.h_frq), by = .(gen, level, ascertainment, model, pop)]

# ===== Population Mean Wavelet Variance











