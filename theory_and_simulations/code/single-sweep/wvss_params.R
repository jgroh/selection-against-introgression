library(tidyr)
params <- NULL
for(g in c(1000)){
  for(i in 1:10){
   params <- c(params,paste(g,i,1:(2^((11-i)-1)),sep = "_"))
  }
}
pars <- params %>% as.data.table() %>%
  separate(1, into = c("gen","scale","shift"), sep = "_") %>%
  .[, lapply(.SD, as.numeric)]

# use sparse wavelets to approximate the integral for speed here
pars2 <- rbindlist(list(pars[scale %in% 8:10 & shift %% 2 != 0,],
                        pars[scale %in% 6:7 & shift %% 4 == 0,],
                        pars[scale %in% 4:5 & shift %% 8 == 0,],
                        pars[scale %in% 2:3 & shift %% 16 == 0,],
                        pars[scale ==1 & shift %% 64 == 0,]
                        ))

write.table(pars2,
            file = "~/workspace/selection-against-introgression/theory_and_simulations/code/single-sweep/wvss_params.txt",
            quote = F,row.names = F,col.names=F)

