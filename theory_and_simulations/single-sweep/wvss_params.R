library(tidyr)
library(data.table)
params <- NULL
for(g in c(10,100,1000)){
  for(i in 1:10){
    for(s in c(0.01,0.1,0.5)){
   params <- c(params,paste(g,i,1:(2^((11-i)-1)),s,sep = "_"))
    }
  }
}
pars <- params %>% as.data.table() %>%
  separate(1, into = c("gen","scale","shift","s"), sep = "_") %>%
  .[, lapply(.SD, as.numeric)]

# use sparse wavelets to approximate the avg. over all wavelet placements at a scale
spars <- rbindlist(list(pars[scale %in% 8:10],
                        pars[scale %in% 6:7 & shift %% 2 == 0,],
                        pars[scale %in% 4:5 & shift %% 4 == 0,],
                        pars[scale %in% 2:3 & shift %% 8 == 0,],
                        pars[scale ==1 & shift %% 16 == 0,]
                        ))
spars.out <- as.matrix(spars[,paste(gen,scale,shift,s, sep = "_")])
write.table(spars.out,
            file = "wvss_params.txt",
            quote = F,row.names = F,col.names=F)

