


params <- NULL
for(g in c(5,10,50,100,500,1000)){
  for(i in 1:10){
   params <- c(params,paste(g,i,1:(2^((11-i)-1))))
  }
}
tail(params)
pars <- as.matrix(params)
write.table(pars,
            file = "~/workspace/selection-against-introgression/simulations/code/wv_batch_params.txt",
            quote = F,row.names = F,col.names=F)

