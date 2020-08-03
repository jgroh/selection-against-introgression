library(data.table)



for(i in 1:replicates){
  
  filename <- "out1.txt"
  slimOutLines <- readLines(filename)
  skip1 <- grep("output", slimOutLines)
  
  dat <- fread("out1.txt", skip = skip1, col.names = c("proportion_introgressed"))
  
  y <- dat$proportion_introgressed
  x <- 1:10
  
  if (i == 1){
    plot(y ~ x , ylim = c(0,.1), type = "n",
         xlab = "generation", ylab = "Avg. genomic proportion introgressed")
  }
  
  lines(y ~ x, add = TRUE, col = "gray")
}
  

