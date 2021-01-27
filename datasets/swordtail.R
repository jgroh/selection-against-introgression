library(data.table)
a1 <- fread("datasets/Population1_TOTO/ancestry-probs-par1_TOTO_allgroups.tsv", nrows = 1)
a2 <- fread("Population1_TOTO/ancestry-probs-par2_TOTO_allgroups.tsv", nrows = 1)
names(a1) == names(a2)
smpl <- sample(names(a1)[-1], size = 100000, replace = F)
a1_smpl <- a1[, ..smpl]
a2_smpl <- a2[, ..smpl]
x <- as.numeric(a1_smpl[1,])
y <- as.numeric(a2_smpl[1,])
par(mar = c(5,4,4,1))
plot(y = y, x = x)
hist(x+y, xlab = "Sum of ancestry posterior probabilities",
     freq = F, main = "", col = "red")

# look at tract lengths on chr 1
cols <- names(a1)[grep("group1:", names(a1))]
a1.c1 <- as.data.frame(a1[, ..cols])
a1.c1.pos <-  as.numeric(sub(".*:", "", x= names(a1.c1)))
a1.c1  <- as.numeric(a1.c1[1,])
plot(y = a1.c1,x = a1.c1.pos, xlab = "Position on chr1", ylab = "Ancestry posterior probability")

mean(a1.c1, na.rm = T)
