source('../ALRA/alra.R')
library(SAVER)
library(Rmagic)
#x <- readRDS("SAVER-data/hrvatin.rds")
if (!file.exists('data/hrvatin_full.rds')) {
	x <- read.csv("SAVER-data/GSE102827_merged_all_raw.csv", header = TRUE, row.names = 1, check.names = FALSE)
	x <- as.matrix(x)
	# Filter out genes with expression less than 0.00003
	x1 <- x[rowMeans(x) >= 0.00003, ]
	# look at non-zero expression
	# Filter out genes with non-zero expression in less than 4 cells
	x2 <- x1[rowSums(x1 != 0) >= 4, ]
	x <- x2
	saveRDS(x,'data/hrvatin_full.rds')
}else{
	print("Reading existing data")
	x <- readRDS('data/hrvatin_full.rds')
}

sizes <- c(1E3, 5E3, 10E3, 20E3, 30E3, 40E3, 50E3)
slra_times <- rep(NA,length(sizes))
saver_times <- rep(NA,length(sizes))
saver_fast_times <- rep(NA,length(sizes))
magic_fast_times <- rep(NA,length(sizes))
for (i in 1:length(sizes)) {
	print(sprintf("Size: %d", sizes[i]))
	set.seed(5)
	samp.cells <- sample(1:ncol(x), sizes[i])
	x.sub <- x[,samp.cells]
	x.sub.t <- t(x.sub)

	starttime <- proc.time()
		x_norm <- normalize_data(x.sub.t)
		k_choice <- choose_k(x_norm)
		print(sprintf("k=%d was chosen", k_choice$k))
		result.completed <- alra(x_norm,k_choice$k)
	endtime <- proc.time()-starttime
	slra_times[i] <- endtime['elapsed']
	save(slra_times,file="data/slra_fast_times.RData")
}


for (i in 1:length(sizes)) {
	print(sprintf("Size: %d", sizes[i]))
	set.seed(5)
	samp.cells <- sample(1:ncol(x), sizes[i])
	x.sub <- x[,samp.cells]
	x.sub.t <- t(x.sub)
	starttime <- proc.time()
	out <- magic(x.sub.t)
	endtime <- proc.time()-starttime
	magic_fast_times[i] <- endtime['elapsed']
	save(sizes,magic_fast_times,file="data/magic_fast_times.RData")
}

library(R.utils)

for (i in 1:2) {
	print(sprintf("Size: %d", sizes[i]))
	print(Sys.time())
	set.seed(5)
	samp.cells <- sample(1:ncol(x), sizes[i])
	x.sub <- x[,samp.cells]
	starttime <- proc.time()
	out <- saver(x.sub, do.fast = TRUE)
	endtime <- proc.time()-starttime
	saver_fast_times[i] <- endtime['elapsed']
	save(sizes,saver_fast_times,file="data/saver_fast_times.RData")
}



#### Load and plot
load("data/saver_fast_times.RData")
load("data/magic_fast_times.RData")
load("data/slra_fast_times.RData")

library(ggplot2)
library(cowplot)
slra_df <- data.frame(sizes=sizes, times=slra_times/60/60)
magic_df <- data.frame(sizes=sizes, times=magic_fast_times/60/60)
saver_df <- data.frame(sizes=sizes[1:2], times=saver_fast_times[1:2]/60/60)
lwidth = 0.7; psize=2;
g <- ggplot() + geom_line(data=slra_df, aes(x=sizes,y=times,color="SLRA"),size=lwidth) + 
	 geom_point(data=slra_df, aes(x=sizes,y=times,color="SLRA"),size=psize) + 
	geom_line(data=magic_df, aes(x=sizes,y=times,color="MAGIC"),size=lwidth) +
	geom_point(data=magic_df, aes(x=sizes,y=times,color="MAGIC"),size=psize) +
	geom_line(data=saver_df, aes(x=sizes,y=times,color="SAVER"),size=lwidth) +
	geom_point(data=saver_df, aes(x=sizes,y=times,color="SAVER"),size=psize) +
	theme(legend.title=element_blank(), legend.position=c(1,0))+
	scale_x_continuous(labels=scales::comma,breaks=sizes)+
	labs(x="N",y="Runtime (hours)")

ggsave(g,filename="figs/runtime.pdf",width=6,height=4)


