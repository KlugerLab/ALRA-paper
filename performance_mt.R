source('../ALRA/alra.R')
library(SAVER)
library(Rmagic)
library(scImpute)
if ( !file.exists( 'data/1M_neurons.RData')) {
	A_norm <- read.csv('data/1M_neurons.csv',header=F)
	save(A_norm,file='data/1M_neurons.RData')
}else{
	load('data/1M_neurons.RData')
}


sizes <- c(1E3, 5E3, 10E3, 20E3, 30E3, 40E3, 50E3)
alra_times <- rep(NA,length(sizes))
saver_times <- rep(NA,length(sizes))
saver_fast_times <- rep(NA,length(sizes))
magic_fast_times <- rep(NA,length(sizes))
dca_fast_times <- rep(NA,length(sizes))
scimpute_fast_times <- rep(NA,length(sizes))


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
	alra_times[i] <- endtime['elapsed']
	save(alra_times,file="data/alra_fast_times2.RData")
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
	save(sizes,magic_fast_times,file="data/magic_fast_times2.RData")
}

library(R.utils)

for (i in 1:3) {
	print(sprintf("Size: %d", sizes[i]))
	print(Sys.time())
	set.seed(5)
	samp.cells <- sample(1:ncol(x), sizes[i])
	x.sub <- x[,samp.cells]
	starttime <- proc.time()
	out <- saver(x.sub, do.fast = TRUE)
	endtime <- proc.time()-starttime
	saver_fast_times[i] <- endtime['elapsed']
	save(sizes,saver_fast_times,file="data/saver_fast_times2.RData")
}


for (i in 1:length(sizes)) {
	print(sprintf("Size: %d", sizes[i]))
	set.seed(5)
	samp.cells <- sample(1:ncol(x), sizes[i])
	x.sub <- x[,samp.cells]
        write.csv(x.sub, 'data/temp_csv.csv') 
	starttime <- proc.time()
        system('export CUDA_VISIBLE_DEVICES=;dca data/temp_csv.csv data/dca_output --threads=1')
	endtime <- proc.time()-starttime
	dca_fast_times[i] <- endtime['elapsed']
	save(sizes,dca_fast_times,file="data/dca_fast_times2.RData")
}


for (i in 1:3) {
        print(sprintf("Size: %d", sizes[i]))
        set.seed(5)
        samp.cells <- sample(1:ncol(x), sizes[i])
        x.sub <- x[,samp.cells]
        fni =  "data/scimpute_bench/scimpute_count_bench.csv"
        write.table( x.sub, fni, sep = ",", quote = F, row.names = T, col.names = T)
	starttime <- proc.time()
        scimpute_result <- scimpute(
          count_path =fni, 
          out_dir = "data/scimpute_bench/", 
          Kcluster = 36, #Hrvatin has 36 different clusters 
          drop_thre = 0.3,
        ncores = 1)
	endtime <- proc.time()-starttime
	scimpute_fast_times[i] <- endtime['elapsed']
	save(sizes,scimpute_fast_times,file="data/scimpute_fast_times2.RData")
}


#### Load and plot
load("data/saver_fast_times2.RData")
load("data/magic_fast_times.RData")
load("data/alra_fast_times2.RData")
load("data/dca_fast_times2.RData")
load("data/scimpute_fast_times2.RData")

library(ggplot2)
library(cowplot
alra_df         <- data.frame(sizes=sizes, times=alra_times/60/60)
magic_df           <- data.frame(sizes=sizes, times=magic_fast_times/60/60)
saver_df              <- data.frame(sizes=sizes[1:2], times=saver_fast_times[1:2]/60/60)
scimpute_df            <- data.frame(sizes=sizes[1:3], times=scimpute_fast_times[1:3]/60/60)
dca_df          <- data.frame(sizes=sizes, times=dca_fast_times/60/60)
lwidth = 0.7; psize=5;

library(reshape2)
df <- melt(list("ALRA"= alra_df,"MAGIC"= magic_df,"SAVER"=saver_df,"scImpute" = scimpute_df,"DCA" = dca_df), id.vars="times")

g<- ggplot(df, aes(value,times, color=L1, shape=L1)) + geom_point(size=psize) +geom_line()+theme(legend.position="bottom",legend.title=element_blank(), legend.justification = "center", panel.grid.major=element_line(colour = "grey90"), legend.spacing.x=unit(0.4, "cm"))+
scale_x_continuous(labels=scales::comma,breaks=seq(0,50000,10000) )+
labs(x="N",y="Runtime (hours)")
g
ggsave(g,filename="figs/runtime.pdf",width=7,height=4)

hello
