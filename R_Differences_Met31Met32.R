require("qvalue")
require("gplots")


#Allegra <- read.table("Allgera_M31M32Different.txt",header=FALSE)

#row.names(Allegra) = Allegra$V2

WD = "/Users/smcisaac/Downloads"
setwd(WD)

y <-read.table("DataMetLimitedMerged_2011_1122_pruned_knn_clustered_SVDsubtract_RECLUSTER.txt",header =TRUE, row.names = 1, sep ="\t")

#only get Met31 and Met32 data
y <- y[,9:24]

#remove genes that are 0 at each timepoint
Check <- apply(y,1,function(x){sum(x==0)==16}) #indices of zero rows

y <- y[!Check,]

Input <- y

#we want there to be variation between the initial time points (0-5 minutes) and the later timepoints when we expect real action from the TF. Require that there is at least a 1.5-fold change between timepoints at(15-90) and either 2.5 and 5. 

Dim <- dim(y)[1]
Cut <- log(1.5,base=2)

y2 <-y

for(i in seq(Dim,1)){
	
	vec1 = y[i,c(2,4,5,6,7,8)]
	vec2 = y[i,3:8]
	
	vec3 = y[i,c(10,12,13,14,15,16)]
	vec4 = y[i,3:8]
	
	vec1 <- vec1 - vec1[[1]]
	vec2 <- vec2 - vec2[[1]]
	vec3 <- vec3 - vec3[[1]]
	vec4 <- vec4 - vec4[[1]]
	
	val1 <- max(abs(vec1 - vec1[[1]]))
	val2 <- max(abs(vec2 - vec2[[1]]))
	val3 <- max(abs(vec3 - vec3[[1]]))
	val4 <- max(abs(vec4 - vec4[[1]]))
	
	if(val1 < Cut && val2 < Cut && val3 < Cut && val4 < Cut){
		y <- y[-i,]
		}
	
	
	}





#do regression on just Met31 and Met32 to check for genes that have significant time dependence.
#Significance is determined by computing the p-value of the F-statistic. If the p-value is less than 0.05
#in at least 1 of the time courses keep it. 
 
Met31 <- y[,1:8]
Met32 <- y[,9:16]

Dim = dim(Met31)[1]
t1 <- c(0,2.5,5,15,30,45,60,90)

model.Met31 <- apply(Met31,1,function(x){lm(x ~ -1 + t1 + I(t1^2))})
model.Met32 <- apply(Met32,1,function(x){lm(x ~ -1 + t1 + I(t1^2))})

fpval.Met31 <- rep(0,dim(Met31)[1]) #store pvalues of the f-statistic
fpval.Met32 <- rep(0,dim(Met31)[1]) #store pvalues of the f-statistic
fpval.keep <- rep(0,dim(Met31)[1])

for(i in seq(1,Dim)){
	
	
	temp.Met31 <- summary(model.Met31[[i]])
	temp.Met32 <- summary(model.Met32[[i]])
	
	fpval.Met31[i] <- pf(temp.Met31$fstatistic[1],temp.Met31$fstatistic[2],temp.Met31$fstatistic[3],lower.tail=FALSE)
	fpval.Met32[i] <- pf(temp.Met32$fstatistic[1],temp.Met32$fstatistic[2],temp.Met32$fstatistic[3],lower.tail=FALSE)
	
	if(is.nan(fpval.Met31[i])){
		fpval.Met31[i] = .1
	}
	
	if(is.nan(fpval.Met32[i])){
		fpval.Met32[i] = .1
	}

	if(fpval.Met31[i]<0.05 | fpval.Met32[i]<0.05){
		fpval.keep[i] = 1
		}
	
	
	}

y <- y[!!fpval.keep,]
Met31 <- Met31[!!fpval.keep,]
Met32 <- Met32[!!fpval.keep,]


#time points
t <- c(0,2.5,5,15,30,45,60,90,0,2.5,5,15,30,45,60,90)

#genes that are different between Met31 and Met32
Ds <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
#Model A
A <- apply(y,1,function(x){lm(x ~ -1 + t + I(t^2) + (Ds:t) + (Ds:I(t^2)))})

#Model B
B <- apply(y,1,function(x){lm(x ~ -1 + t + I(t^2))})

#Log-likelihoods
LA <- sapply(A,function(x){logLik(x)})

LB <- sapply(B,function(x){logLik(x)})

#log-likelihood test statistic
D <- -2*LB + 2*LA

Pval <- pchisq(D,df = 2, lower.tail = F) #is it 3 or 2? 

Q <- qvalue(Pval,lambda = 0)

#get the qvalues for each row 
Q2 <- Q$qvalues

DistStore <- rep(0,dim(Met31)[1])

#find maximum distance between genes by timepoint
for(i in seq(1,dim(Met31)[1])){
	
	temp <- dist(rbind(as.vector(data.matrix(Met31[i,])),as.vector(data.matrix(Met32[i,]))),method="maximum")[1]
	DistStore[i] = temp
	
	
	}

names(DistStore) = row.names(y)

#Store Q-values + distances, and find if Q < 0.05 + dist >0.585 (1.5-fold)
QD <- t(rbind(DistStore,Q2))
QD2 <- cbind(QD,rep(0,dim(Met31)[1]))

for(i in seq(1,dim(Met31)[1])){
	
	if(QD2[i,1]>log(1.5,base=2) && QD2[i,2]<0.05){
		QD2[i,3] = 1	
		}
	
	
	}


Output <- y[!!as.vector(QD2[,3]),]


#compare to genes that are different from Allegra's experiments. 
Row <- row.names(Output)
Z <- sapply(Row,function(x){strsplit(x," ")})
Z2 <- as.vector(sapply(Z,function(x){x[[3]]}))
intersect(Z2,row.names(Allegra))

notchosen <- Input[!(row.names(Input) %in% row.names(Output)),]
View(notchosen)

pairs.breaks <- c(seq(-2, 0, length.out=50),seq(0, 2, length.out=50))
pairs.breaks[51] = 0.01
mycol <- colorpanel(n=99,low="green",mid="black",high="red")
ClusterInfo <- heatmap.2(data.matrix(Output),breaks=pairs.breaks,Colv = FALSE,trace="none",scale="none",col=mycol,distfun=function(x) as.dist((1-cor(t(x)))/2),dendrogram = "none",density.info = "none",labCol="",cexRow=0.15)

ClusterInfo <- heatmap.2(data.matrix(notchosen),breaks=pairs.breaks,Colv = FALSE,trace="none",scale="none",col=mycol,distfun=function(x) as.dist((1-cor(t(x)))/2),dendrogram = "none",density.info = "none",labCol="",cexRow=0.15)





write.table(Output, file=paste(WD,'LinearModeling_Met31Met32_Big.txt', sep='/'), sep='\t', col.names=NA, row.names = TRUE, quote=FALSE)






########### two-sided t-test to compare Met31 and Met32 
Dim = dim(Met31)[1]
Pval.ttest = rep(0,Dim)

for(i in seq(1,Dim)){
	
	temp <- t.test(as.vector(data.matrix(Met31[i,])),as.vector(data.matrix(Met32[i,])),alternative="two.sided",paired = TRUE, conf.level = 0.95)
	Pval.ttest[i] <- temp$p.value
	
	}

Q.ttest <- qvalue(Pval.ttest)
Q2.ttest <- Q.ttest$qvalues



####what if we just took the means?

Met31 = Input[1:8,]
Met32 = Input[9:16,]


A <- as.numeric(apply(Input,1,function(x) {t.test(x[1:8],x[9:16])$p.value}))

mychoose <- A < 0.01

Input_thresholded = Input[mychoose,]

ClusterInfo <- heatmap.2(data.matrix(Input_thresholded),breaks=pairs.breaks,Colv = FALSE,trace="none",scale="none",col=mycol,distfun=function(x) as.dist((1-cor(t(x)))/2),dendrogram = "none",density.info = "none",labCol="",cexRow=0.15)



###do regression in tidy format!

g <- as.numeric(as.vector(Output[20,]))


df <- tibble(gene = rep(strsplit(row.names(Output[20,])," ")[[1]][3],16),
                 g = g,
                 t=t,
                 TF = c(rep("MET31",8),rep("MET32",8)),
                 Ds = Ds)

all_gene_fits = df %>%
  group_by(gene) %>%
  do(fit = lm(g ~ -1 + t + I(t^2) + (Ds:t) + (Ds:I(t^2)), data = .)) %>%
  tidy(fit)


df %>%
  group_by(gene) %>%
  do(fit = lm(g ~ -1 + t + I(t^2) + (Ds:t) + (Ds:I(t^2)), data = .)) %>%
  tidy(fit)

df %>%
  group_by(gene) %>%
  do(fit = lm(g ~ -1 + t + I(t^2) + (Ds:t) + (Ds:I(t^2)), data = .)) %>%
  glance(fit)
#return f-statistic?

# can we do nested model?

plot(lm(g ~ -1 + t + I(t^2) + (Ds:t) + (Ds:I(t^2))))



       