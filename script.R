#######################################
#
# Program to analyze distance between
# Reddit subreddits using the cooccurrence
# of commentors across subreddits. 
# Also implements "subreddit algebra"
# by adding and subtracting subreddit
# vectors. 
# By @martintrevor_ for FiveThirtyEight
#
#######################################

library(reshape2)
library(lsa)
library(ggtern)
library(data.table)

##### Part 1: Load in the data

# This CSV file was created by running the SQL code in processData.sql in Google's BigQuery
rawsubredditvecs = fread("temp",header=TRUE)
rawsubredditvecs <- as.data.frame(rawsubredditvecs)

##### Part 2: Format and clean data for analysis

castsubredditvecs = dcast(rawsubredditvecs,t1_subreddit~t2_subreddit,FUN="identity",fill=0)
subredditvecst = as.matrix(castsubredditvecs[,-1])
rownames(subredditvecst) = castsubredditvecs[,1]
subredditvecs = t(subredditvecst)
subredditvecssums = apply(subredditvecs,1,sum)
subredditvecsnorm = sweep(subredditvecs,1,subredditvecssums,"/")
subredditvecssumscontext = apply(subredditvecs,2,sum)
contextprobs = subredditvecssumscontext/sum(subredditvecssumscontext)
subredditvecspmi = log(sweep(subredditvecsnorm,2,contextprobs,"/")) # PMI version
subredditvecsppmi = subredditvecspmi
subredditvecsppmi[subredditvecspmi<0] = 0 # PPMI version
scalar1 <- function(x) {x / sqrt(sum(x^2))} # Function to normalize vectors to unit length
subredditvecsppminorm = t(apply(subredditvecsppmi,1,scalar1))

##### Part 3: Analysis of subreddit similarities

## Looking at which subreddits are closest to each other (and combinations of subreddits)
cursubmat = subredditvecsppminorm
cursubmatt = t(cursubmat)
currownameslc = tolower(rownames(cursubmat))
# Function to calculate subreddit similarities and perform algebra
# Note that curops always has a leading "+"
findrelsubreddit <- function(cursubs,curops,numret=20) {
    cursubs = tolower(cursubs)
    curvec = 0
    for(i in 1:length(cursubs)) {
	    curvec = ifelse(curops[i]=="+",list(curvec + cursubmat[which(currownameslc==cursubs[i]),]),list(curvec - cursubmat[which(currownameslc==cursubs[i]),]))[[1]]
    }
    curclosesubs = cosine(x=curvec,y=cursubmatt)
    curclosesubso = order(curclosesubs,decreasing=TRUE)
    curclosesubsorder = curclosesubs[curclosesubso]
    curclosesubsorderc = curclosesubsorder[-which(tolower(names(curclosesubsorder))%in%cursubs)]
return(head(curclosesubsorderc,numret))
}

## Political examples

# /r/The_Donald
cursubs = c("the_donald")
curops = c("+")
findrelsubreddit(cursubs,curops,5)

# /r/The_Donald - /r/politics
cursubs = c("the_donald","politics")
curops = c("+","-")
findrelsubreddit(cursubs,curops,5)

# /r/hillaryclinton - /r/politics
cursubs = c("hillaryclinton","politics")
curops = c("+","-")
findrelsubreddit(cursubs,curops,5)

# /r/The_Donald - /r/SandersforPresident
cursubs = c("the_donald","sandersforpresident")
curops = c("+","-")
findrelsubreddit(cursubs,curops,5)

# /r/SandersforPresident - /r/The_Donald
cursubs = c("sandersforpresident","the_donald")
curops = c("+","-")
findrelsubreddit(cursubs,curops,5)

# /r/fatpeoplehate + /r/CoonTown + /r/politics
cursubs = c("fatpeoplehate","coontown","politics")
curops = c("+","+","+")
findrelsubreddit(cursubs,curops,5)

## Validation examples

# /r/nba + /r/minnesota
cursubs = c("nba","minnesota")
curops = c("+","+")
findrelsubreddit(cursubs,curops,5)

# /r/personalfinance - /r/Frugal
cursubs = c("personalfinance","frugal")
curops = c("+","-")
findrelsubreddit(cursubs,curops,5)

# /r/Fitness + /r/TwoXChromosomes
cursubs = c("fitness","twoxchromosomes")
curops = c("+","+")
findrelsubreddit(cursubs,curops,5)

## Creating the ternary plot

# Similatrity to /r/The_Donald
cursubs = c("the_donald")
curops = c("+")
Dsubsims = findrelsubreddit(cursubs,curops,nrow(cursubmat))
# Similarity to /r/SandersforPresident
cursubs = c("sandersforpresident")
curops = c("+")
Ssubsims = findrelsubreddit(cursubs,curops,nrow(cursubmat))
# Similarity to /r/hillaryclinton
cursubs = c("hillaryclinton")
curops = c("+")
Hsubsims = findrelsubreddit(cursubs,curops,nrow(cursubmat))
# List of subreddits we're interested in
ternarysubs = c("theredpill","coontown","fatpeoplehate","politics","worldnews","news","sjwhate","thebluepill","feminism","books","political_revolution","basicincome")
Dternarysubsims = Dsubsims[tolower(names(Dsubsims))%in%ternarysubs]
Sternarysubsims = Ssubsims[tolower(names(Ssubsims))%in%ternarysubs]
Hternarysubsims = Hsubsims[tolower(names(Hsubsims))%in%ternarysubs]
# Normalizing the matrix
allternarysubsims = transform(merge(transform(merge(Sternarysubsims,Dternarysubsims,by="row.names"),row.names=Row.names,Row.names=NULL),Hternarysubsims,by="row.names"),row.names=Row.names,Row.names=NULL)
colnames(allternarysubsims) = c("S","D","H")
allternarysubsimssums = apply(allternarysubsims,1,sum)
allternarysubsimsnorm = sweep(allternarysubsims,1,allternarysubsimssums,"/")
# Creating the plot
pdf("./ternaryplotanno.pdf",height=10,width=10)
ggtern(data=allternarysubsimsnorm,aes(S,D,H)) + geom_point() + geom_text(label=rownames(allternarysubsimsnorm),hjust=0,vjust=0)
dev.off()
pdf("./ternaryplot.pdf",height=10,width=10)
ggtern(data=allternarysubsimsnorm,aes(S,D,H)) + geom_point() + theme_classic()
dev.off()

# Find subreddits that are particularly biased towards any of the three main candidate subreddits
allsubsims = transform(merge(transform(merge(Ssubsims,Dsubsims,by="row.names"),row.names=Row.names,Row.names=NULL),Hsubsims,by="row.names"),row.names=Row.names,Row.names=NULL)
colnames(allsubsims) = c("S","D","H")
chooseunique = c("H") # Set candidate subreddit of interest
curunique = 1/(allsubsims[,(!(colnames(allsubsims)==chooseunique))]/allsubsims[,chooseunique]) # Calculate fold enrichment of target candidate subreddit over other candidate subreddits for all other subreddits
allsubsimsmin = apply(allsubsims,1,min)
curuniquemin = apply(curunique,1,min)
curuniqueminc = curuniquemin[-which(allsubsimsmin==0)]
curuniquemat = data.frame(enrich=curuniqueminc,allsubsims[match(names(curuniqueminc),rownames(allsubsims)),])
curuniquemato = curuniquemat[order(curuniquemat$enrich,decreasing=TRUE),]
curuniquematoc = curuniquemato[which(curuniquemato[,chooseunique]>=0.25),] # Threshold for high enrichment and high raw similarity
head(curuniquematoc,20)


load("z_score_matrix.RData")
nom <- colnames(z_score_matrix)
nom <- do.call(rbind, strsplit(nom, split="\\|"))
colnames(z_score_matrix) <- trimws(nom[,1])
z <- t(z_score_matrix)

mis1 <- apply(z, 1, function(x) sum(is.na(x))/length(x))
mis2 <- apply(z, 2, function(x) sum(is.na(x))/length(x))

summary(mis1)
summary(mis2)

z <- z[, mis2 < 0.75]
z <- apply(z, 2, function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x) })

library(LSAfun)

za <- abs(z)
zl <- lsa(z)

zs <- t(apply(z,1,scale))
dim(zs)
Cosine("Coronary heart disease", "Body mass index", tvectors=z)
Cosine("Age at menarche", "Body mass index", tvectors=z)



get_cosine <- function(v1, v2, mat)
{
	i1 <- which(rownames(mat) == v1)
	i2 <- which(rownames(mat) == v2)
	cosine(mat[i1,], mat[i2,])
}

get_cosine_perm <- function(v1, v2, mat, nperm=1000, scalevec=TRUE)
{
	i1 <- which(rownames(mat) == v1)
	i2 <- which(rownames(mat) == v2)
	v1 <- mat[i1,]
	v2 <- mat[i2,]

	if(scalevec)
	{
		v1 <- as.numeric(scale(v1))
		v2 <- as.numeric(scale(v2))
	}

	emp <- cosine(v1, v2)

	res <- rep(0, nperm)
	for(i in 1:nperm)
	{
		res[i] <- cosine(sample(v1), sample(v2))[1]
	}
	if(emp >= 0)
	{
		r <- sum(res > emp[1]) / nperm
	} else {
		r <- sum(res < emp[1]) / nperm
	}
	return(list(cos=emp, p=r))
}

get_cor_perm <- function(v1, v2, mat, nperm=1000)
{
	i1 <- which(rownames(mat) == v1)
	i2 <- which(rownames(mat) == v2)
	v1 <- mat[i1,]
	v2 <- mat[i2,]

	emp <- cor(v1, v2, use="pair")

	res <- rep(0, nperm)
	for(i in 1:nperm)
	{
		res[i] <- cor(sample(v1), sample(v2), use="pair")
	}
	if(emp >= 0)
	{
		r <- sum(res > emp[1]) / nperm
	} else {
		r <- sum(res < emp[1]) / nperm
	}
	return(list(cos=emp, p=r, dist=res))	
}

cosine(z[3,], z[8,])

cor(scale(z[3,]), scale(z[8,]))

get_cosine("Age at menarche", "Body mass index", z)

get_cor_perm("Age at menarche", "Body mass index", z, nperm=1000)

plot(z[rownames(z) == "Age at menarche",], z[rownames(z) == "Body mass index",])
plot(z[rownames(z) == "Thalamus volume",], z[rownames(z) == "Body mass index",])

summary(lm(z[rownames(z) == "Thalamus volume",] ~ z[rownames(z) == "Body mass index",]))

get_cosine_perm("Thalamus volume", "Body mass index", z)
get_cosine_perm("Thalamus volume", "Thyroid cancer", zs, scalevec=FALSE)

a <- get_cor_perm("Thalamus volume", "Thyroid cancer", z)


cosall <- matrix(0, nrow(z), nrow(z))
cosallp <- matrix(0, nrow(z), nrow(z))
for(i in 1:nrow(z))
{
	message(i)
	for(j in i:nrow(z))
	{
		out <- get_cosine(rownames(z)[i], rownames(z)[j], z)
		cosall[i,j] <- out
		cosall[j,i] <- out
	}
}

heatmap(cosall)

