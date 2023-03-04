library("ggcorrplot")
library("ggplot2")
library("GGally")
library("matrixStats")
library("FactoMineR")
library("cluster")
library("mclust")
library("factoextra")
setwd("~/Documents/WW/BIOE530/final project/HAPT Data Set")

# Load data
X_train = read.table("Train/X_train.txt")
Y_train = unlist(read.table("Train/Y_train.txt"))
subj_id_train = unlist(read.table("Train/subject_id_train.txt"))

X_test = read.table("Test/X_test.txt")
Y_test = unlist(read.table("Test/Y_test.txt"))
subj_id_test = unlist(read.table("Test/subject_id_test.txt"))

activityLabels = read.table("activity_labels.txt")$V2
featureLabels = unlist(read.table("features.txt"))

colnames(X_train) = featureLabels
colnames(X_test) = featureLabels

qqnorm(X_train$`tBodyAcc-Mean-1`)
qqnorm(X_train$`tBodyAcc-Mean-2`)

## Exploratory analysis
Y_all = c(Y_train,Y_test)
Y_all = factor(Y_all)
levels(Y_all) = activityLabels
Y_train = factor(Y_train)
levels(Y_train) = activityLabels
df_allData = data.frame(
  Subject = c(subj_id_train,subj_id_test),
  Activity = Y_all,
  TrainVsTest = factor(c(rep("Train",length(Y_train)),rep("Test",length(Y_test))),levels=c("Train","Test")),
  rbind(X_train,X_test)
)
#correlation matrix
cor.mat <- round(cor(X_train),2)
ggcorrplot(cor.mat, outline.col = "white", tl.cex=0.7)
cor.mat[10,13]


# number of data points per subject
ggplot(df_allData,aes(x=Subject,fill=TrainVsTest)) +
  geom_histogram(color="white") +
  ggtitle("Number of data points per subject")

# number of data points per activity, per subject
ggplot(df_allData,aes(x=Subject,fill=TrainVsTest)) +
  geom_histogram(color="white") +
  ggtitle("Number of data points per subject") +
  facet_wrap(~ Activity, scales="free", nrow = 4)

# look at tBodyAccMean features across subjects and activities
ggplot(df_allData,aes(x=Subject,y=tBodyAccMag.Mean.1,group=Subject,fill=TrainVsTest)) +
  geom_boxplot() +
  ggtitle("tBodyAccMag.Mean.1 per subject") +
  facet_wrap(~ Activity, nrow = 3) #facet_wrap(~ Activity, scales="free", nrow = 4)

## PCA
res.pca = PCA(X_train, ncp=30, graph = FALSE)
print(res.pca)
eig <- res.pca$eig

# individual coordinates
df_trainPCA = data.frame(
  Subject = subj_id_train,
  Activity = Y_train,
  res.pca$ind$coord
)

# Scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# PC1 vs PC2
fviz_pca_ind(res.pca, label="none",habillage=df_trainPCA$Activity,
             addEllipses=TRUE, ellipse.level=0.95)

# pairwise plot of PCs
ggpairs(df_trainPCA[,3:7],aes(color=df_trainPCA$Activity))

### k-means clustering
# perform clustering for k=1:12, with # pcs from 2:2:10, calculating silhouette each time
nCompVec = c(1:5)
kMax = 12 
sil = matrix(nrow = kMax, ncol = length(nCompVec))
for (i in 2:kMax) {
  for (j in 1:length(nCompVec)) {
    kclusterTemp = kmeans(res.pca$ind$coord[,1:nCompVec[j]],centers=i,iter.max = 20, nstart = 25)
    ss = silhouette(kclusterTemp$cluster, dist(res.pca$ind$coord[,1:nCompVec[j]]))
    sil[i,j] <- mean(ss[, 3])
  }
}

sil[1,] = 0
df_silhouette = data.frame(
  Mean_Silhouette = array(sil),
  K = rep(1:kMax,length(nCompVec)),
  Number_of_PCs = factor(rep(nCompVec,each=kMax),levels = nCompVec)
)

ggplot(df_silhouette,aes(x=K,y=Mean_Silhouette,color=Number_of_PCs)) +
  geom_point() +
  geom_line()


# Sublabel activities
staticLabels = c("SITTING","STANDING","LAYING")
dynamicLabels = c("WALKING","WALKING_UPSTAIRS","WALKING_DOWNSTAIRS")
postureLabels = c("STAND_TO_SIT","SIT_TO_STAND","SIT_TO_LIE","LIE_TO_SIT","STAND_TO_LIE","LIE_TO_STAND")
df_trainPCA$Activity3 = array()
df_trainPCA$Activity3[df_trainPCA$Activity %in% staticLabels] = "STATIC"
df_trainPCA$Activity3[df_trainPCA$Activity %in% dynamicLabels] = "DYNAMIC"
df_trainPCA$Activity3[df_trainPCA$Activity %in% postureLabels] = "POSTURE_TRANS"
df_trainPCA$Activity3 = factor(df_trainPCA$Activity3)
fviz_pca_ind(res.pca, label="none",habillage=df_trainPCA$Activity3,
             addEllipses=TRUE, ellipse.level=0.95)


## Perform clustering for k=3, k=12, # pcs from 1:15, calculating ARI each time
nCompMax = 15
k = 12 
ARI_k12 = c()
for (i in 1:nCompMax) {
  kclusterTemp = kmeans(res.pca$ind$coord[,1:i],centers=k,iter.max = 20, nstart = 25)
  ARI_k12[i] <- adjustedRandIndex(kclusterTemp$cluster,df_trainPCA$Activity) 
}

plot(1:nCompMax, ARI_k12, type = "b", pch = 19,
     frame = FALSE, xlab = "Number of princpal components")
abline(v = which.max(ARI), lty = 2)


nCompMax = 15
k = 3 
ARI_k3 = c()
for (i in 1:nCompMax) {
  kclusterTemp = kmeans(res.pca$ind$coord[,1:i],centers=k,iter.max = 20, nstart = 25)
  ARI_k3[i] <- adjustedRandIndex(kclusterTemp$cluster,df_trainPCA$Activity3) 
}

df_ARI = data.frame(
  ARI = c(ARI_k3,ARI_k12),
  K = factor(c(rep(3,nCompMax),rep(12,nCompMax))),
  Number_of_PCs = rep(1:nCompMax,2)
)

ggplot(df_ARI,aes(x=Number_of_PCs,y=ARI,color=K)) +
  geom_point() +
  geom_line()

kcluster1pc = kmeans(res.pca$ind$coord[,1],centers=3,iter.max = 20, nstart = 25)
fviz_pca_ind(res.pca, label="none",habillage=factor(kcluster1pc$cluster),
             addEllipses=TRUE, ellipse.level=0.95)

kcluster2pc = kmeans(res.pca$ind$coord[,1:2],centers=3,iter.max = 20, nstart = 25)
fviz_pca_ind(res.pca, label="none",habillage=factor(kcluster2pc$cluster),
             addEllipses=TRUE, ellipse.level=0.95)

kcluster13pc = kmeans(res.pca$ind$coord[,1:13],centers=12,iter.max = 20, nstart = 25)
fviz_pca_ind(res.pca, label="none",habillage=factor(kcluster13pc$cluster),
             addEllipses=TRUE, ellipse.level=0.95)

# External validation: adjusted Rand Index
df_kmeans = data.frame(
  ARI = array(ARI),
  K = rep(1:kMax,length(nCompVec)),
  Number_of_PCs = factor(rep(nCompVec,each=kMax),levels = nCompVec)
)

ggplot(df_kmeans,aes(x=K,y=ARI,color=Number_of_PCs)) +
  geom_point() +
  geom_line()


# performing clustering for k from 1:12, 15 principal components calculating mean silhouette width EACH TIME
kmax = 12
nComp = 15
sil = rep(0, kmax)
for (i in 2:kmax) {
  kclusterTemp = kmeans(res.pca$ind$coord[,1:nComp],centers=i,iter.max = 20, nstart = 25)
  ss = silhouette(kclusterTemp$cluster, dist(res.pca$ind$coord[,1:nComp]))
  sil[i] <- mean(ss[, 3])
}
# Internal validation: silhouette width
plot(1:kmax, sil, type = "b", pch = 19, 
     frame = FALSE, xlab = "Number of clusters k")
abline(v = which.max(sil), lty = 2)
windows()
plot(ss)

# plot clustering results with k=12, nComp=15
# fviz_pca_ind(res.pca, label="none",habillage=df_trainPCA$Activity,
             # addEllipses=TRUE, ellipse.level=0.95)
# fviz_pca_ind(res.pca, label="none",habillage=kclusterTemp$cluster,
             # addEllipses=TRUE, ellipse.level=0.95)

# # Pairwise plot of PCs with k means clustering
# ggpairs(df_trainPCA[,3:7],aes(color=df_trainPCA$KClusterID))

### Hierarchical clustering

plot(hicluster)

# PCA(X_train, scale.unit = TRUE, ncp = 5, graph = TRUE)
# eig.val <- get_eigenvalue(res.pca)
# head(eig.val)

# hierarchical clustering with 12 clusters, with different linkage methods, varying # pcs
linkageMethods = c("ward.D","single","complete","average")
nCompMax = 15
nClusters = 12

ARI12_hi = matrix(nrow = length(linkageMethods),ncol = nCompMax)
for (i in 1:length(linkageMethods)) {
  for (j in 1:nCompMax) {
    hiclusterTemp <- hclust(dist(res.pca$ind$coord[,1:j]),method = linkageMethods[i])
    clusterCutTemp <- cutree(hiclusterTemp, nClusters)
    ARI12_hi[i,j] <- adjustedRandIndex(clusterCutTemp,df_trainPCA$Activity)
  }
}

df_ARI12_hi = data.frame(
  ARI = array(ARI12_hi),
  Linkage_Method = factor(rep(linkageMethods,nCompMax),levels=linkageMethods),
  Number_of_PCs = rep(1:nCompMax,each=length(linkageMethods))
)

ggplot(df_ARI12_hi,aes(x=Number_of_PCs,y=ARI,color=Linkage_Method)) +
  geom_point() +
  geom_line()

# same thing, external validation with 3 clusters
nClusters = 3
ARI3_hi = matrix(nrow = length(linkageMethods),ncol = nCompMax)
for (i in 1:length(linkageMethods)) {
  for (j in 1:nCompMax) {
    hiclusterTemp <- hclust(dist(res.pca$ind$coord[,1:j]),method = linkageMethods[i])
    clusterCutTemp <- cutree(hiclusterTemp, nClusters)
    ARI3_hi[i,j] <- adjustedRandIndex(clusterCutTemp,df_trainPCA$Activity3)
  }
}

df_ARI3_hi = data.frame(
  ARI = array(ARI3_hi),
  Linkage_Method = factor(rep(linkageMethods,nCompMax),levels=linkageMethods),
  Number_of_PCs = rep(1:nCompMax,each=length(linkageMethods))
)

ggplot(df_ARI3_hi,aes(x=Number_of_PCs,y=ARI,color=Linkage_Method)) +
  geom_point() +
  geom_line()


# plot clustering results for 12 clusters, ward linkage, 8 pcs
hicluster12 <- hclust(dist(res.pca$ind$coord[,1:8]),method="ward.D")
clusterCut12 <- cutree(hicluster12, 12)
fviz_pca_ind(res.pca, label="none",habillage=factor(clusterCut12),
             addEllipses=TRUE, ellipse.level=0.95)

# plot clustering results for 3 clusters, ward linkage, 9 pcs
hicluster3 <- hclust(dist(res.pca$ind$coord[,1:9]),method="ward.D")
clusterCut3 <- cutree(hicluster3, 3)
fviz_pca_ind(res.pca, label="none",habillage=factor(clusterCut3),
             addEllipses=TRUE, ellipse.level=0.95)


# Cophenetic coefficient across PCs for different linkage methods
linkageMethods = c("ward.D","single","complete","average")
nCompMax = 15
cophCorr = matrix(nrow = length(linkageMethods),ncol = nCompMax)
for (i in 1:length(linkageMethods)) {
  for (j in 1:nCompMax) {
    hiclusterTemp <- hclust(dist(res.pca$ind$coord[,1:j]),method = linkageMethods[i])
    coph <- cophenetic(hiclusterTemp)
    cophCorr[i,j] = cor(dist(res.pca$ind$coord[,1:j]),coph)
  }
}

df_Cophenetic = data.frame(
  Cophenetic_Correlation = array(cophCorr),
  Linkage_Method = factor(rep(linkageMethods,nCompMax),levels=linkageMethods),
  Number_of_PCs = rep(1:nCompMax,each=length(linkageMethods))
)

ggplot(df_Cophenetic,aes(x=Number_of_PCs,y=Cophenetic_Correlation,color=Linkage_Method)) +
  geom_point() +
  geom_line()

# Get contributions
coord <- res.pca$var$coord
coord <- as.data.frame(coord)
dim1Contributions <- rownames(coord)[order(coord$Dim.1, decreasing = TRUE)][1:5]
dim2Contributions <- rownames(coord)[order(coord$Dim.2, decreasing = TRUE)][1:5]
dim3Contributions <- rownames(coord)[order(coord$Dim.3, decreasing = TRUE)][1:5]
# order.coord <- coord[order,]

