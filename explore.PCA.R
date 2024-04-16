# PCA with new restrictions on the data:

# ID and then map sites that are outliers:

out.pos <- lapply(df.setup[,v.resp.cols], FindOutliers) %>% unlist %>% unique
outlier.sites <- df.setup$Name[out.pos]

fun.map.DA(outlier.sites)

# set up a df for PCA:

df.PCA <- df.setup
rownames(df.PCA) <- df.PCA[,1] # set the row names as the site names

# determine the top predictors based on variance:

preds.var <- df.PCA[,v.pred.cols] %>% pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  dplyr::summarize(var = var(value)) %>% arrange(desc(var)) %>% ungroup() %>% top_n(16)

# also do variances based on standarized 0-1 variables:

preds.var.01 <- df.PCA[,v.pred.cols] %>%
  lapply(., range01) %>%
  bind_cols() %>%
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  dplyr::summarize(var = var(value)) %>% arrange(desc(var)) %>% ungroup() %>% top_n(16)

# lets look at the differences between the predictor sets for regular and standardized (0-1) top variance:

(x <- preds.var$name)
(y <- preds.var.01$name)

setdiff(x,y) # "mean.daily.Q" "Elev_Avg"     "Elev_Median"  "Sum.WWTP.Q"   "HSG_C" 

setdiff(names(df.setup[,v.pred.cols]),y)

setdiff(y,x) # "R_CROPSNLCD06" "Corn" "R_RIP100_PLANT" "RIP.CSA.100" "BFI"  

# ready for PCA

# there are lots of iterations to be had here:

# 1) all sites with all predictors
# 2) all sites with top half of predictors based on real variance
# 3) all sites with top half of predictors based on 0-1 variance
# 4) remove outlier sites with all predictors
# 5) remove outlier sites with top half of predictors based on real variance
# 6) remove outlier sites with top half of predictors based on 0-1 variance

pc1 <- prcomp(df.PCA[,v.pred.cols], center = TRUE, scale. = TRUE)
pc2 <- prcomp(df.PCA[,v.pred.cols] %>% select(preds.var$name), center = TRUE, scale. = TRUE)
pc3 <- prcomp(df.PCA[,v.pred.cols] %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)
pc4 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name), center = TRUE, scale. = TRUE)
pc5 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var$name), center = TRUE, scale. = TRUE)
pc6 <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)

# make individuals biplots:

l.pca <- list(pc1,pc2,pc3,pc4,pc5,pc6)

l.pca.plot.biplot.obs <- lapply(l.pca, fun.PCA.biplot.individuals)

# look at individuals biplots for each PCA done above,
# as well as the cummulative proportion of explained variance at the second PC
# as well as the PC at which 95% of the variability is explained:

l.pca.plot.biplot.obs[[1]]
as.data.frame(summary(pc1)$importance)$PC2[3]
min(which(summary(pc1)$importance[3,]>.95))

l.pca.plot.biplot.obs[[2]]
as.data.frame(summary(pc2)$importance)$PC2[3]
min(which(summary(pc2)$importance[3,]>.95))

l.pca.plot.biplot.obs[[3]]
as.data.frame(summary(pc3)$importance)$PC2[3]
min(which(summary(pc3)$importance[3,]>.95))

l.pca.plot.biplot.obs[[4]]
as.data.frame(summary(pc4)$importance)$PC2[3]
min(which(summary(pc4)$importance[3,]>.95))

# l.pca.plot.biplot.obs[[5]]
as.data.frame(summary(pc5)$importance)$PC2[3]
min(which(summary(pc5)$importance[3,]>.95))

l.pca.plot.biplot.obs[[6]]
as.data.frame(summary(pc6)$importance)$PC2[3]
min(which(summary(pc6)$importance[3,]>.95))


# i want to make a map based on the three groups I am seeing in pc3 (all sites, restrctions on predictors):

# the first group will be greater than 0 for PC1 and greater than -0.5 for pc2
# the second group will be greater than 0 for PC1 and less than -0.5 for PC2
# the third group will be all sites less than 0 for pc1

df.pc3 <- as.data.frame(pc3$x)

g1<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2>-0.5,])
g2<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2<=-0.5,])
g3<- rownames(df.pc3[df.pc3$PC1<0,])

v.g <- c(g1,g2,g3)
v.g.names <- c(rep('g1', length(g1)),rep('g2', length(g2)),rep('g3', length(g3)))

df.g <- data.frame(v.g, v.g.names) %>% rename(Name = 1, Group = 2)

fun.map.DA(df.g$Name, zcol.key = df.g)

# what if I remove cattargus and genese river and rerun pc3?
# this is basically pc6 with only 2 outlier sites now:

outlier.sites <- c("04231600","04213500")

pc3a <- prcomp(df.PCA[,c(1,v.pred.cols)] %>% filter(!Name %in% outlier.sites) %>% select(-Name) %>% select(preds.var.01$name), center = TRUE, scale. = TRUE)

as.data.frame(summary(pc3a)$importance)$PC2[3]
min(which(summary(pc3a)$importance[3,]>.95))

fun.PCA.biplot.individuals(pc3a)

# the groupings look much different

df.pc3 <- as.data.frame(pc3a$x)

g1<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2>-1,])
g2<- rownames(df.pc3[df.pc3$PC1>0 & df.pc3$PC2<=-1,])
g3<- rownames(df.pc3[df.pc3$PC1<0,])

v.g <- c(g1,g2,g3)
v.g.names <- c(rep('g1', length(g1)),rep('g2', length(g2)),rep('g3', length(g3)))

df.g <- data.frame(v.g, v.g.names) %>% rename(Name = 1, Group = 2)

fun.map.DA(df.g$Name, zcol.key = df.g)

# I dont see much difference between 1 and 2 (I dont know why just three are in 2)
# I suspect it could be something with near stream CSA,
# so I want to look at that map:

df.rip.csa <- df.setup %>% 
  filter(Name %in% df.g$Name) %>% 
  select(Name, R_PLANTNLCD06) %>% 
  rename(Group = 2)

fun.map.DA(df.rip.csa$Name, df.rip.csa)

# so those three stand out because they have much higher watershed percent agriculture

# the issue now becomes we only have three agricultural sites

# lets look at map of all 62 sites (minus the largest ones) or whatever based on agriculture:

outlier.sites <- c(outlier.sites,  c('04260500', '04231600', '01357500', '04249000'))

df.map <- df.datalayers.62 %>% 
  filter(!Name %in% outlier.sites) %>% 
  select(Name, R_PLANTNLCD06) %>% 
  rename(Group = 2)

fun.map.DA(df.map$Name, df.map)

# it doesnt look there are many more sites to add back in that have similar ag to the three I kept...

# lets look at the biplot of predictor variables:

fviz_pca_var(pc3a,
             axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# if we look at the contributions, forest is strongest in the negative direction for Dim1, and plant is strongest in the positive direction
# so can we say that Dim1 is a good contender variable for high correlation with positie nutrient loads?

# lets look at the actual ranks of importance:

df.importance <- as.data.frame(pc3a$rotation) %>% arrange(PC1, desc = F)

