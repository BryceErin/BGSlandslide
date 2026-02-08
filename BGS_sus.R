library(sf)
library(INLA)
library(ggplot2)
library(inlabru)
library(splancs)

#### Data ####
sbdry       <- st_read("Scottishbound.shp")                 # scottish boundary 
covs        <- st_read("SU_Scotland_with_attributes.shp")   # scottish SU covariates with bedrock etc.
DFcount     <- st_read("SU_DFcount.shp")                    # debris flow shapefile


DFcount <- as.data.frame(cbind(covs$ScotsID, DFcount$Join_Count))
colnames(DFcount) <- c("ID", "DFcount")

#### Data prep ####
DFcount$presence     <- ifelse(DFcount$DFcount>0, 1, 0)
data                 <- covs
data$presence        <- DFcount$presence

data$y.count         <- as.numeric(DFcount$DFcount)
y                    <- as.numeric(data$presence)
data$y               <- as.numeric(data$presence)
data$Intercept       <- 1
data$e.pp            <- (covs$area)/1000

data$area.s        <- scale(data$area)
data$LR.s          <- scale(data$LR)
data$slp_avg.s     <- scale(data$slp_avg)
data$slp_std.s     <- scale(data$slp_std)
data$prec_max.s    <- scale(data$prec_max)
data$prec_mean.s    <- scale(data$prec_mean)
data$profC_avg.s   <- scale(data$profC_avg)
data$profC_std.s   <- scale(data$profC_std)


data$quaternary    <- as.factor(data$quaternary)
data$superficia    <- as.factor(data$superficia)
data$bedrockmaj    <- as.factor(data$bedrockmaj)

data$area.g        <- inla.group(data$area.s, method = "quantile", n =20)
data$LR.g          <- inla.group(data$LR.s, method = "quantile", n =20)
data$slp_avg.g     <- inla.group(data$slp_avg.s, method = "quantile", n =20)
data$prec_max.g    <- inla.group(data$prec_max.s, method = "quantile", n =20)
data$prec_mean.g   <- inla.group(data$prec_mean.s, method = "quantile", n =20)
data$profC_avg.g   <- inla.group(data$profC_avg.s, method = "quantile", n = 20)


#### Model components ####
n                      <- rep(1, 153282) 
control.fixed          <- control.fixed(list(prec=0.01), prec.intercept=0.01)
control.family         <- list(link='logit')
control.compute        <- list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE)

dcoords <- data.frame(cbind(data$Longitude, data$Latitude))
dcoords <- dcoords/1000
colnames(dcoords) <- c("x", "y")

# mesh for susceptibility 
dmax.edge    <- diff(range(dcoords$x))/15
dbound.outer <- diff(range(dcoords$x))/100
dmesh        <- inla.mesh.2d(boundary = sbdry,
                             loc=cbind(dcoords$x, dcoords$y),
                             max.edge = c(1,5)*dmax.edge,
                             cutoff = dmax.edge/5,
                             offset = c(dmax.edge, dbound.outer),
                             crs = t.proj)

ggplot() + gg(dmesh)

nv   <- dmesh$n
A    <- inla.spde.make.A(dmesh, loc = as.matrix(dcoords))

spde <- inla.spde2.pcmatern(mesh        = dmesh, 
                            #alpha       = 2, # alpha = 2 corresponds to smoothness = 1 (see Equation 1 in Lindgren at al., 2011, JRSS-B)
                            prior.range = c(1, 0.01), # 10000, 0.1
                            prior.sigma = c(0.01, 0.01)) # 0.05, 0.1 works as well

mesh.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

# effects for susceptibility
effects <- list(c(mesh.index, list(Intercept = 1)),
                list(LR.g       = data$LR.g,
                     slp_avg.g  = data$slp_avg.g,
                     slp_std.s  = data$slp_std.s,
                     prec_max.s = data$prec_max.s, 
                     quaternary = data$quaternary))


# stack for susceptibility
stack.pres <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data$y), 
                         A       = list(A, 1),
                         effects = effects)

A.pred     <- inla.spde.make.A(dmesh, loc = as.matrix(dcoords))

stack.pres.pred <- inla.stack(tag     = 'pred.stack',
                              data    = list(y = NA),
                              A       = list(A.pred, 1),
                              effects = effects)

join.stack.pres <- inla.stack(stack.pres, stack.pres.pred)

# priors for susceptibility
u.LR       <- sd(data$LR.g) 
alpha.LR   <- 0.01
hyper.LR   <- list(theta = list(prior = "pc.prec", param = c(u.LR, alpha.LR)))

u.slp_avg       <- 0.5 
alpha.slp_avg   <- 0.01
hyper.slp_avg   <- list(theta = list(prior = "pc.prec", param = c(u.slp_avg, alpha.slp_avg)))

#### Model fit ####
formula          <- y ~ -1 + Intercept +
                    f(LR.g, model = "rw1", hyper = hyper.LR, cyclic = FALSE, scale.model = TRUE) +
                    f(slp_avg.g, model = "rw1", hyper = hyper.slp_avg, cyclic = FALSE, scale.model = TRUE) +
                    slp_std.s +
                    prec_max.s +
                    f(quaternary, model = "iid", constr = F) +
                    f(spatial.field, model = spde)

inla.setOption(inla.mode = "experimental")
fit              <- inla(formula, 
                         family = "binomial",
                         data = inla.stack.data(join.stack.pres),
                         #control.fixed = control.fixed,
                         control.family = control.family,
                         control.predictor = list(A = inla.stack.A(join.stack.pres), compute = TRUE),
                         control.compute = control.compute, 
                         control.inla  = list(strategy = "simplified.laplace", int.strategy = "eb"),
                         verbose = TRUE)

summary(fit)
save(fit, file = "Fit_BGS_Suscept_DF_Feb23.Rdata")

# checks
plot(fit$summary.random$LR.g$ID, fit$summary.random$LR.g$mean, type = "l")
plot(fit$summary.random$slp_avg$ID, fit$summary.random$slp_avg$mean, type = "l")

# susceptibility fitted values
ID    <- data$ScotsID
igr   <- inla.stack.index(join.stack.pres, 'data.stack')$data

range(fit$summary.fitted.values$mean[igr])
range(inla.link.invlogit(fit$summary.linear.predictor$mean)[igr])


data.pres       <- cbind(ID, fit$summary.fitted.values[igr,])
quantile(data.pres$mean)
data.pres$quant <- data.pres$`0.975quant`- data.pres$`0.025quant`
data.pres       <- data.pres[,-c(4,5, 6, 7)]
colnames(data.pres) <- c("SU_ID", "Mean", "SD", "95CI")

write.csv(data.pres, file = "bgs2_May23_suscept_probs.csv")

#### Fixed effects ####
pdf()
fixed.pres <- fit$summary.fixed[-1,]
inla.object.pres <- fixed.pres
ry <- c(min(inla.object.pres$'0.025quant'), max(inla.object.pres$'0.975quant'))
op <- par(mfrow = c(1,1), mar=c(6,4,3,3))
plot(inla.object.pres$mean, ylim=c(-0.5, 1.2), pch=16, col=2, xlab='', ylab="Regression coefficients",
     main="Fixed Effects", axes=F)
abline(h=0, lty=2, col='grey60')
segments(x0 = 1:nrow(inla.object.pres), y0 = inla.object.pres$'0.025quant', 
         x1 = 1:nrow(inla.object.pres), y1 = inla.object.pres$'0.975quant', col = 4)
axis(2)
axis(1, at = 1:nrow(inla.object.pres), labels = F)
text(x = 1:nrow(inla.object.pres), y = -0.56, labels = c("SLO_SD", "Prec_max"), 
     srt = 90, adj=1.4, xpd=TRUE, cex = 1)
par(op)
dev.off()

#### Random effects ####
xlabs.random = data.frame(LR.g       = sort(unique(inla.group(data$LR, n = 20, method = 'quantile'))), 
                          slp_avg.g  = sort(unique(inla.group(data$slp_avg, n = 20, method = 'quantile'))))

random = fit$summary.random
length(random)
names(random)

random.Bernoulli = c(expression('LR'),
                     expression('SLO_avg'))

pdf(width = 9, height = 11)
op = par(mfrow = c(2,1), mar = c(1,5,6,1), oma = c(4, 4, 1, 1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.1)
for(i in 1:(length(random)-2)){
  inla.object = random[[i]]
  ry          = c(-2, 6)
  id.xlabs    = which(names(random)[i] == names(xlabs.random))
  x           = xlabs.random[, id.xlabs]
# yy = cbind(x,inla.object$mean)
  # which.max(yy[,2])
  # yy[which(yy[,2]>0),1]
  # par(mfrow = c(1,2))
  # plot(inla.object$ID, inla.object$mean, type = 'l', ylim = ry, ylab = 'Regression coefficients', xlab = '', main = '', axes = T)
  # lines(inla.object$ID, inla.object$`0.025quant`, lty = 2)
  # lines(inla.object$ID, inla.object$`0.975quant`, lty = 2)
  plot(x, inla.object$mean, type = 'l', ylim = ry, ylab = 'Regression coefficients', xlab = '', main = '', axes = F)
  polygon(x = c(rev(x), x), y = c(rev(inla.object$`0.025quant`), inla.object$`0.975quant`), col = 'lightblue', border = 'lightblue')
  abline(h = 0, lty = 2, col = 'grey60')
  lines(x, inla.object$mean, col = 'darkblue', lwd = 3)
  title(random.Bernoulli[i], line = 1)
  mtext(paste('Random effects'), side = 3, line = -1.5, outer = TRUE, font = 1, cex = 2)
  axis(2)
  axis(1)
}
par(op)
dev.off()

#### Categorical effects ####
pdf()
cate.pres <- fit$summary.random$quaternary
inla.object.pres <- cate.pres
ry <- c(-4, 4)
op <- par(mfrow = c(1,1), mar=c(8,4,3,3))
plot(inla.object.pres$mean, ylim=ry, pch=16, col=2, xlab='', ylab="Regression coefficients",
     main="Quaternary Effects", axes=F)
abline(h=0, lty=2, col='grey60')
segments(x0 = 1:nrow(inla.object.pres), y0 = inla.object.pres$'0.025quant', 
         x1 = 1:nrow(inla.object.pres), y1 = inla.object.pres$'0.975quant', col = 4)
axis(2)
axis(1, 1:nrow(inla.object.pres), labels = F)
text(x = 1:nrow(inla.object.pres), y = -4.5, labels = c("No_dominant", "Coastal&Estuary", "IceScoured", "Minimal_Till", "Montane&Valley", "Till_Dominant"), 
     srt = 90, adj=1, xpd=TRUE, cex = 1)
par(op)
dev.off()


#### Spatial effects ####
ch                = chull(dcoords)
study.area        = dcoords[c(ch, ch[1]), ]
idx.in.study.area = which(inout(dmesh$loc, study.area)) # mesh nodes in the study area

field = fit$summary.random[['spatial.field']][['mean']]
df    = data.frame(x = dmesh$loc[,1], y= dmesh$loc[,2], mean = field)
df    = df[idx.in.study.area, ]

dev.new()
pmean.nodes = ggplot(data = df, aes(x = x, y = y)) + geom_point(aes(colour = mean)) + 
  scale_colour_gradientn(colours = terrain.colors(10)) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + theme_classic()

field = fit$summary.random[['spatial.field']][['sd']]
df    = data.frame(x = dmesh$loc[,1], y= dmesh$loc[,2], sd = field)
df    = df[idx.in.study.area, ]

psd.nodes = ggplot(data = df, aes(x = x, y = y)) + geom_point(aes(colour = sd)) + 
  scale_colour_gradientn(colours = terrain.colors(10)) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + theme_classic()

pdf()
gridExtra::grid.arrange(pmean.nodes, psd.nodes, nrow = 1, top = 'Posterior mean and SD spatial field')
dev.off()

# for ArcMap plotting
SusLSE <- fit$summary.random$spatial.field
SusLSE$CI <- SusLSE$`0.975quant`-SusLSE$`0.025quant`
SusLSE <- SusLSE[,c(2,9)]
SusLSE$x <- dmesh$loc[,1]
SusLSE$y <- dmesh$loc[,2]

write.csv(SusLSE, file = "SusLSE.csv")

#### ROC ####
library(pROC)
out <- plot.roc(data$y, inla.link.invlogit(fit$summary.linear.predictor$mean)[igr])

dev.new(width = 7, height = 7)
plot((1-out$specificities), out$sensitivities, type='l', lwd=5, main="ROC curve", 
     ylab="Sensitivity (TPR)", xlab="1-Specificity (FPR)", xaxt='n', yaxt='n',
     col="darkblue", font.main=1)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)
polygon(c(1,(1-out$specificities)), c(0,out$sensitivities), col="lightblue")
text(0.8, 0.2, "AUC = 0.9763")
dev.off()


#### CV ####

df <- data.frame(cbind(covs$ScotsID, dcoords))
colnames(df) <- c("ScotsID", "x", "y")

library(sperrorest)
samp  <- partition_kmeans(df, coords = c(2,3), nfold = 2, repetition = 1, seed1 = 1)
folds <- samp[[1]]


# Fit for folds
fold1     <- folds$`1`
train1ID  <- fold1$train
test1ID   <- fold1$test
test1NA   <- data[which(data$ScotsID %in% test1ID), ]
test1NA$y <- NA

train1    <- data[which(data$ScotsID %in% train1ID), ]
data1     <- rbind(train1, test1NA)
data1     <- data1[order(data1$ScotsID),]

stack.pres <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data1$y),
                         A       = list(A, 1),
                         effects = effects)

join.stack.pres <- inla.stack(stack.pres, stack.pres.pred)


inla.setOption(inla.mode = "experimental")
fitFold1          <- inla(formula, 
                         family = "binomial",
                         data = inla.stack.data(join.stack.pres),
                         #control.fixed = control.fixed,
                         control.family = control.family,
                         control.predictor = list(A = inla.stack.A(join.stack.pres), compute = TRUE),
                         control.compute = control.compute, 
                         control.inla  = list(strategy = "simplified.laplace", int.strategy = "eb"),
                         verbose = TRUE)
summary(fitFold1)

# checks
plot(fitFold1$summary.random$LR.g$ID, fitFold1$summary.random$LR.g$mean, type = "l")
plot(fitFold1$summary.random$slp_avg.g$ID, fitFold1$summary.random$slp_avg.g$mean, type = "l")

ID    <- data1$ScotsID
igr   <- inla.stack.index(join.stack.pres, 'data.stack')$data
igrp  <- inla.stack.index(join.stack.pres, 'pred.stack')$data

range(fitFold1$summary.fitted.values$mean[igr])
range(fitFold1$summary.fitted.values$mean[igrp])
range(inla.link.invlogit(fitFold1$summary.linear.predictor$mean)[igr])
range(inla.link.invlogit(fitFold1$summary.linear.predictor$mean)[igrp])

lp   <- fitFold1$summary.linear.predictor[igr,]
p    <- inla.link.invlogit(lp)
p    <- cbind(data1$ScotsID, p)

p1   <- p[which(p$`data1$ScotsID` %in% test1ID),]   # test SUs ID from fold predicted DF probability of occ.
quantile(p1$mean) # range seems slightly different from original but consistent

# compare with model fit 
lp   <- fit$summary.linear.predictor[igr,]
p    <- inla.link.invlogit(lp)
p    <- cbind(ID, p)

p_fit <- p[which(p$ID %in% test1ID),] # test SUs ID from fitted model DF probability of occ.
quantile(p_fit$mean)

plot(p_fit$mean, p1$mean, main = "Fold 1 IDs", xlab = "P(occ) from fitted model", 
     ylab = "P(occ) from fold 1 NA model", xlim = c(0, 0.6))
abline(0,1)

## ROC curves
p    <- data.frame(ID, p1 = fitFold1$summary.fitted.values$mean[igr])
p1   <- p[which(p$ID %in% test1ID),] 

out1      <- plot.roc(data$y, inla.link.invlogit(fitFold1$summary.linear.predictor$mean)[igr])
out1      <- plot.roc(data$y, fitFold1$summary.fitted.values$mean[igr])
out1      <- plot.roc(test1NA$presence, p1$p1)
## repeat lines 302-371 for all folds 

pdf()
dev.new(width=7, height=7)
plot((1-out1$specificities), out1$sensitivities, type='l', lwd=2, main="10 Fold: ROC curves", 
     ylab="Sensitivity (TPR)", xlab="1-Specificity (FPR)", xaxt='n', yaxt='n',
     col="darkblue", font.main=1)
lines((1-out2$specificities), out2$sensitivities, lwd=2, col="red")
lines((1-out3$specificities), out3$sensitivities, lwd=2, col="black")
lines((1-out4$specificities), out4$sensitivities, lwd=2, col="pink")
lines((1-out5$specificities), out5$sensitivities, lwd=2, col="darkgreen")
lines((1-out6$specificities), out6$sensitivities, lwd=2, col="orange")
lines((1-out7$specificities), out7$sensitivities, lwd=2, col="seagreen")
lines((1-out8$specificities), out8$sensitivities, lwd=2, col="grey")
lines((1-out9$specificities), out9$sensitivities, lwd=2, col="purple")
lines((1-out10$specificities), out10$sensitivities, lwd=2, col="turquoise")

axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)

text(0.8, 0.8, "Fold 1 AUC = 0.9765", col = "darkblue")
text(0.8, 0.75, "Fold 2 AUC = 0.9704", col = "red")
text(0.8, 0.7, "Fold 3 AUC = 0.9763", col = "black")
text(0.8, 0.65, "Fold 4 AUC = 0.9764", col = "pink")
text(0.8, 0.6, "Fold 5 AUC = 0.9764", col = "darkgreen")
text(0.8, 0.55, "Fold 6 AUC = 0.9766", col = "orange")
text(0.8, 0.5, "Fold 7 AUC = 0.9765", col = "seagreen")
text(0.8, 0.45, "Fold 8 AUC = 0.9764", col = "grey")
text(0.8, 0.4, "Fold 9 AUC = 0.9760", col = "purple")
text(0.8, 0.35, "Fold 10 AUC = 0.9764", col = "turquoise")
dev.off()

auc10 <- out10$auc

## group CV
set.seed(123)
n <- 100
x <- rnorm(n)
y <- 2 * x + rnorm(n)
dat <- data.frame(ID = c(1:n), x,y)

linear_model <- inla(y ~ x, data = dat, family = "gaussian")

groups <- list(1:10)  # In this simple example, we use individual observations as groups

# Perform Leave-Group-Out Cross-Validation
cv_results <- inla.group.cv(
  result = linear_model,
  groups = groups)



?inla.group.cv
loocv_res <- inla.group.cv(result = fit)
ULOOCV    <- mean(loocv_res$cv[1:153282])

lgocv_res <- inla.group.cv(result = fit, num.level.sets = 1)  # groups as 10 fold # try 10
ULGOCV    <- mean(lgocv_res$cv[1:153282])

groups <- lapply(1:n,FUN = function(i){which(id == id[i])})
lgocv_res = inla.group.cv(result = res,groups = groups)

# ROC for LGOCV
library(pROC)
eeta.bern <- exp(lgocv_res$mean)[1:153282]
pr.bern   <- eeta.bern/(1 + eeta.bern)
out       <- plot.roc(data$presence, pr.bern, Show.labels = F, returnSensitivityMat = T)
out$auc
dev.new(width = 7, height = 7)
plot((1-out$specificities), out$sensitivities, type='l', lwd=5, main="ROC curve", 
     ylab="Sensitivity (TPR)", xlab="1-Specificity (FPR)", xaxt='n', yaxt='n',
     col="darkblue", font.main=1)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)
axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels = c(0.0,0.2,0.4,0.6,0.8,1.0),
     tick=T, lwd=1)
polygon(c(1,(1-out$specificities)), c(0,out$sensitivities), col="lightblue")
text(0.8, 0.2, "AUC = 0.9753")
dev.off()






#### residuals ####
## if say test SUs from whole model fit are 'observed' then we can say test SUs from foldNA models are 'predicted'
## then can compute residuals = observed - predicted
## to test model

# scaling residual = raw residual
# pearson residual = raw residual / variance

#fold residuals
s_resid <- p_fold10 - p10
s_resid$ID <- p_fold10$ID
s_resid10 <- s_resid[,c(1,2,3,4,5,6)]
colnames(s_resid10) <- c("ID", "s_mean", "s_sd", "s_q0.025", "s_q0.5", "s_q0.975")

p_resid <- s_resid10 / p10$sd  
p_resid10 <- p_resid[,c(2,3,4,5,6)]
colnames(p_resid10) <- c("p_mean", "p_sd", "p_q0.025", "p_q0.5", "p_q0.975")

fold10_resid <- cbind(s_resid10, p_resid10)
write.csv(fold10_resid, file = "fold10_resid.csv")

#fold10_resid <- read.csv("fold10_resid.csv", header = T, sep = ",")

#summary stats for residuals
## scaling
s_mean10 <- mean(fold10_resid$s_mean)
scaling_mean <- mean(s_mean1, s_mean2, s_mean3, s_mean4, s_mean5, s_mean6, s_mean7, s_mean8, s_mean9, s_mean10)

p_mean10 <- mean(fold10_resid$p_mean)
pearson_mean <- mean(p_mean1, p_mean2, p_mean3, p_mean4, p_mean5, p_mean6, p_mean7, p_mean8, p_mean9, p_mean10)

s_sd10 <- mean(fold10_resid$s_sd)
scaling_sd <- mean(s_sd1, s_sd2, s_sd3, s_sd4, s_sd5, s_sd6, s_sd7, s_sd8, s_sd9, s_sd10)

p_sd10 <- mean(fold10_resid$p_sd)
pearson_sd <- mean(p_sd1, p_sd2, p_sd3, p_sd4, p_sd5, p_sd6, p_sd7, p_sd8, p_sd9, p_sd10)

s_q0.025_10 <- mean(fold10_resid$s_q0.025)
sacling_q0.025 <- mean(s_q0.025_1, s_q0.025_2, s_q0.025_3, s_q0.025_4, s_q0.025_5, s_q0.025_6, s_q0.025_7, s_q0.025_8,
                       s_q0.025_9, s_q0.025_10)

p_q0.025_10 <- mean(fold10_resid$p_q0.025)
pearson_q0.025 <- mean(p_q0.025_1, p_q0.025_2, p_q0.025_3, p_q0.025_4, p_q0.025_5, p_q0.025_6, p_q0.025_7, p_q0.025_8,
                       p_q0.025_9, p_q0.025_10)

s_q0.5_10 <- mean(fold10_resid$s_q0.5)
scaling_q0.5 <- mean(s_q0.5_1, s_q0.5_2, s_q0.5_3, s_q0.5_4, s_q0.5_5, s_q0.5_6, s_q0.5_7, s_q0.5_8,
                       s_q0.5_9, s_q0.5_10)

p_q0.5_10 <- mean(fold10_resid$p_q0.5)
pearson_q0.5 <- mean(p_q0.5_1, p_q0.5_2, p_q0.5_3, p_q0.5_4, p_q0.5_5, p_q0.5_6, p_q0.5_7, p_q0.5_8,
                     p_q0.5_9, p_q0.5_10)

s_q0.975_10 <- mean(fold10_resid$s_q0.975)
scaling_q0.975 <- mean(s_q0.975_1, s_q0.975_2, s_q0.975_3, s_q0.975_4, s_q0.975_5, s_q0.975_6, s_q0.975_7, s_q0.975_8,
                     s_q0.975_9, s_q0.975_10)

p_q0.975_10 <- mean(fold10_resid$p_q0.975)
pearson_q0.975 <- mean(p_q0.975_1, p_q0.975_2, p_q0.975_3, p_q0.975_4, p_q0.975_5, p_q0.975_6, p_q0.975_7, p_q0.975_8,
                       p_q0.975_9, p_q0.975_10)


fixed <- mbarrier$summary.fixed[3:8,]
inla.object <- fixed
ry <- c(min(inla.object$'0.025quant'), max(inla.object$'0.975quant'))
op <- par(mfrow = c(1,1), mar=c(8,4,1,1))
pdf()
plot(inla.object$mean, ylim = c(-0.5, 0.5), pch=16, cex = 2, col=2, xlab='', ylab="Regression coefficients",
     main="Fixed Effects", axes=F)
abline(h=0, lty=2, col='grey60')
segments(x0 = 1:nrow(inla.object), y0 = inla.object$'0.025quant', 
         x1 = 1:nrow(inla.object), y1 = inla.object$'0.975quant', col = 4, lwd = 4)
axis(2)
axis(1, at = 1:nrow(inla.object), labels = F)
text(x = 1:nrow(inla.object), y = -0.6, labels = c("SI_SD", "SI_SD (Copy)", "Temp_SD", "Temp_SD (Copy)", "Rain_SD", "Rain_SD (Copy)"), 
     srt = 90, adj=1, xpd=TRUE, cex = 1)
par(op)
dev.off()


