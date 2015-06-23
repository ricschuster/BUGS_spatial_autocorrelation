#######################################################################
######################## START OF DEJU ################################
#######################################################################

#shared code
library(unmarked)
library(MuMIn)
library(R2OpenBUGS)
library(ncf)
library(ROCR)

DEJU1<- read.csv(file="DEJU.csv")
covar <- read.csv((file="Covariates_PntCntmaster_2005_11_SHR_WET_HRB_covs.csv"))
surcov <- read.csv((file="Cov_survey_obs_2011.csv"))
reform <- read.csv((file="Reformatted_Birds_OCODES.csv"))
index <- read.csv((file="Sitecov_Index.csv"))
avgcov <- read.csv((file="Avg_Cov_PntCntmaster_2005_11_SHR_WET_HRB_covs.csv"))
UTM <- read.csv((file="UTM.csv"))

for (i in 1:length(names(covar))){
  assign(names(covar)[i],covar[,i])
}
#Spatial Autocorrelation Distances
#dist.indx.20k <- read.csv((file="Dist_Index_20km.csv"))
#distance.20k <- read.csv((file="Dist_Data_20km.csv"))
#dist.indx.10k <- read.csv((file="Dist_Index_10km.csv"))
#distance.10k <- read.csv((file="Dist_Data_10km.csv"))
#dist.indx.5k <- read.csv((file="Dist_Index_5km.csv"))
#distance.5k <- read.csv((file="Dist_Data_5km.csv"))
#dist.indx.2.5k <- read.csv((file="Dist_Index_2_5km.csv"))
#distance.2.5k <- read.csv((file="Dist_Data_2_5km.csv"))
#dist.indx.1k <- read.csv((file="Dist_Index_1km.csv"))
#distance.1k <- read.csv((file="Dist_Data_1km.csv"))
#dist.indx.150m <- read.csv((file="Dist_Index_150m.csv"))
#distance.150m <- read.csv((file="Dist_Data_150m.csv"))
#dist.indx.200m <- read.csv((file="Dist_Index_200m.csv"))
#distance.200m <- read.csv((file="Dist_Data_200m.csv"))
dist.indx.250m <- read.csv((file="Dist_Index_250m.csv"))
distance.250m <- read.csv((file="Dist_Data_250m.csv"))
dist.indx.300m <- read.csv((file="Dist_Index_300m.csv"))
distance.300m <- read.csv((file="Dist_Data_300m.csv"))
#dist.indx.400m <- read.csv((file="Dist_Index_400m.csv"))
#distance.400m <- read.csv((file="Dist_Data_400m.csv"))
#dist.indx.500m <- read.csv((file="Dist_Index_500m.csv"))
#distance.500m <- read.csv((file="Dist_Data_500m.csv"))

TIME <- reform$TIME
OCODES <- reform$OCODES
DATE_JUL <- reform$DATE_JUL
index <- as.matrix(index)

###############################################################################
################ NO Autocorrelation BUGS model ################################
###############################################################################

sink("DEJU_ocodes_avg_res.txt")
cat("
model {

# Priors for p
# values come from unmarked output
beta.p.CR_CL~dnorm(5.38078624908733E-02,154.548897221227)
beta.p.DATE_JUL~dnorm(0.14641908343854,277.84019599501)
beta.p.TIME~dnorm(-9.77688378966269E-03,228.498508306325)
alpha.p[1]~dnorm(-0.805393236528632,30.7283570729493)
alpha.p[2]~dnorm(0.243279206919558,19.1543223352186)
alpha.p[3]~dnorm(-0.15496550916313,18.5349093071377)
alpha.p[4]~dnorm(0.789188376455311,22.7681418280675)
alpha.p[5]~dnorm(0.934257974198359,18.5741355499449)
alpha.p[6]~dnorm(0.159591498787112,16.4741347175178)

# Priors for psi
alpha.psi~dnorm(2.95254976514922,2.03110938047235)
beta.psi.near_roadZ~dnorm(-1.10366263303045,5.17436380822411)
beta.psi.near_saltwZ~dnorm(0.900995243643107,7.60647348295563)
beta.psi.RUR_1KMZ~dnorm(1.42836931069341,1.53582717243546)
beta.psi.FOR1_1KMZ~dnorm(-1.07044673987616,7.25170249185306)
beta.psi.FOR2_1KMZ~dnorm(0.423902077324346,4.89046418279917)
beta.psi.HRB_1KMZ~dnorm(1.68232587513669,1.08435081620085)
beta.psi.SAV_1KMZ~dnorm(0.829813212073497,5.38815738900914)
beta.psi.SHR_1KMZ~dnorm(-1.37596366747437,5.47937625316429)
beta.psi.WET_1KMZ~dnorm(-0.677456826482116,15.0618919262386)
beta.psi.URB_100Z~dnorm(-0.594299947241002,17.9811672763755)
beta.psi.HRB_100Z~dnorm(-0.432572236959316,8.83143012696912)
beta.psi.SHR_100Z~dnorm(0.189870862915399,10.776327188821)
beta.psi.Near_UrbZ~dnorm(1.65013017680579,4.17971177071632)
beta.psi.Is_sizeZ~dnorm(1.04213418448955,6.50561013022006)
beta.psi.rd_up_100Z~dnorm(1.18264351819413,6.4517744006783)
beta.psi.CR_CLZ~dnorm(0.875218861094942,10.7662324891414)

# Likelihood
for (i in 1:R) {
   # True state model for the partially observed true state
   z[i] ~ dbern(psi[i])             # True occupancy z at site i

   logit(psi[i]) <- alpha.psi + beta.psi.near_roadZ * near_roadZ[i] + beta.psi.near_saltwZ * near_saltwZ[i] + beta.psi.RUR_1KMZ * RUR_1KMZ[i] + beta.psi.FOR1_1KMZ * FOR1_1KMZ[i] + beta.psi.FOR2_1KMZ * FOR2_1KMZ[i] + beta.psi.HRB_1KMZ * HRB_1KMZ[i] + beta.psi.SAV_1KMZ * SAV_1KMZ[i] + beta.psi.SHR_1KMZ * SHR_1KMZ[i] + beta.psi.WET_1KMZ * WET_1KMZ[i] + beta.psi.URB_100Z * URB_100Z[i] + beta.psi.HRB_100Z * HRB_100Z[i] + beta.psi.SHR_100Z * SHR_100Z[i] + beta.psi.Near_UrbZ * Near_UrbZ[i] + beta.psi.Is_sizeZ * Is_sizeZ[i] + beta.psi.rd_up_100Z * rd_up_100Z[i] + beta.psi.CR_CLZ * CR_CLZ[i]

   for (j in indx[i,3]:(indx[i,3]+indx[i,2]-1)) {
      # Observation model for the actual observations
      y[j] ~ dbern(p.eff[j])    # Detection-nondetection at i and j
      p.eff[j] <- z[i] * p[j]
      
      logit(p[j]) <- alpha.p[OCODES[j]] + beta.p.CR_CL * CR_CLZ[i] + beta.p.DATE_JUL * DATE_JUL[j] + beta.p.TIME * TIME[j]
      } #j

   yres[i] <- psi[i] - z[i]
   } #i

# Derived quantities
occ.fs <- sum(z[])       # Number of occupied sites among those studied

}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = reform$DEJU, indx = index[,1:3], CR_CLZ = CR_CLZ, DATE_JUL = reform$DATE_JUL, TIME = reform$TIME, OCODES = as.integer(OCODES), near_roadZ=near_roadZ, near_saltwZ=near_saltwZ, RUR_1KMZ=RUR_1KMZ, FOR1_1KMZ=FOR1_1KMZ, FOR2_1KMZ=FOR2_1KMZ, HRB_1KMZ=HRB_1KMZ, SAV_1KMZ=SAV_1KMZ, SHR_1KMZ=SHR_1KMZ, WET_1KMZ=WET_1KMZ, URB_100Z=URB_100Z, HRB_100Z=HRB_100Z, SHR_100Z=SHR_100Z, Near_UrbZ=Near_UrbZ, Is_sizeZ=Is_sizeZ, rd_up_100Z=rd_up_100Z, R = length(CR_CLZ))
    
# Initial values
#zst <- index[,4]   #Good inits for latent states essential
zst <- apply(DEJU1, 1, max, na.rm=T) #Good inits for latent states essential
inits <- function(){list(beta.p.CR_CL=runif(1,-3,3), beta.p.DATE_JUL=runif(1,-3,3), beta.p.TIME=runif(1,-3,3), alpha.p=runif(6,-3,3),alpha.psi=runif(1,-3,3), beta.psi.near_roadZ=runif(1,-3,3), beta.psi.near_saltwZ=runif(1,-3,3), beta.psi.RUR_1KMZ=runif(1,-3,3), beta.psi.FOR1_1KMZ=runif(1,-3,3), beta.psi.FOR2_1KMZ=runif(1,-3,3), beta.psi.HRB_1KMZ=runif(1,-3,3), beta.psi.SAV_1KMZ=runif(1,-3,3), beta.psi.SHR_1KMZ=runif(1,-3,3), beta.psi.WET_1KMZ=runif(1,-3,3), beta.psi.URB_100Z=runif(1,-3,3), beta.psi.HRB_100Z=runif(1,-3,3), beta.psi.SHR_100Z=runif(1,-3,3), beta.psi.Near_UrbZ=runif(1,-3,3), beta.psi.Is_sizeZ=runif(1,-3,3), beta.psi.rd_up_100Z=runif(1,-3,3), beta.psi.CR_CLZ=runif(1,-3,3),z = zst)}

# Parameters monitored
params <- c("yres", "psi", "z", "beta.p.CR_CL", "beta.p.DATE_JUL", "beta.p.TIME", "alpha.p","alpha.psi", "beta.psi.near_roadZ", "beta.psi.near_saltwZ", "beta.psi.RUR_1KMZ", "beta.psi.FOR1_1KMZ", "beta.psi.FOR2_1KMZ", "beta.psi.HRB_1KMZ", "beta.psi.SAV_1KMZ", "beta.psi.SHR_1KMZ", "beta.psi.WET_1KMZ", "beta.psi.URB_100Z", "beta.psi.HRB_100Z", "beta.psi.SHR_100Z", "beta.psi.Near_UrbZ", "beta.psi.Is_sizeZ", "beta.psi.rd_up_100Z", "beta.psi.CR_CLZ",  "occ.fs")

# MCMC settings
ni <- 1100
nt <- 3
nb <- 100
nc <- 3

# Call OpenBUGS from R (BRT 1 min)
out <- bugs(win.data, inits, params, "DEJU_ocodes_avg_res.txt", n.chains = nc,
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())

print(out, dig = 3)
# NO autcov 


###############################################################################
################ Autocorrelation BUGS model ###################################
###############################################################################

sink("DEJU_ocodes_avg_autcov_res.txt")
cat("
model {

# Priors for p
beta.p.CR_CL~dnorm(5.38078624908733E-02,154.548897221227)
beta.p.DATE_JUL~dnorm(0.14641908343854,277.84019599501)
beta.p.TIME~dnorm(-9.77688378966269E-03,228.498508306325)
alpha.p[1]~dnorm(-0.805393236528632,30.7283570729493)
alpha.p[2]~dnorm(0.243279206919558,19.1543223352186)
alpha.p[3]~dnorm(-0.15496550916313,18.5349093071377)
alpha.p[4]~dnorm(0.789188376455311,22.7681418280675)
alpha.p[5]~dnorm(0.934257974198359,18.5741355499449)
alpha.p[6]~dnorm(0.159591498787112,16.4741347175178)

# Priors for psi
alpha.psi~dnorm(2.95254976514922,2.03110938047235)
beta.psi.near_roadZ~dnorm(-1.10366263303045,5.17436380822411)
beta.psi.near_saltwZ~dnorm(0.900995243643107,7.60647348295563)
beta.psi.RUR_1KMZ~dnorm(1.42836931069341,1.53582717243546)
beta.psi.FOR1_1KMZ~dnorm(-1.07044673987616,7.25170249185306)
beta.psi.FOR2_1KMZ~dnorm(0.423902077324346,4.89046418279917)
beta.psi.HRB_1KMZ~dnorm(1.68232587513669,1.08435081620085)
beta.psi.SAV_1KMZ~dnorm(0.829813212073497,5.38815738900914)
beta.psi.SHR_1KMZ~dnorm(-1.37596366747437,5.47937625316429)
beta.psi.WET_1KMZ~dnorm(-0.677456826482116,15.0618919262386)
beta.psi.URB_100Z~dnorm(-0.594299947241002,17.9811672763755)
beta.psi.HRB_100Z~dnorm(-0.432572236959316,8.83143012696912)
beta.psi.SHR_100Z~dnorm(0.189870862915399,10.776327188821)
beta.psi.Near_UrbZ~dnorm(1.65013017680579,4.17971177071632)
beta.psi.Is_sizeZ~dnorm(1.04213418448955,6.50561013022006)
beta.psi.rd_up_100Z~dnorm(1.18264351819413,6.4517744006783)
beta.psi.CR_CLZ~dnorm(0.875218861094942,10.7662324891414)

beta.occ.spat ~ dunif(-10, 10)
# Likelihood
for (i in 1:R) {
   # True state model for the partially observed true state
   z[i] ~ dbern(psi[i])             # True occupancy z at site i

   x[i,1] <- 0
   w[i,1] <- 0
   for (k in dist.ind[i,3]:(dist.ind[i,3] + dist.ind[i,2]-1)) {
      # weigth's for spatial autocorrelation
      x[i,(k-dist.ind[i,3]+2)] <- x[i,(k-dist.ind[i,3]+1)] + z[distance[k,2]]/distance[k,3]
      w[i,(k-dist.ind[i,3]+2)] <- w[i,(k-dist.ind[i,3]+1)] + 1/distance[k,3]
      } #k
   autcov[i] <- x[i,(dist.ind[i,2]+1)]/w[i,(dist.ind[i,2]+1)]

   logit(psi[i]) <- alpha.psi + beta.psi.near_roadZ * near_roadZ[i] + beta.psi.near_saltwZ * near_saltwZ[i] + beta.psi.RUR_1KMZ * RUR_1KMZ[i] + beta.psi.FOR1_1KMZ * FOR1_1KMZ[i] + beta.psi.FOR2_1KMZ * FOR2_1KMZ[i] + beta.psi.HRB_1KMZ * HRB_1KMZ[i] + beta.psi.SAV_1KMZ * SAV_1KMZ[i] + beta.psi.SHR_1KMZ * SHR_1KMZ[i] + beta.psi.WET_1KMZ * WET_1KMZ[i] + beta.psi.URB_100Z * URB_100Z[i] + beta.psi.HRB_100Z * HRB_100Z[i] + beta.psi.SHR_100Z * SHR_100Z[i] + beta.psi.Near_UrbZ * Near_UrbZ[i] + beta.psi.Is_sizeZ * Is_sizeZ[i] + beta.psi.rd_up_100Z * rd_up_100Z[i] + beta.psi.CR_CLZ * CR_CLZ[i] + beta.occ.spat * autcov[i]

   for (j in indx[i,3]:(indx[i,3]+indx[i,2]-1)) {
      # Observation model for the actual observations
      y[j] ~ dbern(p.eff[j])    # Detection-nondetection at i and j
      p.eff[j] <- z[i] * p[j]
      
      logit(p[j]) <- alpha.p[OCODES[j]] + beta.p.CR_CL * CR_CLZ[i] + beta.p.DATE_JUL * DATE_JUL[j] + beta.p.TIME * TIME[j]
      } #j

   yres[i] <- psi[i] - z[i]
   } #i

# Derived quantities
occ.fs <- sum(z[])       # Number of occupied sites among those studied

}
",fill = TRUE)
sink()

# Choose Distance values to use in WinBUGS
# 20k 10k 5k 2.5k 1k 500m 400m 300m 250m 200m 150m
dist.indx <- as.matrix(dist.indx.250m[2:4])
distance <- as.matrix(distance.250m)

# Bundle data
win.data <- list(y = reform$DEJU, indx = index[,1:3], CR_CLZ = CR_CLZ, DATE_JUL = reform$DATE_JUL, TIME = reform$TIME, OCODES = as.integer(OCODES), near_roadZ=near_roadZ, near_saltwZ=near_saltwZ, RUR_1KMZ=RUR_1KMZ, FOR1_1KMZ=FOR1_1KMZ, FOR2_1KMZ=FOR2_1KMZ, HRB_1KMZ=HRB_1KMZ, SAV_1KMZ=SAV_1KMZ, SHR_1KMZ=SHR_1KMZ, WET_1KMZ=WET_1KMZ, URB_100Z=URB_100Z, HRB_100Z=HRB_100Z, SHR_100Z=SHR_100Z, Near_UrbZ=Near_UrbZ, Is_sizeZ=Is_sizeZ, rd_up_100Z=rd_up_100Z, R = length(CR_CLZ), dist.ind = dist.indx, distance = distance)

# Initial values
#zst <- index[,4]   #Good inits for latent states essential
zst <- apply(DEJU1, 1, max, na.rm=T) #Good inits for latent states essential
inits <- function(){list(beta.p.CR_CL=runif(1,-3,3), beta.p.DATE_JUL=runif(1,-3,3), beta.p.TIME=runif(1,-3,3), alpha.p=runif(6,-3,3),alpha.psi=runif(1,-3,3), beta.psi.near_roadZ=runif(1,-3,3), beta.psi.near_saltwZ=runif(1,-3,3), beta.psi.RUR_1KMZ=runif(1,-3,3), beta.psi.FOR1_1KMZ=runif(1,-3,3), beta.psi.FOR2_1KMZ=runif(1,-3,3), beta.psi.HRB_1KMZ=runif(1,-3,3), beta.psi.SAV_1KMZ=runif(1,-3,3), beta.psi.SHR_1KMZ=runif(1,-3,3), beta.psi.WET_1KMZ=runif(1,-3,3), beta.psi.URB_100Z=runif(1,-3,3), beta.psi.HRB_100Z=runif(1,-3,3), beta.psi.SHR_100Z=runif(1,-3,3), beta.psi.Near_UrbZ=runif(1,-3,3), beta.psi.Is_sizeZ=runif(1,-3,3), beta.psi.rd_up_100Z=runif(1,-3,3), beta.psi.CR_CLZ=runif(1,-3,3),z = zst, beta.occ.spat = runif(1, -3, 3))}

# Parameters monitored
params <- c("yres", "psi", "z",  "beta.p.CR_CL", "beta.p.DATE_JUL", "beta.p.TIME", "alpha.p","alpha.psi", "beta.psi.near_roadZ", "beta.psi.near_saltwZ", "beta.psi.RUR_1KMZ", "beta.psi.FOR1_1KMZ", "beta.psi.FOR2_1KMZ", "beta.psi.HRB_1KMZ", "beta.psi.SAV_1KMZ", "beta.psi.SHR_1KMZ", "beta.psi.WET_1KMZ", "beta.psi.URB_100Z", "beta.psi.HRB_100Z", "beta.psi.SHR_100Z", "beta.psi.Near_UrbZ", "beta.psi.Is_sizeZ", "beta.psi.rd_up_100Z", "beta.psi.CR_CLZ", "beta.occ.spat", "occ.fs")

# MCMC settings
ni <- 1100
nt <- 3
nb <- 100
nc <- 3

# Call OpenBUGS from R (BRT 1 min)
out_spat <- bugs(win.data, inits, params, "DEJU_ocodes_avg_autcov_res.txt", n.chains = nc,
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())

print(out_spat, dig = 3)
#250m autcov 

############################################################################
###################### Moran's I ###########################################
############################################################################
res <- as.data.frame(UTM[,3:4])
bayres<- correlog(x = res$UTM_E, y = res$UTM_N, z = as.vector(out$mean$yres), increment=300, resamp=10, quiet = TRUE)
sbayres<- correlog(x = res$UTM_E, y = res$UTM_N, z = as.vector(out_spat$mean$yres), increment=300, resamp=10, quiet = TRUE)

par(mfrow = c(2,1))

plot(bayres$mean.of.class,bayres$correlation, main="No Spatial Autocovariate", 
    xlab="Distance [m]", ylab="Moran's I",ylim=c(-0.1,0.4),xlim=c(0, 5000))
lines(bayres$mean.of.class,bayres$correlation)
bayres$correlation[1:10]
bayres$n[1:10]
bayres$mean.of.class[1:10]

plot(sbayres$mean.of.class,sbayres$correlation, main="Spatial Autocovariate", 
    xlab="Distance [m]", ylab="Moran's I", ylim=c(-0.1,0.4),xlim=c(0, 5000))
lines(sbayres$mean.of.class,sbayres$correlation)
sbayres$correlation[1:10]
sbayres$n[1:10]
sbayres$mean.of.class[1:10]

#######################################################################
############## AUC / Zipkin 2012 adaptations ##########################
#######################################################################
AVG = avgcov$DEJUcov

############################# No Autocov ##############################
#######################################################################
a1=out$sims.list$alpha.psi
b1=out$sims.list$beta.psi.near_roadZ
b3=out$sims.list$beta.psi.near_saltwZ
b4=out$sims.list$beta.psi.RUR_1KMZ
b5=out$sims.list$beta.psi.FOR1_1KMZ
b6=out$sims.list$beta.psi.FOR2_1KMZ
b7=out$sims.list$beta.psi.HRB_1KMZ
b8=out$sims.list$beta.psi.SAV_1KMZ
b9=out$sims.list$beta.psi.SHR_1KMZ
b10=out$sims.list$beta.psi.WET_1KMZ
b11=out$sims.list$beta.psi.URB_100Z
b15=out$sims.list$beta.psi.HRB_100Z
b17=out$sims.list$beta.psi.SHR_100Z
b19=out$sims.list$beta.psi.Near_UrbZ
b20=out$sims.list$beta.psi.Is_sizeZ
b21=out$sims.list$beta.psi.rd_up_100Z
b25=out$sims.list$beta.psi.CR_CLZ

#b1=DEJUout$sims.list$beta.occ                                
ZZ=out$sims.list$z 

site=length(out$mean$z)
iter=length(ZZ[,1])  #iterations/estimates

occ.auto = array(0, dim=c(iter,site))  ##  This time we estimate only for the old sites

 for (i in 1:site) {
  for (z in 1:iter) {
   occ.auto[z,i] = plogis(a1[z] + b1[z] * near_roadZ[i] + b3[z] * near_saltwZ[i] + b4[z] * RUR_1KMZ[i] + b5[z] * FOR1_1KMZ[i] + b6[z] * FOR2_1KMZ[i] + b7[z] * HRB_1KMZ[i] + b8[z] * SAV_1KMZ[i] + b9[z] * SHR_1KMZ[i] + b10[z] * WET_1KMZ[i] + b11[z] * URB_100Z[i] + b15[z] * HRB_100Z[i] + b17[z] * SHR_100Z[i] + b19[z] * Near_UrbZ[i] + b20[z] * Is_sizeZ[i] + b21[z] * rd_up_100Z[i] + b25[z] * CR_CLZ[i]) #a1[z] + b1[z] * AVG[i])
}   }     

auc.auto = rep(0, iter)

  for (z in 1:iter) {
auc.auto[z] = as.numeric(performance(prediction(occ.auto[z,],ZZ[z,]), "auc")@y.values)
}   


summary(auc.auto)
mean(auc.auto)
quantile(auc.auto,prob=0.025)
quantile(auc.auto,prob=0.975)

############################# Autocov #################################
#######################################################################
dist.ind = dist.indx
distance = distance

as1=out_spat$sims.list$alpha.psi
bs1=out_spat$sims.list$beta.psi.near_roadZ
bs3=out_spat$sims.list$beta.psi.near_saltwZ
bs4=out_spat$sims.list$beta.psi.RUR_1KMZ
bs5=out_spat$sims.list$beta.psi.FOR1_1KMZ
bs6=out_spat$sims.list$beta.psi.FOR2_1KMZ
bs7=out_spat$sims.list$beta.psi.HRB_1KMZ
bs8=out_spat$sims.list$beta.psi.SAV_1KMZ
bs9=out_spat$sims.list$beta.psi.SHR_1KMZ
bs10=out_spat$sims.list$beta.psi.WET_1KMZ
bs11=out_spat$sims.list$beta.psi.URB_100Z
bs15=out_spat$sims.list$beta.psi.HRB_100Z
bs17=out_spat$sims.list$beta.psi.SHR_100Z
bs19=out_spat$sims.list$beta.psi.Near_UrbZ
bs20=out_spat$sims.list$beta.psi.Is_sizeZ
bs21=out_spat$sims.list$beta.psi.rd_up_100Z
bs25=out_spat$sims.list$beta.psi.CR_CLZ

bss=out_spat$sims.list$beta.occ.spat
ZZs=out_spat$sims.list$z 

site=length(out_spat$mean$z)
iter=length(ZZs[,1])  #iterations/estimates

occs.auto = array(0, dim=c(iter,site))  ##  This time we estimate only for the old sites
autcov = array(0, dim=c(iter,site))  
#max.rep depends on max. distance file Nbrs; 65 works for up to 2.5km
max.rep = 65
x = array(0, dim=c(iter,site,max.rep)) 
w = array(0, dim=c(iter,site,max.rep)) 

 for (i in 1:site) {
  for (z in 1:iter) {

   x[z,i,1] <- 0
   w[z,i,1] <- 0
  if ((dist.ind[i,3] + dist.ind[i,2]-1) >  dist.ind[i,3]) {
   for (k in dist.ind[i,3]:(dist.ind[i,3] + dist.ind[i,2]-1)) {
      # weigth's for spatial autocorrelation
      x[z,i,(k-dist.ind[i,3]+2)] <- x[z,i,(k-dist.ind[i,3]+1)] + ZZs[z,distance[k,2]]/distance[k,3]
      w[z,i,(k-dist.ind[i,3]+2)] <- w[z,i,(k-dist.ind[i,3]+1)] + 1/distance[k,3]
      } #k
   autcov[z,i] <- x[z,i,(dist.ind[i,2]+1)]/w[z,i,(dist.ind[i,2]+1)]
   }
   else{
      autcov[z,i] <- 0
   }

   occs.auto[z,i] = plogis(as1[z] + bs1[z] * near_roadZ[i] + bs3[z] * near_saltwZ[i] + bs4[z] * RUR_1KMZ[i] + bs5[z] * FOR1_1KMZ[i] + bs6[z] * FOR2_1KMZ[i] + bs7[z] * HRB_1KMZ[i] + bs8[z] * SAV_1KMZ[i] + bs9[z] * SHR_1KMZ[i] + bs10[z] * WET_1KMZ[i] + bs11[z] * URB_100Z[i] + bs15[z] * HRB_100Z[i] + bs17[z] * SHR_100Z[i] + bs19[z] * Near_UrbZ[i] + bs20[z] * Is_sizeZ[i] + bs21[z] * rd_up_100Z[i] + bs25[z] * CR_CLZ[i]+ bss[z] * autcov[z,i])
}   }     



aucs.auto = rep(0, iter)


  for (z in 1:iter) {
aucs.auto[z] = as.numeric(performance(prediction(occs.auto[z,],ZZs[z,]), "auc")@y.values)
} 

summary(aucs.auto)
mean(aucs.auto)
quantile(aucs.auto,prob=0.025)
quantile(aucs.auto,prob=0.975)


#######################################################################
######################## END OF DEJU ##################################
#######################################################################
