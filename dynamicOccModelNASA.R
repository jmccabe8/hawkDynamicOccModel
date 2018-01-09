############################################################################
####################### define MCMC settings ##############################
library(jagsUI)
ni <- 120000; nt <- 30; nb <- 15000; nc <- 3 #iterations, thinning, burnin, chains

############## Specfy model in bugs language:  #####################
sink( "dynOcc.txt" )
cat( "
     model{
     
     # specify priors
     psi1 ~ dunif( 0, 1 ) #occupancy probability for season 1
     for ( k in 1:( K - 1 ) ){
       beta0.gam[ k ] ~ dnorm( logitmu1, var1 ) # random intercept of seasonal colonization probability
       alpha0.phi[ k ] ~ dnorm( logitmu2, var2 ) # random intercept of seasonal persistence probability
     } #close loop k
     
     for ( i in 1:M ){  #loop over sites
       epsGam[ i ] ~ dnorm( 0, var3) # random site error for colonization
       epsPhi[ i ] ~ dnorm( 0, var4) # random site error for persistence
       } #close loop i
     
     beta1.gam  ~ dnorm( 0, 0.01 ) #colonization - impervious surface 
     alpha1.phi  ~ dnorm( 0, 0.01 ) #persistence - impervious surface 
     beta2.gam ~ dnorm( 0, 0.01 ) #colonization - year 
     alpha2.phi ~ dnorm( 0, 0.01 ) #persistence - year 
     beta3.gam ~ dnorm( 0, 0.01 ) #colonization - interaction impervious * year
     alpha3.phi ~ dnorm( 0, 0.01 ) #persistence - interaction impervious * year
     beta4.gam  ~ dnorm( 0, 0.01 ) #colonization - min temperature 
     alpha4.phi  ~ dnorm( 0, 0.01 ) #persistence - min temperature 
     beta5.gam ~ dnorm( 0, 0.01 ) #colonization - prey availability
     alpha5.phi ~ dnorm( 0, 0.01 ) #persistence - prey availability
     beta6.gam  ~ dnorm( 0, 0.01 ) #colonization - tree canopy cover 
     alpha6.phi  ~ dnorm( 0, 0.01 ) #persistence - tree canopy cover 
     beta7.gam ~ dnorm( 0, 0.01 ) #colonization - tree * year effect
     alpha7.phi ~ dnorm( 0, 0.01 ) #persistence - tree * year effect
     tau1.p ~  dnorm( 0, 0.01 ) #detection - min temperature

     #### Missing values in prey availavility covariate ####
     ### Kery and Royle 2016 p. 175 #################################
     ## Covariate mean as model for missing values #################
     for ( i in 1:M ){ #loop over sites
       for ( k in 1:K ){ #loop over seasons
       prey[ i, k ] ~ dnorm(mu.p_prey, tau.prey) # Assume p_obs normally distributed
       } #close loop i
     } #close loop k

     mu.prey ~ dnorm(0, 0.01) 
     tau.prey <- pow(sd.prey, -2)
     sd.prey ~ dunif(0, 1) 
     ###################################

     logitmu1 <- log( mu1/(1-mu1)) #intercept for colonization 
     mu1 ~ dunif( 0, 1 ) #mean colonization
     var1 <- 1/( sigma1 * sigma1 )
     sigma1 ~ dunif( 0, 3 ) #error colonization (standard deviation)
     
     logitmu2 <- log( mu2/(1-mu2)) #intercept for persistence 
     mu2 ~ dunif( 0, 1 ) #mean persistence
     var2 <- 1/( sigma2 * sigma2 )
     sigma2 ~ dunif( 0, 3 ) #error persistence (standard deviation)
     
     var3 <- 1/( sigma3 * sigma3 )
     sigma3 ~ dunif( 0, 3 ) #error colonization (standard deviation): site
     var4 <- 1/( sigma4 * sigma4 )
     sigma4 ~ dunif( 0, 3 ) #error persistence (standard deviation): site
     
     ### ECOLOGICAL SUBMODEL ##########
     for ( i in 1:M ){  #loop over sites
       z[ i, 1 ] ~ dbern( psi1 ) #latent, true occupancy in season 1
         for ( k in 2:K ) { #loop over seasons
         logit(gam[i, k-1]) <- beta0.gam[k-1] + epsGam[ i ] + beta1.gam * imperv[ i ] + beta2.gam * year[ k ] + beta3.gam * year[ k ] * imperv[ i ] + beta4.gam * minT[ k ] + beta5.gam * prey[ i, k ] + beta6.gam * tree[ i ] + beta7.gam * tree[ i ] * year[ k ]
         logit(phi[i, k-1]) <- alpha0.phi[k-1] + epsPhi[ i ] + alpha1.phi * imperv[ i ] + alpha2.phi * year[ k ] + alpha3.phi * year[ k ] * imperv[ i ] + alpha4.phi * minT[ k ] + alpha5.phi * prey[ i, k ] + alpha6.phi * tree[ i ] + alpha7.phi * tree[ i ] * year[ k ]
         muZ[ i, k ] <- z[ i, k-1 ] * phi[i, k-1 ] + ( 1 - z[ i, k-1 ] ) * gam[i, k-1 ] 
         z[ i, k ] ~ dbern( muZ[ i, k ] ) #latent, true occupancy 
         } #close loop k
     } #close loop i
     
     ####### OBSERVATION SUBMODEL ############## 
     for( n in 1:5 ){ #effort catagories
     # define intercepts for observation model as mean probabilities:
     tau0.p[ n ] <-log(mu.p[ n ] / (1-mu.p[ n ])) 
     mu.p[ n ] ~ dunif( 0, 1 )
     } #close loop n

     for ( i in 1:M ){  #loop over sites
       for( j in 1:J ){  #loop over replicate surveys
         for ( k in 1:K ) { #loop over seasons
         muY[ i, j, k ] <- z[ i, k ] * p[ i, j, k ] 
         y_obs[ i, j, k ] ~ dbern( muY[ i, j, k ] )  
         logit(p[ i, j, k ]) <- tau0.p[eff_mat[ i, j, k ]] + tau1.p * minTobs[i, j, k]
         } #close k loop
       } #close j loop
     } #close i loop
     
     ###### DERIVED PARAMETERS ######
     psi[ 1 ] <- psi1 
     n_occ[ 1 ] <- sum( z[ 1:M, 1 ] ) 
     for ( k in 2:K ){
       psi[ k ] <- psi[ k-1 ] * mean(phi[1:M, k-1]) + ( 1 - psi[ k-1 ] ) * mean(gam[1:M, k-1 ])
       n_occ[ k ] <- sum( z[ 1:M, k ] )
       meanGam[ k-1 ] <- mean(gam[1:M, k-1 ])
       meanPhi[ k-1 ] <- mean(phi[1:M, k-1 ])
       meanP[ k-1 ] <- mean(p[1:M, 1:J, k-1 ])
       } #close loop k
     }  
     
     ", fill = TRUE )

sink()
################ end of model specification  ######

modelname <- "dynOcc.txt"
str( win.data <- list( y_obs = y_obs, M = M, J = J, K = K
                       , eff_mat = eff_mat
                       , minTobs=minTobs
                       , imperv = imperv
                       , tree = tree
                       , minT = minT
                       , prey = prey
                       , year = year
) )
#initial values
zst <- ysum
inits <- function(){ list( z = zst ) }

#parameters monitored
params <- c( "psi", "p","n_occ","meanPhi", "meanGam", "sigma1", "sigma2", "sigma3", "sigma4",
             , "beta0.gam", "aplah0.phi", "beta1.gam", "alpha.phi", "beta2.gam", "alpha2.phi",
             , "beta3.gam", "alpha3.phi", "beta4.gam", "alpha4.phi", "beta5.gam", "alpha5.phi",
             , "beta6.gam", "alpha6.phi", "beta7.gam", "alpha7.phi", "tau0.p", "tau1.p"
             
)

fm <- jags( win.data, inits, params, modelname, n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=nb, parallel=TRUE ) # mean p coeff is alpha0.p (sites/means/seasons)

