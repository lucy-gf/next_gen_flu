# testing prior distributions



# want a shifted gamma distribution

# first use it as gamma to calculate. 

limits <- c(1.6, 3)
limits_tester <- limits - 1

r0_gamma_pars <- get.gamma.par(p = c(0.025, 0.975), q = limits_tester, 
                               show.output = F, plot = F)
# r0_gamma_pars <- c(11.082591, 9.248767)
# names(r0_gamma_pars) = c('shape', 'rate')

# plot(dgamma(seq(0,4,0.01), shape = r0_gamma_pars[1], rate = r0_gamma_pars[2]), type='l')
# r0_gamma_pars[1]/r0_gamma_pars[2] #mean

# susceptibility, from Bens paper. 
# taking 95% to be lowest and median value in 2+ year olds (0.19, 0.5) <- immunity
# taking mean of other medians as our median = 0.3925
limits <- c(0.5,0.3925,0.19)
limits_sus <- 1 - limits
# fit to a beta
sus_beta_pars <- get.beta.par(p = c(0.025, 0.5, 0.975), 
                              q = limits_sus, 
                              show.output = F, 
                              plot = F)

# sus_beta_pars <- c(50.19925, 32.55043)
# names(sus_beta_pars) <- c('shape1', 'shape2')

#  x <- c(1:1000)
#  y <- rbeta(n=1000, shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2])
#  plot(x, y)
# hist(y, breaks=20)

# plot(dbeta(seq(0,1,0.01), shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2]), type='l')
# sus_beta_pars[1]/(sus_beta_pars[1] + sus_beta_pars[2]) #mean






  
  
  

