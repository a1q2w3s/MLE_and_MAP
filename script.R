library(ggplot2) #import package

# Import_data --------------------------------------------------------------

setwd('F:/RProjects/221019_HW')
theta_t = as.numeric(as.matrix(read.csv('true_theta.csv'))) #true theta
a_t = as.numeric(as.matrix(read.csv('item_parameter.csv')[1])) #true a parameter
b_t = as.numeric(as.matrix(read.csv('item_parameter.csv')[2])) #true b parameter
answer = as.matrix(read.csv('response_matrix.csv')) #response of subjects


# Initialization ----------------------------------------------------------

n = length(theta_t) #num of subjects
m = length(a_t) #num of items

D=1.702
lambda = 1
episilon = 1e-3

theta_zero = c(1:n)*0
for(i in c(1:n)){
  theta_zero[i]=
    log(
      max(sum(answer[i,]),1)
          /max((m-sum(answer[i,])),1))
}
theta_MLE = theta_zero #to save the estimation by MLE
theta_MAP = theta_zero #to save the estimation by MAP

# FUNC_MLE ---------------------------------------------------------------------

MLE = function(answer,theta,b,a,episilon=1e-3,D=1.702,lambda=1){
  delta = 1
  theta_est = theta
  
  #Iterate following the steps of N-R iteration
  while(abs(delta)>episilon){
    p_est = 1/
      (1+exp(-D*a*(theta_est-b)))
    
    g = sum(D*a*(answer-p_est)) #first derivative of likelihood function
    g_ = -D^2*sum(a^2*p_est*(1-p_est)) #second derivative of likelihood function
    
    if(g_==0){
      break
    }
    
    delta = g/g_
    
    theta_est = theta_est-lambda*delta
    
    #When theta is out of the boundary[-4,4], set it starting value
    if(theta_est>=4 || theta_est<=-4){
      theta_est = theta
      break
    }
  }
  return(theta_est)
}


# FUNC_MAP ---------------------------------------------------------------------

MAP = function(answer,theta,b,a,episilon=1e-3,D=1.702){
  delta = 1
  theta_est = theta
  while(abs(delta)>episilon){
    p_est = 1/
      (1+exp(-D*a*(theta_est-b)))
    
    g = sum(D*a*(answer-p_est))-theta_est
    g_ = -D^2*sum(a^2*p_est*(1-p_est))-1
    
    if(g_==0){
      break
    }
    
    delta = g/g_
    theta_est = theta_est-delta
  }
  return(theta_est)
}


# main --------------------------------------------------------------------

for(i in c(1:n)){
  theta_MLE[i] = MLE(answer[i,],theta_zero[i],b_t,a_t)
  theta_MAP[i] = MAP(answer[i,],theta_zero[i],b_t,a_t)
}


# plot --------------------------------------------------------------------

###
#draw plot
# pic = ggplot(data=NULL,aes(x=theta_t,y=theta_MAP)) +
#   geom_point(alpha=0.5,size=1,color='darkred') +
#   geom_smooth(method = 'lm',se=FALSE, fullrange=TRUE,formula = y~x) +
#   geom_abline(intercept = 0, slope = 1,linetype='longdash')
# print(pic)
# pic + coord_cartesian(xlim=c(-3,3),ylim=c(-3,3))
# 
# ABS_MAP = sum(abs(theta_MAP-theta_t))/n
# RMSE_MAP = sqrt(sum((theta_MAP-theta_t)^2)/n)

###
#Initialize the data-frame for plot
estimation = c(theta_MLE,theta_MAP)
true = c(theta_t,theta_t)
method = c(rep('MLE',2000),rep('MAP',2000))
df = data.frame(estimation,true,method)

#draw plot
pic = ggplot(data = df,aes(x=estimation,y=true,color = method)) +
  geom_point(alpha=0.5,size=1) +
  geom_smooth(method = 'lm',se=FALSE, fullrange=TRUE,formula = y~x)
print(pic)
pic + coord_cartesian(xlim=c(-4,4),ylim=c(-4,4))

#calculate ABS and RMSE
ABS_MLE = sum(abs(theta_MLE-theta_t))/n
RMSE_MLE = sqrt(sum((theta_MLE-theta_t)^2)/n)
ABS_MAP = sum(abs(theta_MAP-theta_t))/n
RMSE_MAP = sqrt(sum((theta_MAP-theta_t)^2)/n)

cat('MLE:','\n','abs=',ABS_MLE,'\n','RMSE=',RMSE_MLE,'\n',
  'MAP:','\n','abs=',ABS_MAP,'\n','RMSE=',RMSE_MAP)


