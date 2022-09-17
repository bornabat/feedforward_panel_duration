# DGP for panel duration data using mixed proportional hazards and grouped data

databuilder = function(n.person, n.incid, beta = c(-1.05,-.09175,.1835), baseline = 'default',prior=list('gamma',1,1),plots=T,cns=F){
  
  ##########################
  
  n = n.person * n.incid
  if(baseline == 'default'){
    lambda.0 = 0.3*c(0,exp(c(-5,-4.8,-3.9,-3.5,-3,-2.881998, -1.8718475, -1.7865933, -1.7411807, -1.2875784, -1.3771235, -0.3480384, -0.5138887, 0.04541635, 0.04156955, 0.06670596, 0.06700147, 0.1252126, .4717, .354523, .456011,
                             0.9990935, 1.056143,1.45870,1.658721,1.0483721,0.993284, 0.8041152, 0.8244722, 0.1739866, 0.1646879, -0.2309252, -0.2567487, -0.01087248, 0.04445797, -0.08015117, -0.4000353,-1.439640042,-1.9200249,
                             -2.8133245,-3.8643217,-2.169146, -3.128616, -2.088019, -1.302666, -3.142478, -3.050000, -10.559315, -10.125574, -12.546020, -13.589619)),0,0)
    lambda.0 = rep(lambda.0,each=2)
  }
  if(baseline != 'default'){
    lambda.0 = c(baseline)
  }
  
  lengam = length(lambda.0)
  lenbet = length(beta)
  if(prior[[1]]=='gamma'){
    thetas = rgamma(n.person,prior[[2]],prior[[3]])
  }
  if(prior[[1]]=='lognorm'){
    thetas = rlnorm(n.person,prior[[2]],prior[[3]])
  }
  if(prior[[1]]=='discrete'){
    thetas = sample(prior[[2]],n.person,replace = T,prob = prior[[3]])
  }
  X = cbind(1,sample(c(0,0,1,1,1),n,replace = T),runif(n,0.1,2.3))
  
  if(lenbet>3){
    for(j in 1:(lenbet-3)){
      X = cbind(X,ifelse(j%%2==0,sample(c(0,0,1,1,1),n,replace=T),runif(n,0.1,2.3)))
    }
  }
  
  lambdas = (rep(thetas,each=n.incid)*exp(X%*%beta))   %*%   t(as.matrix(lambda.0))
  
  cdfs = cbind(lambdas,0)
  for(i in ncol(cdfs):3){
    cdfs[,i] = rowSums(cdfs[,1:i-1])
  }
  cdfs[,2]=cdfs[,1]
  cdfs=cdfs[,2:ncol(cdfs)]
  
  cdfs = 1- exp(-cdfs)
  
  pdfs = cdfs
  for(i in ncol(pdfs):2){
    pdfs[,i] = pdfs[,i] - pdfs[,i-1]
  } 
  
  y = rep(NA,n)
  for(i in 1:n){
    y[i] = which(rmultinom(1,1,pdfs[i,])==1)
    
  }
  ### generate data
  
  
  
  ##########################
  
  ### relevant plots of the data
  
  ##########################
  
  if(plots){
    
    for(lambda.0.plot in 1:1){
      plot(seq(0,length(lambda.0)-0.1,0.1),rep(lambda.0,each=10),type='l',lwd=2.5,col=adjustcolor('steelblue',.7),xlab='days',ylab='baseline hazard function')
    }  ## baseline hazard plot
    
    for(data.gen in 1:1){
      barplot(table(y))
      lines(lambda.0*max(table(y))/max(lambda.0),lwd = 2, col = adjustcolor('steelblue',.7))
    }       ## generated data plots
    
  }
  ##########################
  
  ### prepare for the correlated and independent HHM
  
  #########################
  
  censor = rep(0, length(y))
  censor[which(y>=(max(y)))]=1
  if(cns){
    censor[which(y>=quantile(y,.60))] = 1
    y[which(y>=quantile(y,.75))] = quantile(y,.75)
  }
  
  Xz = array(NA, dim = c(n.incid,ncol(X),n.person))
  dur = array(0,dim = c(n.person,n.incid))
  cens = array(0,dim = c(n.person,n.incid))
  
  for(i in 1:n.person){
    Xz[,,i]=X[(n.incid*i-(n.incid-1)):(n.incid*i),]
    dur[i,] = y[(n.incid*i-(n.incid-1)):(n.incid*i)]
    cens[i,] = censor[(n.incid*i-(n.incid-1)):(n.incid*i)]
  }   
  
  resp = cbind(y,censor)
  
  return(list('Xz' = Xz, 'X' = X, 'dur'=dur, 'cens'=cens, 'resp' = resp, 'lenbet' = lenbet, 'lengam'=lengam))
  
  
  
  ######################################################################################
}

samp = databuilder(50,10,beta = c(1,1,1))
X = samp$X
Xz = samp$Xz
resp = samp$resp
dur = samp$dur
cens = samp$cens
lenbet = samp$lenbet
lengam = samp$lengam
###

# The function that gives the estimated covariate effects using Bateni (2022) model.

BB22 = function(plotshow=T, mtd=7, iter=5e3, first_0=F){
  
  ######################################################################################
  
  ## longitudinal correlated model [Bateni 2022]:
  
  l.corr.universal = function(Xz,dur,cens,params,magn=1e2){
    
    beta = params[1:lenbet]
    gamma = params[(lenbet+1):(length(params)-2)]
    egamma = exp(gamma)
    nu.0 = exp(params[length(params)-1])
    kappa.0 = exp(params[length(params)])
    
    n.person = nrow(dur)
    n.incid = ncol(dur)
    
    egsum = egamma						
    for(i in 2:length(egamma)){egsum[i] = egsum[i] + egsum[i-1]}
    
    exbt = exp(t(apply(Xz,3,FUN = '%*%',beta)))  ## X is a 3d matrix with j's on the rows, X1,X2,X2 on the cols, and i's on the 3rd D.
    
    edelt = array( egsum[(dur)] , dim=c(n.person,n.incid) )		
    edelt.fwd = array( egsum[(dur)+1] , dim=c(n.person,n.incid) )		
    
    nu = matrix( rep(nu.0,n.person*n.incid) , nrow = n.person )
    kappa = matrix( rep(kappa.0,n.person*n.incid) , nrow = n.person )   
    
    for(i in 2:n.incid){
      nu[,i] = nu[,i-1] + (1-cens[,i-1])		
      
      kappa[,i] = kappa[,(i-1)] + 			
        exbt[,(i-1)]*(edelt[,(i-1)] +				
                        .5*egamma[dur[,(i-1)]+1]*(1-cens[,(i-1)]) )
    }
    
    kaxbt = kappa^(-1) * exbt
    x1 = 1 + kaxbt * edelt			
    x2 = 1 + kaxbt * edelt.fwd			
    
    expbt =  (( x1 )^(-nu) - ( x2 )^(-nu)*(1-cens)) 
    approx = (expbt==0)
    approx = ifelse(is.na(approx),FALSE,approx)
    
    if(sum(approx)){
      
      a = exbt*edelt
      b = exbt*edelt.fwd*(1-cens)
      
      #Laurent series
      exp.approx = (nu*( b - a ) / kappa) + (nu*(nu+1)*(a^2 - b^2))/(2*kappa^2) - (nu*(nu+1)*(nu+2)*(a^3-b^3))/(6*kappa^3)
      
      expbt[approx] = exp.approx[approx]
      
    }
    
    expbt = log(expbt)
    
    box  =  sum(  expbt ,na.rm = T )
    box = ifelse(box==0, sum(expbt), box)
    
    return(-magn*box)
  }
  
  
  l.corr.derivative = function(Xz,dur,cens,params,magn=1e2){
    
    ############################################
    
    
    beta = params[1:lenbet]
    gamma = params[(lenbet+1):(length(params)-2)]
    egamma = exp(gamma)
    nu.0 = exp(params[length(params)-1])
    kappa.0 = exp(params[length(params)])
    
    
    n.person = nrow(dur)
    n.incid = ncol(dur)
    
    egsum = egamma						
    for(i in 2:length(egamma)){egsum[i] = egsum[i] + egsum[i-1]}
    
    exbt = exp(t(apply(Xz,3,FUN = '%*%',beta)))  ## X is a 3d matrix with j's on the rows, X1,X2,X2 on the cols, and i's on the 3rd D.
    
    edelt = array(egsum[(dur)],dim=c(n.person,n.incid))		
    edelt.fwd = array(egsum[(dur)+1],dim=c(n.person,n.incid))		
    
    nu = matrix(rep(nu.0,n.person*n.incid),nrow = n.person)
    kappa = matrix(rep(kappa.0,n.person*n.incid),nrow = n.person)   
    
    for(i in 2:n.incid){
      nu[,i] = nu[,i-1] + (1-cens[,i-1])		
      
      kappa[,i] = kappa[,i-1] + 			
        exbt[,i-1]*(edelt[,i-1] +				
                      .5*egamma[dur[,i-1]+1]*(1-cens[,i-1]) )
    }
    
    kaxbt = kappa^(-1) * exbt
    x1 = 1 + kaxbt * edelt			
    x2 = 1 + kaxbt * edelt.fwd
    
    XZ.new = aperm(Xz,c(3,1,2))
    configs  =  c(n.person,n.incid,dim(Xz)[2])
    
    ############################################
    
    ##d/d.beta
    
    x1.beta =  array(apply(XZ.new,3,FUN='*',(exbt*edelt)),configs)
    x2.beta =  array(apply(XZ.new,3,FUN='*',(exbt*edelt.fwd)),configs)
    
    kap.beta = (exbt*(edelt*(1+cens) + edelt.fwd*(1-cens))/2)
    kap.beta = array(apply(XZ.new,3,FUN='*',kap.beta),configs)
    
    for(i in 2:n.incid){
      kap.beta[,i,] = kap.beta[,i,] + kap.beta[,i-1,]
    }
    kap.beta[,2:ncol(kap.beta),] = kap.beta[,1:(ncol(kap.beta)-1),]
    kap.beta[,1,] = matrix(0,nrow=n.person, ncol = dim(Xz)[2])
    
    ##d/d.delta
    
    conf.delta = c(n.person,n.incid,length(egamma))
    
    x1.delta = array( 0, conf.delta )
    x2.delta = array( 0, conf.delta )
    
    for(i in 1:length(egamma)){
      x1.delta[,,i] = egamma[i]*(dur>=i)
      x2.delta[,,i] = egamma[i]*((dur+1)>=i)
    }
    
    kap.delta = (array(apply(x1.delta,3,FUN='*',1+cens),conf.delta) + array(apply(x2.delta,3,FUN='*',1-cens),conf.delta))/2
    
    x1.delta = array(apply(x1.delta,3,FUN='*',exbt),conf.delta)
    x2.delta = array(apply(x2.delta,3,FUN='*',exbt),conf.delta)
    kap.delta = array(apply(kap.delta,3,FUN='*',exbt),conf.delta)
    
    for(i in 2:n.incid){
      kap.delta[,i,] = kap.delta[,i-1,] + kap.delta[,i,]
    }
    kap.delta[,2:ncol(kap.delta),] = kap.delta[,1:(ncol(kap.delta)-1),]
    kap.delta[,1,] = matrix(0,nrow=n.person, ncol = length(egamma))  
    
    ############################################
    
    ##derivative of the log-likelihood
    
    fg.der = function(g,f,g.p,f.p,nu,configs){       ##Gives the derivative of the matrix (1+(f(x)^-1)g(x))^(-nu) over the vector x. The output is a 3D derivative array.
      
      term1 = array(apply(f.p,3,FUN='*',-g/(f^2)),configs) + array(apply(g.p,3,FUN='*',f^(-1)),configs)
      
      term2 = (-nu) * ( 1 + (f^(-1)) * g ) ^ (-nu-1)
      
      return( list(  'org' = array( apply(term1,3,FUN='*',term2) , configs ), 'ext' = term1  )  )
    }
    
    
    ##beta:
    
    fg.1 = fg.der(exbt*edelt,kappa,x1.beta,kap.beta,nu,configs)
    fg.2 = fg.der(exbt*edelt.fwd,kappa,x2.beta,kap.beta,nu,configs)
    
    dbeta.1 = fg.1$org
    dbeta.2 = array(apply(fg.2$org,3,FUN='*',1-cens),configs)
    
    ##delta:
    
    fg.1.delta = fg.der(exbt*edelt,kappa,x1.delta,kap.delta,nu,conf.delta)
    fg.2.delta = fg.der(exbt*edelt.fwd,kappa,x2.delta,kap.delta,nu,conf.delta)
    
    ddelta.1 = fg.1.delta$org
    ddelta.2 = array(apply(fg.2.delta$org,3,FUN='*',1-cens),conf.delta)
    
    ## kappa.0:
    
    U.kap.1 = x1^(-nu-1) * ( (exbt*edelt*nu) / (kappa^2) )
    U.kap.2 = x2^(-nu-1) * ( (exbt*edelt.fwd*nu) / (kappa^2) ) * (1-cens)
    
    ## nu.0:
    
    U.nu.1 = - x1^(-nu) * log(x1)
    U.nu.2 = - x2^(-nu) * log(x2) * (1-cens)
    
    ##
    
    expbt =  (( x1 )^(-nu) - ( x2 )^(-nu)*(1-cens))
    
    approx = (expbt==0)
    approx = ifelse(is.na(approx),FALSE,approx)
    
    if(sum(approx)){
      a = exbt*edelt
      b = exbt*edelt.fwd*(1-cens)
      
      #Laurent series
      exp.approx = (nu*( b - a ) / kappa) + (nu*(nu+1)*(a^2 - b^2))/(2*kappa^2) - (nu*(nu+1)*(nu+2)*(a^3-b^3))/(6*kappa^3)
      
      expbt[approx] = exp.approx[approx]
      
    }
    
    U.beta = array(apply((dbeta.1),3,FUN='/',expbt),configs) - array(apply((dbeta.2),3,FUN='/',expbt),configs)
    U.delta = array(apply((ddelta.1),3,FUN='/',expbt),conf.delta) - array(apply((ddelta.2),3,FUN='/',expbt),conf.delta)
    U.kap = (U.kap.1/expbt  -  U.kap.2/expbt) *  kappa.0
    U.nu = (U.nu.1/expbt  -  U.nu.2/expbt) *  nu.0
    
    expbt = log(expbt)
    
    
    
    U.beta = apply(U.beta,3,sum)
    U.delta = apply(U.delta,3,sum)
    U.kap = sum(U.kap)
    U.nu = sum(U.nu)
    expbt = sum(expbt)
    
    if(first_0){U.delta[1]=0}
    
    return(list('ell'=-magn*expbt,'Ub' = -magn*U.beta, 'Ud' = -magn*U.delta, 'Un' = -magn*U.nu , 'Uk' = -magn*U.kap))
  }
  
  # Xz is a 3D array of observed explanatory variables. #spells (J) on the rows, #covariates (P) on the columns,
  #and #individuals (n) on the layers (Z axis).
  
  ###############################
  
  ### Optimization
  
  mthds = c('grad-descent','fletcher-reeves','polak-ribiere','adagrad','DFP','BFGS','Adam')
  
  
  #############################################
  
  
  
  for(assumptions in 1:1){
    
    param = c(rep(0,lenbet),rep(0,lengam),0,0)
    
    L = l.corr.derivative(Xz,dur,cens,param)
    ell = L$ell
    Utot  = c(L$Ub,L$Ud,L$Un,L$Uk)
    Utot.old = Utot
    
    Q =   diag(length(param))
    
    dir = -Utot / norm(Utot, type = '2')
    
    alpha = 1
    num=0
    TT = iter
    clpnts = c('orange','red','maroon','purple','navy','steelblue','black')
    
    
  }
  
  if(plotshow){
    plot(c(0,1),col='white',xlim = c(1,TT),y=c(ell/2,ell),xlab='iterations',ylab='value function')
    legend('topright',legend = mthds, col = clpnts,lwd=2)
    grid()
  }
  
  
  #############################################
  
  method = mthds[mtd]
  
  ### 1st order methods
  
  if(method %in% c('grad-descent','fletcher-reeves','polak-ribiere')){
    for(tt in 1:TT){
      
      param.candid = param + dir * alpha
      
      ell.candid = l.corr.universal(Xz, dur, cens, param.candid)
      
      if(ell.candid<=(ell+1e-4 * alpha *(dir%*%Utot))){  #1st Wolfe condition
        
        num=0
        param = param.candid
        L = l.corr.derivative(Xz, dur, cens, param.candid)
        ell = L$ell
        Utot.new = c(L$Ub,L$Ud,L$Un,L$Uk)
        alpha = 1
        
        if(method == 'grad-descent'){
          dir = - Utot.new / norm(Utot.new, type = '2')
        }
        if(method == 'fletcher-reeves'){
          beto = (Utot%*%Utot)/(Utot.old%*%Utot.old)
          #dir = - Utot.new + as.vector(beto) * dir
          dir = -Utot.new/norm(Utot.new, type='2') + as.vector(beto) * dir
        }
        if(method == 'polak-ribiere'){
          beto = (Utot%*%(Utot-Utot.old))/(Utot.old%*%Utot.old)
          # dir = - Utot.new + as.vector(beto) * dir
          dir = -Utot.new/norm(Utot.new, type='2') + as.vector(beto) * dir      
        }
        
        Utot.old = Utot
        Utot = Utot.new
      }
      
      else{
        num = num+1
        alpha = alpha*.9
      }
      if(plotshow){
        points(tt,L$ell,col=clpnts[mtd])
      }
    }
  }
  
  if(method %in% c('adagrad')){
    if(method == 'adagrad'){
      Usum = Utot^2
      eps = 1e-8
      dir = -Utot / (eps + sqrt(Usum))}
    for(tt in 1:TT){
      
      param.candid = param + dir * alpha
      
      ell.candid = l.corr.universal(Xz, dur, cens, param.candid)
      
      if(ell.candid<=(ell+1e-4 * alpha *(dir%*%Utot))){  #1st Wolfe condition
        
        num=0
        param = param.candid
        L = l.corr.derivative(Xz, dur, cens, param.candid)
        ell = L$ell
        Utot.new = c(L$Ub,L$Ud,L$Un,L$Uk)
        alpha = 1
        
        if(method == 'adagrad'){
          Usum = Usum + Utot.new^2
          dir = - Utot.new / (eps + sqrt(Usum))
        }
        
        Utot.old = Utot
        Utot = Utot.new
      }
      
      else{
        num = num+1
        alpha = alpha*.9
      }
      if(plotshow){
        points(tt,L$ell,col=clpnts[mtd])
      }
    }
  }
  
  if(method %in% c('Adam')){
    if(method == 'Adam'){
      b1 = 0.9
      b2 = 0.999
      
      m.0 = (1-b1)* (Utot)
      v.0 = (1-b2)* (Utot^2)
      t.0 = 1
      eps = 1e-8
      
      mhat = m.0/(1-(b1^t.0))
      vhat = v.0/(1-(b2^t.0))
      dir = -mhat / (eps + sqrt(vhat))}
    for(tt in 1:TT){
      
      param.candid = param + dir * alpha
      
      ell.candid = l.corr.universal(Xz, dur, cens, param.candid)
      
      if(ell.candid<=(ell-abs(1e-4 * alpha *as.numeric(dir%*%Utot)))){  #1st Wolfe condition
        
        num=0
        param = param.candid
        
        L = l.corr.derivative(Xz, dur, cens, param.candid)
        
        ell = L$ell
        Utot.new = c(L$Ub,L$Ud,L$Un,L$Uk)
        alpha = 1
        
        if(method == 'Adam'){
          m.0 = b1*m.0 + (1-b1)*Utot.new
          v.0 = b2*v.0 + (1-b2)*(Utot.new^2)
          t.0 = t.0+1
          mhat = m.0/(1-(b1^t.0))
          vhat = v.0/(1-(b2^t.0))
          dir = -mhat / (eps + sqrt(vhat))
        }
        
        Utot.old = Utot
        Utot = Utot.new
      }
      
      else{
        num = num+1
        alpha = alpha*.9
      }
      if(plotshow){
        points(tt,L$ell,col=clpnts[mtd])
      }
    }
  }
  
  
  ### 2nd order methods
  
  if(method %in% c('DFP','BFGS')){
    Q =   diag(length(param))
    for(tt in 1:TT){
      
      param.candid = param + dir * alpha
      
      
      ell.candid = l.corr.universal(Xz, dur, cens, param.candid)
      
      while(is.nan(ell.candid) | ell.candid==0 | is.na(ell.candid)){
        alpha = alpha*.9
        
        param.candid = param + dir * alpha
        
        ell.candid = l.corr.universal(Xz, dur, cens, param.candid)
      }
      
      if(ell.candid<=(ell+1e-4 * alpha *(dir%*%Utot))){  #1st Wolfe condition
        
        num=0
        
        L = l.corr.derivative(Xz, dur, cens, param.candid)
        ell = L$ell
        Utot.new = c(L$Ub,L$Ud,L$Un,L$Uk)
        alpha = 1
        delta = Utot.new - Utot
        zeta = param.candid - param
        
        
        if(method == 'DFP'){
          Q = Q - (Q%*%(delta)%*%(t(delta))%*%Q)/as.numeric(t(delta)%*%Q%*%(delta)) + (zeta%*%t(zeta))/as.numeric(t(zeta)%*%delta)
        }
        if(method == 'BFGS'){
          Q = Q - (zeta%*%t(delta)%*%Q + Q%*%(delta)%*%t(zeta))/as.numeric(t(zeta)%*%delta) + (1 + as.numeric(t(delta)%*%Q%*%delta)/as.numeric(t(zeta)%*%delta))*((zeta%*%t(zeta))/as.numeric(t(zeta)%*%delta))
        }
        
        if(sum(is.na(Q))){
          print(c(zeta, delta))
          print(Utot)
          print(Utot.new)
          break()
        }
        
        dir = as.vector(-Q%*%Utot.new)
        param = param.candid
        Utot.old = Utot
        Utot = Utot.new
      }
      
      else{
        num = num+1
        alpha = alpha*.9
      }
      if(plotshow){
        points(tt,L$ell,col=clpnts[mtd])
      }
    }
  }
  
  
  ############################################
  
  ### printing calculated parameters of interest
  
  #############################################
  
  box = matrix(NA, nrow = (lenbet+2), ncol=1)
  
  box[1:lenbet,1] = param[1:lenbet]
  box[(nrow(box)-1):nrow(box)] = c(exp(param[length(param)-1]),exp(param[length(param)]))
  rownames(box) = c('intercept',paste('beta',1:(lenbet-1)),'nu.0','kappa.0')
  
  knitr::kable(box)
  
  
  
  
  ############################################
}

BB22()