## Some preparation;
logit = function (x) { log(x / (1 - x)) }
expit = function (x) { exp(x) / (1 + exp(x)) }

quan = function(tau, dat, F.est){
  ord = order(dat)
  y = dat[ord]
  cum.p = F.est[ord]
  xi = sapply(tau, function(t){min(y[cum.p >= t], max(y))})
  return(xi)
}

## Main function;
## input: 
## x0 negative group data; x1 positive group data; s is the pre-specified FPR in ROC(s), default s=0.3;
## outputs: 
## point estimates of ROC(s), AUC, Youden index and associated cutoff point along with their 95% confidence intervals;

drm.est = function(x0, x1, s = 0.3) { 
  
  n0 = length(x0)
  n1 = length(x1)
  n = n0 + n1
  x01 = x0[x0 != 0]
  x11 = x1[x1 != 0]
  xij = c(x01, x11) # pooled sample;
  m = length(xij)
  n01 = length(x01)
  n11 =  length(x11)
  rho = n11 / m
  pi0 = 1 - n01 / n0
  pi1 = 1 - n11 / n1
  
  if (pi0 == 0 | pi1 == 0) {
    warning("either pi_0 or pi_1 is zero, use with caution")
  }
  
  # Fit logistic regression to estimate theta;
  label = factor(c(rep(0, n01), rep(1, n11)))
  reg = glm(formula = label~log(xij), family = binomial(link = "logit")) # basis function: qx=log(x);
  
  alpha = unname(reg$coefficients[1]) - log(n11 / n01)
  beta = unname(reg$coefficients[2])
  
  exp.func = function(t) { exp(alpha + beta * log(t)) } 
  
  xij = sort(xij)
  g0 = 1 / ((1 + rho * (exp.func(xij) - 1)) * m)
  g1 = exp.func(xij) * g0
  G0 = cumsum(g0)
  G1 = cumsum(g1)
  G0.fun = function(x){ sum(g0*ifelse(xij <= x, 1, 0)) }
  G1.fun = function(x){ sum(g1*ifelse(xij <= x, 1, 0)) }
  F0 = c(pi0, pi0+ (1 - pi0) * G0)
  F1 = c(pi1, pi1+ (1 - pi1) * G1)
  F0.fun = function(x){ sum(c(pi0, (1 - pi0) * g0) * ifelse(c(0, xij) <= x, 1, 0)) }
  F1.fun = function(x){ sum(c(pi1, (1 - pi1) * g1) * ifelse(c(0, xij) <= x, 1, 0)) }
  
  # Cut-off point and Youden index;
  cutpoint = exp((log((1 - pi0) / (1 - pi1)) - alpha) / beta)
  index = findInterval(cutpoint, c(0, xij)) ####
  dif = F0 - F1
  youden = dif[index]
  
  # If the cut-off point does not lie within the sample range, use alternative estimator;
  if (cutpoint < 0 | cutpoint > max(xij)) {
    youden = max(dif)
    cutpoint = mean(c(0, xij)[dif == youden])
  }
  
  # AUC
  S.posi = sum(G0 * g1)
  auc = pi0 * (1 - 0.5 * pi1) + (1 - pi0) * (1 - pi1) * S.posi
  
  # ROC(s)
  c.quan = quan(1 - s, c(0,xij), F0)
  roc = 1 - F1.fun(c.quan)
  
  est = c(roc, auc, youden, cutpoint)
  names(est) = c('ROC', 'AUC', 'Youden', 'cutpoint')
  
  ## Calculate asymptotic variance;
  ## Some preparation;
  w0 = n0 / n
  w1 = 1 - w0
  delta = w0 * (1 - pi0) + w1 * (1 - pi1)
  h.func = function (t) { return (1 + rho * (exp.func(t) - 1)) }
  
  ## matrix A_theta;
  A.theta0 = delta * rho * (1 - rho)* sum(g0 * (exp.func(xij) / h.func(xij)))
  A.theta1 = delta * rho * (1 - rho)* sum(g0 * (exp.func(xij) * log(xij) / h.func(xij)))
  A.theta2 = delta * rho * (1 - rho)* sum(g0 * (exp.func(xij) * log(xij) * log(xij) / h.func(xij)))
  A.theta = matrix(c(A.theta0, A.theta1, A.theta1, A.theta2), 2, 2)
  A.theta.inv = solve(A.theta)
  
  # Sigma matrix;
  Sigma.matrix = function (t) {
    A0 = sum(g0 * ifelse(xij > t,1,0) / h.func(xij))
    A1 = sum(g0 * ifelse(xij > t,1,0) * exp.func(xij) / h.func(xij))
    A2 = sum(g0 * ifelse(xij > t,1,0) * exp.func(xij)^2 / h.func(xij))
    A3 = c(A1, sum(g0 * ifelse(xij > t,1,0) * exp.func(xij) * log(xij) / h.func(xij)))
    
    sigma.tilde = matrix(0,2,2)
    
    sigma.tilde[1,1] = ((1 - pi0)^2 * A0 - (F0.fun(t) - 1)^2) / delta +
      (1 - G0.fun(t))^2 * pi0 * (1 - pi0) / w0 - rho * (F0.fun(t) - 1)^2 / (delta * (1 - rho))
    
    sigma.tilde[1,2] = ((1 - pi0) * (1 - pi1) * A1 - (F0.fun(t) - 1) * (F1.fun(t) - 1)) / delta -
      (F0.fun(t) - 1) * ((1 - pi1) * (1 - G1.fun(t)) + rho * (F1.fun(t) - 1)) / (delta * (1 - rho))
    
    sigma.tilde[2,1] = sigma.tilde[1,2]
    
    sigma.tilde[2,2] = ((1 - pi1)^2 * A2 - (F1.fun(t) - 1)^2) / delta +
      (1 - G1.fun(t))^2 * pi1 * (1 - pi1) / w1 - ((1 - pi1) * (1 - G1.fun(t)) + rho * (F1.fun(t) - 1))^2 /
      (delta * rho * (1 - rho))
    
    sigma.tilde + c((1 - pi0) * rho, (1 - pi1) * (rho - 1)) %*% t(A3) %*%
      A.theta.inv %*% A3 %*% t(c((1 - pi0) * rho, (1 - pi1) * (rho - 1)))
  }
  
  # Asymptotic variance of ROC(s);
  b.vec = c(- (1 - pi1) * exp.func(c.quan) / (1 - pi0), 1)
  var.roc = t(b.vec) %*% Sigma.matrix(c.quan) %*% b.vec
  
  # Asymptotic variance of AUC;
  u.A = pi0 * (1 - 0.5 * pi1) + (1 - pi0) * (1 - pi1) * (exp.func(xij) * G0 + 1 - G1 - S.posi)
  M1 = sum(g0 * ((1 - pi0) * (1 - pi1) * G0 * exp.func(xij) -
                   rho * exp.func(xij) * u.A / h.func(xij)) )
  M2 = sum(g0 * ((1 - pi0) * (1 - pi1) * G0 * log(xij) * exp.func(xij) -
                   rho * u.A * log(xij) * exp.func(xij) / h.func(xij)) )
  
  # calculate two expectations in the asymptotic variance of AUC;
  auc.e1 = sum(g0 * (u.A^2 / h.func(xij)))
  M = c(M1, M2)
  
  var.auc = (auc.e1 - auc^2) / delta + pi0 * (1 - pi0) * (1 - 0.5 * pi1 - (1 - pi1) * S.posi)^2 / w0 +
    pi1 * (1 - pi1) * (0.5 * pi0 + (1 - pi0) * S.posi)^2 / w1 -
    ((1 - pi0) * (1 - pi1) * S.posi - rho * auc)^2 / (delta * rho * (1 - rho)) +
    t(M) %*% A.theta.inv %*% M
  
  # Asymptotic variance of Youden index;
  Sigma.mat.youden = Sigma.matrix(cutpoint)
  var.youden = Sigma.mat.youden[1,1] - Sigma.mat.youden[1,2] - Sigma.mat.youden[2,1] + Sigma.mat.youden[2,2]
  
  # Asymptotic variance of associated cut-off point;
  var.cutpoint = (pi0 / ((1 - pi0) * w0) + pi1 / ((1 - pi1) * w1) +
                    t(c(1, log(cutpoint))) %*%
                    (A.theta.inv - matrix(c(1 / (rho * (1 - rho) * delta), 0, 0, 0), 2, 2)) %*%
                    c(1, log(cutpoint))) / ((beta/cutpoint)^2)
  
  std.err = sqrt(c(var.roc, var.auc, var.youden, var.cutpoint) / n)
  
  ## CI for ROC(s) and AUC; use logit transformation to improve CI;
  CI.lower = expit(logit(est[1:3]) - 1.96 * std.err[1:3] / (est[1:3] * (1 - est[1:3])))
  CI.upper = expit(logit(est[1:3]) + 1.96 * std.err[1:3] / (est[1:3] * (1 - est[1:3])))  
  
  ## CI for optimal cutoff point \in R+; do not need to use logit-transformed CI;
  CI.lower[4] = cutpoint - 1.96 * std.err[4]
  CI.upper[4] = cutpoint + 1.96 * std.err[4]
  names(CI.lower) = names(CI.upper) = c('ROC', 'AUC', 'Youden', 'cutpoint')
  
  list(est = est, CI.lower = CI.lower, CI.upper = CI.upper)
}


####################################################
## An example for illustration;
####################################################

# Generate semi-continuous data from gamma mixture model;
# a0 = 1; a1 = 1.718  # shape parameters;
# b0 = b1 = 0.25  # rate parameters;
# pi0 = 0.2; pi1 = 0.1
# n0 = n1 = 100
# 
# set.seed(1)
# n00 = rbinom(1, n0, pi0)
# n10 = rbinom(1, n1, pi1)
# x0 = c(rep(0, n00), rgamma(n0 - n00, shape = a0, rate= b0))
# x1 = c(rep(0, n10), rgamma(n1 - n10, shape = a1, rate= b1))
# drm.est(x0, x1, s = 0.3)
