## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment='#>', fig.width=8, fig.height=5,
                     dpi=120, warning=FALSE, message=FALSE)


## ----libraries----------------------------------------------------------------
library(spares)
library(ggplot2)
theme_set(theme_minimal(base_size = 12))


## ----data---------------------------------------------------------------------
set.seed(2026)
n <- 100
p <- 12
groups <- rep(1:4, each = 3)
Sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.35^abs(i-j))
X <- matrix(rnorm(n*p), n, p) %*% chol(Sigma)
colnames(X) <- paste0('x', seq_len(p))
beta_true <- c(1.5, -1, 0.8, 1, -0.7, 0, rep(0, 6))
intercept_true <- 1
sigma <- 0.8
mu <- as.numeric(intercept_true + X %*% beta_true)
y <- mu + rnorm(n, sd = sigma)
knitr::kable(data.frame(variable=colnames(X), group=groups,
                       coefficient=beta_true), digits=2)


## ----grid---------------------------------------------------------------------
lambda <- 10^seq(1, -2.5, length.out = 60)
control <- list(maxit=20000, abstol=1e-10, reltol=1e-10, print.out=FALSE)


## ----scalar-fits--------------------------------------------------------------
fit_lasso <- lmreg(X, y, penalty='lasso', intercept=TRUE,
                   standardize=FALSE, lambda=lambda, control=control)
fit_adalasso <- lmreg(X, y, penalty='adalasso', profile='exponential',
                      exponent=1, intercept=TRUE, standardize=FALSE,
                      lambda=lambda, control=control)
fits <- list('Lasso'=fit_lasso, 'Adaptive Lasso'=fit_adalasso)


## ----helpers------------------------------------------------------------------
path_data <- function(fit, method) {
  data.frame(method=method, lambda=rep(fit$lambda, each=nrow(fit$beta)),
             variable=factor(rep(colnames(X), times=length(fit$lambda)), levels=colnames(X)),
             estimate=as.vector(fit$beta),
             group=factor(rep(groups, times=length(fit$lambda))))
}
coefficient_plot <- function(fits) {
  dat <- do.call(rbind, Map(path_data, fits, names(fits)))
  ggplot(dat, aes(lambda, estimate, colour=variable, group=variable)) +
    geom_hline(yintercept=0, colour='grey80') + geom_line(linewidth=0.6) +
    scale_x_log10() + facet_wrap(~method, ncol=2) +
    labs(x=expression(lambda), y=expression(hat(beta)[j]), colour='Variable') +
    theme(legend.position='bottom', legend.box='vertical') +
    guides(colour=guide_legend(nrow=2, byrow=TRUE, order=1))
}
get_diagnostics <- function(fits) {
  do.call(rbind, lapply(names(fits), function(m) {
    d <- dof(fits[[m]])
    data.frame(method=m, lambda=d$lambda, df=d$df,
               status=d$status, nzero=fits[[m]]$nzero)
  }))
}
diagnostics <- get_diagnostics(fits)
knitr::kable(as.data.frame(table(diagnostics$method, diagnostics$status)),
             col.names=c('Method','Derivative status','Grid points'))
stopifnot(all(diagnostics$status == 'ok'))
knitr::kable(data.frame(method=names(fits),
  nonzero_at_largest_lambda=vapply(fits, function(f) f$nzero[which.max(f$lambda)], numeric(1))), row.names=FALSE)


## ----scalar-paths, fig.cap="Coefficient paths for Lasso and Adaptive Lasso. Colours identify the same variables in both panels; the intercept is omitted."----
coefficient_plot(fits)


## ----slices-------------------------------------------------------------------
coefficient_slices <- function(fits, indices=c(1, 25, 45, 60)) {
  do.call(rbind, lapply(names(fits), function(m) {
    f <- fits[[m]]
    b <- coef(f, s=f$lambda[indices])
    z <- data.frame(variable=rownames(b), truth=c(intercept_true, beta_true), b,
                    check.names=FALSE)
    names(z)[-(1:2)] <- paste0('lambda=', format(signif(f$lambda[indices],3), trim=TRUE))
    cbind(method=m, z)
  }))
}
knitr::kable(coefficient_slices(fits), digits=3, row.names=FALSE)


## ----scalar-dof, fig.cap="Divergence and active-coordinate count, including one unpenalized intercept. Points are computed grid values; connecting lines are visual guides."----
plot_dof <- function(d) {
  z <- rbind(data.frame(d, quantity='Divergence', value=d$df),
             data.frame(d, quantity='1 + nonzero coefficients', value=1+d$nzero))
  ggplot(z, aes(lambda, value, colour=quantity, linetype=quantity)) +
    geom_line(linewidth=0.7) + geom_point(size=0.7) + scale_x_log10() +
    facet_wrap(~method, ncol=2) +
    labs(x=expression(lambda), y='Complexity', colour=NULL, linetype=NULL) +
    theme(legend.position='bottom')
}
plot_dof(diagnostics)


## ----scalar-dof-summary-------------------------------------------------------
knitr::kable(do.call(rbind, lapply(split(diagnostics, diagnostics$method), function(d) {
 data.frame(method=d$method[1], min_correction=min(d$df-1-d$nzero),
            max_correction=max(d$df-1-d$nzero))
})), digits=4, row.names=FALSE)


## ----group-fits---------------------------------------------------------------
fit_glasso <- lmreg(X, y, penalty='glasso', groups=groups,
                    intercept=TRUE, standardize=FALSE,
                    lambda=lambda, control=control)
fit_adaglasso <- lmreg(X, y, penalty='adaglasso', groups=groups,
                       profile='exponential', exponent=1,
                       intercept=TRUE, standardize=FALSE,
                       lambda=lambda, control=control)
group_fits <- list('Group Lasso'=fit_glasso, 'Adaptive Group Lasso'=fit_adaglasso)
group_diagnostics <- get_diagnostics(group_fits)
stopifnot(all(group_diagnostics$status == 'ok'))
knitr::kable(as.data.frame(table(group_diagnostics$method, group_diagnostics$status)),
             col.names=c('Method','Derivative status','Grid points'))


## ----group-paths, fig.cap="Individual coefficient estimates along the group-penalized paths. These are coefficients, not group norms."----
coefficient_plot(group_fits)


## ----group-norms, fig.cap="Estimated group norms along the path. The dashed horizontal levels are the true group norms; colours identify groups consistently."----
group_path <- function(fit, method) {
  do.call(rbind, lapply(unique(groups), function(g) {
    data.frame(method=method, lambda=fit$lambda, group=factor(g),
               norm=sqrt(colSums(fit$beta[groups==g,,drop=FALSE]^2)),
               truth=sqrt(sum(beta_true[groups==g]^2)))
  }))
}
gp <- do.call(rbind, Map(group_path, group_fits, names(group_fits)))
truths <- unique(gp[c('method','group','truth')])
ggplot(gp, aes(lambda, norm, colour=group)) + geom_line(linewidth=0.7) +
  geom_hline(data=truths, aes(yintercept=truth, colour=group), linetype=3) +
  scale_x_log10() + facet_wrap(~method, ncol=2) +
  labs(x=expression(lambda), y=expression(group~norm~~rho[g]), colour='Group') +
  theme(legend.position='bottom')


## ----group-slices-------------------------------------------------------------
knitr::kable(coefficient_slices(group_fits), digits=3, row.names=FALSE)


## ----group-dof, fig.cap="Total divergence, number of nonzero coordinates, active-group count and active-group dimension. Each count includes the intercept contribution."----
group_counts <- do.call(rbind, lapply(names(group_fits), function(m) {
  f <- group_fits[[m]]; d <- dof(f)
  active <- attr(d, 'active_groups')
  ng <- d$active_groups
  dimension <- vapply(active, function(a) sum(table(groups)[a]), numeric(1))
  rbind(data.frame(method=m,lambda=f$lambda,quantity='Divergence',value=d$df),
        data.frame(method=m,lambda=f$lambda,quantity='1 + nonzero coefficients',value=1+f$nzero),
        data.frame(method=m,lambda=f$lambda,quantity='1 + active groups',value=1+ng),
        data.frame(method=m,lambda=f$lambda,quantity='1 + active-group dimension',value=1+dimension))
}))
ggplot(group_counts, aes(lambda,value,colour=quantity,linetype=quantity)) +
  geom_line(linewidth=0.7) + scale_x_log10() + facet_wrap(~method,ncol=2) +
  labs(x=expression(lambda), y='Complexity', colour=NULL,linetype=NULL) +
  theme(legend.position='bottom') +
  guides(colour=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE))


## ----decomposition, fig.cap="Decomposition at the adaptive solutions: total divergence equals the fixed-weight contribution plus the adaptive correction."----
all_fits <- c(fits, group_fits)
adaptive_fits <- all_fits[c('Adaptive Lasso','Adaptive Group Lasso')]
components <- do.call(rbind, lapply(names(adaptive_fits), function(m) {
 d <- dof(adaptive_fits[[m]])
 rbind(data.frame(method=m,lambda=d$lambda,component='Total',value=d$df),
       data.frame(method=m,lambda=d$lambda,component='Fixed-weight contribution',value=d$df_fixed),
       data.frame(method=m,lambda=d$lambda,component='Adaptive correction',value=d$adaptive_correction))
}))
ggplot(components,aes(lambda,value,colour=component,linetype=component)) +
  geom_hline(yintercept=0,colour='grey80') + geom_line(linewidth=0.7) +
  scale_x_log10() + facet_wrap(~method,ncol=2) +
  labs(x=expression(lambda),y='Divergence contribution',colour=NULL,linetype=NULL) +
  theme(legend.position='bottom')


## ----selection----------------------------------------------------------------
set.seed(2027)
foldid <- sample(rep(1:5, length.out=n))
penalties <- c('lasso','adalasso','glasso','adaglasso')
selected <- do.call(rbind, lapply(seq_along(all_fits), function(i) {
 m <- names(all_fits)[i]; f <- all_fits[[i]]
 bic <- select(f, criterion='BIC')
 sure <- select(f, criterion='SURE', sigma2=sigma^2)
 args <- list(X=X,y=y,penalty=penalties[i],intercept=TRUE,
              standardize=FALSE,lambda=lambda,foldid=foldid,control=control)
 if(i>=3) args$groups <- groups
 if(i %in% c(2,4)) {args$profile <- 'exponential'; args$exponent <- 1}
 cv <- do.call(cv.lmreg,args)
 data.frame(method=m,rule=c('BIC','SURE','CV minimum','CV one SE'),
            lambda=c(bic$lambda,sure$lambda,cv$lambda.min,cv$lambda.1se))
}))
knitr::kable(selected,digits=5,row.names=FALSE)


## ----selected-coefficients, fig.height=7, fig.cap="Coefficient paths with vertical lines at each selected regularization value. Coincident choices produce overlapping lines."----
coefficient_plot(all_fits) +
  geom_vline(data=selected,aes(xintercept=lambda,linetype=rule),
             colour='grey25',linewidth=0.55) +
  labs(linetype='Selection rule') +
  guides(linetype=guide_legend(nrow=1,order=2)) +
  theme(legend.box='vertical')


## ----selected-dof, fig.height=6, fig.cap="The same tuning choices located on the divergence paths. These are fixed-lambda divergences evaluated at selected knots, not DoF estimates for the entire selection procedure."----
dall <- get_diagnostics(all_fits)
ggplot(dall,aes(lambda,df)) + geom_line(colour='#b5121b',linewidth=0.8) +
  geom_vline(data=selected,aes(xintercept=lambda,linetype=rule),colour='grey30') +
  scale_x_log10() + facet_wrap(~method,ncol=2) +
  labs(x=expression(lambda),y='Divergence',linetype='Selection rule') +
  theme(legend.position='bottom')


## ----test-data----------------------------------------------------------------
set.seed(2028)
ntest <- 2000
Xtest <- matrix(rnorm(ntest*p),ntest,p) %*% chol(Sigma)
mutest <- as.numeric(intercept_true + Xtest %*% beta_true)
ytest <- mutest + rnorm(ntest,sd=sigma)
performance <- do.call(rbind,lapply(seq_len(nrow(selected)),function(i) {
 s <- selected[i,]; f <- all_fits[[s$method]]
 pred <- as.numeric(predict(f,newx=Xtest,s=s$lambda))
 k <- which.min(abs(f$lambda-s$lambda))
 data.frame(s,nonzero=f$nzero[k],df=dof(f)$df[k],
            mean_MSE=mean((pred-mutest)^2),prediction_MSE=mean((pred-ytest)^2))
}))
knitr::kable(performance,digits=4,row.names=FALSE)


## ----interpretation, results='asis', echo=FALSE-------------------------------
for(m in names(all_fits)) {
 z <- performance[performance$method==m,]
 cat(sprintf('* **%s:** the selected penalties range from %.4f to %.4f, with %d to %d nonzero coefficients. Independent test mean MSE ranges from %.4f to %.4f.\n',
 m,min(z$lambda),max(z$lambda),min(z$nonzero),max(z$nonzero),min(z$mean_MSE),max(z$mean_MSE)))
}


## ----finite-differences-------------------------------------------------------
f <- fit_adaglasso
k <- which.min(abs(f$lambda-0.1))
lambda_check <- f$lambda[k]
refit_mean <- function(response) {
 fit <- lmreg(X,response,penalty='adaglasso',groups=groups,
              profile='exponential',exponent=1,intercept=TRUE,
              standardize=FALSE,lambda=lambda_check,control=control)
 stopifnot(all(as.numeric(fit$convergence)==1))
 as.numeric(predict(fit,newx=X))
}
finite_difference <- function(h) {
 sum(vapply(seq_len(n),function(i) {
   yp <- ym <- y; yp[i] <- yp[i]+h; ym[i] <- ym[i]-h
   (refit_mean(yp)[i]-refit_mean(ym)[i])/(2*h)
 },numeric(1)))
}
h <- c(1e-3,3e-4)
analytic <- dof(f)$df[k]
fd <- vapply(h,finite_difference,numeric(1))
validation <- data.frame(lambda=lambda_check,step=h,analytic=analytic,
                         numerical=fd,absolute_error=abs(fd-analytic))
validation_display <- validation
validation_display$absolute_error <- formatC(validation$absolute_error, format='e', digits=2)
knitr::kable(validation_display,digits=7,row.names=FALSE)
stopifnot(all(validation$absolute_error < 0.01))


## ----session------------------------------------------------------------------
sessionInfo()

