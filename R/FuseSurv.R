#' FuseSurv: Targeted data fusion for causal survival analysis
#'
#' Estimate treatment-specific survival curves for a target site using
#' (i) target-only, (ii) federated weighting, and (iii) CCOD estimators.
#' Optionally compute additional contrasts (RD, SR, RMST).
#'
#' @param data A data.frame containing all sites.
#' @param covar.name Character vector of covariate column names.
#' @param site.var Column name of site indicator.
#' @param tgt.name Target site label (must match values in data[, site.var]).
#' @param trt.name Column name of treatment indicator (0/1).
#' @param time.var Column name of observed time.
#' @param event Column name of event indicator (1=event, 0=censor).
#' @param fit.times Numeric vector of time points used for model fitting/prediction.
#' @param eval.times Numeric vector of time points for reporting results; must be subset of fit.times.
#' @param prop.SL.library SuperLearner library for propensity models.
#' @param event.SL.library survSuperLearner library for event hazard.
#' @param cens.SL.library survSuperLearner library for censoring hazard.
#' @param n.folds Number of folds for cross-fitting.
#' @param s Random seed.
#' @param return_contrasts Logical; whether to return additional contrasts.
#' @param contrasts Character vector indicating which contrasts to return.
#'   Any subset of c("RD","SR","RMST"). Ignored if return_contrasts=FALSE.
#'
#' @return A list with components:
#' \itemize{
#'   \item curves: list(df.TGT, df.FED, df.CCOD) treatment-specific survival curves
#'   \item weights: (optional) weights matrix used in federated step (if you choose to return)
#'   \item contrasts: (optional) list of requested contrasts for TGT/FED/CCOD
#'   \item internals: (optional) internal objects for debugging/reproducibility
#' }
#' @export
FuseSurv <- function(data,
                     covar.name,
                     site.var,
                     tgt.name,
                     trt.name,
                     time.var,
                     event,
                     fit.times,
                     eval.times,
                     prop.SL.library,
                     event.SL.library,
                     cens.SL.library,
                     n.folds = 5,
                     s = 1,
                     return_contrasts = FALSE,
                     contrasts = c("RD", "SR", "RMST"),
                     return_internals = FALSE) {

  # ---- basic checks (keep light; package can add more) ----
  contrasts <- unique(contrasts)
  allowed_contrasts <- c("RD", "SR", "RMST")
  if (any(!contrasts %in% allowed_contrasts)) {
    stop("contrasts must be subset of c('RD','SR','RMST').")
  }

  # ---- main survival curves ----
  fit0 <- TrtSurvCurve(
    data = data,
    covar.name = covar.name,
    site.var = site.var,
    tgt.name = tgt.name,
    trt.name = trt.name,
    time.var = time.var,
    event = event,
    fit.times = fit.times,
    eval.times = eval.times,
    prop.SL.library = prop.SL.library,
    event.SL.library = event.SL.library,
    cens.SL.library = cens.SL.library,
    n.folds = n.folds,
    s = s
  )

  out <- list(
    curves = list(
      df.TGT  = fit0$df.TGT,
      df.FED  = fit0$df.FED,
      df.CCOD = fit0$df.CCOD
    )
  )

  # ---- optional contrasts ----
  if (isTRUE(return_contrasts)) {

    con <- OtherContrasts(
      site = fit0$site,
      eval.times = fit0$eval.times,
      IF.00 = fit0$IF.00, IF.01 = fit0$IF.01,
      S.00 = fit0$S.00,   S.01 = fit0$S.01,
      Aug.00.mean = fit0$Aug.00.mean, Aug.01.mean = fit0$Aug.01.mean,
      Aug.R0.mean = fit0$Aug.R0.mean, Aug.R1.mean = fit0$Aug.R1.mean,
      Aug.R0.mean.sour = fit0$Aug.R0.mean.sour, Aug.R1.mean.sour = fit0$Aug.R1.mean.sour,
      IF.R0 = fit0$IF.R0, IF.R1 = fit0$IF.R1,
      IF.CCOD.0 = fit0$IF.CCOD.0, IF.CCOD.1 = fit0$IF.CCOD.1,
      ind.R1.ccod = fit0$ind.R1.ccod,
      s = s + 4399
    )

    contrasts_out <- list()

    if ("RD" %in% contrasts) {
      contrasts_out$RD <- list(
        TGT  = con$df.RD.TGT,
        FED  = con$df.RD.FED,
        CCOD = con$df.RD.CCOD
      )
    }

    if ("SR" %in% contrasts) {
      contrasts_out$SR <- list(
        TGT  = con$df.SR.TGT,
        FED  = con$df.SR.FED,
        CCOD = con$df.SR.CCOD
      )
    }

    if ("RMST" %in% contrasts) {
      contrasts_out$RMST <- list(
        TGT  = list(arm0 = con$df.RMST.0.TGT,  arm1 = con$df.RMST.1.TGT,  diff = con$df.RMST.diff.TGT),
        FED  = list(arm0 = con$df.RMST.0.FED,  arm1 = con$df.RMST.1.FED,  diff = con$df.RMST.diff.FED),
        CCOD = list(arm0 = con$df.RMST.0.CCOD, arm1 = con$df.RMST.1.CCOD, diff = con$df.RMST.diff.CCOD)
      )
    }

    out$contrasts <- contrasts_out
  }

  # ---- optional internals ----
  if (isTRUE(return_internals)) {
    out$internals <- fit0
  }

  return(out)
}


TrtSurvCurve <- function(data,
                         covar.name,
                         site.var,
                         tgt.name,
                         trt.name,
                         time.var,
                         event,
                         fit.times,
                         eval.times,
                         prop.SL.library,
                         event.SL.library,
                         cens.SL.library,
                         n.folds=5,
                         s=1) {

  site <- as.character(data[, site.var])
  tgt.name <- as.character(tgt.name)
  site.names <- unique(site)
  non.tgt.site.names <- site.names[site.names!=tgt.name]
  K <- length(site.names)
  site[site==tgt.name] <- "0"
  for(r in 1:(K-1)) {
    site[site==non.tgt.site.names[r]] <- as.character(r)
  }
  site <- as.numeric(site)
  n.site <- table(site)
  prop.site <- n.site/sum(n.site)

  fit.times <- fit.times[fit.times>0]
  N.time <- length(eval.times)

  if(!all(eval.times%in%fit.times)) stop("eval.times must be a subset of fit.times")
  eval.ind <- which(fit.times%in%eval.times)

  data <- data[order(site),]
  site <- site[order(site)]

  set.seed(seed=s)
  seeds <- round(runif(20*K, 0, 20e5))

  ################################################################
  ## ~~~~~~~~~~~~~~~ Target-site-only estimator ~~~~~~~~~~~~~~~ ##
  ################################################################
  dat0 <- data[site==0, ]
  A <- dat0[, trt.name]
  Y <- dat0[, time.var]
  Delta <- dat0[, event]
  X <- dat0[, covar.name]
  n <- length(Y)

  #### data splitting
  set.seed(seeds[1])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  IF.01 <- IF.00 <- S.00 <- S.01 <- NULL
  theta.01 <- theta.00 <- theta.00.sd <- theta.01.sd <- Aug.00.mean <- Aug.01.mean <- rep(0, N.time)

  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]

    #### fit nuisance functions
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,],
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred

    surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0],
                                event=Delta[train.ind][A.train==0],
                                X=X[train.ind,][A.train==0,],
                                new.times=fit.times,
                                event.SL.library=event.SL.library,
                                cens.SL.library=cens.SL.library)

    surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1],
                                event=Delta[train.ind][A.train==1],
                                X=X[train.ind,][A.train==1,],
                                new.times=fit.times,
                                event.SL.library=event.SL.library,
                                cens.SL.library=cens.SL.library)

    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)

    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict

    #### calculate counterfactual survivals
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)

    theta.00 <- theta.00 + S0$surv[eval.ind]
    theta.01 <- theta.01 + S1$surv[eval.ind]
    theta.00.sd <- theta.00.sd + S0$surv.sd[eval.ind]
    theta.01.sd <- theta.01.sd + S1$surv.sd[eval.ind]
    IF.01 <- rbind(IF.01, S1$IF.vals[,eval.ind])
    IF.00 <- rbind(IF.00, S0$IF.vals[,eval.ind])
    S.00 <- rbind(S.00, S.hats.0[,eval.ind])
    S.01 <- rbind(S.01, S.hats.1[,eval.ind])
    Aug.01.mean <- Aug.01.mean + S1$AUG.means[eval.ind]
    Aug.00.mean <- Aug.00.mean + S0$AUG.means[eval.ind]
  }
  Aug.01.mean <- Aug.01.mean / n.folds
  Aug.00.mean <- Aug.00.mean / n.folds
  theta.00 <- theta.00 / n.folds
  theta.01 <- theta.01 / n.folds
  theta.00.sd <- theta.00.sd / (n.folds*sqrt(n.folds))
  theta.01.sd <- theta.01.sd / (n.folds*sqrt(n.folds))

  df.TGT <- data.frame(time=eval.times,
                       surv1=theta.01, surv1.sd=theta.01.sd,
                       surv0=theta.00, surv0.sd=theta.00.sd )

  ### train models from the target site
  surv.fit.0.tgt=survSuperLearner(time=Y[A==0],
                                  event=Delta[A==0],
                                  X=X[A==0,],
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library,
                                  cens.SL.library=cens.SL.library)

  surv.fit.1.tgt=survSuperLearner(time=Y[A==1],
                                  event=Delta[A==1],
                                  X=X[A==1,],
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library,
                                  cens.SL.library=cens.SL.library)

  ##################################################################
  ## ~~~~~~~~~ Density-ratio adjusted local estimates ~~~~~~~~~~~ ##
  ##################################################################
  X0 <- as.matrix(dat0[, covar.name])
  Aug.R0.mean <- Aug.R1.mean <- Aug.R0.mean.sour <- Aug.R1.mean.sour <-
    matrix(0, nrow=N.time, ncol=K-1)
  IF.R0 <- IF.R1 <- df.SOUR <- list()
  for(r in 1:(K-1)) {
    dat.r <- data[site==r, ]
    A <- dat.r[, trt.name]
    Y <- dat.r[, time.var]
    Delta <- dat.r[, event]
    X <- dat.r[, covar.name]
    n <- length(Y)

    #### data splitting
    set.seed(seeds[r+1])
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    IF.R0[[r]] <- IF.R1[[r]] <- matrix(NA, nrow=1, ncol=N.time)
    theta.R1 <- theta.R0 <- theta.R0.sd <- theta.R1.sd <- rep(0, N.time)

    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]
      A.train <- A[train.ind]

      ### fit density ratio and propensity scores
      omega.hats <- estimate_omega_np(x=as.matrix(X)[train.ind,],
                                      x_target=X0,
                                      x.pred=as.matrix(X)[pred.ind,],
                                      method="logistic")

      ps.fit <- SuperLearner(Y=A[train.ind], X=X[train.ind,],
                             family=binomial(), SL.library=prop.SL.library)
      g.hats <- predict(ps.fit, X[pred.ind, ])$pred

      ### predict conditional event survival from the target site model
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0.tgt, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1.tgt, newdata=X[pred.ind,], new.times=fit.times)
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict

      ### for censoring, use the source site model
      surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0],
                                  event=Delta[train.ind][A.train==0],
                                  X=X[train.ind,][A.train==0,],
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library,
                                  cens.SL.library=cens.SL.library)

      surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1],
                                  event=Delta[train.ind][A.train==1],
                                  X=X[train.ind,][A.train==1,],
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library,
                                  cens.SL.library=cens.SL.library)

      surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
      G.hats.0 <- surv.pred.0$cens.SL.predict
      G.hats.1 <- surv.pred.1$cens.SL.predict

      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1,
                         g.hats=g.hats, omega.hats=omega.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0,
                         g.hats=1-g.hats, omega.hats=omega.hats)

      Aug.R1.mean[,r] <- Aug.R1.mean[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean[,r] <- Aug.R0.mean[,r] + S0$AUG.means[eval.ind]
      IF.R1[[r]] <- rbind(IF.R1[[r]], S1$IF.vals[,eval.ind])
      IF.R0[[r]] <- rbind(IF.R0[[r]], S0$IF.vals[,eval.ind])

      ### naive source site augmented term estimates
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict

      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind],
                         fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind],
                         fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)

      Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] + S0$AUG.means[eval.ind]
    }
    theta.R0 <- theta.R0 / n.folds
    theta.R1 <- theta.R1 / n.folds
    theta.R0.sd <- theta.R0.sd / (n.folds*sqrt(n.folds))
    theta.R1.sd <- theta.R1.sd / (n.folds*sqrt(n.folds))

    IF.R1[[r]] <- IF.R1[[r]][-1,]
    IF.R0[[r]] <- IF.R0[[r]][-1,]
    Aug.R1.mean[,r] <- Aug.R1.mean[,r] / n.folds
    Aug.R0.mean[,r] <- Aug.R0.mean[,r] / n.folds
    Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] / n.folds
    Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] / n.folds
  }

  ##################################################################
  ## ~~~~~~~~~~~~~~~ Federated weighting estimate ~~~~~~~~~~~~~~~ ##
  ##################################################################
  set.seed(seeds[K+5])
  wt1 <- wt0 <- chi0 <- chi1 <- augdiff0 <- augdiff1 <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF1.tgt=c(IF.01[,i]-mean(IF.01[,i]), rep(0, length(site[site!=0])))
    IF0.tgt=c(IF.00[,i]-mean(IF.00[,i]), rep(0, length(site[site!=0])))
    IF0.diff <- IF1.diff <- matrix(0, ncol=K-1, nrow=length(IF0.tgt))
    ind0 <- which(site==0)
    for(r in 1:(K-1)) {
      IF0.diff[,r][ind0] <- IF0.tgt[ind0]
      IF1.diff[,r][ind0] <- IF1.tgt[ind0]

      ind <- which(site==r)
      IF0.diff[,r][ind] <- -IF.R0[[r]][,i]
      IF1.diff[,r][ind] <- -IF.R1[[r]][,i]

      chi0[i,r] <- Aug.00.mean[i]-Aug.R0.mean[i,r]
      chi1[i,r] <- Aug.01.mean[i]-Aug.R1.mean[i,r]

      augdiff0[i,r] <- Aug.00.mean[i]-Aug.R0.mean.sour[i,r]
      augdiff1[i,r] <- Aug.01.mean[i]-Aug.R1.mean.sour[i,r]
    }

    cvfit0=try(cv.glmnet(x=IF0.diff, y=IF0.tgt))
    if(class(cvfit0)[1]!="try-error") {
      fit0=try(glmnet(x=IF0.diff, y=IF0.tgt,
                      penalty.factor=chi0[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cvfit0$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") {
        wt0[i,]=coef(fit0, s=cvfit0$lambda.1se)[-1] } else { wt0[i,]=rep(0, K-1) }
    } else { wt0[i,]=rep(0, K-1) }

    cvfit1=try(cv.glmnet(x=IF1.diff, y=IF1.tgt))
    if(class(cvfit1)[1]!="try-error") {
      fit1=try(glmnet(x=IF1.diff, y=IF1.tgt,
                      penalty.factor=chi1[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cvfit1$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit1)[1]!="try-error") {
        wt1[i,]=coef(fit1, s=cvfit1$lambda.1se)[-1] } else { wt1[i,]=rep(0, K-1) }
    } else { wt1[i,]=rep(0, K-1) }
  }
  wt0.tgt <- 1-apply(wt0,1,sum)
  wt1.tgt <- 1-apply(wt1,1,sum)
  weights <- cbind(wt0.tgt, wt0, wt1.tgt, wt1)

  theta0.fed <- apply(augdiff0*wt0, 1, sum) + theta.00
  theta1.fed <- apply(augdiff1*wt1, 1, sum) + theta.01

  all.var0 <- (apply(IF.00,2,var)*(wt0.tgt^2+2*wt0.tgt*(1-wt0.tgt)) + apply(S.00,2,var)*(1-wt0.tgt)^2) / n.site[1]
  all.var1 <- (apply(IF.01,2,var)*(wt1.tgt^2+2*wt1.tgt*(1-wt1.tgt)) + apply(S.01,2,var)*(1-wt1.tgt)^2) / n.site[1]
  for(k in 1:(K-1)) {
    all.var0 <- all.var0 + apply(IF.R0[[k]],2,var)*wt0[,k]^2 / n.site[k+1]
    all.var1 <- all.var1 + apply(IF.R1[[k]],2,var)*wt1[,k]^2 / n.site[k+1]
  }
  theta0.fed.sd <- sqrt(all.var0)
  theta1.fed.sd <- sqrt(all.var1)

  df.FED <- data.frame(time=eval.times,
                       surv1=theta1.fed, surv1.sd=theta1.fed.sd,
                       surv0=theta0.fed, surv0.sd=theta0.fed.sd )

  ##################################################################
  ## ~~~~~~~~~~~~~~~~~~~~~~~ CCOD estimate ~~~~~~~~~~~~~~~~~~~~~~ ##
  ##################################################################
  A <- data[, trt.name]
  Y <- data[, time.var]
  Delta <- data[, event]
  X <- data[, covar.name]
  R <- as.numeric(site==0)
  n <- length(Y)

  #### data splitting
  set.seed(seeds[10*K])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  IF.CCOD.0 <- IF.CCOD.1 <- NULL
  theta.ccod.1 <- theta.ccod.0 <- theta.ccod.0.sd <- theta.ccod.1.sd <- rep(0, N.time)
  ind.R1.ccod <- list()
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    ind.R1.ccod[[i]] <- c(sum(R[pred.ind]==1),length(pred.ind))

    #### fit nuisance functions
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,],
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred

    # propensity score of the target site R=0
    eta0.fit=SuperLearner(Y=R[train.ind], X=X[train.ind,],
                          family=binomial(), SL.library=prop.SL.library)
    eta0.hats=predict(eta0.fit, X[pred.ind, ])$pred

    surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0],
                                event=Delta[train.ind][A.train==0],
                                X=X[train.ind,][A.train==0,],
                                new.times=fit.times,
                                event.SL.library=event.SL.library,
                                cens.SL.library=cens.SL.library)

    surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1],
                                event=Delta[train.ind][A.train==1],
                                X=X[train.ind,][A.train==1,],
                                new.times=fit.times,
                                event.SL.library=event.SL.library,
                                cens.SL.library=cens.SL.library)

    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)

    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict

    #### calculate counterfactual survivals (CCOD)
    S1 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats, eta0.hats=eta0.hats)
    S0 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats, eta0.hats=eta0.hats)

    IF.CCOD.1 <- rbind(IF.CCOD.1, S1$IF.vals[,eval.ind])
    IF.CCOD.0 <- rbind(IF.CCOD.0, S0$IF.vals[,eval.ind])

    theta.ccod.0 <- theta.ccod.0 + S0$surv[eval.ind]
    theta.ccod.1 <- theta.ccod.1 + S1$surv[eval.ind]
    theta.ccod.0.sd <- theta.ccod.0.sd + S0$surv.sd[eval.ind]
    theta.ccod.1.sd <- theta.ccod.1.sd + S1$surv.sd[eval.ind]
  }
  theta.ccod.0 <- theta.ccod.0 / n.folds
  theta.ccod.1 <- theta.ccod.1 / n.folds
  theta.ccod.0.sd <- theta.ccod.0.sd / (n.folds*sqrt(n.folds))
  theta.ccod.1.sd <- theta.ccod.1.sd / (n.folds*sqrt(n.folds))

  df.CCOD <- data.frame(time=eval.times,
                        surv1=theta.ccod.1, surv1.sd=theta.ccod.1.sd,
                        surv0=theta.ccod.0, surv0.sd=theta.ccod.0.sd )

  ind.R1.ccod <- unlist(lapply(seq_along(ind.R1.ccod), function(j) {
    k <- ind.R1.ccod[[j]][1]
    n <- ind.R1.ccod[[j]][2]
    start <- (j - 1) * n + 1
    end   <- (j - 1) * n + k
    start:end
  }))
  return(list(df.TGT=df.TGT, df.FED=df.FED, df.CCOD=df.CCOD,
              IF.00=IF.00, IF.01=IF.01, S.00=S.00, S.01=S.01,
              IF.R0=IF.R0, IF.R1=IF.R1, IF.CCOD.0=IF.CCOD.0, IF.CCOD.1=IF.CCOD.1,
              Aug.00.mean=Aug.00.mean, Aug.01.mean=Aug.01.mean,
              Aug.R0.mean=Aug.R0.mean, Aug.R1.mean=Aug.R1.mean,
              Aug.R0.mean.sour=Aug.R0.mean.sour, Aug.R1.mean.sour=Aug.R1.mean.sour,
              site=site, ind.R1.ccod=ind.R1.ccod, eval.times=eval.times))
}


OtherContrasts <- function(site,
                           eval.times,
                           IF.00, IF.01,
                           S.00, S.01,
                           Aug.00.mean, Aug.01.mean,
                           Aug.R0.mean, Aug.R1.mean,
                           Aug.R0.mean.sour, Aug.R1.mean.sour,
                           IF.R0, IF.R1,
                           IF.CCOD.0, IF.CCOD.1,
                           ind.R1.ccod,
                           s=4399) {

  set.seed(seed=s)

  ##################################################################
  ## ~~~~~~ Target-only and CCOD estimates for RD and RMST ~~~~~~ ##
  ##################################################################

  ### Risk difference
  n0 <- nrow(IF.00)
  n.all <- nrow(IF.CCOD.0)
  prop.R1 <- n0/n.all
  K <- length(unique(site))
  n.site <- table(site)

  IF.TGT.RD <- IF.01-IF.00
  RD.TGT <- apply(IF.TGT.RD, 2, mean)
  RD.TGT.sd <- apply(IF.TGT.RD, 2, sd)/sqrt(n0)
  df.RD.TGT <- data.frame(time=eval.times, RD=RD.TGT, sd=RD.TGT.sd)

  IF.CCOD.RD <- IF.CCOD.1-IF.CCOD.0
  RD.CCOD <- apply(IF.CCOD.RD, 2, mean)
  RD.CCOD.sd <- sqrt((apply(IF.CCOD.RD[ind.R1.ccod,],2,var)*prop.R1 +
                        apply(IF.CCOD.RD[-ind.R1.ccod,],2,var)*(1-prop.R1)) /n.all)
  df.RD.CCOD <- data.frame(time=eval.times, RD=RD.CCOD, sd=RD.CCOD.sd)

  ### Survival ratio
  IF.00.center <- sweep(IF.00, 2, colMeans(IF.00), `-`)
  IF.01.center <- sweep(IF.01, 2, colMeans(IF.01), `-`)
  S0_hat <- colMeans(IF.00)
  S1_hat <- colMeans(IF.01)

  IF.TGT.SR <- sweep(IF.01.center,2,1/S0_hat, `*`)-sweep(IF.00.center,2,(S1_hat/S0_hat^2), `*`)
  SR.TGT <- S1_hat/S0_hat
  SR.TGT.sd <- apply(IF.TGT.SR, 2, sd) / sqrt(n0)
  df.SR.TGT <- data.frame(time=eval.times, SR=SR.TGT, sd=SR.TGT.sd)

  IF.CCOD.0.center <- sweep(IF.CCOD.0, 2, colMeans(IF.CCOD.0), `-`)
  IF.CCOD.1.center <- sweep(IF.CCOD.1, 2, colMeans(IF.CCOD.1), `-`)
  S0.ccod_hat <- colMeans(IF.CCOD.0)
  S1.ccod_hat <- colMeans(IF.CCOD.1)

  IF.CCOD.SR <- sweep(IF.CCOD.1.center,2,1/S0.ccod_hat, `*`) - sweep(IF.CCOD.0.center,2,(S1.ccod_hat/S0.ccod_hat^2), `*`)
  SR.CCOD <- S1.ccod_hat/S0.ccod_hat
  SR.CCOD.sd <- sqrt((apply(IF.CCOD.SR[ind.R1.ccod,],2,var)*prop.R1 +
                        apply(IF.CCOD.SR[-ind.R1.ccod,],2,var)*(1-prop.R1)) /n.all)
  df.SR.CCOD <- data.frame(time=eval.times, SR=SR.CCOD, sd=SR.CCOD.sd)

  ### RMST (by treatment arm and difference)
  dt <- diff(c(0, eval.times))

  IF.TGT.RMST.0 <- IF.00 %*% dt
  RMST.0.TGT <- mean(IF.TGT.RMST.0)
  RMST.0.TGT.sd <- sd(IF.TGT.RMST.0) / sqrt(n0)
  df.RMST.0.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.0.TGT, sd=RMST.0.TGT.sd)

  IF.TGT.RMST.1 <- IF.01 %*% dt
  RMST.1.TGT <- mean(IF.TGT.RMST.1)
  RMST.1.TGT.sd <- sd(IF.TGT.RMST.1) / sqrt(n0)
  df.RMST.1.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.1.TGT, sd=RMST.1.TGT.sd)

  IF.TGT.RMST.diff <- IF.TGT.RD %*% dt
  RMST.diff.TGT <- mean(IF.TGT.RMST.diff)
  RMST.diff.TGT.sd <- sd(IF.TGT.RMST.diff) / sqrt(n0)
  df.RMST.diff.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.diff.TGT, sd=RMST.diff.TGT.sd)

  IF.CCOD.RMST.0 <- IF.CCOD.0 %*% dt
  RMST.0.CCOD <- mean(IF.CCOD.RMST.0)
  RMST.0.CCOD.sd <-  sqrt((var(IF.CCOD.RMST.0[ind.R1.ccod])*prop.R1 +
                             var(IF.CCOD.RMST.0[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.0.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.0.CCOD, sd=RMST.0.CCOD.sd)

  IF.CCOD.RMST.1 <- IF.CCOD.1 %*% dt
  RMST.1.CCOD <- mean(IF.CCOD.RMST.1)
  RMST.1.CCOD.sd <- sqrt((var(IF.CCOD.RMST.1[ind.R1.ccod])*prop.R1 +
                            var(IF.CCOD.RMST.1[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.1.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.1.CCOD, sd=RMST.1.CCOD.sd)

  IF.CCOD.RMST.diff <- IF.CCOD.RD %*% dt
  RMST.diff.CCOD <- mean(IF.CCOD.RMST.diff)
  RMST.diff.CCOD.sd <- sqrt((var(IF.CCOD.RMST.diff[ind.R1.ccod])*prop.R1 +
                               var(IF.CCOD.RMST.diff[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.diff.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.diff.CCOD, sd=RMST.diff.CCOD.sd)

  ##################################################################
  ## ~~~~~~ Federated weighting estimate for SD, SR, RMST ~~~~~~~ ##
  ##################################################################

  ### Survival difference (SD)
  N.time <- length(eval.times)
  Aug.TGT.mean <- Aug.01.mean-Aug.00.mean
  Aug.R.mean <- Aug.R1.mean-Aug.R0.mean
  Aug.R.mean.sour <- Aug.R1.mean.sour-Aug.R0.mean.sour
  wt.RD <- chi.RD <- augdiff.RD <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF.TGT.RD.center <- c(IF.TGT.RD[,i]-mean(IF.TGT.RD[,i]), rep(0, n.all-n0))
    IF.RD.diff <- matrix(0, ncol=K-1, nrow=n.all)
    ind0 <- which(site==0)

    for(r in 1:(K-1)) {
      IF.RD.diff[,r][ind0] <- IF.TGT.RD.center[ind0]
      ind <- which(site==r)
      IF.RD.diff[,r][ind] <- -(IF.R1[[r]][,i]-IF.R0[[r]][,i])
      chi.RD[i,r] <- Aug.TGT.mean[i]-Aug.R.mean[i,r]
      augdiff.RD[i,r] <- Aug.TGT.mean[i]-Aug.R.mean.sour[i,r]
    }
    cv.fit=try(cv.glmnet(x=IF.RD.diff, y=IF.TGT.RD.center))
    if(class(cv.fit)[1]!="try-error") {
      fit0=try(glmnet(x=IF.RD.diff, y=IF.TGT.RD.center,
                      penalty.factor=chi.RD[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cv.fit$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") {
        wt.RD[i,]=coef(fit0, s=cv.fit$lambda.1se)[-1] } else { wt.RD[i,]=rep(0, K-1) }
    } else { wt.RD[i,]=rep(0, K-1) }
  }
  wt.RD.tgt <- 1-apply(wt.RD,1,sum)
  RD.FED <- apply(augdiff.RD*wt.RD, 1, sum) + RD.TGT
  all.var <- (apply(IF.TGT.RD,2,var)*(wt.RD.tgt^2+2*wt.RD.tgt*(1-wt.RD.tgt)) +
                apply(S.01-S.00,2,var)*(1-wt.RD.tgt)^2) / n.site[1]
  for(k in 1:(K-1)) {
    all.var <- all.var + apply(IF.R1[[k]]-IF.R0[[k]], 2, var)*wt.RD[,k]^2 / n.site[k+1]
  }
  RD.FED.sd <- sqrt(all.var)
  df.RD.FED <- data.frame(time=eval.times, RD=RD.FED, sd=RD.FED.sd)

  ### Survival ratio (SR)
  N.time <- length(eval.times)
  Aug.SR.TGT.mean <- (Aug.01.mean/S0_hat) - Aug.00.mean*(S1_hat/S0_hat^2)
  Aug.SR.R.mean <- sweep(Aug.R1.mean,1,1/S0_hat,`*`)-sweep(Aug.R0.mean,1,S1_hat/(S0_hat^2),`*`)
  Aug.SR.R.mean.sour <-sweep(Aug.R1.mean.sour,1,1/S0_hat,`*`)-sweep(Aug.R0.mean.sour,1,S1_hat/(S0_hat^2),`*`)
  wt.SR <- chi.SR <- augdiff.SR <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF.TGT.SR.center <- c(IF.TGT.SR[,i]-mean(IF.TGT.SR[,i]), rep(0, n.all-n0))
    IF.SR.diff <- matrix(0, ncol=K-1, nrow=n.all)
    ind0 <- which(site==0)

    for(r in 1:(K-1)) {
      IF.SR.diff[,r][ind0] <- IF.TGT.SR.center[ind0]
      ind <- which(site==r)
      IF.SR.diff[,r][ind] <- -(IF.R1[[r]][,i]/S0_hat[i]-IF.R0[[r]][,i]*(S1_hat[i]/S0_hat[i]^2))
      chi.SR[i,r] <- Aug.SR.TGT.mean[i]-Aug.SR.R.mean[i,r]
      augdiff.SR[i,r] <- Aug.SR.TGT.mean[i]-Aug.SR.R.mean.sour[i,r]
    }
    cv.fit=try(cv.glmnet(x=IF.SR.diff, y=IF.TGT.SR.center))
    if(class(cv.fit)[1]!="try-error") {
      fit0=try(glmnet(x=IF.SR.diff, y=IF.TGT.SR.center,
                      penalty.factor=chi.SR[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cv.fit$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") {
        wt.SR[i,]=coef(fit0, s=cv.fit$lambda.1se)[-1] } else { wt.SR[i,]=rep(0, K-1) }
    } else { wt.SR[i,]=rep(0, K-1) }
  }
  wt.SR.tgt <- 1-apply(wt.SR,1,sum)
  SR.FED <- apply(augdiff.SR*wt.SR, 1, sum) + SR.TGT
  all.var <- (apply(IF.TGT.SR,2,var)*(wt.SR.tgt^2+2*wt.SR.tgt*(1-wt.SR.tgt)) +
                apply(sweep(Aug.R1.mean,1,1/S0_hat,`*`)-sweep(Aug.R0.mean,1,S1_hat/(S0_hat^2),`*`),2,var)*(1-wt.SR.tgt)^2) / n.site[1]
  for(k in 1:(K-1)) {
    all.var <- all.var + apply(IF.R1[[k]]-IF.R0[[k]], 2, var)*wt.SR[,k]^2 / n.site[k+1]
  }
  SR.FED.sd <- sqrt(all.var)
  df.SR.FED <- data.frame(time=eval.times, SR=SR.FED, sd=SR.FED.sd)


  ### Restrict mean survival times (RMST)
  Aug.RMST.0.TGT.mean <- sum(Aug.00.mean*dt)
  Aug.RMST.1.TGT.mean <- sum(Aug.01.mean*dt)
  Aug.RMST.0.R.mean <- t(t(Aug.R0.mean) %*% dt)
  Aug.RMST.1.R.mean <- t(t(Aug.R1.mean) %*% dt)
  Aug.RMST.0.R.mean.sour <- t(t(Aug.R0.mean.sour) %*% dt)
  Aug.RMST.1.R.mean.sour <- t(t(Aug.R1.mean.sour) %*% dt)
  wt.RMST.0 <- chi.RMST.0 <- augdiff.RMST.0 <- wt.RMST.1 <- chi.RMST.1 <- augdiff.RMST.1 <-
    wt.RMST <- chi.RMST <- augdiff.RMST <- matrix(NA, nrow=1, ncol=K-1)
  IF.TGT.RMST.0.center <- c(IF.TGT.RMST.0-mean(IF.TGT.RMST.0), rep(0, n.all-n0))
  IF.TGT.RMST.1.center <- c(IF.TGT.RMST.1-mean(IF.TGT.RMST.1), rep(0, n.all-n0))
  IF.RMST.0.diff <- IF.RMST.1.diff <- matrix(0, ncol=K-1, nrow=n.all)

  ind0 <- which(site==0)
  for(r in 1:(K-1)) {
    IF.RMST.0.diff[,r][ind0] <- IF.TGT.RMST.0.center[ind0]
    IF.RMST.1.diff[,r][ind0] <- IF.TGT.RMST.1.center[ind0]

    ind <- which(site==r)
    IF.RMST.0.diff[,r][ind] <- -IF.R0[[r]] %*% dt
    IF.RMST.1.diff[,r][ind] <- -IF.R1[[r]] %*% dt

    chi.RMST.0[,r] <- Aug.RMST.0.TGT.mean-Aug.RMST.0.R.mean[,r]
    chi.RMST.1[,r] <- Aug.RMST.1.TGT.mean-Aug.RMST.1.R.mean[,r]
    chi.RMST[,r] <- chi.RMST.1[,r] - chi.RMST.0[,r]

    augdiff.RMST.0[,r] <- Aug.RMST.0.TGT.mean-Aug.RMST.0.R.mean.sour[,r]
    augdiff.RMST.1[,r] <- Aug.RMST.1.TGT.mean-Aug.RMST.1.R.mean.sour[,r]
    augdiff.RMST[,r] <- augdiff.RMST.1[,r] - augdiff.RMST.0[,r]
  }
  IF.TGT.RMST.center <- IF.TGT.RMST.1.center - IF.TGT.RMST.0.center
  IF.RMST.diff <- IF.RMST.1.diff - IF.RMST.0.diff

  cv.fit0=try(cv.glmnet(x=IF.RMST.0.diff, y=IF.TGT.RMST.0.center))
  if(class(cv.fit0)[1]!="try-error") {
    fit0=try(glmnet(x=IF.RMST.0.diff, y=IF.TGT.RMST.0.center,
                    penalty.factor=chi.RMST.0^2,
                    intercept=FALSE,
                    alpha=1,
                    lambda=cv.fit0$lambda.1se,
                    lower.limits=0,
                    upper.limits=1))
    if(class(fit0)[1]!="try-error") {
      wt.RMST.0=coef(fit0, s=cv.fit0$lambda.1se)[-1] } else { wt.RMST.0=rep(0, K-1) }
  } else { wt.RMST.0=rep(0, K-1) }

  cv.fit1=try(cv.glmnet(x=IF.RMST.1.diff, y=IF.TGT.RMST.1.center))
  if(class(cv.fit1)[1]!="try-error") {
    fit1=try(glmnet(x=IF.RMST.1.diff, y=IF.TGT.RMST.1.center,
                    penalty.factor=chi.RMST.1^2,
                    intercept=FALSE,
                    alpha=1,
                    lambda=cv.fit1$lambda.1se,
                    lower.limits=0,
                    upper.limits=1))
    if(class(fit1)[1]!="try-error") {
      wt.RMST.1=coef(fit1, s=cv.fit1$lambda.1se)[-1] } else { wt.RMST.1=rep(0, K-1) }
  } else { wt.RMST.1=rep(0, K-1) }

  cv.fit=try(cv.glmnet(x=IF.RMST.diff, y=IF.TGT.RMST.center))
  if(class(cv.fit)[1]!="try-error") {
    fit=try(glmnet(x=IF.RMST.diff, y=IF.TGT.RMST.center,
                   penalty.factor=chi.RMST^2,
                   intercept=FALSE,
                   alpha=1,
                   lambda=cv.fit$lambda.1se,
                   lower.limits=0,
                   upper.limits=1))
    if(class(fit)[1]!="try-error") {
      wt.RMST=coef(fit, s=cv.fit$lambda.1se)[-1] } else { wt.RMST=rep(0, K-1) }
  } else { wt.RMST=rep(0, K-1) }

  wt.RMST.0.tgt <- 1-sum(wt.RMST.0)
  wt.RMST.1.tgt <- 1-sum(wt.RMST.1)
  wt.RMST.tgt <- 1-sum(wt.RMST)

  RMST.0.FED <- sum(augdiff.RMST.0*wt.RMST.0) + RMST.0.TGT
  RMST.1.FED <- sum(augdiff.RMST.1*wt.RMST.1) + RMST.1.TGT
  RMST.diff.FED <- sum(augdiff.RMST*wt.RMST) + RMST.diff.TGT

  all.var0 <- (var(IF.TGT.RMST.0)*(wt.RMST.0.tgt^2+2*wt.RMST.0.tgt*(1-wt.RMST.0.tgt)) +
                 var(S.00%*%dt)*(1-wt.RMST.0.tgt)^2) / n.site[1]
  all.var1 <- (var(IF.TGT.RMST.1)*(wt.RMST.1.tgt^2+2*wt.RMST.1.tgt*(1-wt.RMST.1.tgt)) +
                 var(S.01%*%dt)*(1-wt.RMST.1.tgt)^2) / n.site[1]
  all.var <- (var(IF.TGT.RMST.diff)*(wt.RMST.tgt^2+2*wt.RMST.tgt*(1-wt.RMST.tgt)) +
                var((S.01-S.00)%*%dt)*(1-wt.RMST.tgt)^2) / n.site[1]
  for(k in 1:(K-1)) {
    all.var0 <- all.var0 + var(IF.R0[[k]]%*%dt)*wt.RMST.0[k]^2 / n.site[k+1]
    all.var1 <- all.var1 + var(IF.R1[[k]]%*%dt)*wt.RMST.1[k]^2 / n.site[k+1]
    all.var <- all.var + var((IF.R1[[k]]-IF.R0[[k]])%*%dt)*wt.RMST[k]^2 / n.site[k+1]
  }
  RMST.0.FED.sd <- sqrt(all.var0)
  RMST.1.FED.sd <- sqrt(all.var1)
  RMST.diff.FED.sd <- sqrt(all.var)

  df.RMST.0.FED <- data.frame(time.max=max(eval.times), RMST=RMST.0.FED, sd=RMST.0.FED.sd)
  df.RMST.1.FED <- data.frame(time.max=max(eval.times), RMST=RMST.1.FED, sd=RMST.1.FED.sd)
  df.RMST.diff.FED <- data.frame(time.max=max(eval.times), RMST=RMST.diff.FED, sd=RMST.diff.FED.sd)

  return(list(df.RD.TGT=df.RD.TGT,
              df.RD.CCOD=df.RD.CCOD,
              df.RD.FED=df.RD.FED,

              df.SR.TGT=df.SR.TGT,
              df.SR.CCOD=df.SR.CCOD,
              df.SR.FED=df.SR.FED,

              df.RMST.0.TGT=df.RMST.0.TGT,
              df.RMST.1.TGT=df.RMST.1.TGT,
              df.RMST.diff.TGT=df.RMST.diff.TGT,

              df.RMST.0.CCOD=df.RMST.0.CCOD,
              df.RMST.1.CCOD=df.RMST.1.CCOD,
              df.RMST.diff.CCOD=df.RMST.diff.CCOD,

              df.RMST.0.FED=df.RMST.0.FED,
              df.RMST.1.FED=df.RMST.1.FED,
              df.RMST.diff.FED=df.RMST.diff.FED))
}


estimate_omega_np = function(x, x_target, x.pred, method="logistic") {
  x.all = rbind(x, x_target)
  z.all = c(rep(1, nrow(x)), rep(0, nrow(x_target)))

  if(method=="logistic") {
    colnames(x.all) <- colnames(x.pred) <- paste0("X", 1:ncol(x.all))
    fit <- glm(z.all~.-z.all, data=data.frame(z.all, x.all), family=binomial(link="logit"))
    src.predict = predict(fit, newdata=data.frame(x.pred), type="response")
  }
  if(method=="glmnet") {
    fit <- cv.glmnet(x=as.matrix(x.all), y=z.all, nfolds=5, family='binomial')
    src.predict = predict(fit, x.pred, type="response", lambda=fit$lambda.min)
  }
  omega = (1-src.predict)/src.predict * nrow(x)/nrow(x_target)
  omega = pmax(pmin(omega, 20), 0.05)
  return(omega)
}


get.survival <- function(Y, Delta, A, R=0,
                         fit.times,
                         S.hats, G.hats, g.hats, omega.hats=1) {

  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]

  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))

  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- AUG.means <- rep(NA, length(fit.times))

  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    augment <- omega.hats*S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    if.func <- S.hats.t0 - augment
    surv[k1] <- mean(if.func)
    # IF.vals[,k1] <- (if.func - mean(if.func))*I(R==0) - augment*I(R!=0)
    IF.vals[,k1] <- if.func*I(R==0) - augment*I(R!=0)
    AUG.means[k1] <- mean(augment)
  }
  surv = pmin(1, pmax(0, surv))
  # surv.sd <- sqrt(colMeans(IF.vals^2, na.rm=T)/n)
  surv.sd <- sqrt(apply(IF.vals, 2, var, na.rm=T)/n)
  return(list(IF.vals=IF.vals, AUG.means=AUG.means, surv=surv, surv.sd=surv.sd))
}


get.survival.CCOD <- function(Y, Delta, A, R,
                              fit.times,
                              S.hats, G.hats, g.hats, eta0.hats) {

  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]

  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))

  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- surv.sd <- rep(NA, length(fit.times))

  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    augment <- S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    if.func <- (S.hats.t0*I(R==1) - eta0.hats*augment) / mean(R==1)
    surv[k1] <- mean(if.func)
    IF.vals[,k1] <- if.func # - mean(if.func)
    surv.sd[k1] <- sqrt(var(if.func[R==1])*mean(R==1) + var(if.func[R!=1])*mean(R!=1)) / sqrt(n)
  }
  surv = pmin(1, pmax(0, surv))
  return(list(IF.vals=IF.vals, surv=surv, surv.sd=surv.sd))
}
