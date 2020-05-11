#' Forward Variable Selection for Random Forests.
#'
#' @param X A data frame or matrix with the for each covariate a column .
#' @param y A vector of response.
#' @param nsteps The number of selection steps that need to be done if NULL new variables are added until no more variables are available
#' @param model.vars The variables that the model will start with
#' @param taus The sequence of probability levels, the corresponding quantiles are used to evaluate the performance of the model
#' @param alpha The confidence level for the test default set to 95 percent
#' @param ntrees The number of trees estimated in one step
#' @param sf The subsample fraction that is used in the estimation
#' @param output Indicating wheter every forward step the function prints which variable is selected
#' @return The selected variables in model_vars and the selection results in selection.results.
#' @export
varSel <- function(X,y,
                   nsteps = NULL,model_vars=c(),
                   taus=seq(0.1,0.9,0.1),alpha = 0.95,
                   ntrees = 500,sf=0.5,mtry=10,min_node_size = 1,
                   select_sf = F,output=T,num.threads=2){

  if(select_sf) sf <- select_sample_fraction(X,y)$sf
  if(is.null(nsteps)) nsteps = ncol(X)  # if nsteps is not defined set nsteps to the nr of columns in X
  if(ncol(X) < nsteps){ # if nsteps is larger than ncol(X) set equal to ncol(X) and issue a warning
    nsteps <- ncol(X)
    warning("nsteps was larger than the number of columns of X")
  }
  if(is.null(colnames(X))) colnames(X) <- paste0("V",1:ncol(X)) # if X does not have column names give it column names
  if(!is.null(model_vars)){  # if there are no model_vars
    if(!all(model_vars %in% colnames(X))){ # if the model_vars are not in the column names issue a warning and only use the variables in model_vars
      model_vars <- intersect(model_vars,colnames(X))
      warning("Variables dropped from model_vars that were not contained in X")
    }
  }

  ## Initialize the forward variable selection
  steps_done <- 0
  converged <- F

  ## Do a first step of variable selection
  stepPrev <- varSel_step(X,y,model_vars,taus,ntrees,sf,mtry,min_node_size,num.threads)
  selection_results <- list(stepPrev)

  ## Add steppwise variables to the model until steps_done == nsteps
  while(steps_done < nsteps & !converged){
    ## Set the potential selected variables as the already selected variables plus the best performing from previous step
    potSel <- c(model_vars,stepPrev[1,2])

    ## Take another variable selection step
    step_model <- varSel_step(X,y,potSel,taus,ntrees,sf,mtry,min_node_size,num.threads)
    selection_results[[length(selection_results)+1]] <- step_model

    ## Check whether the previous variable did improve the model
    loss_diff <- merge(stepPrev,step_model,"var") %>%
      dplyr::mutate(x=quantileLoss.x,y=quantileLoss.y) %>%
      dplyr::mutate(diff = x-y) %>% .$diff

    ## Use a sign test to determine if the chosen variable has any influence
    p_value <- binom.test(sum(loss_diff>0),length(loss_diff),alternative="greater")$p.value

    ## Selection/Stoppinf Criterion
    if(p_value<(1-alpha)){
      model_vars <- potSel
      stepPrev <- step_model

    } else{
      converged = T
    }

    steps_done <- steps_done + 1

    if (output & !converged) cat(paste(model_vars[steps_done],"added to the model.\n"))
    if (output & converged) cat("converged \n")
  }

  ## Output the selected model variables and the stepwise selection results
  list(model_vars = model_vars,selection_results = selection_results)
}

#' One foreward step for variable selection
#'
#' @param X A data frame or matrix with the for each covariate a column
#' @param y A vector of response.
#' @param model_vars The variables that the model will start with
#' @param taus The sequence of probability levels, the corresponding quantiles are used to evaluate the performance of the model
#' @param ntrees The number of trees estimated in one step
#' @param sf The subsample fraction that is used in the estimation
#' @param mtry the number of variables to consider at each split
#' @param min_node_size the minimum number of vairables in each leaf node
#' @return Data frame with a row with the calculated quantile loss and a row indicating which variable. The data frame is sorted on the quantileLoss
#' @export
varSel_step <- function(X,y,model_vars,taus,ntrees,sf,mtry,min_node_size,num.threads){
  ## Identify the variables in X that can be added to the model
  var_pot <- setdiff(colnames(X),model_vars)

  ## Fit the stepwise models parallised for each instance of var_par
  result <- parallel::mclapply(var_pot,function(v){
    var_mod <- c(model_vars,v) # The model variables plus one added
    fit <-  grf::quantile_forest(as.matrix(X[,var_mod]),y,
                                 quantiles=taus,
                                 min.node.size=min_node_size,
                                 num.trees=ntrees,
                                 sample.fraction = sf,
                                 mtry=mtry,
                                 num.threads = 1) # fit of the model
    return(quant_error(fit,y,taus)) # calculate the quantile error
  },mc.cores=num.threads) %>%
    unlist() %>% as.data.frame() %>% dplyr::rename("quantileLoss"=".") %>% # set in a data.frame
    dplyr::mutate(var = var_pot) %>% # add variable names
    dplyr::arrange(quantileLoss) # order such that best improving variable is on top
  return(result)
}


#' calculation of the quantile error for random forest fit
#'
#' @param fit A fitted gradient random forest object
#' @param y A vector of observations used for fit
#' @param taus A vector of probability levels used for calculation of the quantile error
#' @return The sum of the quantile error for each of the probability levels from taus calculated from the Out-Of-Bag samples
#' @export
quant_error <- function(fit,y,taus){
  diff_tau = predict(fit,quantiles=taus) %>% as.data.frame() %>%     # Predict with OOB samples the performance
    sapply(function(q){y-q}) %>%                          # calculate differences
    t() %>% as.data.frame() %>% dplyr::mutate(tau = taus)   # add for each quantile the probability level as a variable
  error = diff_tau %>% apply(1,function(x){                                  # calculate the quantile error for each row
      tau <- x[length(x)]
      diff <- x[1:(length(x)-1)]
      mean(diff*(tau - as.numeric(diff<0)))
    }) %>% sum()
  return(error)
}

#' Selection of the optimal sub sample size for the variable selection
#'
#' @param X The data matrix of covariates
#' @param y A vector of observations used for fit
#' @param sf_seq a sequence of subsample fractions to try
#' @param ntrees Number of trees for the variable selection procedure
#' @return The optimal subsample fraction
#' @export
select_sample_fraction <- function(X,y,sf_seq = seq(0.05,0.25,0.05),ntrees= 50){
  result <- sapply(sf_seq, function(s){
    print(s)
    vars <- varSel(X,y,sf=s,ntrees=250)$model_vars
    model <- grf::quantile_forest(as.matrix(X[,vars]),y,mtry=1)
    error <- quant_error(model,y,seq(0.1,0.9,0.1))
    cat(error)
    return(error)
  })
  selection_df = data.frame(sf= sf_seq,error=result)
  return(list(sf=selection_df$sf[selection_df$error == min(selection_df$error)][1],selection_df = selection_df))
}
