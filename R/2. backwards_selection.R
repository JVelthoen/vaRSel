#' Backwards Variable Selection for Random Forests.
#'
#' @param X A data frame or matrix with the for each covariate a column .
#' @param y A vector of response.
#' @param ntrees The number of trees estimated in one step
#' @param plot boolean indicating if a plot should be made of the entire procedure
#' @param output boolean whether to print progress of the selection
#' @return vars the chosen variables of the optimal model, the prediction error of the chosen model, selection_results of each step
#' @export
backward_selection <- function(X,y,ntrees = 1000,plot=F,output=T){
  if(is.null(colnames(X))) colnames(X) <- paste0("V",1:ncol(X))
  data <- data.frame(X,y)
  back_sel <- data.frame()

  while(ncol(data) > 1){
    back_step <- backwards_step(data,ntrees)
    back_sel <- rbind(back_sel,back_step)
    data <- data[-match(back_step$var_drop,colnames(data))]
    if(output){
      cat(paste("Variable",back_step$var_drop,"dropped.\n"))
    }
  }

  min_pred_error <- min(back_sel$pred_error)
  best_set <- back_sel$var_drop[which(back_sel$pred_error == min_pred_error):nrow(back_sel)] %>% sort()
  if(plot){
    plot_data = back_sel %>% dplyr::mutate(nvars = dplyr::n():1)
    g <- ggplot2::ggplot(plot_data,aes(x= nvars, y=pred_error)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::geom_errorbar(aes(ymin=pred_error - pred_error_sd, ymax=pred_error + pred_error_sd), width=.2) +
      ggplot2::theme_minimal()+
      ggplot2::labs(x="Prediction Error", y = "Number of Variables")
    print(g)
  }
  result <- list(vars = best_set, pred_error = min_pred_error,selection_results = back_sel)
  return(result)
}

#' One Backwards Variable Selection Step.
#'
#' @param data a dataframe with y and X
#' @param ntrees The number of trees estimated in one step
#' @return var_drop the variable with the lowest importance score, pred_error the prediction error of the model, pred_error_sd the prediction error standard deviation
#' @export
backwards_step <- function(data,ntrees){
  step <- 1:20 %>% lapply(function(i){
    fit <- ranger::ranger(y~.,data=data,
                  num.trees = ntrees,
                  importance = "permutation")
    result <- list(importance = fit$variable.importance, error = fit$prediction.error)
    return(result)
  })

  importance_scores <- step %>%
    lapply(function(x){x$importance}) %>%
    do.call('rbind',.) %>%
    colMeans()
  var_drop <- names(importance_scores)[which(importance_scores == min(importance_scores))]
  pred_error <- step %>% lapply(function(x){x$error}) %>% unlist() %>% mean()
  pred_error_sd <- step %>% lapply(function(x){x$error}) %>% unlist() %>% sd()
  return(data.frame(var_drop = var_drop,pred_error=pred_error,pred_error_sd = pred_error_sd,stringsAsFactors = F))
}
