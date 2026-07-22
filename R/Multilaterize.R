## Estimating animal location from UHF and VHF telemetry data

Multilaterize<-function(data,RSSIpar,Start=NULL,precision=1/2,
                        max_iter=100000,...){
  UseMethod("Multilaterize")
}


Multilaterize.data.frame<-function(data,RSSIpar,Start=NULL,precision=1/2,
                                   max_iter=100000,...){

  if(!inherits(data,"data.frame")){
    stop(" 'data' must be a data.frame objects")
  }

  if(!inherits(RSSIpar,"data.frame")){
    stop(" 'RSSIpar' must be a data.frame objects")
  }

  if(nrow(data)<3){
    stop(" 'data' must have at least 3 detections")

  }

  required_cols <- c("UTM.x", "UTM.y", "tag_rssi", "tag.ID", "timestamp")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns for 'data': ", paste(missing_cols, collapse = ", "))
  }

  required_cols <- c("RSSI", "Shape", "Lambda")
  missing_cols <- required_cols[!required_cols %in% names(RSSIpar)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns for 'RSSIpar': ", paste(missing_cols, collapse = ", "))
  }


  r<- if (is.null(Start)) {
    c(mean(data$UTM.x, na.rm = TRUE),
      mean(data$UTM.y, na.rm = TRUE))
  } else {
    if (length(Start) != 2) stop("'Start' must be a vector of length 2")
    Start
  }


  ## Take only the RSS value and create a dataframe
  Pred_df<-data.frame(RSSI=data$tag_rssi)

  ## Get shape parameter for each RSSI observation
  Shape<-RSSIpar$Shape[match(Pred_df$RSSI,RSSIpar$RSSI)]

  ## Get lambda parameter for each RSSI observation
  lambda<-RSSIpar$Lambda[match(Pred_df$RSSI,RSSIpar$RSSI)]


  ## Location matrix of nodes
  r_i<-as.matrix(data[,c("UTM.x","UTM.y")])

  ## Small constant to avoid division by zero
  eps<-1e-10

  Converged<-FALSE
  hess_inv <- NULL

  for(t in 1:max_iter){

    ## Deviation from current estimates to node locations
    Diff<-sweep(r_i,2,r)
    Diff<--Diff

    ## Calculate squared distance
    D2<-rowSums(Diff^2) + eps

    ## Number of rows
    n_nodes<-nrow(r_i)

    ## Number of parameters
    n_params<-length(r)


    ## Calculate Gradient
    Grad<-2*colSums((lambda-(Shape/D2))*Diff)


    ## Calculate Hessian matrix
    hess<-matrix(0,n_params,n_params)

    for(i in 1:n_nodes){

      V<-Diff[i,]

      OuterV<-tcrossprod(V)

      T1<-(lambda[i]-(Shape[i]/D2[i]))*diag(n_params)

      T2<-2*(Shape[i]/D2[i]^2)*OuterV


      hess<-hess + 2*(T1 + T2)
    }

    ## Calculate the inverse of the matrix with robust matrix inversion
    hess_inv <- tryCatch(
      solve(hess),
      error = function(e) {
        if (requireNamespace("MASS", quietly = TRUE)) {
          MASS::ginv(hess)
        } else {
          stop("Hessian is singular and MASS package is not available")
        }
      }
    )

    ## The step for updating the location
    Step<-hess_inv%*%Grad

    ## Update location
    r_new<-r-Step


    ## Convergence check
    if(sqrt(sum((r_new-r)^2))<.Machine$double.eps^precision){

      message("Convergence achieved after: ",t, " iterations\n")

      Converged<-TRUE

      r<-r_new

      break
    }

    r<-r_new



  }

  if(!Converged){

    warning("The estimates has not converged after ",
            max_iter,"iterations")
  }
  ## Return results
  return(data.frame(
    tag.local.identifier = data$tag.ID[1],
    timestamp = data$timestamp[1],
    UTM.x= r[1],
    UTM.y= r[2],
    COV.xx=hess_inv[1,1],
    COV.yy=hess_inv[2,2],
    COV.xy=hess_inv[2,1],
    COV.yx=hess_inv[1,2]
  ))

}


