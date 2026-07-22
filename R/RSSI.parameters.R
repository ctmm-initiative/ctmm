# Generate gamma distribution parameters (shape and scale) across an RSSI range

RSSI.parameters<- function(model) {
  
  mf<- model.frame(model)
  
  if(ncol(mf) != 2) {
    stop("The model must contain exactly one predictor.")
  }
  
  rssi<- mf[[2]]
  
  if(!is.numeric(rssi)) {
    stop("The predictor must be numeric.")
  }
  
  predictor_name<- names(mf)[2]
  
  param_df<- data.frame(seq(min(rssi, na.rm = TRUE) - 30,
                 max(rssi, na.rm = TRUE) + 30,by = 1))
  
  names(param_df)<- predictor_name
  
  Shape<- fitted(model, newdata=param_df, dpar ="shape")[,"Estimate"]
  Mu<- fitted(model, newdata = param_df)[,"Estimate"]
  
  param_df$Shape<- Shape
  param_df$Lambda<- Mu / Shape
  
  param_df
}
