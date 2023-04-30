# ------------------- Data simulation -------------------

#' CreateSyntheticData
#' function used to simulated data as per the lasso paper
#' According to my notes the seeds used were:
#'   data set 1: 23359 (except it doesn't seem to reproduce the data...)
#'   data set 2: 7
#'   data set 3: 130  (this is confirmed - the lasso example I have online does replicate paper data set 3)
#'   data set 4: 2872945

#'
#' @param whichsim which data set to simulate - 1/2/3/4
#' @param numperiods the number of periods, eg 40
#'
#' @return
#' @export
#'
#' @examples
CreateSyntheticData<-function(whichsim, numperiods)
{
  
  # create the acc/dev/cal parameters
  kk <- rep(1:numperiods, each = numperiods) #AQ
  jj <- rep(1:numperiods, times= numperiods) #DQ
  tt <- kk+jj-1 # PQ
  
  # originla beta
  beta  <- (16/3 - 1)*log(jj)- (1/3)*jj  # original
  # revised
  #beta  <- (16/3 - 1)*log(jj)- (1/2)*jj
  
  # function to generate calendar period effects
  gammafunc <- function(t){
    gg <-
      ifelse( t<=12, gg <- 0.0075*LinearSpline(t,1,12),
              ifelse(t<=24,  gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (t-12)*(t-11)/2,
                     ifelse(t<=32, gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2,
                            ifelse(t<=40, gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2 + 0.002*(t-32)*(t-31)/2,
                                   0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2 + 0.002*(40-32)*(40-31)/2
                            ))))
    1*gg
  }
  
  
  # set alpha/beta/gamma - hard-code up the sim values
  lvl <- 100000
  
  if (whichsim == 1){
    alpha <- log(lvl)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40)
    gamma <- 0
    mu <- exp( alpha + beta + gamma)
  }
  else if (whichsim == 2){
    alpha <- log(lvl)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40)
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma)
  }
  else if (whichsim == 3){
    alpha <- log(lvl)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40)
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma + 0.3*beta*ifelse(kk>16 & jj>20,1,0))
  }
  else if (whichsim == 4){
    alpha <- log(lvl)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40)
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma*((numperiods-1)-LinearSpline(jj,1,numperiods))/(numperiods-1) )  # need to check
  }
  
  varbase <- (0.3 * mu[  kk==1 & jj ==16] )^2 # can scale variance up and down here
  CC  <-  varbase / mu[  kk==1 & jj ==16]
  
  vars   <- CC*mu
  tausq  <- log (vars / (mu^2) + 1)
  
  pmts <- exp( rnorm( numperiods^2, mean = log(mu)-0.5*tausq , sd = sqrt(tausq)  ) )
  
  # indicator for past/future
  past_ind<-(tt<=numperiods)
  
  ### data fram for output
  full<-data.table(pmts, acc=as.integer(kk), dev=as.integer(jj), cal=as.integer(tt), mu, past_ind )
  full
}






#' LinearSpline - spline generator function
#'
#' @param var
#' @param start
#' @param stop
#'
#' @return
#' @export
#'
#' @examples
LinearSpline <- function(var, start, stop){
  pmin(stop - start, pmax(0, var - start))
}



# ------------------- Basis functions -------------------



#' GetScaling - calculate the rho factor scaling for a vector
#'
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
GetScaling <- function(vec) {
  fn <- length(vec)
  fm <- mean(vec)
  fc <- vec - fm
  rho_factor <- ((sum(fc^2))/fn)^0.5
}



#' GetRamps - gets ramps for a vector (both discrete and continuous)
#'
#' @param vec fundamental regressor
#' @param vecname name of regressor
#' @param scaling scaling factor to use
#' @param min_var minimum value of the variable. Required
#' @param max_var maximum value of the variable. Required
#' @param ramp_inc size of ramp increment - not needed if np is defined
#' @param np number of different periods in the data - not needed if ramp_inc is defined
#'
#' @return
#' @export
#'
#' @examples
GetRamps <- function(vec, vecname, scaling, min_var=NULL, max_var=NULL, ramp_inc=NULL, np=NULL){
  
  # add safety checks sometime
  
  # set min and max values if not defined
  if (is.null(min_var)) min_var <- min(vec)
  if (is.null(max_var)) max_var <- max(vec)
  
  
  # calculate np if not defined
  if (is.null(np)) np <- (max_var - min_var) / ramp_inc + 1
  
  # calculate ramp_inc if not defined
  if (is.null(ramp_inc)) ramp_inc <- (max_var - min_var) / (np - 1)
  
  
  # get the number of ramp functions
  nramps <- (np-1)
  
  
  # pre-allocate the matrix to hold the results for speed/efficiency
  n <- length(vec)
  mat <- matrix(data=NA, nrow=n, ncol=nramps)
  cnames <- vector(mode="character", length=nramps)
  
  
  col_indx <- 0
  
  for (i in 1:nramps){
    col_indx <- col_indx + 1
    
    i_start <- (i-1) * ramp_inc + min_var
    
    mat[, col_indx] <- LinearSpline(vec, i_start, 999) / scaling
    cnames[col_indx] <- paste0("L_", i_start, "_999_", vecname)
  }
  
  colnames(mat) <- cnames
  
  return(mat)
}





#' GetInteractions - function to do an arbitrary level of interaction
#'
#' @param vec_list list of vectors to interact
#' @param vecname_list list of the names of these vectors
#' @param scaling_list list of the scaling factors to use
#' @param np_list list of the number of periods to use. Not required if ramp_inc_list is defined
#' @param min_var_list list of the min values of the variables
#' @param max_var_list list of the max values of the variables
#' @param ramp_inc_list list of the ramp increments to use, ie ge this number. Not required if np_list is defined
#'
#' @return
#' @export
#'
#' @examples \dontrun{m <- GetInteractions(vec_list=list(c(1,1,1,1,2,2,2,3,3,4), c(1,2,3,4,1,2,3,1,2,1), c(0.1,0.2,0.3,0.4,0.2,0.3,0.4,0.3,0.4,0.4)),
#'                                         vecname_list=list("acc", "dev", "opt"),
#'                                         scaling_list=list(1,1,2),
#'                                         np_list=list(NULL, NULL, NULL),
#'                                         min_var_list=list(1,1,0),
#'                                         max_var_list=list(4,4,1),
#'                                         ramp_inc_list=list(1,1,0.05) )}
GetInteractions <- function(vec_list, vecname_list, scaling_list, np_list=NULL, min_var_list, max_var_list, ramp_inc_list=NULL ){
  
  # missing lots of safety checks!
  
  # one of np_list and ramp_inc_list must be defined
  # the other one needs to be set within the function
  
  define_np_list <- FALSE
  define_ramp_inc_list <- FALSE
  
  if (is.null(np_list) & is.null(ramp_inc_list)) stop("One of np_list or ramp_inc_list must be defined")
  else if (is.null(np_list)){
    
    define_np_list <- TRUE
    np_list <- vector(mode="list", length=length(ramp_inc_list))
    
  }else if (is.null(ramp_inc_list)) {
    
    define_ramp_inc_list <- TRUE
    ramp_inc_list <- vector(mode="list", length=length(ramp_inc_list))
    
  }
  
  
  # nvars gives the interaction level - usually 2-way but potentially more. NOT tested for higher than 2 though
  nvars <- length(vec_list)
  
  # set values for min, max, np, ramp where none provided
  # also get nramps for each vector
  # finally maintain a count of the total number of interactions
  nramps_list <- vector(mode="list", length=nvars)
  nints <- 1
  
  for (i in 1:nvars){
    if (is.null(min_var_list[[i]])) min_var_list[[i]] <- min(vec_list[[i]])
    if (is.null(max_var_list[[i]])) max_var_list[[i]] <- max(vec_list[[i]])
    
    if (define_np_list) np_list[[i]] <- (max_var_list[[i]] - min_var_list[[i]]) / ramp_inc_list[[i]] + 1
    if (define_ramp_inc_list) ramp_inc_list[[i]] <- (max_var_list[[i]] - min_var_list[[i]]) / (np_list[[i]] - 1)
    
    if (is.null(scaling_list[[i]])) scaling_list[[i]] <- 1
    
    nramps_list[[i]] <- np_list[[i]] - 1
    
    nints <- nints * nramps_list[[i]]
    
  }
  
  # print("number of ramps")
  # print(nramps_list)
  
  
  # preallocate a matrix to hold the interactions
  n <- length(vec_list[[1]])
  mat <- matrix(data=1, nrow=n, ncol=nints)    # data=1 here since will be building interactions step by step
  
  cnames <- vector(mode="character", length=nints)
  
  
  # now calculate the interactions
  # bit of footwork needed to allow an arbitrary number of ints
  
  col_indx <- 1:nints
  
  for (i in 1:nvars){
    
    # print(paste("doing variable", i))
    
    # first get the index vector of where the columns need to be slotted in - use the rep command to generate column indices
    arg_times <- 1
    arg_each <- 1
    
    if (i>1){
      for (j in 1:(i-1)) arg_times <- arg_times * nramps_list[[j]]
    }
    
    if (i<nvars){
      for (j in (i+1):nvars) arg_each <- arg_each * nramps_list[[j]]
    }
    
    vec_col_indx <- rep(1:nramps_list[[i]], each=arg_each, times=arg_times)
    # print("vec_col_indx")
    # print(vec_col_indx)
    
    # make names vector - compile this at the end
    incs <- seq(from=min_var_list[[i]] + ramp_inc_list[[i]], by=ramp_inc_list[[i]], length.out=nramps_list[[i]])
    
    # print(incs)
    
    base_names <- paste("I", vecname_list[[i]], "ge", incs, sep="_")
    names_indx <- rep(base_names, each=arg_each, times=arg_times)
    
    
    # now get the interaction term and slot it into the right place
    for (k in 1:nramps_list[[i]]){
      # print("working out bug")
      # print(vec_list[[i]])
      # print(ramp_inc_list[[i]])
      # print(min_var_list[[i]])
      
      int_vec <- as.integer( vec_list[[i]] >= (k * ramp_inc_list[[i]] + min_var_list[[i]]) ) / scaling_list[[i]]
      # print("int_vec")
      # print(int_vec)
      
      # work out where this goes - depends what level of the hierarchy it is at
      insert_vec <- col_indx[vec_col_indx==k]
      # print("where to insert")
      # print(insert_vec)
      mat[, insert_vec] <- mat[, insert_vec] * int_vec
      
    }   # getting the interaction contribution from a particular variable
    
    # combine names now
    if (i==1) cnames <- names_indx
    else cnames <- paste(cnames, names_indx, sep="_")
    
    # print("cnames")
    # print(cnames)
    
  }   # end of loop over the variables
  
  colnames(mat) <- cnames
  
  # print("cnames at end")
  # print(cnames)
  return(mat)
  
}






#' GenerateVarset
#'
#' @param data = data.frame containing the data
#' @param predictor_object hold the raw predictors and information about these as a list of list. Entries must include
#'                         $num_predictors (integer), $predictors (list of predictor names), $min (list of predictor min values),
#'                         $max (list of predictor max values), $ramp_inc(list of predictor ramp increment to use), $scaling (scaling for the predictor)
#' @param main_effects_list list of the main effects to include in the model
#' @param int_effects_list list of interactions to include in model in form \dontrun{list(c("a", "b), c("a", "b", "c))}
#'
#' @return matrix with all the basis functions including names
#' @export
#'
#' @examples \dontrun{me_main_effects <- c("cal", "opt")
#' me_list_int_effects <- list(c("cal", "opt"), c("acc", "opt"))
#'
#'
#' varset <- GenerateVarset(data = me_data_agg,
#'                          predictor_object = me_predictor_object,
#'                          main_effects_list = me_main_effects,
#'                          int_effects_list = me_list_int_effects)
#' }

GenerateVarset<-function(data,
                         predictor_object=NULL,
                         main_effects_list=NULL,
                         int_effects_list=NULL){
  
  
  # predictor object should contain
  #  $num_predictors
  #  $predictors
  #  $min
  #  $max
  #  $ramp_inc
  #  $scaling
  # It should be a list of lists containing this information - since lists used in args to GetInteractions
  
  
  #----------------------------------------------------------;
  # add main effects
  
  main_mat <- NULL
  
  if (!is.null(main_effects_list)){
    
    for (i in 1:predictor_object$num_predictors){
      pred_name <- predictor_object$predictors[[i]]
      
      if (pred_name %in% main_effects_list){
        
        tem <- GetRamps(vec = data[[pred_name]],
                        vecname = pred_name,
                        scaling = predictor_object$scaling[[i]],
                        min_var = predictor_object$min[[i]],
                        max_var = predictor_object$max[[i]],
                        ramp_inc = predictor_object$ramp_inc[[i]])
        
        main_mat <- if (is.null(main_mat)) tem else cbind(main_mat, tem)
      }
    }
  }
  
  
  #----------------------------------------------------------;
  # add interactions
  
  int_mat <- NULL
  
  if (!is.null(int_effects_list)){
    
    for (i in 1:length(int_effects_list)){
      
      # extract the interaction
      this_int <- int_effects_list[[i]]
      
      int_degree <- length(this_int)
      
      ivec_list      <- vector(mode="list", length=int_degree)
      ivecname_list  <- vector(mode="list", length=int_degree)
      iscaling_list  <- vector(mode="list", length=int_degree)
      imin_var_list  <- vector(mode="list", length=int_degree)
      imax_var_list  <- vector(mode="list", length=int_degree)
      iramp_inc_list <- vector(mode="list", length=int_degree)
      
      
      # extract the predictors in the interaction and the information on them needed for GetInteractions
      for (j in 1:int_degree){
        
        #work out which predictor is the jth member of the interaction - call this k
        for (k in 1:predictor_object$num_predictors){
          if (this_int[[j]] == predictor_object$predictors[[k]]) break
        }
        
        pred_name <- predictor_object$predictors[[k]]
        ivec_list[[j]] <- data[[pred_name]]
        ivecname_list[[j]]  <- pred_name
        iscaling_list[[j]]  <- predictor_object$scaling[[k]]
        imin_var_list[[j]]  <- predictor_object$min[[k]]
        imax_var_list[[j]]  <- predictor_object$max[[k]]
        iramp_inc_list[[j]] <- predictor_object$ramp_inc[[k]]
      }
      
      # make the interactions
      
      tem <- GetInteractions(vec_list = ivec_list,
                             vecname_list = ivecname_list,
                             scaling_list = iscaling_list,
                             min_var_list = imin_var_list,
                             max_var_list = imax_var_list,
                             ramp_inc_list = iramp_inc_list )
      
      
      int_mat <- if (is.null(int_mat)) tem else cbind(int_mat, tem)
      
    }
  }
  
  
  # bind together main and interaction effects
  mat <- if(is.null(main_mat)) int_mat else if(is.null(int_mat)) main_mat else cbind(main_mat, int_mat)
  
  return(mat)
  
  
}


# ------------------- Plotting functions -------------------



#' GraphHeatMap
#'
#' @param dt data.table of input data
#' @param x variable to put across the columns (x-axis)
#' @param y variable to put across the rows (y-axis)
#' @param facet facet variable for generating multiple plots
#' @param actual actual value
#' @param fitted fitted or predicted value
#' @param lims floor and ceilting for actual/fitted values. Default is (25%, 400%)
#' @param xlab columns label (x-axis)
#' @param ylab row labels (y-axis)
#' @param past_line set to TRUE for acc/dev to show diagonal last line of past data. Default=TRUE
#'
#' @return
#' @export
#'
#' @examples
GraphHeatMap <- function(dt, x, y, facet, actual, fitted, lims=c(0.25, 4), xlab, ylab, past_line=TRUE, ncol=2){
  
  
  # copy data to avoid modifying original
  dt <- copy(dt)
  
  # get fails if there is a variable with the same name so make local copies
  local_x <- x
  local_y <- y
  local_actual <- actual
  local_fitted <- fitted
  
  # make restricted Avs F for heatmap and set up past/future split line
  np <- max(dt[[y]])
  
  dt[, .avsf := get(local_actual) / get(local_fitted)]
  dt[, .avsf_restrict_log := log(pmax(min(lims), pmin(max(lims), .avsf)))]
  dt[, .past_line := np + 1 - get(local_y)]
  
  x_var <- "foo"
  y_var <- "bar"
  aes(.data[[x_var]], .data[[y_var]])
  
  
  g <- ggplot(data=dt, aes(x=.data[[local_x]], y=.data[[local_y]])) +
    geom_tile(aes(fill = .avsf_restrict_log))+scale_y_reverse()+
    facet_wrap(~get(facet), ncol=ncol)+
    scale_fill_gradient2(name="AvF_min", low="royalblue", mid="white", high="red3", midpoint=0, space="Lab", na.value="grey50", guide="colourbar")+
    labs(x=xlab, y=ylab)+
    theme(strip.text = element_text(size=8,colour="grey30"), strip.background = element_rect(colour="white", fill="white"))+
    theme(axis.title.x = element_text(size=8), axis.text.x  = element_text(size=7))+
    theme(axis.title.y = element_text(size=8), axis.text.y  = element_text(size=7))+
    theme(element_line(linewidth = 0.25, colour="grey30"))+
    #theme(legend.position=c(0.75,0.1), legend.direction = "horizontal", legend.title=element_blank(), legend.text=element_text(size=8))+
    theme(legend.position="none")+
    NULL
  
  if(past_line) g <- g + geom_line(aes(x=.data[[".past_line"]], y=.data[[local_y]]), colour="grey30", linewidth = 2)
  
  invisible(list(data=dt, graph=g))
  
  
}


# ------------------- Posterior calculations -------------------



#' GetPredAndCoef
#' Returns scaled and filtered predictions and coefficients
#'
#' @param dt 
#' @param glmnet_obj 
#' @param varset 
#' @param yvar name of y variable - need to keep for later calcs
#' @param scaling_factor 
#' @param list_filter 
#'
#' @return
#' @export
#'
#' @examples
GetPredAndCoef <- function(dt, glmnet_obj, varset, yvar, scaling_factor = NULL, list_filter = NULL){
  
  # get all fitted values ------;
  fitted <- predict(glmnet_obj, newx = varset, s = glmnet_obj$lambda, type = "response")
  
  # get all coefficients ------;
  coef <- as.matrix(predict(glmnet_obj, newx = varset, s = glmnet_obj$lambda, type = "coefficients"))
  
  # apply scaling if applicable ------;
  if(!is.null(scaling_factor)){ 
    # fitted values
    fitted <- fitted * scaling_factor
    
    # intercept of coefficients (reqd for priors later)
    coef[1,] <- coef[1,] + log(scaling_factor)
    
    # past data must also be scaled
    dt <- copy(dt)
    dt[!is.na(get(yvar)), (yvar) := get(yvar) * scaling_factor]
  }
  
  # convert to data.table ------;
  coef <- setDT(as.data.frame(coef))
  fitted <- setDT(as.data.frame(fitted))
  fnames <- names(fitted)
  
  # add fitted to dt and make a long file
  dt_wide <- cbind(dt, fitted)
  
  id_vars <- c("acc", "dev", "cal", "past_ind", yvar, "pred")
  #if("boot" %in% names(dt)) id_vars <- c("boot", id_vars)
  
  dt_long <- melt(dt_wide, 
                  id.vars = id_vars, 
                  measure.vars = fnames,
                  value.name = "value",
                  variable.name = "variable",
                  variable.factor = FALSE)
  
  dt_long[, .drop := FALSE]
  
  
  
  # filter results ------;
  if(!is.null(list_filter)){
    
    # summarise observed and .fitted
    fvals <- sort(unique(list_filter$acc$period, list_filter$cal$period))
    
    last_cal <- dt[past_ind == FALSE, min(cal)] - 1   
    last_acc <- dt[past_ind == FALSE, max(acc)]
    
    for(v in fvals){
      
      # calculate stats now  (had this condition in - past_ind == FALSE, but N/A for CAL and not used for ACC in original work and filters calibrated on past vals being in there)
      dt_summ <- rbindlist(
        list(
          acc = dt_wide[acc %between% c(last_acc - v + 1, last_acc), lapply(.SD, sum), .SDcols = c("pred", fnames)],
          cal = dt_wide[cal %between% c(last_cal + 1, last_cal + v), lapply(.SD, sum), .SDcols = c("pred", fnames)]
        ), idcol = "stat_type"
      )
      
      
      # get ratios
      dt_summ[, (fnames) := lapply(.SD, function(x) x / pred), .SDcols = fnames]
      
      # get cutoffs - return TRUE if x falls outside bounds
      for(vv in c("acc", "cal")){
        # handle different periods in acc and cal filtering
        if(!(v %in% list_filter[[vv]]$period)) next
        
        .cutoff <- list_filter[[vv]][period == v, cutoff]
        dt_summ[stat_type == vv, (fnames) := lapply(.SD, function(x) !(x %between% c(1/.cutoff, .cutoff))), .SDcols = fnames]
      } 
      
      # convert to long file to collapse over acc and cal stats
      # add add and if > 0 then .drop = TRUE
      dt_summ_long <- melt(dt_summ, 
                           id.vars = "stat_type",
                           measure.vars = fnames,
                           value.name = ".drop",
                           variable.name = "variable", 
                           variable.factor = FALSE)
      
      dt_summ_long <- dt_summ_long[, .(.drop = sum(.drop)), keyby = .(variable)]
      
      # convert to logicals
      dt_summ_long[, .drop := as.logical(.drop)]
      
      # apply filtering - if new value is FALSE, take that otherwise keep prev (which may be T or F)
      dt_long[ dt_summ_long, on =.(variable), .drop := fifelse(i..drop == TRUE, TRUE, .drop)]
      
    }
  }
  
  # retain only values that pass filtering (all if no filtering)
  dt_long <- dt_long[.drop == FALSE,][, .drop := NULL]
  
  
  # long file for coefficients with matching lambda values ----;
  dt_coef_long <- melt(coef, 
                       measure.vars = fnames,
                       value.name = "coefficient",
                       variable.name = "model_num",
                       variable.factor = FALSE)
  
  # remove the filtered models
  dt_coef_long <- dt_coef_long[model_num %in% dt_long$variable]
  
  
  # finally replace variable with model_num and add in actual lambda value ----;
  # Also rename pred to primary_pred and rename value to pred
  dt_long[, model_num := as.integer(sub("s", "", variable))][, variable := NULL]
  dt_long[, lambda := glmnet_obj$lambda[model_num]]
  setnames(dt_long, c("pred", "value"), c("primary_pred", "pred"))
  
  dt_coef_long[, model_num := as.integer(sub("s", "", model_num))]
  dt_coef_long[, lambda := glmnet_obj$lambda[model_num]]
  
  
  # return the filtered values and coefficients ---->
  return(list(data = dt_long, coefficients = dt_coef_long))
  
}
#debugonce(GetPredAndCoef)
# x <- GetPredAndCoef(list_primary$data, glmnet_obj = list_primary$glmnet_obj, varset = list_primary$varset$all,
#                              scaling_factor = NULL, list_filter = list_dt_filter_params$main, yvar = "pmts")



#' Title
#'
#' @param dt 
#' @param dt_coefficients 
#' @param list_prior_params 
#' @param yvar 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
IntCalculatePosterior <- function(dt, dt_coefficients, list_prior_params, yvar, scale){
  
  dt <- copy(dt)
  dt_coefficients <- copy(dt_coefficients)
  
  # calculate prior -----;
  # add params - note intercept is different to rest
  # could use rowid to do this, dt_coefficients[, laplace_mean := fifelse(rowid(model_num) == 1, list_prior_params$laplace_m_intercept, list_prior_params$laplace_m_parameters)]
  # but will be more efficient to just replicate the vector
  ncoef <- dt_coefficients[model_num == dt_coefficients[1, model_num], .N]
  nmod <- nrow(dt_coefficients[, .N, by = .(model_num)])
  vec_mean <- rep( c(list_prior_params$laplace_m_intercept, rep(list_prior_params$laplace_m_parameters, times = ncoef - 1)), times = nmod)
  vec_scale <- rep( c(list_prior_params$laplace_s_intercept, rep(list_prior_params$laplace_s_parameters, times = ncoef - 1)), times = nmod)
  
  dt_coefficients[, laplace_mean := vec_mean][, laplace_scale := vec_scale]
  dt_coefficients[, log_prior := rmutil::dlaplace(coefficient, laplace_mean, laplace_scale, log=TRUE)]
  
  dt_prior <- dt_coefficients[, .(log_prior = sum(log_prior)), keyby = .(model_num)]  
  dt_prior[, prior_name := list_prior_params$prior_name]
  
  
  # calculate likelihood -----;
  dt[past_ind == TRUE, ll := dgamma(x = get(yvar), shape = 1/scale, scale = pred*scale, log = TRUE)]
  dt_ll <- dt[past_ind == TRUE, .(log_likelihood = sum(ll)), keyby = .(model_num)]
  
  
  # add process error and sum reserves -----;
  g_shape <- 1 / scale  
  dt[past_ind == FALSE, g_scale := pred / g_shape]
  dt[past_ind == FALSE, pred_pe := rgamma(.N, shape = g_shape, scale = g_scale)]
  dt_res <- dt[past_ind == FALSE, .(reserve = sum(pred), reserve_pe = sum(pred_pe)), keyby = .(model_num)]
  
  # combine results into a data.table ----;
  dt_res[ dt_ll, on = .(model_num), log_likelihood := i.log_likelihood]
  dt_res[dt_prior, on = .(model_num), `:=`(log_prior = i.log_prior, prior_name = i.prior_name)]
  
  # calculate the posterior ------;
  # need to do some scaling to avoid underflow errors
  dt_res[, log_unscaled_posterior := log_prior + log_likelihood]
  # scale now
  dt_res[, log_unscaled_posterior := log_unscaled_posterior - max(dt_res$log_unscaled_posterior)]
  # exponentiate and scale
  dt_res[, unscaled_posterior := exp(log_unscaled_posterior)]
  dt_res[, posterior := unscaled_posterior / sum(dt_res$unscaled_posterior)]
  
  # return -----;
  return(dt_res)
}

# list_pp <- list(laplace_m_intercept = 19.39,
#                 laplace_s_intercept = 2,
#                 laplace_m_parameters = 0,
#                 laplace_s_parameters = 0.06,
#                 prior_name = "lambda.1se")
# 
# debugonce(IntCalculatePosterior)
# x_pos <- IntCalculatePosterior(x$data, x$coefficients, list_prior_params = list_pp, yvar = "pmts", scale = list_primary$scale)


#' Title
#'
#' @param list_prior_params 
#' @param prior_name 
#'
#' @return
#' @export
#'
#' @examples
GetPriorParams <- function(list_prior_params, prior_name){
  
  pnum <- which(list_prior_params$prior_name == prior_name)
  
  return( lapply(list_prior_params, function(x){ if(length(x) == 1) return(x) else return(x[pnum]) }) )
  
}



#' Title
#'
#' @param dt 
#' @param glmnet_obj 
#' @param varset 
#' @param yvar 
#' @param scaling_factor 
#' @param list_filter 
#' @param list_prior_params 
#' @param prior_name 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
CalculatePosterior <- function(dt, 
                               glmnet_obj, 
                               varset, 
                               yvar, 
                               dt_scaling_factor = NULL, 
                               list_filter = NULL, 
                               list_prior_params, 
                               prior_name = NULL,
                               scale,
                               return_fitted = FALSE){
  
  # extract prior parameters if name given -----;
  # otherwise we assume list is ready to go
  if(!is.null(prior_name)){
    list_prior_params <- GetPriorParams(list_prior_params, prior_name)
    .prior_name <- prior_name
    scaling_factor <- if(!is.null(dt_scaling_factor)) dt_scaling_factor[prior_name == .prior_name, scaling_factor] else NULL
  } else scaling_factor <- if(!is.null(dt_scaling_factor)) dt_scaling_factor$scaling_factor else NULL
  
  
  # get fitted values, coefficients, filter and scale the values (including past values) if required ----;;
  list_fits <- GetPredAndCoef(dt, glmnet_obj, varset, yvar, scaling_factor = scaling_factor, list_filter = list_filter)
  
  # trap cases where all models removed by filtering ----;
  if(nrow(list_fits$data) == 0){
    list_ret <- list(probs = NULL, fitted = NULL)
    return(list_ret)
  }
  
  # calculate posterior probabilities ----;
  
  dt_post <- IntCalculatePosterior(list_fits$data, list_fits$coefficients, list_prior_params = list_prior_params, yvar = yvar, scale = scale)
  
  # return the posterior probs and future fitted if required
  list_ret <- list(probs = dt_post)
  if(return_fitted){
    list_ret$fitted <- list_fits$data[past_ind == FALSE]
    if(!is.null(prior_name)) list_ret$fitted[, .prior_name := prior_name]
    setnames(list_ret$fitted, ".prior_name", "prior_name")
  } 
  
  return(list_ret)
}  



#' Title
#'
#' @param dt 
#'
#' @return
#' @export
#'
#' @examples
SummariseBayesianAveraging <- function(dt, use_prior = TRUE, use_boot = FALSE){
  
  # handle missing prior_name ----;
  if(!use_prior) dt[, prior_name := "temp_dummy_prior"]
  
  # handle missing boot variable ----;
  if(!use_boot) dt[, boot := 999]
  
  
  # get the summary quantities ----;
  dt <- copy(dt)
  
  # add means in for convenience
  dt[, `:=`(post_mean = sum(reserve * posterior), post_mean_pe = sum(reserve_pe * posterior)), by=.(prior_name, boot)]
  
  dt_res <- dt[, .(post_mean = sum(reserve * posterior),
                   post_sd = sqrt(sum(((reserve - post_mean)^2) * posterior)),
                   post_mean_pe = sum(reserve_pe * posterior),
                   post_sd_pe = sqrt(sum(((reserve_pe - post_mean_pe)^2) * posterior))), keyby = .(prior_name, boot)]
  
  # add in cov
  dt_res[, `:=`(post_cov = post_sd / post_mean, post_cov_pe = post_sd_pe / post_mean_pe)]
  
  # column ordering ----;
  setcolorder(dt_res, c("prior_name", "post_mean", "post_sd", "post_cov", "post_mean_pe", "post_sd_pe", "post_cov_pe"))
  
  # clean up ----;
  if(!use_prior) dt_res[, prior_name := NULL]
  if(!use_boot) dt_res[, boot := NULL]
  
  
  # return ----;
  return(dt_res)
}



#' GetPosteriorForecasts
#' Calculates the posterior mean cashflows for use in trimming
#'
#' @param obj 
#'
#' @return
#' @export
#'
#' @examples
GetPosteriorForecasts <- function(dt, dt_probs){
  
  dt <- copy(dt)
  
  # merge on posteriors
  dt[ dt_probs, on = .(prior_name, model_num), posterior := i.posterior]
  
  # get contribution to mean
  dt[, posterior_mean_fitted := posterior * pred]
  
  # summarise down to one set of cashflows
  key_cols <- setdiff(names(dt), c("pred", "model_num", "lambda", "posterior_mean_fitted", "posterior"))
  dt_res <- dt[, .(pred = sum(posterior_mean_fitted)), keyby = key_cols]
  
  return(dt_res)
}


# ------------------- Bootstrapping -------------------


#' Title
#'
#' @param dt 
#' @param nboot 
#' @param b1 
#' @param scale 
#' @param resid_name 
#' @param min_val 
#'
#' @return
#' @export
#'
#' @examples
SemiParametricSample <- function(dt, nboot, b1 = 1, resid_name, var_pred_name, min_val, yvar){
  
  # get vector of resisduals to sample
  resids <- dt[!is.na(get(resid_name))][[resid_name]]
  
  # get number of samples to draw
  n <- dt[past_ind == TRUE, .N] * nboot
  
  # sample with replacement  
  samp_resids <- sample(resids, n, replace = TRUE)
  
  # create framework - past and future vals, and get the bootstrapped values
  dt <- copy(dt)
  dt_boot <- dt[rep(1:.N, times = nboot)][, boot := rep(1:nboot, each = nrow(dt)) + b1 - 1]
  
  # add the residuals and generate the pseudo fitted values now
  dt_boot[past_ind == TRUE, pseudo_residuals := samp_resids]
  
  # calculate the pseudo values
  dt_boot[, (yvar) := pmax(pred + pseudo_residuals*sqrt(get(var_pred_name)), min_val)]
  
  # clean-up
  drop_vars <- intersect(unique(c("lambda.1se", "lambda.min", var_pred_name, "resid", resid_name, "pseudo_residuals")), names(dt_boot))
  dt_boot[, c(drop_vars) := NULL]
  
  # return
  return(dt_boot)
  
}


#' Title
#'
#' @param dt 
#' @param nboot 
#' @param b1 
#' @param scale 
#' @param resid_name 
#' @param var_pred_name 
#' @param min_val 
#' @param vec_lambda 
#' @param list_varset 
#' @param yvar 
#' @param dfmax 
#' @param pmax 
#' @param thresh 
#' @param list_prior_params 
#' @param dt_scaling_factor 
#' @param list_filter 
#' @param list_glmnet NULL if we need to fit the lasso models
#'
#' @return
#' @export
#'
#' @examples
BootstrapBatch <- function(dt, nboot, b1 = 1, scale, resid_name, var_pred_name, min_val = 10, 
                           vec_lambda, list_varset, yvar, dfmax, pmax, thresh = 1e-7,
                           list_prior_params, dt_scaling_factor = NULL, list_filter,
                           list_glmnet = NULL, retain_bs_data = TRUE, seed = NULL){
  
  # generate the pseudo samples if required----;
  # set a seed if required based on seed and b1 (e.g. when process called via parallel workers)
  if(!is.null(seed)) set.seed(seed + b1)
  
  # generate the pseudo samples if required----;
  if(is.null(list_glmnet)){
    dt_boot <- SemiParametricSample(dt = dt, 
                                    nboot = nboot, 
                                    b1 = 1, 
                                    resid_name = resid_name, 
                                    var_pred_name = var_pred_name, 
                                    min_val = min_val, 
                                    yvar = yvar)
  } else{
    # dt here should be the bootstrapped data sets
    if(!"boot" %in% names(dt))
      stop(paste("When running bootstrap in rescale mode, dt = first arg must be the bootstrapped data\n"))
    
    # extract the relevant bootstraps from dt_bs_data - needed for likelihood calcs
    dt_boot <- dt[boot %between% c(b1, min(b1 + nboot - 1, max(dt$boot)))]
    # renumber bootstraps from 1 to b
    dt_boot[, boot := boot - b1 + 1]
  } 
  
  # run the bootstrapping process ----;
  
  # first set up some storage
  list_probs <- vector(mode = "list", length = nboot)
  
  if(is.null(list_glmnet)){
    list_glmnet <- vector(mode = "list", length = nboot) 
  } else{
    # pull out the correct elements of list_glmnet  -ordering is important since must match data
    list_glmnet <- list_glmnet[as.character(b1:min(b1+nboot-1, length(list_glmnet)))]
  } 
  
  # now bootstrap, running each bootstrap separately
  for(b in 1:nboot){
    # fit the lassos if required
    if(is.null(list_glmnet[[b]])){
      suppressWarnings(list_glmnet[[b]] <- cv.glmnet(x = list_varset$past,
                                                     y = dt_boot[past_ind == TRUE & boot == b,][[yvar]],
                                                     lambda = vec_lambda,
                                                     family = "poisson",
                                                     dfmax = dfmax,
                                                     pmax = pmax,
                                                     thresh = thresh,	
                                                     alpha = 1,     
                                                     standardize = FALSE,
                                                     nfolds = 8, 
                                                     parallel=FALSE,
                                                     relax=FALSE,
                                                     trace.it = FALSE,  # TRUE in notebooks leads to unpleasant output
                                                     maxit = 1000000))
    }
    
    # calculate the posterior probs now for each prior
    list_tem <- vector(mode = "list", length = length(list_prior_params$prior_name))
    names(list_tem) <- list_prior_params$prior_name
    
    for(v in list_prior_params$prior_name){
      list_tem[[v]] <- CalculatePosterior(dt_boot[boot == b], 
                                          glmnet_obj = list_glmnet[[b]],
                                          varset = list_varset$all,
                                          yvar = yvar, 
                                          dt_scaling_factor = dt_scaling_factor, 
                                          list_filter = list_filter, 
                                          list_prior_params = list_prior_params, 
                                          prior_name = v,
                                          scale = scale,
                                          return_fitted = FALSE)$probs
    }
    # still all tables across the priors together
    list_probs[[b]] <- rbindlist(list_tem)
  }
  
  # join probs into tables ----;
  dt_probs <- rbindlist(list_probs, idcol = "boot")
  
  # fix up bootstrap numbers
  if(b1 > 1){ 
    dt_probs[, boot := boot + b1 - 1]
    dt_boot[,  boot := boot + b1 - 1]
  }
  
  # make glmnet names characters no matter what for now
  names(list_glmnet) <- as.character(b1:(nboot + b1 - 1))
  
  # drop data if not needed (on a rescaling run)
  if(!retain_bs_data) dt_boot <- NULL
  
  gc()
  
  # return
  return(list(list_glmnet = list_glmnet, dt_probs = dt_probs, dt_bs_data = dt_boot))
  
}

# debugonce(CalculatePosterior)
# debugonce(IntCalculatePosterior)
# debugonce(BootstrapBatch)
# x <- BootstrapBatch(dt = list_primary$data, 
#                     nboot = 2, 
#                     b1 = 1, 
#                     scale = list_primary$scale, 
#                     resid_name = "adj_resid", 
#                     var_pred_name = "var_pred",
#                     min_val = 10, 
#                     vec_lambda = lambdavec, 
#                     list_varset = list_primary$varset, 
#                     yvar = "pmts", 
#                     dfmax = my_dfmax, 
#                     pmax = my_pmax, 
#                     thresh = 1e-7,
#                     list_prior_params = list_primary$prior$list_prior_params, 
#                     dt_scaling_factor = NULL, 
#                     list_filter = list_dt_filter_params$initial)

# need a seed in parallel mode because otherwise we repeat bootstraps

ParallelBootstrap <- function(dt, nboot, scale, resid_name, var_pred_name, min_val = 10, 
                              vec_lambda, list_varset, yvar, dfmax, pmax, thresh = 1e-7,
                              list_prior_params, dt_scaling_factor = NULL, list_filter,
                              parallel = TRUE, nworkers = max(1, parallel::detectCores()/2 -1), batchsize = 10, seed = NULL,
                              list_glmnet = NULL){
  
  # timing
  ptm <- proc.time()
  
  # some parallel checks                                   
  if(parallel == TRUE){
    # only works on linux for now
    # to run on windows, drop the type="FORK" and also export required libraries and code via parallel::clusterEvalQ/clusterExport
    if(Sys.info()['sysname'] != "Linux"){
      warning("Parallel code only set up to run on linux\n Code will run sequentially")
      parallel <- FALSE
    }
    
    if(nworkers <= 1){
      warning("Not enough workers to run parallel code\n Code will run sequentially")
      parallel <- FALSE
    }
    
  }
  
  # adjust batchsize for small nboot
  if(nboot <= batchsize){
    if(parallel) {
      parallel <- FALSE
      warning("Batchsize greater than total number of bootstraps\n Code will run sequentially")
    }
    batchsize <- nboot
  }
  
  if(nboot %% batchsize != 0) warning("nboot will be increased so it is a multiple of batchsize\n")
  
  
  # setup up parallel infrastructure ---;
  if(parallel == TRUE){
    
    # error out if seed not set
    # due to forking, bootstraps will repeat in each batch otherwise....
    if(is.null(seed)) stop("Seed must be set when running in parallel to prevent identical bootstraps between concurrent batches\n")
    
    if(is.null(nworkers)) nworkers <- parallel::detectCores()/2
    workers <- parallel::makeCluster(nworkers, outfile = "parallel.txt", type = "FORK")
    on.exit(parallel::stopCluster(workers))
    doParallel::registerDoParallel(workers)
    
  } else{
    foreach::registerDoSEQ()
  }
  
  
  # batch up code now -----;
  nbatches <- ceiling(nboot/batchsize)
  
  list_res <- foreach(b = 1:nbatches, .combine = CustomCombine, .inorder = FALSE) %dopar%{
    
    
    x <- BootstrapBatch(dt = dt,
                        nboot = batchsize,
                        b1 = (b-1) * batchsize + 1,
                        scale = scale,
                        resid_name = resid_name,
                        var_pred_name = var_pred_name,
                        min_val = min_val,
                        vec_lambda = vec_lambda,
                        list_varset = list_varset,
                        yvar = yvar,
                        dfmax = dfmax,
                        pmax = pmax,
                        thresh = thresh,
                        list_prior_params = list_prior_params,
                        dt_scaling_factor = dt_scaling_factor,
                        list_filter = list_filter,
                        list_glmnet = list_glmnet,
                        retain_bs_data = is.null(list_glmnet),
                        seed = seed)
    
  }
  
  
  #browser()
  
  # fix up orderings in list_res ----;
  list_res$dt_probs[, prior_name := factor(prior_name, levels = list_prior_params$prior_name)]
  setkeyv(list_res$dt_probs, c("prior_name", "boot", "model_num"))
  
  if(!is.null(list_res$dt_bs_data)) setkeyv(list_res$dt_bs_data, c("boot", "acc", "dev", "cal"))
  
  list_res$list_glmnet <- list_res$list_glmnet[as.character(1:length(list_res$list_glmnet))]
  
  # timing info ----;
  print(paste0("Time taken: ", round((proc.time()[3] - ptm[3])/60, 4), " mins"))
  
  # return ----;    
  return(list_res)
}


# used by foreach loop above
CustomCombine <- function(a, b){
  
  ret_val <- vector(mode = "list", length = 3)
  names(ret_val) <- c("dt_probs", "list_glmnet", "dt_bs_data")
  ret_val[["dt_probs"]] <- rbindlist(list(a[["dt_probs"]], b[["dt_probs"]]))  # data.table
  ret_val[["list_glmnet"]] <- c(a[["list_glmnet"]], b[["list_glmnet"]])  # c works with list objects
  ret_val[["dt_bs_data"]] <- if(!is.null(a[["dt_bs_data"]])) rbindlist(list(a[["dt_bs_data"]], b[["dt_bs_data"]])) else NULL
  invisible(ret_val)
  
}



#debug(BootstrapBatch)

# start_time <- Sys.time()
# x <- ParallelBootstrap(dt = list_primary$data,
#                     nboot = 100,
#                     scale = list_primary$scale,
#                     resid_name = "adj_resid",
#                     var_pred_name = "var_pred",
#                     min_val = 10,
#                     vec_lambda = lambdavec,
#                     list_varset = list_primary$varset,
#                     yvar = "pmts",
#                     dfmax = my_dfmax,
#                     pmax = my_pmax,
#                     thresh = 1e-7,
#                     list_prior_params = list_primary$prior$list_prior_params,
#                     dt_scaling_factor = NULL,
#                     list_filter = list_dt_filter_params$initial,
#                     parallel = TRUE,
#                     nworkers = 5,
#                     batchsize = 20
#                     )
# 
# Sys.time()- start_time 



# assumes variable names
CalculateModelError <- function(dt, min_model_number = 5, min_model_prob_threshold = 1e-4){
  
  dt <- copy(dt)
  
  # add bootstrap mean to each obs to make calculating rest easier
  dt[, post_mean := sum(reserve * posterior), by=.(boot, prior_name)]
  dt[, post_mean_pe := sum(reserve_pe * posterior), by=.(boot, prior_name)]
  
  # now get variances values
  dt_sum <- dt[, .(post_var = sum( ((reserve - post_mean)^2)*posterior ),
                   post_var_pe = sum( ((reserve_pe - post_mean_pe)^2)*posterior )),
               keyby = .(post_mean, post_mean_pe, prior_name, boot)]
  
  # sd and cov
  dt_sum[, `:=`(post_sd = sqrt(post_var), post_sd_pe = sqrt(post_var_pe))]
  dt_sum[, `:=`(post_cov = post_sd / post_mean, post_cov_pe = post_sd_pe / post_mean_pe)]
  
  # add standalone process error now - just need means
  dt[, diff_pe := reserve_pe - reserve]
  dt_sum2 <- dt[, .(post_mean_diff_pe = sum(diff_pe * posterior)), keyby=.(prior_name, boot)]
  
  # join up all results now
  dt_sum[ dt_sum2, on = .(prior_name, boot), post_mean_diff_pe := i.post_mean_diff_pe]
  
  
  # filter out models with massive probs
  num_retained <- NULL
  if(min_model_number > 1 && min_model_prob_threshold > 0){
    bs_retain <- dt[posterior > min_model_prob_threshold, .(n = .N), keyby = .(prior_name, boot)][n >= min_model_number]
    dt_sum <- dt_sum[bs_retain, on=.(prior_name, boot), nomatch = 0]
    num_retained <- bs_retain[, .(num_bs = .N), keyby = .(prior_name)]
  }
  
  
  # calculate the grand means - with and without process error
  dt_res <- dt_sum[, .(bs_means_mean = mean(post_mean),
                       bs_variances_mean = mean(post_sd^2),
                       bs_means_variance = var(post_mean),
                       bs_covs_mean = mean(post_cov),
                       bs_means_mean_pe = mean(post_mean_pe),
                       bs_variances_mean_pe = mean(post_sd_pe^2),
                       bs_means_variance_pe = var(post_mean_pe),
                       bs_covs_mean_pe = mean(post_cov_pe),
                       bs_means_mean_diff_pe = mean(post_mean_diff_pe),
                       bs_means_variance_diff_pe = var(post_mean_diff_pe)
  ), keyby = .(prior_name)]
  
  dt_res[, `:=`(bs_variances_cov = sqrt(bs_variances_mean) / bs_means_mean,
                bs_means_cov = sqrt(bs_means_variance) / bs_means_mean,
                bs_variances_cov_pe = sqrt(bs_variances_mean_pe) / bs_means_mean_pe,
                bs_means_cov_pe = sqrt(bs_means_variance_pe) / bs_means_mean_pe,
                bs_means_cov_diff_pe = sqrt(bs_means_variance_diff_pe) / bs_means_mean)] # denom is deliberate here   
  
  # now do final calculations around model, parameter and process error
  # Bs_covs_mean gives us pure model error
  # Bs_means_cov gives us parameter error (and maybe a bit of model error)
  # Bs_means_cov_pe gives us parameter and process error 
  # We use bs_means_cov and this to get a standalone estimate of the parameter and process error variances.
  # We then get the model error variance from bs_covs_mean and add all these up to get the grand total variance (assuming independence of the 3 – probably reasonable since we think dependency between parameter and model error might already by included in bs_means_cov).
  
  dt_res[, sd_model_error := bs_covs_mean * bs_means_mean]
  dt_res[, sd_param_error := bs_means_cov * bs_means_mean]   # or just sqrt(bs_variances_mean)
  dt_res[, sd_process_error := bs_means_cov_diff_pe * bs_means_mean]   # or just sqrt(bs_means_variance_diff_pe)
  
  # get grand estimate of all errors
  dt_res[, sd_model_param_process_error := sqrt(sd_model_error^2 + sd_param_error^2 + sd_process_error^2)]
  
  # get CoVs - some will repeat other entries in the table
  dt_res[, cov_model_error := sd_model_error / bs_means_mean]
  dt_res[, cov_param_error := sd_param_error / bs_means_mean]
  dt_res[, cov_process_error := sd_process_error / bs_means_mean]
  dt_res[, cov_model_param_process_error := sd_model_param_process_error / bs_means_mean]
  
  if(!is.null(num_retained)){
    dt_res[ num_retained, on=.(prior_name), num_bs := i.num_bs]
    return(list(results = dt_res, retained_bs = bs_retain))
  } 
  
  return(list(results = res))
}






