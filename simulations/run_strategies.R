# ==============================================================================
# Analysis of public debt management strategies
# ==============================================================================


outputs <- c("mean_d","stdv_d","DaR95","mean_rr","stdv_rr",
             "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
# Convert the names of output into "expressions" for charts labels & titles:
outputs4charts <- c(expression(paste(E(d),sep="")),
                    expression(paste(sqrt(V(d)),sep="")),
                    expression(paste(q[95](d),sep="")),
                    expression(paste(E(r),sep="")),
                    expression(paste(sqrt(V(r)),sep="")),
                    expression(paste(sqrt(Delta(d)),sep="")),
                    expression(paste(E(PD[10]),sep="")),
                    expression(paste(E(spd[10]),sep="")))
# For Latex table:
latex.column.names <- c("$\\mathbb{E}(d)$","$\\sqrt{\\mathbb{V}(d)}$",
                        "$q_{95}(d)$",
                        "$\\mathbb{E}(r)$","$\\sqrt{\\mathbb{V}(r)}$",
                        "$\\sqrt{\\mathbb{V}(\\Delta d)}$",
                        "$\\mathbb{E}(PD)$","$\\mathbb{E}(spd)$")


grids <- make_grid(nb_grid = nb_grid,
                   min_d  = min_d,
                   max_d  = max_d,
                   min_rr = min_rr,
                   max_rr = max_rr,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

strategy_grid <- expand.grid(chi = values_of_chi,
                             kappa_pi = values_of_kappa_pi,
                             kappa_y = values_of_kappa_y)
if(exists("additional_strategies")){
  additional_strategies <- additional_strategies[, c("chi", "kappa_pi", "kappa_y")]
  strategy_grid <- unique(rbind(strategy_grid, additional_strategies))
}
nb_strategies <- nrow(strategy_grid)
use_parallel <- isTRUE(indic_parallel_strategies) &&
  .Platform$OS.type != "windows" &&
  nb_cores_strategies > 1
nb_cores_to_use <- if(use_parallel){
  min(nb_cores_strategies, nb_strategies)
}else{
  1
}
start_time <- Sys.time()

format_clock_time <- function(x){
  format(x, "%Y-%m-%d %H:%M:%S")
}

format_minutes <- function(x){
  if(is.na(x)){
    return("unknown")
  }
  if(x < 60){
    return(sprintf("%.1f minutes", x))
  }
  sprintf("%.1f hours", x / 60)
}

previous_elapsed_minutes <- NA_real_
if(file.exists("results/results_strategies.Rda")){
  previous_results <- new.env(parent = emptyenv())
  load("results/results_strategies.Rda", envir = previous_results)
  if(exists("strategy_run_info", envir = previous_results)){
    previous_info <- get("strategy_run_info", envir = previous_results)
    if(isTRUE(previous_info$nb_strategies == nb_strategies) &&
       isTRUE(previous_info$nb_cores == nb_cores_to_use)){
      previous_elapsed_minutes <- previous_info$elapsed_minutes
    }
  }
}

message("Starting strategy run: ", nb_strategies, " strategies.")
message("Result matrix will be ", nb_strategies, " x ", length(outputs),
        " (strategies x performance metrics).")
message("Start time: ", format_clock_time(start_time))
message("Worker processes: ", nb_cores_to_use)
if(is.finite(previous_elapsed_minutes)){
  estimated_end_time <- start_time + previous_elapsed_minutes * 60
  message("Estimated duration from previous comparable run: ",
          format_minutes(previous_elapsed_minutes),
          " (estimated finish: ", format_clock_time(estimated_end_time), ").")
}else{
  message("No previous comparable timing found; ETA will be updated after the first batch.")
}


get_strategy_output <- function(strategy, output_name){
  if(output_name == "avg_PD[maxH]"){
    return(strategy$avg_PD[maxH])
  }
  if(output_name == "avg_spreads[maxH]"){
    return(strategy$avg_spreads[maxH])
  }
  strategy[[output_name]]
}

run_one_strategy <- function(row_id){
  chi <- strategy_grid$chi[row_id]
  kappa_pi <- strategy_grid$kappa_pi[row_id]
  kappa_y <- strategy_grid$kappa_y[row_id]
  
  if(!isTRUE(indic_parallel_strategies)){
    message(sprintf("[%d/%d] chi = %.4f, kappa_pi = %.4f, kappa_y = %.4f",
                    row_id, nb_strategies, chi, kappa_pi, kappa_y))
  }
  Model_i <- Model
  Model_i$kappa_pi <- kappa_pi
  Model_i$kappa_y  <- kappa_y
  Model_i$chi      <- chi
  
  # ============================================
  # To model liquidity risks:
  Model_i$delta <- Model$delta - .0 * kappa_pi
  # ============================================
  
  Model_solved_i <- solve_ToyModel(Model_i,grids,
                                   nb_iter = nb_iter,
                                   nb_iter_sdf = nb_iter_sdf)
  
  # -----------------------------------
  p <- compute_uncond_distri(Model_solved_i$indicators_x,
                             Model_solved_i$Probas,nb_iter4probas)
  if(isTRUE(indic_strategy_diagnostic_plots)){
    distri_d  <- compute_distri_x(grids$all_d,Model_solved_i$d,p)
    plot(grids$all_d,distri_d,type="l")
    distri_rr  <- compute_distri_x(grids$all_rr,Model_solved_i$rr,p)
    plot(grids$all_rr,distri_rr,type="l")
  }
  # -----------------------------------
  
  strat_i <- run_strategy(Model_solved_i,maxH=maxH,
                          nb_iter4probas = nb_iter4probas,
                          nb_iter_sdf = nb_iter_sdf,
                          p = p)
  
  thisLine <- vapply(outputs, function(output){
    get_strategy_output(strat_i, output)
  }, numeric(1))
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  if(!isTRUE(indic_parallel_strategies)){
    message(sprintf("Finished strategy %d/%d (elapsed %.1f minutes).",
                    row_id, nb_strategies, elapsed))
  }
  
  list(parameters = c(chi,kappa_pi,kappa_y),
       outputs = thisLine)
}

if(use_parallel){
  message("Parallel strategy run using ", nb_cores_to_use, " worker processes.")
  strategy_batches <- split(seq_len(nb_strategies),
                            ceiling(seq_len(nb_strategies) / nb_cores_to_use))
  all_results <- vector("list", nb_strategies)
  for(batch_id in seq_along(strategy_batches)){
    rows <- strategy_batches[[batch_id]]
    message(sprintf("Starting batch %d/%d: strategies %d-%d.",
                    batch_id, length(strategy_batches), min(rows), max(rows)))
    batch_results <- parallel::mclapply(rows,
                                        run_one_strategy,
                                        mc.cores = min(nb_cores_to_use, length(rows)),
                                        mc.preschedule = FALSE)
    if(any(vapply(batch_results, inherits, logical(1), "try-error"))){
      stop("At least one parallel worker failed; rerun with indic_parallel_strategies <- FALSE for a detailed traceback.")
    }
    all_results[rows] <- batch_results
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(sprintf("Completed %d/%d strategies (elapsed %.1f minutes).",
                    max(rows), nb_strategies, elapsed))
    if(batch_id == 1){
      estimated_total_minutes <- elapsed * length(strategy_batches)
      estimated_end_time <- start_time + estimated_total_minutes * 60
      message("Estimated total duration after first batch: ",
              format_minutes(estimated_total_minutes),
              " (estimated finish: ", format_clock_time(estimated_end_time), ").")
    }
  }
}else{
  message("Sequential strategy run.")
  all_results <- vector("list", nb_strategies)
  for(row_id in seq_len(nb_strategies)){
    all_results[[row_id]] <- run_one_strategy(row_id)
    if(row_id == 1){
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      estimated_total_minutes <- elapsed * nb_strategies
      estimated_end_time <- start_time + estimated_total_minutes * 60
      message("Estimated total duration after first strategy: ",
              format_minutes(estimated_total_minutes),
              " (estimated finish: ", format_clock_time(estimated_end_time), ").")
    }
  }
}

parameters <- do.call(rbind, lapply(all_results, `[[`, "parameters"))
M <- do.call(rbind, lapply(all_results, `[[`, "outputs"))

colnames(parameters) <- c("chi","kappa_pi","kappa_y")
colnames(M)          <- outputs

if(isTRUE(indic_strategy_diagnostic_plots)){
  plot(M[,"mean_d"],M[,"stdv_d"])
  plot(M[,"mean_rr"],M[,"avg_PD[maxH]"])
  plot(M[,"mean_rr"],M[,"DaR95"])
  
  plot(M[(parameters[,"kappa_y"]==0),"mean_d"],M[(parameters[,"kappa_y"]==0),"stdv_d"])
  plot(M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"mean_rr"],
       M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"avg_PD[maxH]"])
}

end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
strategy_run_info <- list(start_time = start_time,
                          end_time = end_time,
                          elapsed_minutes = elapsed_minutes,
                          nb_strategies = nb_strategies,
                          nb_outputs = length(outputs),
                          nb_cores = nb_cores_to_use,
                          nb_grid = nb_grid,
                          nb_iter = nb_iter,
                          nb_iter_sdf = nb_iter_sdf,
                          nb_iter4probas = nb_iter4probas,
                          maxH = maxH)

# Save results:
save(parameters,M,Model,grids,strategy_run_info,
     file="results/results_strategies.Rda")
message("End time: ", format_clock_time(end_time))
message("Total strategy-run duration: ", format_minutes(elapsed_minutes), ".")
message("Saved strategy results to results/results_strategies.Rda.")
