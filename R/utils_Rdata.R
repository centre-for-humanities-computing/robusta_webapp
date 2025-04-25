require(spatstat)
require(sf) # change to sf/stars/terra after 2023
require(terra)
require(doParallel)
require(foreach)
require(tidyverse)
library(gridExtra)
library(scales)
library(ggspatial)
library(ggridges)
library(shiny)

extract_path_and_name <- function(path){
  # Remove the last 5 characters
  path_directory <- substr(path, 1, nchar(path) - 5)
  
  # Extract the 5th last character
  path_filename <- substr(path, nchar(path) - 4, nchar(path) - 4)
  
  return(list(path_directory = path_directory,
              path_filename = path_filename))
}

load_data <- function(file_shp, file_poly){
  sites <- file_shp
  sites.um <- unmark(as.ppp(sites))
  
  p <- file_poly
  pwc <- as.owin(p)
  
  sp_unmark <- sites.um[pwc]
  
  sites_asppp <- as.ppp(sites)
  
  sp_mark <- sites_asppp[pwc]
  
  return(list(sp_unmark = sp_unmark,
              sp_mark = sp_mark))
}

generate_mc_envelopes <- function(sp_unmark, clusters, nsim, portions = seq(0.5,1.0,0.1)){
  cl = makeCluster(clusters)
  registerDoParallel(cl)
  
  PO = length(portions)
  
  simulations_list = list()
  
  for (j in 1:PO){
    set.seed(j)
    indToRemove = sample(x = 1:sp_unmark$n, size = floor((1-portions[j])*sp_unmark$n), replace = FALSE)
    indToKeep = setdiff(1:sp_unmark$n, indToRemove)
    spModified = sp_unmark[indToKeep]
    simulations_list[[paste0(100*(portions[j]),"/02-_MCenveloppe1E4-")]] <- foreach(i = 1:clusters, .packages = c("spatstat")) %dopar% {
      attr(envelope(spModified, fun=pcfinhom, nsim = nsim, divisor="d", correction = "iso", savepatterns = FALSE, savefuns = TRUE), "simfuns")
    }
  }
  stopCluster(cl)
  return(simulations_list = simulations_list)
}

compute_mc_quantiles <- function(simulations_list, clusters, nsim, quantiles, quantile_50){
  
  folders = c("50","60","70","80","90","100")
  levels_low = rev(1-quantiles)
  if (quantile_50){
    levels = c(levels_low,0.5,quantiles)
  }
  else{
    levels = c(levels_low,quantiles)
  }
  
  allQuantiles_list <- list()
  allQuantilesT_list <- list()
  plot_list <- list()
  
  for (f in 1:length(folders)){
    
    allSim = c()
    simulations_ <- simulations_list[[paste0(folders[f],"/02-_MCenveloppe1E4-")]]
    for (i in 1:clusters){
      simulations <- simulations_[[i]]
      allSim = cbind(allSim,as.matrix(simulations)[,2:nsim+1])
      print(paste0(f," -- ",i))
    }
    
    allQuantiles = apply(allSim, 1, quantile, probs = levels)
    colnames(allQuantiles) = as.matrix(simulations)[,1]
    rownames(allQuantiles) = levels
    allQuantilesT = as_tibble(as.data.frame.table(allQuantiles))
    colnames(allQuantilesT) = c("level","r","value")
    allQuantilesT$r = as.numeric(as.character(allQuantilesT$r))
    allQuantiles_list[[paste0(folders[f],"-quantiles")]] <- allQuantiles
    allQuantilesT_list[[paste0(folders[f],"-quantiles")]] <- allQuantilesT
    
    p = ggplot(allQuantilesT) +
      geom_hline(yintercept = 1, lty = "dotted") +   
      geom_line(aes(x = r, y = value, col = level)) + 
      scale_y_continuous(limits = c(0,9), breaks = 0:9) + 
      ggtitle(paste0("MC quantiles for the robustness framework ",folders[f]," (", nsim*5, " MC scenarios)")) + ylab("PCF value") + 
      theme_light() + guides(col=guide_legend(nrow=3,byrow=TRUE)) + theme(legend.position = "bottom")
    
    plot_list[[paste0("temp_data/plots/",folders[f],"-quantiles")]] <- p
  }
  
  return(list(allQuantiles_list = allQuantiles_list,
              allQuantilesT_list = allQuantilesT_list,
              plot_list = plot_list))
}

generate_robustness_scenarios_a <- function(sp_unmark, nbRobSc){
  
  folders = c("50","60","70","80","90","100")
  dataPortionsToKeep = seq(from = 0.5, to = 1.0, by = 0.1)
  
  PT = length(dataPortionsToKeep)
  
  pctOriginal = pcfinhom(sp_unmark)
  RV = length(pctOriginal$r)
  
  pctShockedData_list <- list()
  shockedSp_list <- list()
  
  #t1 = Sys.time()
  for (i in 1:PT){
    set.seed(i)
    # First we generate all random samples
    indToKeep = replicate(n = nbRobSc, expr = sample(x = 1:sp_unmark$n, size = floor(dataPortionsToKeep[i]*sp_unmark$n), replace = FALSE))
    
    # Then for each random sample we perform and store the calculations
    pctShockedData = array(NA, c(RV,nbRobSc))
    
    shockedSp = array(list(),nbRobSc)
    for (j in 1:nbRobSc){
      shockedSp[[j]] = sp_unmark[indToKeep[,j]]
      pctShockedData[,j] = pcfinhom(shockedSp[[j]])$iso
      print(paste0(i," -- ",j/nbRobSc))
    }
    pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]] <- pctShockedData
    shockedSp_list[[paste0(folders[i],"-robustnessScenariosA")]] <- shockedSp
  }
  
  return(list(pctShockedData_list = pctShockedData_list,
              shockedSp_list = shockedSp_list))
}

assess_robustness_scenarios_a <- function(sp_unmark, nbRobSc, pctShockedData_list, allQuantilesT_list, quantiles){
  
  experiment = "A"
  labexperiment = "E1: uniform sampling"
  
  folders = c("50","60","70","80","90","100")
  dataPortionsToKeep = seq(0.5, 0.9,0.1)
  PT = length(dataPortionsToKeep)
  
  levelsHigh = quantiles
  LE = length(levelsHigh)
  
  pctOriginal = pcfinhom(sp_unmark)
  RV = length(pctOriginal$r)
  
  # 1. Graphique binaire ----
  allBinary = array(NA, c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    pctShockedData <- pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]]
    allQuantilesT <- allQuantilesT_list[[paste0(folders[i],"-quantiles")]]
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indFinite = which(is.finite(pctShockedData[,k]))
        allBinary[i,j,k] = any(pctShockedData[indFinite,k] > thisQuantile$value[indFinite])
      }
    }
  }
  
  proportions = apply(allBinary,c(1,2), mean)
  results = tibble(experiment = NA, removedPortion = NA, level = NA, value = NA, .rows = 0)
  colnames(proportions) = levelsHigh
  rownames(proportions) = dataPortionsToKeep
  ajout = as.data.frame(proportions)
  ajoutL = gather(ajout)
  colnames(ajoutL) = c("level","value")
  ajoutL = transform(ajoutL, experiment = experiment, keptPortion = rep(dataPortionsToKeep,LE))
  results = ajoutL[,c(3,4,1,2)]
  
  p1 = ggplot(results) + 
    geom_point(aes(x = keptPortion, y = value, col = level)) + 
    geom_line(aes(x = keptPortion, y = value, col = level)) + 
    scale_y_continuous(limits = c(0.49,1)) + 
    theme_light() + guides(col=guide_legend(nrow=1, title = "Quantile level"))+
    xlab("% of sites that are kept in each robustness scenario") + 
    ylab("% of robustness scenarios\nin which the conclusion is similar") + 
    theme(legend.position = "bottom", 
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 9, margin = margin(r = 10, unit = "pt")),
          legend.title = element_text(
            size = 9,
            margin = margin(r = 10)
          )) +
    ggtitle("First comparison tool")
  
  robustness_experiment_a01_plot <- p1
  
  # 2. Distribution des milieux ----
  allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    pctShockedData <- pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]]
    allQuantilesT <- allQuantilesT_list[[paste0(folders[i],"-quantiles")]]
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indExceed = which(pctShockedData[,k] > thisQuantile$value)
        if (length(indExceed)>1){
          difExceed = c(0,diff(indExceed))
          indJump = which(difExceed>2) # we define jumps in such a way that a single index below the threshold is not taken into account. Exemple : (3,4,5,7,8,9) is the same as (3,4,5,6,7,8,9)
          if (length(indJump)>0){
            if (length(indJump)>1){
              iStart = indExceed[indJump]
              rStart = thisQuantile$r[iStart]
              iEnd = indExceed[c(indJump[2:length(indJump)]-1,length(indExceed))]
              rEnd = thisQuantile$r[iEnd]
            } else {
              rStart = thisQuantile$r[indExceed[indJump]]
              rEnd = thisQuantile$r[indExceed[length(indExceed)]]
            }
            rIntv = cbind(rStart,rEnd)
            allFirstIntervals[[i,j,k]] = rIntv[1,]
          }
        }
      }
    }
  }
  
  originalInterval = tibble(inf = NA, sup = NA, level = levelsHigh)
  for (j in 1:LE){
    thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
    indExceed = which(pctOriginal$iso > thisQuantile$value)
    difExceed = c(0,diff(indExceed))
    indJump = which(difExceed>2)
    rStart = thisQuantile$r[indExceed[indJump]]
    rEnd = thisQuantile$r[indExceed[length(indExceed)]]
    originalInterval$inf[j] <- rStart
    originalInterval$sup[j] = rEnd
  }
  
  stats = tibble(portion = NA, level = NA, center = NA, .rows = 0)
  for (i in 1:PT){
    for (j in 1:LE){
      theseFirstIntervals = allFirstIntervals[i,j,]
      existence = sapply(theseFirstIntervals, length)
      tous = do.call(rbind, theseFirstIntervals[existence == 2])
      ajout = tibble(exp = experiment, portion = dataPortionsToKeep[i], 
                     level = levelsHigh[j], center = (tous[,1]+tous[,2])/2)
      stats = rbind(stats, ajout)
    }
  }
  
  p2 = ggplot(stats) +
    geom_density_ridges(aes(x = center, y = paste0(portion*100,"%"), height = after_stat(density), fill = as.factor(level)), linewidth = 0.25) +
    geom_vline(aes(xintercept = (inf+sup)/2, lty = "Center of original interval"), linewidth = 0.75, alpha = 0.5,
               data = originalInterval) +
    scale_linetype_manual(values = 2) + 
    scale_x_continuous(limits = c(0,8000)) +
    labs(y = "Portion of kept sites in the simulation", x = "Center of distance interval", lty = "", fill = NULL) + 
    guides(fill="none") + 
    ggtitle("Second comparison tool") +
    theme_light() + 
    facet_wrap(~paste0("Quantile level: ",level)) +
    theme(legend.position = "bottom", axis.title = element_text(size = 8), strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"))
  
  robustness_experiment_a02_plot <- p2
  
  plot_compare_exp_a <- grid.arrange(p1,p2, nrow=2, heights = c(2, 2))
  
  return(list(plot_compare_exp_a = plot_compare_exp_a))
}

comparison_tool_a <- function(sp_unmark, nbRobSc, pctShockedData_list, shockedSp_list, allQuantilesT_list, quantiles, quantile_50){
  
  folders = c("50","60","70","80","90","100")
  dataPortionsToKeep = seq(0.5, 1.0,0.1)
  PT = length(dataPortionsToKeep)
  
  pctOriginal = pcfinhom(sp_unmark)

  levelsLow = rev(1-quantiles)
  if (quantile_50){
    levels <- c(0.5,quantiles)
  }
  else{
    levels <- quantiles
  }
  LE <- length(levelsLow)
  levelsHigh <- rev(quantiles)
  
  mycolors = heat.colors(LE)
  
  example_plots_a_1 <- list()
  example_plots_a_2 <- list()
  all_plots_a <- list()
  map_plots_a <- list()
  example_map_plots_a_1 <- list()
  example_map_plots_a_2 <- list()
  combined_plots <- list()
  
  for (i in 1:PT){
    pctShockedData <- pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]]
    shockedSp <- shockedSp_list[[paste0(folders[i],"-robustnessScenariosA")]]
    allQuantilesT <- allQuantilesT_list[[paste0(folders[i],"-quantiles")]]
    
    explications = tibble(x = rep(-50, length(levels)), 
                          y = allQuantilesT$value[allQuantilesT$level %in% as.character(levels) & allQuantilesT$r == 0],
                          lab = paste0(levels * 100, "%"))
    
    p0 = ggplot(allQuantilesT) + 
      geom_blank() + theme_light()
    for (j in 1:LE){
      thisAdd = tibble(x = pctOriginal$r, ymin = subset(allQuantilesT, level == levelsLow[j])$value,
                       ymax = subset(allQuantilesT, level == levelsHigh[j])$value)
      p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                            fill = mycolors[j], alpha = 0.5, lty =  "dotted", 
                            col=alpha("black", 0.2))
    }
    
    if (quantile_50){
      p0 = p0 + geom_line(aes(x=r, y = value), data =  subset(allQuantilesT, level == 0.5), col = "blue", alpha = 0.5, lty = "dotted")
    }
    
    for (k in 1:nbRobSc){
      thisSp = shockedSp[[k]]
      p1 = ggplot() +
        layer_spatial(data = sp_unmark, aes(color = "Selected Sites"), fill = "antiquewhite", size = 0.85) +
        layer_spatial(data = thisSp, aes(color = "Sites"), fill = NA, size = 0.9) +
        scale_color_manual(
          values = c("Selected Sites" = "red", "Sites" = "black"),
          guide = guide_legend(
            override.aes = list(
              fill = c("red","black"),
              color = c("red","black"),
              shape = 21
            )
          )
        ) +
        labs(color = NULL) +
        theme_void() +
        theme(
          legend.position = "bottom",
          legend.text = element_text(size = 9, margin = margin(l = 5, r = 10, unit = "pt"))
        )
      # if folders 100 and nbRobSc 1 save plot that is slightly different
      if (folders[i] == "100" & k == 1){
        p100 = ggplot() +
          layer_spatial(data = thisSp, aes(color = "Sites"), fill = "antiquewhite", size = 1.25) +
          scale_color_manual(
            values = c("Sites" = "black"),
            guide = guide_legend(
              override.aes = list(
                fill = "black",
                color = "black",
                shape = 21
              ))) +
          labs(color = NULL) +
          theme_void() +
          theme(
            legend.position = "bottom",
            legend.text = element_text(size = 9, margin = margin(l = 5, r = 10, unit = "pt"))
          )
        example_map_plot <- p100
      }
      thisData = tibble(y = pctShockedData[,k], x = pctOriginal$r)
      p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
        scale_y_continuous(limits = c(0,4.5)) + 
        geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
        ggtitle(paste0("Robustness scenario #",k," for ",folders[i],"% of points")) + 
        xlab("r") + ylab("PCF value")
      p_combine = grid.arrange(p3,p1,ncol = 2)

      map_plots_a[[folders[i]]][[k]] <- p1
      
      all_plots_a[[folders[i]]][[k]] <- p3
      combined_plots[[folders[i]]][[k]] <- p_combine
      print(paste0(i," -- ",k))
    }

  }
  
  return(list(combined_plots = combined_plots,
              example_map_plot = example_map_plot))
}

original_pcf_function <- function(sp_unmark, nsim, quantiles){
  
  s100 <- envelope(sp_unmark, fun=pcfinhom, nsim = nsim, divisor="d", correction = "iso")
  
  
  p0 = ggplot() + 
    geom_blank() + theme_light()
  
  thisAdd = tibble(x = s100$r, ymin = s100$lo,
                   ymax = s100$hi)
  p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                        fill = "#8d8d8d", alpha = 0.5, lty =  "dotted", 
                        col=alpha("black", 0.2))
  
  p3 = p0 + geom_line(aes(x=r, y=obs), data = s100, linewidth = 0.75) + 
    scale_y_continuous(limits = c(0,4.5)) +
    ggtitle("Monte Carlo Envelope of the PCF - 100% of sites") +
    xlab("r") + ylab("PCF value")
  
  # plot(s100, ylim=c(0,4), legend=FALSE, main="Monte Carlo Envelope of the PCF - 100% of sites")
  
  # Repeat this analysis 100 times to assess its variability.
  # initialize s100l.env as an empty list
  # s100lf.env <- list()
  # 
  # # set.seed(100) # For reproducibility
  # nsim <- 100
  # for (i in 1:nsim) {
  #   # generate a new random sample each time through the loop
  #   s100lf.env[[i]] <- envelope(sp, fun=pcfinhom, nsim = 999, divisor="d", correction = "iso") # run envelope function and store result
  #   print(s100lf.env[[i]]) # print envelope result for current iteration
  # }
  # 
  # # Plot the models
  # par(mfrow=c(2,5))
  # par(mar=c(5.1,5.1,5.1,2.1))
  # 
  # for (i in 1:100) {
  #   plot(s100lf.env[[i]], ylim=c(0,4), legend=FALSE, main=paste("PCF 100% of sites -", i))
  # }
  # 
  # dev.off()
  return(p3)
}

big_processing_func <- function(file_shp, file_poly, nsim, clusters, nbRobSc, quantiles, quantile_50){
  pdf(file = NULL)
  
  withProgress(message = 'Processing Data', value = 0, {
    
    incprogress_number <- 1/6
    
    
    incProgress(incprogress_number, detail = "Loading Data")
    data <- load_data(file_shp, file_poly)
    original_pcf_plot <- original_pcf_function(sp_unmark = data[["sp_unmark"]], nsim = nsim*5, quantiles = quantiles)
    
    incProgress(incprogress_number, detail = "Generating MC envelopes")
    simulations_list <- generate_mc_envelopes(data[["sp_unmark"]], clusters = clusters, nsim = nsim, portions=seq(0.5,1.0,0.1))
    
    incProgress(incprogress_number, detail = "Computing MC quantiles")
    print("start of script 04_computeMC....R")
    output <- compute_mc_quantiles(simulations_list, clusters = clusters, nsim = nsim, quantiles = quantiles, quantile_50 = quantile_50)
    allQuantiles_list <- output[["allQuantiles_list"]]
    allQuantilesT_list <- output[["allQuantilesT_list"]]
    plot_list <- output[["plot_list"]]
    
    incProgress(incprogress_number, detail = "Generating Robustness Scenarios A")
    print("start of script 04_generateRobust...A.R")
    output <- generate_robustness_scenarios_a(data[["sp_unmark"]], nbRobSc)
    
    pctShockedData_list <- output[["pctShockedData_list"]]
    shockedSp_list <- output[["shockedSp_list"]]
    
    incProgress(incprogress_number, detail = "Robustness Experiment A")
    print("start of script 07_assessRobust...A.R")
    output_plots <- assess_robustness_scenarios_a(data[["sp_unmark"]], 
                                                            nbRobSc, 
                                                            pctShockedData_list, 
                                                            allQuantilesT_list, 
                                                            quantiles = quantiles)
    
    plot_compare_exp_a <- output_plots[["plot_compare_exp_a"]]

    incProgress(incprogress_number, detail = "Comparison Tools")
    print("start of script 09_Comparison...E1.R")
    output <- comparison_tool_a(data[["sp_unmark"]], nbRobSc, pctShockedData_list, shockedSp_list, allQuantilesT_list, quantiles = quantiles, quantile_50 = quantile_50)
    combined_plots <- output[["combined_plots"]]
    example_map_plot <- output[["example_map_plot"]]
    original_pcf_plot_combined <- grid.arrange(original_pcf_plot, example_map_plot, ncol=2)
    
  })
  
  return(list(original_pcf_plot_combined, combined_plots, plot_compare_exp_a, nbRobSc))
}

