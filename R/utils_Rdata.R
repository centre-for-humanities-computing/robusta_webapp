require(spatstat)
require(rgdal) # change to sf/stars/terra after 2023
require(maptools)
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
  #sites <-readOGR(dsn = filepath_shp, layer = file_shp)
  sites <- file_shp
  sites.um <- unmark(as.ppp(sites))
  
  #p <-readOGR(dsn = filepath_poly, layer = file_poly)
  p <- file_poly
  pwc <- as.owin(p)
  
  sp_unmark <- sites.um[pwc]
  
  # save(sp_unmark, file = "temp_data/RData/01_SpObject.RData") ###
  
  sites_asppp <- as.ppp(sites)
  
  sp_mark <- sites_asppp[pwc]
  
  # save(sp_mark, file = "temp_data/RData/01_SpObject_withmarks.RData") ###
  
  return(list(sp_unmark = sp_unmark,
              sp_mark = sp_mark))
}

generate_mc_envelopes <- function(sp_unmark, clusters, nsim, portions = seq(0.5,1.0,0.1)){
  #cl = makeCluster(10)
  cl = makeCluster(clusters)
  registerDoParallel(cl)
  
  # load(file = "temp_data/RData/01_SpObject.RData")
  
  #cat(paste0("Debut : ",Sys.time(),"\n"),file = paste0("output/manRun",numero,".txt"))
  
  # portions = seq(0.5,1.0,0.1)
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
  # levels = c(0.001,0.002,0.005,0.01,0.02,0.05,0.1)
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
    #allSim <- do.call(cbind, lapply(simulations_list[[f]], function(x) as.matrix(x)[,1:nsim+1]))
    
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
    # save(allQuantiles, allQuantilesT, file = paste0("temp_data/RData/",folders[f],"-quantiles.RData"))
    allQuantiles_list[[paste0(folders[f],"-quantiles")]] <- allQuantiles
    allQuantilesT_list[[paste0(folders[f],"-quantiles")]] <- allQuantilesT
    
    pdf(paste0("temp_data/plots/",folders[f],"-quantiles.pdf"))
    p = ggplot(allQuantilesT) +
      geom_hline(yintercept = 1, lty = "dotted") +   
      geom_line(aes(x = r, y = value, col = level)) + 
      scale_y_continuous(limits = c(0,9), breaks = 0:9) + 
      # geom_text(aes(x = r+1, y = value, col = level, label = level), data = subset(allQuantilesT, r == max(allQuantilesT$r))) + 
      ggtitle(paste0("MC quantiles for the robustness framework ",folders[f]," (", nsim*5, " MC scenarios)")) + ylab("PCF value") + 
      theme_light() + guides(col=guide_legend(nrow=3,byrow=TRUE)) + theme(legend.position = "bottom")
    print(p)
    dev.off()
    
    plot_list[[paste0("temp_data/plots/",folders[f],"-quantiles")]] <- p
  }
  
  return(list(allQuantiles_list = allQuantiles_list,
              allQuantilesT_list = allQuantilesT_list,
              plot_list = plot_list))
}

generate_robustness_scenarios_a <- function(sp_unmark, nbRobSc){
  # load(file = "temp_data/RData/01_SpObject.RData")
  
  folders = c("50","60","70","80","90","100")
  dataPortionsToKeep = seq(from = 0.5, to = 1.0, by = 0.1)
  
  PT = length(dataPortionsToKeep)
  #nbRobSc = 10000
  #nbRobSc = 10
  
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
    # save(pctShockedData, shockedSp, file = paste0("temp_data/RData/",folders[i],"-robustnessScenariosA.RData"))
    pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]] <- pctShockedData
    shockedSp_list[[paste0(folders[i],"-robustnessScenariosA")]] <- shockedSp
  }
  
  return(list(pctShockedData_list = pctShockedData_list,
              shockedSp_list = shockedSp_list))
}

assess_robustness_scenarios_a <- function(sp_unmark, nbRobSc, pctShockedData_list, allQuantilesT_list, quantiles){

  # load(file = "temp_data/RData/01_SpObject.RData")
  
  experiment = "A"
  labexperiment = "E1: uniform sampling"
  
  folders = c("50","60","70","80","90","100")
  dataPortionsToKeep = seq(0.5, 0.9,0.1)
  PT = length(dataPortionsToKeep)
  
  levelsHigh = quantiles
  LE = length(levelsHigh)
  
  #nbRobSc = 10000   
  # nbRobSc = 10
  
  pctOriginal = pcfinhom(sp_unmark)
  RV = length(pctOriginal$r)
  
  # 1. Graphique binaire ----
  allBinary = array(NA, c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    # load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    # load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
    pctShockedData <- pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]]
    allQuantilesT <- allQuantilesT_list[[paste0(folders[i],"-quantiles")]]
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indFinite = which(is.finite(pctShockedData[,k]))
        allBinary[i,j,k] = any(pctShockedData[indFinite,k] > thisQuantile$value[indFinite])
        # if (k %% 1000 == 0){
        #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
        # }
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
  
  pdf("temp_data/plots/07_robustnessExperimentA01.pdf", width = 7, height = 5)
  p1 = ggplot(subset(results, level != 0.98)) + 
    geom_point(aes(x = keptPortion, y = value, col = level)) + 
    geom_line(aes(x = keptPortion, y = value, col = level)) + 
    scale_y_continuous(limits = c(0.49,1)) + 
    theme_light() + guides(col=guide_legend(nrow=1, title = "quantile level"))+
    xlab("% of sites that are kept in each robustness scenario") + ylab("% of robustness scenarios\nin which the conclusion is similar") + 
    theme(legend.position = "bottom") + ggtitle("First comparison tool")
  print(p1)
  dev.off()
  
  robustness_experiment_a01_plot <- p1
  
  # 2. Distribution des milieux ----
  allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    # load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    # load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
    pctShockedData <- pctShockedData_list[[paste0(folders[i],"-robustnessScenariosA")]]
    allQuantilesT <- allQuantilesT_list[[paste0(folders[i],"-quantiles")]]
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indExceed = which(pctShockedData[,k] > thisQuantile$value)
        if (length(indExceed)>1){
          difExceed = c(0,diff(indExceed))
          # indJump = which(difExceed>1)
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
        # if (k %% 100 == 0){
        #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
        # }
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
    originalInterval$inf[j] = rStart
    originalInterval$sup[j] = rEnd
  }
  # cbind(1:513,pctOriginal$r,pctOriginal$iso - thisQuantile$value) 
  
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
  
  pdf("temp_data/plots/07_robustnessExperimentA02.pdf", width = 7, height = 5)
  p2 = ggplot(subset(stats, level != 0.98)) +
    geom_density_ridges(aes(x = center, y = paste0(portion*100,"%"), height = after_stat(density), fill = as.factor(level)), linewidth = 0.25) +
    geom_vline(aes(xintercept = (inf+sup)/2, lty = "Center of original interval"), linewidth = 0.75, alpha = 0.5,
               data = subset(originalInterval, level != 0.98)) +
    scale_linetype_manual(values = 2) + 
    scale_x_continuous(limits = c(0,8000)) + 
    labs(y = "Portion of kept sites in the simulation", x = "Center of distance interval", lty = "", fill = NULL) + 
    guides(fill="none") + 
    ggtitle("Second comparison tool") +
    theme_light() + 
    facet_wrap(~paste0("Quantile level: ",level)) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"))
  print(p2)
  dev.off()
  
  robustness_experiment_a02_plot <- p2
  
  # Et on met tout ensemble : 
  pdf("temp_data/plots/07_robustnessExperimentA99.pdf", width = 7, height = 10)
  grid.arrange(p1,p2)
  dev.off()
  
  # plot_compare_exp_a <- arrangeGrob(p1,p2,nrow=2)
  plot_compare_exp_a <- list(p1, p2)
  
  return(list(robustness_experiment_a01_plot = robustness_experiment_a01_plot,
              robustness_experiment_a02_plot = robustness_experiment_a02_plot,
              plot_compare_exp_a = plot_compare_exp_a))
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
    
    pdf(paste0("temp_data/plots/11_examplePlots",folders[i],"-01.pdf"), width = 7, height = 5)
    for (k in 1:nbRobSc){
      thisSp = shockedSp[[k]]
      p1 = ggplot() +
        layer_spatial(data = sp_unmark, col = "red", fill="antiquewhite") +
        layer_spatial(data = thisSp, fill = NA, col = "black") +
        theme_void()
      thisData = tibble(y = pctShockedData[,k], x = pctOriginal$r)
      p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
        scale_y_continuous(limits = c(0,4.5)) + 
        geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
        ggtitle(paste0("Scenario #",k," for ",folders[i],"% of points")) + 
        xlab("r") + ylab("PCF value")
      # grid.arrange(p1,p3,nrow = 2)
      # ptot = p3 + annotation_custom(ggplotGrob(p1), xmin = 1750, xmax = 8000, 
      #                               ymin = 2, ymax = 4.5)
      # print(ptot)
      # if (k == 1){
      #   example_plots_a[[i]] <- ptot
      # }
      map_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[k]] <- p1
      all_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[k]] <- p3
      print(paste0(i," -- ",k))
    }
    dev.off()
    
    # Randomly select one plot from all_plots_a and assign it to example_plots_a
    random_index = sample(seq_along(all_plots_a[[paste0("11_examplePlots",folders[i],"-01")]]), 2)
    example_plots_a_1[[i]] <- all_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[random_index[[1]]]]
    example_map_plots_a_1[[i]] <- map_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[random_index[[1]]]]
    example_plots_a_2[[i]] <- all_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[random_index[[2]]]]
    example_map_plots_a_2[[i]] <- map_plots_a[[paste0("11_examplePlots",folders[i],"-01")]][[random_index[[2]]]]
    
  }
  
  return(list(example_plots_a_1 = example_plots_a_1,
              example_map_plots_a_1 = example_map_plots_a_1,
              example_plots_a_2 = example_plots_a_2,
              example_map_plots_a_2 = example_map_plots_a_2,
              all_plots_a = all_plots_a))
}

original_pcf_function <- function(sp_unmark, nsim){
  
  s100 <- envelope(sp_unmark, fun=pcfinhom, nsim = nsim, divisor="d", correction = "iso")
  
  
  p0 = ggplot() + 
    geom_blank() + theme_light()
  
  thisAdd = tibble(x = s100$r, ymin = s100$lo,
                   ymax = s100$hi)
  p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                        fill = heat.colors(1), alpha = 0.5, lty =  "dotted", 
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
  withProgress(message = 'Processing Data', value = 0, {
    
    incprogress_number <- 1/6
    
    
    incProgress(incprogress_number, detail = "Loading Data")
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    
    data <- load_data(file_shp, file_poly)
    
    ########## end of scripts 02_...R. Output two RData objects #########
    
    
    
    
    
    ########## start of scripts 01_...R. Output two RData objects #########
    original_pcf_plot <- original_pcf_function(data[["sp_unmark"]], nsim*5)
    
    
    
    ########## start of script 02_generateMC....R from Sébastien ##########
    incProgress(incprogress_number, detail = "Generating MC envelopes")
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    
    simulations_list <- generate_mc_envelopes(data[["sp_unmark"]], clusters = clusters, nsim = nsim, portions=seq(0.5,1.0,0.1))
    
    
    ########## end of script 02_...R (no parallel right now) ##########
    
    
    
    
    
    ########## start of script 04_computeMC....R from Eduardo ##########
    incProgress(incprogress_number, detail = "Computing MC quantiles")
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 04_computeMC....R")
    
    output <- compute_mc_quantiles(simulations_list, clusters = clusters, nsim = nsim, quantiles = quantiles, quantile_50 = quantile_50)
    allQuantiles_list <- output[["allQuantiles_list"]]
    allQuantilesT_list <- output[["allQuantilesT_list"]]
    plot_list <- output[["plot_list"]]
    
    ########## end of script 04_computeMC....R from Eduardo ##########
    
    
    
    
    
    
    ########## start of script 04_generateRobust...A.R from Sébastien ##########
    incProgress(incprogress_number, detail = "Generating Robustness Scenarios A")
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 04_generateRobust...A.R")
    
    output <- generate_robustness_scenarios_a(data[["sp_unmark"]], nbRobSc)
    pctShockedData_list <- output[["pctShockedData_list"]]
    shockedSp_list <- output[["shockedSp_list"]]
    
    
    ########## end of script 04_generateRobust...A.R from Sébastien ##########
    
    
    
    ########## start of script 07_assessRobust...A.R from Sébastien ##########
    incProgress(incprogress_number, detail = "Robustness Experiment A")
    print("start of script 07_assessRobust...A.R")
    
    output_plots <- assess_robustness_scenarios_a(data[["sp_unmark"]], nbRobSc, pctShockedData_list, allQuantilesT_list, quantiles = quantiles)
    plot_compare_exp_a <- output_plots[["plot_compare_exp_a"]]
    
    ########## end of script 07_assessRobust...A.R from Sébastien ##########
    
    
    
    
    
    
    

    
    
    ########## start of script 09_Comparison...E1.R from Eduardo ##########
    incProgress(incprogress_number, detail = "Comparison Tools")
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b", "plot_test")))
    print("start of script 09_Comparison...E1.R")
    
    output <- comparison_tool_a(data[["sp_unmark"]], nbRobSc, pctShockedData_list, shockedSp_list, allQuantilesT_list, quantiles = quantiles, quantile_50 = quantile_50)
    example_plots_a_1 <- output[["example_plots_a_1"]]
    example_map_plots_a_1 <- output[["example_map_plots_a_1"]]
    example_plots_a_2 <- output[["example_plots_a_2"]]
    example_map_plots_a_2 <- output[["example_map_plots_a_2"]]
    
    ########## end of script 09_Comparison...E1.R from Eduardo ##########
  })
  
  return(list(plot_compare_exp_a, example_plots_a_1, example_map_plots_a_1, example_plots_a_2, example_map_plots_a_2, original_pcf_plot))
}











big_processing_func_original <- function(filepath_shp, file_shp, filepath_poly, file_poly, nsim, n_iter, nbRobSc){
  withProgress(message = 'Processing Data', value = 0, {
    incProgress(1/9, detail = "Loading Data")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    
    #sites <-readOGR(dsn = filepath_shp, layer = file_shp)
    sites <- file_shp
    sites.um <- unmark(as.ppp(sites))
    
    #p <-readOGR(dsn = filepath_poly, layer = file_poly)
    p <- file_poly
    pwc <- as.owin(p)
    
    sp_unmark <- sites.um[pwc]
    
    save(sp_unmark, file = "temp_data/RData/01_SpObject.RData") ###
    
    sites_asppp <- as.ppp(sites)
    
    sp_mark <- sites_asppp[pwc]
    
    save(sp_mark, file = "temp_data/RData/01_SpObject_withmarks.RData") ###
    
    ########## end of scripts 02_...R. Output two RData objects #########
    
    
    
    
    
    ########## start of script 02_generateMC....R from Sébastien ##########
    incProgress(1/9, detail = "Generating MC envelopes")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    
    #cl = makeCluster(10)
    cl = makeCluster(5)
    registerDoParallel(cl)
    
    load(file = "temp_data/RData/01_SpObject.RData")
    
    #cat(paste0("Debut : ",Sys.time(),"\n"),file = paste0("output/manRun",numero,".txt"))
    
    portions = seq(0.5,1.0,0.1)
    PO = length(portions)
    
    for (j in 1:PO){
      set.seed(j)
      indToRemove = sample(x = 1:sp_unmark$n, size = floor((1-portions[j])*sp_unmark$n), replace = FALSE)
      indToKeep = setdiff(1:sp_unmark$n, indToRemove)
      spModified = sp_unmark[indToKeep]
      #print(paste0('Debut ronde ',j," (",as.character(Sys.time()),")\n"))
      foreach(i = 1:n_iter, .packages = c("spatstat")) %dopar% {
        MCenvs = envelope(spModified, fun=pcfinhom, nsim = nsim, divisor="d", correction = "iso", savepatterns = FALSE, savefuns = TRUE)
        simulations = attr(MCenvs, "simfuns")
        save(simulations, file = paste0("temp_data/CECI/RData/0",100*(portions[j]),"/02-_MCenveloppe1E4-",i,".RData"))
      }
      #print(paste0('Ronde ',j," terminee (",as.character(Sys.time()),")\n"))
    }
    
    #cat(paste0("Fin : ",Sys.time(),"\n"),file = paste0("output/manRun",numero,".txt"), append = TRUE)
    
    stopCluster(cl)
    
    
    ########## end of script 02_...R (no parallel right now) ##########
    
    
    
    
    
    ########## start of script 04_computeMC....R from Eduardo ##########
    incProgress(1/9, detail = "Computing MC quantiles")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 04_computeMC....R")
    
    folders = c("050","060","070","080","090","0100")
    #folders = c("temp")
    #temp_list = c("50","60","70","80","90","100")
    levels = c(0.001,0.002,0.005,0.01,0.02,0.05,0.1)
    levels = c(levels,0.5,rev(1-levels))
    
    for (f in 1:length(folders)){
      #allSim <- do.call(cbind, lapply(simulations_list[[f]], function(x) as.matrix(x)[,1:nsim+1]))
      
      allSim = c()
      for (i in 1:n_iter){
        load(paste0("temp_data/CECI/RData/",folders[f],"/02-_MCenveloppe1E4-",i,".RData"))
        #load(paste0("temp/RData/02_MCenveloppe100000-",temp_list[i],".RData"))
        allSim = cbind(allSim,as.matrix(simulations)[,1:nsim+1])
        print(paste0(f," -- ",i))
      }
      
      allQuantiles = apply(allSim, 1, quantile, probs = levels)
      colnames(allQuantiles) = as.matrix(simulations)[,1]
      rownames(allQuantiles) = levels
      allQuantilesT = as_tibble(as.data.frame.table(allQuantiles))
      colnames(allQuantilesT) = c("level","r","value")
      allQuantilesT$r = as.numeric(as.character(allQuantilesT$r))
      save(allQuantiles, allQuantilesT, file = paste0("temp_data/RData/",folders[f],"-quantiles.RData"))
      
      pdf(paste0("temp_data/plots/",folders[f],"-quantiles.pdf"))
      p = ggplot(allQuantilesT) +
        geom_hline(yintercept = 1, lty = "dotted") +   
        geom_line(aes(x = r, y = value, col = level)) + 
        scale_y_continuous(limits = c(0,9), breaks = 0:9) + 
        # geom_text(aes(x = r+1, y = value, col = level, label = level), data = subset(allQuantilesT, r == max(allQuantilesT$r))) + 
        ggtitle(paste0("MC quantiles for the robustness framework ",folders[f]," (", nsim*5, " MC scenarios)")) + ylab("PCF value") + 
        theme_light() + guides(col=guide_legend(nrow=3,byrow=TRUE)) + theme(legend.position = "bottom")
      print(p)
      dev.off()
    }
    
    ########## end of script 04_computeMC....R from Eduardo ##########
    
    
    
    
    
    
    
    
    
    ########## start of script 03_computeMC....R from Sébastien ##########
    # rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    # print("start of script 03_computeMC....R")
    # #nbScMC = c(100,1000,10000,100000)
    # nbScMC = c(2,5,10,nsim,nsim*n_iter)
    # NS = length(nbScMC)
    # 
    # levels = c(0.005,0.01,0.05,0.1)
    # levels = c(levels,0.5,rev(1-levels))
    # 
    # allSim = c()
    # for (i in 1:n_iter){
    #   load(paste0("temp_data/CECI/RData/0100/02-_MCenveloppe1E4-",i,".RData"))
    #   allSim = cbind(allSim,as.matrix(simulations)[,2:nsim+1])
    #   print(paste0(" -- ",i))
    # }
    # 
    # allQuantilesT = tibble(level = NA, r = NA, value = NA, nberSc = NA, repet = NA, .rows = 0)
    # # ????? why a minus 1?????
    # for (i in 1:(NS-1)){
    #   set.seed(2)
    #   indAvailable = 1:dim(allSim)[2]
    #   for (j in 1:n_iter){
    #     indSampled = sample(indAvailable, size = nbScMC[i], replace = FALSE)
    #     allQuantiles = apply(allSim[,indSampled], 1, quantile, probs = levels)
    #     colnames(allQuantiles) = as.matrix(simulations)[,1]
    #     rownames(allQuantiles) = levels
    #     ajout = as_tibble(as.data.frame.table(allQuantiles))
    #     colnames(ajout) = c("level","r","value")
    #     ajout = transform(ajout, nberSc = nbScMC[i], repet = j)
    #     allQuantilesT = rbind(allQuantilesT, ajout)
    #     indAvailable = setdiff(indAvailable, indSampled)
    #     print(paste0(i,"--",j))
    #   }
    # }
    # 
    # allQuantilesT$r = as.numeric(as.character(allQuantilesT$r))
    # allQuantilesT$nberSc = as.character(allQuantilesT$nberSc)
    # allQuantilesT$repet = as.character(allQuantilesT$repet)
    # save(allQuantilesT, file = paste0("temp_data/RData/full-quantiles.RData"))
    # 
    # 
    # essai = subset(allQuantilesT, level == 0.995 & repet == 1)
    # 
    # ggplot(essai) + 
    #   geom_line(aes(x = r, y = value, col = nberSc))
    # 
    # essai2a = subset(allQuantilesT, level == 0.995)
    # essai2b = subset(allQuantilesT, level == 0.005)
    # 
    # load(file = "temp_data/RData/01_SpObject.RData")
    # pctOriginal = pcfinhom(sp_unmark)
    # pctOriginalT = tibble(r = pctOriginal$r, value = pctOriginal$trans)
    # 
    # pdf("temp_data/plots/03_quantileComparison01.pdf", height = 10, width = 7)
    # ggplot(essai2a) + 
    #   geom_line(aes(x = r, y = value, col = repet)) +
    #   geom_line(aes(x = r, y = value, col = repet), data = essai2b) + 
    #   geom_line(aes(x = r, y = value), data = pctOriginalT, linewidth = 0.75) +
    #   coord_cartesian(ylim = c(0,3)) + 
    #   theme_light() + theme(legend.position = "none") + ggtitle(label = "MC enveloppe when varying the number of scenarios\nused to compute the 99.5% quantile", subtitle = "Each color curve corresponds to one of the 10 repetitions of the quantile computations") +  
    #   facet_wrap(~paste0("Number of simulations = ",nberSc),nrow = 3)
    # dev.off()
    
    ########## end of script 03_computeMC....R from Sébastien ##########
    
    
    
    
    
    
    
    ########## start of script 04_generateRobust...A.R from Sébastien ##########
    incProgress(1/9, detail = "Generating Robustness Scenarios A")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 04_generateRobust...A.R")
    load(file = "temp_data/RData/01_SpObject.RData")
    
    folders = c("050","060","070","080","090","0100")
    dataPortionsToKeep = seq(from = 0.5, to = 1.0, by = 0.1)
    
    PT = length(dataPortionsToKeep)
    #nbRobSc = 10000
    #nbRobSc = 10
    
    pctOriginal = pcfinhom(sp_unmark)
    RV = length(pctOriginal$r)
    
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
      save(pctShockedData, shockedSp, file = paste0("temp_data/RData/",folders[i],"-robustnessScenariosA.RData"))
    }
    #t2 = Sys.time()
    
    
    
    
    ########## end of script 04_generateRobust...A.R from Sébastien ##########
    
    
    
    
    
    
    
    
    
    ########## start of script 04_generateRobust...B.R from Sébastien ##########
    incProgress(1/9, detail = "Generating Robustness Scenarios B")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 04_generateRobust...B.R")
    
    load(file = "temp_data/RData/01_SpObject_withmarks.RData")
    
    folders = c("050","060","070","080","090","0100")
    dataPortionsToKeep = seq(from = 0.5, to = 1.0, by = 0.1)
    
    PT = length(dataPortionsToKeep)
    # ???????? this was 10000 but why? Is it connected to number of simulations??????+
    #nbRobSc = 10000
    #nbRobSc = 10
    
    pctOriginal = pcfinhom(sp_mark)
    RV = length(pctOriginal$r)
    
    #t1 = Sys.time()
    for (i in 1:PT){
      set.seed(i)
      
      # Probabilities of begin kept proportional to area of site
      probaOfBeingKept = sp_mark$marks$Area.sq.m/sum(sp_mark$marks$Area.sq.m)
      # First we generate all random samples
      indToKeep = replicate(n = nbRobSc, 
                            expr = sample(x = 1:sp_mark$n, 
                                          size = floor(dataPortionsToKeep[i]*sp_mark$n), 
                                          prob = probaOfBeingKept,
                                          replace = FALSE))
      
      # Then for each random sample we perform and store the calculations
      pctShockedData = array(NA, c(RV,nbRobSc))
      
      shockedSp = array(list(),nbRobSc)
      for (j in 1:nbRobSc){
        shockedSp[[j]] = sp_mark[indToKeep[,j]]
        pctShockedData[,j] = pcfinhom(shockedSp[[j]])$iso
        print(paste0(i," -- ",j/nbRobSc))
      }
      save(pctShockedData, shockedSp, file = paste0("temp_data/RData/",folders[i],"-robustnessScenariosB.RData"))
    }
    #t2 = Sys.time()
    
    ########## end of script 04_generateRobust...B.R from Sébastien ##########
    
    
    
    
    
    
    
    
    
    ########## start of script 07_assessRobust...A.R from Sébastien ##########
    incProgress(1/9, detail = "Robustness Experiment A")
    print("start of script 07_assessRobust...A.R")
    
    load(file = "temp_data/RData/01_SpObject.RData")
    
    experiment = "A"
    labexperiment = "E1: uniform sampling"
    
    folders = c("050","060","070","080","090","100")
    dataPortionsToKeep = seq(0.5, 0.9,0.1)
    PT = length(dataPortionsToKeep)
    
    levelsLow = c(0.005,0.01,0.02,0.05,0.1)
    LE = length(levelsLow)
    levelsHigh = 1-levelsLow
    
    #nbRobSc = 10000   
    nbRobSc = 10
    
    pctOriginal = pcfinhom(sp_mark)
    RV = length(pctOriginal$r)
    
    # 1. Graphique binaire ----
    allBinary = array(NA, c(PT,LE,nbRobSc))
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      for (j in 1:LE){
        thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
        for (k in 1:nbRobSc){
          indFinite = which(is.finite(pctShockedData[,k]))
          allBinary[i,j,k] = any(pctShockedData[indFinite,k] > thisQuantile$value[indFinite])
          # if (k %% 1000 == 0){
          #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
          # }
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
    
    pdf("temp_data/plots/07_robustnessExperimentA01.pdf", width = 7, height = 5)
    p1 = ggplot(subset(results, level != 0.98)) + 
      geom_point(aes(x = keptPortion, y = value, col = level)) + 
      geom_line(aes(x = keptPortion, y = value, col = level)) + 
      scale_y_continuous(limits = c(0.49,1)) + 
      theme_light() + guides(col=guide_legend(nrow=1, title = "quantile level"))+
      xlab("% of sites that are kept in each robustness scenario") + ylab("% of robustness scenarios\nin which the conclusion is similar") + 
      theme(legend.position = "bottom") + ggtitle("First comparison tool")
    print(p1)
    dev.off()
    
    # 2. Distribution des milieux ----
    allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      for (j in 1:LE){
        thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
        for (k in 1:nbRobSc){
          indExceed = which(pctShockedData[,k] > thisQuantile$value)
          if (length(indExceed)>1){
            difExceed = c(0,diff(indExceed))
            # indJump = which(difExceed>1)
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
          # if (k %% 100 == 0){
          #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
          # }
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
      originalInterval$inf[j] = rStart
      originalInterval$sup[j] = rEnd
    }
    # cbind(1:513,pctOriginal$r,pctOriginal$iso - thisQuantile$value) 
    
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
    
    pdf("temp_data/plots/07_robustnessExperimentA02.pdf", width = 7, height = 5)
    p2 = ggplot(subset(stats, level != 0.98)) +
      geom_density_ridges(aes(x = center, y = paste0(portion*100,"%"), height = ..density.., fill = as.factor(level)), size = 0.25) +
      geom_vline(aes(xintercept = (inf+sup)/2, lty = "Center of original interval"), size = 0.75, alpha = 0.5,
                 data = subset(originalInterval, level != 0.98)) +
      scale_linetype_manual(values = 2) + 
      scale_x_continuous(limits = c(0,8000)) + 
      labs(y = "Portion of kept sites in the simulation", x = "Center of distance interval", lty = "", fill = NULL) + 
      guides(fill="none") + 
      ggtitle("Second comparison tool") +
      theme_light() + 
      facet_wrap(~paste0("Quantile level: ",level)) +
      theme(legend.position = "bottom", strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"))
    print(p2)
    dev.off()
    
    
    # Et on met tout ensemble : 
    pdf("temp_data/plots/07_robustnessExperimentA99.pdf", width = 7, height = 10)
    grid.arrange(p1,p2)
    dev.off()
    
    plot_compare_exp_b <- arrangeGrob(p1,p2,nrow=2)
    
    ########## end of script 07_assessRobust...A.R from Sébastien ##########
    
    
    
    
    
    
    ########## start of script 07_assessRobust...B.R from Sébastien ##########
    incProgress(1/9, detail = "Robustness Experiment B")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b")))
    print("start of script 07_assessRobust...B.R")
    
    
    load(file = "temp_data/RData/01_SpObject.RData")
    
    experiment = "B"
    labexperiment = "E1: uniform sampling"
    
    folders = c("050","060","070","080","090","100")
    dataPortionsToKeep = seq(0.5, 0.9,0.1)
    PT = length(dataPortionsToKeep)
    
    levelsLow = c(0.005,0.01,0.02,0.05,0.1)
    LE = length(levelsLow)
    levelsHigh = 1-levelsLow
    
    #nbRobSc = 10000   
    nbRobSc = 10
    
    pctOriginal = pcfinhom(sp_unmark)
    RV = length(pctOriginal$r)
    
    # 1. Graphique binaire ----
    allBinary = array(NA, c(PT,LE,nbRobSc))
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      for (j in 1:LE){
        thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
        for (k in 1:nbRobSc){
          indFinite = which(is.finite(pctShockedData[,k]))
          allBinary[i,j,k] = any(pctShockedData[indFinite,k] > thisQuantile$value[indFinite])
          # if (k %% 1000 == 0){
          #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
          # }
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
    
    pdf("temp_data/plots/07_robustnessExperimentB01.pdf", width = 7, height = 5)
    p1 = ggplot(subset(results, level != 0.98)) + 
      geom_point(aes(x = keptPortion, y = value, col = level)) + 
      geom_line(aes(x = keptPortion, y = value, col = level)) + 
      scale_y_continuous(limits = c(0.49,1)) + 
      theme_light() + guides(col=guide_legend(nrow=1, title = "quantile level"))+
      xlab("% of sites that are kept in each robustness scenario") + ylab("% of robustness scenarios\nin which the conclusion is similar") + 
      theme(legend.position = "bottom") + ggtitle("First comparison tool")
    print(p1)
    dev.off()
    
    # 2. Distribution des milieux ----
    allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      for (j in 1:LE){
        thisQuantile = subset(allQuantilesT, level == levelsHigh[j])
        for (k in 1:nbRobSc){
          indExceed = which(pctShockedData[,k] > thisQuantile$value)
          if (length(indExceed)>1){
            difExceed = c(0,diff(indExceed))
            # indJump = which(difExceed>1)
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
          # if (k %% 100 == 0){
          #   print(paste0("Rob framework ",i," ; level ",j," ; scenario ",k)) 
          # }
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
      originalInterval$inf[j] = rStart
      originalInterval$sup[j] = rEnd
    }
    # cbind(1:513,pctOriginal$r,pctOriginal$iso - thisQuantile$value) 
    
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
    
    pdf("temp_data/plots/07_robustnessExperimentB02.pdf", width = 7, height = 5)
    p2 = ggplot(subset(stats, level != 0.98)) +
      geom_density_ridges(aes(x = center, y = paste0(portion*100,"%"), height = ..density.., fill = as.factor(level)), size = 0.25) +
      geom_vline(aes(xintercept = (inf+sup)/2, lty = "Center of original interval"), size = 0.75, alpha = 0.5,
                 data = subset(originalInterval, level != 0.98)) +
      scale_linetype_manual(values = 2) + 
      scale_x_continuous(limits = c(0,8000)) + 
      labs(y = "Portion of kept sites in the simulation", x = "Center of distance interval", lty = "", fill = NULL) + 
      guides(fill="none") + 
      ggtitle("Second comparison tool") +
      theme_light() + 
      facet_wrap(~paste0("Quantile level: ",level)) +
      theme(legend.position = "bottom", strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black"))
    print(p2)
    dev.off()
    
    
    # Et on met tout ensemble : 
    pdf("temp_data/plots/07_robustnessExperimentB99.pdf", width = 7, height = 10)
    grid.arrange(p1,p2)
    dev.off()
    
    plot_compare_exp_a <- arrangeGrob(p1,p2,nrow=2)
    
    
    ########## end of script 07_assessRobust...B.R from Sébastien ##########
    
    
    
    
    
    
    ########## start of script 08_Otherplots.R from Eduardo ########
    
    
    
    ########## end of script 08_Otherplots.R from Eduardo ########
    
    
    
    
    
    
    
    ########## start of script 09_Comparison...E1.R from Eduardo ##########
    incProgress(1/9, detail = "Comparison Tools")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b", "plot_test")))
    print("start of script 09_Comparison...E1.R")
    
    folders = c("050","060","070","080","090","0100")
    dataPortionsToKeep = seq(0.5, 1.0,0.1)
    PT = length(dataPortionsToKeep)
    
    load(file = "temp_data/RData/01_SpObject.RData")
    pctOriginal = pcfinhom(sp_unmark)
    
    load("temp_data/RData/0100-quantiles.RData")
    levelsLow = c(0.005,0.01,0.02,0.05,0.1)
    LE = length(levelsLow)
    levelsHigh = 1-levelsLow
    
    mycolors = (heat.colors(LE))
    
    example_plots_a <- list()
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenariosA.RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      
      explications = tibble(x = rep(-50, 6), 
                            y = allQuantilesT$value[allQuantilesT$level %in% as.character(c(0.5,0.9,0.95,0.98,0.99,0.995)) & allQuantilesT$r == 0],
                            lab = c('50%','90%','95%','98%','99%','99.5%'))
      
      p0 = ggplot(allQuantilesT) + 
        geom_blank() + theme_light()
      for (j in 1:LE){
        thisAdd = tibble(x = pctOriginal$r, ymin = subset(allQuantilesT, level == levelsLow[j])$value,
                         ymax = subset(allQuantilesT, level == levelsHigh[j])$value)
        p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                              fill = mycolors[j], alpha = 0.5, lty =  "dotted", 
                              col=alpha("black", 0.2))
      }
      p0 = p0 + geom_line(aes(x=r, y = value), data =  subset(allQuantilesT, level == 0.5), col = "blue", alpha = 0.5, lty = "dotted")
      
      pdf(paste0("temp_data/plots/11_examplePlots",folders[i],"-01.pdf"), width = 7, height = 5)
      for (k in 1:nbRobSc){
        thisSp = shockedSp[[k]]
        p1 = ggplot() +
          layer_spatial(data = sp_unmark, col = "red", fill="antiquewhite") +
          layer_spatial(data = thisSp, fill = NA, col = "black") +
          theme_void()
        thisData = tibble(y = pctShockedData[,k], x = pctOriginal$r)
        p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
          scale_y_continuous(limits = c(0,4.5)) + 
          geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
          ggtitle(paste("Scenario #",k,"for",folders[i],"% of points")) + 
          xlab("r") + ylab("PCF value")
        # grid.arrange(p1,p3,nrow = 2)
        ptot = p3 + annotation_custom(ggplotGrob(p1), xmin = 1750, xmax = 8000, 
                                      ymin = 2, ymax = 4.5)
        print(ptot)
        if (k == 1){
          example_plots_a[[i]] <- ptot
        }
        print(paste0(i," -- ",k))
      }
      dev.off()
      
    }
    
    ########## end of script 09_Comparison...E1.R from Eduardo ##########
    
    
    
    
    
    
    
    ########## start of script 09_Comparison...E2.R from Eduardo ##########
    incProgress(1/9, detail = "Additional Plots")
    rm(list=setdiff(ls(), c("filepath_shp", "file_shp", "filepath_poly", "file_poly", "nsim", "n_iter", "nbRobSc", "plot_compare_exp_a", "plot_compare_exp_b", "plot_test", "example_plots_a")))
    print("start of script 09_Comparison...E2.R")
    
    folders = c("050","060","070","080","090","0100")
    dataPortionsToKeep = seq(0.5, 1.0,0.1)
    PT = length(dataPortionsToKeep)
    
    load(file = "temp_data/RData/01_SpObject.RData")
    pctOriginal = pcfinhom(sp_unmark)
    
    load("temp_data/RData/0100-quantiles.RData")
    levelsLow = c(0.005,0.01,0.02,0.05,0.1)
    LE = length(levelsLow)
    levelsHigh = 1-levelsLow
    
    mycolors = (heat.colors(LE))
    
    example_plots_b <- list()
    
    for (i in 1:PT){
      load(paste0("temp_data/RData/",folders[i],"-robustnessScenariosB.RData"))
      load(paste0("temp_data/RData/",folders[i],"-quantiles.RData"))
      
      explications = tibble(x = rep(-50, 6), 
                            y = allQuantilesT$value[allQuantilesT$level %in% as.character(c(0.5,0.9,0.95,0.98,0.99,0.995)) & allQuantilesT$r == 0],
                            lab = c('50%','90%','95%','98%','99%','99.5%'))
      
      p0 = ggplot(allQuantilesT) + 
        geom_blank() + theme_light()
      for (j in 1:LE){
        thisAdd = tibble(x = pctOriginal$r, ymin = subset(allQuantilesT, level == levelsLow[j])$value,
                         ymax = subset(allQuantilesT, level == levelsHigh[j])$value)
        p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                              fill = mycolors[j], alpha = 0.5, lty =  "dotted", 
                              col=alpha("black", 0.2))
      }
      p0 = p0 + geom_line(aes(x=r, y = value), data =  subset(allQuantilesT, level == 0.5), col = "blue", alpha = 0.5, lty = "dotted")
      
      pdf(paste0("temp_data/plots/11E2_examplePlots",folders[i],"-01.pdf"), width = 7, height = 5)
      for (k in 1:nbRobSc){
        thisSp = shockedSp[[k]]
        p1 = ggplot() +
          layer_spatial(data = sp_unmark, col = "red", fill="antiquewhite") +
          layer_spatial(data = thisSp, fill = NA, col = "black") +
          theme_void()
        thisData = tibble(y = pctShockedData[,k], x = pctOriginal$r)
        p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
          scale_y_continuous(limits = c(0,4.5)) + 
          geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
          ggtitle(paste("Scenario #",k,"for",folders[i],"% of points")) + 
          xlab("r") + ylab("PCF value")
        # grid.arrange(p1,p3,nrow = 2)
        ptot = p3 + annotation_custom(ggplotGrob(p1), xmin = 1750, xmax = 8000, 
                                      ymin = 2, ymax = 4.5)
        print(ptot)
        if (k == 1){
          example_plots_b[[i]] <- ptot
        }
        print(paste0(i," -- ",k))
      }
      dev.off()
    }
    
    ########## end of script 09_Comparison...E2.R from Eduardo ##########
  })
  
  return(list(plot_compare_exp_a, plot_compare_exp_b, example_plots_a, example_plots_b))
}












