temp_function <- function(x1, x2){
  return(x1+x2)
}

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

function_1 <- function(filepath_shp, file_shp, filepath_poly, file_poly, n_simulations_1,n_simulations_2, proportion){
  # for the functions to work, the files of the user need to be uploaded and temporarily saved. I need to make sure that stuff are deleted afterwards
  
  ## Load and prepare vector data
  # To avoid complications along the way, it is better to load a database that only contain
  # site IDs and coordinates.
  sites <-readOGR(dsn = filepath_shp, layer = file_shp)
  
  # Import the verctor file of the study region as a SpatialPolygons
  p <-readOGR(dsn = filepath_poly, layer = file_poly)
  
  # Convert the study region polygon to an Object Windon Class
  pwc <- as.owin(p)
  
  # Convert the groups of archaeological sites to a Spatial Point Pattern class and remove the 'marks'
  sites.um <- unmark(as.ppp(sites)) #Some functions only work with unmark sites
  
  # Combine the Spatial Point Pattern class groups of sites with the Object Window of the study polygon
  sp <- sites.um[pwc]
  
  ## Random sampling - Database with _x_% of the sites
  # Extract a 10% random sample of the 95 sites, i.e. 9.5 and round up to 10 sites, n = 85.
  # Then, reproduce the same analysis 100 times to determine the site structure variance.
  
  # initialize s_x_rl as an empty list
  s_x_rlf <- list()
  # initialize s_x_rlf.env as an empty list
  s_x_rlf.env <- list()
  
  set.seed(100) # For reproducibility
  for (i in 1:n_simulations_2) {
    # generate a new random sample each time through the loop
    s_x_rlf[[i]] <- sites[sample(nrow(sites), size = nrow(sites)-round(nrow(sites)*(1-proportion)), replace = FALSE, prob = NULL), ]
    s_x_f.uml <- unmark(as.ppp(s_x_rlf[[i]])) # apply unmark function
    s_x_f.spl <- s_x_f.uml[pwc] # create Object Window polygon
    # run envelope function with a Pair Correlation Function for Inhomogeneous process function and store result
    s_x_rlf.env[[i]] <- envelope(s_x_f.spl, fun=pcfinhom, nsim = n_simulations_1, divisor="d", correction = "iso") 
  }
  
  # Extract the significant values for each model
  # Create an empty list to store the points outside the envelope for each simulation
  outside_upper_x_f <- list()
  outside_lower_x_f <- list()
  
  # Loop through each simulation and extract the upper and lower envelope values
  for (i in 1:n_simulations_2) {
    upper_env_x_f <- s_x_rlf.env[[i]]$hi
    lower_env_x_f <- s_x_rlf.env[[i]]$lo
    
    # Subset the data to only include points outside the envelope
    outside_upper_x_f[[i]] <- subset(s_x_rlf.env[[i]]$obs, s_x_rlf.env[[i]]$obs > upper_env_x_f)
    outside_lower_x_f[[i]] <- subset(s_x_rlf.env[[i]]$obs, s_x_rlf.env[[i]]$obs < lower_env_x_f)
  }
  
  
  
  
  return(list(s_x_rlf, s_x_rlf.env, pwc, p))
}



big_processing_func <- function(filepath_shp, file_shp, filepath_poly, file_poly, nsim){
  
  sites <-readOGR(dsn = filepath_shp, layer = file_shp)
  sites.um <- unmark(as.ppp(sites))
  
  p <-readOGR(dsn = filepath_poly, layer = file_poly)
  pwc <- as.owin(p)
  
  sp_unmark <- sites.um[pwc]
  
  ### save(sp, file = "RData/01_SpObject.RData") ###
  
  sites_asppp <- as.ppp(sites)

  sp_mark <- sites_asppp[pwc]
  
  ### save(sp, file = "RData/01_SpObject_withmarks.RData") ###
  
  ########## end of scripts 02_...R. Output two RData objects #########
  
  
  
  
  
  ########## start of script 02_generateMC....R from Sébastien ##########
  
  # library(spatstat)
  # library(doParallel)
  
  #cl = makeCluster(10)
  cl = makeCluster(5)
  registerDoParallel(cl)
  
  #load(file = "RData/01_SpObject.RData")
  
  #cat(paste0("Debut : ",Sys.time(),"\n"),file = paste0("output/manRun",numero,".txt"))
  
  portions = seq(0.5,1.0,0.1)
  PO = length(portions)
  
  simulations_list <- list() # saving it in a list is new, involving also the .combine = 'c'
  
  for (j in 1:PO){
    set.seed(j)
    indToRemove = sample(x = 1:sp_unmark$n, size = floor((1-portions[j])*sp_unmark$n), replace = FALSE)
    indToKeep = setdiff(1:sp_unmark$n, indToRemove)
    spModified = sp_unmark[indToKeep]
    #print(paste0('Debut ronde ',j," (",as.character(Sys.time()),")\n"))
    simulations_list[[j]] <- foreach(i = 1:5, .packages = c("spatstat")) %dopar% {
      MCenvs = envelope(spModified, fun=pcfinhom, nsim = nsim, divisor="d", correction = "iso", savepatterns = FALSE, savefuns = TRUE)
      simulations = attr(MCenvs, "simfuns")
      #save(simulations, file = paste0("RData/0",100*(portions[j]),"/02-_MCenveloppe1E4-",i,".RData"))
    }
    #print(paste0('Ronde ',j," terminee (",as.character(Sys.time()),")\n"))
  }
  
  #cat(paste0("Fin : ",Sys.time(),"\n"),file = paste0("output/manRun",numero,".txt"), append = TRUE)
  
  stopCluster(cl)
  
  
  ########## end of script 02_...R (no parallel right now) ##########
  
  
  
  
  
  ########## start of script 04_computeMC....R from Eduardo ##########
  print("start of script 04_computeMC....R")
  
  folders = c("050","060","070","080","090","100")
  #folders = c("temp")
  #temp_list = c("50","60","70","80","90","100")
  levels = c(0.001,0.002,0.005,0.01,0.02,0.05,0.1)
  levels = c(levels,0.5,rev(1-levels))
  
  allQuantiles_list <- list()
  allQuantilesT_list <- list()
  
  for (f in 1:length(folders)){
    allSim <- do.call(cbind, lapply(simulations_list[[f]], function(x) as.matrix(x)[,1:nsim+1]))
    
    # allSim = c()
    # for (i in 1:10){
    #   #load(paste0("CECI/RData/",folders[f],"/02-_MCenveloppe1E4-",i,".RData"))
    #   #load(paste0("RData/02_MCenveloppe100000-",temp_list[i],".RData"))
    #   allSim = cbind(allSim,as.matrix(simulations)[,2:nsim+1])
    #   print(paste0(f," -- ",i))
    # }
    
    allQuantiles = apply(allSim, 1, quantile, probs = levels)
    colnames(allQuantiles) = as.matrix(simulations_list[[f]][[1]])[,1]
    rownames(allQuantiles) = levels
    allQuantilesT = as_tibble(as.data.frame.table(allQuantiles))
    colnames(allQuantilesT) = c("level","r","value")
    allQuantilesT$r = as.numeric(as.character(allQuantilesT$r))
    #save(allQuantiles, allQuantilesT, file = paste0("RData/",folders[f],"-quantiles.RData"))
    allQuantiles_list[[f]] <- allQuantiles
    allQuantilesT_list[[f]] <- allQuantilesT
    
    #pdf(paste0("plots/",folders[f],"-quantiles.pdf"))
    p = ggplot(allQuantilesT) +
      geom_hline(yintercept = 1, lty = "dotted") +   
      geom_line(aes(x = r, y = value, col = level)) + 
      scale_y_continuous(limits = c(0,9), breaks = 0:9) + 
      # geom_text(aes(x = r+1, y = value, col = level, label = level), data = subset(allQuantilesT, r == max(allQuantilesT$r))) + 
      ggtitle(paste0("MC quantiles for the robustness framework ",folders[f]," (", nsim*5, " MC scenarios)")) + ylab("PCF value") + 
      theme_light() + guides(col=guide_legend(nrow=3,byrow=TRUE)) + theme(legend.position = "bottom")
    print(p)
    #dev.off()
  }
  
  ########## end of script 04_computeMC....R from Eduardo ##########
  
  
  
  
  
  
  
  
  
  ########## start of script 03_computeMC....R from Sébastien ##########
  print("start of script 03_computeMC....R")
  #nbScMC = c(100,1000,10000,100000)
  nbScMC = c(2,5,10)
  NS = length(nbScMC)
  
  levels = c(0.005,0.01,0.05,0.1)
  levels = c(levels,0.5,rev(1-levels))
  
  # Apply as.matrix(), index specific columns (e.g., only the first column), and then cbind
  allSim <- do.call(cbind, lapply(simulations_list[[length(simulations_list)]], function(x) as.matrix(x)[,1:nsim+1]))

  # allSim = c()
  # for (i in 1:10){
  #   allSim = cbind(allSim,as.matrix(simulations)[,2:10001])
  #   print(paste0(" -- ",i))
  # }
  
  allQuantilesT = tibble(level = NA, r = NA, value = NA, nberSc = NA, repet = NA, .rows = 0)
  # ????? why a minus 1?????
  for (i in 1:(NS-1)){
    set.seed(2)
    indAvailable = 1:dim(allSim)[2]
    for (j in 1:10){
      indSampled = sample(indAvailable, size = nbScMC[i], replace = FALSE)
      allQuantiles = apply(allSim[,indSampled], 1, quantile, probs = levels)
      colnames(allQuantiles) = as.matrix(simulations_list[[length(simulations_list)]][[1]])[,1]
      rownames(allQuantiles) = levels
      ajout = as_tibble(as.data.frame.table(allQuantiles))
      colnames(ajout) = c("level","r","value")
      ajout = transform(ajout, nberSc = nbScMC[i], repet = j)
      allQuantilesT = rbind(allQuantilesT, ajout)
      indAvailable = setdiff(indAvailable, indSampled)
      print(paste0(i,"--",j))
    }
  }
  
  allQuantilesT$r = as.numeric(as.character(allQuantilesT$r))
  allQuantilesT$nberSc = as.character(allQuantilesT$nberSc)
  allQuantilesT$repet = as.character(allQuantilesT$repet)
  #save(allQuantilesT, file = paste0("RData/full-quantiles.RData"))
  
  
  essai = subset(allQuantilesT, level == 0.995 & repet == 1)
  
  ggplot(essai) + 
    geom_line(aes(x = r, y = value, col = nberSc))
  
  essai2a = subset(allQuantilesT, level == 0.995)
  essai2b = subset(allQuantilesT, level == 0.005)
  
  #load(file = "RData/01_SpObject.RData")
  pctOriginal = pcfinhom(sp_unmark)
  pctOriginalT = tibble(r = pctOriginal$r, value = pctOriginal$trans)
  
  #pdf("plots/03_quantileComparison01.pdf", height = 10, width = 7)
  ggplot(essai2a) + 
    geom_line(aes(x = r, y = value, col = repet)) +
    geom_line(aes(x = r, y = value, col = repet), data = essai2b) + 
    geom_line(aes(x = r, y = value), data = pctOriginalT, linewidth = 0.75) +
    coord_cartesian(ylim = c(0,3)) + 
    theme_light() + theme(legend.position = "none") + ggtitle(label = "MC enveloppe when varying the number of scenarios\nused to compute the 99.5% quantile", subtitle = "Each color curve corresponds to one of the 10 repetitions of the quantile computations") +  
    facet_wrap(~paste0("Number of simulations = ",nberSc),nrow = 3)
  #dev.off()
  
  ########## end of script 03_computeMC....R from Sébastien ##########
  
  
  
  
  
  
  
  ########## start of script 04_generateRobust...A.R from Sébastien ##########
  print("start of script 04_generateRobust...A.R")
  #load(file = "RData/01_SpObject.RData")
  
  folders = c("050","060","070","080","090","100")
  dataPortionsToKeep = seq(from = 0.5, to = 0.9, by = 0.1)
  
  PT = length(dataPortionsToKeep)
  #nbRobSc = 10000
  nbRobSc = 10
  
  pctOriginal = pcfinhom(sp_unmark)
  RV = length(pctOriginal$r)
  
  pctShockedData_a_list <- list()
  shockedSp_a_list <- list()
  
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
    pctShockedData_a_list[[i]] <- pctShockedData
    shockedSp_a_list[[i]] <- shockedSp
    #save(pctShockedData, shockedSp, file = paste0("RData/",folders[i],"-robustnessScenariosA.RData"))
  }
  #t2 = Sys.time()
  
  
  
  
  ########## end of script 04_generateRobust...A.R from Sébastien ##########
  
  
  
  
  
  
  
  
  
  ########## start of script 04_generateRobust...B.R from Sébastien ##########
  print("start of script 04_generateRobust...B.R")
  
  #load(file = "RData/01_SpObject_withmarks.RData")
  
  #folders = c("050","060","070","080","090","100")
  dataPortionsToKeep = seq(from = 0.5, to = 0.9, by = 0.1)
  
  PT = length(dataPortionsToKeep)
  # ???????? this was 10000 but why? Is it connected to number of simulations??????+
  #nbRobSc = 10000
  nbRobSc = 10
  
  pctOriginal = pcfinhom(sp_mark)
  RV = length(pctOriginal$r)
  
  pctShockedData_b_list <- list()
  shockedSp_b_list <- list()
  
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
    pctShockedData_b_list[[i]] <- pctShockedData
    shockedSp_b_list[[i]] <- shockedSp
    #save(pctShockedData, shockedSp, file = paste0("RData/",folders[i],"-robustnessScenariosB.RData"))
  }
  #t2 = Sys.time()
  
  ########## end of script 04_generateRobust...B.R from Sébastien ##########
  
  
  
  
  
  
  
  
  
  ########## start of script 07_assessRobust...A.R from Sébastien ##########
  print("start of script 07_assessRobust...A.R")
  
  #load(file = "RData/01_SpObject.RData")
  
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
    #load(paste0("CECI/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    #load(paste0("RData/",folders[i],"-quantiles.RData"))
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT_list[[i]], level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indFinite = which(is.finite(pctShockedData_a_list[[i]][,k]))
        allBinary[i,j,k] = any(pctShockedData_a_list[[i]][indFinite,k] > thisQuantile$value[indFinite])
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
  
  #pdf("plots/07_robustnessExperimentA01.pdf", width = 7, height = 5)
  p1 = ggplot(subset(results, level != 0.98)) + 
    geom_point(aes(x = keptPortion, y = value, col = level)) + 
    geom_line(aes(x = keptPortion, y = value, col = level)) + 
    scale_y_continuous(limits = c(0.49,1)) + 
    theme_light() + guides(col=guide_legend(nrow=1, title = "quantile level"))+
    xlab("% of sites that are kept in each robustness scenario") + ylab("% of robustness scenarios\nin which the conclusion is similar") + 
    theme(legend.position = "bottom") + ggtitle("First comparison tool")
  print(p1)
  #dev.off()
  
  # 2. Distribution des milieux ----
  allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    # load(paste0("../R/CECI/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    # load(paste0("RData/",folders[i],"-quantiles.RData"))
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT_list[[i]], level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indExceed = which(pctShockedData_a_list[[i]][,k] > thisQuantile$value)
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
    thisQuantile = subset(allQuantilesT_list[[length(allQuantilesT_list)]], level == levelsHigh[j])
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
  
  #pdf("plots/07_robustnessExperimentA02.pdf", width = 7, height = 5)
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
  #dev.off()
  
  
  # Et on met tout ensemble : 
  #pdf("plots/07_robustnessExperimentA99.pdf", width = 7, height = 10)
  grid.arrange(p1,p2)
  #dev.off()
  
  ########## end of script 07_assessRobust...A.R from Sébastien ##########
  
  
  
  
  
  
  ########## start of script 07_assessRobust...B.R from Sébastien ##########
  print("start of script 07_assessRobust...B.R")
  
  
  #load(file = "RData/01_SpObject.RData")
  
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
    #load(paste0("CECI/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    #load(paste0("RData/",folders[i],"-quantiles.RData"))
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT_list[[i]], level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indFinite = which(is.finite(pctShockedData_b_list[[i]][,k]))
        allBinary[i,j,k] = any(pctShockedData_b_list[[i]][indFinite,k] > thisQuantile$value[indFinite])
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
  
  #pdf("plots/07_robustnessExperimentB01.pdf", width = 7, height = 5)
  p1 = ggplot(subset(results, level != 0.98)) + 
    geom_point(aes(x = keptPortion, y = value, col = level)) + 
    geom_line(aes(x = keptPortion, y = value, col = level)) + 
    scale_y_continuous(limits = c(0.49,1)) + 
    theme_light() + guides(col=guide_legend(nrow=1, title = "quantile level"))+
    xlab("% of sites that are kept in each robustness scenario") + ylab("% of robustness scenarios\nin which the conclusion is similar") + 
    theme(legend.position = "bottom") + ggtitle("First comparison tool")
  print(p1)
  #dev.off()
  
  # 2. Distribution des milieux ----
  allFirstIntervals = array(list(), c(PT,LE,nbRobSc))
  
  for (i in 1:PT){
    #load(paste0("../R/CECI/RData/",folders[i],"-robustnessScenarios",experiment,".RData"))
    #load(paste0("RData/",folders[i],"-quantiles.RData"))
    for (j in 1:LE){
      thisQuantile = subset(allQuantilesT_list[[i]], level == levelsHigh[j])
      for (k in 1:nbRobSc){
        indExceed = which(pctShockedData_b_list[[i]][,k] > thisQuantile$value)
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
    thisQuantile = subset(allQuantilesT_list[[length(allQuantilesT_list)]], level == levelsHigh[j])
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
  
  #pdf("plots/07_robustnessExperimentB02.pdf", width = 7, height = 5)
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
  #dev.off()
  
  
  # Et on met tout ensemble : 
  #pdf("plots/07_robustnessExperimentB99.pdf", width = 7, height = 10)
  grid.arrange(p1,p2)
  #dev.off()
  
  
  
  ########## end of script 07_assessRobust...B.R from Sébastien ##########
  
  
  
  
  
  
  ########## start of script 08_Otherplots.R from Eduardo ########
  
  
  
  ########## end of script 08_Otherplots.R from Eduardo ########
  
  
  
  
  
  
  
  ########## start of script 09_Comparison...E1.R from Eduardo ##########
  print("start of script 09_Comparison...E1.R")
  
  folders = c("050","060","070","080","090","100")
  dataPortionsToKeep = seq(0.5, 0.9,0.1)
  PT = length(dataPortionsToKeep)
  
  #load(file = "RData/01_SpObject.RData")
  pctOriginal = pcfinhom(sp_unmark)
  
  #load("RData/100-quantiles.RData")
  levelsLow = c(0.005,0.01,0.02,0.05,0.1)
  LE = length(levelsLow)
  levelsHigh = 1-levelsLow
  
  mycolors = (heat.colors(LE))
  
  for (i in 1:PT){
    #load(paste0("RData/",folders[i],"-robustnessScenariosA.RData"))
    #load(paste0("RData/",folders[i],"-quantiles.RData"))
    
    explications = tibble(x = rep(-50, 6), 
                          y = allQuantilesT_list[[i]]$value[allQuantilesT_list[[i]]$level %in% as.character(c(0.5,0.9,0.95,0.98,0.99,0.995)) & allQuantilesT_list[[i]]$r == 0],
                          lab = c('50%','90%','95%','98%','99%','99.5%'))
    
    p0 = ggplot(allQuantilesT_list[[i]]) + 
      geom_blank() + theme_light()
    for (j in 1:LE){
      thisAdd = tibble(x = pctOriginal$r, ymin = subset(allQuantilesT_list[[i]], level == levelsLow[j])$value,
                       ymax = subset(allQuantilesT_list[[i]], level == levelsHigh[j])$value)
      p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                            fill = mycolors[j], alpha = 0.5, lty =  "dotted", 
                            col=alpha("black", 0.2))
    }
    p0 = p0 + geom_line(aes(x=r, y = value), data =  subset(allQuantilesT_list[[i]], level == 0.5), col = "blue", alpha = 0.5, lty = "dotted")
    
    #pdf(paste0("plots/11_examplePlots",folders[i],"-01.pdf"), width = 7, height = 5)
    for (k in 1:nbRobSc){
      thisSp = shockedSp_a_list[[i]][[k]]
      p1 = ggplot() +
        layer_spatial(data = sp_unmark, col = "red", fill="antiquewhite") +
        layer_spatial(data = thisSp, fill = NA, col = "black") +
        theme_void()
      thisData = tibble(y = pctShockedData_a_list[[i]][,k], x = pctOriginal$r)
      p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
        scale_y_continuous(limits = c(0,4.5)) + 
        geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
        ggtitle(paste("Scenario #",k,"for",folders[i],"% of points")) + 
        xlab("r") + ylab("PCF value")
      # grid.arrange(p1,p3,nrow = 2)
      ptot = p3 + annotation_custom(ggplotGrob(p1), xmin = 1750, xmax = 8000, 
                                    ymin = 2, ymax = 4.5)
      print(ptot)
      print(paste0(i," -- ",k))
    }
    #dev.off()
  }
  
  ########## end of script 09_Comparison...E1.R from Eduardo ##########
  
  
  
  
  
  
  
  ########## start of script 09_Comparison...E2.R from Eduardo ##########
  print("start of script 09_Comparison...E2.R")
  
  folders = c("050","060","070","080","090","100")
  dataPortionsToKeep = seq(0.5, 0.9,0.1)
  PT = length(dataPortionsToKeep)
  
  #load(file = "RData/01_SpObject.RData")
  pctOriginal = pcfinhom(sp_unmark)
  
  #load("RData/100-quantiles.RData")
  levelsLow = c(0.005,0.01,0.02,0.05,0.1)
  LE = length(levelsLow)
  levelsHigh = 1-levelsLow
  
  mycolors = (heat.colors(LE))
  
  for (i in 1:PT){
    #load(paste0("RData/",folders[i],"-robustnessScenariosB.RData"))
    #load(paste0("RData/",folders[i],"-quantiles.RData"))
    
    explications = tibble(x = rep(-50, 6), 
                          y = allQuantilesT_list[[i]]$value[allQuantilesT_list[[i]]$level %in% as.character(c(0.5,0.9,0.95,0.98,0.99,0.995)) & allQuantilesT_list[[i]]$r == 0],
                          lab = c('50%','90%','95%','98%','99%','99.5%'))
    
    p0 = ggplot(allQuantilesT_list[[i]]) + 
      geom_blank() + theme_light()
    for (j in 1:LE){
      thisAdd = tibble(x = pctOriginal$r, ymin = subset(allQuantilesT_list[[i]], level == levelsLow[j])$value,
                       ymax = subset(allQuantilesT_list[[i]], level == levelsHigh[j])$value)
      p0 = p0 + geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), data = thisAdd, 
                            fill = mycolors[j], alpha = 0.5, lty =  "dotted", 
                            col=alpha("black", 0.2))
    }
    p0 = p0 + geom_line(aes(x=r, y = value), data =  subset(allQuantilesT_list[[i]], level == 0.5), col = "blue", alpha = 0.5, lty = "dotted")
    
    #pdf(paste0("plots/11E2_examplePlots",folders[i],"-01.pdf"), width = 7, height = 5)
    for (k in 1:nbRobSc){
      thisSp = shockedSp_b_list[[i]][[k]]
      p1 = ggplot() +
        layer_spatial(data = sp_unmark, col = "red", fill="antiquewhite") +
        layer_spatial(data = thisSp, fill = NA, col = "black") +
        theme_void()
      thisData = tibble(y = pctShockedData_b_list[[i]][,k], x = pctOriginal$r)
      p3 = p0 + geom_line(aes(x=x, y=y), data = thisData, linewidth = 0.75) + 
        scale_y_continuous(limits = c(0,4.5)) + 
        geom_text(data = explications, aes(x=x, y=y, label = lab), size = 2.5, hjust = 1) + 
        ggtitle(paste("Scenario #",k,"for",folders[i],"% of points")) + 
        xlab("r") + ylab("PCF value")
      # grid.arrange(p1,p3,nrow = 2)
      ptot = p3 + annotation_custom(ggplotGrob(p1), xmin = 1750, xmax = 8000, 
                                    ymin = 2, ymax = 4.5)
      print(ptot)
      print(paste0(i," -- ",k))
    }
    #dev.off()
  }
  
  ########## end of script 09_Comparison...E2.R from Eduardo ##########
  
  
  return()
}












