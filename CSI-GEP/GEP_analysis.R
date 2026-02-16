args = (commandArgs(trailingOnly = TRUE))

library(jaccard)
library(dplyr)
library(igraph)
library(vroom)
library(elbow)
library(ggplot2)

set.seed(1234)
calcSplitStats <- function(split0, split1, jval){
  print("calcSplitStats")
  # for each row of split0, compare it to every row of split1
  jaccardInd <- numeric()
  jaccardIndPs <- numeric()
  
  for(i in 1:ncol(split0)){
    for(j in 1:ncol(split1)){
      GEP_split0 = as.numeric(rank(-split0[, i]) < jval)
      GEP_split1 = as.numeric(rank(-split1[, j]) < jval)
      jaccardTest <- jaccard.test.mca(GEP_split0, GEP_split1)
      jaccardInd <- c(jaccardInd, jaccardTest$statistics)
      jaccardIndPs <- c(jaccardIndPs, jaccardTest$pvalue)
    }
  }
  
  # Calculate FWERs
  fwerJaccard <- p.adjust(jaccardIndPs, method="bonferroni")
  fwerJaccard <- t(matrix(fwerJaccard, nrow = dim(split1)[2], ncol = dim(split0)[2]))
  
  # Calculate number of significant findings.
  sigFwerJaccard <- which(fwerJaccard < 0.05, arr.ind = T)
  
  #Calculate community number
  if(nrow(sigFwerJaccard) > 0){
    community.df <- as.data.frame(sigFwerJaccard)
    community.df$row <- paste0("Split0_GEP", community.df$row)
    community.df$col <- paste0("Split1_GEP", community.df$col)
    community.graph <- graph_from_data_frame(community.df, directed = F)
    ceb <- cluster_edge_betweenness(community.graph)
    community_number <- length(communities(ceb))
  }else{
    community_number <- 0
  }
  return(list(jaccardIndPs=jaccardIndPs, jaccardInd=jaccardInd, sigFwerJaccard=sigFwerJaccard, jaccard.com=community_number))
}
#############################################################################################3
#Interpret the Jaccard Test results for individual dataset benchmarking study
#Obtain the Rank and Jaccard length that maximize the number of communitites in network graph
#Once these paramaters are identified, extract GEPs from cNMF output (in next code block)
###############################################################################################
home.out.dir <- args[1] 
dataset.name <- args[2]

output.dir <- args[3]
ranks = scan(text=args[4], sep = ",")
print(ranks)

check_for_rescue = as.logical(as.numeric(args[5]))

jaccard.val <- c(seq(10,100,10))
    
jaccardtest_results <- list()

for(k in ranks){
    #These are the results of the calcSplitStats function, which calculates the FWER-corrected jaccard significance test 
    splitstats <- readRDS(file = paste0(home.out.dir, "results/jaccardtest_results_K", k, ".RDS"))
    jaccardtest_results[[k]] <- splitstats
}

jaccard.val <- as.character(jaccard.val)

com.df <- NULL

for(x in ranks){
    u <- NULL
    #print(paste0("Rank: ", x))
    for(y in jaccard.val){
        u <- c(u, jaccardtest_results[[x]][[y]][["jaccard.com"]]) 
    }
    
    u <- cbind(u, x)
    u <- data.frame(u, jaccard.val)
    com.df <- rbind(com.df, u)
}   

colnames(com.df) <- c("Community_Number", "Rank", "Jaccard_Length")
saveRDS(com.df, file = paste0(output.dir, dataset.name, "_jaccardtest_communities.RDS"))


#Insert plotting function here
#This block of code evaluates which curve has the best inflection point
jl.eval <- data.frame()

for(x in jaccard.val){
    x <- as.numeric(x)  
    jl.df <- com.df %>% filter(Jaccard_Length == x)
    #Run a loess regression on JL data to smooth the curve
    sl <- jl.df %>% ggplot(aes(Rank, Community_Number)) + geom_point() + geom_smooth(method = "loess", span = 0.3, method.args = list(degree = 1))
    gbuild <- ggplot_build(sl)
    elbow.df <- data.frame(X = gbuild$data[[2]]$x, Y = gbuild$data[[2]]$y)
    #Find inflection point
    ipoint <- elbow(data = elbow.df, plot = F)
    #Since the data was run on smoothed data, retreive the closest rank to the x intercept calculated by elbow package
    optimum.rank <- ranks[which.min(abs(ranks - ipoint$X_selected))]
    #Run linear regression on cost values (X values after inflection point) and determine slope
    ipoint.filtered <- ipoint$data %>% filter(X >= ipoint$X_selected)
    reg <- lm(ipoint.filtered$benefits ~ ipoint.filtered$X)
    jl.eval <- rbind(jl.eval, c(x, optimum.rank, reg$coefficients[2]))
}

colnames(jl.eval) <- c("JL", "Rank", "Slope")
#This will determine which JL curve gives the steepest slope for the cost values
max.slope.params <- jl.eval %>% filter(JL %in% jl.eval$JL[jl.eval$Slope == min(jl.eval$Slope)])
comnum <- com.df %>% filter(Jaccard_Length %in% max.slope.params$JL & Rank %in% max.slope.params$Rank)


params <- as.data.frame(cbind(comnum$Community_Number, max.slope.params$Rank, max.slope.params$JL))
colnames(params) <- c("Community_Number", "Rank", "Jaccard_Length")

writeLines(paste0(names(params)[1], ": ", as.character(params$Community_Number)), 
    file(paste0(output.dir, "output_summary.txt")))
write(paste0(names(params)[2], ": ", as.character(params$Rank)),
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)
write(paste0(names(params)[3], ": ", as.character(params$Jaccard_Length)),
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)

if(!(check_for_rescue) | (nrow(params)>1)){
    rescue = FALSE
    params = params[1,]
} else{
    rescue = TRUE
    n_1 = params$Community_Number
}


pdf(paste0(output.dir, dataset.name, "_Curves.pdf"))
ggplot(com.df, aes(x=Rank,y=Community_Number, group = Jaccard_Length, color = Jaccard_Length)) +
stat_smooth(method="loess", span=0.2, se=TRUE, aes(fill=Jaccard_Length), alpha=0.3) +
theme_linedraw() + geom_point(aes(x = params$Rank, y = params$Community_Number), color = "black") +
ggtitle(paste0("Community curves of GPU-based torchnmf")) + theme(plot.title = element_text(size=10))
dev.off()


#Get the significant matches(indeces) from the Jaccard test
sig.labels <- as.data.frame(jaccardtest_results[[params$Rank]][[as.character(params$Jaccard_Length)]]$sigFwerJaccard)
sig.labels$row <- paste0("Split0_GEP", sig.labels$row)
sig.labels$col <- paste0("Split1_GEP", sig.labels$col)
    

#####second round of programs rescuing
if(rescue){
    k <- params$Rank
    jval = params$Jaccard_Length                                                          

    #Read in NMF gene spectra scores  
    print(paste0("Rank: ", k))
    split0 <- t(vroom(file = paste0(home.out.dir, "cNMF_split0/cNMF_split0.gene_spectra_score.k_", k, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))
    split1 <- t(vroom(file = paste0(home.out.dir, "cNMF_split1/cNMF_split1.gene_spectra_score.k_", k, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))
    print(paste0("split0 Filename: cNMF_split0.gene_spectra_score.k_", k, ".dt_0_1.txt"))
    print(paste0("split1 Filename: cNMF_split1.gene_spectra_score.k_", k, ".dt_0_1.txt"))

    # remove columns already forming the network
    row_ind_keep = setdiff(1:k, unique(as.data.frame(jaccardtest_results[[k]][[as.character(jval)]]$sigFwerJaccard)$row))
    col_ind_keep = setdiff(1:k, unique(as.data.frame(jaccardtest_results[[k]][[as.character(jval)]]$sigFwerJaccard)$col))

    if (length(row_ind_keep)>0 & length(col_ind_keep)>0){
        split0_rescue = data.frame(split0[,row_ind_keep])
        split1_rescue = data.frame(split1[,col_ind_keep])
        #Compute jaccard index for fixed JL 
        new_jval = min(as.numeric(jaccard.val)[as.numeric(jaccard.val)>jval])

        print(paste0("Jaccard Length: ", new_jval))
        splitstatlist <- list()
        splitstatlist <- calcSplitStats(split0_rescue, split1_rescue, new_jval) 

        rescued.labels = as.data.frame(splitstatlist$sigFwerJaccard)

        if (nrow(rescued.labels) > 0){
            rescued.labels$row <- paste0("Split0_GEP", row_ind_keep[rescued.labels$row])
            rescued.labels$col <- paste0("Split1_GEP", col_ind_keep[rescued.labels$col])
            sig.labels = rbind(sig.labels, rescued.labels)
            } else{
                rescue = FALSE
            }
    } else {
        rescue = FALSE
    }
}

write(paste0("Rescued more programs: ", rescue), 
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)

saveRDS(sig.labels, file = paste0(output.dir, dataset.name,"_labels.RDS"))

###############################################################################
#Use network (igraph) results to obtain mean spectra score for each community 
##############################################################################3

dataset.com.graphs = list()
community_geps = list()
print(paste0("Calculating ", dataset.name, " averaged community-wide GEPs"))
#Build network graph from splits
dat.graph <- graph_from_data_frame(sig.labels, directed = F)
ceb <- cluster_edge_betweenness(dat.graph)
membership <- sort(membership(ceb))

dataset.com.graphs <- list(Graph = dat.graph, CEB = ceb)

#Obtain the labels for the community membership

gep.index <- list()
x <- numeric()
y <- numeric()

for(i in 1:max(membership)){
    z <- names(membership[membership == i])
    for(j in z){
        if(grepl("Split0", j)){
            temp <- as.numeric(substr(j,11,nchar(j)))
            x <- c(x, temp)
            }else{
                temp <- as.numeric(substr(j,11,nchar(j)))
                y <- c(temp, y)
            }
        }
        gep.index[[i]] <- list("Split0" = x, "Split1" = y)
        x <- NULL
        y <- NULL
}

saveRDS(gep.index, file = paste0(output.dir, dataset.name, "_gepindex.RDS"))
params$Community_Number = length(gep.index)
saveRDS(params, file = paste0(output.dir, dataset.name,"_optimized_parameters.RDS"))

#cNMF outpu

K = params$Rank
max.com.num = length(gep.index)


split0 <- t(vroom(file = paste0(home.out.dir, "cNMF_split0/cNMF_split0.gene_spectra_score.k_", K, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))
split1 <- t(vroom(file = paste0(home.out.dir, "cNMF_split1/cNMF_split1.gene_spectra_score.k_", K, ".dt_0_1.txt"), col_select = c(-...1), show_col_types = FALSE))

community_vectors <- list()

#Average Z-score ranked GEP for each community
for(i in 1:length(gep.index)){
    
    #These are the node names in each community
    #These names correspond to columns in the cNMF output
    s0 <- unlist(gep.index[i][[1]][1])
    s1 <- unlist(gep.index[i][[1]][2])

    if(length(s0) == 0 | length(s1) == 0){
        if(length(s0) > length(s1)){
            community_vectors[[i]] <- split0[, s0]
        }else{
            community_vectors[[i]] <- split1[, s1]
        }
    }else{
    
    #Compute the average score for the community
    mean.df.temp <- NULL
    mean.df.temp2 <- NULL
    if(length(s0) == 1){
        #print("Length Split0 = 1")
        if(length(s1) == 1){
            #print("Length Split1 = 1")
            community_vectors[[i]] <- rowMeans(cbind(split0[, s0], split1[, s1]), na.rm = T)
        }else{
            #print("Length Split1 > 1")
            for(j in s1){
                mean.df.temp <- cbind(mean.df.temp, split1[, j])
            }
            mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
            community_vectors[[i]] <- rowMeans(cbind(split0[, s0], mean.df.temp), na.rm = T)
        }
    }else if(length(s1) == 1){
        #print("Length Split0 > 1")
        #print("Length Split1 = 1")
        for(j in s0){
            mean.df.temp <- cbind(mean.df.temp, split0[, j])
        }
        mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
        community_vectors[[i]] <- rowMeans(cbind(mean.df.temp, split1[, s1]), na.rm = T)
    }else{
        #print("Length Split1 > 1")
        for(j in s0){
            mean.df.temp <- cbind(mean.df.temp, split0[, j])
        }
        for(k in s1){
            mean.df.temp2 <- cbind(mean.df.temp2, split1[, k])
        }

    mean.df.temp <- rowMeans(mean.df.temp, na.rm = T)
    mean.df.temp2 <- rowMeans(mean.df.temp2, na.rm = T)
    community_vectors[[i]] <- rowMeans(cbind(mean.df.temp, mean.df.temp2), na.rm = T)
    }
    }
}

community_vector_mat <- matrix(as.numeric(unlist(community_vectors)), ncol = max.com.num, nrow = nrow(split0))
rownames(community_vector_mat) <- rownames(split0)
community_vector_mat_filtered <- community_vector_mat[complete.cases(community_vector_mat), ]


saveRDS(community_vector_mat_filtered, file = paste0(output.dir, dataset.name, "_community_geps.RDS"))  

if (rescue){
    write(paste0("Number of programs rescued: ", max.com.num-n_1), 
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)
}

write(paste0("Automatically selected rank k = ",K), 
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)

write("CSIGEP finished!", 
    file = paste0(output.dir, "output_summary.txt"), append = TRUE)


