
#' Function doing the actual analysis before calling the plotting functions.
#'
#' @title run() : Invokes a routine SPACEmapX analysis to Infer CNV changes given a matrix of RNASeq counts.
#'
#' @param SPACEmapX_obj An SPACEmapX object populated with raw count data
#'
#' @param cutoff Cut-off for the min average read counts per gene among reference cells. (default: 1)
#'
#' @param min_cells_per_gene minimum number of reference cells requiring expression measurements to include the corresponding gene.
#'                           default: 3
#'
#' @param out_dir path to directory to deposit outputs (default: NULL, required to provide non NULL)
#'
#' ## Smoothing params
#' @param window_length Length of the window for the moving average
#'                          (smoothing). Should be an odd integer. (default: 101)#'
#'
#' @param smooth_method  Method to use for smoothing: c(runmeans,pyramidinal,coordinates)  default: pyramidinal
#'
#' #####
#'
#' @param num_ref_groups The number of reference groups or a list of
#'                           indices for each group of reference indices in
#'                           relation to reference_obs. (default: NULL)
#'
#' @param ref_subtract_use_mean_bounds   Determine means separately for each ref group, then remove intensities within bounds of means (default: TRUE)
#'                                       Otherwise, uses mean of the means across groups.
#'
#' #############################
#'
#' @param cluster_by_groups   If observations are defined according to groups (ie. patients), each group
#'                            of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
#'
#' @param cluster_references Whether to cluster references within their annotations or not. (dendrogram not displayed)
#'                             (default: TRUE)
#'
#'
#' @param k_obs_groups Number of groups in which to break the observations. (default: 1)
#'
#'
#'
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' default("ward.D2")
#'
#' @param max_centered_threshold The maximum value a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value. (default: 3),
#'                               can set to a numeric value or "auto" to bound by the mean bounds across cells.
#'                               Set to NA to turn off.
#'
#' @param scale_data  perform Z-scaling of logtransformed data (default: FALSE).  This may be turned on if you have
#'                    very different kinds of data for your normal and tumor samples. For example, you need to use GTEx
#'                    representative normal expression profiles rather than being able to leverage normal single cell data
#'                    that goes with your experiment.
#'
#' #########################################################################
#' ## Downstream Analyses (HMM or non-DE-masking) based on tumor subclusters
#'
#' @param HMM  when set to True, runs HMM to predict CNV level (default: FALSE)
#'
#' @param HMM_transition_prob   transition probability in HMM (default: 1e-6)
#'
#'
#' @param HMM_report_by   cell, consensus, subcluster (default: subcluster)  Note, reporting is performed entirely separately from the HMM prediction.  So, you can predict on subclusters, but get per-cell level reporting (more voluminous output).
#'
#'
#' @param HMM_type  HMM model type. Options: (i6 or i3):
#'                          i6: SPACEmapX 6-state model (0, 0.5, 1, 1.5, 2, >2) where state emissions are calibrated based on simulated CNV levels.
#'                          i3: SPACEmapX 3-state model (del, neutral, amp) configured based on normal cells and HMM_i3_pval
#'
#' @param HMM_i3_pval     p-value for HMM i3 state overlap (default: 0.05)
#'
#'
#' @param HMM_i3_use_KS   boolean: use the KS test statistic to estimate mean of amp/del distributions (ala HoneyBadger). (default=TRUE)
#'
#'
#' ## Filtering low-conf HMM preds via BayesNet P(Normal)
#'
#' @param BayesMaxPNormal  maximum P(Normal) allowed for a CNV prediction according to BayesNet. (default=0.5, note zero turns it off)
#' 
#' @param reassignCNVs (boolean) Given the CNV associated probability of belonging to each possible state, 
#'                      reassign the state assignments made by the HMM to the state that has the highest probability. (default: TRUE)
#'
#' ######################
#' ## Tumor subclustering
#'
#' @param analysis_mode options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions.
#'                                                          default: samples (fastest, but subclusters is ideal)
#'
#' @param tumor_subcluster_partition_method  method for defining tumor subclusters. Options('leiden', 'random_trees', 'qnorm')
#'                                           leiden: Runs a nearest neighbor search, where communities are then partitionned with the Leiden algorithm.
#'                                           random_trees: Slow, uses permutation statistics w/ tree construction.
#'                                           qnorm: defines tree height based on the quantile defined by the tumor_subcluster_pval
#'
#' @param tumor_subcluster_pval max p-value for defining a significant tumor subcluster (default: 0.1)
#'
#' @param k_nn number k of nearest neighbors to search for when using the Leiden partition method for subclustering (default: 20)
#'
#' @param leiden_method Method used to generate the graph on which the Leiden algorithm is applied, one of "PCA" or "simple". (default: "PCA")
#'
#' @param leiden_function Whether to use the Constant Potts Model (CPM) or modularity in igraph. Must be either "CPM" or "modularity". (default: "CPM")
#'
#' @param leiden_resolution resolution parameter for the Leiden algorithm using the CPM quality score (default: 0.05)
#'
#' @param leiden_method_per_chr Method used to generate the graph on which the Leiden algorithm is applied for the per chromosome subclustering, one of "PCA" or "simple". (default: "simple")
#'
#' @param leiden_function_per_chr Whether to use the Constant Potts Model (CPM) or modularity in igraph for the per chromosome subclustering. Must be either "CPM" or "modularity". (default: "modularity")
#'
#' @param leiden_resolution_per_chr resolution parameter for the Leiden algorithm for the per chromosome subclustering (default: 1)
#'
#' @param per_chr_hmm_subclusters Run subclustering per chromosome over all cells combined to run the HMM on those subclusters instead. Only applicable when using Leiden subclustering.
#'                                This should provide enough definition in the predictions while avoiding subclusters that are too small thus providing less evidence to work with. (default: FALSE)
#'
#' @param per_chr_hmm_subclusters_references Whether the per chromosome subclustering should also be done on references, which should not have as much variation as observations. (default = FALSE)
#'
#' 
#' @param z_score_filter Z-score used as a treshold to filter genes used for subclustering. 
#'                       Applied based on reference genes to automatically ignore genes with high expression variability such as MHC genes. (default: 0.8)
#'
#'
#' #############################
#' ## de-noising parameters ####
#'
#' @param denoise       If True, turns on denoising according to options below
#'
#' @param noise_filter  Values +- from the reference cell mean will be set to zero (whitening effect)
#'                      default(NA, instead will use sd_amplifier below.
#'
#' @param sd_amplifier  Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier
#'                      default: 1.5
#'
#' @param noise_logistic use the noise_filter or sd_amplifier based threshold (whichever is invoked) as the midpoint in a
#'                       logistic model for downscaling values close to the mean. (default: FALSE)
#'
#'
#' ##################
#' ## Outlier pruning
#'
#' @param outlier_method_bound Method to use for bounding outlier values. (default: "average_bound")
#'                             Will preferentially use outlier_lower_bounda and outlier_upper_bound if set.
#' @param outlier_lower_bound  Outliers below this lower bound will be set to this value.
#' @param outlier_upper_bound  Outliers above this upper bound will be set to this value.
#'
#'
#' ##########################
#' ## Misc options
#'
#' @param final_scale_limits The scale limits for the final heatmap output by the run() method. Default "auto". Alt, c(low,high)
#'
#' @param final_center_val   Center value for final heatmap output by the run() method.
#'
#' @param debug If true, output debug level logging.
#'
#' @param num_threads (int) number of threads for parallel steps (default: 4)
#'
#' @param plot_steps If true, saves SPACEmapX objects and plots data at the intermediate steps.
#'
#' @param inspect_subclusters If true, plot subclusters as annotations after the subclustering step to easily see if the subclustering options are good. (default = TRUE)
#'
#' @param resume_mode  leverage pre-computed and stored SPACEmapX objects where possible. (default=TRUE)
#'
#' @param png_res Resolution for png output.
#'
#' @param no_plot   don't make any of the images. Instead, generate all non-image outputs as part of the run. (default: FALSE)
#'
#' @param no_prelim_plot  don't make the preliminary SPACEmapX image (default: FALSE)
#'
#' @param write_expr_matrix Whether to write text files with the content of matrices when generating plots (default: FALSE)
#'
#' @param write_phylo Whether to write newick strings of the dendrograms displayed on the left side of the heatmap to file (default: FALSE)
#'
#' @param output_format Output format for the figure. Choose between "png", "pdf" and NA. NA means to only write the text outputs without generating the figure itself. (default: "png")
#'
#' @param plot_chr_scale Whether to scale the chromosme width on the heatmap based on their actual size rather than just the number of expressed genes.
#'
#' @param chr_lengths A named list of chromsomes lengths to use when plot_chr_scale=TRUE, or else chromosome size is assumed to be the last chromosome's stop position + 10k bp
#'
#' @param useRaster Whether to use rasterization for drawing heatmap. Only disable if it produces an error as it is much faster than not using it. (default: TRUE)
#'
#' @param plot_probabilities option to plot posterior probabilities (default: TRUE)
#'
#' @param save_rds Whether to save the current step object results as an .rds file (default: TRUE)
#'
#' @param save_final_rds Whether to save the final object results as an .rds file (default: TRUE)
#'
#' @param diagnostics option to create diagnostic plots after running the Bayesian model (default: FALSE)
#'
#' #######################
#' ## Experimental options
#'
#' @param remove_genes_at_chr_ends experimental option: If true, removes the window_length/2 genes at both ends of the chromosome.
#'
#' @param prune_outliers  Define outliers loosely as those that exceed the mean boundaries among all cells.  These are set to the bounds.
#'
#' ## experimental opts involving DE analysis
#'
#' @param mask_nonDE_genes If true, sets genes not significantly differentially expressed between tumor/normal to
#'                          the mean value for the complete data set (default: 0.05)
#'
#' @param mask_nonDE_pval  p-value threshold for defining statistically significant DE genes between tumor/normal
#
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.'
#'
#' @param require_DE_all_normals If mask_nonDE_genes is set, those genes will be masked only if they are are found as DE according to test.use and mask_nonDE_pval in each of the comparisons to normal cells options: {"any", "most", "all"} (default: "any")
#'
#' other experimental opts
#'
#' @param sim_method   method for calibrating CNV levels in the i6 HMM (default: 'meanvar')
#'
#' @param sim_foreground  don't use... for debugging, developer option.
#'
#' @param hspike_aggregate_normals  instead of trying to model the different normal groupings individually, just merge them in the hspike.
#'
#' @param up_to_step run() only up to this exact step number (default: 100 >> 23 steps currently in the process)
#'
#' @return SPACEmapX_obj containing filtered and transformed data
#'
#' @export
#'
#' @examples
#' data(SPACEmapX_data_example)
#' data(SPACEmapX_annots_example)
#' data(SPACEmapX_genes_example)
#'
#' SPACEmapX_object_example <- SPACEmapX::CreateSPACEmapXObject(raw_counts_matrix=SPACEmapX_data_example, 
#'                                                           gene_order_file=SPACEmapX_genes_example,
#'                                                           annotations_file=SPACEmapX_annots_example,
#'                                                           ref_group_names=c("normal"))
#'
#' SPACEmapX_object_example <- SPACEmapX::run(SPACEmapX_object_example,
#'                                          cutoff=1,
#'                                          out_dir=tempfile(), 
#'                                          cluster_by_groups=TRUE, 
#'                                          denoise=TRUE,
#'                                          HMM=FALSE,
#'                                          num_threads=2,
#'                                          analysis_mode="samples",
#'                                          no_plot=TRUE)
#'


#' Subtracting the mean of the reference expr distributions from the observed cells.
#'
#' @title subtract_ref_expr_from_obs()
#'
#' @description Remove the average of the genes of the reference observations from all
#' observations' expression. In the case there are multiple reference groupings,
#' the averages are computed separately for each reference grouping, and the min|max
#' of the averages are subtracted from the observation expression levels. Any values within the range
#' of the min,max of the group are set to zero.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param inv_log mean values will be determined based on (2^x -1)
#'
#' @param use_bounds if multiple normal data sets are used, it takes the bounds of the means from each set for subtraction.
#'                   Alternatively, will use the mean( mean(normal) for each normal)    default: TRUE
#'
#' @return SPACEmapX_obj containing the reference subtracted values.
#'
#' @keywords internal
#' @noRd
#'

subtract_ref_expr_from_obs <- function(SPACEmapX_obj, inv_log=FALSE, use_bounds=TRUE) {
    ## r = genes, c = cells
    flog.info(sprintf("::subtract_ref_expr_from_obs:Start inv_log=%s, use_bounds=%s", inv_log, use_bounds))
    
    
    if (has_reference_cells(SPACEmapX_obj)) {
        ref_groups = SPACEmapX_obj@reference_grouped_cell_indices
        flog.info("subtracting mean(normal) per gene per cell across all data")
    } else {
        ref_groups = list('proxyNormal' = unlist(SPACEmapX_obj@observation_grouped_cell_indices))
        flog.info("-no reference cells specified... using mean of all cells as proxy")
    }
    
    ref_grp_gene_means <- .get_normal_gene_mean_bounds(SPACEmapX_obj@expr.data, ref_groups, inv_log=inv_log)
    
    SPACEmapX_obj@expr.data <- .subtract_expr(SPACEmapX_obj@expr.data, ref_grp_gene_means, use_bounds=use_bounds)
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- subtract_ref_expr_from_obs(SPACEmapX_obj@.hspike, inv_log=inv_log, use_bounds=use_bounds)
    }
    
    return(SPACEmapX_obj)
    
}

#' @keywords internal
#' @noRd
#'

.get_normal_gene_mean_bounds <- function(expr.data, ref_groups, inv_log=FALSE) {
    
    get_indiv_gene_group_means_bounds_fun <- function(x) {
        
        grp_means = c()
        
        if (inv_log) {
            grp_means <- vapply(ref_groups, function(ref_group) {
                log2(mean(2^x[ref_group] - 1) + 1)
            }, double(1))
        }
        else {
            grp_means <- vapply(ref_groups, function(ref_group) {
                mean(x[ref_group])
            }, double(1))
        }

        names(grp_means) <- names(ref_groups)
        
        return(as.data.frame(t(data.frame(grp_means))))
    }
    
    gene_ref_grp_means <- do.call(rbind, apply(expr.data, 1, get_indiv_gene_group_means_bounds_fun))
    
    rownames(gene_ref_grp_means) <- rownames(expr.data)
    
    return(gene_ref_grp_means)
}


#' @keywords internal
#' @noRd
#'

.subtract_expr <- function(expr_matrix, ref_grp_gene_means, use_bounds=FALSE) {
    
    my.rownames = rownames(expr_matrix)
    my.colnames = colnames(expr_matrix)
    
    flog.info(sprintf("-subtracting expr per gene, use_bounds=%s", use_bounds))
    
    subtract_normal_expr_fun <- function(row_idx) {
        
        gene_means <- as.numeric(ref_grp_gene_means[row_idx, , drop=TRUE])
        
        gene_means_mean <- mean(gene_means)
        
        x <- as.numeric(expr_matrix[row_idx, , drop=TRUE])
        
        row_init = rep(0, length(x))
        
        if (use_bounds) {
            
            grp_min = min(gene_means)
            grp_max = max(gene_means)
            
            above_max = which(x>grp_max)
            below_min = which(x<grp_min)
            
            row_init[above_max] <- x[above_max] - grp_max
            row_init[below_min] <- x[below_min] - grp_min
            ## note, in-between values are left at zero!
            
        } else {
            
            row_init <- x - gene_means_mean
            
        }
        
        return(row_init)
    }
    
    subtr_data <- do.call(rbind, lapply(seq_len(nrow(expr_matrix)), subtract_normal_expr_fun))
    rownames(subtr_data) <- my.rownames
    colnames(subtr_data) <- my.colnames
    
    return(subtr_data)
    
}


#' @description Helper function allowing greater control over the steps in a color palette.
#'              Source: http://menugget.blogspot.com/2011/11/define-color-steps-for-
#'              colorramppalette.html#more
#'
#' @title Helper function allowing greater control over the steps in a color palette.
# 
#' @param steps Vector of colors to change use in the palette
#' @param between Steps where gradients change
#' @param ... Additional arguments of colorRampPalette
#
#' @return Color palette
#'
#' @export
#'
#' @examples
#' color.palette(c("darkblue", "white", "darkred"),
#'               c(2, 2))
#'

color.palette <- function(steps,
                          between=NULL, ...){
    
    if (is.null(between)){
        between <- rep(0, (length(steps) - 1))
    }
    if (length(between) != length(steps) - 1){
        stop("Must have one less \"between\" value than steps")
    }
    
    fill.steps <- cumsum(rep(1, length(steps)) + c(0, between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[, fill.steps] <- col2rgb(steps)
    for(i in which(between > 0)){
        col.start <- RGB[, fill.steps[i]]
        col.end <- RGB[, fill.steps[i + 1]]
        for(j in seq(3)){
            vals <- seq(col.start[j],
                        col.end[j],
                        length.out=between[i] + 2)
            vals <- vals[2:(2 + between[i] - 1)]
            RGB[j, (fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue=255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Create a sepList for the heatmap.3 plotting function given integer vectors
# of rows and columns where seperation should take place.
# The expected input to the heatmap function is a list of 2 lists.
# The first list are column based rectangles, and the second row.
# To define a rectangle the index of the row or column where the line of the rectagle
# should be placed is done with a vector of integers, left, bottom, right and top line.
# Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
#
# Args:
# row_count Total number of rows
# col_count Total number of columns
# row_seps Vector of integers indices for row breaks
# col_seps Vector of integer indices for column breaks
#
# Returns
# List of lists of vectors

#' @keywords internal
#' @noRd
#'

create_sep_list <- function(row_count,
                            col_count,
                            row_seps=NULL,
                            col_seps=NULL){
    sepList <- list()
                                        # Heatmap.3 wants a list of boxes for seperating columns
                                        # Column data
    if(!is.null(col_seps) &&
       !any(is.na(col_seps)) &&
       (length(col_seps)>0) &&
       col_count > 0){
        colList <- list()
        for(sep in seq_along(col_seps)){
            colList[[sep]] <- c(col_seps[sep],0,col_seps[sep],row_count)
        }
        sepList[[1]] <- colList
    } else {
        sepList[[1]] <- list()
        sepList[[1]][[1]] <- c(0,0,0,0)
    }
    
                                        # Row data
                                        # This is measured from bottom to top
                                        # So you have to adjust the values of the data
    row_seps <- row_count-row_seps
    if(!is.null(row_seps) &&
       !any(is.na(row_seps)) &&
       (length(row_seps)>0) &&
       row_count > 0){
        rowList <- list()
        for(sep in seq_along(row_seps)){
            rowList[[sep]] <- c(0,row_seps[sep],col_count,row_seps[sep])
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}


#' @title split_references()
#'
#' @description Split up reference observations in to k groups based on hierarchical clustering.
#'
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param num_groups (default: 2)
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

split_references <- function(SPACEmapX_obj,
                             num_groups=2,
                             hclust_method='complete') {
    
    
    flog.info(paste("::split_references:Start", sep=""))
    
    if ("sparseMatrix" %in% is(SPACEmapX_obj@expr.data)) {
        SPACEmapX_obj@expr.data = as.matrix(SPACEmapX_obj@expr.data)
    }

    ref_expr_matrix = SPACEmapX_obj@expr.data[ , get_reference_grouped_cell_indices(SPACEmapX_obj) ]
    
    hc <- hclust(parallelDist(t(ref_expr_matrix), threads=SPACEmapX.env$GLOBAL_NUM_THREADS), method=hclust_method)
    
    split_groups <- cutree(hc, k=num_groups)
    
    ref_groups <- list()
    
    grp_counter = 0
    for (cut_group in unique(split_groups)) {
        grp_counter = grp_counter + 1
        grp_name = sprintf("refgrp-%d", grp_counter)
        cell_names <- names(split_groups[split_groups == cut_group])
        ref_groups[[grp_name]] <- which(colnames(SPACEmapX_obj@expr.data) %in% cell_names)
    }
    
    SPACEmapX_obj@reference_grouped_cell_indices <- ref_groups
    
    return(SPACEmapX_obj)
}


#' @title remove_outliers_norm()
#'
#' @description Set outliers to some upper or lower bound.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param out_method method for computing the outlier bounds (default: "average_bound", involving
#'                   determining the range of values for each cell, and then taking the mean of those bounds.)
#'
#' @param lower_bound setting the lower bound for the data (default: NA, uses out_method above)
#'
#' @param upper_bound setting the upper bound for the data (default: NA, uses out_method above)
#'
#' @return SPACEmapX_obj with data bounds set accordingly.
#'
#' @keywords internal
#' @noRd
#'

remove_outliers_norm <- function(SPACEmapX_obj,
                                 out_method="average_bound",
                                 lower_bound=NA,
                                 upper_bound=NA) {
    
    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))
    
    
    SPACEmapX_obj@expr.data <- .remove_outliers_norm(data=SPACEmapX_obj@expr.data,
                                                    out_method=out_method,
                                                    lower_bound=lower_bound,
                                                    upper_bound=upper_bound)
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- remove_outliers_norm(SPACEmapX_obj@.hspike, out_method, lower_bound, upper_bound)
    }
    
    return(SPACEmapX_obj)
    
}

#' @keywords internal
#' @noRd
#'

.remove_outliers_norm <- function(data,
                                  out_method="average_bound",
                                  lower_bound=NA,
                                  upper_bound=NA) {
    
    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))
    
    
    
    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        flog.error("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }
    
    
    
    if ( (! is.na(lower_bound)) & (! is.na(upper_bound)) ) {
        
        ## using the values as specified.
        flog.info(paste("::remove_outlier_norm: using hard thresholds: ",
                        "lower_bound:" , lower_bound,
                        "upper_bound:", upper_bound) )
        
    } else if (! is.na(out_method)) {
        
                                        # using out_method instead of specified bounds.
        flog.info(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))
        
        if (out_method == "average_bound"){
            
            bounds = .get_average_bounds(data)
            lower_bound = bounds[1]
            upper_bound = bounds[2]
            
            flog.info(sprintf("outlier bounds defined between: %g - %g", lower_bound, upper_bound))
            
        } else {
            flog.error(paste("::remove_outlier_norm:Error, please",
                             "provide an approved method for outlier",
                             "removal for visualization."))
            stop(991)
        }
    } else {
        flog.error("::remove_outlier_norm:Error, must specify outmethod or define exact bounds")
        stop(992)
    }
    
                                        # apply bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound
    
    
    return(data)
}




#' @title center_cell_expr_across_chromosome()
#'
#' @description Centers expression data across all genes for each cell, using the cell mean or median expression
#' value as the center.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param method method to select the center of the cell expression value. (default: 'mean', options: 'mean,'median')
#'
#' @return SPACEmapX_object
#'
#' @keywords internal
#' @noRd
#'

center_cell_expr_across_chromosome <- function(SPACEmapX_obj, method="mean") { # or median
    
    flog.info(paste("::center_smooth across chromosomes per cell"))
    
    SPACEmapX_obj@expr.data <- .center_columns(SPACEmapX_obj@expr.data, method)
    
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- center_cell_expr_across_chromosome(SPACEmapX_obj@.hspike, method)
    }
    
    
    return(SPACEmapX_obj)
}

#' @keywords internal
#' @noRd
#'

.center_columns <- function(expr_data, method) {
    
                                        # Center within columns (cells)
    if (method == "median") {
        row_median <- apply(expr_data, 2, function(x) { median(x, na.rm=TRUE) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_median))
    }
    else {
                                        # by mean
        row_means <- apply(expr_data, 2, function(x) { mean(x, na.rm=TRUE) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_means))
    }
    return(expr_data)
}




#' @title require_above_min_mean_expr_cutoff ()
#'
#' @description Filters out genes that have fewer than the corresponding mean value across all cell values.
#'
#' @param SPACEmapX_obj  SPACEmapX_object
#'
#' @param min_mean_expr_cutoff  the minimum mean value allowed for a gene to be retained in the expression matrix.
#'
#' @return SPACEmapX_obj  the SPACEmapX_object with lowly or unexpressed genes removed.
#'
#' @keywords internal
#' @noRd
#'

require_above_min_mean_expr_cutoff <- function(SPACEmapX_obj, min_mean_expr_cutoff) {
    
    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))
    
    
    indices <-.below_min_mean_expr_cutoff(SPACEmapX_obj@expr.data, min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(sprintf("Removing %d genes from matrix as below mean expr threshold: %g",
                          length(indices), min_mean_expr_cutoff))
        
        SPACEmapX_obj <- remove_genes(SPACEmapX_obj, indices)
        
        expr_dim = dim(SPACEmapX_obj@expr.data)
        flog.info(sprintf("There are %d genes and %d cells remaining in the expr matrix.",
                          expr_dim[1], expr_dim[2]))
        
    }
    
    return(SPACEmapX_obj)
    
}

#' @keywords internal
#' @noRd
#'

.below_min_mean_expr_cutoff <- function(expr_data, min_mean_expr) {
    
    average_gene <- rowMeans(expr_data)
    
                                        # Find averages above a certain threshold
    indices <- which(average_gene < min_mean_expr)
    
    return(indices)
    
}




#' @title require_above_min_cells_ref()
#'
#' @description Filters out genes that have fewer than specified number of cells expressing them.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#'
#' @return SPACEmapX_obj SPACEmapX_object with corresponding genes removed.
#'
#' @keywords internal
#' @noRd
#'

require_above_min_cells_ref <- function(SPACEmapX_obj, min_cells_per_gene) {
    
    genes_passed = which(apply(SPACEmapX_obj@expr.data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))
    
    num_genes_total = dim(SPACEmapX_obj@expr.data)[1]
    num_removed = num_genes_total - length(genes_passed)
    if (num_removed > 0) {
        
        flog.info(sprintf("Removed %d genes having fewer than %d min cells per gene = %g %% genes removed here",
                          num_removed, min_cells_per_gene, num_removed / num_genes_total * 100))
        
        
        if (num_removed == num_genes_total) {
            
            flog.warn(paste("::All genes removed! Must revisit your data..., cannot continue here."))
            stop(998)
        }
        
        
        SPACEmapX_obj <- remove_genes(SPACEmapX_obj, -1 * genes_passed)
        
        
    }
    else {
        
        flog.info("no genes removed due to min cells/gene filter")
        
    }
    
    return(SPACEmapX_obj)
    
}


#' @title clear_noise()
#'
#' @description Remove values that are too close to the reference cell expr average and are considered noise.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param threshold values within reference mean +- threshold are set to zero.
#'
#' @param noise_logistic uses a logistic (sigmoidal) function to noise removal.
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

clear_noise <- function(SPACEmapX_obj, threshold, noise_logistic=FALSE) {
    
    flog.info(paste("********* ::clear_noise:Start. threshold: ", threshold, sep=""))
    
    if (threshold == 0) {
        return(SPACEmapX_obj); # nothing to do
    }
    
    if (has_reference_cells(SPACEmapX_obj)) {
        ref_idx = get_reference_grouped_cell_indices(SPACEmapX_obj)
        mean_ref_vals = mean(SPACEmapX_obj@expr.data[,ref_idx])
    } else {
        ## no reference
        ## use mean of all data
        mean_ref_vals = mean(SPACEmapX_obj@expr.data)
    }
    
    if (noise_logistic) {
        
        SPACEmapX_obj <- depress_log_signal_midpt_val(SPACEmapX_obj, mean_ref_vals, threshold)
        
    } else {
        
        SPACEmapX_obj@expr.data <- .clear_noise(SPACEmapX_obj@expr.data, threshold, center_pos=mean_ref_vals)
    }
    
    # if (! is.null(SPACEmapX_obj@.hspike)) {
    #     flog.info("-mirroring for hspike")
    #     SPACEmapX_obj@.hspike <- clear_noise(SPACEmapX_obj@.hspike, threshold, noise_logistic)
    # }
    
    return(SPACEmapX_obj)
}

#' @keywords internal
#' @noRd
#'

.clear_noise <- function(expr_data, threshold, center_pos=0) {
    
    upper_bound = center_pos + threshold
    lower_bound = center_pos - threshold
    
    expr_data[expr_data > lower_bound & expr_data < upper_bound] = center_pos
    
    return(expr_data)
}



#' @title clear_noise_via_ref_mean_sd()
#'
#' @description Define noise based on the standard deviation of the reference cell expression data.
#' The range to remove noise would be mean +- sdev * sd_amplifier
#' where sd_amplifier expands the range around the mean to be removed as noise.
#' Data points defined as noise are set to zero.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise
#'                     range (default: 1.5)
#'
#' @param noise_logistic uses a logistic (sigmoidal) function to noise removal.
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

clear_noise_via_ref_mean_sd <- function(SPACEmapX_obj, sd_amplifier=1.5, noise_logistic=FALSE) {
    
    if (has_reference_cells(SPACEmapX_obj)) {
        ref_idx = get_reference_grouped_cell_indices(SPACEmapX_obj)
        flog.info("denoising using mean(normal) +- sd_amplifier * sd(normal) per gene per cell across all data")
    }
    else {
        ref_idx = unlist(SPACEmapX_obj@observation_grouped_cell_indices)
        flog.info("-no reference cells specified... using mean and sd of all cells as proxy for denoising")
    }
    vals = SPACEmapX_obj@expr.data[,ref_idx]
    
    mean_ref_vals = mean(vals)
    
    mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=TRUE))) * sd_amplifier
    
    upper_bound = mean_ref_vals + mean_ref_sd
    lower_bound = mean_ref_vals - mean_ref_sd
    
    flog.info(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                    lower_bound, "-", upper_bound, sep=" "))
    
    
    
    if (noise_logistic) {
        
        threshold = mean_ref_sd
        SPACEmapX_obj <- depress_log_signal_midpt_val(SPACEmapX_obj, mean_ref_vals, threshold)
        
    } else {
        # smooth_matrix <- SPACEmapX_obj@expr.data
        
        # smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = mean_ref_vals
        SPACEmapX_obj@expr.data[SPACEmapX_obj@expr.data > lower_bound & SPACEmapX_obj@expr.data < upper_bound] = mean_ref_vals
        
        # SPACEmapX_obj@expr.data <- smooth_matrix
    }
    
    # if (! is.null(SPACEmapX_obj@.hspike)) {
    #     flog.info("-mirroring for hspike")
    #     SPACEmapX_obj@.hspike <- clear_noise_via_ref_mean_sd(SPACEmapX_obj@.hspike, sd_amplifier, noise_logistic)
    # }
    
    return(SPACEmapX_obj)
}




# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr Indices of the chr in which the tails are to be removed.
# tail_length Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.

#' @keywords internal
#' @noRd
#'


.remove_tails <- function(smooth_matrix, chr, tail_length){
    
                                        #flog.info(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
        tail_length <- floor(chr_length / 3)
    }
    remove_indices <- chr[seq_len(tail_length)]
    remove_indices <- c(remove_indices,
                        chr[ ( (chr_length + 1) - tail_length):
                             chr_length])
    
    return(remove_indices)
}


#' @title smooth_by_chromosome()
#'
#' @description Smooth expression values for each cell across each chromosome by using a
#' moving average with a window of specified length.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param window_length length of window (number of genes) for the moving average
#'
#' @param smooth_ends perform smoothing at the ends of the chromosomes (default: TRUE)
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

smooth_by_chromosome <- function(SPACEmapX_obj, window_length, smooth_ends=TRUE) {
    
    gene_chr_listing = SPACEmapX_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=SPACEmapX_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window(data=chr_data,
                                               window_length=window_length)
            
            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))
            
            SPACEmapX_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- smooth_by_chromosome(SPACEmapX_obj@.hspike, window_length, smooth_ends)
    }
    
    
    return(SPACEmapX_obj)
}

#' @keywords internal
#' @noRd
#'

.smooth_window <- function(data, window_length) {
    
    flog.debug(paste("::smooth_window:Start.", sep=""))
    
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }
    
    ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
    
    data_sm <- apply(data,
                     2,
                     .smooth_helper,
                     window_length=window_length)
    
    
    ## Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    
    
    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))
    
    
    return(data_sm)
}


# Helper function for smoothing with a moving average.
#
# Args:
# obs_data: Data to smooth
# window_length:  Length of the window to smooth.
#
# Returns:
# Data smoothed.
##.smooth_helper <- function(obs_data, tail_length) {

#' @keywords internal
#' @noRd
#'

.smooth_helper <- function(obs_data, window_length) {
                                        # strip NAs out and replace after smoothing
    orig_obs_data = obs_data
    
    nas = is.na(obs_data)
    
    obs_data = obs_data[!nas]
    
    obs_length <- length(obs_data)
    end_data <- obs_data
    
    tail_length = (window_length - 1)/2
    if (obs_length >= window_length) {
        end_data <- .smooth_center_helper(obs_data, window_length)
    }
    
                                        # end_data will have the end positions replaced with mean values, smoothing just at the ends.
    
    obs_count <- length(obs_data)
    
    numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))
    
                                        # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
    iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))
    
    for (tail_end in seq_len(iteration_range)) {
        end_tail = obs_count - tail_end + 1
        
        d_left = tail_end - 1
        d_right = obs_count - tail_end
        d_right = ifelse(d_right > tail_length, tail_length, d_right)
        
        r_left = tail_length - d_left
        r_right = tail_length - d_right
        
        denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)
        
        left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
        right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]
        
        numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
        
        end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
        end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
    }
    
    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data
    
    return(orig_obs_data)
}

smooth_by_chromosome_coordinates <- function(SPACEmapX_obj, window_length) {
    
    gene_chr_listing = SPACEmapX_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=SPACEmapX_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window_by_coordinates(data=chr_data,
                                                              window_length=window_length,
                                                              genes=SPACEmapX_obj@gene_order[which(SPACEmapX_obj@gene_order$chr == chr) ,])
            
            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))
            
            SPACEmapX_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }
    
    if (! is.null(SPACEmapX_obj@.hspike)) {  ## for now not supported.
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- smooth_by_chromosome_coordinates(SPACEmapX_obj@.hspike, window_length=51)
    }
    
    
    return(SPACEmapX_obj)
}

.smooth_window_by_coordinates <- function(data, window_length, genes) {
    
    flog.debug(paste("::smooth_window:Start.", sep=""))
    
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }
    
    ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
    
    data_sm <- apply(data,
                     2,
                     .smooth_helper_by_coordinates,
                     window_length=window_length,
                     genes=genes)
    
    
    ## Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    
    
    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))
    
    
    return(data_sm)
}

.smooth_helper_by_coordinates <- function(obs_data, window_length=10000000, genes) {

    ### No handling of NAs for this method of smoothing.

    end_data <- obs_data

    for (i in seq_along(obs_data)) {

        current_pos = (genes$start[i] + genes$stop[i]) / 2

        around_indices = which((genes$start > (current_pos - window_length)) & (genes$stop < (current_pos + window_length)))
        if (length(around_indices) == 0) {  ## in case the window is smaller than the current gene's size
            around_indices = i
        }
        weights = 1 - abs(((genes$stop[around_indices] + genes$start[around_indices])/2) - current_pos)/window_length
        if (length(around_indices < 10)) {
            to_add = floor(length(around_indices - 10)/2)

            new_low = max(1, (min(around_indices) - to_add))
            new_high = min(length(obs_data), (max(around_indices) + to_add))
            weights = c(rep(0.1, (min(around_indices) - new_low)), weights, rep(0.1, (new_high - max(around_indices))))
            around_indices = seq(new_low, new_high)
        }

        end_data[i] = sum(obs_data[around_indices] * weights) / sum(weights)
    }

    return(end_data)
}



# Smooth vector of values over the given window length.
#
# Args:
# obs_data Vector of data to smooth with a moving average.
# window_length Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.

#' @keywords internal
#' @noRd
#'

.smooth_center_helper <- function(obs_data, window_length){
    
    nas = is.na(obs_data)
    vals = obs_data[! nas]
    
    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))
    
    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)
    
                                        #    flog.info(paste("custom filter = ", custom_filter, "\n and window_length =", window_length, "\n and nrow(data) = ", nrow(obs_data), sep=""))
    
                                        #smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)
    smoothed = filter(vals, custom_filter, sides=2)
    
    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]
    
    obs_data[! nas] = vals
    
    return(obs_data)
}


#' @title smooth_by_chromosome_runmeans
#'
#' @description uses the simpler caTools:runmeans() to perform smoothing operations.
#'
#' @param SPACEmapX_obj SPACEmapX object
#'
#' @param window_length window length to use for smoothing.
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'


smooth_by_chromosome_runmeans <- function(SPACEmapX_obj, window_length) {
    
    gene_chr_listing = SPACEmapX_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=SPACEmapX_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            chr_data = apply(chr_data, 2, caTools::runmean, k=window_length)
            
            SPACEmapX_obj@expr.data[chr_genes_indices, ] <- chr_data
        }
    }
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- smooth_by_chromosome_runmeans(SPACEmapX_obj@.hspike, window_length)
    }
    
    
    return(SPACEmapX_obj)
}






#' @title get_average_bounds()
#'
#' @description Computes the mean of the upper and lower bound for the data across all cells.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return (lower_bound, upper_bound)
#'
#' @keywords internal
#' @noRd
#'

get_average_bounds <- function (SPACEmapX_obj) {
    
    return(.get_average_bounds(SPACEmapX_obj@expr.data))
    
}


#' @keywords internal
#' @noRd
#'

.get_average_bounds <- function(expr_matrix) {
    
    lower_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))
    
    return(c(lower_bound, upper_bound))
}

#' @title log2xplus1()
#'
#' @description Computes log(x+1), updates SPACEmapX_obj@expr.data
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

log2xplus1 <- function(SPACEmapX_obj) {
    
    flog.info("transforming log2xplus1()")
    
    SPACEmapX_obj@expr.data <- log2(SPACEmapX_obj@expr.data + 1)
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- log2xplus1(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
    
}



#' @title invert_log2xplus1()
#'
#' @description Computes 2^x - 1
#' Updates SPACEmapX_obj@expr.data
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

invert_log2xplus1 <- function(SPACEmapX_obj) {
    
    flog.info("inverting log2xplus1()")
    
    SPACEmapX_obj@expr.data <- 2^SPACEmapX_obj@expr.data - 1
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- invert_log2xplus1(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
}


#' @title invert_log2()
#'
#' @description Computes 2^x
#' Updates SPACEmapX_obj@expr.data
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

invert_log2 <- function(SPACEmapX_obj) {
    
    flog.info("invert_log2(), computing 2^x")
    
    SPACEmapX_obj@expr.data <- 2^SPACEmapX_obj@expr.data
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- invert_log2(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
}


#' @title make_zero_NA()
#'
#' @description Converts zero to NA
#' Updates SPACEmapX_obj@expr.data
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

make_zero_NA <- function(SPACEmapX_obj) {
    
    flog.info("make_zero_NA()")
    
    SPACEmapX_obj@expr.data <- SPACEmapX_obj@expr.data[SPACEmapX_obj@expr.data == 0] <- NA
    
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- make_zero_NA(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
    
}


#' @title transform_to_reference_based_Zscores()
#'
#' @description Computes mean and standard deviation for the reference cells, then uses these values
#' to compute gene Z-scores for both the reference and observation cells.
#'
#' Note, reference cell gene expression values will be centered at zero after this operation.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

transform_to_reference_based_Zscores <- function(SPACEmapX_obj) {

    ## center and convert to z-scores
    flog.info(paste("::center_and_Zscore_conversion", sep=""))
    
    ## remember, genes are rows, cells are cols
    
    ## centering and z-scores based on the reference (normal) cells:

    ## ref data represent the null distribution
    ref_idx = get_reference_grouped_cell_indices(SPACEmapX_obj)
    
    ref_data = SPACEmapX_obj@expr.data[,ref_idx]
    
    gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=TRUE)})
    gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=TRUE)})
    
    ## assume at least Poisson level variation
    gene_ref_sd = pmax(gene_ref_sd, gene_ref_mean)
    
    ## center all genes at the ref (normal) center:
    SPACEmapX_obj@expr.data = sweep(SPACEmapX_obj@expr.data, 1, gene_ref_mean, FUN="-")
    
    ## convert to z-scores based on the ref (normal) distribution
    SPACEmapX_obj@expr.data = sweep(SPACEmapX_obj@expr.data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- transform_to_reference_based_Zscores(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
    
}


#' @title mean_center_gene_expr()
#'
#' @description mean-center all gene expression values across all cells
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

mean_center_gene_expr <- function(SPACEmapX_obj) {
    
    flog.info(paste("::centering", sep=""))
    
    SPACEmapX_obj@expr.data <- sweep(SPACEmapX_obj@expr.data, 1, rowMeans(SPACEmapX_obj@expr.data, na.rm=TRUE), FUN="-")
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- mean_center_gene_expr(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
}


#' @title get_reference_grouped_cell_indices()
#'
#' @description Retrieves the matrix indices for the columns correspoinding to the reference cells.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return vector of column indices
#'
#' @keywords internal
#' @noRd
#'

get_reference_grouped_cell_indices <- function(SPACEmapX_obj) {
    
    return( unlist(SPACEmapX_obj@reference_grouped_cell_indices) )
    
}


#' @title apply_max_threshold_bounds()
#'
#' @description Assumes centered at zero and sets bounds to +- threshold value.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param threshold value to threshold the data
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

apply_max_threshold_bounds <- function(SPACEmapX_obj, threshold) {
    
    flog.info(paste("::process_data:setting max centered expr, threshold set to: +/-: ", threshold))
    
    SPACEmapX_obj@expr.data[SPACEmapX_obj@expr.data > threshold] <- threshold
    SPACEmapX_obj@expr.data[SPACEmapX_obj@expr.data < (-1 * threshold)] <- -1 * threshold
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- apply_max_threshold_bounds(SPACEmapX_obj@.hspike, threshold)
    }
    
    return(SPACEmapX_obj)
}


#' @title remove_genes_at_ends_of_chromosomes()
#'
#' @description Removes genes that are within window_length/2 of the ends of each chromosome.
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param window_length length of the window to use.
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

remove_genes_at_ends_of_chromosomes <- function(SPACEmapX_obj, window_length) {
    
    contig_tail = (window_length - 1) / 2
    
    remove_indices <- c()
    gene_chr_listing = SPACEmapX_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    for (chr in chrs){
                                        #flog.info(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- .remove_tails(SPACEmapX_obj@expr.data,
                                    which(gene_chr_listing == chr),
                                    contig_tail)
        
                                        #flog.debug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))
        
        remove_indices <- c(remove_indices, remove_chr)
        
    }
    if (length(remove_indices) > 0){
        
        SPACEmapX_obj = remove_genes(SPACEmapX_obj, remove_indices)
        
        flog.info(paste("::process_data:Remove genes at chr ends, ",
                        "new dimensions (r,c) = ",
                        paste(dim(SPACEmapX_obj@expr.data), collapse=","),
                        " Total=", sum(SPACEmapX_obj@expr.data, na.rm=TRUE),
                        " Min=", min(SPACEmapX_obj@expr.data, na.rm=TRUE),
                        " Max=", max(SPACEmapX_obj@expr.data, na.rm=TRUE),
                        ".", sep=""))
        
    }
    else {
        flog.error("No genes removed at chr ends.... something wrong here")
        stop(1234)
    }
    
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- remove_genes_at_ends_of_chromosomes(SPACEmapX_obj@.hspike, window_length)
    }
    
    return(SPACEmapX_obj)
    
}


#' @title normalize_counts_by_seq_depth()
#'
#' @description Normalizes count data by total sum scaling
#'
#' For single cell data, a typical normalization factor is 1e5, providing counts per 100k total counts.
#' If a normalization factor is not provided, the median lib size is used.:
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param normalize_factor  total counts to scale the normalization to (default: NA, computed as described above)
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

normalize_counts_by_seq_depth <- function(SPACEmapX_obj, normalize_factor=NA) {
    
    data <- SPACEmapX_obj@expr.data
    
    normalized_data <- .normalize_data_matrix_by_seq_depth(data, normalize_factor)
    
    
    SPACEmapX_obj@expr.data <- normalized_data
    
    return(SPACEmapX_obj)
    
}


#' @keywords internal
#' @noRd
#'

.normalize_data_matrix_by_seq_depth <- function(counts.matrix, normalize_factor=NA) {
    
    flog.info("normalizing counts matrix by depth")
    
    data <- counts.matrix
    
    cs = colSums(data)
    
    ## make fraction of total counts:
    data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
    
    if (is.na(normalize_factor)) {

        normalize_factor = median(cs)
        
        flog.info(sprintf("Computed total sum normalization factor as median libsize: %f", normalize_factor))
        
    } else {
        flog.info(sprintf("Using specified normalization factor: %f", normalize_factor))
    }
    
    if (is.na(normalize_factor)) {
        stop("Error, normalize factor not estimated")
    }
    
    data <- data * normalize_factor
    
    return(data)
    
}



#' @title anscombe_transform()
#'
#' @description Performs Anscombe's transformation:
#'    y = 2 * sqrt(x + 3/8)
#' as per
#' https://en.wikipedia.org/wiki/Anscombe_transform
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

anscombe_transform <- function(SPACEmapX_obj) {
    
    SPACEmapX_obj@expr.data <- 2 * sqrt(SPACEmapX_obj@expr.data + 3/8)
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- anscombe_transform(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
    
}

#' @keywords internal
#' @noRd
#'
add_pseudocount <- function(SPACEmapX_obj, pseudocount) {
    
    flog.info(sprintf("Adding pseudocount: %g", pseudocount))
    
    SPACEmapX_obj@expr.data = SPACEmapX_obj@expr.data + pseudocount
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- add_pseudocount(SPACEmapX_obj@.hspike, pseudocount)
    }
    
    return(SPACEmapX_obj)
}


#' @title scale_SPACEmapX_expr
#'
#' @description performs scaling to expression values for each cell,
#' assigning all values to a standard normal centered at zero.
#'
#' @param SPACEmapX_obj SPACEmapX object
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

scale_SPACEmapX_expr <- function(SPACEmapX_obj) {
    
    flog.info("-scaling expr data")
    SPACEmapX_obj@expr.data = t(scale(t(SPACEmapX_obj@expr.data)))
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- scale_SPACEmapX_expr(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
}



#' @keywords internal
#' @noRd
#'

cross_cell_normalize <- function(SPACEmapX_obj) {
    
    ## using upper quartile normalization
    
    flog.info("-cross cell normalization")
    
    upper_quart = apply(SPACEmapX_obj@expr.data, 2, quantile, probs=0.75)
    mean_upper_quart = mean(upper_quart)
    SPACEmapX_obj@expr.data = sweep(SPACEmapX_obj@expr.data, 2, mean_upper_quart/upper_quart, "*")
    
    
    if (! is.null(SPACEmapX_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        SPACEmapX_obj@.hspike <- cross_cell_normalize(SPACEmapX_obj@.hspike)
    }
    
    return(SPACEmapX_obj)
    
    
}



#' @title .redo_hierarchical_clustering()
#'
#' @description Recomputes hierarchical clustering for subclusters
#'
#'
#' @param SPACEmapX_obj SPACEmapX_object
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return SPACEmapX_obj
#'
#' @keywords internal
#' @noRd
#'

.redo_hierarchical_clustering <- function(SPACEmapX_obj, hclust_method) {

    subclusters = SPACEmapX_obj@tumor_subclusters$subclusters

    for (grp_name in names(subclusters)) {

        grp_cell_idx = unlist(subclusters[[grp_name]], recursive=TRUE)
        grp_cell_idx = grp_cell_idx[order(grp_cell_idx)] # retain relative ordering in the expr matrix
        
        grp_expr_data = SPACEmapX_obj@expr.data[, grp_cell_idx, drop=FALSE]
        
        hc <- hclust(parallelDist(t(grp_expr_data), threads=SPACEmapX.env$GLOBAL_NUM_THREADS), method=hclust_method)

        SPACEmapX_obj@tumor_subclusters$hc[[grp_name]] <- hc
    }

    return(SPACEmapX_obj)
    
}


# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
compareNA <- function(v1,v2) {
    if ((is.null(v1) != is.null(v2)) || (length(v1) != length(v2))) {
        return(FALSE)
    }
    if (all((v1 == v2) | (is.na(v1) & is.na(v2)))) {
        return(TRUE)
    }
    return(FALSE)
}

#' @return FALSE if relevant args are not identical, TRUE if they are
#'
#' @keywords internal
#' @noRd
#'

.compare_args <- function(current_args, relevant_args, loaded_options) {

    diff_a = any(relevant_args %in% names(current_args[!(names(current_args) %in% names(loaded_options))]))
    diff_b = any(relevant_args %in% names(loaded_options[!(names(loaded_options) %in% names(current_args))]))

    if (diff_a || diff_b) {
        return(FALSE)
    }

    return(all(sapply(relevant_args[relevant_args %in% names(current_args)], function(arg_n) {
        all(compareNA(current_args[[arg_n]], loaded_options[[arg_n]]))
    })))
}


#' @keywords internal
#' @noRd
#'

.get_relevant_args_list <- function(
# .get_relevant_args_list <- function(out_dir=NULL,
#                                     HMM=FALSE,
#                                     HMM_type='i6',
#                                     num_ref_groups=NULL,
#                                     analysis_mode='samples',
#                                     tumor_subcluster_partition_method,
#                                     smooth_method='pyramidinal',
#                                     max_centered_threshold=3,
#                                     remove_genes_at_chr_ends=FALSE,
#                                     prune_outliers=FALSE,
#                                     BayesMaxPNormal=0.5,
#                                     mask_nonDE_genes=FALSE,
#                                     denoise=FALSE,
#                                     noise_filter=NA,
#                                     sd_amplifier=1.5,
#                                     noise_logistic=FALSE,
                                    ...) {
    run_arguments = list(...)

    # creation args  ## add check for matrix size, col/row names, and ref/obs indices
    relevant_args = vector("list", 21) # 21 = max steps
    expected_file_names = vector("character", 21)
    step_i = 1

    # 1 _incoming_data.SPACEmapX_obj
    relevant_args[[step_i]] = c("chr_exclude", "max_cells_per_group", "min_max_counts_per_cell", "counts_md5")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_incoming_data.SPACEmapX_obj", step_i))
    step_i = step_i + 1

    # 2 _reduced_by_cutoff.SPACEmapX_obj
    relevant_args[[step_i]] = c("cutoff", "min_cells_per_gene")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_reduced_by_cutoff.SPACEmapX_obj", step_i))
    step_i = step_i + 1

    # 3 _normalized_by_depth%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type") 
        if (run_arguments$HMM_type == "i6") {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "sim_method", "hspike_aggregate_normals", "sim_foreground")
        }
    }
    resume_file_token = ifelse((run_arguments$HMM), paste0("HMM", run_arguments$HMM_type), "")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_normalized_by_depth%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 4 _logtransformed%s.SPACEmapX_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_logtransformed%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 5 _scaled%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("scale_data")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_scaled%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 6 _split_%sf_refs%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("num_ref_groups")
    if (!is.null(run_arguments$num_ref_groups)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_split_%sf_refs%s.SPACEmapX_obj", step_i, resume_file_token, run_arguments$num_ref_groups))
    }
    step_i = step_i + 1

    # 7 _tumor_subclusters%s.%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method == 'random_trees') {
        resume_file_token = paste0(resume_file_token, ".rand_trees")
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method", "tumor_subcluster_pval", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.%s.SPACEmapX_obj", step_i, resume_file_token, run_arguments$tumor_subcluster_partition_method))
    }
    step_i = step_i + 1

    # 8 _remove_ref_avg_from_obs_logFC%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_logFC%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 9 _apply_max_centered_expr_threshold%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("max_centered_threshold")
    if(!is.na(run_arguments$max_centered_threshold)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "max_centered_threshold")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_apply_max_centered_expr_threshold%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 10 _smoothed_by_chr%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("smooth_method", "window_length")
    if (run_arguments$smooth_method == 'pyramidinal') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "smooth_ends")
    }
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_smoothed_by_chr%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 11 _recentered_cells_by_chr%s.SPACEmapX_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_recentered_cells_by_chr%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 12 _remove_ref_avg_from_obs_adjust%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_adjust%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 13 _remove_gene_at_chr_ends%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("remove_genes_at_chr_ends")
    if (run_arguments$remove_genes_at_chr_ends == TRUE) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "window_length")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_gene_at_chr_ends%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 14 _invert_log_transform%s.SPACEmapX_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_invert_log_transform%s.SPACEmapX_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 15 _tumor_subclusters%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method", "z_score_filter")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method != 'random_trees') {
        resume_file_token = paste0(resume_file_token, '.', run_arguments$tumor_subcluster_partition_method)
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "k_nn", "leiden_resolution", "tumor_subcluster_pval", "leiden_method", "leiden_function", "leiden_method_per_chr", "leiden_function_per_chr", "leiden_resolution_per_chr", "per_chr_hmm_subclusters", "per_chr_hmm_subclusters_references", "hclust_method", "cluster_by_groups") 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    else if (run_arguments$analysis_mode != 'subclusters') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "tumor_subcluster_pval", "hclust_method", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_no_subclustering%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 16 _remove_outlier%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("prune_outliers")
    if (run_arguments$prune_outliers) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "outlier_method_bound", "outlier_lower_bound", "outlier_upper_bound")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_outlier%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 17 _HMM_pred%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type", "analysis_mode", "HMM_transition_prob", "HMM_report_by")
        hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", run_arguments$analysis_mode)
        # hmm.SPACEmapX_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.SPACEmapX_obj", step_i, hmm_resume_file_token))  ########
        if (run_arguments$analysis_mode == 'subclusters' && run_arguments$tumor_subcluster_partition_method == 'random_trees') {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        }

        if (run_arguments$HMM_type == "i6") {
            # relevant_args[[step_i]] = c(relevant_args[[step_i]], )

        }
        else {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_i3_pval", "HMM_i3_use_KS")
        }
        expected_file_names[[step_i]] = hmm.SPACEmapX_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.SPACEmapX_obj", step_i, hmm_resume_file_token))
    }
    step_i = step_i + 1

    # 18 _HMM_pred.Bayes_Net%s.mcmc_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.SPACEmapX_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.SPACEmapX_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
    }
    step_i = step_i + 1

    # 19 _HMM_pred.Bayes_Net%s.Pnorm_%g.SPACEmapX_obj
    relevant_args[[step_i]] = c("HMM", "BayesMaxPNormal")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        # relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.SPACEmapX_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.SPACEmapX_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.SPACEmapX_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 20 _HMM_pred.repr_intensities%s.Pnorm_%g.SPACEmapX_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "BayesMaxPNormal")
        # hmm.SPACEmapX_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.SPACEmapX_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.SPACEmapX_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 21 _mask_nonDE%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("mask_nonDE_genes")
    if (run_arguments$mask_nonDE_genes) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "require_DE_all_normals", "test.use", "mask_nonDE_pval")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_mask_nonDE%s.SPACEmapX_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 22 _denoise%s.NF_%s.SD_%g.NL_%s.SPACEmapX_obj
    relevant_args[[step_i]] = c("denoise")
    if (run_arguments$denoise) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "noise_filter", "sd_amplifier", "noise_logistic")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_denoise%s.NF_%s.SD_%g.NL_%s.SPACEmapX_obj", step_i, resume_file_token, run_arguments$noise_filter, run_arguments$sd_amplifier, run_arguments$noise_logistic))
    }
    step_i = step_i + 1

    return(list(relevant_args=relevant_args, expected_file_names=expected_file_names))
}
