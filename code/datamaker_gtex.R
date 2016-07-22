################################
## Filename: datamaker_gtex.R
## Created by: Lei Sun
## Created on: 07/12/2016
## Synopsis: 
################################


library(grex)
library(limma)
library(edgeR)


## sample Nsample samples from tissue 1, Nsample samples from tissue 2, choose fully expressed genes only, order the genes from low to high expression, to form a count matrix
datamaker_gtex <- function(tissue, Nsample) {
    rawdata1 <- read.csv(paste0("../data/", tissue[1], ".csv"), header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
    gene_names1 = grex::cleanid(rawdata1[, 1])
    if (length(tissue) > 1) {
        rawdata2 <- read.csv(paste0("../data/", tissue[2], ".csv"), header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
        gene_names2 = grex::cleanid(rawdata2[, 1])
        gene_names = intersect(gene_names1, gene_names2)
        rawdata1 = rawdata1[is.element(gene_names1, gene_names), ]
        rawdata1 = rawdata1[order(rawdata1[, 1]), ]
        rawdata1 = rawdata1[, -c(1, 2)]
        rawdata2 = rawdata2[is.element(gene_names2, gene_names), ]
        rawdata2 = rawdata2[order(rawdata2[, 1]), ]
        rawdata1 = rawdata2[, -c(1, 2)]
        if (is.null(Nsample)) {
            Nsample <- min(dim(rawdata1)[2], dim(rawdata2)[2])
        }
        if (dim(rawdata1)[2] < Nsample | dim(rawdata2)[2] < Nsample) {
            stop("Not enough samples in the raw dataset!")
        }

            ## All genes are alternatives
            subsample1 <- sample(colnames(rawdata1), Nsample)
            counts1 <- rawdata1[, subsample1]
            subsample2 = sample(setdiff(colnames(rawdata2), subsample1), Nsample)
            counts2 <- rawdata2[, subsample2]
            counts <- cbind(counts1, counts2)
            subsample <- list(subsample1, subsample2)
     } else {
     	rawdata1 = rawdata1[, -c(1, 2)]
        if (is.null(Nsample)) {
            Nsample <- floor(dim(rawdata1)[2] / 2)
        }
        if (dim(rawdata1)[2] < 2 * Nsample) {
            stop("Not enough samples in the raw dataset!")
        }


            subsample <- sample(1 : (dim(rawdata1)[2]), 2 * Nsample)
            counts <- rawdata1[, subsample]
        gene_names <- gene_names1
    }

    ## Remove unexpressed genes (min expression > 0)
    geneInd = apply(counts, 1, min) > 0
    gene_names <- gene_names[geneInd]
    counts <- counts[geneInd, ]
    rm(geneInd)
    
    ## Remove genes without an entrez ID
    geneInd = !is.na(grex::grex(gene_names)[, 2])
    gene_names <- gene_names[geneInd]
    counts <- counts[geneInd, ]
    rm(geneInd)
    
    ## Remove genes on XY chromosomes
	
# 	## Separate all expressed genes to a designated gene set and the rest, take genes in gene set as Group 0
# 	gene_set <- intersect(dfargs$gene_set, gene_names)
# ##	if (length(gene_set) < length(dfargs$gene_set)*0.95) {
# ##		stop("Gene set not expressed: more than 5% of genes in this gene set not expressed at all")
# ##	}
# 	counts_gene_set = counts[is.element(gene_names,gene_set),]
# 	subsample_gene_set = subsample[is.element(gene_names,gene_set),]
# 	gene_names_gene_set = gene_names[is.element(gene_names,gene_set)]
# 	Ngene_gene_set = length(gene_names_gene_set)
# 
# 	counts = counts[!is.element(gene_names,gene_set),]
# 	subsample = subsample[!is.element(gene_names,gene_set),]
# 	gene_names = gene_names[!is.element(gene_names,gene_set)]
# 
#     ## Take another top Ngene high-expressed genes from the rest as Group 0
#     if (!is.null(dfargs$Ngene)) {
#         dfargs$Ngene <- min(dfargs$Ngene, dim(counts)[1])
#         subsample <- subsample[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
#         gene_names <- gene_names[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)])]
#         counts <- counts[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
#     }
#     dfargs$Ngene <- dim(counts)[1]
#     
#     ## Combine Group 0 and Group 1 genes
#     counts <- rbind(counts, counts_gene_set)
#     subsample <- rbind(subsample, subsample_gene_set)
#     Ngene_total = dfargs$Ngene + Ngene_gene_set

    ## order genes from high to low expression
    gene_order = order(rowSums(counts))
    gene_names <- gene_names[gene_order]
    counts <- counts[gene_order, ]
    rm(gene_order)
    
    
    ## Model's design: Nsamp samples for group A and Nsamp samples for group B
    condition <- factor(rep(1:2, each = Nsample))



    


    # ## Ground truth of null hypotheses: beta_g=0
    # null <- rep(0, Ngene_total)
    # null[sample(Ngene_total, round(Ngene_total * dfargs$nullpi))] <- 1
    # ## default is that dfargs$nullpi is 1, so all genes are null genes -- dcg

    # ## Poisson thinning (optional)
    # pois_out <- pois_thinning(counts, dfargs, null) ## does nothing if args$poisthin == FALSE -- dcg
    # counts <- pois_out$counts
    # true_log2foldchange <- pois_out$log2foldchanges

    # ## Mix null and alternative genes from different samples (optional)
    # mix <- mix_sample(counts, dfargs, null, subsample) ## only matters if not all genes are null genes -- dcg
    # counts <- mix$counts
    # subsample <- mix$subsample



	# input <- list(counts = counts, condition = condition, Group_0_gene_names = gene_names, Group_1_gene_names = gene_names_gene_set, Group_0_size = dfargs$Ngene, Group_1_size = Ngene_gene_set, expressed_gene_proportion = length(gene_set)/max(1,length(dfargs$gene_set)))
    # meta <- list(null = null, true_log2foldchange = true_log2foldchange, dfargs = dfargs)

    return(list(counts = counts, condition = condition, gene_names = gene_names, subsample = subsample))
}



# # default_datamaker_args <- function(args) {
    # ## poisthin: flag of Poisson thinning
    # if (is.null(args$poisthin)) {
        # args$poisthin <- FALSE
    # }

    # ## number of top genes to skip
    # if (is.null(args$skip_genes)) {
        # args$skip_genes <- 0
    # }

    # ## log2foldmean, log2foldsd: Poisson thinning params
    # if (args$poisthin == TRUE) {
        # if (is.null(args$log2foldmean)) {
            # args$log2foldmean <- 0
        # }
        # if (is.null(args$log2foldsd)) {
            # args$log2foldsd <- 1
        # }
    # }

    # ## breaksample: flag of each gene randomly select samples
    # if (is.null(args$breaksample)) {
        # args$breaksample <- FALSE
    # }


    # ## nullpi: proportion of null genes
    # if (is.null(args$nullpi)) {
        # if (args$poisthin == TRUE) {
            # args$nullpi <- 0.9
        # } else if (length(args$tissue) == 1) {
            # args$nullpi <- 1
        # } else if (length(args$tissue) > 1) {
            # args$nullpi <- 0
        # } else if (is.null(args$tissue)) {
            # args$nullpi <- 1
        # }
    # }

    # ## pseudocounts: add pseudocounts to count matrix
    # if (is.null(args$pseudocounts)) {
        # args$pseudocounts <- 1
    # }

    # if(is.null(args$sig_diag)) {
        # args$sig_diag <- rep(1, args$Ngene)
    # }

    # if(is.null(args$sig_alpha)) {
        # args$sig_alpha <- 1
    # }

    # if(is.null(args$beta0)) {
         # args$beta0 <- 10
    # }

    # if(is.null(args$get_null)) {
        # args$get_null <- TRUE
    # }

    # if (is.null(args$alt_type)) {
        # args$alt_type <- "normal"
    # }

    # if (is.null(args$log2fold_inflate_beta)) {
        # args$log2fold_inflate_beta <- 1
    # }

    # return(args)
# }

# selectsample <- function(counts, Nsample) {
	# subsample <- sample(1:dim(counts)[2], Nsample)
    # counts <- counts[, subsample]
    # subsample <- t(matrix(rep(subsample, dim(counts)[1]), ncol = dim(counts)[1]))
    # return(list(counts = counts, subsample = subsample))
# }

# pois_thinning <- function(counts, args, null) {
    # if (args$poisthin == TRUE) {
        # if (args$alt_type == "normal") {
            # log2foldchanges <- rnorm(sum(!null), mean = args$log2foldmean, sd = args$log2foldsd)
        # } else if (args$alt_type == "mixnorm") {
            # if (is.null(args$pi_vals)) {
                # stop("args$alt_type = \"mixnorm\" but pi_vals not specified")
            # } else if (is.null(args$tau_seq)) {
                # stop("args$alt_type = \"mixnorm\" but tau_seq not specified")
            # } else if (is.null(args$mu_seq)) {
                # stop("args$alt_type = \"mixnorm\" but mu_seq not specified")
            # }
            # log2foldchanges <- rmixnorm(pi_vals = args$pi_vals,
                                        # mu_seq = args$mu_seq,
                                        # tau_seq = args$tau_seq,
                                        # p = sum(!null))
            # cat("yay!\n")
        # }
        # log2foldchanges <- log2foldchanges * args$log2fold_inflate_beta ## inflation defaults to 1.

        # foldchanges <- 2 ^ log2foldchanges

        # ## thin group A
        # counts[which(!null)[log2foldchanges > 0], 1:args$Nsamp] <-
            # matrix(rbinom(sum(log2foldchanges >
            # 0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges >
            # 0], 1:args$Nsamp])), prob = rep(1 / foldchanges[log2foldchanges > 0], args$Nsamp)),
            # ncol = args$Nsamp)
        # ## thin group B
        # counts[which(!null)[log2foldchanges < 0], (args$Nsamp + 1):(2 * args$Nsamp)] <-
            # matrix(rbinom(sum(log2foldchanges <
            # 0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges <
            # 0], (args$Nsamp + 1):(2 * args$Nsamp)])), prob = rep(foldchanges[log2foldchanges <
            # 0], args$Nsamp)), ncol = args$Nsamp)

    # } else {
        # log2foldchanges <- rep(0, length = dim(counts)[1])
    # }
    # return(list(counts = counts, log2foldchanges = log2foldchanges))
# }

# mix_sample <- function(counts, args, null, subsample) {
    # if (args$nullpi < 1 & args$nullpi > 0 & args$breaksample == TRUE) {
        # newcounts <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        # newcounts[as.logical(null), ] <- counts[as.logical(null), 1:(2 * args$Nsamp)]
        # newcounts[!null, ] <- counts[!null, c(1:args$Nsamp, (2 * args$Nsamp + 1):(3 *
            # args$Nsamp))]
        # counts <- newcounts
        # newsubsample <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        # newsubsample[as.logical(null), ] <- subsample[as.logical(null), 1:(2 * args$Nsamp)]
        # newsubsample[!null, ] <- subsample[!null, c(1:args$Nsamp, (2 * args$Nsamp +
            # 1):(3 * args$Nsamp))]
        # subsample <- newsubsample
        # rm(newcounts)
        # rm(newsubsample)
    # }
    # return(list(counts = counts, subsample = subsample))
# }

# fit_ash <- function(out_obj) {
    # if(is.null(out_obj$sebetahat)) { ## some methods do not return sebetahat.
        # return()
    # } else {
        # ash_out <- ashr::ash(betahat = out_obj$betahat, sebetahat = out_obj$sebetahat,
                             # df = out_obj$df)
        # return(ash_out)
    # }
# }

# #' Get wald p-values and adjust by benjamini and hochberg.
# fit_freq_methods <- function(out_obj) {
    # if (is.null(out_obj$pvalue)) { ## if given pvalues, use those
        # if (is.null(out_obj$df)) {
            # p_values <- 2 * (1 - pnorm(abs(out_obj$betahat / out_obj$sebetahat)))
        # } else {
            # p_values <- 2 * (1 - pt(abs(out_obj$betahat / out_obj$sebetahat), df = out_obj$df))
        # }
    # } else {
        # p_values <- out_obj$pvalue
    # }

    # q_bh <- stats::p.adjust(p_values, method = "BH")

    # if (all(is.na(p_values))) {
        # q_storey <- NULL
    # } else {
        # q_storey <- qvalue::qvalue(p_values)
    # }
    # return(list(p_values = p_values, q_bh = q_bh, q_storey = q_storey))
# }

# extract_ashpi0 <- function(ash_obj) {
    # return(ash_obj$fitted.g$pi[1])
# }

# extract_ashlfdr <- function(ash_obj) {
    # return(ash_obj$lfdr)
# }

# extract_qvaluepi0 <- function(freq_obj) {
    # return(freq_obj$q_storey$pi0)
# }

# extract_pvalues <- function(freq_obj) {
    # return(freq_obj$p_values)
# }


# ##' Simple ordinary least squares.
# ##'
# ##' @param log_counts A matrix of numerics. The responses.
# ##' @param condition A matrix of numerics. The predictors.
# get_ols <- function(log_counts, condition) {
    # limma_out <- limma::lmFit(log_counts, model.matrix(~condition))
    # betahat <- limma_out$coefficients[, 2]
    # sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    # df <- limma_out$df.residual
    # return(list(betahat = betahat, sebetahat = sebetahat, df = df))
# }

# #' Generate a random sample from a mixture of normals.
# #'
# #' @param pi_vals The mixing proportions.
# #' @param mu_seq The mixing means.
# #' @param tau_seq The mixing standard deviations.
# #' @param p The number of samples.
# rmixnorm <- function (pi_vals, mu_seq, tau_seq, p)
# {
    # M <- length(pi_vals)
    # beta <- rep(NA, length = p)
    # which.mix <- sample(1:M, size = p, replace = TRUE, prob = pi_vals)
    # for (index in 1:M) {
        # current.ind <- which.mix == index
        # n_m <- sum(current.ind)
        # if (n_m > 0) {
            # beta[current.ind] <- stats::rnorm(n = n_m, mean = mu_seq[index],
                # sd = tau_seq[index])
        # }
    # }
    # return(beta)
# }