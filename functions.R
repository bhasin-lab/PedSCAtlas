# getGenes(text)
# extract gene names from text entry in BulkExpression Gene Set section
#' @param text text entry from PedSCAtlas
#' @return list of genes 
getGenes = function(text) {
    # find locations of new lines in text entry, separate each line as genes
    raw_gs = text
    newlines = gregexpr('\n',raw_gs)[[1]]
    if (newlines == -1) {
        if (any(grep(" ", raw_gs))) {
            return(gsub(" ","", raw_gs))
        } else {
            return(raw_gs)
        }
    }
    gene_i = toupper(substr(raw_gs, 1, newlines[1]))
    genes = gsub("\n", "", gene_i)
    for (i in 1:length(newlines)) {
        if (i < length(newlines)) {
            gene_i = toupper(substr(raw_gs, newlines[i], newlines[i+1]))
            gene_i = gsub("\n", "", gene_i)
            genes = c(genes, gene_i)
        } else {
            gene_i = toupper(substr(raw_gs, newlines[i], nchar(raw_gs)))
            gene_i = gsub("\n", "", gene_i)
            genes = c(genes, gene_i)
        }
    }    
    
    # make sure there are no spaces in gene names
    if (any(grep(" ", genes))) {
        space_genes = grep(" ", genes)
        for (i in 1:length(space_genes)) {
            genes[i] = gsub(" ","", genes[i])
        }
    }

    return(genes)
}

# getComparisons(dataset, comp)
#' get comparison list for stat_compare_means
#' @param comp
#' @param groups
#' @return list of comparisons
getComparisons = function(comp, groups) {
    other = groups[groups != comp]
    out = vector(mode = "list", length = length(other))
    for (i in 1:length(other)) {
        out[[i]] = c(comp, other[i])
    }
    return(out)
}

# ribbonPaths()
#' generate ribbon plot for plotting enrichment of multiple pathways (adapted from script by Chenbin Huang)
#' @param en
#' @param coords
#' @param group_ribbon
#' @param g1
#' @param g2
#' @param paths
#' @return ggplot object for ribbon plot of multiple pathways
ribbonPaths = function(en, coords, group_ribbon, g1, g2, paths) { 
    # if leukemia dataset, change all blast cell types to same label
    if ("B MPAL blast" %in% coords$CellType) {
        coords$CellType[grep("blast",coords$CellType)] = "blast"
    }

    # reformat dataframe
    en.table = en
    en.table = en.table[colnames(en.table) %in% paths]
    en.table$Group = unlist(coords[group_ribbon])
    en.table$CellType = unlist(coords$CellType)
    en.table = en.table %>% filter(Group == g1 | Group == g2)
    en.group = en.table$Group

    # create dataframe to fill in values
    en.summary = en.table %>% select(!Group) %>% pivot_longer(!CellType, values_to = "Enrichment", names_to = "Pathway") %>% group_by(CellType,Pathway) %>% summarize(avgEnrichment = mean(Enrichment))
    en.summary["diff"] = NA

    # calculate significance between g1 and g2, for each cell type and format table for ribbon plot
    ct = unique(en.table$CellType)
    for (i in 1:length(ct)) {
        en.filt = en.table %>% filter(CellType == ct[i])
        group.filt = en.filt$Group
        if (length(unique(group.filt)) == 2) {
            # rows as pathways, columns as cells
            en.lm = en.filt[,1:length(paths)] 
            en.lm = t(en.lm) 

            # calculate significance
            mod = model.matrix(~ factor(group.filt))
            colnames(mod) = c(g1, paste0(g2,"vs",g1))
            fit = lmFit(en.lm, mod)
            fit = eBayes(fit)
            res = decideTests(fit, p.value=0.01)
            table = topTable(fit, coef=2, n=Inf)

            for (j in 1:length(paths)) {
                path_row = table[rownames(table) == paths[j],]
                if (path_row$logFC>0 & path_row$adj.P.Val < 0.01) {
                    en.summary$diff[en.summary$Pathway == paths[j] & en.summary$CellType == ct[i]] = paste0("High_",g2)
                } else if (path_row$logFC<0 & path_row$adj.P.Val < 0.01) {
                    en.summary$diff[en.summary$Pathway == paths[j] & en.summary$CellType == ct[i]] = paste0("High_",g1)
                } 
            }
        }
    }

    # generate ribbon plot
    plot = ggplot(en.summary, 
        aes(y=avgEnrichment, axis1=CellType, axis2=Pathway)) + 
        geom_alluvium(aes(fill = diff)) + geom_stratum() +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
        scale_x_discrete(limits = c("CellType", "Pathway"), expand = c(.05, .05))
    plot = plot + ggtitle(paste("Enrichment Differences between", g2, "and", g1))

    return(plot)
}

# userGSEA(genes, exp)
#' perform ssgsea on user input biomarker set to get enrichment in chosen dataset
#' @param genes input genes, from getGenes(text)
#' @param exp loaded expression matrix 
#' @return encrichmet values of input biomarker set for each cell in exp
userGSEA = function(genes, exp) {
    genes = list(genes)
    names(genes) = "userBiomarkerSet_Enrichment"
    exp = selectMethod("summary","dgCMatrix")(exp)
    exp = as.matrix(exp)
    en = GSVA::gsva(exp, genes, method="ssgsea",verbose=TRUE)
    return(en)
}

# Get Classifier Predictions From User-Input 
#' @param model supervised classifier
#' @param data user-input data frame or matrix with colnames as genes and rownames as samples
#' @param genes character vector of features used to train classifier
#' @return data.frame of predicted probabilities columns represent classes, rows represent samples

cancerPred = function(model, data, genes) {
    # correct NME1-NME2 dash problem with classifiers
    if ("NME1-NME2" %in% colnames(data) | "NME1.NME2" %in% colnames(data)) {
        colnames(data)[which(colnames(data) == "NME1-NME2" | colnames(data) == "NME1.NME2")] = "NME1_NME2"
    }
    genes[which(genes == "NME1-NME2")] = "NME1_NME2"

    # correct for missing features
    if (!(all(genes %in% colnames(data)))) {
        print("Using Median Sample Values for Missing Genes ...")
        missing.genes = genes[-which(genes %in% colnames(data))]
        missing.data = matrix(rep(0, length(missing.genes)*nrow(data)), nrow = nrow(data))
        rownames(missing.data) = rownames(data)
        colnames(missing.data) = missing.genes
        for (i in 1:length(missing.genes)) {
            print(paste("Estimating for", missing.genes[i]))
            for (j in 1:nrow(missing.data)) {
                missing.data[j, missing.genes[i]] = median(data[j,])
            }
        }
        data = cbind(data, missing.data)
    }

    # generate cancer type predictions
    data.preds = data.frame(predict(model, data))

    return(data.preds)
}

#' provide cutoff threshold and gene expression values for specific gene
#'
#' @param gene the gene we want to use to classify samples based on expression
#' @param mat expression values (FPKM) tibble with samples as columns and symbol column containing gene names
#' @param os tibble with samples, Overall.Survival, and Vital.Status values
#' @param cut method to use to separate into high and low groups, can be either median (defualt) or cutP
#' @param type if there is a column "class" in the os object, then the disease subtype (class) to anlayze
#' @return the cutoff threshold and expression values (log2(FPKM+1)) as a list
#' @export
#' @examples
#' data = abs(data.frame(data = matrix(rnorm(100),ncol=10)))
#' expr = dplyr::tibble(symbol = paste0(rep("gene",10),1:10), data)
#' colnames(expr) = c("symbol", letters[1:10])
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' mCut("gene1", mat = expr, os = os, cut = "cutP")
mCut = function(gene, mat, os, cut = "median", type = "All") {
   # geneExp = mat %>% dplyr::filter(symbol == gene) %>% dplyr::select(!"symbol") %>% unlist(., use.names = FALSE) %>% as.numeric() # extract all values for gene    
    # split mat and os based on type input
    if (type != "All") {
        os = os %>% filter(class == type)
    }   
    mat = mat[,colnames(mat) %in% os$sample]
    geneExp = mat[rownames(mat) == gene,]
    if (cut == "median") {
        cutoff = median(geneExp[2:length(geneExp)]) # find median for gene (removing NA at beginning due to symbol column)
    } else if (cut == "cutP") {
        os["Exp"] = geneExp
        cox.os = suppressWarnings(survival::coxph(survival::Surv(time=os$Overall.Survival, event=os$Vital.Status=="Dead")~Exp, data = os)) # get survival based on Expression
        c = survMisc::cutp(cox.os)$Exp # get cutp output
        data.table::setorder(c, "Exp") # create table with cut points
        percentile <- ecdf(c$Exp) # find cdf of Expression points
        cutoff <- as.numeric(c[order(c$Q), ][nrow(c), 1]) # find optimal expression CutOff
        print(paste0("Cutoff = ", signif(cutoff,3)))
    } else {
        print("invalid cut method, please insert either median or cutP ... returning os without separation")
        os$group = "Low"
        return(os)
    }
    out = list(signif(cutoff,3), geneExp)
    names(out) = c("cutoff","geneExp")
    return(out)
}

#' provide cutoff threshold and GSEA values for gene set of two or more genes
#'
#' @param gs the genes we want to use to classify samples based on expression
#' @param exp expression values (FPKM) tibble with samples as columns and ENSMBL ids as rows
#' @param os tibble with samples, overall survival (OS) values, and Vital.Status values
#' @param ref reference tibble with ensmble ids and gene symbols as columns, converts to ENS using this
#' @param cut method to use to separate into high and low groups, can be either median (defualt) or cutP
#' @return the cutoff threshold and GSEA scores as a list
#' @export
#' @examples
#' expr = abs(matrix(rnorm(100),ncol=10))
#' rownames(expr) = paste0(rep("ENSG",10),1:10)
#' colnames(expr) = letters[1:10]
#' ref = dplyr::tibble(ensg = paste0(rep("ENSG",10), 1:10), symbol = paste0(rep("gene",10), 10:1))
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' mSet(c("gene1","gene4","gene5"), ENS = expr, os = os, ref = ref, cut = "cutP")
mSet = function(gs, exp, os, cut = "median", type) {
    #gs = findENS(gs, ref) # convert to ensmbl ids
    # filter based on type
    if (type != "All") {
        os = os %>% filter(class == type)
    }
    exp.filt = exp[rownames(exp) %in% unlist(gs),] # extract only genes of interest
    exp.filt = exp.filt[,colnames(exp.filt) %in% os$sample]
    #exp.filt = log2(exp.filt + 1) # retrieve log2(FPKM+1)
    scores = callGSVA(exp.filt,gs)
    if (cut == "median") {
        cutoff = median(scores) # find median score for gene set
    } else if (cut == "cutP") {
        os["Scores"] = 0
        # add scores into os
        for (i in 1:length(scores)) {
            os[os$sample == rownames(scores)[i],"Scores"] = scores[i]
        }
        cox.os = suppressWarnings(survival::coxph(survival::Surv(time=os$Overall.Survival, event=os$Vital.Status=="Dead")~Scores, data = os)) # get survival based on GSEA scores
        c = survMisc::cutp(cox.os)$Scores # get cutp output
        data.table::setorder(c, Scores) # create table with cut points
        percentile <- ecdf(c$Scores) # find cdf of GSEA score points
        cutoff <- as.numeric(c[order(c$Q), ][nrow(c), 1]) # find optimal GSEA score CutOff
        print(paste0("Cutoff = ", signif(cutoff,3)))
    } else {
        print("invalid cut method, please insert either median or cutP ... returning os without separation")
        os$group = "Low"
        return(os)
    }
    out = list(signif(cutoff,3),scores)
    names(out) = c("cutoff","scores")
    return(out) # return cutoff point and Scores
}

#' convert gene symbols to ENSMBL ids
#'
#' @param genes gene symbols
#' @param ref reference tibble with ensmble ids and gene symbols as columns
#' @return character vector of ENSMBL ids for genes
#' @examples
#' genes = c("gene3","gene4","gene10")
#' ref = dplyr::tibble(ensg = paste0(rep("ENSG",10), 1:10), symbol = paste0(rep("gene",10), 10:1))
#' findENS(genes, ref)
findENS = function(genes, ref) {
    ens = "" # create empty vector to add onto
    for (i in 1:length(genes)) {
        ind = which(ref$symbol == genes[i]) # find index or indices of gene symbol
        if (length(ind) > 1) {
            ens = c(ens, ref$ensg[ind[1]]) # take first index
        } else if (length(ind) == 1) {
            ens = c(ens, ref$ensg[ind])
        } else {
            print(paste(genes[i], "not found in ENSEMBL reference."))
        }
    }
    ens = ens[2:length(ens)] # remove "" at beginning
    return(ens)
}

# function survStats
    # inputs:
    # os -- tibble with samples, OS values, and filled in column "group"
    # outputs:
    # stats -- tibble filled in with test statistics (Comparison, survdiffP, coxHR, coxP)
#' calculate survival statistics for split survival data
#' 
#' @param os tibble with sample, OS, Vital.Status, and group (High/Low) columns
#' @return stats tibble filled with test statistics for survival data: Comparison, survdiffP, coxHR, coxP
#' @export
survStats = function(os) {
    stats = data.frame("Comparison" = "", "survdiffP" = "", "coxHR" = "", "coxP" = "")
    surv = survival::Surv(time=os$OS, event=os$Vital.Status=="Dead")
    sdOut = survival::survdiff(surv~group, data = os) # save survdiff output
    cpOut = survival::coxph(surv~group, data = os) # save coxph output
    stats$survdiffP = 1 - pchisq(sdOut$chisq, length(sdOut$n) - 1) # fill in survdiff p-value
    stats$coxHR = summary(cpOut)$coefficients[2] # fill in coxph Hazard Ratio
    stats$coxP = summary(cpOut)$coefficients[5] # fill in coxph p-value
    return(stats)
}

#' split survival data into High and Low groups based on cutoff value and scores/expression values
#'
#' @param os tibble with columns: samples, OS values, and Vital.Status
#' @param out output from mCut or mSet. It contains the cutoff threshold and either GSEA scores ("scores") or Gene Expression ("geneExp") values for each sample
#' @return the os["group"] column with samples split into High and Low groups
#' @export
#' @examples
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' out = list("cutoff" = 0.44, "scores" = rnorm(10))
#' os["group"] = cut(os,out)
#' os
cutGroups = function(os, out) {
    if (names(out)[2] == "geneExp") {
        os["Exp"] = out[[2]]
        os[os$Exp > out[[1]],"group"] = "High"
        os[os$Exp <= out[[1]],"group"] = "Low"
    } else {
        os["Scores"] = out[[2]]
        os[os$Scores > out[[1]],"group"] = "High"
        os[os$Scores <= out[[1]],"group"] = "Low"
    }
    os$group = factor(os$group, levels = c("Low","High"))# change levels so that the reference group is Low expression
    return(os$group)
}

# FROM SURVIVALGENIE
#'@name callGSVA
#'@aliases callGSVA
#'@title GSVA enrichment analysis
#'@description Estimates GSVA enrichment zscores (from SurvivalGenie).
#'@usage callGSVA(x,y)
#'@param x A data frame or matrix of gene or probe expression values where rows corrospond to genes and columns corrospond to samples
#'@param y A list of genes as data frame or vector
#'@details This function uses "zscore" gene-set enrichment method in the estimation of gene-set enrichment scores per sample.
#'@return A gene-set by sample matrix of GSVA enrichment zscores.
#'@import GSVA
#'@examples 
#'g <- 10 ## number of genes
#'s <- 30 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'rnames <- paste("g", 1:g, sep="")
#'cnames <- paste("s", 1:s, sep="")
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(rnames, cnames))
#'## genes of interest
#'genes <- paste("g", 1:g, sep="")
#'## Estimates GSVA enrichment zscores.
#'callGSVA(expr,genes)
#'@seealso GSVA
callGSVA = function(x,y) {
    if(missing(x)){
    stop("input expression data missing!")
    }
    if(missing(y)){
    stop("input gene set missing!")
    }
    #genes <- list(set1=y)
    genes = list(y)
    gsva.results <- GSVA::gsva(x, genes, method="zscore",verbose=FALSE, parallel.sz=2)
    tr_gsva.results <- t(gsva.results)
    colnames(tr_gsva.results) <- c("GSVA score")
    return (tr_gsva.results)
}

#' plot the survival Kaplan-Meier curve
#' 
#' @param os overall survival tibble with High/Low groups formed
#' @param name name of the gene/geneSet for the title, default is "gene"
#' @return plotting object
#' @export
#' @examples 
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' out = list("cutoff" = 0.44, "scores" = rnorm(10))
#' os["group"] = cut(os,out)
#' plotSurv(os, "gene1")
plotSurv = function(os, name = "gene", type = "All") {
    os$group = factor(os$group, levels = c("Low","High"))
    # change survival time to months from days
    os$Overall.Survival = os$Overall.Survival * 0.0328767

    surv = survival::survfit(survival::Surv(time=Overall.Survival, event=Vital.Status=="Dead")~group, data = os, conf.type = "log-log") # calculate survival
    cpOut = survival::coxph(survival::Surv(time=Overall.Survival,event=Vital.Status=="Dead")~group, data = os) # save coxph output
    p = survminer::ggsurvplot(surv, data = os,
        palette = c("#4f5379","#d90429"), 
        pval = TRUE,
        xlab = "Time (Months)"
    )$plot + ggplot2::annotate("text", x = Inf, y = Inf, vjust = 2, hjust = 1, label = paste0("Cox HR (p-val): ", 
    round(summary(cpOut)$coefficients[2],3), " (", round(summary(cpOut)$coefficients[5],3),")"), size = 5) 

    if (type != "All") {
       p + ggplot2::ggtitle(paste("Survival in", type, "based on", name, "Expression"))
    } else {
        p + ggplot2::ggtitle(paste("Survival based on", name, "Expression"))
    }
}