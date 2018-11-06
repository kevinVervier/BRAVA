# install packages if needed
if(!require(biomaRt)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}

if(!require(GenomicRanges)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GenomicRanges", version = "3.8")
}


if(!require(RColorBrewer)){
  install.packages('RColorBrewer')
  library(RColorBrewer)
}

if(!require(gplots)){
  install.packages('gplots')
  library(gplots)
}
######################################################################
#paste function
"%&%" = function(a,b) paste(a,b,sep="")

# load random forest model weights for each feature
load('../data/randomForest_model_relativeimportance.Rdata')

##############################################################################################
# faster version of getContext
#chr = c(4,1,3,20,9,7,'X','X',6,6)
#pos = c(55643311,4429872,132424172,51948402,25552267,114055091,66765627,66748355,33057511,33057835)

getContextMat_opti <- function(chr=NULL,pos=NULL,snpids=NULL,weights = cur_weights){
  require(biomaRt)
  require(GenomicRanges)
  
    load('../data/BRAINSPAN/Brain_span.rdata') 
    bgd = apply(X_all[1:dim(X_all)[1],,],2:3,function(x)sum(x,na.rm=TRUE))/(dim(X_all)[1] - length(which(is.na(X_all[,1,1]) | is.nan(X_all[,1,1]))))
    row.names(bgd) = row.names(X_all[1,,])
    # hic_average
    load('../data/HiC/GSE77565_FBD_IC-heatmap-res-100k_average_connected_bins_Mat_corrected.Rdata') 
    # hic_average_weighted
    load('../data/HiC/GSE77565_FBD_IC-heatmap-res-100k_weighted_average_connected_bins_Mat_corrected.Rdata')


 
    # if only position are provided, look for orresponding SNP IDs on Ensembl database
  if(!is.null(snpids)){
    # map ID on locations
    snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="grch37.ensembl.org") 
    snp_attributes = c("refsnp_id", "chr_name", "chrom_start")
    snp_loc = getBM(attributes=snp_attributes, filters="snp_filter", 
                    values=snpids, mart=snp_mart)
    # re-order snps
    snp_loc = snp_loc[match(snpids,snp_loc$refsnp_id),]
  }else{
    snp_loc = cbind(chr,pos)
    colnames(snp_loc) = c('chr_name','chrom_start')
    snp_loc = as.data.frame(snp_loc,stringsAsFactors = FALSE)
    snp_loc$chrom_start = as.numeric(snp_loc$chrom_start)
  }
  
  # remove unmapped SNPs
  snp_loc <- snp_loc[which(!is.na(snp_loc[,1])),]
  #mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")
  mart.hs <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  # get loci
  if(!is.null(snpids)){
    target = snp_loc[,1]#snpids
  }else{target = paste(chr,pos,sep='-')}
  
  loci = GRanges(
    seqnames = Rle(paste('chr',snp_loc$chr_name,sep='')), 
    ranges = IRanges(start = snp_loc$chrom_start,end = snp_loc$chrom_start)
  )
  
  # init output matrix
  mat = c(rep(list(matrix(0,nrow=dim(X_all[1,,])[1],ncol=dim(X_all[1,,])[2])),length(loci)))
  # find closest gene(s)
  g = getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position"), filters = "chromosomal_region", values = paste(snp_loc$chr_name,snp_loc$chrom_start,snp_loc$chrom_start,sep=':'), mart = mart.hs)
  # if any gene is found
  if (nrow(g)  > 1){
    g.gr = GRanges(
      seqnames = Rle(paste('chr',g$chromosome_name,sep='')), 
      ranges = IRanges(start = g$start_position,end = g$end_position),
      ens = g$ensembl_gene_id
    )
    
    # locate which loci are in genic regions
    map = findOverlaps(loci, g.gr) 
    idx.genic = unique(queryHits(map))#@queryHits)
    #if(length(idx.genic) > 0){
    idx.intergenic = c(1:length(loci))[-idx.genic]
    
    ###################
    # intragenic case #
    # for loci found in gene body we estimate its context only as expression level in brain.
    
    map@metadata$ens = g.gr@elementMetadata$ens[subjectHits(map)]#@subjectHits]
    mat[idx.genic] = lapply(idx.genic,function(i){
      idx = match(map@metadata$ens[queryHits(map)==i],Meta_data_brainSpan$ensembl_gene_id)
      if(length(which(is.na(idx))) > 0) idx = idx[-which(is.na(idx))]
      # if at least one gene is found as a match between ENSEMBL and BRAINSPAN
      if(length(idx) == 1) return(X_all[idx,,] - bgd)
      # if more than one gene is found at a given location, we sum their expression
      if(length(idx) > 1){
        return(colSums(X_all[idx,,])/length(idx) - bgd)
      }
      # no gene found --> intergenic
      if(length(idx) == 0){
        return(matrix(-1000,nrow=dim(X_all[1,,])[1],ncol=dim(X_all[1,,])[2])) # distinct value
      }
    })
    
    idx.intergenic = c(idx.intergenic, idx.genic[which(sapply(mat[idx.genic],function(x)x[1,1]) == -1000)])
    # reboot 
    mat[idx.genic[which(sapply(mat[idx.genic],function(x)x[1,1]) == -1000)]] = matrix(0,nrow=dim(X_all[1,,])[1],ncol=dim(X_all[1,,])[2])
    
  }else{
    idx.intergenic = c(1:length(loci))
  }
  
  ##################
  #intergenic case #
  # for intergenic loci, we combine different features ti estimate its context
  if(length(idx.intergenic)>0){
    ###################  
    # get closest eQTL --> will not work on chrX
    close.eqtl = fetch_eQTL(paste('chr',snp_loc$chr_name[idx.intergenic],sep=''),snp_loc$chrom_start[idx.intergenic],closest=TRUE,brain.only=TRUE,db.brain = db.brain)$prox.brain # just need to do it once, given we consider all brain eQTLs as important, not just region specific
    seqlevelsStyle(close.eqtl) <- "UCSC"
    # check for not found eQTLs (especially for chr X and Y)
    idx = idx.intergenic[nearest(close.eqtl,loci[idx.intergenic])]
    if(length(idx) < length(idx.intergenic)) mis.idx = idx.intergenic[which(!idx.intergenic %in% idx)]
    # remove the eQTLs further than 100,000 bp
    mis.idx = idx.intergenic[which(close.eqtl@elementMetadata$dist.brain > 100000)]
    idx = idx.intergenic[which(!idx.intergenic %in% mis.idx)]
    
    # get expression
    if(length(idx) > 0){
      closest = nearest(loci[idx],close.eqtl)
      if(length(test <- which(is.na(closest))) > 0) idx = idx[-test]
      tmp = get_exp_from_ens_matrix(close.eqtl[nearest(loci[idx],close.eqtl)],rescaled = rescaled) 
      if(length(dim(tmp)) == 3){ # more than one SNP
        mat[idx] = lapply(1:length(idx),function(i){
          if(!is.na(tmp[i,1,1])) mat[[idx[i]]] =  weights[2,]*(tmp[i,,] - bgd)
        })
      }else{
        mat[[idx]] =   weights[2,]*(tmp - bgd)
      }
    }
    
    ###################
    # get closest gene
    close.gene = closest_gene(paste('chr',snp_loc$chr_name[idx.intergenic],sep=''),snp_loc$chrom_start[idx.intergenic])
    mis.idx = idx.intergenic[which(close.gene@elementMetadata$distance > 100000)]
    idx = idx.intergenic[which(!idx.intergenic %in% mis.idx)]
    
    # get expression
    if(length(idx) > 0){
      tmp = get_exp_from_ens_matrix(close.gene[nearest(loci[idx],close.gene)],rescaled = rescaled) 
      if(length(dim(tmp)) == 3){ # more than one SNP
        mat[idx] = lapply(1:length(idx),function(i){
          if(!is.na(tmp[i,1,1])) mat[[idx[i]]] = mat[[idx[i]]] + weights[3,]*(tmp[i,,] - bgd)
        })
      }else{
        mat[[idx]] =  mat[[idx]] + weights[3,]*(tmp - bgd)
      }
    }
    
    ########################
    # get average TADD_ll_cp
    no_tad_idx =  which(close.gene@elementMetadata$distance > 100000)
    # which gene matrices we need
    if(length(no_tad_idx)>0){
      tmp.idx = idx.intergenic[-no_tad_idx]
      
    }else{
      tmp.idx = idx.intergenic
    }
    if(length(tmp.idx) > 0){
      map <- in_tad(paste('chr',snp_loc$chr_name[tmp.idx],sep=''),snp_loc$chrom_start[tmp.idx],option='cp')
      
      gene_set = ll_cp[subjectHits(map)]#@subjectHits]
      
      mat[idx.intergenic[-no_tad_idx]] = lapply(1:length(gene_set),function(i){
        idx = which(Meta_data_brainSpan$entrez_id %in% gene_set[[i]])
        if(length(idx)>1){
          mat[[tmp.idx[i]]] = mat[[tmp.idx[i]]]  + weights[1,]*(colSums(X_all[idx,,])/length(idx) - bgd)
        }else{
          if(length(idx == 1)){
            mat[[tmp.idx[i]]]  = mat[[tmp.idx[i]]] + weights[1,]*(X_all[idx,,] - bgd)
          }else{
            mat[[tmp.idx[i]]] = mat[[tmp.idx[i]]] 
          }
        }
      })
    }
    ################################
    # get HiC_unweighted_average
    # have Granges from bin location
    tmp.str = strsplit(names(hic_average_bin),split = '-')
    tmp.str = do.call('rbind',tmp.str)
    
    bins = GRanges(
      seqnames = Rle(paste('chr',tmp.str[,1],sep='')),
      ranges = IRanges(start = as.numeric(tmp.str[,2])+1,end = as.numeric(tmp.str[,2])+100000)
    )
    
    loci = GRanges(
      seqnames = Rle(paste('chr',snp_loc$chr_name[idx.intergenic],sep='')), 
      ranges = IRanges(start = snp_loc$chrom_start[idx.intergenic],end = snp_loc$chrom_start[idx.intergenic])
    )
    
    map = findOverlaps(loci, bins) 
    
    if(length(map) == length(idx.intergenic)){
      mat[idx.intergenic] = lapply(1:length(idx.intergenic),function(i) mat[[idx.intergenic[i]]] + weights[4,]*(hic_average_bin[subjectHits(map)][[i]]))#@subjectHits][[i]]))# not needed anymore, given we used the corrected version for HiC (removed after WCPG poster): - bgd)
    }else{
      # missing HiC ???
    }
    
    
    
    ###############################
    #get HiC_weighted_average
    #have Granges from bin location
    tmp.str = strsplit(names(hic_weighted_average_bin),split = '-')
    tmp.str = do.call('rbind',tmp.str)
    
    bins = GRanges(
      seqnames = Rle(paste('chr',tmp.str[,1],sep='')),
      ranges = IRanges(start = as.numeric(tmp.str[,2])+1,end = as.numeric(tmp.str[,2])+100000)
    )
    
    map = findOverlaps(loci, bins) 
    
    if(length(map) == length(idx.intergenic)){
      mat[idx.intergenic] = lapply(1:length(idx.intergenic),function(i) mat[[idx.intergenic[i]]] + weights[5,]*(hic_weighted_average_bin[subjectHits(map)][[i]]))#@subjectHits][[i]]))# not needed anymore, given we used the corrected version for HiC (removed after WCPG poster): - bgd)
    }else{
      # missing HiC ???
    }
    
    
    ##########################################
    # get weighted_genes_intraTAD_ll_cp
    # find relative distances 
    # get corresponding gene locations
    
    if(length(no_tad_idx)>0){
      idx.in.tad = idx.intergenic[-no_tad_idx]
    }else{
      idx.in.tad = idx.intergenic
    }
    
    #idx.in.tad = idx.intergenic[-no_tad_idx]
    if(length(idx.in.tad) > 0){
      mat[idx.in.tad] = lapply(1:length(idx.in.tad),function(i){
        idx = which(Meta_data_brainSpan$entrez_id %in% gene_set[[i]])
        idx2 = which(genes_loc$entrezgene %in% Meta_data_brainSpan$entrez_id[idx])
        if(length(idx2)>1){
          ref = snp_loc$chrom_start[idx.in.tad[i]]
          tmp = sapply(idx2,function(j) min(abs(ref-genes_loc$start_position[j]), abs(ref-genes_loc$end_position[j]))) # no need for an intragenic case here
          tmp = (1/tmp)/sum(1/tmp)
          return(mat[[idx.in.tad[i]]] + weights[6,]*(colSums(tmp*X_all[idx2,,],na.rm = TRUE)- bgd))
        }else{
          if(length(idx2) == 1){ # /!\ corrected typo
            # no need to weight genes (only one)
            return(mat[[idx.in.tad[i]]] + cur_weights[6,]*(X_all[idx2,,]- bgd))
          }else{
            mat[[idx.in.tad[i]]] = mat[[idx.in.tad[i]]]
          }
        }
      })
    }
  }
  
  
  names(mat) = target
  
  return(mat)
}


# Function for closest GTEx eQTL (brain or not)
# find the eQTLs in +/-offset regions
# if closest = TRUE: return only the closest eQTL for each input position (disregard of the offset value)
fetch_eQTL <- function(chr,position,offset=50000,closest=FALSE){
  # load brainspan for gen names matching (TODO: pre-filter GTEX data)
  load('../data/BRAINSPAN/Brain_span.rdata')
  #GTEx eQTLS files for brain tissues
  en.dir = '../data/GTEX/'
  db.brain = read.table(file=en.dir %&% 'Brain_snps.db',header=T,as.is = TRUE)
  # create Granges for loci
  loci = GRanges(
    seqnames = Rle(gsub(pattern = 'chr',replacement = '',chr)), 
    ranges = IRanges(start = position,end = position)
  )
  # filter db.brain based on genes available in BRAINSPAN
  idx = which(db.brain$ens.id %in% Meta_data_brainSpan$ensembl_gene_id) # 1,099,899 
  db.brain = db.brain[idx,]
  
  # check if the numbers of chr and pos match
  if(length(chr) == length(position)){
    if(closest == FALSE){
      idx.brain = sapply(1:length(chr), function(i) which(db.brain$snp_chrom == chr[i] & db.brain$snp_pos >= position[i]-offset & db.brain$snp_pos <= position[i]+offset)) 
      return(list('brain'=lapply(1:length(position),function(i) cbind(db.brain[idx.brain[[i]],],abs(db.brain$snp_pos[idx.brain[[i]]]-position[i])))))

    }else{
      brain.gr <- GRanges(
        seqnames = Rle(db.brain$snp_chrom), 
        ranges = IRanges(start = db.brain$snp_pos,end = db.brain$snp_pos),
        ens = db.brain$ens.id
      )
      idx.map.brain = nearest(loci, brain.gr)
      dist.brain = distanceToNearest(loci, brain.gr)
      tmp = brain.gr[subjectHits(dist.brain)]
      tmp@elementMetadata$dist.brain = dist.brain@elementMetadata$distance
      dist.brain = tmp
      
      return(list('prox.brain'=dist.brain))
    }
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}


# Function extracting spatiotemp. gene expression matrix based on Granges with 'ens' field
get_exp_from_ens <- function(gr,rescaled=FALSE){
  load('../data/BRAINSPAN/Brain_span.rdata') 
  idx = match(gr@elementMetadata$ens,Meta_data_brainSpan$ensembl_gene_id)
  # only return the cells used in the analysis (oldest, 5 brain regions)
  return(X_all[idx,,dim(X_all)[3]])
}

# Function extracting spatiotemp. gene expression matrix based on Granges with 'ens' field
get_exp_from_ens_matrix <- function(gr,rescaled=FALSE){

  load('../data/BRAINSPAN/Brain_span.rdata') 

  idx = match(gr@elementMetadata$ens,Meta_data_brainSpan$ensembl_gene_id)
  tmp = X_all[idx,,]
  nan.detect = apply(tmp,1, function(x) is.nan(x[1]))
  if(any(nan.detect)) tmp[which(nan.detect),,] = 0
  return(tmp)
}

#######################################################################
# closest gene function
closest_gene <- function(chr,position){
  library(GenomicRanges)
  load('../data/BRAINSPAN/Brain_span.rdata') 
  gene = read.table('../data/genes/gencode',header=FALSE,sep='\t')
  # remove the genes not found in BRAINSPAN, to make sure that the closest gene we selected has a matrix
  gene$V5 = gsub(pattern = '\\..*',replacement = '',gene$V5)
  idx = which(gene$V5 %in% Meta_data_brainSpan$ensembl_gene_id) # 18,741 / 22,327
  gene = gene[idx,]
  gene = GRanges(
    seqnames = Rle(gene$V1), 
    ranges = IRanges(start = gene$V3,end = gene$V4),
    ens = gene$V5
  )
  
  loci = GRanges(
    seqnames = Rle(chr), 
    ranges = IRanges(start = position,end = position)
  )
  tmp = distanceToNearest(loci, gene, ignore.strand=T)
  elementMetadata(loci) = c(elementMetadata(loci),mcols(gene[subjectHits(tmp)]))
  elementMetadata(loci) = c(elementMetadata(loci),mcols(tmp))
  return(loci)
}

##############################################
# in TAD function
in_tad <- function(chr,position,option='cp'){
  library(GenomicRanges)
  load('/sdata/Hi-C/geschwind/jake/genes_to_TADs.Rdata')
  load('/sdata/Hi-C/geschwind/jake/developing_brain_TADs.Rdata')
  
  loci = GRanges(
    seqnames = Rle(chr), 
    ranges = IRanges(start = position,end = position)
  )
  if(option == 'cp'){
    cp = as(cp,"GRanges")
  }else{
    cp = as(gz,"GRanges")
  }
  tmp = distanceToNearest(loci, cp, ignore.strand=T)
  return(tmp)
}


#####################################
# plot function
# input: context matrix
# region.order: brain regions order (vector of length 16), default order is alphabetical order
# output: figure

plotContext <- function(mat=NULL,outputFile=NULL,fix.quantile=FALSE,region.order=NULL){
  # get dimension names
  load('../data/BRAINSPAN/BrainSpan2cerebroViz.Rda')
  if(!is.null(outputFile)) pdf(outputFile,width = 100,height = 75)
  
  # column labels
  age_tmp = round(exp(seq(2,log(2118),length.out=50))) 
  f = function(x) {
    if(x <= 38) out=paste(x,'pcw',sep='')
    if(x > 38 & x < (38+52)) out =paste((x-38),"wks",sep='') 
    if(x > (38+52)) out =paste(round(x/52),"yrs",sep='') 
    return(out)
  }
  age_wk = sapply(age_tmp,f)
  # remove duplicates
  for(i in length(age_wk):2){
    if(age_wk[i-1] == age_wk[i]) age_wk[i] = ''
  }
  # potentially reorder the region names
  lab.region = rownames(expr[[1]])
  if(!is.null(region.order)) lab.region = lab.region[region.order]
  if(!is.list(mat)) mat = list(mat)
  for(i in 1:length(mat)){
    if(fix.quantile){
      quantile.range = quantile(as.vector(unlist(mat)), probs = seq(0, 1, 0.05))
    }else{
      quantile.range <- quantile(as.vector(mat[[i]]), probs = seq(0, 1, 0.05))
    }
    
    if(length(unique(quantile.range)) > 1){
      # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
      library(RColorBrewer)
      library(gplots)
      color.palette  <- colorRampPalette(c("#67a9cf", "#f7f7f7","#ef8a62"))(length(unique(quantile.range))-1)
      

        heatmap.2(mat[[i]],col=color.palette,Rowv=NA,Colv=NA,trace="none",density.info ='none',main = names(mat)[i], # re-ordering brain regions !!
                  labCol= age_wk, 
                  labRow= lab.region,
                  margins =c(30,36),
                  lmat=rbind(c(0,3),c(2,1),c(0,4)), 
                  lhei=c(0.1,9,0.1), 
                  lwid=c(0.5,9),
                  cexRow = 15,
                  cexCol = 15,
                  cex.main = 25,
                  breaks = unique(quantile.range),key = FALSE,dendrogram='none',scale='none')
      }
  }
  if(!is.null(outputFile)){
    dev.off()
  }else{
    return()
  }
}

# example
#plotContext(getContextMat(snpids = c('rs7266390','rs17177122'))[[2]])
