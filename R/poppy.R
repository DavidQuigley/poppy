#library(RColorBrewer)
#library(VariantAnnotation)
#library(CLONETv2)
#library(TPES)
#library(deconstructSigs)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(plyr)

#' Find copy number of two sets of CNA segments, e.g. pre vs post
#'
#' @description
#' This function would be used to generate a plot showing copy number across the
#' genome at two samples from the same individual obtained at two different timepoints.
#' It uses GenomicRanges to find the segments that are in both pre and post
#' and writes the value of copies_corr into these segments
#'
#' @param sample_id sample ID to write into output dataframe
#' @param segments_pre segment file (e.g. from a somatic object)
#' @param segments_post segment file (e.g. from a somatic object)
#' @return A data frame with copy number of overlapping segments in pre and post-event samples
#' @export
relative_CN = function( sample_id, segments_pre, segments_post ){
    #https://blog.dominodatalab.com/genomic-ranges-an-introduction-to-working-with-genomic-data/
    #test_pre = data.frame( chrom=c("chr1", "chr1"), start = c(100, 201), end=c(200,300), copies_corr=c(2.3, 1.8))
    #test_post = data.frame( chrom=c("chr1"), start = c(150), end=c(250), copies_corr=c(5))
    gr_pre = GenomicRanges::makeGRangesFromDataFrame( segments_pre, keep.extra.columns = TRUE )
    gr_post = GenomicRanges::makeGRangesFromDataFrame( segments_post, keep.extra.columns = TRUE )

    # find union of all segments in either pre or post
    gr_segs = GenomicRanges::disjoin( unlist( GenomicRanges::GenomicRangesList( gr_pre, gr_post ) ) )

    # segments may be only in pre or in post.
    # intersect to find the segments in both pre and post
    idx_segs_in_both =
        intersect( S4Vectors::queryHits( GenomicRanges::findOverlaps( gr_segs, gr_post ) ),
                   S4Vectors::queryHits( GenomicRanges::findOverlaps( gr_segs, gr_pre ) ) )
    gr_segs = gr_segs[idx_segs_in_both]
    # find the pre and post copy value for each segment in gr_segs
    copies_pre = gr_pre$copies_corr[ S4Vectors::subjectHits( GenomicRanges::findOverlaps( gr_segs, gr_pre ) ) ]
    copies_post = gr_post$copies_corr[ S4Vectors::subjectHits( GenomicRanges::findOverlaps( gr_segs, gr_post ) ) ]
    relative_CN=data.frame(
        sample_id = rep(sample_id, length(gr_segs)),
        chrom=as.character( GenomeInfoDb::seqnames(gr_segs)  ),
        start = GenomicRanges::start( gr_segs ),
        end = GenomicRanges::end( gr_segs ),
        copies_pre,
        copies_post, stringsAsFactors=FALSE)
    relative_CN
}

#' Constuctor for object that holds metadata about the genome itself
#' @param fn_gene_locs TSV, chrom start end NM (refseq_id) score (ignored) strand
#' @param fn_cyto TSV, Chromosome ChromStart ChromEnd Band Stain
#' @param fn_chrom_lengths TSV, (no header) chrom length
#' @param fn_centromeres TSV, chrom chromStart (pos) chromEnd (pos) chrom_idx
#' @param fn_ensembl TSV, generated from GTF, ensembl_id symbol chrom start end strand gene_type
#' @param fn_exons TSV, symbol (no individual exon identifier) chrom start end strand
#' @return A new Genome object
#' @export
Genome = function( fn_gene_locs, fn_cyto, fn_chrom_lengths, fn_centromeres, fn_ensembl, fn_exons){
    # create genome object for metadata about the genome itself

    genome = list("version"="hg38")
    file_check( "Genome", fn_cyto )
    file_check( "Genome", fn_centromeres )
    file_check( "Genome", fn_exons )
    file_check( "Genome", fn_gene_locs )

    genome = c(genome, list( "cyto_info" = utils::read.table(fn_cyto, stringsAsFactors=FALSE,header=TRUE, sep='\t') ) )
    chrom_lengths = utils::read.table(fn_chrom_lengths, row.names=1, stringsAsFactors=FALSE)
    names(chrom_lengths) = "length"
    chrom_valid = rownames(chrom_lengths)
    genome = c(genome, list( "chrom_lengths" = chrom_lengths ) )
    genome = c(genome, list( "centromeres" = utils::read.table( fn_centromeres, header=TRUE, stringsAsFactors = FALSE) ) )
    genome = c(genome, list( "exons" = utils::read.table( fn_exons, sep='\t', header=TRUE ) ) )
    gene_locs = utils::read.table(fn_gene_locs, header=FALSE, sep='\t', stringsAsFactors=FALSE)
    names(gene_locs) = c("chrom", "start", "end", "NM", "score", "strand")
    gene_locs$symbol = get.split.col(gene_locs$NM, "~", first=TRUE)
    gene_locs$NM = get.split.col(gene_locs$NM, "~", last=TRUE)
    gene_locs = gene_locs[order(gene_locs$symbol),]
    gene_locs = gene_locs[ match.idx( chrom_valid, gene_locs$chrom, allow.multiple.B = TRUE )$idx.B,]
    gene_locs_full = gene_locs

    # check for a few genes that harbor annotations on more than one chromosome
    # one example: IL9R~NM_176786.1-2 (chrY), IL9R~NM_176786.1 (chrX)

    symbols = unique(gene_locs$symbol)
    n_chr=rep( 0, length(symbols))
    for(i in  1:length(n_chr)){
        n_chr[i] = length( unique( gene_locs$chrom[ gene_locs$symbol==symbols[i] ] ))
    }
    gene_locs = gene_locs[ ! (gene_locs$symbol %in% symbols[ which(n_chr>1) ] & gene_locs$chrom=="chrY"),]

    gene_locs = plyr::ddply(gene_locs, "symbol", plyr::summarize, start = min(GenomicRanges::start), end = max(GenomicRanges::end) )
    gene_locs$strand = gene_locs_full$strand[ match.idx( gene_locs_full$symbol, gene_locs$symbol  )$idx.A ]
    gene_locs$chrom = gene_locs_full$chrom[ match.idx( gene_locs_full$symbol, gene_locs$symbol  )$idx.A ]
    chrom_idx = get.split.col( gene_locs$chrom, "chr", last=TRUE)
    chrom_idx[chrom_idx=="X"]=23
    chrom_idx[chrom_idx=="Y"]=24
    chrom_idx = as.numeric(chrom_idx)
    gene_locs$chrom_idx = chrom_idx
    si = GenomeInfoDb::Seqinfo(chrom_valid, seqlengths=chrom_lengths$length, isCircular=rep(FALSE, dim(chrom_lengths)[1]), genome=NA)

    ensembl = utils::read.table(fn_ensembl, row.names=1, stringsAsFactors=FALSE, header=TRUE)
    ensembl = ensembl[order(dimnames(ensembl)[[1]] ),]
    idx_chr_alt = grep( "_alt", ensembl$chrom)
    ensembl = ensembl[ setdiff( 1:dim(ensembl)[1], idx_chr_alt ), ]

    genome = c(genome, "seqinfo" = si)
    genome = c(genome, list( "gene_locs" = gene_locs ) )
    genome = c(genome, list( "gene_locs_full" = gene_locs_full ) )
    genome = c(genome, list( "ensembl" = ensembl ) )

    genome
}


# ----------------------------------
# somatic object functions
# ----------------------------------
#' Constuctor for object that holds somatic data for a sample
#'
#' @param sample_base donor identifier, might not be unique
#' @param sample_id unique sample identifier
#' @param library_type One of {WES,WGS}
#' @export
#' @return A new Somatic object
Somatic = function(sample_base, sample_id, library_type){
    if( ! (library_type=="WES" | library_type=="WGS") ){
        stop("library_type must be one of WES,WGS")
    }
    list( "sample_base" = sample_base,
          "sample_id" = sample_id,
          "library_type" = library_type)
}


# ------------------------------------------------------------------------------------------------------
# Quality Control functions
# ------------------------------------------------------------------------------------------------------

#' parses the Picard alignment summary to obtain the number of paired reads
#' @param fn path to file written by Picard
#' @export
read_Picard_alignmentsummary = function( fn ){
    as = utils::read.table( fn, sep='\t', stringsAsFactors = FALSE, header=TRUE)
    as = as[ which(as$CATEGORY=="PAIR"),]
    as
}

#' parses the BamQC summmary file to obtain n_mapped_paired, n_duplicated, median_insert_size, mean_mapping_quality
#' @param fn path to file written by BamQC
#' @export
read_bamqc_summary = function( fn ){
    con=file(fn, open="r")
    lines=readLines(con)
    n_mapped_paired=NA
    n_duplicated=NA
    median_insert_size=NA
    mean_mapping_quality=NA
    for(i in 1:length( lines ) ){
        line = lines[i]
        if( length(line)>0 ){
            a = strsplit(line, " = ", fixed=TRUE)
            if( !is.na( a[[1]][1] )  ){
                if( stringr::str_trim( a[[1]][1], side="both" ) == "number of mapped paired reads (both in pair)") {
                    x = stringr::str_trim( a[[1]][2], side="both" )
                    n_mapped_paired = as.numeric( as.numeric( paste( Biostrings::strsplit( x, ",")[[1]], collapse="" ) ) )
                }
                if( stringr::str_trim( a[[1]][1], side="both" ) == "number of duplicated reads (flagged)") {
                    x = stringr::str_trim( a[[1]][2], side="both" )
                    n_duplicated = as.numeric( as.numeric( paste( Biostrings::strsplit( x, ",")[[1]], collapse="" ) ) )
                }
                if( stringr::str_trim( a[[1]][1], side="both" ) == "median insert size") {
                    x = stringr::str_trim( a[[1]][2], side="both" )
                    median_insert_size = as.numeric( as.numeric( paste( Biostrings::strsplit( x, ",")[[1]], collapse="" ) ) )
                }
                if( stringr::str_trim( a[[1]][1], side="both" ) == "mean mapping quality") {
                    x = stringr::str_trim( a[[1]][2], side="both" )
                    mean_mapping_quality = as.numeric( as.numeric( paste( Biostrings::strsplit( x, ",")[[1]], collapse="" ) ) )
                }
            }
        }
    }
    close(con)
    data.frame( n_mapped_paired, n_duplicated, median_insert_size, mean_mapping_quality, stringsAsFactors=FALSE)
}

# --------------------------------------------------------------------
# plotting functions
# --------------------------------------------------------------------

#' without putting any data up, create the background for a genome-wide plot of CNA data
#'
#' @param genome Genome object
#' @param ylim ylim for plot
#' @param main title
#' @param ylab y axis label
#' @export
plot_CNA_genome_background = function( genome, ylim=c(0,5), main="", ylab="copies"  ){
    chrom_lengths = genome$chrom_lengths$length
    ss = data.frame(chrom_int=1:24, start=chrom_lengths)

    xs = convert.genome.to.scaled.x( ss$chrom_int, ss$start )
    graphics::plot( -100, -100, pch='.', ylim=ylim, xlim=c(0,max(xs, na.rm=TRUE)),
          col="lightgrey", xaxs="i", las=1, xlab="",
          ylab=ylab, main=main, axes=FALSE )

    graphics::axis(2, las=1)
    for(i in seq(from=ylim[1], to=ylim[2], by=1)){
        graphics::lines(c(0,max(xs,na.rm=TRUE) ), c(i,i), col="#eeeeee")
    }
    graphics::text( mean(c(0, xs[ 1 ] ) ), ylim[1], "1", cex=1)
    for(i in 1:23){
        graphics::lines( c( xs[ i ], xs[ i ]),
               c( ylim[1], ylim[2] ),
               col="lightgrey")
        chrname=i+1
        if(chrname==23){
            chrname = "X"
        }else if(chrname==24){
            chrname = "Y"
        }
        graphics::text( mean(c(xs[i], xs[i+1] ) ), ylim[1], chrname, cex=1)
    }
}

#' Plot a track of CNA data, called after plot_CNA_genome_background()
#'
#' @description
#' Can be called multiple times
#'
#' @param chrom_ints chromosomes expressed as integers, so X = 23 and chr5 = 5
#' @param starts positions to start segments
#' @param ends positions to end segments
#' @param values heights
#' @param color for line
#' @export
plot_CNA_add_segments = function( chrom_ints, starts, ends, values, color="black" ){

    seg_xs_start = convert.genome.to.scaled.x( chrom_ints, starts )
    seg_xs_end = convert.genome.to.scaled.x( chrom_ints, ends)
    for(i in 1:length(seg_xs_start) ){
        graphics::lines( c( seg_xs_start[i], seg_xs_end[i] ),
               c( values[i], values[i]) , col=color, lwd=2)
    }
}


#' Plot a uncorrected bins of copy number data in a genomic region
#' @param genome genome object
#' @param somatic somatic object
#' @param chrom chromosome to start segments
#' @param loc_start position to start bins
#' @param loc_end positions to end bins
#' @param symbols symbols to mark with rectangles at the bottom of the plot
#' @param ylim ylim values for plot
#' @return data frame of bins plotted
plot_CNA_bins=function( genome, somatic, chrom, loc_start, loc_end, symbols, ylim){
    df=copynumber_range( genome, somatic, chrom, loc_start, loc_end )
    #if(chrom=="chrX"){
    #    df$copies = df$copies/2
    #}
    plot(df$start, df$copies, pch=19, cex=0.2, ylim=ylim,
         las=1, ylab="uncorrected copies", xlab="",
         main=paste(somatic$sample_id, "purity",somatic$purity) )
    #text( df$start[1], ylim[1]*1.1,
    #      paste(somatic$sample_id, "purity",somatic$purity), font=2, adj=0)
    for(i in 1:length(symbols)){
        symbol=symbols[i]
        idx = which(genome$gene_locs$symbol==symbol)
        graphics::lines( c( genome$gene_locs$start[idx], genome$gene_locs$end[idx] ),
               c(1.1*ylim[1], 1.1*ylim[1]), lwd=3, col="orange")
        #text( genome$gene_locs$start[idx], 1.05*ylim[1], "AR", font=3, adj=0  )
        if(symbol=="AR"){
            graphics::lines( c( 66920031,66920031 ),ylim, lwd=3, col="#66666666")
        }
    }
    df
}

#' Plot a uncorrected bins of copy number data in a genomic region
#' @param genome genome object
#' @param SO list of somatic objects
#' @param chrom chromosome to start segments
#' @param loc_start position to start bins
#' @param loc_end positions to end bins
#' @param symbols symbols to mark with rectangles at the bottom of the plot
#' @param ylim ylim values for plot
#' @return data frame of bins plotted
plot_CNA_bins_across_samples = function( genome, SO, chrom, loc_start, loc_end, symbols, ylim ){
    odds = seq(from=1, to=length(SO)-1, by=2 )
    evens = odds+1
    graphics::layout(matrix(1:2,2,1))
    graphics::par(mar=c(3,3,2,1))
    for(i in 1:length(odds)){

        df_pre=copynumber_range( genome, SO[[ odds[i] ]], chrom, loc_start, loc_end )
        df_post=copynumber_range( genome, SO[[ evens[i] ]], chrom, loc_start, loc_end )
        if( chrom == "chrX" ){
            df_pre$copies=df_pre$copies/2
            df_post$copies=df_post$copies/2
        }
        max_pre=1
        max_post=1
        for(mm in 1:100){
            if( sum( df_pre$copies > mm, na.rm=TRUE ) > 10 ){
                max_pre=mm
            }
            if( sum( df_post$copies > mm, na.rm=TRUE ) > 10 ){
                max_post=mm
            }
        }
        ylim[2] = max(max_pre, max_post)+1
        df_pre=plot_CNA_bins( genome, SO[[ odds[i] ]], chrom, loc_start, loc_end, symbols, ylim )
        df_post=plot_CNA_bins( genome, SO[[ evens[i] ]], chrom, loc_start, loc_end, symbols, ylim )
    }
}


#' Plot a log10 barplot of SNV frequency
#' @param SO list of somatic objects
#' @param ylim ylim
#' @param x_labels labels for objects on X axis
#' @param denominator divide counts by this value, 3000 would give muts / Mb
#' @param min_alt_count minimum reads for alternate allele to count
#' @export
plot_SNV_frequency = function( SO, ylim=c(0,5), x_labels=NULL, denominator=3000, min_alt_count=20 ){

    lbls = c("0","10","100","1,000","10,000","100,000", "1,000,000")
    snv_counts = rep(0, length(SO))
    for(i in 1:length(SO)){

        alts = as.character( unlist( VariantAnnotation::alt( SO[[i]]$SNV_strelka ) ) )
        As = VariantAnnotation::geno( SO[[i]]$SNV_strelka  )$AU[,,1]
        Cs = VariantAnnotation::geno( SO[[i]]$SNV_strelka  )$CU[,,1]
        Gs = VariantAnnotation::geno( SO[[i]]$SNV_strelka  )$GU[,,1]
        Ts = VariantAnnotation::geno( SO[[i]]$SNV_strelka  )$TU[,,1]
        n_alt = rep(0, dim(As)[1])
        idx_alt_A = which( alts=="A" )
        n_alt[idx_alt_A] = As[idx_alt_A, 2]
        idx_alt_C = which( alts=="C" )
        n_alt[idx_alt_C] = Cs[idx_alt_C, 2]
        idx_alt_G = which( alts=="G" )
        n_alt[idx_alt_G] = Gs[idx_alt_G, 2]
        idx_alt_T = which( alts=="T" )
        n_alt[idx_alt_T] = Ts[idx_alt_T, 2]
        snv_counts[i] = sum(n_alt>=min_alt_count)
    }
    b=graphics::barplot( log10(snv_counts / denominator), col='cornflowerblue',axes=FALSE, xaxs="i",
             ylim=ylim, xlab="")
    ats = seq(from=ylim[1],to=ylim[2],by=1)
    labels = lbls[1:length(ats)]
    graphics::axis(2, labels = labels, at=ats, las=1)
    if( !is.null(x_labels)){
        graphics::axis(1, labels=x_labels, at=b[1:length(x_labels)], las=2)
    }
    snv_counts
}


# ------------------------------------------------------------------------------------------------------
# RNA functions
# ------------------------------------------------------------------------------------------------------

#' parse the output of a STAR mapping call to extract transcript counts and metadata of alignment
#' writes N_unmapped, N_multimapping, N_noFeature, N_ambiguous, and STAR_counts
#' @param fn_path path to STAR output
#' @param somatic incoming somatic object to modify
#' @param genome genome data, used to map transcripts
#' @export
read_STAR_counts = function( fn_path, somatic, genome ){
    cnt = utils::read.table( fn_path, stringsAsFactors = FALSE, header=FALSE, row.names=1)
    names(cnt) = c("total", "strand1", "strand2")
    somatic = c(somatic, list( "STAR_unmapped" = cnt$total[ which(dimnames(cnt)[[1]] == "N_unmapped" ) ] ) )
    somatic = c(somatic, list( "STAR_multimapping" = cnt$total[ which(dimnames(cnt)[[1]] == "N_multimapping" ) ] ) )
    somatic = c(somatic, list( "STAR_nofeature" = cnt$total[ which(dimnames(cnt)[[1]] == "N_noFeature" ) ] ) )
    somatic = c(somatic, list( "STAR_ambiguous" = cnt$total[ which(dimnames(cnt)[[1]] == "N_ambiguous" ) ] ) )
    cnt = cnt[5:dim(cnt)[1],]
    cnt_std = rep(0, dim(cnt)[1])
    m = match.idx( dimnames( genome$ensembl )[[1]],  dimnames(cnt)[[1]] )
    cnt_std[ m$idx.A ] = cnt$total[ m$idx.B ]
    somatic = c(somatic, list( "STAR_counts" = cnt_std ) )
    somatic
}

#' Parse the output of a STAR metadata output to extract mean insert length
#'
#' @description
#' input file has loose text format. Writes STAR_mean_insert_length.
#'
#' @param fn_path path to star metadata
#' @param somatic incoming somatic object to modify
#' @export
read_STAR_mean_insert_length = function( fn_path, somatic ){
    con=file(fn_path, open="r")
    lines=readLines(con)
    for(i in 1:length( lines ) ){
        line = lines[i]
        if( length(line)>0 ){
            a = strsplit(line, "|", fixed=TRUE)
            if( !is.na( a[[1]][1] )  ){
                if( stringr::str_trim( a[[1]][1], side="both" ) == "Average mapped length") {
                    aml = as.numeric( stringr::str_trim( a[[1]][2], side="both" ) )
                }
            }
        }
    }
    close(con)
    somatic = c(somatic, list("STAR_mean_insert_length" = aml ) )
    somatic
}

#' Calculate TPM from feature counts
#'
#' @description
#' Uses effective length calculation and writes STAR_tpm to somatic
#'
#' @param somatic incoming somatic object to modify
#' @param genome genome data, used to map transcripts
#' @export
calculate_STAR_TPM = function( somatic, genome  ){
    counts = somatic$STAR_counts
    meanFragmentLength = somatic$STAR_mean_insert_length
    # this is imperfect because it include intron regions. The GTF file contains exons
    # for each transcript rather than each gene (since genes are not really a thing)
    # so each gene has very large numbers of exons assigned to overlapping transcripts.
    # imperfect solution for now
    featureLength = genome$ensembl$end-genome$ensembl$start
    stopifnot(length(featureLength) == length(counts))

    # Compute effective lengths of features in each library.
    effLen = featureLength - meanFragmentLength + 1

    # Exclude genes with length less than the mean fragment length.
    #idx <- which( featureLength >= meanFragmentLength )
    #counts <- counts[idx,]
    #effLen <- effLen[idx,]
    #featureLength <- featureLength[idx]

    # Process one column at a time.
    rate = log(counts) - log(effLen)
    denom = log(sum(exp(rate), na.rm=TRUE))
    tpm = exp(rate - denom + log(1e6))

    somatic = c(somatic, list( "STAR_tpm" = tpm ) )
    somatic
}

#' Aggregates TPM counts from all samples into a matrix for protein-coding genes
#' defined in genome$ensembl
#' @param genome genome data, used to map transcripts
#' @param SO incoming list of somatic objects
#' @param max_TPM maximum allowed TPM value, exceed this and get set to 0
#' @export
calculate_STAR_TPM_matrix = function(genome, SO, max_TPM = 1e9){
    # generate protein-coding TPM matrix
    #
    idx = which( genome$ensembl$gene_type=="protein_coding" )
    matrix_tpm_full = matrix( 0, nrow=length(idx), ncol=length(SO))
    for(i in 1:length(SO)){
        matrix_tpm_full[, i ] = SO[[i]]$STAR_tpm[idx]
    }
    matrix_tpm_full[is.na(matrix_tpm_full)] = 0
    matrix_tpm_full[ matrix_tpm_full > max_TPM ] = 0
    dimnames(matrix_tpm_full)[[1]] = genome$ensembl$symbol[idx]
    dimnames(matrix_tpm_full)[[2]] = names(SO)

    # contains multiple entries for genes with more than one transcript
    symbols = sort( unique( dimnames(matrix_tpm_full)[[1]] ) )
    matrix_tpm = matrix(NA, ncol=dim(matrix_tpm_full)[2], nrow=length(symbols))
    dimnames(matrix_tpm)[[1]] = symbols
    for(i in 1:dim(matrix_tpm)[1]){
        idx = which( dimnames(matrix_tpm_full)[[1]] == dimnames(matrix_tpm)[[1]][i] )
        if( length(idx)==1 ){
            matrix_tpm[i,] = matrix_tpm_full[idx,]
        }else{
            matrix_tpm[i,] = colSums( matrix_tpm_full[idx,], na.rm=TRUE )
        }
    }
    matrix_tpm
}




# ------------------------------------------------------------------------------------------------------
# Structural variations
# ------------------------------------------------------------------------------------------------------

#' Parse a Manta VCF file to calculate read depth in the WT and alternate
#' @param vcf_f VCF file read by VariantAnnotation package
parse_manta_PR_SR = function( vcf_f ){
    reads = data.frame(
        pr_n_wt=rep(0, length(vcf_f)),
        pr_n_alt=rep(0, length(vcf_f)),
        pr_t_wt=rep(0, length(vcf_f)),
        pr_t_alt=rep(0, length(vcf_f)),
        sr_n_wt=rep(0, length(vcf_f)),
        sr_n_alt=rep(0, length(vcf_f)),
        sr_t_wt=rep(0, length(vcf_f)),
        sr_t_alt=rep(0, length(vcf_f)),
        stringsAsFactors = FALSE
    )
    if( length(vcf_f) > 0 ){
        # had to add this check because a manta vcf because a vcf with no PASS entries threw errors
        prs_n = VariantAnnotation::geno(vcf_f)$PR[,1]
        prs_t = VariantAnnotation::geno(vcf_f)$PR[,2]
        srs_n = VariantAnnotation::geno(vcf_f)$SR[,1]
        srs_t = VariantAnnotation::geno(vcf_f)$SR[,2]

        for(i in 1:length(vcf_f)){
            reads$pr_n_wt[i] = prs_n[[i]][1]
            reads$pr_n_alt[i] = prs_n[[i]][2]
            reads$pr_t_wt[i] = prs_t[[i]][1]
            reads$pr_t_alt[i] = prs_t[[i]][2]
            reads$sr_n_wt[i] = srs_n[[i]][1]
            reads$sr_n_alt[i] = srs_n[[i]][2]
            reads$sr_t_wt[i] = srs_t[[i]][1]
            reads$sr_t_alt[i] = srs_t[[i]][2]
        }
        reads$sr_n_wt[ is.na(reads$sr_n_wt) ] = 0
        reads$sr_n_alt[ is.na(reads$sr_n_alt) ] = 0
        reads$sr_t_wt[ is.na(reads$sr_t_wt) ] = 0
        reads$sr_t_alt[ is.na(reads$sr_t_alt) ] = 0
    }
    reads
}

#' Add Manta output to somatic object

#' @param somatic incoming somatic object to modify, writes n_bnd, n_del, n_dup, n_ins, n_inv, manta
#' where manta is the vcf file read by VariantAnnotation
#' @param fn_manta manta output file
#' @param MIN_READS discard Manta output with fewer than this number of reads in alternate
#' @export
add_manta = function( somatic, fn_manta, MIN_READS=2 ){
    sample_id_tumor = somatic$sample_id
    file_check( "add_manta", fn_manta )
    vcf = VariantAnnotation::readVcf(fn_manta, "hg38")
    vcf_f = vcf[VariantAnnotation::fixed(vcf)$FILTER =="PASS"]

    reads = parse_manta_PR_SR( vcf_f )
    keep = reads$pr_n_alt==0 & reads$sr_n_alt == 0 &
        reads$pr_t_alt + reads$sr_t_alt >= MIN_READS
    vcf_f = vcf_f[keep]
    somatic$n_bnd = sum( VariantAnnotation::info(vcf_f)$SVTYPE=="BND") / 2 # otherwise both ends are counted
    somatic$n_del = sum( VariantAnnotation::info(vcf_f)$SVTYPE=="DEL")
    somatic$n_dup = sum( VariantAnnotation::info(vcf_f)$SVTYPE=="DUP")
    somatic$n_ins = sum( VariantAnnotation::info(vcf_f)$SVTYPE=="INS")
    somatic$n_inv = sum( VariantAnnotation::info(vcf_f)$SVTYPE=="INV")
    somatic$manta = vcf_f
    somatic
}

#' Collate manta data into list_sv, including start/stop of each variant type
#' and the evidence for microhomology at deletions
#' @param somatic incoming somatic object to modify
#' @param nt_for_microhomology number of nucleotides to search on shoulder of deletions
#' @export
calculate_manta_summary = function( somatic, nt_for_microhomology = 15){
    data.frame( IRanges::ranges( somatic$manta ) )$start

    df = data.frame(
        svtype = as.character( VariantAnnotation::info( somatic$manta)$SVTYPE),
        chrom_1 = as.character( GenomeInfoDb::seqnames(somatic$manta )),
        start_1 = data.frame( IRanges::ranges( somatic$manta ) )$start, #start( somatic$manta ),
        end_1 = data.frame( IRanges::ranges( somatic$manta ) )$end, # end( somatic$manta ),
        chrom_2 = as.character( GenomeInfoDb::seqnames(somatic$manta )),
        start_2 = data.frame( IRanges::ranges( somatic$manta ) )$start, #start( somatic$manta ),
        end_2 = data.frame( IRanges::ranges( somatic$manta ) )$end, # end( somatic$manta ),
        filter = as.character( SummarizedExperiment::rowRanges( somatic$manta)$FILTER),
        row.names=names( SummarizedExperiment::rowRanges( somatic$manta)),
        stringsAsFactors=FALSE
    )

    if( dim(df)[1] > 0 ){
        df$sample_id = rep( somatic$sample_id, dim(df)[1])
        df$sv_length = rep(NA, dim(df)[1])
        mate_id = rep(NA, dim(df)[1])
        for(i in 1:length( df$sv_length )){
            ll = VariantAnnotation::info( somatic$manta)$SVLEN[[i]]
            if(  length(ll)>0 ){
                df$sv_length[i]=ll
            }
            mid = VariantAnnotation::info( somatic$manta)$MATEID[[i]]
            if( length(mid)>0 ){
                mate_id[i] = mid
            }
        }
        idx_del = which( df$svtype=="DEL" )
        idx_ins = which( df$svtype=="INS" )
        idx_dup = which( df$svtype=="DUP" )
        idx_inv = which( df$svtype=="INV" )
        idx_bnd = which( df$svtype=="BND" )
        if( length(idx_bnd)>0 ){
            str2 = unlist( VariantAnnotation::alt( somatic$manta )[ idx_bnd ] )
            is_left = stringr::str_detect(string=str2, pattern = stringr::fixed("[") )
            is_right = stringr::str_detect(string=str2, pattern = stringr::fixed("]") )
            chrom_dest = rep( "", length(str2) )
            pos_dest = rep( 0, length(str2) )
            for( x in 1:length( str2 ) ){
                if( is_left[x] ){
                    slug = strsplit( str2[x], "[", fixed=TRUE )[[1]][2]
                }else{
                    slug = strsplit( str2[x], "]", fixed=TRUE )[[1]][2]
                }
                chr_pos = strsplit( slug, ":")[[ 1 ]]
                chrom_dest[x] = chr_pos[ 1 ]
                pos_dest[x] = chr_pos[ 2 ]
            }
            df$chrom_2[ idx_bnd ] = chrom_dest
            df$start_2[ idx_bnd ] = as.numeric( pos_dest )
            df$end_2[ idx_bnd ] = df$start_2[ idx_bnd ]+1
        }
        if( length(idx_del)>0 ){
            df$end_1[ idx_del ] = df$start_1[ idx_del ] + 1
            df$start_2[ idx_del ] = df$start_1[ idx_del ] - df$sv_length[ idx_del ]
            df$end_2[ idx_del ] = df$start_2[ idx_del ] + 1
        }
        if( length(idx_ins)>0 ){
            df$end_1[ idx_ins ] = df$start_1[ idx_ins ] + 1
            df$start_2[ idx_ins ] = df$start_1[ idx_ins ] + df$sv_length[ idx_ins ]
            df$end_2[ idx_ins ] = df$start_2[ idx_ins ] + 1
        }
        if( length(idx_dup)>0 ){
            df$end_1[ idx_dup ] = df$start_1[ idx_dup ] + 1
            df$start_2[ idx_dup ] = df$start_1[ idx_dup ] + df$sv_length[ idx_dup ]
            df$end_2[ idx_dup ] = df$start_2[ idx_dup ] + 1
        }
        if( length( idx_inv ) > 0 ){
            df$end_1[ idx_inv ] = df$start_1[ idx_inv ] + 1
            df$start_2[ idx_inv ] = df$start_1[ idx_inv ] + df$sv_length[ idx_inv ]
            df$end_2[ idx_inv ] = df$start_2[ idx_inv ] + 1
        }
        df$mate_id = mate_id
        df = cbind(df, parse_manta_PR_SR( somatic$manta ) )

        df$n_microhomology = rep(NA, dim( df )[1] )
        for(i in 1:dim( df )[1] ){
            if( df$svtype[i]=="DEL" ){
                df$n_microhomology[i] = identify_microhomology( df$chrom_1[i],
                                                                df$start_1[i],
                                                                df$end_1[i],
                                                                nt_for_microhomology )
            }
        }
    }
    somatic$list_sv = df
    somatic$nt_for_microhomology = nt_for_microhomology
    somatic
}


# ------------------------------------------------------------------------------------------------------
# SNVs
# ------------------------------------------------------------------------------------------------------

#' Parse Strelka output from SNV and Indel text files
#' @param fn_strelka_clinvar_snv path to SNV file
#' @param fn_strelka_clinvar_indel path to indel file
read_SNV_gene_summary_strelka = function( fn_strelka_clinvar_snv, fn_strelka_clinvar_indel ){
    snv = utils::read.table(fn_strelka_clinvar_snv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    has_allele_counts = dim(snv)[2] > 9
    for(i in 1:4){
        snv = cbind( snv, rep(NA, dim(snv)[1]))
    }
    if( has_allele_counts ){
        names(snv) = c("chrom", "pos", "ref", "alt", "gene", "transcript", "effect", "aa", "clinsig",
                       "A_N", "C_N", "G_N", "T_N", "A_T", "C_T", "G_T", "T_T",
                       "DP_indel_ref_N", "DP_indel_alt_N", "DP_indel_ref_T", "DP_indel_alt_T" )
    }else{
        names(snv) = c("chrom", "pos", "ref", "alt", "gene", "transcript", "effect", "aa", "clinsig",
                       "DP_indel_ref_N", "DP_indel_alt_N", "DP_indel_ref_T", "DP_indel_alt_T")
    }
    indel = utils::read.table(fn_strelka_clinvar_indel, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    names(indel) = c("chrom", "pos", "ref", "alt", "gene", "transcript", "effect", "aa", "clinsig",
                     "DP_indel_ref_N", "DP_indel_alt_N", "DP_indel_ref_T", "DP_indel_alt_T")
    if( has_allele_counts ){
        for(i in 1:8){
            indel = cbind( indel, rep(NA, dim(indel)[1]))
        }
        indel = indel[, c(1:9, 14:21, 10:13)]
        names(indel) = names(snv)
    }
    rbind(snv, indel)
}

#' add_SNV_strelka_germline
#'
#' Deprecated; strelka germline calling does not appear reliable in our hands
#' @param somatic incoming somatic object to modify
#' @param fn_strelka_germline path to germline VCF file
#' @export
add_SNV_strelka_germline=function( somatic, fn_strelka_germline ){
    file_check( "add_SNV_strelka_germline", fn_strelka_germline )
    snv = utils::read.table(fn_strelka_germline, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    names(snv) = c("chrom", "pos", "ref", "alt", "gene", "transcript", "effect", "aa", "clinsig",
                   "DP_ref", "DP_alt" )
    somatic = c(somatic, list( "strelka_germline" = snv ) )
    somatic
}

#' Wrapper around readVcf
#' @param fn_vcf path to VCF file
#' @param build default to hg38
#' @export
read_VCF = function(fn_vcf, build="hg38"){
    if( !file.exists( fn_vcf ) ){
        stop( paste("read_VCF cannot find file", fn_vcf) )
    }
    VariantAnnotation::readVcf(fn_vcf, build)
}


#' Add Strelka results to somatic object
#'
#' @description
#' Reads files that have been processed by snpSift extractFields in the annotate_mutation_calls
#' step of the somatic workflow
#'
#' @param somatic incoming somatic object to modify
#' @param fn_strelka_clinvar_snv path to SNVs
#' @param fn_strelka_clinvar_indel path to indels
#' @export
add_SNV_strelka = function( somatic, fn_strelka_clinvar_snv, fn_strelka_clinvar_indel ){
    file_check( "add_SNV_strelka", fn_strelka_clinvar_snv )
    file_check( "add_SNV_strelka", fn_strelka_clinvar_indel )
    somatic = c(somatic, list( "strelka" = read_SNV_gene_summary_strelka( fn_strelka_clinvar_snv, fn_strelka_clinvar_indel ) ) )
    somatic
}

#' Add Mutect results to somatic object
#'
#' @description
#' Reads files that have been processed by snpSift extractFields in the annotate_mutation_calls
#' step of the somatic workflow. Currently only handles Mutect version 1.
#'
#' @param somatic incoming somatic object to modify
#' @param fn_mutect_clinvar_snv path to SNVs
#' @export
add_SNV_mutect = function(somatic, fn_mutect_clinvar_snv ){
    file_check( "add_SNV_mutect", fn_mutect_clinvar_snv )
    mutect = utils::read.table(fn_mutect_clinvar_snv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    names(mutect) = c("chrom", "pos", "ref", "alt", "gene", "transcript", "effect", "aa", "clinsig",
                    "DP_SNV_ref_T", "DP_SNV_alt_T", "DP_SNV_ref_N", "DP_SNV_alt_N")
    somatic = c(somatic, list( "mutect" = mutect ) )

    if(sum( which( names(somatic) =="strelka") )>0 ){
        slug_strelka = paste(somatic$strelka$chrom, somatic$strelka$pos, sep=":")
        slug_mutect = paste(somatic$mutect$chrom, somatic$mutect$pos, sep=":")
        m=match.idx( slug_strelka, slug_mutect )
        df = somatic$strelka[ m$idx.A,]
        idx_alt_a = which( df$alt == "A" )
        idx_alt_c = which( df$alt == "C" )
        idx_alt_g = which( df$alt == "G" )
        idx_alt_t = which( df$alt == "T" )
        df$vaf_obs = rep( NA, dim(df)[1] )
        T_sum = df$A_T + df$C_T + df$G_T + df$T_T
        df$n_alt = rep(NA, dim(df)[1])
        df$n_ref = rep(NA, dim(df)[1])
        df$n_alt[ idx_alt_a ] = df$A_T[ idx_alt_a ]
        df$n_ref[ idx_alt_a ] = T_sum[ idx_alt_a ]-df$n_alt[ idx_alt_a ]

        df$n_alt[ idx_alt_c ] = df$C_T[ idx_alt_c ]
        df$n_ref[ idx_alt_c ] = T_sum[ idx_alt_c ]-df$n_alt[ idx_alt_c ]

        df$n_alt[ idx_alt_g ] = df$G_T[ idx_alt_g ]
        df$n_ref[ idx_alt_g ] = T_sum[ idx_alt_g ]-df$n_alt[ idx_alt_g ]

        df$n_alt[ idx_alt_t ] = df$T_T[ idx_alt_t ]
        df$n_ref[ idx_alt_t ] = T_sum[ idx_alt_t ]-df$n_alt[ idx_alt_t ]
        df$vaf_obs = round( df$n_alt / (df$n_alt+df$n_ref), 4 )
        somatic = c( somatic, list( "strelka_mutect" = df ) )
    }
    somatic
}



#' Calculate number of strelka calls where the alt allele has at least min_depth
#' @param somatic incoming somatic object to modify
#' @param min_depth minimum depth to count a call
#' @export
calculate_strelka_properties = function(somatic, min_depth=5 ){
    if( length( which( names(somatic) == "strelka" ) )== 0 ){
        stop( "somatic does not have strelka data loaded" )
    }
    somatic$n_strelka_pass = 0
    somatic$n_strelka_mindepth = 0

    if( dim(somatic$strelka)[1] > 0 ){
        slug = paste( somatic$strelka$chrom, somatic$strelka$pos, sep=":" )
        idx_keep = match.idx( unique(slug), slug )$idx.B
        strelka_unique = somatic$strelka[idx_keep,]
        n_alt = rep(0, dim(strelka_unique)[1])
        n_ref = n_alt
        n_ref[which( strelka_unique$ref=="A")] = strelka_unique$A_N[ which( strelka_unique$ref=="A") ]
        n_ref[which( strelka_unique$ref=="C")] = strelka_unique$C_N[ which( strelka_unique$ref=="C") ]
        n_ref[which( strelka_unique$ref=="G")] = strelka_unique$G_N[ which( strelka_unique$ref=="G") ]
        n_ref[which( strelka_unique$ref=="T")] = strelka_unique$T_N[ which( strelka_unique$ref=="T") ]
        n_alt[which( strelka_unique$alt=="A")] = strelka_unique$A_T[ which( strelka_unique$alt=="A") ]
        n_alt[which( strelka_unique$alt=="C")] = strelka_unique$C_T[ which( strelka_unique$alt=="C") ]
        n_alt[which( strelka_unique$alt=="G")] = strelka_unique$G_T[ which( strelka_unique$alt=="G") ]
        n_alt[which( strelka_unique$alt=="T")] = strelka_unique$T_T[ which( strelka_unique$alt=="T") ]

        somatic$n_strelka_pass = dim(strelka_unique)[1]
        somatic$n_strelka_mindepth = sum( n_alt >= min_depth, na.rm=TRUE )
    }
    somatic
}


#' Calculate mutation signatures with ReconstructSigs methods
#'
#' Writes out mutsig to somatic object
#'
#' @param somatic incoming somatic object to modify
#' @param tri.counts.method passed through to whichSignatures
#' @param min_depth minimum alt depth for mutation signature analysis
calculate_SNV_signatures = function(somatic, tri.counts.method = "genome", min_depth=10){

    slug = paste( somatic$strelka$chrom, somatic$strelka$pos, sep=":" )
    idx_keep = match.idx( unique(slug), slug )$idx.B
    strelka_unique = somatic$strelka[idx_keep,]
    n_alt = rep(0, dim(strelka_unique)[1])
    n_ref = n_alt
    n_ref[which( strelka_unique$ref=="A")] = strelka_unique$A_N[ which( strelka_unique$ref=="A") ]
    n_ref[which( strelka_unique$ref=="C")] = strelka_unique$C_N[ which( strelka_unique$ref=="C") ]
    n_ref[which( strelka_unique$ref=="G")] = strelka_unique$G_N[ which( strelka_unique$ref=="G") ]
    n_ref[which( strelka_unique$ref=="T")] = strelka_unique$T_N[ which( strelka_unique$ref=="T") ]
    n_alt[which( strelka_unique$alt=="A")] = strelka_unique$A_T[ which( strelka_unique$alt=="A") ]
    n_alt[which( strelka_unique$alt=="C")] = strelka_unique$C_T[ which( strelka_unique$alt=="C") ]
    n_alt[which( strelka_unique$alt=="G")] = strelka_unique$G_T[ which( strelka_unique$alt=="G") ]
    n_alt[which( strelka_unique$alt=="T")] = strelka_unique$T_T[ which( strelka_unique$alt=="T") ]
    muts = data.frame( chrom=strelka_unique$chrom,
                       pos=strelka_unique$pos,
                       ref=strelka_unique$ref,
                       alt=strelka_unique$alt,
                       n_alt, stringsAsFactors=FALSE)
    muts$sample_id = rep(somatic$sample_id, dim(muts)[1])


    sig_input = deconstructSigs::mut.to.sigs.input( muts[muts$n_alt >= min_depth,],
                                   sample.id = "sample_id",
                                   chr="chrom",
                                   pos="pos",
                                   ref="ref",
                                   alt="alt",
                                   bsg=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

    sig=deconstructSigs::whichSignatures(tumor.ref=sig_input, sample.id=somatic$sample_id,
                        signatures.ref = deconstructSigs::signatures.cosmic,
                        associated = c(), signatures.limit = NA,
                        signature.cutoff = 0.02, contexts.needed = TRUE,
                        tri.counts.method = tri.counts.method)
    somatic = c( somatic, list( mutsig = sig ) )
    somatic
}



# ------------------------------------------------------------------------------------------------------
# gene by gene analysis
# ------------------------------------------------------------------------------------------------------
#' Completely characterize one gene in one sample
#'
#' @description
#' uses any information that has been loaded into SO to characterize a symbol in one sample
#' for copy alterations, SNVs, inactivating SVs, and number of predicted functional alleles
#'
#' @param S somatic object, a list of data measurements such as SNVs and CNA data
#' @param genome genome metadata generated by Genome() function
#' @param symbol gene symbol to analyze
#' @param sex M or F
#' @param CNA_max_deepdel threshold for calling deep deltion
#' @param CNA_max_shallowdel threshold for calling shallow deletion
#' @param CNA_min_amp threshold for calling amplification
#' @param SNV_method which SNV summary to use, one of {strelka,strelka_mutect,mutect}
#' @param effects_missense list of effects to count as missense, defaults "missense_variant", "missense_variant&splice_region_variant"
#' @param effects_inactivate list of effects to count as inactivating, defaults "stop_gained", "stop_gained&splice_region_variant", "frameshift_variant"
#' @param effects_splice list of effects to count as splice inactivating, defaults "splice_acceptor_variant&intron_variant", "splice_donor_variant&intron_variant"
#' @export
symbol_summary = function( S, genome, symbol, sex="M",
                           CNA_max_deepdel=0.7, CNA_max_shallowdel=1.25, CNA_min_amp=3,
                           SNV_method="strelka_mutect",
                           effects_missense=c("missense_variant","missense_variant&splice_region_variant"),
                           effects_inactivate=c( "stop_gained", "stop_gained&splice_region_variant", "frameshift_variant"),
                           effects_splice=c("splice_acceptor_variant&intron_variant","splice_donor_variant&intron_variant")){
    ret = list(
        snv_missense = NA,
        snv_inactivate = NA,
        snv_splice = NA,
        snv_germline_missense = NA,
        snv_germline_inactivate = NA,
        snv_germline_splice = NA,
        sv_inactivate = NA,
        DNA_copy_mean = NA,
        DNA_copy_max = NA,
        DNA_copy_min = NA,
        has_DNA_copy_amp = NA,
        has_DNA_copy_deepdel = NA,
        has_DNA_copy_shallowdel = NA,
        has_SV_inactivate = NA,
        has_inactivating_germline = NA,
        tpm = NA,
        allele_count = 2, # expected diploid
        n_inactivate = NA,
        n_missense = NA,
        n_splice = NA
    )
    gene_loc = genome$gene_locs[ genome$gene_locs$symbol==symbol,]

    if( dim(gene_loc )[1] == 0 ){
        stop( paste("symbol not found:",symbol ) )
    }
    gr_gene = GenomicRanges::makeGRangesFromDataFrame(
        data.frame( chrom=gene_loc$chrom,
                    start=gene_loc$start,
                    end=gene_loc$end, strand=gene_loc$strand) )
    if( sex=="M" & (gene_loc$chrom=="chrX" | gene_loc$chrom=="chrY") ){
        ret$allele_count=1
    }

    if( ! sum(names(S)==SNV_method)==1 ){
        stop( paste0( "SNV_method parameter ",SNV_method," is not present in S"))
    }
    if( SNV_method=="strelka"){
        cv = S$strelka[which( S$strelka$gene == symbol ),]
    }else if( SNV_method=="strelka_mutect"){
        cv = S$strelka_mutect[which( S$strelka_mutect$gene == symbol ),]
    }else if( SNV_method=="mutect"){
        cv = S$mutect[which( S$mutect$gene == symbol ),]
    }
    if( sum( names( S ) == SNV_method ) == 1 ){
        idx_missense = match.idx( effects_missense, cv$effect, allow.multiple.B=TRUE)$idx.B
        idx_inactivate = match.idx( effects_inactivate, cv$effect, allow.multiple.B=TRUE)$idx.B
        idx_splice = match.idx( effects_splice, cv$effect, allow.multiple.B=TRUE)$idx.B
        if( length(idx_missense) > 0 ){
            ret$snv_missense = cv[ idx_missense, ]
            ret$n_missense =   length( unique( ret$snv_inactivate$pos ) )
        }else{
            ret$n_missense = 0
        }
        if( length( idx_inactivate )>0 ){
            ret$snv_inactivate = cv[ idx_inactivate,]
            ret$n_inactivate = length( unique( ret$snv_inactivate$pos ) )
        }else{
            ret$n_inactivate = 0
        }
        if( length(idx_splice) > 0 ){
            ret$snv_splice = cv[ idx_splice, ]
            ret$n_splice = length( unique( ret$snv_splice$pos ) )
        }else{
            ret$n_splice = 0
        }
    }

    # germline pathogenic are counted in n_inactivate, n_splice, n_missense
    if( sum( names( S ) == "strelka_germline" ) == 1 ){
        cv = S$strelka_germline[which( S$strelka_germline$gene == symbol ),]

        idx_missense = match.idx( effects_missense, cv$effect, allow.multiple.B=TRUE)
        idx_inactivate = match.idx( effects_inactivate, cv$effect, allow.multiple.B=TRUE)
        idx_splice = match.idx( effects_splice, cv$effect, allow.multiple.B=TRUE)
        ret$snv_germline_missense = cv[ idx_missense, ]
        ret$snv_germline_inactivate = cv[ idx_inactivate,]
        ret$snv_germline_splice = cv[ idx_splice, ]
        ret$n_inactivate = ret$n_inactivate + length( unique( ret$snv_germline_inactivate$pos ) )
        ret$n_missense =   ret$n_missense + length( unique( ret$snv_germline_missense$pos ) )
        ret$n_splice =     ret$n_splice + length( unique( ret$snv_germline_splice$pos ) )
        ret$has_inactivating_germline = dim( ret$snv_germline_inactivate )[1]>0 | dim( ret$snv_germline_splice )[1]>0
    }

    # only report large deletions if they overlap an exon, vs deletions entirely within intron regions
    if( sum( names( S ) == "list_sv" ) == 1 ){
        idx_del = which( S$list_sv$svtype=="DEL" )
        idx_bndinv = which( S$list_sv$svtype=="BND" | S$list_sv$svtype=="INV")
        gr_exons = GenomicRanges::makeGRangesFromDataFrame( genome$exons[ genome$exons$symbol == symbol, ],
                                                                seqnames.field="chrom")
        gr_sv = GenomicRanges::makeGRangesFromDataFrame(
            data.frame( chrom=S$list_sv$chrom_1[idx_del],
                        start=S$list_sv$start_1[idx_del],
                        end=S$list_sv$end_1[idx_del] ) )
        idx_del_in_sv = S4Vectors::subjectHits( GenomicRanges::findOverlaps( gr_exons, gr_sv) )
        idx_bndinv_in_sv = which( (S$list_sv$chrom_1[idx_bndinv]==gene_loc$chrom &
                                       S$list_sv$start_1[idx_bndinv] > gene_loc$start &
                                       S$list_sv$end_1[idx_bndinv] < gene_loc$start ) |
                                      (S$list_sv$chrom_2[idx_bndinv] == gene_loc$chrom &
                                           S$list_sv$start_2[idx_bndinv] > gene_loc$start &
                                           S$list_sv$end_2[idx_bndinv] < gene_loc$start ) )
        ret$sv_inactivate = S$list_sv[ c( idx_del[ idx_del_in_sv ], idx_bndinv[ idx_bndinv_in_sv ] ),]
        ret$has_SV_inactivate = length(idx_del_in_sv)>0 | length(idx_bndinv_in_sv)>0
    }

    if( sum( names( S ) == "segments" ) == 1 ){
        gr_segs = GenomicRanges::makeGRangesFromDataFrame(
            data.frame( chrom=S$segments$chrom,
                        start=S$segments$start,
                        end=S$segments$end ) )
        idx_in_segs = S4Vectors::subjectHits( GenomicRanges::findOverlaps( gr_gene, gr_segs) )
        if( length(idx_in_segs)==0 ){
            print( paste( "No segment data for", symbol ) )
        }else{
            cc = S$segments$copies_corr[ idx_in_segs ]
            ret$DNA_copy_mean = mean( cc )
            ret$DNA_copy_max = max( cc )
            ret$DNA_copy_min = min( cc )
            ret$has_DNA_copy_deepdel = ret$DNA_copy_min <= CNA_max_deepdel
            ret$has_DNA_copy_shallowdel = !ret$has_DNA_copy_deepdel & ret$DNA_copy_min <= CNA_max_shallowdel
            ret$has_DNA_copy_amp = ret$DNA_copy_max >= CNA_min_amp
        }
    }

    if( sum( names( S ) == "STAR_tpm" ) == 1 ){
        ret$tpm = S$STAR_tpm[ which( dimnames(S$STAR_tpm)[[1]] == symbol ) ]
    }

    # amplified CNA means a higher total allele count
    if( !is.na(ret$has_DNA_copy_amp) && ret$has_DNA_copy_amp ){
        ret$allele_count = round( ret$DNA_copy_max )
    }
    # deep deletions are automatically zero alleles
    if( !is.na(ret$has_DNA_copy_deepdel) && ret$has_DNA_copy_deepdel ){
        ret$allele_count=0
    }
    # subtract alleles for inactivation (either stop or frameshift) or splices
    ret$allele_count = ret$allele_count - ret$n_inactivate - ret$n_splice
    if( !is.na(ret$has_SV_inactivate) && ret$has_SV_inactivate ){
        ret$allele_count=ret$allele_count-1
    }

    if( !is.na(ret$has_DNA_copy_shallowdel) && ret$has_DNA_copy_shallowdel ){
        ret$allele_count = ret$allele_count-1
    }

    # lower bound on allele count has to be zero
    if( ret$allele_count < 0 ){
        ret$allele_count=0
    }
    ret
}

#' Generate a symbol by sample grid illustrating somatic alterations
#'
#' @param SO list of somatic objects
#' @param genome genome metadata generated by Genome() function
#' @param symbol_list vector of symbols (in genome) to plot after a call to symbol_summary
#' @param sample_ids names of elements in SO to plot
#' @param block.height grid element height
#' @param block.width grid element width
#' @param cex.y cex Y axis
#' @param show_x_axis show the X axis
#' @param show_y_labels show the Y labels
#' @param sex one of {M,F}, defaults to M
#' @param CNA_amp_threshold threshold to call amplification, default 4
#' @param CNA_deepdel_threshold threshold to call deep deletion, default 0.5
#' @param missense_whitelist only show missense mutations from symbols on this list
#' @param germline_whitelist only show germline alterations from symbols on this list
#' @param LOH_list dataframe with columns symbol,sample_id. Show LOH for these specific cases and no others.
plot_symbol_summary_grid = function( SO, genome, symbol_list, sample_ids,
                                     block.height=6, block.width=10,
                                     cex.y = 1, show_x_axis = TRUE,
                                     show_y_labels = TRUE,
                                     sex="M", CNA_amp_threshold = 4,
                                     CNA_deepdel_threshold=0.5,
                                     missense_whitelist=NA,
                                     germline_whitelist=NA,
                                     LOH_list=NA){

    col.noevent =   "#eeeeee"
    col.CN.amp =    "#e31a1c" # "red"
    col.CN.single = "#a6cee3" # "cornflowerblue"
    col.CN.double = "#1f78b4" # "darkblue"
    col.CN.LOH =    "darkgrey"
    col.germline =  "black"
    col.svbroken =  "#ff7f00" # "orange"
    col.missense =  "#33a02c" # "chartreuse4"
    col.inactive = "#6a3d9a" # "saddlebrown"
    col.missense.active="yellow"
    col.fusion =    "#fb9a99" # "purple"
    col.svactivated="hotpink"

    spacer_labels=0

    n.rows = length(symbol_list) - spacer_labels
    n.cols = length(sample_ids)
    total.width =  ( block.width * n.cols)
    total.height = ( block.height * n.rows )

    # accumulate summary statistics as we draw
    has_alt_sv = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    has_alt_cna = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    has_alt_mut = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    n.alterations = rep(0, length( symbol_list) )

    graphics::plot(0, 0, col = col.noevent, xlim=c(0,total.width), ylim=c(0,total.height),
         axes=F, xlab="", ylab="", bg=col.noevent, xaxs="i", yaxs="i")
    graphics::rect(0, 0, total.width, total.height, col=col.noevent)
    cur.y = total.height
    symbol_ctr=1
    xlab_locs=rep(0, n.cols)
    ylab_locs=rep(0, length(symbol_list) )

    for(rr in 1:n.rows){
        this.x.left = 0
        for(i in 1:n.cols){
            cc = i
            this.x.right = this.x.left + block.width
            this.x.right.halfway = this.x.left + ( block.width / 2 )
            this.y.top = cur.y
            this.y.mid = cur.y - ( ( block.height)/2 )
            this.y.bot = cur.y - block.height
            this.y.upperthird = cur.y - ( ( block.height ) / 3 )
            this.y.lowerthird = cur.y - block.height + ( ( block.height ) / 3 )

            graphics::rect( this.x.left, cur.y - block.height,
                  this.x.right, cur.y, col=col.noevent, border=NA)

            col.top = col.noevent
            col.bot = col.noevent

            idx_s = which( names(SO)==sample_ids[cc] )
            ss = symbol_summary( SO[[ idx_s ]], genome, symbol_list[rr], sex="M", CNA_min_amp = CNA_amp_threshold )
            colors = c()
            xlab_locs[cc] = mean( c(this.x.left, this.x.right) )

            show_CNA_deep =    !is.na( ss$has_DNA_copy_deepdel) && ss$has_DNA_copy_deepdel
            show_CNA_shallow = !is.na( ss$has_DNA_copy_shallowdel) && ss$has_DNA_copy_shallowdel
            show_CNA_amp =     !is.na( ss$has_DNA_copy_amp) && ss$has_DNA_copy_amp
            if( is.na( missense_whitelist )[1] ){
                show_SNV_missense = dim( ss$snv_missense )[1] > 0
            }else{
                show_SNV_missense = FALSE
                for(i in 1:dim(missense_whitelist)[1]){
                    if( sum( ss$snv_missense$pos == missense_whitelist$pos[i] &
                               ss$snv_missense$alt == missense_whitelist$alt[i] ) > 0 ){
                        show_SNV_missense=TRUE
                    }
                }
            }
            show_SNV_inact =    dim( ss$snv_inactivate )[1] > 0
            show_SNV_splice =   dim( ss$snv_splice )[1] > 0
            show_SV_inactivate = !is.na( ss$has_SV_inactivate ) && ss$has_SV_inactivate
            #show_SV_ACT = !is.na( ss$has_fusion) && ss$has_fusion
            #show_mis_INACT =  alleles$inactivating_missense[ cc ]
            if( is.na( germline_whitelist)[1] ){
                show_germ = !is.na( ss$has_inactivating_germline) && ss$has_inactivating_germline
            }else{
                show_germ = !is.na( ss$has_inactivating_germline) && ss$has_inactivating_germline &&
                            sum( germline_whitelist$sample_id==sample_ids[cc] &
                                 germline_whitelist$symbol_id==symbol_list[rr] ) > 0
            }
            #show_mis_ACT = alleles$activating_missense[cc]
            if( !is.na( LOH_list )[1] ){
                show_LOH =  sum( LOH_list$symbol==symbol_list[rr] &
                                          LOH_list$sample_id==sample_ids[cc])>0
            }else{
                show_LOH=FALSE
            }
            if( show_SV_inactivate & !show_CNA_deep ){
                colors = c(colors, col.svbroken)
                has_alt_sv[rr,cc]=TRUE
            }
            if( show_germ){
                colors = c(colors, col.germline )
                has_alt_mut[rr,cc]=TRUE
            }
            if( show_CNA_deep ){
                colors = c(colors, col.CN.double )
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_CNA_shallow ){
                colors = c(colors, col.CN.single )
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_LOH ){
                colors = c(colors, col.CN.LOH )
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_CNA_amp ){
                colors = c(colors, col.CN.amp )
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_SNV_missense ){
                colors = c(colors, col.missense )
                has_alt_mut[rr,cc]=TRUE
            }
            if( show_SNV_inact | show_SNV_splice ){
                colors = c(colors, col.inactive )
                has_alt_mut[rr,cc]=TRUE
            }
            if( length(colors)>0 ){
                n.alterations[rr]=n.alterations[rr]+1
            }
            if( length(colors)==0 ){
                colors = c(col.noevent, col.noevent)
            }else if( length(colors)==1){
                if( colors[1]==col.CN.double ){
                    colors = c(col.CN.double,col.CN.double) # draw deep dels as two events
                }else{
                    colors = c(colors, col.noevent)
                }
            }
            if( length( colors) == 2 ){
                graphics::rect( this.x.left, this.y.top, this.x.right, this.y.mid, col=colors[1], border=NA)
                graphics::rect( this.x.left, this.y.mid, this.x.right, this.y.bot, col=colors[2], border=NA)
            }else{
                graphics::rect( this.x.left, this.y.top, this.x.right, this.y.upperthird, col=colors[1], border=NA)
                graphics::rect( this.x.left, this.y.upperthird, this.x.right, this.y.lowerthird, col=colors[2], border=NA)
                graphics::rect( this.x.left, this.y.lowerthird, this.x.right, this.y.bot, col=colors[3], border=NA)
            }
            this.x.left = this.x.left + block.width
        }
        cur.y = cur.y - block.height
        symbol_ctr=symbol_ctr+1
    }
    for( cc in 1:n.cols ){
        line_loc = (cc-1)*block.width
        graphics::lines( c(line_loc,line_loc  ), c(0, total.height), col="white")
    }
    if( spacer_labels==0 ){
        for( rr in seq(from=1, to=n.rows, by=1 ) ){
            line_loc = (rr-1) * block.height
            graphics::abline( line_loc, 0, col="black")
        }
    }
    if( spacer_labels==1 ){
        symbol_list = symbol_list[1]
        n.alterations = n.alterations[1]
    }
    results = list( cna=has_alt_cna,
                    mut=has_alt_mut,
                    sv=has_alt_sv,
                    n_symbols=length(symbol_list),
                    n_samples_altered = n.alterations)
    ylab_locs = 0 : ( length ( symbol_list) - 1 ) * block.height + ( block.height / 2 )
    ylab_locs = ylab_locs[ length( ylab_locs ):1 ]
    if( ! show_y_labels ){
        symbol_list=rep("", length(symbol_list))
    }
    graphics::axis(2, at=ylab_locs, labels=symbol_list, las=2, cex.axis=cex.y,
         tick=FALSE, hadj=1, line=-0.5, font=3)
    if( show_x_axis ){
        graphics::axis(1, at=xlab_locs, labels = sample_ids, cex.axis=0.75, las=2 )
    }
    results
}

#' Generate circos plot of structural and copy number alterations
#'
#' @param SO list of somatic objects
#' @param sample_id sample to plot, must be name of element in SO
#' @param chroms vector of chromosomes to plot
#' @param xlims x limits on chromosomes; 2 element vector if 1, matrix if >1 chromosome
#' @param ylims common y limits, 2 element vector
#' @param min_SV_coverage restrict SV to those with combined ALT at least this depth
#' @param horizontal_axis_lines_at where to draw the horizontal grid lines
#' @param sv_manual data frame of additional SV to plot
#' @param main text to draw at center, NA if leave blank
plot_circlize = function( SO,
                          sample_id,
                          chroms,
                          xlims,
                          ylims,
                          min_SV_coverage=10,
                          horizontal_axis_lines_at,
                          sv_manual=data.frame(),
                          main=NA){
    if( sum( names(SO)==sample_id)==0 ){
        stop(paste("sample",sample_id,"not found in list of somatic data"))
    }
    idx_sample = which(names(SO)==sample_id)
    idx_sv = which(  SO[[ idx_sample ]]$list_sv$chrom_1 %in% chroms &
                         SO[[ idx_sample ]]$list_sv$chrom_2 %in% chroms)
    sv = SO[[ idx_sample ]]$list_sv[ idx_sv,]
    sv =  sv[ sv$pr_t_alt + sv$sr_t_alt >= min_SV_coverage,]
    idx_cna = which( SO[[ idx_sample ]]$segments$chrom %in% chroms )
    cna = SO[[idx_sample]]$segments[idx_cna,]
    # cut CNA down to regions that we'll plot; avoids circlize wraparound bug
    for( i in 1:length(chroms) ) {
        if( is.matrix(xlims) ){  xmin = xlims[i,1] ; xmax = xlims[i,2]
        }else{                   xmin = xlims[1];    xmax = xlims[2]  }
        keep = cna$chrom != chroms[i] |
            (cna$chrom == chroms[i] & (
                (cna$start<=xmin & cna$end >= xmin ) |
                    (cna$start<=xmax & cna$end >= xmax ) |
                    (cna$start>xmin & cna$end <= xmax )
            )
            )
        cna = cna[keep,]
        cna$start[ cna$chrom==chroms[i] & cna$start<xmin] = xmin
        cna$end[ cna$chrom==chroms[i] & cna$end>xmax] = xmax
    }

    circlize::circos.clear()
    circlize::circos.initialize(chroms, xlim=xlims )
    circlize::circos.track(sectors=chroms, ylim = ylims )

    # vertical axis for CNA
    for( i in 1:length(chroms)){
        circlize::circos.yaxis(side="right",
                     sector.index = chroms[i],
                     at = horizontal_axis_lines_at )
    }

    # genomic axis and tick marks
    for(i in 1:length(chroms)){
        if( is.matrix(xlims) ){  xmin = xlims[i,1] ; xmax = xlims[i,2]
        }else{                   xmin = xlims[1];    xmax = xlims[2]  }
        circlize::circos.genomicAxis( h = "top",
                            major.at = seq(from=xmin,to=xmax,by=1000000),
                            labels = round( seq(from=xmin,to=xmax,by=1000000) / 1000000 ),
                            sector.index = chroms[i],
                            labels.cex = 1)
    }

    # CNA
    for(i in 1:length(chroms)){
        circlize::circos.lines( x = c(rbind( cna$start[cna$chrom==chroms[i]],
                                   cna$end[cna$chrom==chroms[i]] ) ),
                      y = rep( cna$copyNumber[cna$chrom==chroms[i] ], each=2  ) ,
                      sector.index = chroms[i],
                      track.index = 1,
                      area=TRUE,
                      col="#feb24c33",border="red")
    }

    # axis lines
    for( cc in 1:length(chroms) ){
        for(i in horizontal_axis_lines_at ){
            if( is.matrix(xlims) ){  xmin = xlims[cc,1] ; xmax = xlims[cc,2]
            }else{                   xmin = xlims[1];    xmax = xlims[2]  }
            circlize:: circos.lines( x = c(xmin,xmax),
                          y = c(i,i),
                          sector.index = chroms[cc],
                          track.index = 1,
                          area=FALSE,
                          col="lightgrey")
        }
    }
    # SVs
    if( dim(sv)[1] > 0 ){
        for(i in 1:dim(sv)[1]){
            if( sv$svtype[i]=="BND" ){
                circlize::circos.link(sector.index1=sv$chrom_1[i],
                            point1=c(sv$start_1[i]),
                            sector.index2=sv$chrom_2[i],
                            point2=c(sv$end_1[i]),
                            directional = 1)
            }
        }
        for(i in 1:dim(sv)[1]){
            if( sv$svtype[i]=="DEL" ){
                circlize::circos.link(sector.index1=sv$chrom_1[i],
                            point1=c(sv$start_1[i]),
                            sector.index2=sv$chrom_2[i],
                            point2=c(sv$end_2[i]), col="blue")
            }
        }
        for(i in 1:dim(sv)[1]){
            if( sv$svtype[i]=="DUP" ){
                circlize::circos.link(sector.index1=sv$chrom_1[i],
                            point1=c(sv$start_1[i]),
                            sector.index2=sv$chrom_2[i],
                            point2=c(sv$end_2[i]), col="#2ca25f", lwd=2)
            }
        }
    }
    for(i in 1:length(chroms)){
        if( is.matrix(xlims)){
            mids = rowMeans(xlims)[i]
        }else{
            mids = mean( xlims )
        }
        circlize::circos.text(x = mids, y=mean(ylims),
                    labels = chroms[i], sector.index = chroms[i], track.index = 1, font=2, niceFacing = TRUE )
    }

    # manually added SV (e.g. not present in manta calls)
    if( dim(sv_manual)[1]>0 ){
        for(i in 1:dim(sv_manual)[1]){
            if( sv_manual$svtype[i]=="BND"){
                circlize::circos.link(sector.index1 = sv_manual$chrom1[i], point1 = sv_manual$point1[i],
                            sector.index2 = sv_manual$chrom2[i], point2 = sv_manual$point2[i],
                            directional=1)
            }else{
                circlize::circos.link(sector.index1 = sv_manual$chrom1[i], point1 = sv_manual$point1[i],
                            sector.index2 = sv_manual$chrom2[i], point2 = sv_manual$point2[i],
                            col="blue")
            }
        }
    }
    if(!is.na(main)){
        graphics::text(0,0,main)  # label in center
    }
    list( sv=sv, cna=cna )
}
