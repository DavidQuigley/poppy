# CNA functions
# ------------------------------------------------------------------------------------------------------

#' Read in a segment file created by CopyCat
#'
#' Output has chrom, start, end, segments, copies, chrom_int
#' where chrom_in is the integer version of chromosomes (X = 23)
#' @param fn_path path to segment file
#' @export
read_CNA_segments_copycat = function( fn_path ){
    s = utils::read.table( fn_path, header=FALSE, sep='\t', stringsAsFactors=FALSE )
    names(s) = c("chrom", "start", "end", "segments", "copies")
    s$chrom_int = get.split.col( s$chrom, "chr", last=TRUE)
    s$chrom_int[s$chrom_int=="X"] = 23
    s$chrom_int[s$chrom_int=="Y"] = 24
    s$chrom_int = as.numeric( s$chrom_int )
    s$copies=round(s$copies, 3)
    s
}

#' Read in a segment file created by CNVkit
#'
#' Output has chrom, start, end, log2, depth, weight, copies, chrom_int
#' where chrom_in is the integer version of chromosomes (X = 23)
#' and copies is converted from log2 using knowledge of whether sample is M or F
#'
#' @param fn_path path to segment file
#' @param sex M or F, used to calculate copies
#' @export
read_CNA_segments_cnvkit = function( fn_path, sex="M" ){
    s = utils::read.table( fn_path, header=TRUE, sep='\t', stringsAsFactors=FALSE )
    s = s[, c(1,2,3,5,6,7)]
    names(s) = c("chrom", "start", "end", "log2", "depth", "weight")
    is_diploid = rep( TRUE, dim(s)[1] )
    is_diploid[ which( s$chrom=="chrY" ) ] = FALSE
    if( sex=="M" ){
        is_diploid[  s$chrom=="chrX" ] = FALSE
    }
    s$copies = s$log2
    s$copies[ is_diploid ] = 2 * ( 2 ^ s$copies[ is_diploid ] )
    s$copies[ !is_diploid ] = ( 2 ^ s$copies[ !is_diploid ] )
    s$chrom_int = get.split.col( s$chrom, "chr", last=TRUE)
    s$chrom_int[s$chrom_int=="X"] = 23
    s$chrom_int[s$chrom_int=="Y"] = 24
    s$chrom_int = as.numeric( s$chrom_int )
    s$copies=round(s$copies, 3)
    s
}


#' Read in binned copy number data created by copycat
#'
#' Output is chrom, start, end, copies
#' Infers interval from the gaps between bins.
#' Note that copies doesn't know whether a chromosome is diploid or not.
#'
#' @param fn_path path to binned data
#' @export
read_CNA_rdbins_copycat = function(fn_path){

    # TODO: should update to adjust for diploid status to be consistent with read_CNA_rdbins_cnvkit
    #       would need to adjust downstream as well
    s = utils::read.table( fn_path, header=FALSE, sep='\t', stringsAsFactors=FALSE )
    names(s) = c("chrom", "start", "copies")
    interval = (s$start[2] - s$start[1])-1
    s$end = s$start + interval
    s = s[,c("chrom","start","end","copies")]
    s$chrom_int = get.split.col( s$chrom, "chr", last=TRUE)
    s$chrom_int[s$chrom_int=="X"] = 23
    s$chrom_int[s$chrom_int=="Y"] = 24
    s$chrom_int = as.numeric( s$chrom_int )
    s$copies=round(s$copies, 3)
    s
}

#' Read in binned copy number data created by CNVkit
#'
#' Output is chrom, start, end, copies
#' Infers interval from the gaps between bins.
#'
#' @param fn_path path to binned data
#' @param sex M or F, used to adjust copies for diploid status
#' @export
read_CNA_rdbins_cnvkit = function( fn_path, sex="M" ){
    s = utils::read.table( fn_path, header=TRUE, sep='\t', stringsAsFactors=FALSE )
    names( s ) = c("chrom",	"start", "end",	"gene", "log2", "depth", "weight")
    is_diploid = rep( TRUE, dim(s)[1] )
    is_diploid[ s$chrom=="chrY" ] = FALSE
    if( sex=="M" ){
        is_diploid[  s$chrom=="chrX" ] = FALSE
    }
    s$copies = s$log2
    s$copies[ is_diploid ] = 2 * ( 2 ^ s$copies[is_diploid] )
    s$copies[ !is_diploid ] = ( 2 ^ s$copies[!is_diploid] )
    s$copies=round(s$copies, 3)
    s = s[,c("chrom", "start", "end", "copies", "gene", "log2", "depth", "weight") ]
    s
}


#' Filter segments for overlap with centromere and minimum size to cut down noise
#'
#' @param genome genome object created by Genome(), contains centromere locations
#' @param somatic somatic object containing segments
#' @param max_percent_centromere_overlap eliminate segments where total overlap exceeds this percent, 0-1
#' @param min_segment_percent minimal size of a segment expressed as a percent of the chromosome size
#' @export
adjust_CNA_segments_filter = function( genome, somatic,
                                       max_percent_centromere_overlap=.20,
                                       min_segment_percent=0.00005 ){

    # record percent for each segment
    somatic$segments_raw$chrom_percent = rep(0, dim(somatic$segments_raw)[1])
    chroms = sort(unique(somatic$segments_raw$chrom))
    for(i in 1:length(chroms)){
        idx = which( somatic$segments_raw == chroms[i] )
        seg_chrom = somatic$segments_raw[idx,]
        somatic$segments_raw$chrom_percent[idx] = ( seg_chrom$end - seg_chrom$start ) / max( seg_chrom$end )
    }
    somatic$segments_raw$chrom_percent = signif( somatic$segments_raw$chrom_percent, 3)
    # if more than max_percent_centromere_overlap of a segment overlaps with a centromere,
    # eliminate it.
    somatic$segments_raw$centro_percent = rep(0, dim(somatic$segments_raw)[1])
    gr_seg = GenomicRanges::makeGRangesFromDataFrame(somatic$segments_raw)
    gr_cent = GenomicRanges::makeGRangesFromDataFrame(genome$centromeres)
    hits = GenomicRanges::findOverlaps( gr_seg, gr_cent )
    overlaps = IRanges::pintersect(gr_seg[S4Vectors::queryHits(hits)], gr_cent[S4Vectors::subjectHits(hits)])
    percentOverlap = BiocGenerics::width(overlaps) / BiocGenerics::width(gr_seg[S4Vectors::queryHits(hits)])
    somatic$segments_raw$centro_percent[ S4Vectors::queryHits(hits) ] = percentOverlap
    somatic$segments_raw$centro_percent=round(somatic$segments_raw$centro_percent,2)

    segments = somatic$segments_raw
    segments = segments[ segments$chrom_percent >= min_segment_percent,]
    segments = segments[ segments$centro_percent < max_percent_centromere_overlap,]

    c(somatic, list( "segments" = segments) )
}

#' used to bin the genome into bins of size bin_size
#'
#' the last bin may be smaller than bin_size
#' returns data frame: chrom, start, end
#'
#' @param genome object created by Genome, source of chromosome sizes
#' @param bin_size nucleotides per bin
calculate_genome_bins = function( genome, bin_size ){
    df_bins = data.frame()
    for( chrom in 1:24 ){
        cl = genome$chrom_lengths$length[chrom]
        starts = seq(from=1, to=cl, by=bin_size)
        if( starts[ length(starts) ] > cl ){
            starts = starts[1: ( length(starts) - 1 ) ]
        }
        ends = starts[ 2:length(starts)]-1
        ends = c(ends, cl)
        chroms = rep(dimnames(genome$chrom_lengths)[[1]][chrom], length(ends))
        df_bins = rbind(df_bins,data.frame(chrom=chroms, start=starts, end=ends,
                                           stringsAsFactors=FALSE) )
    }
    df_bins
}

#' Given a bin size, converts genome-wide segments into values per bin
#'
#' Returns a matrix with one row per bin and one column per somatic object in SO
#' @param SO list of somatic objects, final matrix will include all
#' @param genome object created by Genome(), used for chromosome sizes
#' @param bin_size size of each bin, 3000000 creates approximately 1000 bins
#' @export
calculate_CNA_matrix = function(SO, genome, bin_size=3000000 ){

    df_bins = calculate_genome_bins( genome, bin_size )
    gr_bins = GenomicRanges::makeGRangesFromDataFrame( df_bins, seqnames.field = "chrom")
    matrix_CNA = matrix( NA, nrow=length(gr_bins), ncol=length(SO))

    for(s in 1:length(SO)){
        gr_segs = GenomicRanges::makeGRangesFromDataFrame( SO[[s]]$segments,
                                            seqnames.field = "chrom",
                                            keep.extra.columns = TRUE)
        ol = IRanges::findOverlaps( gr_bins, gr_segs)

        ol_bins = S4Vectors::queryHits( ol )
        ol_segs = S4Vectors::subjectHits( ol )
        for(i in 1: dim(df_bins)[1] ){
            idx = which(ol_bins==i)
            if( length(idx)>0 ){
                seg_lengths = rep(0, length(idx))
                seg_vals = rep(0, length(idx))
                for( j in 1:length(idx)){
                    idx_segment = ol_segs[ idx[j] ]
                    seg_vals[j] = gr_segs[ idx_segment ]$copies
                    seg_start =  SO[[s]]$segments$start[idx_segment]
                    seg_end =  SO[[s]]$segments$end[idx_segment]
                    if( seg_start < df_bins$start[i] ){
                        if( seg_end >= df_bins$end[i] ){
                            seg_lengths[j] = bin_size # completely covers
                        }else{
                            seg_lengths[j] = seg_end - df_bins$start[i] # covers left border not right
                        }
                    }else{
                        if( seg_end < df_bins$end[i] ){
                            seg_lengths[j] = seg_end - seg_start #  interior
                        }else{
                            seg_lengths[j] = df_bins$end[i] - seg_start # covers right border not left
                        }
                    }
                    matrix_CNA[i, s] = stats::weighted.mean(seg_vals, w=seg_lengths / sum(seg_lengths))
                }
            }
        }
    }
    dimnames(matrix_CNA)[[1]] = paste(df_bins$chrom,":",df_bins$start,sep='')
    dimnames(matrix_CNA)[[2]] = names(SO)
    matrix_CNA
}



#' Wrapper function to add CNA segments, filter, and optionally adjust for ploidy/purity using CLONET
#' @param genome object created by call to Genome()
#' @param somatic somatic data object to modify
#' @param fn_segs path to segment file
#' @param fn_bins path to binned copy numbers
#' @param fn_pileup_normal pileup file created by ASEQ program
#' @param fn_pileup_tumor pileup file created by ASEQ program
#' @param min_segment_percent minimal size of a segment expressed as a percent of the chromosome size
#' @param max_percent_centromere_overlap eliminate segments where total overlap exceeds this percent, 0-1
#' @param purity_override use this value instead of estimating purity
#' @param ploidy_override use this value instead of estimating ploidy
#' @param min_required_SNPs minimum number of SNPs in a segment to make a call
#' @param sex M or F
#' @param n_cores allow for multithreaded calculations
#' @param experiment_type WGS or WES
#' @param verbose report progress to stdout, default FALSE
#' @export
somatic_add_CNA_cnvkit = function( genome, somatic, fn_segs, fn_bins, fn_pileup_normal, fn_pileup_tumor,
                            min_segment_percent = 0.005,
                            max_percent_centromere_overlap = 0.2,
                            purity_override = NA,
                            ploidy_override = NA,
                            min_required_SNPs = 10,
                            sex="M",
                            n_cores = 1,
                            experiment_type = "WGS",
                            verbose=FALSE ){

    file_check( "somatic_add_CNA", fn_pileup_normal )
    file_check( "somatic_add_CNA", fn_pileup_tumor )
    file_check( "somatic_add_CNA", fn_segs )
    file_check( "somatic_add_CNA", fn_bins )
    CNA_method = "cnvkit"
    sample_base = somatic$sample_base
    sample_id_tumor = somatic$sample_id
    pileup_tumor <- utils::read.table(fn_pileup_tumor, header = TRUE, stringsAsFactors = FALSE)
    pileup_normal <- utils::read.table(fn_pileup_normal, header = TRUE, stringsAsFactors = FALSE)

    somatic = c( somatic, list( "segments_raw" = read_CNA_segments_cnvkit( fn_segs ) ) )
    somatic = c( somatic, list( "rdbins" = read_CNA_rdbins_cnvkit( fn_bins ) ) )

    somatic = adjust_CNA_segments_filter(genome, somatic, min_segment_percent = min_segment_percent )
    somatic
}

#' Wrapper function to add CNA segments, filter, and adjust for ploidy/purity using CLONET
#' @param genome object created by call to Genome()
#' @param somatic somatic data object to modify
#' @param fn_segs path to segment file
#' @param fn_bins path to binned copy numbers
#' @param fn_pileup_normal pileup file created by ASEQ program
#' @param fn_pileup_tumor pileup file created by ASEQ program
#' @param min_segment_percent minimal size of a segment expressed as a percent of the chromosome size
#' @param max_percent_centromere_overlap eliminate segments where total overlap exceeds this percent, 0-1
#' @param purity_override use this value instead of estimating purity
#' @param ploidy_override use this value instead of estimating ploidy
#' @param min_required_SNPs minimum number of SNPs in a segment to make a call
#' @param sex M or F
#' @param n_cores allow for multithreaded calculations
#' @param verbose TRUE or FALSE
#' @export
somatic_add_CNA_copycat = function( genome, somatic, fn_segs, fn_bins, fn_pileup_normal, fn_pileup_tumor,
                            min_segment_percent = 0.005,
                            max_percent_centromere_overlap = 0.2,
                            purity_override = NA,
                            ploidy_override = NA,
                            min_required_SNPs = 10,
                            sex="M",
                            n_cores = 1,
                            verbose=FALSE ){

    file_check( "somatic_add_CNA", fn_pileup_normal )
    file_check( "somatic_add_CNA", fn_pileup_tumor )
    file_check( "somatic_add_CNA", fn_segs )
    file_check( "somatic_add_CNA", fn_bins )
    CNA_method = "copycat"
    sample_base = somatic$sample_base
    sample_id_tumor = somatic$sample_id
    pileup_tumor <- utils::read.table(fn_pileup_tumor, header = TRUE, stringsAsFactors = FALSE)
    pileup_normal <- utils::read.table(fn_pileup_normal, header = TRUE, stringsAsFactors = FALSE)
    somatic = c( somatic, list( "segments_raw" = read_CNA_segments_copycat( fn_segs ) ) )
    somatic = c( somatic, list( "rdbins" = read_CNA_rdbins_copycat( fn_bins ) ) )
    somatic = adjust_CNA_segments_filter(genome, somatic, min_segment_percent = min_segment_percent )
    somatic
}


#' Wrapper function to add gene level copy number calls from PURPLE
#' @param somatic somatic data object to modify
#' @param fn_cna path to gene level copy number file, required
#' @export
somatic_add_CNA_by_gene_purple = function( somatic, fn_cna){
    file_check("somatic_add_CNA_by_gene_purple", fn_cna)
    CNA_genes = utils::read.table( fn_cna, header=TRUE, sep='\t', stringsAsFactors=FALSE )
    names( CNA_genes )[1] = "chrom"
    dimnames(CNA_genes)[[1]] = CNA_genes$gene
    CNA_genes = CNA_genes[, c("transcriptId", "minCopyNumber","maxCopyNumber", "minMinorAlleleCopyNumber")]
    somatic = c( somatic, list( "CNA_genes" = CNA_genes ) )
    somatic
}

#' calculate gene level copy number calls using segments defined by PURPLE and gene locations
#' @param somatic somatic data object to modify
#' @param genome genome object created by Genome(), contains gene locations
#' @export
somatic_add_CNA_by_gene_from_segments = function( somatic, genome ){

    CN_weighted_mean = rep(NA, N )
    CN_largest_segment = rep(NA, N )
    gr_segments <- GenomicRanges::makeGRangesFromDataFrame(somatic$segments[,c("chrom","start","end","copyNumber")],
                                                           keep.extra.columns = T)
    # identify segments that intersect the start and end of each gene
    # trim intersecting segments to gene bounds
    # calculate 1) weighted average segment CN and 2) CN of largest segment
    for(i in 1:N ){
        gr_gene = GenomicRanges::GRanges(seqnames=genome$gene_locs$chrom[i],
                                         ranges = IRanges(start=genome$gene_locs$start[i],
                                                          end=genome$gene_locs$end[i]))
        df_segments_overlapping_gene = data.frame( subsetByOverlaps( gr_segments, gr_gene) )
        df_segments_overlapping_gene$start[ df_segments_overlapping_gene$start < genome$gene_locs$start[i] ] = genome$gene_locs$start[i]
        df_segments_overlapping_gene$end[ df_segments_overlapping_gene$end > genome$gene_locs$end[i] ] = genome$gene_locs$end[i]
        CN_weighted_mean[i] = weighted.mean( df_segments_overlapping_gene$copyNumber, df_segments_overlapping_gene$width)
        CN_largest_segment[i] = df_segments_overlapping_gene$copyNumber[ which( df_segments_overlapping_gene$width == max(df_segments_overlapping_gene$width)[1] ) ]
    }
    somatic <- c( somatic, list( CN_weighted_mean = CN_weighted_mean, CN_largest_segment = CN_largest_segment ) )
    somatic
}

#' Wrapper function to add CNA segments from PURPLE
#' @param somatic somatic data object to modify
#' @param fn_segs path to segment file, required
#' @param fn_bins path to binned copy numbers, can be null
#' @param fn_purity path to purity summary, can be null
#' @export
somatic_add_CNA_purple = function( somatic,
                                   fn_segs,
                                   fn_purity=NULL,
                                   fn_bins=NULL){

    file_check("somatic_add_CNA_purple", fn_segs)

    sample_base = somatic$sample_base
    sample_id_tumor = somatic$sample_id

    segments_raw = utils::read.table( fn_segs, header=TRUE, sep='\t', stringsAsFactors=FALSE )
    names( segments_raw )[1] = "chrom"
    segments_raw$copies = round( segments_raw$copyNumber, 3 )
    segments_raw$chrom_int = get.split.col( segments_raw$chrom, "chr", last=TRUE )
    segments_raw$chrom_int[ segments_raw$chrom_int=="X" ] = 23
    segments_raw$chrom_int[ segments_raw$chrom_int=="Y" ] = 24
    segments_raw$chrom_int = as.numeric( segments_raw$chrom_int )
    somatic = c( somatic, list( "segments_raw" = segments_raw ) )

    somatic$segments = somatic$segments_raw # already corrected for purity/ploidy

    if(!is.null( fn_bins )){
        file_check("somatic_add_CNA_purple", fn_bins)
        rdbins = utils::read.table( fn_bins, header=TRUE, sep='\t', stringsAsFactors=FALSE )
        names( rdbins )[ 1 ] = c("chrom")
        names( rdbins )[ 2 ] = c("start")
        window_size = rdbins$start[2] - rdbins$start[1] - 1
        rdbins$end = rdbins$start + window_size
        rdbins$copies = rdbins$tumorReadCount
        somatic = c( somatic, list( "rdbins" = rdbins ) )
    }else{
        somatic = c( somatic, list( "rdbins" = NULL ) )
    }

    somatic$purity = NA
    somatic$diploidProportion = NA
    somatic$ploidy = NA
    somatic$wholeGenomeDuplication = NA
    somatic$msStatus = NA
    somatic$tmbPerMb = NA
    somatic$svTumorMutationalBurden = NA
    if( !is.null( fn_purity ) ){
        file_check("somatic_add_CNA_purple", fn_purity)
        puritysummary = utils::read.table( fn_purity, header=TRUE, sep='\t', stringsAsFactors=FALSE )
        somatic$purity = puritysummary$purity
        somatic$diploidProportion = puritysummary$diploidProportion
        somatic$ploidy = puritysummary$ploidy
        somatic$wholeGenomeDuplication = puritysummary$wholeGenomeDuplication
        somatic$msStatus = puritysummary$msStatus
        somatic$tmbPerMb = puritysummary$tmbPerMb
        somatic$svTumorMutationalBurden = puritysummary$svTumorMutationalBurden
    }
    somatic
}

