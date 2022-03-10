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

#' Use CLONET to calculate ploidy
#' @param beta_table calculated by CLONET
#' @param n_cores threads, default 1
#' @param library_type WGS or WES
#' @param max_homo_dels_fraction maximum allowed fraction of segment with homozygous deletion, default 0.01
#' @param min_required_snps minimum required SNPs within a segment, default 5
#' @return the ploidy table calculated by CLONET
compute_ploidy_with_retries = function( beta_table,
                                        n_cores=1,
                                        library_type,
                                        max_homo_dels_fraction=0.01,
                                        min_required_snps = 5){
    if( utils::packageVersion("CLONETv2") == '2.2.0' ){
        pl_table =
            tryCatch( CLONETv2::compute_ploidy(beta_table,
                                     n_cores = n_cores,
                                     library_type=library_type,
                                     max_homo_dels_fraction=max_homo_dels_fraction),
                      error=function(cond){ return( NA ) },
                      warning=function(cond){ return( NA ) }
            )
    }else{
        pl_table = CLONETv2::compute_ploidy(beta_table, n_cores = n_cores)
    }
    if( is.na( pl_table[1] ) ){
        # compute_ploidy can fail to converge and throw an error in 2.2.0
        pl_table = CLONETv2::compute_ploidy(beta_table, n_cores = n_cores)
        #        pl_table = compute_ploidy(beta_table, n_cores = n_cores, library_type = library_type,
        #                                  max_homo_dels_fraction = max_homo_dels_fraction,
        #                                    beta_limit_for_neutral_reads = 0.9, min_coverage = 20,
        #                                    min_required_snps = min_required_snps,
        #                                    n_digits = 3 )
    }
    pl_table
}

#' correct raw segment CNA values using calculated purity & ploidy
#'
#' CLONET's estimate can fail, so this code has some built-in attempts to retry
#' the calculation with parameters chosen to increase the chances of a successful
#' convergence. This is not an ideal solution.
#'
#' If purity and ploidy cannot be estimated but you still want to correct the segment
#' copy number estimates, pass in purity_override and/or ploidy_override parameters
#'
#' @param genome object created by call to Genome()
#' @param somatic somatic data object to modify
#' @param pileup_normal pileup file created by ASEQ program
#' @param pileup_tumor pileup file created by ASEQ program
#' @param segment_source source, default copycat
#' @param sex sample sex, M or F
#' @param experiment_type one of WGS, WES
#' @param min_required_SNPs minimum number of SNPs in a segment to make a call
#' @param n_cores allow for multithreaded calculations
#' @param purity_override use this value instead of estimating purity
#' @param ploidy_override use this value instead of estimating ploidy
#' @export
calculate_CNA_segments_clonality_CLONET = function( genome, somatic,
                                                    pileup_normal,
                                                    pileup_tumor, sex="M",
                                                    segment_source="copycat",
                                                    experiment_type="WGS",
                                                    min_required_SNPs=10,
                                                    n_cores = 1,
                                                    purity_override = purity_override,
                                                    ploidy_override = ploidy_override){

    segmentation = somatic$segments
    segmentation <- cbind( somatic$sample_id, segmentation )
    names(segmentation)[1] = "sample_id"
    segmentation$sample_id <- as.character(segmentation$sample_id)

    if( segment_source=="copycat" ){
        segmentation = segmentation[c("sample_id", "chrom", "start", "end", "segments", "copies")]
        dimnames(segmentation)[[2]] <- c("sample_id","chr","start","end","call","logR")
        # convert estimate from chromosome count to log2
        # CopyCat doesn't know X/Y are haploid on males, so single copy X called as 2 copies by CopyCat
        if( sex=="M"){
            is_dip = segmentation$chr != "chrX" & segmentation$chr != "chrY"
            # copies is not used for calculation, but adjusting to keep consistent.
            somatic$segments$copies[ !is_dip ] = somatic$segments$copies[ !is_dip] / 2
        }else{
            is_dip = rep(TRUE, dim(segmentation)[1])
        }
        segmentation$logR <- log2(segmentation$logR / 2)
    }else{
        if( sex=="M"){
            is_dip = segmentation$chr != "chrX" & segmentation$chr != "chrY"
        }else{
            is_dip = rep(TRUE, dim(segmentation)[1])
        }
        names(segmentation)[ which(names(segmentation)=="chrom") ] = "chr"
        names(segmentation)[ which(names(segmentation)=="log2") ] = "logR"
        segmentation = segmentation[,c("sample_id", "chr", "start", "end", "depth", "logR")]
    }
    beta_table <- CLONETv2::compute_beta_table(seg_tb = segmentation,
                                     pileup_tumor = pileup_tumor,
                                     pileup_normal = pileup_normal,
                                     min_required_snps=min_required_SNPs,
                                     n_cores = n_cores)

    # Compute ploidy and purity
    #
    pl_table = compute_ploidy_with_retries( beta_table, n_cores=n_cores, library_type=experiment_type,
                                            min_required_snps=min_required_SNPs, max_homo_dels_fraction = 0.01)
    if( utils::packageVersion("CLONETv2") == '2.2.0' ){
        adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                       ploidy_table = pl_table,
                                       library_type = experiment_type,
                                       min_required_snps=min_required_SNPs,
                                       n_cores = n_cores)
    }else{
        adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                           ploidy_table = pl_table,
                                           min_required_snps=min_required_SNPs,
                                           n_cores = n_cores)
    }
    if( !is.na(purity_override ) ){
        purity = purity_override
        adm_table[1,"adm"] = 1-purity
    }
    purity = 1-adm_table[1,"adm"]
    #cpt = check_ploidy_and_admixture(beta_table, pl_table, adm_table)
    #ggsave("/data1/projects/WCDT_WGS_paired/results_v2/DTB-137-PRO/purity/cpa_plot_max_01.pdf")
    if( is.na(purity) ){
        pl_table = compute_ploidy_with_retries( beta_table, n_cores=n_cores, library_type=experiment_type,
                                                min_required_snps = min_required_SNPs, max_homo_dels_fraction = 0.15)
        if( utils::packageVersion("CLONETv2") == '2.2.0' ){
            adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                               ploidy_table = pl_table,
                                               library_type = experiment_type,
                                               min_required_snps=min_required_SNPs,
                                               n_cores = n_cores)
        }else{
            adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                               ploidy_table = pl_table,
                                               min_required_snps=min_required_SNPs,
                                               n_cores = n_cores)
        }
        purity = 1-adm_table[1,"adm"]
    }
    if( is.na(purity) ){
        pl_table = compute_ploidy_with_retries( beta_table, n_cores=n_cores, library_type=experiment_type,
                                                min_required_snps = min_required_SNPs, max_homo_dels_fraction = 0.2)
        if( utils::packageVersion("CLONETv2") == '2.2.0' ){
            adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                               ploidy_table = pl_table,
                                               library_type = experiment_type,
                                               min_required_snps=min_required_SNPs,
                                               n_cores = n_cores)
        }else{
            adm_table <- CLONETv2::compute_dna_admixture(beta_table = beta_table,
                                               ploidy_table = pl_table,
                                               min_required_snps=min_required_SNPs,
                                               n_cores = n_cores)
        }
        purity = 1-adm_table[1,"adm"]
    }
    ploidy = pl_table[1,2]
    if( !is.na( ploidy_override )){
        pl_table[1, 2] = ploidy_override

    }
    ast = CLONETv2::compute_allele_specific_scna_table(beta_table,
                                             pl_table,
                                             adm_table,
                                             error_tb = CLONETv2::error_table,
                                             n_cores = n_cores)
    somatic$segments$copies_int_corr = ast$cnA.int + ast$cnB.int
    somatic$segments$logR_corr = round( ast$log2.corr, 2)
    somatic$segments$copies_major_corr = round( ast$cnA, 2)
    somatic$segments$copies_minor_corr = round( ast$cnB, 2)
    somatic$segments$copies_corr[is_dip] = round( 2* 2^( ast$log2.corr[is_dip] ), 2)
    somatic$segments$copies_corr[!is_dip] = round( 2^( ast$log2.corr[!is_dip] ), 2)
    corr_failed = is.na( somatic$segments$logR_corr )
    somatic$segments$copies_corr[corr_failed] = somatic$segments$copies[corr_failed]

    somatic = c(somatic, list( "purity" = purity ) )
    somatic = c(somatic, list( "ploidy" = ploidy ) )
    somatic = c(somatic, list( "beta_table" = beta_table) )
    somatic = c(somatic, list( "pl_table" = pl_table) )
    somatic = c(somatic, list( "adm_table" = adm_table) )
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
    somatic = calculate_CNA_segments_clonality_CLONET( genome, somatic,
                                                       pileup_normal, pileup_tumor, sex=sex,
                                                       segment_source=CNA_method,
                                                       experiment_type=experiment_type,
                                                       n_cores = n_cores,
                                                       min_required_SNPs = min_required_SNPs,
                                                       purity_override = purity_override,
                                                       ploidy_override = ploidy_override)

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
    somatic = calculate_CNA_segments_clonality_CLONET( genome, somatic,
                                                       pileup_normal, pileup_tumor, sex=sex,
                                                       segment_source = CNA_method,
                                                       experiment_type = "WGS",
                                                       n_cores = n_cores,
                                                       min_required_SNPs = min_required_SNPs,
                                                       purity_override = purity_override,
                                                       ploidy_override = ploidy_override)
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

