#' Wrapper function to add structural variants and fusions annotations from Gridss and Linx, respectivelly.
#' @param somatic somatic data object to modify
#' @param fn_linx path to fusion annotations linx text file (e.g. SAMPLE_ID.linx.fusion.tsv)
#' @param fn_gripss path to structural variants filtered gridss vcf file (e.g. SAMPLE_ID_GRIPSS.somatic.filtered.vcf)
#' @export
#' library( VariantAnnotation )
#' library( stringr )
#' library( plyr )
#' library( poppy )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### new function add_SV_gridss
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_SV_gridss <-function (somatic, fn_gripss, MIN_READS = 2) 
{
    sample_id_tumor = somatic$sample_id
    file_check("add_SV_gridss", fn_gripss)
    vcf = readVcf(fn_gripss, "hg38")
    if ( length(rowRanges(vcf)) != 0){
        vcf_f = vcf[VariantAnnotation::fixed(vcf)$FILTER == "PASS" || VariantAnnotation::fixed(vcf)$FILTER == "PON"]
        reads = parse_gridss_vcf(vcf_f)
        keep = reads$rp_n_alt == 0 & reads$vf_n_alt == 0 & reads$refpair_t + 
            reads$vf_t_alt >= MIN_READS
        vcf_f = vcf_f[keep]
        
        somatic$gn_bnd = length(unique(info(vcf_f)[grepl("BND",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gn_del = length(unique(info(vcf_f)[grepl("DEL",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gn_dup = length(unique(info(vcf_f)[grepl("DUP",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gn_ins = length(unique(info(vcf_f)[grepl("INS",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gn_inv = length(unique(info(vcf_f)[grepl("INV",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gn_sgl = length(unique(info(vcf_f)[grepl("SGL",info(vcf_f)$EVENTTYPE),]$EVENT))
        somatic$gridss = vcf_f
        somatic
    }
    else {
        print (paste("None somatic SV detected or passed the filters by GRIPSS on", somatic$sample_id, "sample" ) )
        if(length(rowRanges(vcf)) == 0){
            somatic$gridss = vcf
        }
        else{
            somatic$gridss = GRanges(seqnames = NULL, ranges = NULL)
        }
        somatic$gn_bnd = 0
        somatic$gn_del = 0
        somatic$gn_dup = 0
        somatic$gn_ins = 0
        somatic$gn_inv = 0
        somatic$gn_sgl = 0
        somatic
    }
    
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new function parse_gridss_vcf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parse_gridss_vcf <- function( vcf_f ){
    reads = data.frame(
        rp_n_alt=rep(0, length(vcf_f)),
        rp_t_alt=rep(0, length(vcf_f)),
        pr_ref_n=rep(0, length(vcf_f)),
        vf_n_alt=rep(0, length(vcf_f)),
        pr_ref_t=rep(0, length(vcf_f)),
        vf_t_alt=rep(0, length(vcf_f)),
        refpair_n=rep(0, length(vcf_f)),
        refpair_t=rep(0, length(vcf_f)),
        TAF=rep(0, length(vcf_f)),
        hom_len=rep(NA, length(vcf_f)),
        stringsAsFactors = FALSE
    )
    rp = geno(vcf_f)$RP #equivalent to pr_n_alt and pr_t_alt PR (spanning pair-reads) from manta
    pr_ref = geno(vcf_f)$REF #similar to sr_n_wt and sr_t_wt
    vf = geno(vcf_f)$VF #equivalent to sr_n_alt and sr_t_alt
    refpair = geno(vcf_f)$REFPAIR # read pairs "spanning" this breakend supporting the reference allele and necessary to calculated VAF
    TAF = info(vcf_f)$TAF
    hom_len = info(vcf_f)$HOMLEN
    
    for(i in 1:length(vcf_f)){    
        
        reads$rp_n_alt[i] = rp[i,1]
        reads$rp_t_alt[i] = rp[i,2]
        reads$pr_ref_n[i] = pr_ref[i,1]
        reads$vf_n_alt[i] = vf[i,1]
        reads$pr_ref_t[i] = pr_ref[i,2]
        reads$vf_t_alt[i] = vf[i,2]
        reads$refpair_n[i] = refpair[i,1]
        reads$refpair_t[i] = refpair[i,2]
        reads$TAF[i] = TAF[i]
        if (is.integer(hom_len[[i]]) && length(hom_len[[i]]) != 0L)
        {
            reads$hom_len[i] = hom_len[[i]]   
        }
        else{
            reads$hom_len[i] = NA
        }
        
    }
    
    reads$vf_n_alt[ is.na(reads$vf_n_alt) ] = 0
    reads$pr_ref_t[ is.na(reads$pr_ref_t) ] = 0
    reads$refpair_n[ is.na(reads$refpair_n) ] = 0
    reads$refpair_t[ is.na(reads$refpair_t) ] = 0
    reads
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new function calculate_gridss_summary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' GRIDSS no longer writes the SVLEN header [https://github.com/PapenfussLab/gridss/issues/561]
#' SVLEN is calculated manually


calculate_gridss_summary<- function (somatic, nt_for_microhomology = 15) 
{
    
    if ( dim(somatic$gridss)[1] != 0){
        df = data.frame(svtype = as.character(info(somatic$gridss)$EVENTTYPE), 
                        chrom_1 = as.character(seqnames(somatic$gridss)), start_1 = start(somatic$gridss), 
                        end_1 = end(somatic$gridss), chrom_2 = as.character(seqnames(somatic$gridss)), 
                        start_2 = start(somatic$gridss), end_2 = end(somatic$gridss), 
                        filter = as.character(rowRanges(somatic$gridss)$FILTER), 
                        row.names = names(rowRanges(somatic$gridss)), stringsAsFactors = FALSE)
        if (dim(df)[1] > 0) {
            df$sample_id = rep(somatic$sample_id, dim(df)[1])
            df$sv_length = rep(NA, dim(df)[1])
            mate_id = rep(NA, dim(df)[1])
            
            str2 = unlist(alt(somatic$gridss))
            is_left = str_detect(string = str2, pattern = stringr::fixed("["))
            is_right = str_detect(string = str2, pattern = stringr::fixed("]"))
            chrom_dest = rep("", length(str2))
            pos_dest = rep(0, length(str2))
            for (x in 1:length(str2)) {
                if (is_left[x]) {
                    slug = strsplit(str2[x], "[", fixed = TRUE)[[1]][2]
                }
                else {
                    slug = strsplit(str2[x], "]", fixed = TRUE)[[1]][2]
                }
                chr_pos = strsplit(slug, ":")[[1]]
                chrom_dest[x] = chr_pos[1]
                pos_dest[x] = chr_pos[2]
            }
            df$chrom_2 = chrom_dest
            df$start_2 = as.numeric(pos_dest)
            df$end_2 = df$start_2 + 1
            df$end_1 = df$start_1 + 1
        }
        
        for (i in 1:length(df$sv_length)) {
            
            if (df$svtype[i] != "BND") {
                df$sv_length[i] = df$end_1[i] - df$start_2[i]
            }
            mid = info(somatic$gridss)$MATEID[[i]]
            if (length(mid) > 0) {
                mate_id[i] = mid
            }
        }
        df$mate_id = mate_id
        df = cbind(df, parse_gridss_vcf(somatic$gridss))
        ##df$n_microhomology = rep(NA, dim(df)[1])
        ##for (i in 1:dim(df)[1]) {
        #  if (df$svtype[i] == "DEL") {
        #     df$n_microhomology[i] = identify_microhomology(df$chrom_1[i], 
        #                                                   df$start_1[i], df$end_1[i], nt_for_microhomology)
        #}
        #}
        somatic$list_sv_gripss = df
        somatic$nt_for_microhomology = nt_for_microhomology
        somatic
    }
    else {
        print (paste("None somatic SV detected or passed the filters by GRIPSS on", somatic$sample_id, "sample" ) )
        somatic$list_sv_gripss = ""
        somatic$nt_for_microhomology = nt_for_microhomology
        somatic
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### new function add_fusions_linx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

somatic_add_fusion_linx <- function( somatic, fn_linx ){
    file_check("somatic_add_fusion_linx", fn_linx)
    fusions_linx <- read.table(fn_linx, header=TRUE, sep = "\t",
                               stringsAsFactors=FALSE)
    if (dim(fusions_linx)[1]==0){
        print(paste("None fusions detected or passed the filters by LINX on", somatic$sample_id, "sample"))
    }
    else{
        print(paste(dim(fusions_linx)[1]," fusions detected by LINX on", somatic$sample_id, "sample"))
    }
    somatic <- c( somatic, list( fusions_linx = fusions_linx ) )
    somatic
}
