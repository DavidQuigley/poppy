#' wrapper to check if file exists
#' @param func_name report back which function was called
#' @param fn file path to check
file_check=function(func_name, fn){
    if( !file.exists( fn ) ){
        stop( paste("Missing in ",func_name,": ", fn, sep="" ) )
    }
}


#' Extract rdbins within a genomic range
#' @param genome Genome object
#' @param somatic somatic object from which to extract
#' @param chrom chromosome
#' @param loc_start location to start
#' @param loc_end location to end
#' @return subset of rdbins
copynumber_range = function( genome, somatic, chrom, loc_start, loc_end ){

    idx = which( somatic$rdbins$chrom==chrom &
                     somatic$rdbins$start>= loc_start &
                     somatic$rdbins$end<= loc_end )
    somatic$rdbins[idx,]
}


#' Identify microhomology between two loci in the genome on the same chromosome
#' @param chrom chromosome
#' @param pos_start starting locus
#' @param pos_end ending locus
#' @param n_shoulder number of nucleotides in either direction to scan
identify_microhomology = function( chrom, pos_start, pos_end, n_shoulder ){

    #    getSeq(Hsapiens,"chr3",start=117018070,end=117018089,strand="+")
    #    ATTTTCCTTTTCCTTTTTT
    #             TTCCTT
    #                           ATGTTCCTTCTTCCCAGAAT
    #                                     TTCCCAG

    # chrom="chr13"
    # pos_start = 32000010
    # pos_end = 32000061
    # n_shoulder=6
    # output: 4
    #          |                                                 |
    #AGCGTCTCTGCCCGGCCGCCCATCATCGTCTGAGATGTGGGGAGCGCCTCTGCCCCGCCGCCCCGTCTGGG
    #    TCTCTG                                                   CCCG
    #      sh5                                                     sh3
    #           CCCG                                        CGCCG
    #           inner5                                       inner3

    n_microhomology = 0
    if( pos_start > pos_end ){
        tmp = pos_start
        pos_start = pos_end
        pos_end = tmp
    }
    if( pos_start - 1 - n_shoulder < 0 ){
        n_microhomology = 0
    }else{
        sh5 = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                     chrom,
                     start = pos_start - n_shoulder,
                     end   = pos_start - 1,
                     strand="+" )
        sh3 = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                     chrom,
                     start = pos_end + 1,
                     end   = pos_end + 1 + n_shoulder,
                     strand="+" )
        inner5 = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                        chrom,
                        start = pos_start + 1 ,
                        end =   pos_start + 1 + n_shoulder,
                        strand="+" )
        inner3 = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                        chrom,
                        start = pos_end - n_shoulder,
                        end =   pos_end - 1,
                        strand="+" )
        n_hom_inner5 = 0
        n_hom_inner3 = 0
        for( idx5 in 1 : n_shoulder ){
            if( inner5[idx5] == sh3[idx5] ){
                n_hom_inner5 = n_hom_inner5 + 1
            }else{
                break
            }
        }
        for( idx3 in 1:n_shoulder){
            inv_idx = idx3 * -1
            if( inner3[inv_idx] == sh5[inv_idx] ){
                n_hom_inner3 = n_hom_inner3 + 1
            }else{
                break
            }
        }
        n_microhomology = max( c(n_hom_inner5, n_hom_inner3) )
    }
    n_microhomology
}

#' return dataframe of indices into A and B restricted to perfect matches
#'
#' If there is more than one match between A and B, reports only the first match
#'
#' @param A first vector
#' @param B second vector
#' @param allow.multiple.B if TRUE, allow single hits in A to match multiple hits in B
match.idx = function(A, B, allow.multiple.B=FALSE){
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}

#' splits elemts of a vector
#' @param v vector
#' @param string search string
#' @param col index of the token to return
#' @param last return the last token
#' @param first return the first token
get.split.col = function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")

    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}

#' Convert chromosomes and positions to a scaled proportionate value
#' @param chr.list chromosomes
#' @param loc.list locations on those chromosomes
#' @param x.max highest value for scaling
#' @param chrom_lengths length of each chromosome, defaults to hardcoded hg38 lengths
convert.genome.to.scaled.x=function( chr.list, loc.list, x.max=1000, chrom_lengths=NULL ){

    chrom_lengths_hg38 = c(248956422,242193529,198295559,190214555,181538259,
                           170805979,159345973,145138636,138394717,133797422,
                           135086622,133275309,114364328,107043718,101991189,
                           90338345,83257441,80373285,58617616,64444167,46709983,
                           50818468,156040895,57227415)
    if( is.null(chrom_lengths) ){
        chrom_lengths = chrom_lengths_hg38
    }
    if( !is.numeric( chr.list )){
        chr.list = get.split.col( chr.list, "chr", last=TRUE)
        chr.list[ chr.list=="X" ] = 23
        chr.list[ chr.list=="Y" ] = 24
        chr.list = as.numeric(chr.list)
    }
    size.chroms = data.frame( chr=1:24,
                              loc=chrom_lengths,
                              stringsAsFactors = FALSE)
    size.genome = sum(size.chroms)
    xs = rep(0.0, length(chr.list))
    for(i in 1:length(xs)){
        cur.chrom = chr.list[i]
        cur.loc = loc.list[i]
        if( cur.chrom>1 ){
            xs[i] = sum(size.chroms$loc[size.chroms$chr < cur.chrom]) + as.numeric(cur.loc)
        }
        else{
            xs[i] = cur.loc
        }
    }
    round( xs / size.genome * x.max )
}


#' used to taste files to see what kind they are
#' @param fn path to file to test
#' @return first token in file, converted to numeric
read_first_token_from_file = function(fn){
    as.numeric( utils::read.table(fn, stringsAsFactors = FALSE, header=FALSE)$V1[1] )
}
