
#' make_dt
#'
#' Creates a data.table from a vector of files.
#'
#' @param files character vector of file paths.
#' @param group_lev If provided, should be directory names at equivalent depth
#'   in file paths. Default of NULL disables this functionality and group will
#'   be none and batch will be A.
#' @param max_name_len Names are truncated to this length. Default is 30.
#'
#' @return A data.table with attribute "file" and "name", "group", and "batch" attributes.
#' @export
#'
#' @examples
#' files = c("exp1/file1", "exp1/file2", "exp2/file3")
#' #default no group levels
#' make_dt(files)
#' dt1 = make_dt(files, group_lev = c("exp1", "exp2"))
#' levels(dt1$group)
#' dt1$batch
#' #reversed group_lev
#' dt2 = make_dt(files, group_lev = rev(c("exp1", "exp2")))
#' levels(dt2$group)
#' dt2$batch
make_dt = function(files, group_lev = NULL, max_name_len = 30){
  file_dir = group = batch = name = NULL#global data.table bindings
  p_dt = data.table(file = files)
  p_dt[, sample := sub("\\..+", "", basename(file))]
  p_dt[, sample := sub("(?<=rep[0-9])_.+", "", sample, perl = TRUE)]

  p_dt[, file_dir := file]

  if(is.null(group_lev)){
    p_dt$group = "none"
    p_dt$batch = factor("A")
  }else{
    while(!any(p_dt$file_dir == "/")){
      if(all(basename(p_dt$file_dir) %in% group_lev)){
        p_dt[, group := basename(file_dir)]
        break
      }else{
        p_dt[, file_dir := dirname(file_dir)]
      }
    }
    p_dt$group = factor(p_dt$group, levels = group_lev)
    p_dt[, batch := LETTERS[as.numeric(group)]]
  }
  p_dt[, name := paste0(substr(sample, 0, max_name_len), ".", batch)]
  dupe_num = 1
  if(any(duplicated(p_dt$name))){
    p_dt[, dupe_num := seq(.N)-1, list(name)]
    p_dt[dupe_num > 0, name := paste0(name, ".", dupe_num) ]
    p_dt$dupe_num = NULL
  }
  p_dt[]
}

#' make_anno_grs
#'
#' Loads a GTF file and processes it to list of GRanges suitable for serial
#' annotation by \code{\link{make_feature_as_signal_dt}}.
#'
#' Output is intended for use with \code{\link{make_feature_as_signal_dt}}
#'
#' @param gtf_file A gtf or gtf.gz file from GENCODE  Other sources may work but
#'   have not been tested.
#'
#' @return a named list of GRanges where each corresponds to an annotated
#'   feature type
#' @export
#'
#' @import rtracklayer GenomicRanges
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' make_anno_grs(gtf_file)
make_anno_grs = function(gtf_file){
  type = tag = NULL # global binding for data.table
  ref_gr = rtracklayer::import.gff(gtf_file, format = "gtf")
  gene_gr = reduce(subset(ref_gr, type == "gene"))
  exon_gr = reduce(subset(ref_gr, type == "exon" & tag == "basic"))
  intron_gr = setdiff(gene_gr, exon_gr)
  tx_gr = subset(ref_gr, type == "transcript" & tag == "basic")
  tss_gr = flank(tx_gr, 1e3, start = TRUE, both = TRUE)
  tts_gr = flank(tx_gr, 1e3, start = FALSE, both = TRUE)
  artifact_gr = rtracklayer::import.bed("~/R/qc_cutnrun/reference/blacklist_hg38.bed")
  anno_grs = rev(list(artifact = artifact_gr, tss = tss_gr, tts = tts_gr, exon = exon_gr, intron = intron_gr, genebody = gene_gr))
  anno_grs
}

#' make_anno_dt
#'
#' @param peak_grs a named list of GRanges where each feature set should be
#'   annotated. Use output from \link[seqsetvis]{easyLoad_narrowPeak} or
#'   similar.
#' @param anno_grs a named list of GRanges where each corresponds to an
#'   annotated feature type. Use output from \code{\link{make_anno_grs}}
#' @param name_lev Optional levels to impose order in feature sets.  Default of
#'   NULL uses input order of peak_grs.
#'
#' @return data.table with annotation overlaps for inerval sets in peak_grs.
#' @export
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#'
#' make_anno_dt(peak_grs, anno_grs)
make_anno_dt =  function(peak_grs, anno_grs, name_lev = NULL){
  sample_cnt = count = fraction = type = tag = NULL #data.table global declaration
  if(is.null(names(peak_grs))){
    names(peak_grs) = paste("peaks", LETTERS[seq_along(peak_grs)])
  }
  if(is.null(name_lev)){
    name_lev = names(peak_grs)
  }

  message("features overlap")
  .apply_annotation = function(x){
    x$anno = "intergenic"
    for(i in seq_along(anno_grs)){
      olaps = findOverlaps(x, anno_grs[[i]])
      x$anno[S4Vectors::queryHits(olaps)] = names(anno_grs)[i]

    }
    x
  }

  peak_grs.anno = lapply(peak_grs, .apply_annotation)

  .dt_count_anno = function(x){
    tab = table(x$anno)
    data.table(feature = names(tab), count = as.numeric(tab))
  }

  peak_grs.anno_cnt = lapply(peak_grs.anno, .dt_count_anno)

  anno_cnt = rbindlist(peak_grs.anno_cnt, idcol = "sample")
  anno_cnt[, sample_cnt := paste0(sample, "\n", sum(count)), list(sample)]
  anno_cnt[, fraction := count / sum(count), list(sample)]
  stopifnot(setequal(anno_cnt$sample, name_lev))
  anno_cnt$sample = factor(anno_cnt$sample, levels = name_lev)
  anno_cnt = anno_cnt[order(sample)]
  anno_cnt$sample_cnt = factor(anno_cnt$sample_cnt, levels = unique(anno_cnt$sample_cnt))
  anno_cnt
}

#' make_assign_dt
#'
#' @param clust_dt Output from \link[seqsetvis]{ssvSignalClustering}
#' @param cluster_var variable name with cluster assignment info. Default is "cluster_id".
#' @param id_var variable name with id info. Default is "id".
#'
#' @return data.table with id and cluster_id assignment info.
#' @export
#'
#' @examples
#' bw_files = dir(system.file("extdata", package = "seqqc"), pattern = "^M.+bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 6e2, fix = "center")
#'
#' prof_dt = seqsetvis::ssvFetchBigwig(query_dt, query_gr, return_data.table = TRUE)
#'
#' clust_dt = seqsetvis::ssvSignalClustering(prof_dt, nclust = 3)
#'
#' assign_dt = make_assign_dt(clust_dt)
#' assign_dt
#'
make_assign_dt = function(clust_dt, cluster_var = "cluster_id", id_var = "id"){
  assign_dt = unique(clust_dt[, list(get(cluster_var), get(id_var))])
  setnames(assign_dt, c(cluster_var, id_var))
  assign_dt
}

#' make_feature_as_signal_dt
#'
#' Create a data.table of overlaps between query_gr and anno_grs (from
#' \code{\link{make_anno_grs}} )
#'
#' @param anno_grs named list of GRanges objects where each is an annotation
#'   class.  Use \code{\link{make_anno_grs}} to generate.
#' @param query_gr GRanges of regions to characterize via anno_grs.
#'
#' @return a data.table of counts for each annotation class.
#' @export
#'
#' @import seqsetvis
#'
#' @examples
#' gtf_file = system.file(package = "seqqc", "extdata/gencode.v35.annotation.at_peaks.gtf")
#' anno_grs = make_anno_grs(gtf_file)
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"), pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#' query_gr = resize(seqsetvis::ssvOverlapIntervalSets(peak_grs), 2e4, fix = "center")
#' anno_dt = make_feature_as_signal_dt(anno_grs, query_gr)
#' anno_dt$id = factor(anno_dt$id, levels = rev(unique(anno_dt$id)))
#' ggplot(anno_dt, aes(x = x, fill = y, y = id)) + geom_tile() + facet_wrap(~sample)
make_feature_as_signal_dt = function(anno_grs, query_gr){
  seqsetvis::ssvFetchGRanges(anno_grs, query_gr, return_data.table = TRUE)
}


#' make_fq_dt
#'
#' @param fastq_files paths to fast files (can be gzipped with .gz extension)
#' @param fastq_names optional parallel vector of names for fastq files.
#'   Defaults to basename of fastq_files. Should be unique.
#' @param fastq_treatments optional parallel vector of treatments. Defaults to
#'   fastq_names. May be duplicated.
#' @param n_cores number of cores to use to count lines in fastq files. Defaults
#'   to mc.cores if set or 1.
#' @param cache_counts logical. Should the counts be saved to *.cnt files
#'   alongside the fastq_files?  Default is TRUE
#'
#' @return a data.table countaining fastq, count, name, and treatment attributes
#' @export
#'
#' @examples
#' fq_files = dir("inst/extdata", pattern = "(fq$)|(fq.gz$)|(fastq$)|(fastq.gz$)", full.names = TRUE)
#'  #no idea why this make_fq_dt example won't run
#'  #make_fq_dt(fq_files,
#'  #  fastq_names = c("4_reads_fq", "4_reads_gz", "5_reads_fq", "5_reads_gz"),
#'  #  cache_counts = FALSE)
make_fq_dt = function(fastq_files, fastq_names = basename(fastq_files), fastq_treatments = c("fq", "gz", "fq", "gz"), n_cores = getOption("mc.cores", 1), cache_counts = TRUE){
  message("count fastq reads...")

  .cnt_fq = function(f){
    cnt_f = paste0(f, ".cnt")
    if(!file.exists(cnt_f)){
      if(grepl(".gz$", f)){
        cnt = system(paste0("gunzip -c ", f, "| wc -l"), intern = TRUE)
      }else{
        cnt = system(paste0("cat ", f, "| wc -l"), intern = TRUE)
      }

      cnt_dt = data.table(f, as.numeric(cnt)/4)
      if(cache_counts){
        fwrite(cnt_dt, cnt_f, sep = "\t", col.names = FALSE)
      }
    }else{
      cnt_dt = fread(cnt_f, sep = "\t")
    }
    # message(class(cnt_dt))
    cnt_dt
  }

  fq_dt = data.table::rbindlist(pbmcapply::pbmclapply(fastq_files, mc.cores = n_cores, .cnt_fq))
  setnames(fq_dt, c("fastq", "count"))
  fq_dt$name = fastq_names
  fq_dt$treatment = fastq_treatments
  fq_dt
}

#' make_centered_query_gr
#'
#' Returns version of query_gr GRanges centered on the maximum signal in
#' bams/bigwigs supplied by query_dt.
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges to be centered.
#' @param view_size Size of final regions.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed on the ssvFechBam or ssvFetchBigwig. Do not use, n_cores,
#'   win_size, win_method, return_data.table or n_region_splits.
#'
#' @return GRanges centered on the maximum signal
#' @export
#'
#' @examples
#' library(seqsetvis)
#' #bigwig example with 3 bigwigs
#' bw_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "bw$", full.names = TRUE)
#' query_dt = make_dt(bw_files)
#' query_dt[, sample := sub("_FE_random100.A", "", name)]
#'
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = easyLoad_narrowPeak(peak_files)
#' query_gr = resize(ssvOverlapIntervalSets(peak_grs), 700, fix = "center")
#'
#'
#'
#' #heatmap before centering
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBigwig(query_dt, query_gr), nclust = 4) +
#'   labs(title = "Before centering")
#'
#' query_gr.centered = make_centered_query_gr(query_dt, query_gr)
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBigwig(query_dt, query_gr.centered), nclust = 4)+
#'   labs(title = "After centering")
#'
#' #bam example with 1 bam, query_dt can just be file paths
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr2 = easyLoad_bed(peak_file)[[1]]
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBam(bam_file, query_gr2), nclust = 4) +
#'   labs(title = "Before centering")
#'
#' query_gr2.centered = make_centered_query_gr(bam_file, query_gr2, fragLens = 180)
#' set.seed(0)
#' ssvSignalHeatmap(ssvFetchBam(bam_file, query_gr2.centered), nclust = 4) +
#'   labs(title = "After centering")
make_centered_query_gr = function(query_dt, query_gr, view_size = NULL,
                                  n_cores = getOption("mc.cores", 1), ...){
  if(is.character(query_dt)) query_dt = data.table(file = query_dt)
  stopifnot(is.data.table(query_dt))
  if(!is.null(query_dt$file)){
    files = query_dt$file
  }else{
    files = query_dt[[1]]
  }
  is_bam = grepl(".bam$", files)
  if(all(is_bam)){
    fetch_fun = seqsetvis::ssvFetchBam
  }else if(all(!is_bam)){
    fetch_fun = seqsetvis::ssvFetchBigwig
  }else{
    stop("Files should be either all bams or all bigwigs. No mixing!")
  }
  stopifnot("GRanges" %in% class(query_gr))
  if(is.null(view_size)){
    view_size = median(width(query_gr))
  }
  n_region_splits = max(1, floor(length(query_gr) / 1e3))
  if(n_cores == 1) n_region_splits = 1
  if(max(width(query_gr)) > 100){
    message("seeking coarse maxima...")
    ws = ceiling(max(width(query_gr)) / 100)
    if(ws > min(width(query_gr))){
      stop("Too much variation in width of input query_gr.  Please increase minimum or decrease maximum width and retry.")
    }
    prof_dt.coarse = fetch_fun(query_dt, query_gr, n_cores = n_cores,
                               win_size = ws, win_method = "summary",
                               return_data.table = TRUE,
                               n_region_splits = n_region_splits, ...)
    query_gr = unique(resize(centerGRangesAtMax(prof_dt.coarse, query_gr), 100, fix = 'center'))
  }
  message("seeking fine maxima...")
  prof_dt.fine = fetch_fun(query_dt, query_gr, n_cores = n_cores,
                           win_size = 1, win_method = "sample",
                           return_data.table = TRUE,
                           n_region_splits = n_region_splits, ...)

  centered_gr = unique(resize(centerGRangesAtMax(prof_dt.fine, query_gr), view_size, fix = 'center'))
  centered_gr
}

#' make_frip_dt
#'
#' create a data.table with FRIP data.  To be used with
#' \code{\link{plot_frip_dt}}
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges to be centered.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#'
#' @return a data.table with reads_in_peak, mapped_reads, and frip data for all
#'   bam files in query_dt at each region in query_gr.
#' @export
#' @import Rsamtools
#' @importFrom stats median
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#'
#' make_frip_dt(bam_file, query_gr)
make_frip_dt = function(query_dt, query_gr, n_cores = getOption("mc.cores", 1)){
  name = treatment = qname = id = frip = N = mapped_reads = V1 = NULL#global data.table bindings
  if(is.character(query_dt)){
    query_dt = data.table(file = query_dt)
  }
  if(is.null(query_dt$file)){
    stop("query_dt must contain file variable")
  }
  if(is.null(query_dt$name)){
    query_dt[, name := basename(file)]
  }
  if(is.null(query_dt$treatment)){
    query_dt[, treatment := name]
  }

  n_region_splits = max(1, floor(length(query_gr) / 1e3))
  message("fetch read counts...")
  reads_dt = seqsetvis::ssvFetchBam(query_dt, query_gr, fragLens = NA, return_unprocessed = TRUE, n_region_splits = n_region_splits, n_cores = n_cores)

  frip_dt = reads_dt[, list(N = length(unique(qname))), list(id, name, treatment, sample)]
  unique(frip_dt[, list(name, treatment, sample)])

  frip_dt_filled = melt(dcast(frip_dt, id~name, value.var = "N", fill = 0), id.vars = "id", value.name = "N", variable.name = "name")
  frip_dt = merge(frip_dt_filled, unique(frip_dt[, list(name, treatment, sample)]), by = "name")

  message("fetch total mapped reads...")
  mapped_counts = sapply(query_dt$file, function(f){
    stats = Rsamtools::idxstatsBam(f)
    stats = subset(stats, grepl("chr[0-9XY]+$", seqnames ))
    sum(stats[,3])
  })
  frip_dt$mapped_reads = mapped_counts[frip_dt$sample]
  frip_dt[, frip := N/mapped_reads]

  name_lev = frip_dt[, stats::median(N) , list(name)][rev(order(V1))]$name
  stopifnot(all(frip_dt$name %in% name_lev))
  frip_dt$name = factor(frip_dt$name, levels = name_lev)
  setnames(frip_dt, "N", "reads_in_peak")
  frip_dt

}

#' make_peak_dt
#'
#' @param peak_grs list of GRanges for peak sets
#' @param treatments opional character vector of treatments for each peak set.
#'
#' @return a data.table with peak_count data for each GRanges in peak_grs.
#' @export
#' @import seqsetvis
#' @examples
#' peak_files = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "Peak$", full.names = TRUE)
#' peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files)
#'
#' make_peak_dt(peak_grs)
make_peak_dt = function(peak_grs, treatments = NULL){
  treatment = NULL#global binding for data.table
  peak_dt = seqsetvis::ssvFeatureBars(peak_grs, return_data = TRUE)
  if(is.null(treatments)){
    peak_dt$treatment = peak_dt$group
  }else{
    stopifnot(nrow(peak_dt) == length(treatments))
    peak_dt$treatment = treatment
  }
  if(is.null(names(peak_grs))){
    peak_dt$name = peak_dt$treatment
  }else{
    peak_dt$name = names(peak_grs)
  }
  setnames(peak_dt, "count", "peak_count")
  peak_dt
}



#' make_scc_dt
#'
#' Calculate Strand Cross Correlation (SCC) for several bam files defined in
#' query_dt at regions defined by query_gr and return tidy data.table.
#'
#' @param query_dt data.table with query information. Only really needs file as
#'   first column.
#' @param query_gr GRanges of regions to calculate SCC for
#' @param frag_sizes optional numeric. Fragment sizes to calculate correlation
#'   at.  The higher the resolution the longer calculation will take.  The
#'   default is to count by 10 from 50 to 350.
#' @param fetch_size optional numeric. Size in bp centered around each interval
#'   in query_gr to retrieve.  Should be greater than max frag_size. The default
#'   is 3*max(frag_sizes).
#' @param cache_path path to cache location for BiocFileCache to use.
#' @param cache_version Modifying the cache version will force recalulation of
#'   all results going forward. Default is v1.
#' @param force_overwrite Logical, if TRUE, cache contents will be overwritten.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed to Rsamtools::ScanBamParam()
#'
#' @return list fo tidy data.table of SCC data for every bam file in query_dt
#' @export
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#' query_dt = data.table(file = rep(bam_file, 2))
#' #attributes added here will be carried forward to the final data.tables
#' query_dt$sample = c("A", "B")
#' query_dt$value = c(.3, .7)
#'
#' scc_res = make_scc_dt(query_dt, query_gr)
#'
#' #making the data.table can be skipped for conveinence
#' #also run at higher resolution if a more precise estimate of fragment size is required
#' scc_res = make_scc_dt(bam_file, query_gr, frag_sizes = seq(150, 210, 1))
make_scc_dt = function(query_dt,
                       query_gr,
                       frag_sizes = seq(50, 350, 10),
                       fetch_size = 3*max(frag_sizes),
                       cache_path = "~/.cache_peakrefine",
                       cache_version = "v1",
                       force_overwrite = FALSE,
                       n_cores = getOption("mc.cores", 1L),
                       ...){
  if(is.character(query_dt)) query_dt = data.table(file = query_dt, sample = basename(query_dt))
  if(!is.null(query_dt$file)){
    files = query_dt$file
  }else{
    files = query_dt[[1]]
  }
  if(!is.null(query_dt$sample)){
    uniq_names = query_dt$sample
  }else{
    uniq_names = basename(files)
  }
  stopifnot(!any(duplicated(uniq_names)))
  names(files) = uniq_names
  scc_res_l = lapply(files, function(f){
    make_scc_dt.single(f, query_gr,
                       frag_sizes = frag_sizes,
                       fetch_size = fetch_size,
                       cache_path = cache_path,
                       cache_version = cache_version,
                       force_overwrite = force_overwrite,
                       n_cores = n_cores,
                       ...)
  })

  vnames = names(scc_res_l[[1]])

  scc_res = lapply(vnames, function(nam){
    part = lapply(scc_res_l, function(x){
      xv = x[[nam]]
      # if(!is.data.table(xv)){
      #   xv = data.table(xv)
      #   setnames(xv, nam)
      # }
      xv
    })
    dt = rbindlist(part, idcol = "sample")
    merge(dt, query_dt, by = 'sample')
  })
  names(scc_res) = vnames
  scc_res
}

#' make_scc_dt.single
#'
#' based on peakrefine::calcSCCMetrics
#'
#' @param bam_file a single path to a bam file
#' @param query_gr GRanges of regions to calculate SCC for
#' @param frag_sizes optional numeric. Fragment sizes to calculate correlation
#'   at.  The higher the resolution the longer calculation will take.  The
#'   default is to count by 10 from 50 to 350.
#' @param fetch_size optional numeric. Size in bp centered around each interval
#'   in query_gr to retrieve.  Should be greater than max frag_size. The default
#'   is 3*max(frag_sizes).
#' @param cache_path path to cache location for BiocFileCache to use.
#' @param cache_version Modifying the cache version will force recalulation of
#'   all results going forward. Default is v1.
#' @param force_overwrite Logical, if TRUE, cache contents will be overwritten.
#' @param n_cores Number of cores to use. Defaults to mc.cores if set or 1.
#' @param ... passed to Rsamtools::ScanBamParam()
#'
#' @return list fo tidy data.table of SCC data for bam_file
#'
#' @import tools digest
#'
#' @examples
#' peak_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bed$", full.names = TRUE)
#' bam_file = dir(system.file("extdata", package = "seqqc"),
#'   pattern = "test_peaks.bam$", full.names = TRUE)
#'
#' query_gr = seqsetvis::easyLoad_bed(peak_file)[[1]]
#'
#' scc_res = seqqc:::make_scc_dt.single(bam_file, query_gr, seq(50,300, by = 10))
make_scc_dt.single = function(bam_file,
                              query_gr,
                              frag_sizes,
                              fetch_size = 3*max(frag_sizes),
                              cache_path = "~/.cache_peakrefine",
                              cache_version = "v1",
                              force_overwrite = FALSE,
                              n_cores = getOption("mc.cores", 1L),
                              ...){
  bam_md5 = NULL
  qgr_md5 = NULL
  if(is.null(bam_md5)){
    bam_md5 = tools::md5sum(bam_file)
  }
  if(is.null(qgr_md5)){
    qgr_md5 = digest::digest(query_gr)
  }
  if(fetch_size <= max(frag_sizes)){
    stop("fetch_size (", fetch_size, ") must exceed max of frag_sizes (", max(frag_sizes), ").")
  }
  stopifnot(file.exists(bam_file))
  stopifnot(class(query_gr) == "GRanges")
  if(!file.exists(paste0(bam_file, ".bai"))){
    stop("bam_file index not found. expected at ", paste0(bam_file, ".bai"),
         "\ntry running:\nsamtools index ", bam_file)
  }
  stopifnot(n_cores >= 1)

  query_gr = resize(query_gr, fetch_size, fix = 'center')


  bfc_corr = BiocFileCache::BiocFileCache(cache_path, ask = FALSE)
  corr_key = paste(qgr_md5, bam_md5, digest::digest(frag_sizes), fetch_size, cache_version, sep = "_")
  corr_res = bfcif(bfc_corr, corr_key, function(){
    message("cached results not found, gathering correlation info.")
    nper = ceiling(length(query_gr) / n_cores)
    grps = ceiling(seq_along(query_gr)/ nper)
    table(grps)
    # browser()
    rl = getReadLength(bam_file, query_gr)
    lres = parallel::mclapply(unique(grps), function(g){
      k = grps == g
      crossCorrByRle(bam_file, query_gr[k], fragment_sizes = frag_sizes, read_length = rl, ...)
    })
    peak_strand_corr = rbindlist(lres)
    gather_metrics(peak_strand_corr, rl)
  }, force_overwrite = force_overwrite)

  #make everything a data.table
  vnames = names(corr_res)
  corr_res = lapply(vnames, function(nam){
    xv = corr_res[[nam]]
    if(!is.data.table(xv)){
      xv = data.table(xv)
      setnames(xv, nam)
    }
    xv
  })
  names(corr_res) = vnames
  corr_res
}


