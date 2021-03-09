#' write_bed_frip
#'
#' @param query_gr GRanges of regions used to make frip_dt and assign_dt
#' @param assign_dt data.table with "id" and "cluster_id" that maps region ids to clusters. In output from \code{\link{plot_signals}}.
#' @param frip_dt output from \code{\link{make_frip_dt}}
#' @param file File to write to
#'
#' @return
#' @export
#'
#' @examples
write_bed_frip = function(query_gr, assign_dt, frip_dt, file = "regions_with_FRIP.txt"){
  region_frip_dt = dcast(frip_dt[, .(id, name, N, treatment, frip)], id~name, value.var = "frip")
  setkey(region_frip_dt, id)
  setkey(assign_dt, id)
  bed_frip_towrite = query_gr

  mcols(bed_frip_towrite) = region_frip_dt[.(names(bed_frip_towrite))]
  # bed_frip_towrite$id = NULL

  message("write files...")
  bed_frip_towrite$cluster_id = assign_dt[.(names(bed_frip_towrite))]$cluster_id
  # rtracklayer::export.bed(bed_frip_towrite, paste0("consensus_regions_with_FRIP.", file_tag, ".bed"))

  bed_frip_towrite_dt = as.data.table(bed_frip_towrite)
  bed_frip_towrite_dt$strand = "."
  bed_frip_towrite_dt$score = 0
  extra_cols = setdiff(colnames(mcols(bed_frip_towrite)), "id")
  extra_cols = gsub("-", ".", extra_cols)
  fwrite(bed_frip_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE],
         file = file,
         sep = "\t")
  invisible(bed_frip_towrite_dt)
}

#' write_bed_overlaps
#'
#' @param query_gr GRanges of regions used to make assign_dt
#' @param assign_dt data.table with "id" and "cluster_id" that maps region ids to clusters. In output from \code{\link{plot_signals}}.
#' @param file File to write to
#'
#' @return
#' @export
#'
#' @examples
write_bed_overlaps = function(query_gr, assign_dt, file = "regions_with_overlaps.txt"){
  setkey(assign_dt, id)
  bed_peak_towrite = query_gr
  bed_peak_towrite$cluster_id = assign_dt[.(names(bed_peak_towrite))]$cluster_id
  bed_peak_towrite_dt = as.data.table(bed_peak_towrite)
  bed_peak_towrite_dt$strand = "."
  bed_peak_towrite_dt$score = 0
  extra_cols = setdiff(colnames(mcols(bed_peak_towrite)), "id")
  extra_cols = gsub("-", ".", extra_cols)
  fwrite(bed_peak_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE],
         file = paste0("consensus_regions_with_peak_call.", file_tag, ".txt"),
         sep = "\t")
}
