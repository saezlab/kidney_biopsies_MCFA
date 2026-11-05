# (C) Sergio Alejandro Gómez
# Small helpers shared by the genetics scripts: command-line arguments, table I/O and logging.

suppressPackageStartupMessages(library(data.table))

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (a in args) {
    a <- sub("^--", "", a)
    key <- sub("=.*$", "", a)
    out[[key]] <- if (grepl("=", a)) sub("^[^=]*=", "", a) else TRUE
  }
  out
}

arg <- function(args, key, default = NULL) {
  value <- args[[key]]
  if (is.null(value)) {
    if (is.null(default)) stop("missing argument --", key, call. = FALSE)
    return(default)
  }
  value
}

arg_true <- function(args, key, default = "true") tolower(arg(args, key, default)) %in% c("true", "1", "yes")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

read_tsv <- function(path, ...) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", ...)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  log_msg("wrote ", path, " (", NROW(df), " rows)")
  invisible(path)
}

save_rds <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  log_msg("wrote ", path)
  invisible(path)
}
