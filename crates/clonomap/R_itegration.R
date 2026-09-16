library(tools)

pca_tree <- function(seqs, k = 30, prefix, plots = TRUE) {

    input_file  <- paste0(prefix, "_input.txt")
    coords_file <- paste0(prefix, "_pca.tsv")
    edges_file  <- paste0(prefix, "_tree.tsv")
    pca_plot    <- paste0(prefix, "_pca.png")
    tree_plot   <- paste0(prefix, "_tree.png")

    writeLines(seqs, input_file)

    cmd <- "pca_tree"
    args <- c(
        input_file,
        "--k", k,
        "--coords", coords_file,
        "--edges", edges_file
    )

    if (plots) {
        args <- c(args,
            "--plot-pca", pca_plot,
            "--plot-tree", tree_plot
        )
    }

    # Capture both stdout and stderr
    out <- system2(cmd, args,
                   stdout = TRUE,
                   stderr = TRUE,
                   wait = TRUE)

    status <- attr(out, "status")

    if (!is.null(status) && status != 0) {
        cat("\n❌ ERROR: Rust binary 'pca_tree' exited with code", status, "\n")
        cat("------------ Rust error output ------------\n")
        cat(out, sep="\n")
        cat("------------ end Rust error ---------------\n")
        stop("pca_tree failed — see Rust error above.")
    }

    return(list(
        coords_file = coords_file,
        edges_file  = edges_file,
        pca_plot    = if (plots) pca_plot else NULL,
        tree_plot   = if (plots) tree_plot else NULL
    ))
}

assign_clone_length_groups <- function(df, length_ratio_threshold = 0.20) {
    stopifnot("clone_id" %in% names(df))
    stopifnot("sequence_alignment" %in% names(df))

    len = nchar(df[,'sequence_alignment'])
    df$"clone_id_length"= paste( sep="_len_", df$"clone_id", len)

    return(df)
}


base_dir <- "/home/med-sal/sens05_shared/jyuan/no_backup/GM/Stefans_analysis/ChangeO_Db_2025/DefineClones_2025"

files <- list.files(
  path = base_dir,
  pattern = "2025_ProductiveCloneDfined\\.tsv$",
  full.names = TRUE
)

print(files)

total_start <- Sys.time()   # start timing everything

for (file in files) {

    file_start <- Sys.time()  # start timing this file

    cat("\n------------------------------------------------------------\n")
    cat("Processing file:", file, "\n")

    df <- read.delim(file)

    df <- assign_clone_length_groups(df)

    clones <- table(df[,'clone_id_length'])
    clones <- names(clones[clones > 200])

    prefix_base <- file_path_sans_ext(file)

    for (clone_id in clones) {

        clone_start <- Sys.time()  # start timing clone

        n_cells <- sum(df[,'clone_id_length'] == clone_id)
        cat(sprintf("   - clone %s (n = %d)\n", clone_id, n_cells))

        strs <- df[df[,'clone_id_length'] == clone_id, 'sequence_alignment']

        # prefix for this clone
        clone_prefix <- paste0(prefix_base, "_clone_", clone_id)

        pca_tree(
            strs,
            k = 30,
            prefix = clone_prefix,
            plots = TRUE
        )

        clone_end <- Sys.time()
        cat(sprintf("      clone time: %.2f sec\n",
                    as.numeric(difftime(clone_end, clone_start, units = "secs"))))
    }

    file_end <- Sys.time()
    cat(sprintf("File completed in: %.2f sec\n",
                as.numeric(difftime(file_end, file_start, units = "secs"))))
}

total_end <- Sys.time()

cat("\n============================================================\n")
cat(sprintf("TOTAL RUNTIME: %.2f seconds (%.2f minutes)\n",
            as.numeric(difftime(total_end, total_start, units="secs")),
            as.numeric(difftime(total_end, total_start, units="mins"))))
cat("============================================================\n\n")
