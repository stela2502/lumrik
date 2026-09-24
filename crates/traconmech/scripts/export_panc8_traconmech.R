#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratData)
    library(Matrix)
})

outdir <- "/tmp/panc8_traconmech"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading panc8...\n")
panc8 <- LoadData("panc8")

expression <- panc8[["RNA"]]$data

cat(sprintf(
    "Expression matrix: %d genes x %d cells; %d non-zero entries\n",
    nrow(expression), ncol(expression), length(expression@x)
))

matrix_file <- file.path(outdir, "expression.mtx")
genes_file <- file.path(outdir, "genes.tsv")
cells_file <- file.path(outdir, "cells.tsv")
info_file <- file.path(outdir, "README.txt")

writeMM(expression, matrix_file)

write.table(
    rownames(expression), genes_file,
    quote = FALSE, row.names = FALSE, col.names = FALSE
)

meta <- data.frame(
    cell = colnames(expression),
    celltype = panc8$celltype,
    tech = panc8$tech,
    dataset = panc8$dataset,
    replicate = panc8$replicate,
    stringsAsFactors = FALSE
)

stopifnot(identical(meta$cell, colnames(expression)))

write.table(
    meta, cells_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

info <- c(
    "TraConMech panc8 test dataset",
    "",
    sprintf("Genes: %d", nrow(expression)),
    sprintf("Cells: %d", ncol(expression)),
    sprintf("Non-zero matrix entries: %d", length(expression@x)),
    "",
    "expression.mtx:",
    "  Matrix Market sparse matrix.",
    "  Rows correspond exactly to genes.tsv.",
    "  Columns correspond exactly to cells.tsv.",
    "",
    "IMPORTANT:",
    "  The panc8 SeuratData object was converted to Seurat v5.",
    "  RNA counts and data layers were identical in the inspected object.",
    "  Values include non-integers.",
    "  Therefore expression.mtx must NOT be interpreted as raw UMI counts.",
    "",
    "Cell types:",
    paste(capture.output(print(table(meta$celltype))), collapse = "\n"),
    "",
    "Technologies:",
    paste(capture.output(print(table(meta$tech))), collapse = "\n")
)

writeLines(info, info_file)

cat("Compressing output...\n")
status <- system2("gzip", c("-f", matrix_file, genes_file, cells_file))
if (status != 0) stop("gzip failed")

cat("\nDone. Output directory:\n", outdir, "\n\n", sep = "")
print(file.info(list.files(outdir, full.names = TRUE))[, "size", drop = FALSE])
cat("\nCell types:\n")
print(table(meta$celltype))
cat("\nTechnologies:\n")
print(table(meta$tech))
