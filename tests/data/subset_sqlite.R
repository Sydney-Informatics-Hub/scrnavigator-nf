library(ensembldb)
library(RSQLite)
library(Seurat)

# Download the full EnsDb sqlite using the same script the pipeline uses.
system2("Rscript", args = c("scrnavigator-nf/bin/download_ensdb.R", "human", "v113", ".cache"))
src <- "scrnavigator-nf/tests/data/EnsDb_hsapiens_v113.sqlite"

# Collect every unique gene identifier present across all test RDS files
genes_in_test_data <- list.files("scrnavigator-nf/tests/data/rds", pattern = "\\.Rds$", full.names = TRUE) |>
  lapply(\(f) rownames(readRDS(f))) |>
  unlist() |>
  unique()
message(sprintf("Genes in test data: %d", length(genes_in_test_data)))

# Copy the full sqlite to the test fixture destination so the original is
# left untouched. All subsetting is done on the copy.
dst <- "scrnavigator-nf/tests/data/EnsDb_hsapiens_test.sqlite"
file.copy(src, dst, overwrite = TRUE)

con <- dbConnect(SQLite(), dst)
genes_to_keep <- data.frame(lookup = genes_in_test_data)
dbWriteTable(con, "genes_to_keep", genes_to_keep, temporary = TRUE)

dbExecute(con, "DELETE FROM gene    WHERE gene_id  NOT IN (SELECT lookup FROM genes_to_keep)
                                      AND gene_name NOT IN (SELECT lookup FROM genes_to_keep)")
dbExecute(con, "DELETE FROM tx      WHERE gene_id  NOT IN (SELECT gene_id  FROM gene)")
dbExecute(con, "DELETE FROM tx2exon WHERE tx_id    NOT IN (SELECT tx_id    FROM tx)")
dbExecute(con, "DELETE FROM exon    WHERE exon_id  NOT IN (SELECT exon_id  FROM tx2exon)")

# 5. Reclaim freed pages to shrink the file, then close.
dbExecute(con, "VACUUM")
dbDisconnect(con)

# 6. Validate: every gene in the test data must appear in the subset db.
#    Genes absent from Ensembl (e.g. non-standard contigs) are expected and
#    reported but do not fail the script.
con <- dbConnect(SQLite(), dst)
db_ids   <- dbGetQuery(con, "SELECT gene_id   FROM gene")$gene_id
db_names <- dbGetQuery(con, "SELECT gene_name FROM gene")$gene_name
dbDisconnect(con)

found     <- genes_in_test_data %in% db_ids | genes_in_test_data %in% db_names
n_found   <- sum(found)
n_total   <- length(genes_in_test_data)
n_missing <- n_total - n_found

message(sprintf("Validation: %d / %d genes found in subset db", n_found, n_total))
if (n_missing > 0) {
  message(sprintf("  %d genes not found e.g.:", n_missing))
  message(paste(" ", head(genes_in_test_data[!found], 10), collapse = "\n"))
}

message(sprintf("Done. Subset sqlite: %.1f MB", file.info(dst)$size / 1e6))
