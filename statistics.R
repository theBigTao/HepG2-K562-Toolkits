library(data.table)

# Define command-line arguments
commandArgs <- commandArgs(trailingOnly = TRUE)
root_dir <- commandArgs[1]
file1_name <- commandArgs[2]
file2_name <- commandArgs[3]
out_file <- commandArgs[4]

calculate_stats <- function(directory_path, file1_name, file2_name) {

    # Create an empty data table to store the results
    results <- data.table(cell_line = character(),
                          tf1 = character(),
                          tf2 = character(),
                          log2_fold_change = numeric(),
                          pvalue_paired_t_test = numeric(),
                          pvalue_paired_wilcox_test = numeric())

    # Iterate over all subdirectories
    for (subdir in list.dirs(directory_path, full.names = TRUE, recursive=FALSE)) {

        # Skip the top-level directory
        if (subdir == directory_path) next

        # Get the file paths for the two files in this subdirectory
        file1_path <- file.path(subdir, file1_name)
        file2_path <- file.path(subdir, file2_name)

        # Get TF1, TF2, Cellline
        cell_line <- strsplit(subdir, "_")[[1]][2]
        tf1 <- strsplit(subdir, "_")[[1]][3]
        tf2 <- strsplit(subdir, "_")[[1]][4]

        # Read the two files into data tables
        cat (file1_path, "\n")
        file1 <- fread(file1_path)
        file2 <- fread(file2_path)

        # Calculate the log2-fold change of the medians
        log2_fc <- log2(median(file2$V1) / median(file1$V1))

        # Calculate the paired t-test p-value
        pvalue_t <- t.test(file1$V1, file2$V1, var.equal=FALSE, paired = TRUE)$p.value
        pvalue_t <- ifelse(is.na(pvalue_t), 1, pvalue_t)

        # Calculate the paired Wilcoxon rank sum test p-value
        pvalue_wilcox <- wilcox.test(file1$V1, file2$V1, var.equal=FALSE, paired = TRUE)$p.value
        pvalue_wilcox <- ifelse(is.na(pvalue_wilcox), 1, pvalue_wilcox)

        # Add the results to the data table
        results <- rbind(results, data.table(cell_line = cell_line,
                                             tf1 = tf1,
                                             tf2 = tf2,
                                             log2_fold_change = log2_fc,
                                             pvalue_paired_t_test = pvalue_t,
                                             pvalue_paired_wilcox_test = pvalue_wilcox))

    }

    results[, "qvalues_paired_t_test"] <- p.adjust(results[, pvalue_paired_t_test], method = "BH")
    results[, "qvalue_paired_wilcox_test"] <- p.adjust(results[, pvalue_paired_wilcox_test], method = "BH")

    return (results)
}

results <- calculate_stats(root_dir, file1_name, file2_name)
fwrite(results, out_file, sep="\t")
