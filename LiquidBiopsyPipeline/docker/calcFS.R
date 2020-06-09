# Collect Family Size and other Depth Metrics
# This produces several files that are intended to be used as output
# in cromwell to a Firecloud, or Terra workspace.

library(dplyr)
library(data.table)

calculateMedianFromHistogram <- function(df, size_column, count_column) {
	df$cumulative = cumsum(df[count_column])
	midpoint = sum(df[count_column]) / 2
	median = df[nrow(subset(df, cumulative <= midpoint)) + 1, size_column]
	median
}

getDuplexFamilySizeMetrics <- function(fileName) {
	metrics = read.table(fileName, header = T)
	metrics[, "total_size"] = metrics[, "ab_size"] + metrics[, "ba_size"]

	total_family_sizes = as.data.frame(metrics %>% group_by(total_size) %>% summarise(total_size_count = sum(count)))

	median_total_family_size = calculateMedianFromHistogram(total_family_sizes, "total_size", "total_size_count")

	write(median_total_family_size, "median_total_family_size.txt")

	ab_family_sizes = as.data.frame(metrics %>% group_by(ab_size) %>% summarise(ab_size_count = sum(count)))
	median_ab_family_size = calculateMedianFromHistogram(ab_family_sizes, "ab_size", "ab_size_count")
	write(median_ab_family_size, "median_ab_family_size.txt")

	ba_family_sizes = as.data.frame(metrics %>% group_by(ba_size) %>% summarise(ba_size_count = sum(count)))

	median_ba_family_size = calculateMedianFromHistogram(ba_family_sizes, "ba_size", "ba_size_count")
	write(median_ba_family_size, "median_ba_family_size.txt")

	mean_family_size = sum(metrics$total_size * metrics$fraction)
	write(mean_family_size, "mean_total_family_size.txt")

	mean_ab_size = sum(metrics$ab_size * metrics$fraction)
	write(mean_ab_size, "mean_ab_family_size.txt")

	mean_ba_size = sum(metrics$ba_size * metrics$fraction)
	write(mean_ba_size, "mean_ba_family_size.txt")
}

