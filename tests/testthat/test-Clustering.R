test_that("cluster function performs K-means clustering as expected", {
  # Create an expanded sample dataset with more samples
  datafile <- data.frame(
    Locus = rep(c(1, 2), each = 9),
    Locus_allele = rep(c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2), 2),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0),
    Sample4 = c(NA, 25, 0.2, 0.3, 0.5, NA, 15, 0.1, 0.2, NA, 25, 0.2, 0.3, 0.5, NA, 15, 0.1, 0.2),
    Sample5 = c(NA, 35, 0.3, 0.4, 0.3, NA, 5, 0.3, 0.4, NA, 35, 0.3, 0.4, 0.3, NA, 5, 0.3, 0.4),
    Sample6 = c(NA, 40, 0.1, 0.2, 0.4, NA, 20, 0.2, 0.3, NA, 40, 0.1, 0.2, 0.4, NA, 20, 0.2, 0.3),
    Sample7 = c(NA, 50, 0.2, 0.3, 0.4, NA, 30, 0.3, 0.4, NA, 50, 0.2, 0.3, 0.4, NA, 30, 0.3, 0.4),
    Sample8 = c(NA, 45, 0.1, 0.2, 0.5, NA, 25, 0.2, 0.3, NA, 45, 0.1, 0.2, 0.5, NA, 25, 0.2, 0.3)
  )

  # Call the cluster function with the sample dataset and K=2
  clustering_result <- cluster(datafile, K=2)

  # Test that the function returns a list
  expect_true(is.list(clustering_result))

  # Test that the list contains the expected elements
  expect_true(all(c("PopFactor", "clust", "dat") %in% names(clustering_result)))

  # Test that the PopFactor is a factor
  expect_true(is.factor(clustering_result$PopFactor))

  # Test that the clust is an integer vector with length equal to number of samples
  expect_true(is.integer(clustering_result$clust))
  expect_equal(length(clustering_result$clust), ncol(datafile) - 4) # Number of samples

  # Test specific values in the clustering result
  # Since K-means clustering results can vary, we will not test specific cluster assignments
  expect_true(all(clustering_result$clust > 0 & clustering_result$clust <= 2)) # All clusters should be between 1 and K
})

