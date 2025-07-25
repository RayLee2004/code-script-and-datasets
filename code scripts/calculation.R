# R version: R-4.5.0

########## PART 1: COMMUNITY INDICES CALCULATION ##########

# Read merged datasets from Rdata file
load('merged datasets.Rdata')

# Load required packages
library(vegan)
library(tidyr)
library(dplyr)
library(codyn)
library(openxlsx)
library(tibble)

# Create collection of all taxonomic groups' datasets
all_groups <- list(
  algae        = algae,
  bacteria     = bacteria,
  invertebrate = invertebrate,
  fish         = fish,
  fungi        = fungi,
  insect       = insect,
  zooplankton  = zooplankton
)

# Standardize first column name across all tables
all_groups <- lapply(all_groups, function(df) {
  colnames(df)[1] <- 'species'
  return(df)
})

# Clean environment by removing unnecessary objects
rm(list = setdiff(ls(), 'all_groups'))

# Convert absolute abundance to relative abundance
standardize <- function(x) {
  species <- as.data.frame(x$species)
  colnames(species)[1] <- 'species'
  x <- x[, -1]
  col_sums <- colSums(x)
  x <- sweep(x, 2, col_sums, FUN = '/')
  x <- cbind(species, x)
  return(x)
}
all_groups_standardized <- lapply(all_groups, standardize)

# Transpose data tables
transpose <- function(x) {
  current_df <- as.data.frame(x)
  col_names  <- current_df[, 1]
  current_df <- t(current_df[, -1])
  current_df <- as.data.frame(current_df)
  colnames(current_df) <- col_names
  return(current_df)
}
all_groups_transposed <- lapply(all_groups, transpose)


###### Community Indices For Single Taxonomic Groups ######

# Calculate species richness per group
richness <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    richness_values <- vegan::specnumber(x[[y]])
    temp_df <- data.frame(richness = as.numeric(richness_values))
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}
richness_single_groups <- richness(all_groups_transposed)

# Calculate Shannon diversity per group
diversity <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    diversity_values <- vegan::diversity(x[[y]])
    temp_df <- data.frame(diversity = as.numeric(diversity_values))
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}
diversity_single_groups <- diversity(all_groups_transposed)

# Convert data to long format for synchrony/stability calculations
convert_to_longdata <- function(df) {
  df %>%
    pivot_longer(
      cols      = -1,
      names_to  = c('sampling_site', 'sampling_order'),
      names_sep = '_',
      values_to = 'abundance'
    ) %>%
    rename(species = 1) %>%
    mutate(sampling_order = as.integer(sampling_order))
}
all_groups_longdata <- lapply(all_groups, convert_to_longdata)

# Calculate community synchrony
synchrony <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    synchrony_values <- codyn::synchrony(
      x[[y]],
      species.var   = 'species',
      time.var      = 'sampling_order',
      abundance.var = 'abundance',
      replicate.var = 'sampling_site'
    )
    colnames(synchrony_values)[2] <- substitute(y)
    result_df <- cbind(result_df, synchrony_values[2])
  }
  return(result_df)
}
synchrony_single_groups <- synchrony(all_groups_longdata)

# Calculate community stability
stability <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    stability_values <- codyn::community_stability(
      x[[y]],
      time.var      = 'sampling_order',
      abundance.var = 'abundance',
      replicate.var = 'sampling_site'
    )
    colnames(stability_values)[2] <- substitute(y)
    result_df <- cbind(result_df, stability_values[2])
  }
  return(result_df)
}
stability_single_groups <- stability(all_groups_longdata)


###### Community Indices For Multiple Taxonomic Groups ######

# Merge all taxonomic groups
multi_groups <- bind_rows(all_groups)

multi_groups_transposed <- transpose(multi_groups)
multi_groups_longdata   <- convert_to_longdata(multi_groups)

# Calculate multi-group indices
richness_multi_groups  <- data.frame(multi_groups = vegan::specnumber(multi_groups_transposed))
diversity_multi_groups <- data.frame(multi_groups = vegan::diversity(multi_groups_transposed))
synchrony_multi_groups <- synchrony(list(multi_groups = multi_groups_longdata))
stability_multi_groups <- stability(list(multi_groups = multi_groups_longdata))

# Combine all indices into collections
collection_richness   <- cbind(richness_single_groups, richness_multi_groups)
collection_diversity  <- cbind(diversity_single_groups, diversity_multi_groups)
collection_synchrony  <- cbind(synchrony_single_groups, multi_groups = synchrony_multi_groups[, 2])
collection_stability  <- cbind(stability_single_groups, multi_groups = stability_multi_groups[, 2])

# Export results to Excel
community_indices_collection <- list(
  richness  = collection_richness,
  diversity = collection_diversity,
  synchrony = collection_synchrony,
  stability = collection_stability
)
write.xlsx(community_indices_collection, 'results/community indices.xlsx')


########## PART 2: NETWORK INDICES CALCULATION ##########
###### WARNING: CLEARS ALL OBJECTS IN R ENVIRONMENT
rm(list = ls())

# Load required packages
library(openxlsx)
library(igraph)
library(bipartite)
library(tibble)
library(purrr)
library(Hmisc)

# Load species occurrence data
load('datasets/species occurrences.Rdata')

# Create occurrence list from loaded objects
all_objects <- ls()
occurrences <- list()
for (x in all_objects) {
  if (is.data.frame(get(x))) {
    occurrences[[x]] <- get(x)
  }
}
occurrences <- lapply(occurrences, as.data.frame)
occurrences <- occurrences[!names(occurrences) %in% c('sample2022', 'sample2023', 'sample2024')]
rm(list = setdiff(ls(), 'occurrences'))

# Load interaction network and abundance data
original_network <- read.xlsx('datasets/trophic interaction 0-1 adjacent matrix.xlsx')
total_dataset <- read.xlsx('datasets/total abundance dataset.xlsx')

# Prepare abundance matrix
name_list <- total_dataset$species
total_dataset <- total_dataset[, -(1:10)]
rownames(total_dataset) <- name_list

# Weighted network option (keep when considering relationship weights)
mat <- as.matrix(total_dataset)
res <- rcorr(t(mat), type = "pearson")
R   <- (1 + res$r) / 2  # Scale correlation coefficients to [0,1]
diag(R) <- 0
original_network <- original_network[, -(1:2)] * R

# Prepare final network matrix
diag(original_network) <- 0
original_network <- as.data.frame(original_network)
rownames(original_network) <- colnames(original_network)

# Initialize storage objects
dfs      <- list()
matrices <- list()
graphs   <- list()

# Extract subnetwork data frames
get_df <- function(x) {
  occurrence_list <- x[, ncol(x)]
  current_df <- original_network[occurrence_list, occurrence_list]
  rownames(current_df) <- colnames(current_df)
  dfs <<- c(dfs, list(current_df))
}
lapply(occurrences, get_df)

# Convert to matrices and graphs
get_matrix_graph <- function(x) {
  matrix <- as.matrix(x)
  graph  <- igraph::graph_from_adjacency_matrix(
    matrix, 
    mode = "directed", 
    weighted = TRUE, 
    diag = FALSE
  )
  # Retain only connected nodes
  none_isolated_nodes <- igraph::V(graph)[igraph::degree(graph) > 0]
  graph2 <- igraph::induced_subgraph(graph, vids = none_isolated_nodes)
  matrices <<- c(matrices, list(matrix))
  graphs   <<- c(graphs, list(graph2))
}
lapply(dfs, get_matrix_graph)

# Name the network objects
names(dfs)      <- names(occurrences)
names(matrices) <- names(occurrences)
names(graphs)   <- names(occurrences)

# Calculate chain-based network indices
chain_indices <- function(x) {
  result_df <- data.frame(
    sampling_site_order = character(0),
    connectance         = numeric(0),
    max_path_length     = numeric(0),
    mean_path_length    = numeric(0),
    no_nodes            = numeric(0)
  )
  for (y in names(x)) {
    current_graph <- x[[y]]
    metrics <- data.frame(
      sampling_site_order = y,
      connectance         = edge_density(current_graph),
      max_path_length     = diameter(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      mean_path_length    = mean_distance(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      no_nodes            = sum(igraph::degree(current_graph) > 0)
    )
    result_df <- rbind(result_df, metrics)
  }
  return(result_df)
}
chain_indices <- chain_indices(graphs)

# Calculate structural network indices
structure_indices <- function(x, y) {
  # Modularity calculation
  result_df1 <- data.frame(sampling_site_order = character(0), modularity = numeric(0))
  for (z in names(y)) {
    current_community <- cluster_walktrap(y[[z]])
    mod_df <- data.frame(
      sampling_site_order = z,
      modularity          = modularity(current_community, weights = E(current_community)$weight)
    )
    result_df1 <- rbind(result_df1, mod_df)
  }
  
  # Nestedness and robustness metrics
  result_df2 <- data.frame(nestedness = numeric(0), vulnerability = numeric(0), robustness = numeric(0))
  
  vulnerability <- function(n) {
    gen_vul <- networklevel(n, index = "vulnerability", weighted = TRUE)
    gen_vul[!grepl("generality", names(gen_vul))]
  }
  
  for (w in names(x)) {
    current_matrix <- x[[w]]
    net_metrics <- data.frame(
      nestedness    = nested(current_matrix, "wine"),
      vulnerability = vulnerability(current_matrix),
      robustness    = robustness(
        second.extinct(current_matrix, participant = "higher", method = "random", nrep = 10, details = FALSE)
      )
    )
    result_df2 <- rbind(result_df2, net_metrics)
  }
  return(cbind(result_df1, result_df2))
}
network_indices <- structure_indices(matrices, graphs)

# Combine and export network indices
collection_network_indices <- cbind(chain_indices, network_indices[, -1])
write.xlsx(collection_network_indices, 'results/meta-foodweb indices.xlsx')
