library("mda")
library("scPred")
library("Seurat")
library("magrittr")

reference <- scPred::pbmc_1
query <- scPred::pbmc_2

reference <- reference %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)

get_probabilities(reference) %>% head()
get_scpred(reference)
# plot_probabilities(reference)

reference <- trainModel(reference, model = "mda", reclassify = c("cMono", "ncMono"))
get_scpred(reference)
plot_probabilities(reference)


query <- query %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

query <- RenameAssays(query, "RNA" = "data")
reference <- RenameAssays(reference, "RNA" = "data")

query <- scPredict(query, reference)

DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")

# query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)

# We can compare the results with the original labels
DimPlot(query, group.by = "cell_type", label = TRUE, repel = TRUE)

crossTab(query, "cell_type", "scpred_prediction")

crossTab(query, "cell_type", "scpred_prediction", output = "prop")
