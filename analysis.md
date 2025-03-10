# Single-Cell RNA Sequencing Analysis Workflow

## 1. Data Preprocessing

### 1.1 Organize Data Directory

```bash
# Define the data directory
data_dir="/public/home/99017/02users/07ca/coursework/GSE121638_RAW"

# Get all file names
files=$(ls "${data_dir}")

# Loop through files
for file in ${files}; do
    # Extract sample identifier (e.g., GU0744_P)
    sample=$(echo "${file}" | grep -oP 'GU[0-9]+_[PT]')

    # If sample identifier is successfully extracted
    if [[ -n "${sample}" ]]; then
        # Create corresponding directory (if it doesn't exist)
        mkdir -p "${data_dir}/${sample}"

        # Move files to the corresponding directory
        mv "${data_dir}/${file}" "${data_dir}/${sample}/"
    fi
done

# Print the organized directory structure
echo "Organized directory structure:"
tree "${data_dir}"

```
### 1.2 Data preanalysis
```r
library(Seurat)

# Define the data directory
data_dir <- "/public/home/99017/02users/07ca/coursework/GSE121638_RAW"
samples <- c("GU0700_T", "GU0715_T", "GU0744_T",
             "GU0700_P", "GU0715_P", "GU0744_P")

# Create an empty list to store Seurat objects for each sample
seurat_list <- list()

# Loop through each sample to read data and create Seurat objects
for (sample in samples) {
  sample_dir <- file.path(data_dir, sample)
  data <- Read10X(data.dir = sample_dir)  # Read 10X data
  seurat_obj <- CreateSeuratObject(counts = data, project = sample)  # Create Seurat object
  seurat_list[[sample]] <- seurat_obj  # Store the Seurat object in the list
}

# Add sample names to cell barcodes to avoid conflicts after merging
for (sample_name in names(seurat_list)) {
  seurat_obj <- seurat_list[[sample_name]]
  new_cellnames <- paste0(sample_name, "_", colnames(seurat_obj))  # Append sample name to cell barcodes
  colnames(seurat_obj) <- new_cellnames  # Update cell barcodes
  seurat_list[[sample_name]] <- seurat_obj  # Update the Seurat object
}

# Merge all samples into a single Seurat object
sc <- merge(seurat_list[[1]], seurat_list[-1]])

# Print the merged Seurat object
print(sc)
seurat_list
$GU0700_T

An object of class Seurat 
33694 features across 4859 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

$GU0715_T
An object of class Seurat 
33694 features across 4680 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

$GU0744_T
An object of class Seurat 
33694 features across 2715 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

$GU0700_P
An object of class Seurat 
33694 features across 5850 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

$GU0715_P
An object of class Seurat 
33694 features across 5125 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

$GU0744_P
An object of class Seurat 
33694 features across 2459 samples within 1 assay 
Active assay: RNA (33694 features, 0 variable features)
 2 layers present: counts, data

```
### 1.3 Data visualize
```r
# Add sample information to metadata
sc$sample <- sc$orig.ident

# Calculate mitochondrial gene percentage
mito_genes <- grep(pattern = "^MT-", x = rownames(sc), value = TRUE)  # Identify mitochondrial genes
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")  # Calculate percentage

# Calculate ribosomal gene percentage
rb.genes <- rownames(sc)[grep("^RP[SL]", rownames(sc))]  # Identify ribosomal genes
sc[["percent.ribo"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")  # Calculate percentage

# Filter cells based on quality control metrics
sc <- subset(sc, subset = nFeature_RNA > 300 & nFeature_RNA < 20000 & percent.ribo < 30 & percent.mt < 20)
# Normalize data
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
# Scale data
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
# Run PCA
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
# Run Harmony for batch correction
library(harmony)
sc <- RunHarmony(sc, reduction = "pca", group.by.vars = "sample", reduction.save = "harmony")
# Determine the optimal number of principal components (PCs)
pct <- sc[["harmony"]]@stdev / sum(sc[["harmony"]]@stdev) * 100  # Percentage of variance explained by each PC
cumu <- cumsum(pct)  # Cumulative variance
co1 <- which(cumu > 90 & pct < 5)[1]  # First PC where cumulative variance > 90% and individual variance < 5%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1  # PC where the drop in variance is significant
pcs <- min(co1, co2)  # Select the minimum number of PCs
bestpc <- 1:pcs  # Define the range of PCs to use
# Run UMAP
sc <- RunUMAP(sc, reduction = "harmony", dims = bestpc, reduction.name = "umap", min.dist = 0.4)
# Find neighbors and cluster cells
sc <- FindNeighbors(sc, reduction = "harmony", dims = bestpc)
sc <- FindClusters(sc, resolution = seq(0.1, 0.5, by = 0.1))  # Test multiple resolutions
# Visualize clustering tree
library(clustree)
pdf("Harmony_clustree.pdf")
clustree(sc)  # Plot clustering tree
dev.off()
# Final clustering at resolution 0.5
sc <- FindClusters(sc, resolution = 0.5)
# Visualize clusters using UMAP
library(scCustomize)
library(ggsci)
pal <- pal_d3('category20')(12)  # Define color palette
p3 <- DimPlot_scCustom(seurat_object = sc, group.by = "seurat_clusters", colors_use = pal, figure_plot = TRUE, label = TRUE)
ggsave("Clusters_umap.pdf", p3)
# Visualize sample distribution using UMAP
p3 <- DimPlot_scCustom(seurat_object = sc, group.by = "sample", figure_plot = TRUE, label = TRUE)
ggsave("Clusters_sample.pdf", p3)
```

![相对路径](D:/22222/0课题/output/fig0.png)
### 1.4 Marker Gene featureplot
```r
# Define marker genes for cell type annotation
genelist <- list(
  Monocyte = c("CD14", "LYZ", "S100A8", "S100A9", "FCGR3A"),
  Macrophage = c("CD163", "MSR1", "FCGR1A", "C1QA", "C1QB", "C1QC"),
  DendriticCell = c("CD1C", "FCER1A", "CLEC9A", "CD83", "CCR7"),
  NKCell = c("NKG7", "GNLY", "KLRD1", "KLRB1", "CD16"),
  CD8T = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "CCL5", "CXCR3"),
  CD4T = c("CD4", "IL7R", "FOXP3"),
  Effect = c("PDCD1", "HAVCR2", "TNFRSF9"),
  BCell = c("CD19", "CD79A", "CD79B", "MS4A1"),
  PlasmaCell = c("JCHAIN", "IGHA1", "IGHG1"),
  Mast = c("CCR3", "CPA3"),
  Neutrophil = c("FCGR3B", "CXCR2")
)
# Dot plot of marker genes
p1 <- DotPlot(sc, features = genelist, cols = c("grey", "red")) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin = margin(b = 3, unit = "mm")),
    axis.line = element_blank()
  ) +
  labs(x = "", y = "")
ggsave("dotplot_celltype.pdf", p1, width = 12, height = 5)
# Annotate cell types
celltype <- c("NK cell", "Monocyte", "Cytotoxic CD8+ T cell", "Macrophage", "CD8+ T cell", "Effector T cell", "Monocyte", "DC", "Proliferative", "CD4+ T cell", 
              "NK cell", "DC", "Plasma", "DC", "B cell", "Monocyte", "DC", "Mast cell")
names(celltype) <- levels(sc)
sc <- RenameIdents(sc, celltype)  # Rename clusters based on cell type
sc@meta.data$celltype <- Idents(sc)  # Add cell type information to metadata
Idents(sc) <- sc@meta.data$celltype  # Set cell type as active identity
# Visualize annotated cell types
p3 <- DimPlot_scCustom(seurat_object = sc, group.by = "celltype", colors_use = pal, figure_plot = TRUE, label = TRUE)
ggsave("Celltype_umap.pdf", p3)
# Feature plot of marker genes
gene <- c("CD3E", "CD8A", "CD4", "IL7R", "CD79A", "MS4A1", "CD14", "FCGR3A", "NKG7", "CPA3", "JCHAIN", "CD1C")
plots <- list()
for (i in 1:length(gene)) {
  plots[[i]] <- FeaturePlot_scCustom(seurat_object = sc, colors_use = colorRampPalette(c("#3288BD", "white", "#D53E4F"))(50), features = gene[i]) + NoAxes()
}
```

![相对路径](D:/22222/0课题/output/fig3.png)
![相对路径](D:/22222/0课题/output/fig1.png)
![相对路径](D:/22222/0课题/output/fig4.png)
### 1.5 Sample cell ratio
```r
# Combine feature plots using patchwork
library(patchwork)
p <- wrap_plots(plots, ncol = 3)
ggsave(p, file = "featureplot3.pdf", width = 8, height = 9)
# Define a function to create cell ratio plots
create_cell_ratio_plot <- function(sc, celltype, output_path) {
  cellratio <- table(sc@meta.data[[celltype]], sc@meta.data$sample) %>% as.data.frame()
  colnames(cellratio) <- c("celltype", "sample", "Counts")
  cellratio$celltype <- factor(cellratio$celltype)
  col <- pal
  cell_ratio <- ggplot(data = cellratio, aes(x = Counts, y = sample, fill = celltype)) +
    geom_bar(stat = "identity", width = 0.8, position = "fill") +
    scale_fill_manual(values = col) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Ratio", y = "") +
    theme(axis.text.y = element_text(size = 18, colour = "black")) +
    theme(axis.text.x = element_text(size = 18, colour = "black")) +
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
  # Save the plot
  ggsave(file.path(output_path, "Cell_ratio_sample.pdf"), egg::set_panel_size(cell_ratio, width = unit(6, "in"), height = unit(3, "in")), 
         width = 12, height = 6, units = "in", dpi = 300)
}
# Create and save cell ratio plot
create_cell_ratio_plot(sc, "celltype", "/public/home/99017/02users/07ca/coursework/GSE121638_RAW/output")
```
![相对路径](D:/22222/0课题/output/fig5.png)
### 1.5 Doheatmap plot
```r
names(sc) <- levels(sc)
sc <- RenameIdents(sc, seurat_clusters)
sc.markers.all <- Seurat::FindAllMarkers(sc,
                               only.pos = TRUE,
                               logfc.threshold = 0.5)
sc.markers <- sc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)
head(sc.markers)
pdf('sc1.pdf',height = 6,width = 8,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:12))
dev.off()
```
![相对路径](D:/22222/0课题/output/fig6.png)
### 1.6 CD8+ subtype
```r
# Subset CD8+ T cells for further analysis
sc2 <- subset(x = sc, idents = c("Cytotoxic CD8+ T cell", "CD8+ T cell", "Effector T cell"))
# Re-run analysis on the subset
sc2 <- NormalizeData(sc2, normalization.method = "LogNormalize", scale.factor = 10000)
sc2 <- FindVariableFeatures(sc2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc2)
sc2 <- ScaleData(sc2, features = all.genes)
sc2 <- RunPCA(sc2, features = VariableFeatures(object = sc2))
sc2 <- RunHarmony(sc2, reduction = "pca", group.by.vars = "sample", reduction.save = "harmony")
pct <- sc2[["harmony"]]@stdev / sum(sc2[["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs <- min(co1, co2)
bestpc <- 1:pcs
sc2 <- RunUMAP(sc2, reduction = "harmony", dims = bestpc, reduction.name = "umap", min.dist = 0.4)
sc2 <- FindNeighbors(sc2, reduction = "harmony", dims = bestpc)
sc2 <- FindClusters(sc2, resolution = seq(0.1, 0.5, by = 0.1))
# Visualize clustering tree for the subset
pdf("Harmony_clustree2.pdf")
clustree(sc2)
dev.off()
# Final clustering at resolution 0.5
sc2 <- FindClusters(sc2, resolution = 0.5)
# Visualize clusters for the subset
p3 <- DimPlot_scCustom(seurat_object = sc2, group.by = "seurat_clusters", figure_plot = TRUE, label = TRUE)
ggsave("cd8_cluster_umap.pdf", p3)
```
![相对路径](D:/22222/0课题/output/fig7.png)
### 1.6 CD8+ subtype slingshot
```r
# Perform trajectory analysis using Slingshot
library(slingshot)
sce <- as.SingleCellExperiment(sc2, assay = "RNA")
sce_slingshot1 <- slingshot(sce, reducedDim = "UMAP", clusterLabels = sce$seurat_clusters, approx_points = 150)

# Define color palette for visualization
cell_pal <- function(cell_vars, pal_fun, ...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(sce_slingshot1$seurat_clusters, brewer_pal("qual", "Set2"))

# Plot trajectory
pdf("clusters.pdf", height = 10, width = 10)
plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch = 16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce_slingshot1), lwd = 2, col = "black")

# Add cluster labels
celltype_label <- sc2@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(seurat_clusters = sc2@meta.data$seurat_clusters) %>%
  group_by(seurat_clusters) %>%
  summarise(UMAP1 = median(UMAP_1), UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$seurat_clusters[i], x = celltype_label$UMAP1[i] - 1, y = celltype_label$UMAP2[i])
}
dev.off()
```
![相对路径](D:/22222/0课题/output/fig8.png)
### 1.7 CD8+ cluster1 KEGG pathway
```r
# Find marker genes for a specific cluster
marker1 <- FindMarkers(sc2, ident.1 = "1", ident.2 = NULL, only.pos = TRUE)

# Perform KEGG pathway enrichment analysis
up <- rownames(marker1)
up_ENTREZID <- mapIds(x = org.Hs.eg.db, keys = up, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
KEGG_diff <- enrichKEGG(gene = up_ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
KEGG_diff <- setReadable(KEGG_diff, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Calculate enrichment metrics
KEGG_diff2 <- mutate(KEGG_diff, RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio)))

# Save KEGG results
write.csv(KEGG_diff2@result, file = "KEGG_diff.csv")

# Plot top 10 enriched KEGG pathways
KEGG_top10 <- KEGG_diff2@result[1:10,]
KEGG_top10$pathway <- factor(KEGG_top10$Description, levels = rev(KEGG_top10$Description))
mytheme <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11), plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), legend.title = element_text(size = 13), legend.text = element_text(size = 11))

p <- ggplot(data = KEGG_top10, aes(x = FoldEnrichment, y = pathway)) +
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_color_distiller(palette = "Spectral", direction = -1) +
  labs(x = "Fold Enrichment", y = "", title = "Dotplot of Enriched KEGG Pathways", size = "Gene Number") +
  theme_bw() +
  mytheme

ggsave("test.pdf", p)
```
![相对路径](D:/22222/0课题/output/fig9.png)
