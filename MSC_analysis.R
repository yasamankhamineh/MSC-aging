############################### Visualization of data: GSE68374 ###########################
setwd ()
library (pheatmap)
library (RColorBrewer)

sample_labels <- rep("Fetal MSCs", 3)
sample_labels[4:6] <- "Adult MSCs"
y <- heatmap.raw
y <- y[,-1]
annotation_col <- data.frame(Group = factor(sample_labels, levels = c("Fetal MSCs", "Adult MSCs")))
rownames(annotation_col) <- colnames(y)
any(is.na(y))  # Returns TRUE if there are NAs  
sum(is.na(y))  # Counts the number of NAs  
y <- as.matrix(y)  
if (any(is.na(y))) {  
  cat("Found NA values. Removing rows with NAs.\n")  
  y <- y[complete.cases(y), ]  # Remove rows with any NAs  
}  
if (any(is.nan(y))) {  
  cat("Found NaN values. You may want to handle them.\n")  
  y[is.nan(y)] <- 0  
}  
if (any(is.infinite(y))) {  
  cat("Found Inf values. You may want to handle them.\n")  
  # Example: Replacing Infs with max value  
  y[y == Inf] <- max(y[is.finite(y)], na.rm = TRUE)  
}  
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)  
pheatmap(y,   
         annotation_col = annotation_col,  
         cluster_rows = T,   
         cluster_cols = F,  
         show_rownames = TRUE,   
         show_colnames = TRUE,  
         color = my_colors,      
         scale = "row",          
         fontsize = 10,       
         cellwidth = 30,         
         cellheight = 0.075, 
         fontsize_col = 5)         

pdf("heatmap_sample_groups.pdf", width = 8, height = 6)  
###################################################################### GSEA


Setwd()

library(GEOquery)
library(data.table)
library(dplyr)
library(tidyr)
library(biomaRt)

# Download the dataset
gse <- getGEO("GSE68374", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])
pheno <- pData(gse[[1]])
# Load platform annotation
annot<- annot2
annot <- annot %>%
  dplyr::mutate(GB_LIST = gsub("[,].*", "", GB_LIST)) %>%
  dplyr::filter(GB_LIST != "") %>%
  dplyr::select(ID, GB_LIST)
annot <- annot %>%
  filter(grepl("^NM_", GB_LIST)) # Keep only NM_ RefSeq IDs
annot <- annot %>% filter(GB_LIST != "" & !is.na(GB_LIST))
annot_long <- annot %>%
  separate_rows(GB_LIST, sep = ",") %>%
  mutate(GB_LIST = trimws(GB_LIST))  # Remove extra spaces
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "description")
mart_data <- getBM(attributes = attributes, mart = ensembl)


write.table(mart_data, "mart_data.txt", sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)

mart_data <- mart_data %>%
  mutate(refseq_mrna = gsub("\\..*$", "", refseq_mrna))
merged_data <- annot_long %>%
  left_join(mart_data, by = c("GB_LIST" = "refseq_mrna"))
chip_data <- merged_data %>%
  dplyr::select(ID, external_gene_name, description) %>%
  distinct(ID, .keep_all = TRUE) %>%
  rename(`Probe Set ID` = ID,
         `GeneSymbol` = external_gene_name,
         `Gene Title` = description)

write.table(chip_data, file ="chip_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


joined_data <- left_join(results_with_names, chip_data, by = "ID")
filtered_results <- joined_data %>%
  filter(!is.na(GeneSymbol.y))
ranked_genes <- filtered_results %>%
  group_by(GeneSymbol.x) %>%
  summarise(stat = max(t)) %>%
  arrange(desc(stat))

ranked_genes_clean <- ranked_genes %>%
  mutate(GeneSymbol.x = trimws(gsub('"', '', GeneSymbol.x))) %>%  # remove quotes and spaces
  filter(!is.na(stat)) %>%                                        # remove NA stats
  distinct(GeneSymbol.x, .keep_all = TRUE) %>%                    # remove duplicates
  arrange(desc(stat))

ranks <- setNames(ranked_genes_clean$stat, ranked_genes_clean$GeneSymbol.x)

names(ranks) <- trimws(names(ranks))


library(msigdbr)
library(dplyr)
library(tidyr)
library(fgsea)
all_sets <- msigdbr(species = "Homo sapiens")
kegg_sets_df <- all_sets %>%
  filter(gs_cat == "C2" & grepl("KEGG", gs_name))
hallmark_sets_df <- all_sets %>%
  filter(gs_cat == "H")

kegg_sets <- split(kegg_sets_df$gene_symbol, kegg_sets_df$gs_name)
hallmark_sets <- split(hallmark_sets_df$gene_symbol, hallmark_sets_df$gs_name)

set.seed(123)
fgsea_kegg <- fgseaMultilevel(pathways = kegg_sets, stats = ranks)

set.seed(123)
fgsea_hallmark <- fgseaMultilevel(pathways = hallmark_sets, stats = ranks)

fgsea_kegg <- fgseaMultilevel(pathways = kegg_sets, stats = ranks, minSize=15, maxSize=500)
fgsea_hallmark <- fgseaMultilevel(pathways = hallmark_sets, stats = ranks, minSize=15, maxSize=500)

top_kegg <- fgsea_kegg %>% arrange(padj) %>% head(50)
top_hallmark <- fgsea_hallmark %>% arrange(padj) %>% head(50)


library(fgsea)
library(ggplot2)
library(cowplot)
library(gridExtra)  # for arranging multiple ggplots
library(dplyr)
library(magrittr)
library(enrichplot)
n_plot <- 50
top50_kegg <- top_kegg %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)
top50_hallmark <- top_hallmark %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)

generate_enrichment_plots <- function(pathways, pathway_list, stats = ranks, prefix) {
  plot_list <- vector("list", length(pathway_list))
  for (i in seq_along(pathway_list)) {
    p <- plotEnrichment(pathways[[pathway_list[i]]], stats = stats) +
      ggtitle(paste0(prefix, ": ", pathway_list[i])) +
      xlab("Rank") +
      ylab("Enrichment score") +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 7,          # Larger title
          face = "bold",
          vjust = 2           # Move title up
        ),
        axis.title.x = element_text(
          size = 7,
          margin = margin(t = 12)  # Space above x-axis title
        ),
        axis.title.y = element_text(
          size = 7,
          margin = margin(r = 8)   # Space to the right of y-axis title
        ),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin = margin(t = 18, r = 10, b = 18, l = 10) # Top/Bottom/Right/Left
      )
    plot_list[[i]] <- p
  }
  return(plot_list)
}

kegg_plots <- generate_enrichment_plots(kegg_sets, top50_kegg, stats = ranks, "Enrichment plot")
hallmark_plots <- generate_enrichment_plots(hallmark_sets, top50_hallmark, stats = ranks, "Enrichment plot")

grid.arrange(grobs = kegg_plots[1:4], ncol = 2)

for (i in seq_along(kegg_plots)) {
  ggsave(filename = paste0("Enrichment plot", i, ".png"), plot = kegg_plots[[i]], 
         width = 4, height = 3)
}

for (i in seq_along(hallmark_plots)) {
  ggsave(filename = paste0("Hallmark_enrichment_", i, ".png"), plot = hallmark_plots[[i]], 
         width = 4, height = 3)
}


############################ For Validation: GSE97311###############################
setwd()

library(GEOquery)  
library(limma)  

data = getGEO(GEO = "GSE97311", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)  
eset = data[["GSE97311_series_matrix.txt.gz"]]  
norm = eset@assayData[["exprs"]]  
boxplot(norm)  
norm = normalizeQuantiles(norm)  
exprs(eset) = norm  
norm <- norm[, -c(2, 3, 5, 7:14)]  # Remove unwanted columns  
# Log2 transformation  
log2norm = log2(norm)  
groups = c("Aged_MSCs", "Aged_MSCs", "Fetal_MSCs", "Aged_MSCs", "Aged_MSCs", "Fetal_MSCs", "Fetal_MSCs")  
groups = factor(groups, levels = unique(groups))  
design = model.matrix(object = ~ 0 + groups)  
colnames(design) = levels(groups)  
fit = lmFit(object = log2norm, design = design)  
cm = makeContrasts(contrasts = c("Aged_MSCs - Fetal_MSCs"), levels = design)  
fit2 = contrasts.fit(fit = fit, contrasts = cm)  
fit3 = eBayes(fit = fit2)  
df = data.frame(fit3)  
results=topTable(fit=fit3,coef = "Aged_MSCs - Fetal_MSCs",
                 number = Inf,
                 adjust.method = "fdr",
                 sort.by = "logFC")
# Filter for DEGs according to your criteria |log2(FC)| > 1 and P < 0.05
final_deg = results[abs(results$logFC) > 1 & results$P.Value < 0.05,]

top_genes <- final_deg %>% group_by(GeneSymbol) %>% 
  filter(logFC == max(logFC, na.rm = TRUE)) %>%
  ungroup()

write.table(results, file="results.txt", quote = F, sep="\t")
#########################################################
results$IlluminaID <- rownames(results)  
results <- results[, c("IlluminaID", setdiff(names(results), "IlluminaID"))]  
rownames(results) <- NULL  
print(str(results$IlluminaID))  
print(str(illumina$IlluminaID))  
results$IlluminaID <- as.character(results$IlluminaID)  
illumina$IlluminaID <- as.character(illumina$IlluminaID) 
final= merge(results, illumina , by= "IlluminaID", all = T)

write.table(final, file="final.txt", quote = F, sep="\t")

# Filter for DEGs according to your criteria |log2(FC)| > 1 and P < 0.05  
final_deg = final[abs(final$logFC) > 1 & final$P.Value < 0.05,] 

write.table(final_deg, file="final_deg.txt", quote = F, sep="\t")
####################################### OR find out the DEGs separately:1. merge
final2= merge(DEGs, GeneSymbol, by= "GeneSymbol", all = T)

write.table(final3, file = "final3.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

######################################## :2. DEGs
final_deg = final2[abs(final$logFC) > 1 & final$P.Value < 0.05,] 

write.table(final_deg, file= "final_deg.txt", quote = FALSE, sep = "\t") 

################################################################# Create Heatmap
library(pheatmap)
library(RColorBrewer)
fetal_labels <- rep("Fetal MSCs", 3)  
gsm_labels <- c("GSM2561504", "GSM2561515", "GSM2561516")  
sample_labels <- c(fetal_labels, gsm_labels)  
aged_labels <- rep("Aged MSCs", 4)  
gsm_labels_aged <- c("GSM2561499", "GSM2561502", "GSM2561513", "GSM2561514")  
sample_labels_aged <- c(aged_labels, gsm_labels_aged)  
y <- heatmap97311  
y <- y[,-1]  # Remove the first column  
annotation_col <- data.frame((Group = c(rep("Fetal MSCs", 3), rep("Aged MSCs", 4)))  
                             rownames(annotation_col) <- c("GSM2561504", "GSM2561515", "GSM2561516", "GSM2561499", "GSM2561502", "GSM2561513", "GSM2561514")  
                             if (ncol(y) != nrow(annotation_col)) {  
                               stop("The number of columns in y must match the number of rows in annotation_col.")  
                             }  
                             my_colors <- colorRampPalette(c("blue", "white", "red"))(100) 
                             pheatmap(y,  annotation_col = annotation_col,   cluster_rows = TRUE,   
                                      cluster_cols = T,  
                                      show_rownames = F,   
                                      show_colnames = TRUE,  
                                      color = my_colors,      
                                      scale = "row",          
                                      fontsize = 10,          
                                      cellwidth = 20,         
                                      cellheight = 0.2,   
                                      fontsize_col = 4)    
                             
                             ################################### GSEA
                             
                             setwd()
                             
                             library(GEOquery)
                             library(data.table)
                             library(dplyr)
                             library(tidyr)
                             library(biomaRt)
                             
                             gse <- getGEO("GSE97311", GSEMatrix = TRUE)
                             exprSet <- exprs(gse[[1]])
                             pheno <- pData(gse[[1]])
                             
                             # Load platform annotation
                             annot2<- annot
                             
                             # Extract relevant columns and clean up
                             annot2 <- annot2 [, -2:-7]
                             annot2 <- annot2 [, -3:-24]
                             annot2 <- annot2 %>%
                               filter(grepl("^NM_", RefSeq_ID)) # Keep only NM_ RefSeq IDs
                             
                             
                             
                             annot_long <- annot2 %>%
                               separate_rows(RefSeq_ID, sep = ",") %>%
                               mutate(RefSeq_ID = trimws(RefSeq_ID))  # Remove extra spaces
                             
                             annot_long$RefSeq_ID <- sub("\\..*$", "", annot_long$RefSeq_ID)
                             
                             annot_long<- annot_long[,-2]
                             
                             ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
                             
                             attributes <- c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "description")
                             mart_data <- getBM(attributes = attributes, mart = ensembl)
                             
                             mart_data <- mart_data %>%
                               mutate(refseq_mrna = gsub("\\..*$", "", refseq_mrna))
                             
                             merged_data <- annot_long %>%
                               left_join(mart_data, by = c("RefSeq_ID" = "RefSeq_ID"))
                             
                             chip_data <- merged_data %>%
                               dplyr::select(ID, RefSeq_ID, external_gene_name, description) %>%
                               distinct(ID, .keep_all = TRUE) %>%
                               rename(`Ref_ID` = RefSeq_ID,
                                      `GeneSymbol` = external_gene_name,
                                      `Gene Title` = description)
                             
                             chip_data_clean <- na.omit(chip_data)
                             
                             data = getGEO(GEO = "GSE97311", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)  
                             eset = data[["GSE97311_series_matrix.txt.gz"]]  
                             norm = eset@assayData[["exprs"]]  
                             boxplot(norm)  
                             norm = normalizeQuantiles(norm)  
                             exprs(eset) = norm  
                             norm <- norm[, -c(2, 3, 5, 7:14)]  # Remove unwanted columns  
                             log2norm = log2(norm)  
                             groups = c("Aged_MSCs", "Aged_MSCs", "Fetal_MSCs", "Aged_MSCs", "Aged_MSCs", "Fetal_MSCs", "Fetal_MSCs")  
                             groups = factor(groups, levels = unique(groups))  
                             design = model.matrix(object = ~ 0 + groups)  
                             colnames(design) = levels(groups)  
                             fit = lmFit(object = log2norm, design = design)  
                             cm = makeContrasts(contrasts = c("Aged_MSCs - Fetal_MSCs"), levels = design)  
                             fit2 = contrasts.fit(fit = fit, contrasts = cm)  
                             fit3 = eBayes(fit = fit2)  
                             df = data.frame(fit3)  
                             results_with_names=topTable(fit=fit3,coef = "Aged_MSCs - Fetal_MSCs",
                                                         number = Inf,
                                                         adjust.method = "fdr",
                                                         sort.by = "logFC")
                             ID<- rownames(results_with_names) 
                             results_with_names<- cbind(results_with_names, ID)
                             
                             joined_data <- left_join(ID, chip_data_clean, by = "ID")
                             
                             filtered_results <- joined_data %>%
                               filter(!is.na(GeneSymbol))
                             
                             ranked_genes <- filtered_results %>%
                               group_by(GeneSymbol) %>%
                               summarise(stat = max(t)) %>%
                               arrange(desc(stat))
                             
                             ranks <- setNames(ranked_genes$stat, ranked_genes$`GeneSymbol`)
                             
                             names(ranks) <- trimws(names(ranks))
                             
                             library(msigdbr)
                             library(dplyr)
                             library(tidyr)
                             all_sets <- msigdbr(species = "Homo sapiens")
                             
                             kegg_sets_df <- all_sets %>%
                               filter(gs_cat == "C2" & grepl("KEGG", gs_name))
                             
                             kegg_list <- kegg_sets_df %>%
                               group_by(gs_name) %>%
                               summarise(gene_list = list(gene_symbol), .groups = "drop") %>%
                               pull(gene_list, name = gs_name)
                             
                             hallmark_sets_df <- all_sets %>%
                               filter(gs_cat == "H")
                             
                             kegg_sets <- split(kegg_sets_df$gene_symbol, kegg_sets_df$gs_name)
                             hallmark_sets <- split(hallmark_sets_df$gene_symbol, hallmark_sets_df$gs_name)
                             
                             set.seed(123)
                             fgsea_kegg <- fgseaMultilevel(pathways = kegg_sets, stats = ranks)
                             
                             set.seed(123)
                             fgsea_hallmark <- fgseaMultilevel(pathways = hallmark_sets, stats = ranks)
                             
                             top_kegg <- fgsea_kegg %>% arrange(padj) %>% head(50)
                             top_hallmark <- fgsea_hallmark %>% arrange(padj) %>% head(50)
                             
                             library(fgsea)
                             library(ggplot2)
                             library(cowplot)
                             library(gridExtra)  # for arranging multiple ggplots
                             library(dplyr)
                             library(magrittr)
                             n_plot <- 50
                             top50_kegg <- top_kegg %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)
                             top50_hallmark <- top_hallmark %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)
                             
                             generate_enrichment_plots <- function(pathways, pathway_list, stats = ranks, prefix) {
                               plot_list <- vector("list", length(pathway_list))
                               for (i in seq_along(pathway_list)) {
                                 p <- plotEnrichment(pathways[[pathway_list[i]]], stats = stats) +
                                   ggtitle(paste0(prefix, ": ", pathway_list[i])) +
                                   xlab("Rank") +
                                   ylab("Enrichment score") +
                                   theme_minimal() +
                                   theme(
                                     plot.title = element_text(
                                       hjust = 0.5,
                                       size = 7,          # Larger title
                                       face = "bold",
                                       vjust = 2           # Move title up
                                     ),
                                     axis.title.x = element_text(
                                       size = 7,
                                       margin = margin(t = 12)  # Space above x-axis title
                                     ),
                                     axis.title.y = element_text(
                                       size = 7,
                                       margin = margin(r = 8)   # Space to the right of y-axis title
                                     ),
                                     axis.text.x = element_text(size = 6),
                                     axis.text.y = element_text(size = 6),
                                     plot.margin = margin(t = 18, r = 10, b = 18, l = 10) # Top/Bottom/Right/Left
                                   )
                                 plot_list[[i]] <- p
                               }
                               return(plot_list)
                             }
                             
                             kegg_plots <- generate_enrichment_plots(kegg_sets, top50_kegg, stats = ranks, "Enrichment plot")
                             hallmark_plots <- generate_enrichment_plots(hallmark_sets, top50_hallmark, stats = ranks, "Enrichment plot")
                             grid.arrange(grobs = kegg_plots[1:4], ncol = 2)
                             
                             # To save all KEGG plots as separate files:
                             for (i in seq_along(kegg_plots)) {
                               ggsave(filename = paste0("Enrichment plot", i, ".png"), plot = kegg_plots[[i]], 
                                      width = 4, height = 3)
                             }
                             
                             # Similarly, save Hallmark plots if desired
                             for (i in seq_along(hallmark_plots)) {
                               ggsave(filename = paste0("Hallmark_enrichment_", i, ".png"), plot = hallmark_plots[[i]], 
                                      width = 4, height = 3)
                             }
                             
                             ############################ For Validation: GSE119987###############################
                             setwd()
                             library("Biobase")
                             library("GEOquery")
                             library("Biobase")
                             library("GEOquery")
                             library("limma") 
                             library(dplyr)
                             
                             gse <- getGEO("GSE119987", GSEMatrix = TRUE)
                             eset <- gse[[1]]
                             selected_samples <- c("GSM3389898","GSM3389922","GSM3389899",
                                                   "GSM3389910","GSM3389923","GSM3389944", "GSM3389943","GSM3389942","GSM3389941", "GSM3389921","GSM3389940")
                             eset <- eset[, selected_samples]
                             group <- factor(c(rep("Early", 6), rep("Late", 5)), levels = c("Early", "Late"))
                             design   <- model.matrix(~0 + group)
                             colnames(design) <- levels(group)
                             contrast <- makeContrasts(Late - Early, levels = design)
                             fit      <- lmFit(eset, design)
                             fit2     <- contrasts.fit(fit, contrast) |> eBayes()
                             results  <- topTable(fit2, coef = "Late - Early", number = Inf, adjust.method = "fdr", sort.by = "logFC")
                             get_symbol_map <- function(eset_obj, probes) {
                               fd      <- fData(eset_obj)
                               sym_col <- grep("symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
                               if (!is.na(sym_col)) {
                                 syms <- fd[probes, sym_col]
                                 return(setNames(syms, probes))}
                               ann_pkg <- paste0(annotation(eset_obj), ".db")
                               if (!requireNamespace(ann_pkg, quietly = TRUE)) {
                                 BiocManager::install(ann_pkg, ask = FALSE)  }
                               library(ann_pkg, character.only = TRUE)
                               mapping <- AnnotationDbi::select(get(ann_pkg),
                                                                keys    = probes,
                                                                keytype = "PROBEID",
                                                                columns = "SYMBOL")
                               mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
                               setNames(mapping$SYMBOL, mapping$PROBEID)
                             }
                             symbol_map           <- get_symbol_map(eset, rownames(results))
                             results$GeneSymbol  <- symbol_map[rownames(results)]
                             
                             # Filter DEGs by logFC and raw p-value
                             deg <- results %>% filter(abs(logFC) > 1, P.Value < 0.05, !is.na(GeneSymbol))
                             top_genes <- deg %>% group_by(GeneSymbol) %>% 
                               filter(logFC == max(logFC, na.rm = TRUE)) %>%
                               ungroup()
                             
                             write.table(top_genes, file=" top_genes.txt", quote = F, sep="\t")
                             ################################################################ Create Heatmap
                             setwd ()
                             library(pheatmap)
                             library(RColorBrewer)
                             
                             heatmap_data <- exprs_data[rownames(deg), ]
                             ann_col      <- data.frame(Group = group, row.names = colnames(heatmap_data))
                             pheatmap(heatmap_data, annotation_col = ann_col,
                                      scale = "row", color = colorRampPalette(c("blue", "white", "red"))(100),
                                      cluster_rows = TRUE, cluster_cols = TRUE,  show_rownames = FALSE,  fontsize_col = 8)
                             
                             ######################################################################### GSEA
                             Setwd()
                             
                             library(GEOquery)
                             library(data.table)
                             library(dplyr)
                             library(tidyr)
                             library(biomaRt)
                             gse <- getGEO("GSE119987", GSEMatrix = TRUE)
                             exprSet <- exprs(gse[[1]])
                             pheno <- pData(gse[[1]])
                             
                             annot2<- annot
                             annot2 <- annot2 [, -2:-12]
                             annot2 <- annot2 [, -3:-5]
                             annot2 <- annot2 %>%
                               filter(grepl("^NM_", RefSeq_ID)) # Keep only NM_ RefSeq IDs
                             annot_long <- annot2 %>%
                               separate_rows(RefSeq_ID, sep = ",") %>%
                               mutate(RefSeq_ID = trimws(RefSeq_ID))  # Remove extra spaces
                             
                             annot_long$RefSeq_ID <- sapply(strsplit(annot_long$RefSeq_ID, " /// "), `[`, 1)
                             ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
                             attributes <- c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "description")
                             mart_data <- getBM(attributes = attributes, mart = ensembl)
                             
                             mart_data <- mart_data %>%
                               mutate(refseq_mrna = gsub("\\..*$", "", refseq_mrna))
                             
                             merged_data <- annot_long %>%
                               left_join(mart_data, by = c("RefSeq_ID" = "refseq_mrna"))
                             chip_data <- merged_data %>%
                               dplyr::select(ID, RefSeq_ID, external_gene_name, description) %>%
                               distinct(ID, .keep_all = TRUE) %>%
                               rename(`Ref_ID` = RefSeq_ID,
                                      `GeneSymbol` = external_gene_name,
                                      `Gene Title` = description)
                             
                             chip_data_clean <- na.omit(chip_data)
                             
                             # get the total genes
                             library(GEOquery)  
                             library(limma)  
                             
                             data = getGEO(GEO = "GSE119987", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)  
                             eset = data[["GSE119987_series_matrix.txt.gz"]]  
                             norm = eset@assayData[["exprs"]]  
                             boxplot(norm)  
                             norm = normalizeQuantiles(norm)  
                             exprs(eset) = norm
                             norm <- norm[, -c(3:12, 14:23, 27:42, 48:50)]  # Remove unwanted columns
                             norm <- norm[, -12:-18]  # Remove unwanted columns
                             
                             groups = c("Aged_MSCs", "Aged_MSCs", "Aged_MSCs", "Young_MSCs", "Aged_MSCs", "Aged_MSCs", 
                                        "Young_MSCs","Young_MSCs","Young_MSCs","Young_MSCs","Aged_MSCs") 
                             
                             groups = factor(groups, levels = unique(groups))  
                             design = model.matrix(object = ~ 0 + groups)  
                             colnames(design) = levels(groups)  
                             fit = lmFit(object = norm, design = design)  
                             cm = makeContrasts(contrasts = c("Aged_MSCs - Young_MSCs"), levels = design)  
                             fit2 = contrasts.fit(fit = fit, contrasts = cm)  
                             fit3 = eBayes(fit = fit2)  
                             df = data.frame(fit3)  
                             results_with_names=topTable(fit=fit3,coef = "Aged_MSCs - Young_MSCs",
                                                         number = Inf,
                                                         adjust.method = "fdr",
                                                         sort.by = "logFC")
                             ID<- rownames(results_with_names) 
                             results_with_names<- cbind(results_with_names, ID)
                             
                             joined_data <- left_join(results_with_names, chip_data_clean, by = "ID")
                             
                             filtered_results <- joined_data %>%
                               filter(!is.na(GeneSymbol))
                             
                             ranked_genes <- filtered_results %>%
                               group_by(GeneSymbol) %>%
                               summarise(stat = max(t)) %>%
                               arrange(desc(stat))
                             ranks <- setNames(ranked_genes$stat, ranked_genes$`GeneSymbol`)
                             
                             library(msigdbr)
                             library(dplyr)
                             library(tidyr)
                             library(fgsea)
                             library(ggplot2)
                             library(cowplot)
                             library(gridExtra)  # for arranging multiple ggplots
                             library(magrittr)
                             all_sets <- msigdbr(species = "Homo sapiens")
                             
                             kegg_sets_df <- all_sets %>%
                               filter(gs_cat == "C2" & grepl("KEGG", gs_name))
                             kegg_list <- kegg_sets_df %>%
                               group_by(gs_name) %>%
                               summarise(gene_list = list(gene_symbol), .groups = "drop") %>%
                               pull(gene_list, name = gs_name)
                             hallmark_sets_df <- all_sets %>%
                               filter(gs_cat == "H")
                             kegg_sets <- split(kegg_sets_df$gene_symbol, kegg_sets_df$gs_name)
                             hallmark_sets <- split(hallmark_sets_df$gene_symbol, hallmark_sets_df$gs_name)
                             
                             names(ranks) <- trimws(names(ranks))
                             set.seed(123)
                             fgsea_kegg <- fgseaMultilevel(pathways = kegg_sets, stats = ranks)
                             
                             set.seed(123)
                             fgsea_hallmark <- fgseaMultilevel(pathways = hallmark_sets, stats = ranks)
                             
                             top_kegg <- fgsea_kegg %>% arrange(padj) %>% head(50)
                             top_hallmark <- fgsea_hallmark %>% arrange(padj) %>% head(50)
                             n_plot <- 50
                             top50_kegg <- top_kegg %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)
                             top50_hallmark <- top_hallmark %>% arrange(padj) %>% slice(1:n_plot) %>% pull(pathway)
                             
                             
                             generate_enrichment_plots <- function(pathways, pathway_list, stats = ranks, prefix) {
                               plot_list <- vector("list", length(pathway_list))
                               for (i in seq_along(pathway_list)) {
                                 p <- plotEnrichment(pathways[[pathway_list[i]]], stats = stats) +
                                   ggtitle(paste0(prefix, ": ", pathway_list[i])) +
                                   xlab("Rank") +
                                   ylab("Enrichment score") +
                                   theme_minimal() +
                                   theme(
                                     plot.title = element_text(
                                       hjust = 0.5,
                                       size = 7,          # Larger title
                                       face = "bold",
                                       vjust = 2           # Move title up
                                     ),
                                     axis.title.x = element_text(
                                       size = 7,
                                       margin = margin(t = 12)  # Space above x-axis title
                                     ),
                                     axis.title.y = element_text(
                                       size = 7,
                                       margin = margin(r = 8)   # Space to the right of y-axis title
                                     ),
                                     axis.text.x = element_text(size = 6),
                                     axis.text.y = element_text(size = 6),
                                     plot.margin = margin(t = 18, r = 10, b = 18, l = 10) # Top/Bottom/Right/Left
                                   )
                                 plot_list[[i]] <- p
                               }
                               return(plot_list)
                             }
                             kegg_plots <- generate_enrichment_plots(kegg_sets, top50_kegg, stats = ranks, "Enrichment plot")
                             hallmark_plots <- generate_enrichment_plots(hallmark_sets, top50_hallmark, stats = ranks, "Enrichment plot")
                             grid.arrange(grobs = kegg_plots[1:4], ncol = 2)
                             
                             for (i in seq_along(kegg_plots)) {
                               ggsave(filename = paste0("Enrichment plot", i, ".png"), plot = kegg_plots[[i]], 
                                      width = 4, height = 3)
                             }
                             
                             # Similarly, save Hallmark plots if desired
                             for (i in seq_along(hallmark_plots)) {
                               ggsave(filename = paste0("Hallmark_enrichment_", i, ".png"), plot = hallmark_plots[[i]], 
                                      width = 4, height = 3)
                             }
                             
                             