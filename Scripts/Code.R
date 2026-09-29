library(ape)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggtext)
library(pheatmap)
library(phyloseq)
library(RColorBrewer)
library(readr)
library(reshape2)
library(rstatix)
library(scales)
library(tibble)
library(tidyr)

#---- Fig 1----
SNPS<- read.csv("/Users/irene/OneDrive/TESIS/3Cap/ARBOL/tip_SNP_counts.SNPs_all.parsimony", header =T, sep= "\t")
SNPs_filter <- SNPS %>%
  filter(name_on_tree %in% genomas$Id)
SNPs_filter$combined <- paste(SNPs_filter$name_on_tree, SNPs_filter$SNP_counts, sep = "_")

metadata<- read.csv("metadata_tesis_gephi.csv", header =T, sep= ";")
metadata <- metadata %>%
  left_join(SNPs_filter %>% select(name_on_tree, combined), 
            by = c("Id" = "name_on_tree"))
write.csv(metadata, "/Users/irene/OneDrive/TESIS/3Cap/ARBOL/metadata_parsimony.csv", row.names = FALSE)



#---- GENOMAD Fig 1 & Table 2---- 
# Leer un archivo TSV
datap <- read.table("genomad_ALL_plasmid_summary.tsv", header = TRUE, sep = "\t")
genomad_plasmid <- datap %>%
  filter(folder_name %in% filtered_genomes$Id)
write.csv(genomad_plasmid, "genomad_plasmid_tesis.csv", row.names = FALSE)
datav <- read.table("genomad_ALL_virus_summary.tsv", header = TRUE, sep = "\t")
genomad_virus <- datav %>%
  filter(folder_name %in% filtered_genomes$Id)
write.csv(genomad_virus, "genomad_virus_tesis.csv", row.names = FALSE)

datavi <- read.csv("genomad_virus_humanizado.csv", sep= ",")
datapla <- read.csv("genomad_plasmid_humanizado.csv", sep= ",")
genomad_plasmid_huma <- datapla %>%
  filter(folder_name %in% filtered_genomes$Id)
write.csv(genomad_plasmid_huma, "genomad_plasmid_tesis_humanizado.csv.csv", row.names = FALSE)
datav <- read.table("genomad_ALL_virus_summary.tsv", header = TRUE, sep = "\t")
genomad_virus_huma <- datavi %>%
  filter(folder_name %in% filtered_genomes$Id)
write.csv(genomad_virus_huma, "genomad_virus_tesis_humanizado.csv", row.names = FALSE)

#---- Extended Figure 2d----
busco <- read.table(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/busco_summary_ALL.csv",
  header = TRUE,
  sep = ";",
  stringsAsFactors = FALSE
)
busco$missing_fraction <- busco$Missing.BUSCOs..M./
  busco$total.BUSCO.groups.searched

p_missing <- ggplot(busco, aes(x = missing_fraction)) +
  geom_histogram(
    bins = 25,
    fill = "#7BAF9E",
    color = "white",
    linewidth = 0.3
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Fraction of missing BUSCO genes",
    y = "Genomes"
  )

p_missing
#---- Extended Figure 1----
setwd("N:/PAPERS/NatMicroPaggPan/ppts_Figuras")

tabla <- read.table(text = "
Especie P.agglomerans P.vagans P.pleuroti P.annatis P.alli P.beijingensis P.rodasii
Aegilops_tauschii 58.33333333 25 0 8.333333333 0 4.166666667 4.166666667
Triticum_monococcum 70.83333333 25 0 0 4.166666667 0 0
T.turgidum 50 33.33333333 0 0 16.66666667 0 0
T.dicoccum 78.26086957 21.73913043 0 0 0 0 0
T.durum 60 36 0 0 4 0 0
T.spelta 70 20 10 0 0 0 0
T.aestivum 83.33333333 8.333333333 8.333333333 0 0 0 0
", header = TRUE)

# Convertir a formato largo
tabla_long <- melt(tabla, id.vars = "Especie", variable.name = "Pantoea", value.name = "Porcentaje")

# Ordenar especies como en tu tabla (izq → der)
orden_especies <- c("Aegilops_tauschii", "Triticum_monococcum", "T.turgidum", 
                    "T.dicoccum", "T.durum", "T.aestivum", "T.spelta")
tabla_long$Especie <- factor(tabla_long$Especie, levels = orden_especies)

# Ordenar Pantoea por abundancia global (P.agglomerans primero)
pantoea_orden <- tabla %>%
  select(-Especie) %>%
  colSums() %>%
  sort(decreasing = TRUE) %>%
  names()
tabla_long$Pantoea <- factor(tabla_long$Pantoea, levels = rev(pantoea_orden))

# Colores como en tu figura
colores <- c("P.agglomerans" = "#F8766D",    # Rojo/naranja
             "P.vagans" = "#00BA38",        # Verde
             "P.pleuroti" = "#619CFF",      # Azul
             "P.annatis" = "#F564E3",       # Rosa
             "P.alli" = "#00BFC4",          # Cyan
             "P.beijingensis" = "#FECA57",  # Amarillo
             "P.rodasii" = "#A58CFF")       # Violeta

# Crear el plot
p <- ggplot(tabla_long, aes(x = Especie, y = Porcentaje, fill = Pantoea)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colores) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), 
                     breaks = seq(0, 100, 25)) +
  labs(x = "Wheat species", y = "Percentage (%)", 
       title = "Pantoea species composition") +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 14),  # Cursiva en eje Y
    axis.title.y = element_markdown(size = 16),
    axis.text.x = element_text(size = 14),
    legend.title = element_markdown(size = 14),  # Cursiva en leyenda
    legend.text = element_markdown(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = guide_legend(reverse = FALSE))

p

#---- Extended Table 6----
roary_matrix = read_tsv("gene_presence_absence.Rtab", col_names = TRUE)
grupo <- read_delim('ji_grupo.csv', col_names = TRUE, delim=";")
grupo <- tibble(Id = 'ISP212TRI626322bis', ji_group = 3) %>% full_join(grupo)

ji_grupo1 <- grupo %>%
  filter(ji_group == 1) %>%  # Filtra para el grupo 1
  pull(Id)

ji_grupo2 <- grupo %>%
  filter(ji_group == 2) %>%  # Filtra para el grupo 1
  pull(Id)

ji_grupo3 <- grupo %>%
  filter(ji_group == 3) %>%  # Filtra para el grupo 1
  pull(Id)

C1_matrix <- roary_matrix[, c("Gene", ji_grupo1)]
C2_matrix <- roary_matrix[, c("Gene", ji_grupo2)]
C3_matrix <- roary_matrix[, c("Gene", ji_grupo3)]

C1_matrix <- as.data.frame(C1_matrix)
C2_matrix <- as.data.frame(C2_matrix)
C3_matrix <- as.data.frame(C3_matrix)

rownames(C1_matrix) <- C1_matrix$Gene
C1_matrix <- C1_matrix[, -1]
rownames(C2_matrix) <- C2_matrix$Gene
C2_matrix <- C2_matrix[, -1]
rownames(C3_matrix) <- C3_matrix$Gene
C3_matrix <- C3_matrix[, -1]

rows_to_keep_C1 <- rowSums(C1_matrix != 0) > 0
rows_to_keep_C2 <- rowSums(C2_matrix != 0) > 0
rows_to_keep_C3 <- rowSums(C3_matrix != 0) > 0

final_C1 <- C1_matrix[rows_to_keep_C1, ]
final_C2 <- C2_matrix[rows_to_keep_C2, ]
final_C3 <- C3_matrix[rows_to_keep_C3, ]

notation_matrix = read_csv("gene_presence_absence.csv", col_names = TRUE) %>% select('Gene', 'Annotation')

genes_to_filter_C1_core <- row.names(C1_final_core)
C1_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C1_core, ]
genes_to_filter_C2_core <- row.names(C2_final_core)
C2_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C2_core, ]
genes_to_filter_C3_core <- row.names(C3_final_core)
C3_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C3_core, ]

final_rows_to_keep_C1_core99 <- rowSums(final_C1 != 0) >= 0.99*55  # Reemplaza el 0.99 por el threshold para core que quieras aplicar
final_rows_to_keep_C2_core99 <- rowSums(final_C2 != 0) >= 0.99*252  
final_rows_to_keep_C3_core99 <- rowSums(final_C3 != 0) >= 0.99*41  

C1_final_core99 <- final_C1[final_rows_to_keep_C1_core99, ]
C2_final_core99 <- final_C2[final_rows_to_keep_C2_core99, ]
C3_final_core99 <- final_C3[final_rows_to_keep_C3_core99, ]

final_rows_to_keep_C1_any <- rowSums(final_C1 != 0) >= 1  # Reemplaza el 0 por el numero de genomas del grupo que quieras analizar
final_rows_to_keep_C2_any <- rowSums(final_C2 != 0) >= 1  
final_rows_to_keep_C3_any <- rowSums(final_C3 != 0) >= 1  

# Filtrar la matriz para eliminar filas que tienen solo ceros
C1_final_any <- final_C1[final_rows_to_keep_C1_any, ]
C2_final_any <- final_C2[final_rows_to_keep_C2_any, ]
C3_final_any <- final_C3[final_rows_to_keep_C3_any, ]

#Ahora seleciono los genes de cada cluster y filtro la matriz de notation con esa lista de genes
genes_to_filter_C1_any <- row.names(C1_final_any)
C1_notation_any <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C1_any, ]
genes_to_filter_C2_any <- row.names(C2_final_any)
C2_notation_any <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C2_any, ]
genes_to_filter_C3_any <- row.names(C3_final_any)
C3_notation_any <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C3_any, ]

rownames(C3_final_core)[which(! rownames(C3_final_core) %in% union(rownames(C1_final_any), rownames(C2_final_any)))]

# Eliminar las columnas de esos genomas
C3.2_final_core <- C3_final_core %>%
  select(-GCF_037152265.1, -GCF_029623495.1)
#ver si tiene cosas especificas q no tienen los clusters 2 y 1
rownames(C3.2_final_core)[which(! rownames(C3_final_core) %in% union(rownames(C1_final_any), rownames(C2_final_any)))]

write.csv(C3_notation_any, file = "C3_ALL_genes_99.csv", row.names = FALSE)

#---- Extended Figure 2b----
roary <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/Pangenoma_NO_split_prokka/gene_presence_absence.csv", check.names = FALSE)

mat <- roary[,15:ncol(roary)]

mat_bin <- ifelse(mat == "", 0, 1)
mat_bin <- as.matrix(mat_bin)
storage.mode(mat_bin) <- "numeric"

cat("Matrix size:", dim(mat_bin), "\n")

n <- ncol(mat_bin)
freq <- rowSums(mat_bin) / n

core <- freq == 1
softcore <- freq >= 0.95 & freq < 1
shell <- freq >= 0.15 & freq < 0.95
cloud <- freq < 0.15

category <- ifelse(core,"Core",
                   ifelse(softcore,"Softcore",
                          ifelse(shell,"Shell","Cloud")))

category_order <- factor(category,
                         levels=c("Core","Softcore","Shell","Cloud"))

order_genes <- order(category_order, freq, decreasing=TRUE)

mat_bin <- mat_bin[order_genes, ]
category <- category[order_genes]

mat_bin_small <- mat_bin[seq(1,nrow(mat_bin),by=5), ]
category_small <- category[seq(1,length(category),by=5)]

cat("Reduced matrix size:", dim(mat_bin_small), "\n")

row_anno <- rowAnnotation(
  Category = category_small,
  col=list(Category=c(
    Core="#0E7C7B",
    Softcore="#79C6A3",
    Shell="#55B69A",
    Cloud="#8CC4C4"))
)

pdf("pangenome_heatmap_categories.pdf", width=8, height=10)
col = c(
  "0" = "#F5F5F5",
  "1" = "#4C8577"
)
row_anno + Heatmap(mat_bin_small,
                   name="Gene presence",
                   col=c("0"="grey","1"="#4C8577"),
                   cluster_rows=FALSE,
                   cluster_columns=FALSE,
                   show_row_names=FALSE,
                   show_column_names=FALSE,
                   use_raster=TRUE,
                   heatmap_legend_param=list(
                     at=c(0,1),
                     labels=c("Absence","Presence")
                   ))

dev.off()
#---- Extended Figure 2c----
reference_genome <- "GCF_019048385.1" #primero elegimos el genoma de referencia Para comprobar si el genoma que quieres usar como referencia está en la lista, hacer esto: "grep("GCF_019048385.1", colnames(roary), value=TRUE)"
roary <- read.csv("gene_presence_absence.csv", check.names=FALSE)

mat <- roary[,15:ncol(roary)]

mat_bin <- ifelse(mat=="",0,1)
mat_bin <- as.matrix(mat_bin)
storage.mode(mat_bin) <- "numeric"

n <- ncol(mat_bin)

freq <- rowSums(mat_bin)/n

category <- ifelse(freq==1,"Core",
                   ifelse(freq>=0.95,"Softcore",
                          ifelse(freq>=0.15,"Shell","Cloud")))

ref_presence <- mat_bin[,reference_genome]

df <- data.frame(
  category = category,
  reference = ref_presence
)

summary_df <- df %>%
  group_by(category, reference) %>%
  summarise(n=n()) %>%
  mutate(type = ifelse(reference==1,"Reference","New genes"))

summary_df$category <- factor(summary_df$category, 
                              levels = c("Cloud", "Shell", "Softcore", "Core"))

# Ahora el plot tomará ese orden
p <- ggplot(summary_df,
            aes(x=category,
                y=n,
                fill=type)) +
  
  geom_bar(stat="identity", position="fill") +
  
  scale_fill_manual(values=c(
    "Reference"="darkgreen",
    "New genes"="darkred"
  )) +
  
  ylab("Proportion of gene clusters") +
  xlab("Pangenome category") +
  
  theme_classic()

p

#---- Extended Figure 3 ----
setwd("C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/")

#Read comparisons data from file. Add ANI column based on the Mutational Distance
ji_pag <- read.delim('Filter_ji_filtered.tsv', stringsAsFactors = FALSE) %>%
  mutate(ANI = (1 - Mut_Distance) * 100)
# Add selfcomparisons (not included by bindash)
for (i in union(unique(ji_pag$Source), unique(ji_pag$Target))){
  new_row <- list(Source=i, Target=i, Mut_Distance=0, e_Value=0, JI=1, ANI=100)
  ji_pag <- rbind(ji_pag, new_row, stringsAsFactors = FALSE)
}

# Load JI_Group information. Filter rows with no comparisons
ji_group <- read.delim('PCOA_groups.csv', sep = ';', stringsAsFactors = FALSE) %>%
  filter(Id %in% union(unique(ji_pag$Source), unique(ji_pag$Target))) %>%
  mutate(across(JI_Group, factor))

# Append columns with JI_Group info of reference and query strains to ji_pag dataframe
ji_pag <- left_join(ji_pag, ji_group, by = join_by('Source' == 'Id')) %>%
  rename('Ref.JI_Group' = 'JI_Group')
ji_pag <- left_join(ji_pag, ji_group, by = join_by('Target' == 'Id')) %>% 
  rename('Qry.JI_Group' = 'JI_Group')
# Filter all comparisons with a member of the JI-Group 0
ji_pag <- ji_pag %>%
  filter(Ref.JI_Group != 0) %>%
  filter(Qry.JI_Group != 0)
# Reformat list of comparisons into a matrix
mat <- acast(ji_pag, Source~Target, value.var='ANI')
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    if (is.na(mat[i, j]) == TRUE) {
      mat[i, j] = mat[j, i]
    }
  }
}
# Create dendogram to order heatmap columns and rows
dendr_cls <- hclust(dist(t(mat)), method = 'ward.D')

any(is.na(mat))   # Verifica si hay valores NA
any(is.nan(mat))  # Verifica si hay valores NaN
any(is.infinite(mat))  # Verifica si hay valores Inf
is.numeric(mat) # Asegúrate de que mat sea numérica
mat[is.na(mat)] <- 0
ani_max = 100
ani_thr = 98
ani_low = 97
ani_min = floor(min(97, min(mat, na.rm = TRUE)))
ani_step = 0.1
ani_steps = ((ani_max - ani_min) / ani_step) + 1
ani_breaks = seq(ani_min, ani_max, length.out = ani_steps)
gradient_high = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")[1:4]))(sum(ani_breaks > ani_thr & ani_breaks <= ani_max))
gradient_mid = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")[5:8]))(sum(ani_breaks > ani_low & ani_breaks <= ani_thr))
gradient_bottom = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")[10:11]))(sum(ani_breaks[-1] <= ani_low) + 1)
mat_colors = (c(head(gradient_bottom, -1), gradient_mid, gradient_high))
mat_colors.tiles <- data.frame(xmin = ani_breaks[1:(length(ani_breaks) - 1)],
                               xmax = ani_breaks[2:length(ani_breaks)],
                               ymin = 0,
                               ymax = max(density(ji_pag$ANI, bw = 'SJ')$y),
                               fill = mat_colors,
                               stringsAsFactors = FALSE)
ggplot(mat_colors.tiles) +
  geom_rect(aes(xmin = xmin, xmax = xmax,ymin = ymin, ymax = ymax, fill = fill)) +
  geom_density(data = ji_pag, aes(ANI), bw = 'SJ') +
  scale_fill_identity(aes(fill)) +
  #theme_set(theme_bw(base_size = 16)) +
  theme(legend.position = 'None') +
  ylab('Density')
filename = paste('Pag_ANI_density.png', sep = '')
ggsave(filename)
filename <- 'Pag_ANI_density.pdf'
ggsave(filename, width = 8, height = 6)

filename_ <- 'Pag_ANI_heatmap.pdf'
filename_ = paste('Pag_ANI_heatmap.png', sep = '')
pheatmap(mat,
         cluster_cols = dendr_cls,
         cluster_rows = dendr_cls,
         #treeheight_col = 0,
         #treeheight_row = 0,
         #breaks = mat_breaks,
         #color = mat_colors,
         na_col = '#444444',
         border_color = NA,
         #show_rownames = FALSE,
         angle_col = 315,
         fontsize_col = 3,
         fontsize_row = 3,
         cellwidth = 4,
         cellheight = 4,
         display_numbers = FALSE,
         #number_format = "%.0f",
         fontsize_number = 3,
         #legend = FALSE,
         filename = filename_,
)
#---- Fig2b----
ji_pag %>%
  filter(Source != Target) %>%
  filter(Ref.JI_Group != Qry.JI_Group) %>%
  filter(Ref.JI_Group != 0) %>%
  filter(Qry.JI_Group != 0) %>%
  mutate(Inter = if_else(as.integer(Ref.JI_Group) < as.integer(Qry.JI_Group), paste0(Ref.JI_Group, '<->', Qry.JI_Group), paste0(Qry.JI_Group, '<->', Ref.JI_Group))) %>%
  ggplot(aes(x = Inter, y = ANI)) +
  geom_violin() +
  labs(x = 'Genetic Cluster', y = 'JI Inter-group dispersion')
filename = paste('ANI_inter5group_dispersion.png', sep = '')
ggsave(filename)

ji_pag %>%
  filter(Source != Target) %>%
  filter(Ref.JI_Group != Qry.JI_Group) %>%
  filter(Ref.JI_Group != 4) %>%
  filter(Qry.JI_Group != 4) %>%
  filter(Ref.JI_Group != 2) %>%
  filter(Qry.JI_Group != 2) %>%
  mutate(Inter = if_else(as.integer(Ref.JI_Group) < as.integer(Qry.JI_Group), 
                         paste0(Ref.JI_Group, '<->', Qry.JI_Group), 
                         paste0(Qry.JI_Group, '<->', Ref.JI_Group))) %>%
  ggplot(aes(x = Inter, y = ANI)) +
  geom_violin() +
  scale_x_discrete(labels = c("1<->3" = "C1<->C2", 
                              "1<->5" = "C1<->C3", 
                              "3<->5" = "C2<->C3")) +
  labs(x = 'Genetic Cluster', y = 'JI Inter-group dispersion')+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))  
filename = paste('ANI_inter5group_dispersion.png', sep = '')
#---- Fig2d----
ji_pag %>%
  filter(Source != Target) %>%
  filter(Ref.JI_Group == Qry.JI_Group) %>%
  filter(Ref.JI_Group != 4) %>%
  filter(Ref.JI_Group != 2) %>%
  ggplot(aes(x = Ref.JI_Group, y = JI)) +
  geom_violin() +
  scale_x_discrete(labels = c("1" = "C1", "3" = "C2", "5" = "C3")) +
  labs(x = 'Genetic Cluster', y = 'JI Intra-cluster dispersion')+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
filename = paste('JI_intra5group_dispersion.png', sep = '')
#---- Fig2a----
# Cargar la matriz de similitud (Jaccard Index) en R
jaccard_matrix <- as.matrix(read.csv("jaccard_matrix.csv", row.names=1))
distance_matrix <- 1 - jaccard_matrix  # Convertir matriz de similitud en matriz de distancias

pcoa_result <- pcoa(distance_matrix)

# Visualizar la varianza explicada por los componentes principales
pcoa_result$values$Relative_eig

# Graficar el PCoA con los primeros dos componentes principales
plot(pcoa_result$vectors[,1], pcoa_result$vectors[,2], 
     xlab="PC1", ylab="PC2", main="PCoA basado en Jaccard Index")

# Convertir los resultados de PCoA a un data frame para ggplot2
pcoa_df <- data.frame(PC1 = pcoa_result$vectors[,1], 
                      PC2 = pcoa_result$vectors[,2],
                      Genoma = rownames(pcoa_result$vectors))

# Obtener el porcentaje de varianza explicada por los dos primeros componentes
var_explained <- 100 * pcoa_result$values$Relative_eig[1:2]

# Graficar con puntos en gris oscuro, borde en el gráfico y líneas discontinuas en el origen
p <-ggplot(pcoa_df, aes(x=PC1, y=PC2)) +
  geom_point(size=2, color="#4d4d4d") +  # Puntos en gris oscuro
  labs(
    #title="PCoA basado en Jaccard Index", 
    x=paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y=paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_light() +  # Estilo minimalista
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1), # Borde del gráfico
    #plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 16),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 18)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea horizontal en y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  # Línea vertical en x=0
print(p)
#---- X-square Fig 3f & 4g----
meta <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/gephis/metadatos_gephi2026.csv", sep = ",", stringsAsFactors = FALSE)
clusters <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/Pangenoma_NO_split_prokka/ji_grupo.csv", sep = ";", stringsAsFactors = FALSE)
colnames(clusters)[2] <- "cluster"
data <- merge(meta, clusters, by = "Id")
cluster_counts <- table(data$cluster)
cluster_props <- cluster_counts / sum(cluster_counts)
print(cluster_props)
variables <- c("location", "chi.plant.tissue", "chi.source", "plant_family", "data.source", "sequencing_method","astA", "sitC", "Flagellum", "T1SS","T3SS","T4aP","T4bP","T5aSS","T5bSS","T5cSS","T6SSi", "T4SSi", "T4SSt")
results <- list()

for (var in variables) {
  categories <- unique(data[[var]])
  for (cat in categories) {
    subset_data <- data[data[[var]] == cat, ]
    counts <- table(subset_data$cluster)
    counts_full <- rep(0, length(cluster_props))
    names(counts_full) <- names(cluster_props)
    counts_full[names(counts)] <- counts
    expected <- sum(counts_full) * cluster_props
    std_residuals <- ifelse(expected > 0,
                            (counts_full - expected) / sqrt(expected),
                            0)
    if (sum(counts_full) > 5) {
      test <- chisq.test(counts_full, p = cluster_props, simulate.p.value = TRUE)
      P_val <- test$p.value
    } else {
      P_val <- NA
    }
    results[[paste(var, cat, sep="_")]] <- data.frame(
      Variable = var,
      Category = cat,
      C1_obs = counts_full["1"],
      C2_obs = counts_full["2"],
      C3_obs = counts_full["3"],
      C1_prop = counts_full["1"]/sum(counts_full),
      C2_prop = counts_full["2"]/sum(counts_full),
      C3_prop = counts_full["3"]/sum(counts_full),
      C1_exp = expected["1"],
      C2_exp = expected["2"],
      C3_exp = expected["3"],
      C1_res = std_residuals["1"],
      C2_res = std_residuals["2"],
      C3_res = std_residuals["3"],
      Total = sum(counts_full),
      P_value = P_val
    )
  }
}

results_df <- do.call(rbind, results)
results_df$P_adj_BH <- p.adjust(results_df$P_value, method = "BH")
results_df$P_adj_Bonf <- p.adjust(results_df$P_value, method = "bonferroni")
write.csv(results_df, "results_clusters_all_categories.csv", row.names = FALSE)
#write.csv(meta, "metadatos_gephi2026.csv", row.names = FALSE)

# Preparar datos para heatmap
results_sig <- results_df %>%
  filter(P_adj_Bonf < 0.05)

heatmap_data <- results_sig %>%
  select(Variable, Category, C1_res, C2_res, C3_res)

heatmap_long <- reshape2::melt(
  heatmap_data,
  id.vars = c("Variable", "Category"),
  variable.name = "Cluster",
  value.name = "Residual"
)

heatmap_long$Cluster <- gsub("_res", "", heatmap_long$Cluster)

heatmap_long <- heatmap_long %>%
  mutate(
    VarCat = paste(Variable, Category, sep = " = "),
    VarCat = factor(VarCat, levels = unique(VarCat))
  )

lim <- max(abs(heatmap_long$Residual), na.rm = TRUE)

p <- ggplot(heatmap_long, aes(x = VarCat, y = Cluster, fill = Residual)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = c("#2C3E75", "#CFD1FC", "white", "grey", "#8B1E3F"),
    values  = scales::rescale(c(-lim, -2, 0, 2, lim)),
    limits  = c(-4, 4),
    oob     = scales::squish
  ) +
  labs(x = "Variable / Category", y = "Cluster", fill = "Residual") +
  theme_minimal() +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p
heatmap_long$Category <- factor(
  heatmap_long$Category,
  levels = unique(heatmap_long$Category)
)
heatmap_long <- heatmap_long %>%
  mutate(x = as.numeric(Category))
var_lines <- heatmap_long %>%
  distinct(Variable, Category, x) %>%
  group_by(Variable) %>%
  summarise(x_end = max(x))
p +
  geom_vline(
    data = var_lines,
    aes(xintercept = x_end + 0.5),
    color = "black",
    linewidth = 0.3
  )
#---- Fig5a----
df <- read.table(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/COG__counts.tsv",
  header = FALSE,
  sep = "\t"
)

colnames(df) <- c("COG", "cluster", "count")
wide <- df %>%
  pivot_wider(names_from = cluster,
              values_from = count,
              values_fill = 0)
genomes_per_cluster <- c(
  `1` = 55,   # <-- sustituye por tu valor real
  `2` = 254,
  `3` = 38
)
wide_norm <- wide %>%
  mutate(
    across(c(`1`,`2`,`3`),
           ~ .x / genomes_per_cluster[cur_column()])
  )
COG_fun <- function(x) {
  cog_map <- c(
    "C"="Energy production and conversion",
    "D"="Cell cycle control",
    "E"="Amino acid transport and metabolism",
    "F"="Nucleotide transport and metabolism",
    "G"="Carbohydrate transport and metabolism",
    "H"="Coenzyme transport and metabolism",
    "I"="Lipid transport and metabolism",
    "J"="Translation",
    "K"="Transcription",
    "L"="Replication and repair",
    "M"="Cell wall/membrane/envelope biogenesis",
    "N"="Cell motility",
    "O"="Posttranslational modification",
    "P"="Inorganic ion transport",
    "Q"="Secondary metabolites",
    "T"="Signal transduction",
    "U"="Intracellular trafficking",
    "V"="Defense mechanisms",
    "S"="Function unknown"
  )
  
  cog_map[x]
}

wide_norm$description <- COG_fun(wide_norm$COG)
wide_norm <- wide_norm %>%
  mutate(cluster12 = `1` + `2`)
total12 <- sum(wide$`1` + wide$`2`)
total3  <- sum(wide$`3`)

run_fisher <- function(a, total_a, b, total_b) {
  fisher.test(
    matrix(c(a, total_a - a,
             b, total_b - b),
           nrow = 2)
  )$p.value
}

res <- wide %>%
  mutate(
    cluster12 = `1` + `2`,
    pval = mapply(run_fisher, cluster12, total12, `3`, total3),
    FDR  = p.adjust(pval, method = "fdr"),
    logFDR = -log10(FDR)
  ) %>%
  filter(cluster12 > 1)

# añadir descripción a resultados
res$description <- COG_fun(res$COG)
res_long <- wide_norm %>%
  dplyr::select(description, cluster12, `3`) %>%
  pivot_longer(cols = c(cluster12, `3`),
               names_to = "cluster",
               values_to = "count")

p1 <- ggplot(res_long,
             aes(x = description,
                 y = count,
                 fill = cluster)) +
  geom_col(position = position_dodge()) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Anotaciones por genoma (log10)",
       x = "",
       fill = "Cluster",
       title = "COG - Motility (normalizado)") +
  scale_fill_manual(values = c("3" = "#8dd3c7",
                               "cluster12" = "#bebada")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p1

sig_res <- res %>%
  filter(!is.na(description)) %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice_head(n = 4)

p2 <- ggplot(sig_res %>% arrange(desc(logFDR)),
             aes(x = reorder(description, logFDR),
                 y = logFDR)) +
  geom_col(fill = "#ccebc5") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed") +
  coord_flip() +
  theme_minimal() +
  labs(y = "-log10(FDR)",
       x = "",
       title = "COG enrichment (cluster12 vs 3)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2

counts_long <- wide_norm %>%
  dplyr::select(description, `1`, `2`, `3`) %>%
  tidyr::pivot_longer(
    cols = c(`1`, `2`, `3`),
    names_to = "cluster",
    values_to = "count"
  ) %>%
  group_by(description) %>%
  mutate(mean_count = mean(count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(description = reorder(description, mean_count, FUN = mean, decreasing = TRUE))
p3 <- ggplot(counts_long,
             aes(x = description,
                 y = count,
                 fill = cluster)) +
  geom_col(position = position_dodge()) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Anotaciones por genoma (log10)",
       x = "",
       fill = "Cluster",
       title = "COG categories by cluster (normalizado)") +
  scale_fill_manual(values = c("1" = "#9EC3E6",
                               "2" = "#F8E8A0",
                               "3" = "#C53C2F")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3

#---- Fig5d----
data <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/FW_ Swarming/Summary_swarming_motilit_dataR.csv")
data <- data %>%
  mutate(
    cluster = factor(cluster, levels = c("C1", "C2", "C3")),
    time = factor(time, levels = c("48h", "120h", "160h"))
  )
plot_time <- function(df, t) {
  df %>%
    filter(time == t) %>%
    ggplot(aes(x = cluster, y = radius, fill = cluster)) +
    
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 1.8, alpha = 0.7) +
    
    scale_fill_manual(values = c("gold", "steelblue", "red")) +
    
    theme_classic() +
    theme(legend.position = "none") +
    
    labs(
      title = paste("Swarming at", t),
      x = "Cluster",
      y = "Radius"
    )
}
p48  <- plot_time(data, "48h")
p120 <- plot_time(data, "120h")
p160 <- plot_time(data, "160h")

p48
p120
p160

kruskal_results <- data %>%
  group_by(time) %>%
  kruskal_test(radius ~ cluster)
kruskal_results

dunn_results <- data %>%
  group_by(time) %>%
  dunn_test(radius ~ cluster, p.adjust.method = "bonferroni")
dunn_table <- dunn_results %>%
  mutate(
    signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )
dunn_table
data <- data %>%
  arrange(cluster, strain) %>%
  mutate(strain = factor(strain, levels = unique(strain)))
pall<-ggplot(data, aes(x = strain, y = radius, fill = cluster)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.6) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  facet_wrap(~time, scales = "free_x") +
  scale_fill_manual(values = c("gold", "steelblue", "red")) +
  theme_classic() +
  labs(
    x = "Strain",
    y = "Swarming radius",
    fill = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
data <- data %>%
  mutate(
    time = factor(time, levels = c("48h", "120h", "160h"))
  )
C1_strains <- c("C396", "C8", "C88")
C2_strains <- c("C5", "C136", "C21")
C3_strains <- c("C360", "C385")
C1_colors <- colorRampPalette(c("#FFFFCC", "#FFD700"))(length(C1_strains))  # amarillos
C2_colors <- colorRampPalette(c("#CCE5FF", "#1E90FF"))(length(C2_strains))  # azules
C3_colors <- colorRampPalette(c("#FFCCCC", "#FF0000"))(length(C3_strains))  # rojos

strain_colors <- c(C1_colors, C2_colors, C3_colors)
names(strain_colors) <- c(C1_strains, C2_strains, C3_strains)

p<-ggplot(data, aes(x = time, y = radius, group = strain, color = strain)) +
  geom_line(size = 1) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = strain_colors) +
  theme_classic() +
  labs(
    x = "Time",
    y = "Swarming radius",
    color = "Bacteria"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cluster_stats <- data %>%
  group_by(cluster, time) %>%
  summarise(
    mean_radius = mean(radius, na.rm = TRUE),
    sd_radius = sd(radius, na.rm = TRUE),
    .groups = "drop"
  )
time0 <- data.frame(
  cluster = unique(data$cluster),
  time = "0h",
  mean_radius = 0,
  sd_radius = 0
)
cluster_stats <- bind_rows(time0, cluster_stats)
cluster_stats$time <- factor(cluster_stats$time, levels = c("0h", "48h", "120h", "160h"))
cluster_colors <- c("C1" = "gold", "C2" = "steelblue", "C3" = "red")
p1<-ggplot(cluster_stats, aes(x = time, y = mean_radius, group = cluster, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = mean_radius - sd_radius, ymax = mean_radius + sd_radius),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  labs(
    x = "Time",
    y = "Swarming radius",
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2<-ggplot(cluster_stats, aes(x = time, y = mean_radius, group = cluster, color = cluster)) +
  geom_errorbar(aes(ymin = mean_radius - sd_radius, ymax = mean_radius + sd_radius),
                width = 0.2, size = 1) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  labs(
    x = "Time",
    y = "Swarming radius",
    color = "Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

data <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/FW_ Swarming/Swarming_120h/Swarming_120h/20260213Cluster.csv")

df_long <- pivot_longer(data,
                        cols = everything(),
                        names_to = "Cluster",
                        values_to = "Valor")

p<-ggplot(df_long, aes(x = Cluster, y = Valor, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#9ecae1", "#fddc8c", "#f8766d")) +
  theme_minimal() +
  labs(title = "Comparación de clusters")

#---- Fig5b----
data <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/FW_ Swarming/swimming_24_48_72h/swimming_24_48_72h/all_swimming_3_48h.csv")

df_long <- pivot_longer(data,
                        cols = everything(),
                        names_to = "Cluster",
                        values_to = "Valor")

p<-ggplot(df_long, aes(x = Cluster, y = Valor, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("#9ecae1", "#fddc8c", "#f8766d")) +
  theme_minimal() +
  labs(title = "Comparación de clusters")
p

#---- Fig 5f & 5g----
bact_tab <- read.table("output_table_automatic.txt", sep="\t", dec=".", header=T)
plant_tab <- read.table("plant_reads_per_sample.txt", sep="\t", dec=",", header=T)
bact_genome_length <- read.table("bact_genome_lengthRoot.txt", sep="\t", dec=".", header=T)
plant_genome_length <- read.table("plant_genome_length.txt", sep="\t", dec=".", header=T)

#trato las tablas
rownames(bact_tab) <- bact_tab$Name #esto es para convertir la primera columna
bact_tab<-bact_tab%>% select(-Name) #esto para eliminar la columna con el codigo de OTUs
rownames(plant_tab) <- plant_tab$Name #esto es para convertir la primera columna
plant_tab<-plant_tab%>% select(-Name) #esto para eliminar la columna con el codigo de OTUs
rownames(bact_genome_length) <- bact_genome_length$Genome #esto es para convertir la primera columna
bact_genome_length<-bact_genome_length%>% select(-Genome) #esto para eliminar la columna con el codigo de OTUs
rownames(plant_genome_length) <- plant_genome_length$Genome #esto es para convertir la primera columna
plant_genome_length<-plant_genome_length%>% select(-Genome) #esto para eliminar la columna con el codigo de OTUs

#Calculo log10 de la longitud genómica ---
bact_genome_length$Log10_Length <- log10(bact_genome_length$Length)  # Ajusta 'Length' si tu columna tiene otro nombre
plant_genome_length$Log10_Length <- log10(plant_genome_length$Length)

#Normalizo por longitud genómica ---
bact_norm <- sweep(bact_tab, 1, bact_genome_length$Log10_Length[match(rownames(bact_tab), rownames(bact_genome_length))],"/")
plant_norm <- sweep(plant_tab, 1, plant_genome_length$Log10_Length[match(rownames(plant_tab), rownames(plant_genome_length))], "/")

#Normalizar por total de reads (abundancia relativa) por muestra ---
bact_rel_abund <- sweep(bact_norm, 2, colSums(bact_norm), "/") * 100
plant_rel_abund <- sweep(plant_norm, 2, colSums(plant_norm), "/") * 100

# Combinamos en una sola tabla (asegurando que las filas tengan nombres únicos)
missing_cols <- setdiff(names(plant_norm), names(bact_norm))
bact_norm[missing_cols] <- NA
missing_cols <- setdiff(names(bact_norm), names(plant_norm))
plant_norm[missing_cols] <- NA
combined_norm <- rbind(bact_norm, plant_norm)

#normalizamos por muestra, pero usando TODA la tabla
combined_rel_abund <- sweep(combined_norm, 2, colSums(combined_norm), "/") * 100

#Guardar resultados
write.csv(bact_rel_abund, "bacteria_normalizada.csv", quote=FALSE)
write.csv(plant_rel_abund, "planta_normalizada.csv", quote=FALSE)
write.csv(combined_rel_abund, "combinada_normalizada.csv", quote=FALSE)

#---- Phyloseq Object Combinado
otu_table <- otu_table(combined_rel_abund, taxa_are_rows = TRUE)

taxap <- read.delim("plant_bact_taxonomy.csv", sep = ";", header = TRUE)
rownames(taxap) <- taxap$Genome
tax_table <- tax_table(as.matrix(taxap))

metadata <- read.table("Metadata_SequencedSamples_Root.csv", sep = ";", header = TRUE, row.names = 1)
sample_data <- sample_data(metadata)

ps_plant <- phyloseq(otu_table, tax_table, sample_data)

#---- Taxabarplot plant per sample
ps_plant_exp<- subset_samples(ps_plant, Experiment == "RootCompetence") # como son muchas muestras, filtrar para quedarte solo con las muestras de cada experimento y analizarlo por separado
ps_plant_exp_rel <- transform_sample_counts(ps_plant_exp, function(x) x / sum(x)) #Convertir a abundancia relativa
df <- psmelt(ps_plant_exp_rel) #Pasar a data frame

palette<- c("HV"="lightgreen",
            "TA"="darkgreen",
            "Bacteria"= "grey")

p<-ggplot(df, aes(x = sample_code, y = Abundance, fill = kingdom)) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))+
  scale_fill_manual(values = palette) + 
  labs(
    x = "Sample",
    y = "Relative Abundance",
    title = "Plant abundance in Root Competence samples"
  )
p

#---- Phyloseq Object Bact
otu_table <- otu_table(bact_rel_abund, taxa_are_rows = TRUE)

taxap <- read.delim("bact_genomes_taxonomy_Root.csv", sep = ";", header = TRUE)
rownames(taxap) <- taxap$Genome
tax_table <- tax_table(as.matrix(taxap))

metadata <- read.table("Metadata_SequencedSamplesRoot.csv", sep = ";", header = TRUE, row.names = 1)
sample_data <- sample_data(metadata)

ps_bact <- phyloseq(otu_table, tax_table, sample_data)
ps_bact_RC<- subset_samples(ps_bact, Experiment == "RootCompetence") # como son muchas muestras, filtrar para quedarte solo con las muestras de cada experimento y analizarlo por separado

#---- Fig 5f----
df <- psmelt(ps_bact_RC)  
meta <- data.frame(sample_data(ps_bact_RC))
meta_long <- meta %>%
  select(SampleType,
         G1 = Pantoea_G1,
         G2 = Pantoea_G2,
         G3 = Pantoea_G3) %>%
  mutate(across(c(G1, G2, G3), as.character)) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(
    cols = c(G1, G2, G3),
    names_to = "Group",
    values_to = "Genome"
  )

df_filtered <- inner_join(df, meta_long, by = c("Sample", "OTU" = "Genome"))
genome_group_lookup <- meta_long %>%
  select(Genome, Group) %>%
  distinct()
summary_df <- df_filtered %>%
  group_by(SampleType.x, OTU) %>%
  summarise(
    mean_abundance = mean(Abundance),
    sd_abundance = sd(Abundance),
    n = n(),
    .groups = "drop"
  ) %>%
  left_join(genome_group_lookup, by = c("OTU" = "Genome"))

df_norm <- df %>%
  filter(!is.na(Genome), Abundance > 0) %>%
  mutate(
    Time = case_when(
      SampleType == "Inocula" ~ "T_inicial",
      SampleType == "Sample" ~ "T_final"
    )
  ) %>%
  group_by(SampleName, Genome, Group, Time) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Time,
    values_from = Abundance
  ) %>%
  mutate(
    log2FC = log2((T_final + 1e-6) / (T_inicial + 1e-6))
  )
df_norm <- df_norm %>%
  mutate(
    Genome = factor(Genome),
    SampleName = factor(SampleName)
  )
p<-ggplot(df_norm, aes(
  x = Genome,
  y = SampleName,
  fill = log2FC
)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2C3E75",
    mid = "white",
    high = "#8B1E3F",
    midpoint = 0,
    name = "log2(Final / Inóculo)"
  ) +
  facet_wrap(~ Group, scales = "free_x") +
  labs(
    title = "Cambio relativo de abundancia respecto al inóculo",
    x = "Genome",
    y = "Pantoea combination"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p
#---- Fig 5g----
df_log2fc <- df %>%
  filter(!is.na(Genome), Abundance > 0) %>%
  mutate(
    Time = case_when(
      SampleType == "Inocula" ~ "T_inicial",
      SampleType == "Sample" ~ "T_final"
    )
  ) %>%
  group_by(SampleName, Genome, Group, Time) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Time,
    values_from = Abundance
  ) %>%
  mutate(
    log2FC = log2((T_final + 1e-6) / (T_inicial + 1e-6))
  )
p <- ggplot(
  df_log2fc,
  aes(x = Group, y = log2FC, fill = Group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.6,
    color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.3,
    alpha = 0.6
  ) +
  scale_fill_manual(
    values = c("G1" = "#9ecae1", "G2" = "#74c476", "G3" = "#fd8d3c")
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2,
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    x = "Genomic Cluster",
    y = "log2(Final / Inóculo)",
    title = "Colonization per Genomic Cluster (log2FC)"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p

kruskal.test(log2FC ~ Group, data = df_log2fc)
kw_res <- kruskal.test(log2FC ~ Group, data = df_log2fc)
kw_res$p.value
pairwise_res <- pairwise.wilcox.test(
  x = df_log2fc$log2FC,
  g = df_log2fc$Group,
  p.adjust.method = "BH"
)

pairwise_res

df_log2fc %>%
  group_by(Group) %>%
  summarise(
    p_vs_zero = wilcox.test(log2FC, mu = 0)$p.value,
    median_log2FC = median(log2FC)
  ) %>%
  mutate(p_adj = p.adjust(p_vs_zero, method = "BH"))

#---- Fig4a & 4b & Extended Table 4----
meta <- read.csv("C:/Users/ipuente/Desktop/02_PaggPan/RE_ Script Heatmap/metadata_tesis_gephi.csv",
                 stringsAsFactors = FALSE)

traits <- read.delim("C:/Users/ipuente/Desktop/02_PaggPan/RE_ Script Heatmap/trait_summary_ALL_noC0.tsv",
                     sep = "\t",
                     stringsAsFactors = FALSE)
gephi_merge <- merge(meta, traits, by = "ID")

metadata <- gephi_merge %>%
  select(ID, ji_group, isolation.source, strainnumber)

data <- gephi_merge %>%
  select(-ji_group, -isolation.source)
trait_map <- read.delim(
  "C:/Users/ipuente/Desktop/02_PaggPan/RE_ Script Heatmap/trait_rename.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
trait_map <- trait_map %>% filter(Class == "t")
trait.labels <- setNames(trait_map$NEW, trait_map$OLD)

colnames(data) <- ifelse(
  colnames(data) %in% names(trait.labels),
  trait.labels[colnames(data)],
  colnames(data)
)
data_traits <- data[, intersect(colnames(data), trait_map$NEW)]

data_traits <- as.data.frame(ifelse(data_traits > 0, 1, 0))
data_mat <- as.matrix(data_traits)

metadata <- metadata[match(gephi_merge$ID, metadata$ID), ]
rownames(data_mat) <- metadata$strainnumber
rownames(metadata) <- metadata$strainnumber
order_idx <- order(metadata$ji_group)

data_mat <- data_mat[order_idx, ]
metadata <- metadata[order_idx, ]
annot_row <- data.frame(
  JI_Group = metadata$ji_group
)
rownames(annot_row) <- rownames(data_mat)
group_colors <- c(
  "C1" = "#9ecae1",
  "C2" = "#fddc8c",
  "C3" = "tomato3"
)

annot_colors <- list(
  JI_Group = group_colors
)


col_clust <- hclust(dist(t(data_mat)), method = "complete")

pheatmap(
  data_mat,
  cluster_rows = FALSE,
  cluster_cols = col_clust,
  annotation_row = annot_row,
  annotation_colors = annot_colors,
  color = c("#2C3E75", "#8B1E3F"),
  breaks = c(-0.5, 0.5, 1.5),
  cellwidth = 6,
  cellheight = 4,
  fontsize_row = 4,
  fontsize_col = 8,
  border_color = NA,
  gaps_row = cumsum(table(metadata$ji_group)),
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  filename = "heatmap_ALL_t.pdf"
)

#COLLAPSED OKEY
gephi_merge <- merge(meta, traits, by = "ID")

metadata <- gephi_merge %>%
  select(ID, ji_group, isolation.source, strainnumber)

data <- gephi_merge %>%
  select(-ji_group, -isolation.source)

trait_map <- read.delim(
  "C:/Users/ipuente/Desktop/02_PaggPan/RE_ Script Heatmap/trait_rename.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

trait_map <- trait_map %>% filter(Class == "t")
trait.labels <- setNames(trait_map$NEW, trait_map$OLD)

colnames(data) <- ifelse(
  colnames(data) %in% names(trait.labels),
  trait.labels[colnames(data)],
  colnames(data)
)

data_traits <- data[, intersect(colnames(data), trait_map$NEW)]
data_traits <- as.data.frame(ifelse(data_traits > 0, 1, 0))
data_mat <- as.matrix(data_traits)

metadata <- metadata[match(gephi_merge$ID, metadata$ID), ]

rownames(data_mat) <- metadata$strainnumber
rownames(metadata) <- metadata$strainnumber

order_idx <- order(metadata$ji_group)
data_mat <- data_mat[order_idx, ]
metadata <- metadata[order_idx, ]

pattern_df <- data.frame(
  strainnumber = rownames(data_mat),
  ji_group = metadata$ji_group,
  pattern = apply(data_mat, 1, paste0, collapse = ""),
  stringsAsFactors = FALSE
)

pattern_df$combo <- paste(pattern_df$ji_group, pattern_df$pattern, sep = "_")

data_collapsed <- aggregate(
  data_mat,
  by = list(combo = pattern_df$combo),
  FUN = max
)

rownames(data_collapsed) <- data_collapsed$combo
data_collapsed$combo <- NULL
data_collapsed <- as.matrix(data_collapsed)

strain_groups <- tapply(
  pattern_df$strainnumber,
  pattern_df$combo,
  function(x) paste(x, collapse = ", ")
)

rownames(data_collapsed) <- strain_groups[rownames(data_collapsed)]

row_cluster <- sub("_.*", "", names(strain_groups))
row_cluster <- row_cluster[match(names(strain_groups), names(strain_groups))]

row_cluster <- as.character(row_cluster)
row_cluster <- row_cluster[match(rownames(data_collapsed), strain_groups)]

keep <- !is.na(row_cluster)
data_collapsed <- data_collapsed[keep, , drop = FALSE]
row_cluster <- row_cluster[keep]

order_idx <- order(row_cluster)
data_collapsed <- data_collapsed[order_idx, ]
row_cluster <- row_cluster[order_idx]

group_colors <- c(
  "1" = "#9ecae1",
  "2" = "#fddc8c",
  "3" = "tomato3"
)

row_cluster <- factor(row_cluster, levels = names(group_colors))

annot_row <- data.frame(JI_Group = row_cluster)
rownames(annot_row) <- rownames(data_collapsed)

annot_colors <- list(JI_Group = group_colors)

col_clust <- hclust(dist(t(data_collapsed)), method = "complete")

gaps <- cumsum(table(row_cluster))

pheatmap(
  data_collapsed,
  cluster_rows = FALSE,
  cluster_cols = col_clust,
  annotation_row = annot_row,
  annotation_colors = annot_colors,
  gaps_row = gaps,
  color = c("#2C3E75", "#8B1E3F"),
  breaks = c(-0.5, 0.5, 1.5),
  cellwidth = 6,
  cellheight = 4,
  fontsize_row = 4,
  fontsize_col = 8,
  border_color = NA,
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  filename = "heatmap_PATTERNSt_t.pdf"
)

#NAMES CAMBIADOS
gephi_merge <- merge(meta, traits, by = "ID")

metadata <- gephi_merge %>%
  select(ID, ji_group, isolation.source, strainnumber)

data <- gephi_merge %>%
  select(-ji_group, -isolation.source)

trait_map <- read.delim(
  "C:/Users/ipuente/Desktop/02_PaggPan/RE_ Script Heatmap/trait_rename.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

trait_map <- trait_map %>% filter(Class == "t")

trait.labels <- setNames(trait_map$NEW, trait_map$OLD)

colnames(data) <- ifelse(
  colnames(data) %in% names(trait.labels),
  trait.labels[colnames(data)],
  colnames(data)
)

data_traits <- data[, intersect(colnames(data), trait_map$NEW), drop = FALSE]
data_traits <- as.data.frame(ifelse(data_traits > 0, 1, 0))
data_mat <- as.matrix(data_traits)

metadata <- metadata[match(gephi_merge$ID, metadata$ID), ]

rownames(data_mat) <- metadata$strainnumber
rownames(metadata) <- metadata$strainnumber

order_idx <- order(metadata$ji_group)
data_mat <- data_mat[order_idx, , drop = FALSE]
metadata <- metadata[order_idx, ]

pattern_df <- data.frame(
  strainnumber = rownames(data_mat),
  ji_group = metadata$ji_group,
  pattern = apply(data_mat, 1, paste0, collapse = ""),
  stringsAsFactors = FALSE
)

pattern_df$combo <- paste(pattern_df$ji_group, pattern_df$pattern, sep = "_")
pattern_df <- pattern_df[!is.na(pattern_df$combo), ]

data_mat <- data_mat[pattern_df$strainnumber, , drop = FALSE]

data_collapsed <- aggregate(
  data_mat,
  by = list(combo = pattern_df$combo),
  FUN = max
)

rownames(data_collapsed) <- data_collapsed$combo
data_collapsed$combo <- NULL
data_collapsed <- as.matrix(data_collapsed)

strain_groups <- tapply(
  pattern_df$strainnumber,
  pattern_df$combo,
  function(x) paste(x, collapse = ", ")
)

strain_groups <- strain_groups[!is.na(names(strain_groups))]

common_combos <- intersect(names(strain_groups), rownames(data_collapsed))

data_collapsed <- data_collapsed[common_combos, , drop = FALSE]
strain_groups <- strain_groups[common_combos]

ji_map <- pattern_df %>%
  distinct(combo, ji_group)

ji_map <- ji_map[match(common_combos, ji_map$combo), ]

row_cluster <- ji_map$ji_group
pattern_ids <- paste0("Pattern_", seq_along(common_combos))
#pattern_ids <- paste0("Pattern_", seq(from = 37, length.out = length(common_combos)))
names(pattern_ids) <- common_combos

pattern_n <- sapply(strsplit(strain_groups, ", "), length)
names(pattern_n) <- common_combos

pattern_labels_with_n <- paste0(
  pattern_ids[common_combos],
  " (n=", pattern_n[common_combos], ")"
)
names(pattern_labels_with_n) <- common_combos

rownames(data_collapsed) <- pattern_labels_with_n[common_combos]

pattern_table <- data.frame(
  Pattern = pattern_ids[common_combos],
  Pattern_label = pattern_labels_with_n[common_combos],
  Combo = common_combos,
  JI_Group = row_cluster,
  Strains = unname(strain_groups[common_combos]),
  N_strains = as.integer(pattern_n[common_combos]),
  stringsAsFactors = FALSE
)

pattern_table_long <- data.frame(
  Pattern = rep(pattern_ids[common_combos], pattern_n[common_combos]),
  Pattern_label = rep(pattern_labels_with_n[common_combos], pattern_n[common_combos]),
  JI_Group = rep(row_cluster, pattern_n[common_combos]),
  Strain = unlist(strsplit(strain_groups[common_combos], ", ")),
  stringsAsFactors = FALSE
)

data_for_clust <- apply(data_collapsed, 2, as.numeric)
rownames(data_for_clust) <- rownames(data_collapsed)
data_for_clust[!is.finite(data_for_clust)] <- 0

group_colors <- c(
  "1" = "#9ecae1",
  "2" = "#fddc8c",
  "3" = "tomato3"
)

row_cluster <- factor(row_cluster, levels = names(group_colors))

annot_row <- data.frame(JI_Group = row_cluster)
rownames(annot_row) <- rownames(data_collapsed)

annot_colors <- list(JI_Group = group_colors)

gaps <- cumsum(table(row_cluster))

stopifnot(
  nrow(data_collapsed) > 1,
  ncol(data_collapsed) > 1,
  all(!is.na(rownames(data_collapsed))),
  all(rownames(annot_row) == rownames(data_collapsed))
)

pheatmap(
  data_collapsed,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annot_row,
  annotation_colors = annot_colors,
  gaps_row = gaps,
  color = c("#2C3E75", "#8B1E3F"),
  breaks = c(-0.5, 0.5, 1.5),
  cellwidth = 6,
  cellheight = 4,
  fontsize_row = 4,
  fontsize_col = 8,
  border_color = NA,
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  filename = "heatmap_PATTERNS_t_nameFT.pdf"
)

write.table(
  pattern_table,
  file = "pattern_to_genomes_t.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pattern_table_long,
  file = "pattern_to_genomes_t_long.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#---- Fig4i----
df <- read_csv("DataSet_TESIS/M_C_aa.csv")

df_long <- df %>%
  pivot_longer(
    cols = starts_with("d"),
    names_to = "day",
    values_to = "coti"
  ) %>%
  mutate(
    day = as.numeric(str_remove(day, "d"))
  )
df_long <- df_long %>%
  mutate(
    genotype = factor(
      genotype,
      levels = c("col-0", "abi 5-1", "ein", "jin 1-2", "coi 1-16", "pyr", "ctr1")
    ),
    treatment = factor(treatment)
  )
my_colors <- c(
  "Control" = "black",
  "79" = "#4C72B0",
  "88" = "#F4C56A",
  "143" = "#D65F4A"
)

p_c_aa <- ggplot(
  df_long,
  aes(x = day, y = coti, color = treatment, group = treatment)
) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.5
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.35
  ) +
  
  scale_color_manual(values = my_colors) +
  
  facet_wrap(~ genotype, ncol = 1) +
  
  scale_x_continuous(
    breaks = c(2, 4, 7, 9, 12)
  ) +
  
  theme_bw(base_size = 7, base_family = "Helvetica") +
  
  theme(
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 7
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 7
    ),
    
    axis.text = element_text(
      color = "black",
      size = 7
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 7
    ),
    
    legend.text = element_text(
      size = 7
    ),
    
    panel.spacing = unit(0.15, "cm"),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.2
    ),
    
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Days post inoculation",
    y = "cotiledon development",
    color = "Treatment"
  )


p_c_aa
#---- Extended Figure 2a----

class_counts <- data.frame(
    Class = c("Core", "Soft-core", "Shell", "Cloud"),
    Percent = c(10.2, 0.9, 3.4, 85.5)
  )

class_counts$ymax <- cumsum(class_counts$Percent)
class_counts$ymin <- c(0, head(class_counts$ymax, n = -1))

class_counts$label_pos <- (class_counts$ymax + class_counts$ymin) / 2

p <- ggplot(class_counts,
            aes(ymax = ymax,
                ymin = ymin,
                xmax = 4,
                xmin = 2,
                fill = Class)) +
  
  geom_rect(color = "white", linewidth = 0.7) +
  
  coord_polar(theta = "y") +
  
  xlim(c(0, 4)) +
  
  geom_text(
    aes(
      x = 3,
      y = label_pos,
      label = paste0(Percent, "%")
    ),
    size = 4
  ) +
  
  scale_fill_manual(values = c(
    "Core" = "#3B6FB6",
    "Soft-core" = "#F2C14E",
    "Shell" = "#B8B8B8",
    "Cloud" = "#D96C6C"
  )) +
  
  theme_void(base_size = 12) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.margin = margin(10,10,10,10)
  )

p





#---- Extended Figure 2h----
gpa <- read_csv(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/Pangenoma_NO_split_prokka/gene_presence_absence.csv",
  show_col_types = FALSE
)
genome_cols <- colnames(gpa)[15:ncol(gpa)]
n_genomes <- length(genome_cols)

gpa_class <- gpa %>%
  mutate(
    presence_count = rowSums(!is.na(across(all_of(genome_cols))) & across(all_of(genome_cols)) != ""),
    frequency = presence_count / n_genomes,
    Class = case_when(
      frequency >= 0.99 ~ "Core",
      frequency >= 0.95 ~ "Soft-core",
      frequency >= 0.15 ~ "Shell",
      TRUE ~ "Cloud"
    )
  )
table(gpa_class$Class)
prop.table(table(gpa_class$Class)) * 100
faa <- readAAStringSet(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/all_proteins.faa"
)

head(names(faa))
head(gpa_class[1, genome_cols[1:5]])
protein_df <- data.frame(
  fasta_header = names(faa),
  protein_id = sub("^[^|]+\\|", "", names(faa)),
  length_aa = width(faa),
  stringsAsFactors = FALSE
)
protein_df$protein_id <- sub(" .*", "", protein_df$protein_id)

head(protein_df)
gpa_long <- gpa_class %>%
  dplyr::select(Gene, Class, all_of(genome_cols)) %>%
  pivot_longer(
    cols = all_of(genome_cols),
    names_to = "Genome",
    values_to = "Protein_ids"
  ) %>%
  filter(!is.na(Protein_ids), Protein_ids != "")
gpa_long <- gpa_long %>%
  separate_rows(Protein_ids, sep = "\\s+|;|,") %>%
  filter(Protein_ids != "")
protein_lengths_class <- gpa_long %>%
  left_join(protein_df, by = c("Protein_ids" = "protein_id"))
sum(!is.na(protein_lengths_class$length_aa))
sum(is.na(protein_lengths_class$length_aa))
mean(!is.na(protein_lengths_class$length_aa)) * 100
plot_df <- protein_lengths_class %>%
  filter(!is.na(length_aa))
plot_df$Class <- factor(
  plot_df$Class,
  levels = c("Core", "Soft-core", "Shell", "Cloud")
)

p <- ggplot(plot_df, aes(x = Class, y = length_aa, fill = Class)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    color = NA,
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.3,
    alpha = 0.7
  ) +
  coord_cartesian(
    ylim = c(0, quantile(plot_df$length_aa, 0.99, na.rm = TRUE))
  ) +
  scale_fill_manual(values = c(
    "Core" = "#3B6FB6",
    "Soft-core" = "#F2C14E",
    "Shell" = "#B8B8B8",
    "Cloud" = "#D96C6C"
  )) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Protein length (aa)"
  )

p
kruskal_res <- plot_df %>%
  kruskal_test(length_aa ~ Class)

kruskal_res
pairwise_res <- plot_df %>%
  pairwise_wilcox_test(
    length_aa ~ Class,
    p.adjust.method = "BH"
  )

pairwise_res
summary_lengths <- plot_df %>%
  group_by(Class) %>%
  summarise(
    n = n(),
    median_length = median(length_aa, na.rm = TRUE),
    mean_length = mean(length_aa, na.rm = TRUE),
    q25 = quantile(length_aa, 0.25, na.rm = TRUE),
    q75 = quantile(length_aa, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

summary_lengths
#---- Extended Figure 2g----
faa <- readAAStringSet(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/all_proteins.faa"
)

headers <- names(faa)

genome_ids <- str_extract(headers, "^[^|]+")

cds_df <- data.frame(
  Genome = genome_ids
) %>%
  count(Genome, name = "Coding_genes")
p_cds <- ggplot(cds_df, aes(x = Coding_genes)) +
  
  geom_histogram(
    bins = 25,
    fill = "#7BAF9E",
    color = "white",
    linewidth = 0.3
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Number of coding genes",
    y = "Genomes"
  )

p_cds
#---- Fig2c----
aln <- read.dna(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/core_gene_alignment_pv.aln",
  format = "fasta"
)

genomes <- rownames(aln)

ji <- read.csv(
  "C:/Users/ipuente/Desktop/02_PaggPan/DataSet_TESIS/Pangenoma_NO_split_prokka/ji_grupo.csv",
  sep = ";",
  stringsAsFactors = FALSE
)

ji <- ji[match(genomes, ji$Id), ]

valid <- !is.na(ji$JI_Group)

aln <- aln[valid, ]
ji <- ji[valid, ]

pop_factor <- factor(paste0("C", ji$JI_Group))
window_size <- 10000
alignment_length <- ncol(aln)
windows <- seq(1, alignment_length, by = window_size)

pi_window_df <- list()

for (pop in levels(pop_factor)) {
  
  dna_pop <- aln[pop_factor == pop, ]
  
  vals <- data.frame()
  
  for (start in windows) {
    end <- min(start + window_size - 1, alignment_length)
    window_data <- dna_pop[, start:end]
    
    pi_val <- tryCatch(
      nuc.div(window_data, quiet = TRUE),
      error = function(e) NA
    )
    
    vals <- rbind(vals, data.frame(
      Cluster = pop,
      Start = start,
      End = end,
      Midpoint = (start + end) / 2,
      Pi = pi_val
    ))
  }
  
  pi_window_df[[pop]] <- vals
}

pi_window_df <- bind_rows(pi_window_df)
p_pi <- ggplot(pi_window_df, aes(x = Midpoint / 1000, y = Pi, color = Cluster)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~Cluster, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = "Posición en genoma core (kb)",
    y = expression(pi~"nucleotide diversity")
  )

p_pi
#---- sitios segregantes
seg_window_df <- list()

for (pop in levels(pop_factor)) {
  
  dna_pop <- aln[pop_factor == pop, ]
  dna_char <- as.character(dna_pop)
  
  vals <- data.frame()
  
  for (start in windows) {
    end <- min(start + window_size - 1, alignment_length)
    
    window_char <- dna_char[, start:end, drop = FALSE]
    
    seg_count <- sum(apply(window_char, 2, function(site) {
      alleles <- unique(site[site %in% c("a", "c", "g", "t")])
      length(alleles) > 1
    }))
    
    seg_fraction <- seg_count / (end - start + 1)
    
    vals <- rbind(vals, data.frame(
      Cluster = pop,
      Start = start,
      End = end,
      Midpoint = (start + end) / 2,
      Seg_sites = seg_count,
      Seg_fraction = seg_fraction
    ))
  }
  
  seg_window_df[[pop]] <- vals
}

seg_window_df <- bind_rows(seg_window_df)
p_seg <- ggplot(seg_window_df, aes(x = Midpoint / 1000, y = Seg_fraction, color = Cluster)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  theme_classic(base_size = 11) +
  labs(
    x = "Posición en genoma core (kb)",
    y = "Fracción sitios segregantes"
  )

p_seg

#----Cotiledones----
                                            file <- "COTILEDONES_ISP.xlsx"
sheets <- excel_sheets(file)
read_cotyledon_sheet <- function(sheet_name) {
  
  day <- str_extract(sheet_name, "\\d+") |> as.integer()
  
  read_excel(
    file,
    sheet = sheet_name,
    skip = 1
  ) %>%
    rename(
      genotype = 1,
      replicate = 2
    ) %>%
    pivot_longer(
      cols = -c(genotype, replicate),
      names_to = "strain",
      values_to = "cotyledons"
    ) %>%
    mutate(
      day = day
    )
}


cotyledons <- map_dfr(
  sheets,
  read_cotyledon_sheet
)
strain_info <- tibble(
  strain = c(
    "C396", "C8", "C88",
    "C136", "C5", "C21",
    "C385", "C143", "C360"
  ),
  cluster = c(
    "C1", "C1", "C1",
    "C2", "C2", "C2",
    "C3", "C3", "C3"
  )
)


cotyledons <- cotyledons %>%
  mutate(
    strain = recode(
      strain,
      "C121" = "C21"
    )
  ) %>%
  left_join(
    strain_info,
    by = "strain"
  )
cotyledons <- cotyledons %>%
  mutate(
    treatment = case_when(
      strain == "PBS" ~ "Control",
      strain == "C79" ~ "C79",
      cluster == "C1" ~ "C1",
      cluster == "C2" ~ "C2",
      cluster == "C3" ~ "C3",
      TRUE ~ NA_character_
    )
  )
summary_cotyledons <- cotyledons %>%
  filter(!is.na(treatment)) %>%
  group_by(genotype, treatment, day) %>%
  summarise(
    mean = mean(cotyledons, na.rm = TRUE),
    se = sd(cotyledons, na.rm = TRUE) /
      sqrt(sum(!is.na(cotyledons))),
    n = sum(!is.na(cotyledons)),
    .groups = "drop"
  )
my_colors <- c(
  "Control" = "black",
  "C79" = "grey",
  "C1" = "#F4C56A",
  "C2" = "#4C72B0",
  "C3" = "#D65F4A"
)


p_cotyledons <- ggplot(
  summary_cotyledons,
  aes(
    x = day,
    y = mean,
    color = treatment,
    group = treatment
  )
) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.25,
    linewidth = 0.35
  ) +
  scale_color_manual(
    values = my_colors
  ) +
  facet_wrap(
    ~ genotype,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = c(2, 5, 7, 9, 12)
  ) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, 2)
  ) +
  theme_bw(
    base_size = 7,
    base_family = "Helvetica"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 7
    ),
    axis.title = element_text(
      face = "bold",
      size = 7
    ),
    axis.text = element_text(
      color = "black",
      size = 7
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 7
    ),
    legend.text = element_text(
      size = 7
    ),
    panel.spacing = unit(
      0.15,
      "cm"
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Days post inoculation",
    y = "Cotyledon-bearing seedlings",
    color = "Treatment"
  )


p_cotyledons
cotyledon_stats <- cotyledons %>%
  filter(!is.na(treatment)) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "Control",
        "C79",
        "C1",
        "C2",
        "C3"
      )
    ),
    day = as.numeric(day)
  )
model_cotyledons <- lmer(
  cotyledons ~ treatment * day * genotype + (1 | strain),
  data = cotyledon_stats
)
anova(model_cotyledons)
emm_cotyledons <- emmeans(
  model_cotyledons,
  ~ treatment | genotype * day
)
contrasts_cotyledons <- contrast(
  emm_cotyledons,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "holm"
)
summary(contrasts_cotyledons)
#----germination----
file <- "GERMINACIÓN_ISP.xlsx"
sheets <- excel_sheets(file)
read_germination_sheet <- function(sheet_name) {
  
  day <- str_extract(sheet_name, "\\d+") |> as.integer()
  
  read_excel(
    file,
    sheet = sheet_name,
    skip = 1
  ) %>%
    rename(
      genotype = 1,
      replicate = 2
    ) %>%
    pivot_longer(
      cols = -c(genotype, replicate),
      names_to = "strain",
      values_to = "germinated"
    ) %>%
    mutate(
      day = day
    )
}

germination_new <- map_dfr(
  sheets,
  read_germination_sheet
)
strain_info <- tibble(
  strain = c(
    "C396", "C8", "C88",
    "C136", "C5", "C121",
    "C385", "C143", "C360"
  ),
  cluster = c(
    "C1", "C1", "C1",
    "C2", "C2", "C2",
    "C3", "C3", "C3"
  )
)

germination_new <- germination_new %>%
  mutate(
    strain = recode(
      strain,
      "C121" = "C21"
    )
  ) %>%
  left_join(strain_info, by = "strain")
strain_info <- tibble(
  strain = c(
    "C396", "C8", "C88",
    "C136", "C5", "C21",
    "C385", "C143", "C360"
  ),
  cluster = c(
    "C1", "C1", "C1",
    "C2", "C2", "C2",
    "C3", "C3", "C3"
  )
)

# Read all data and assign clusters
germination_new <- map_dfr(
  sheets,
  read_germination_sheet
) %>%
  mutate(
    strain = recode(strain, "C121" = "C21")
  ) %>%
  left_join(strain_info, by = "strain") %>%
  mutate(
    treatment = case_when(
      strain == "PBS" ~ "Control",
      strain == "C79" ~ "C79",
      cluster == "C1" ~ "C1",
      cluster == "C2" ~ "C2",
      cluster == "C3" ~ "C3",
      TRUE ~ NA_character_
    )
  )
summary_germination <- germination_new %>%
  filter(!is.na(treatment)) %>%
  group_by(genotype, treatment, day) %>%
  summarise(
    mean = mean(germinated, na.rm = TRUE),
    se = sd(germinated, na.rm = TRUE) / sqrt(sum(!is.na(germinated))),
    n = sum(!is.na(germinated)),
    .groups = "drop"
  )
my_colors <- c(
  "Control" = "black",
  "C79" = "grey",
  "C1" = "#F4C56A",
  "C2" = "#4C72B0",
  "C3" = "#D65F4A"
)
p_germination <- ggplot(
  summary_germination,
  aes(
    x = day,
    y = mean,
    color = treatment,
    group = treatment
  )
) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.25,
    linewidth = 0.35
  ) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ genotype, ncol = 1) +
  scale_x_continuous(
    breaks = c(2, 5, 7, 9, 12)
  ) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, 2)
  ) +
  theme_bw(
    base_size = 7,
    base_family = "Helvetica"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7),
    axis.title = element_text(face = "bold", size = 7),
    axis.text = element_text(color = "black", size = 7),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 7),
    panel.spacing = unit(0.15, "cm"),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Days post inoculation",
    y = "Germinated seeds",
    color = "Treatment"
  )

p_germination
germination_stats <- germination_new %>%
  filter(!is.na(treatment)) %>%
  mutate(
    treatment = factor(treatment,
                       levels = c("Control", "C79", "C1", "C2", "C3")),
    day = as.numeric(day)
  )
model <- lmer(
  germinated ~ treatment * day * genotype + (1 | strain),
  data = germination_stats
)

anova(model)
emm <- emmeans(
  model,
  ~ treatment | genotype * day
)

contrasts <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "holm"
)

summary(contrasts)        
#----Rarefraccion----
pantoea <- data.frame(
  agglomerans = c(33, 35, 7, 28, 57, 78, 48),
  ananatis = c(2, 0, 2, 0, 0, 0, 0),
  pleuroti = c(0, 0, 0, 0, 1, 1, 1),
  allii = c(0, 1, 1, 0, 1, 0, 0),
  vagans = c(5, 2, 0, 4, 7, 2, 2),
  beijingensis = c(1, 0, 0, 0, 0, 0, 0),
  rodasii = c(1, 0, 0, 0, 0, 0, 0)
)          
rownames(pantoea) <- c(
  "Aegilops tauschii",
  "T. monococcum",
  "T. turgidum",
  "T. dicoccum",
  "T. durum",
  "T. aestivum",
  "T. spelta"
)

observed_richness <- specnumber(pantoea)
observed_richness

rare_12 <- rarefy(pantoea, sample = 12)
rare_12

results <- data.frame(
  Host = rownames(pantoea),
  Total_isolates = rowSums(pantoea),
  Observed_species = observed_richness,
  Rarefied_species_12 = as.numeric(rare_12)
)
results

rare_curves <- rarecurve(
  pantoea,
  step = 1,
  sample = 12,
  label = TRUE #Poner FALSE si no las quiero
)
rare_df <- map_dfr(
  rownames(pantoea),
  function(host) {
    
    max_n <- sum(pantoea[host, ])
    
    data.frame(
      Host = host,
      Isolates = 1:max_n,
      Species = as.numeric(
        rarefy(
          pantoea[host, , drop = FALSE],
          sample = 1:max_n
        )
      )
    )
  }
)
p<-ggplot(
  rare_df,
  aes(
    x = Isolates,
    y = Species,
    color = Host
  )
) +
  geom_line(linewidth = 0.8) +
  theme_bw(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  labs(
    x = "Number of Pantoea isolates",
    y = "Species richness",
    color = "Wheat host"
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

