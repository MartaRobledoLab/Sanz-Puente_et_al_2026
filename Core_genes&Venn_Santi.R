setwd("C:/Users/irene/OneDrive/IBBTEC/Paper_Genomas")

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

#----LISTA DE GENES POR CLUSTERS----
# CARGA LA GENE PRESENCE.RTAB
roary_matrix = read_tsv("gene_presence_absence.Rtab", col_names = TRUE)
grupo <- read_delim('ji_grupo.csv', col_names = TRUE, delim=";")
# Añadimos el grupo del genoma que falta. Otra posibilidad es eliminarlo de roary_matrix
grupo <- tibble(Id = 'ISP212TRI626322bis', ji_group = 3) %>% full_join(grupo)

# Crear una lista de genomas POR GRUPO
ji_grupo1 <- grupo %>%
  filter(ji_group == 1) %>%  # Filtra para el grupo 1
  pull(Id)

ji_grupo2 <- grupo %>%
  filter(ji_group == 2) %>%  # Filtra para el grupo 1
  pull(Id)

ji_grupo3 <- grupo %>%
  filter(ji_group == 3) %>%  # Filtra para el grupo 1
  pull(Id)

### Separar los datos de presencia/ausencia SEGUN EL GRUPO AL QUE PERTENEZCAN
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

# Eliminar las filas de Genes que sean 0 para cada grupo de JI
rows_to_keep_C1 <- rowSums(C1_matrix != 0) > 0
rows_to_keep_C2 <- rowSums(C2_matrix != 0) > 0
rows_to_keep_C3 <- rowSums(C3_matrix != 0) > 0

# Filtrar la matriz para eliminar filas que tienen solo ceros
final_C1 <- C1_matrix[rows_to_keep_C1, ]
final_C2 <- C2_matrix[rows_to_keep_C2, ]
final_C3 <- C3_matrix[rows_to_keep_C3, ]

# SI AHORA QUIERES VER QUE GENES CORE AL 100% HAY EN CADA GRUPO
final_rows_to_keep_C1_core <- rowSums(final_C1 != 0) >= 55  # Reemplaza el 0 por el numero de genomas del grupo que quieras analizar
final_rows_to_keep_C2_core <- rowSums(final_C2 != 0) >= 252  
final_rows_to_keep_C3_core <- rowSums(final_C3 != 0) >= 41  

# Filtrar la matriz para eliminar filas que tienen solo ceros
C1_final_core <- final_C1[final_rows_to_keep_C1_core, ]
C2_final_core <- final_C2[final_rows_to_keep_C2_core, ]
C3_final_core <- final_C3[final_rows_to_keep_C3_core, ]

#Ahora quiero añadir la notacion genica a cada gen en cada matriz final por cluster
#Cargo los datos y quito la info de los genomas para que no pese mucho y sea más facil manejarla
notation_matrix = read_csv("gene_presence_absence.csv", col_names = TRUE) %>% select('Gene', 'Annotation')

#Ahora seleciono los genes core de cada cluster y filtro la matriz de notation con esa lista de genes
genes_to_filter_C1_core <- row.names(C1_final_core)
C1_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C1_core, ]
genes_to_filter_C2_core <- row.names(C2_final_core)
C2_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C2_core, ]
genes_to_filter_C3_core <- row.names(C3_final_core)
C3_notation_core <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C3_core, ]

# SI AHORA QUIERES VER QUE GENES CORE AL 99% HAY EN CADA GRUPO
final_rows_to_keep_C1_core99 <- rowSums(final_C1 != 0) >= 0.99*55  # Reemplaza el 0.99 por el threshold para core que quieras aplicar
final_rows_to_keep_C2_core99 <- rowSums(final_C2 != 0) >= 0.99*252  
final_rows_to_keep_C3_core99 <- rowSums(final_C3 != 0) >= 0.99*41  

# Filtrar la matriz para eliminar filas que tienen solo ceros
C1_final_core99 <- final_C1[final_rows_to_keep_C1_core99, ]
C2_final_core99 <- final_C2[final_rows_to_keep_C2_core99, ]
C3_final_core99 <- final_C3[final_rows_to_keep_C3_core99, ]

#Ahora seleciono los genes core de cada cluster y filtro la matriz de notation con esa lista de genes
genes_to_filter_C1_core99 <- row.names(C1_final_core99)
C1_notation_core99 <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C1_core99, ]
genes_to_filter_C2_core99 <- row.names(C2_final_core99)
C2_notation_core99 <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C2_core99, ]
genes_to_filter_C3_core99 <- row.names(C3_final_core99)
C3_notation_core99 <- notation_matrix[notation_matrix$Gene %in% genes_to_filter_C3_core99, ]

# SI AHORA QUIERES VER QUE GENES APARECEN ALGUNA VEZ EN CADA GRUPO
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

# Los grupos 1 y 2 no tienen genes core que no estén en alguno de los otros grupos (en roary_nosplit_prokka)
#rownames(C1_final_core)[which(! rownames(C1_final_core) %in% union(rownames(C2_final_any), rownames(C3_final_any)))]
#rownames(C2_final_core)[which(! rownames(C2_final_core) %in% union(rownames(C1_final_any), rownames(C3_final_any)))]

# El grupo 3 tiene 2 proteinas específicas que no están en los otros 2 grupos (en roary_nosplit_prokka)
# Lo que son 2 proteinas hipotéticas: "group_22169" "group_22244"
rownames(C3_final_core)[which(! rownames(C3_final_core) %in% union(rownames(C1_final_any), rownames(C2_final_any)))]


# Primero voy a eliminar GCF_037152265.1 y GCF_029623495.1 (que son los genoms del C3.1)
library(dplyr)

# Eliminar las columnas de esos genomas
C3.2_final_core <- C3_final_core %>%
  select(-GCF_037152265.1, -GCF_029623495.1)
#ver si tiene cosas especificas q no tienen los clusters 2 y 1
rownames(C3.2_final_core)[which(! rownames(C3_final_core) %in% union(rownames(C1_final_any), rownames(C2_final_any)))]



write.csv(C3_notation_any, file = "C3_ALL_genes_99.csv", row.names = FALSE)

#----DIAGRAMA DE VENN----
library(VennDiagram)
genes_cluster1 <- C1_notation_core99$Gene  
genes_cluster2 <- C2_notation_core99$Gene
genes_cluster3 <- C3_notation_core99$Gene
write.csv(data.frame(Gene = genes_cluster1),
          "genes_cluster1.csv", row.names = FALSE)

write.csv(data.frame(Gene = genes_cluster2),
          "genes_cluster2.csv", row.names = FALSE)

write.csv(data.frame(Gene = genes_cluster3),
          "genes_cluster3.csv", row.names = FALSE)

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(Cluster1 = genes_cluster1, Cluster2 = genes_cluster2, Cluster3 = genes_cluster3),
  category.names = c("Cluster 1", "Cluster 2", "Cluster 3"),
  filename = NULL,
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.5), alpha('#21908dff',0.5), alpha('#fde725ff',0.5)),
  cex = 0.7,
  fontfamily = "sans",
  cat.cex = 0.7,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = "black",
  rotation = 1)
grid.draw(venn.plot)  
png("VennDiagram_clusters99.png", width = 800, height = 800, res = 300)

intersection_all <- intersect(intersect(genes_cluster1, genes_cluster2), genes_cluster3)
print(intersection_all)
genes_shared_1_2_not_3 <- intersect(genes_cluster1, genes_cluster2)
genes_shared_1_2_not_3 <- setdiff(genes_shared_1_2_not_3, genes_cluster3)
# Combinar anotaciones de ambos clústeres (por si hay diferencias en la descripción)
combined_annotation <- bind_rows(C1_notation_core99, C2_notation_core99) %>%
  distinct(Gene, .keep_all = TRUE)  # evitar duplicados
# Combinar anotaciones de ambos clústeres (por si hay diferencias en la descripción)
combined_annotation <- bind_rows(C1_notation_core99, C2_notation_core99) %>%
  distinct(Gene, .keep_all = TRUE)  # evitar duplicados
genes_shared_1_2_not_3_data <- combined_annotation %>%
  filter(Gene %in% genes_shared_1_2_not_3)

write.csv(genes_shared_1_2_not_3_data, "genes_shared_1_2_not_3.csv", row.names = FALSE)

  


# Dibujar el diagrama de Venn y guardarlo en el archivo
grid.draw(venn.plot)

# Cerrar el dispositivo gráfico
dev.off()

#----genes solo core por cada cluster----

# Identificar los genes únicos en cada clúster
unique_genes_cluster1 <- setdiff(genes_cluster1, union(genes_cluster2, genes_cluster3))
unique_genes_cluster2 <- setdiff(genes_cluster2, union(genes_cluster1, genes_cluster3))
unique_genes_cluster3 <- setdiff(genes_cluster3, union(genes_cluster1, genes_cluster2))

# Filtrar las filas originales para obtener los genes únicos y sus anotaciones en cada clúster
unique_cluster1_data <- C1_notation_core99 %>% filter(Gene %in% unique_genes_cluster1)
unique_cluster2_data <- C2_notation_core99 %>% filter(Gene %in% unique_genes_cluster2)
unique_cluster3_data <- C3_notation_core99 %>% filter(Gene %in% unique_genes_cluster3)

# Guardar los datos únicos de cada clúster en archivos CSV separados, conservando ambas columnas
write.csv(unique_cluster1_data, "unique_cluster1.csv", row.names = FALSE)
write.csv(unique_cluster2_data, "unique_cluster2.csv", row.names = FALSE)
write.csv(unique_cluster3_data, "unique_cluster3.csv", row.names = FALSE)

#----Ronnie suggestion----
library(clusterProfiler)
library(KEGGREST)

# Enriquecimiento de KEGG pathways
enrich_c1 <- enrichKEGG(
  gene         = genes_cluster1,
  organism     = 'pag',  # KEGG code de Pantoea agglomerans
  pvalueCutoff = 0.05
)
enrich_c2 <- enrichKEGG(
  gene         = genes_cluster2,
  organism     = 'pag',  # KEGG code de Pantoea agglomerans
  pvalueCutoff = 0.05
)
enrich_c3 <- enrichKEGG(
  gene         = genes_cluster3,
  organism     = 'pag',  # KEGG code de Pantoea agglomerans
  pvalueCutoff = 0.05
)