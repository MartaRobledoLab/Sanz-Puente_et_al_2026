## Árboles del Cluster 3 (PAgg) anotados con metadatos de MetadatosC3.xlsx
## Uso:  Rscript plot_arboles_C3.R [ruta_excel] [carpeta_arboles] [carpeta_salida]
## Paquetes: ape, ggtree, patchwork, ggplot2, dplyr, readxl, svglite

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(patchwork)
  library(dplyr); library(readxl)
})

args     <- commandArgs(trailingOnly = TRUE)
xlsx     <- if (length(args) >= 1) args[1] else "MetadatosC3.xlsx"
tree_dir <- if (length(args) >= 2) args[2] else "."
out_dir  <- if (length(args) >= 3) args[3] else file.path(tree_dir, "figuras")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

outgroup <- "GCF_019048385.1"   # prefijo; la etiqueta lleva además _n

## Truncado de ramas: se truncan las que miden más de trunc_factor veces la
## mediana, y se dibujan con trunc_frac veces la profundidad del resto del árbol.
trunc_factor <- 100
trunc_frac   <- 0.3

## Quita las ramas indicadas (pone su longitud a 0) para medir el resto del árbol
drop_edges <- function(tr, e) { tr$edge.length[e] <- 0; tr }

## Columnas del Excel que se dibujan (nombre en el Excel = título en la figura)
cols <- c("Host Species"          = "Host",
          "Genome Ploidy"         = "Ploidy",
          "Domestication status"  = "Domestication",
          "Continent of Origin"   = "Continent",
          "Regeneration location" = "Regeneration site",
          "Strain Isolation Year" = "Isolation year",
          "Provider"              = "Provider",
          "Provider group"        = "Provider group",
          "Sequencing"            = "Sequencing")

## Variables cuya leyenda va en cursiva (nombres de especie)
italic_cols <- c("Host")

## Color del sufijo _n (número de SNPs de cada hoja)
snp_col <- "#B2182B"

palettes <- list(
  "Host"              = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"),
  "Ploidy"            = c("#FDE0A6", "#F4A259", "#B5541F"),
  "Domestication"     = c("#8C6BB1", "#9EBC9F"),
  "Continent"         = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"),
  "Regeneration site" = c("#4E79A7", "#F28E2B", "#59A14F"),
  "Isolation year"    = c("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"),
  "Provider"          = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                          "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                          "#8C564B", "#C49C94", "#E377C2", "#7F7F7F", "#BCBD22",
                          "#17BECF"),
  "Provider group"    = c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
                          "#76B7B2"),
  "Sequencing"        = c("#3B6E8F", "#E3A33B")
)

## --- Metadatos: limpiar espacios duros y espacios sobrantes ---
clean <- function(x) {
  x <- gsub("[  ]", " ", as.character(x))
  x <- gsub("\\s+", " ", trimws(x))
  x <- sub("^Leon$", "León (Spain)", x)
  x <- sub("^Asturias.*$", "Asturias (Spain)", x)
  x
}
meta <- read_excel(xlsx)
names(meta) <- clean(names(meta))
meta <- meta %>% mutate(across(everything(), clean)) %>%
  filter(Cluster == "3") %>%
  ## Tecnología de secuenciación: Seq Run ID que empieza por "Run" = Nanopore
  ## Proveedores agrupados: todos los "IPK (...)" comparten categoría
  mutate(`Provider group` = ifelse(grepl("^IPK", Provider), "IPK", Provider)) %>%
  mutate(Sequencing = ifelse(grepl("^Run", `Seq Run ID`), "Nanopore", "Illumina")) %>%
  select(`Strain ID`, all_of(names(cols)))
names(meta) <- c("label", unname(cols))

plot_tree <- function(tree_file, title) {
  tr <- read.tree(tree_file)
  ## Etiquetas tipo "C-239_0": el sufijo _n es el número de SNPs de la hoja.
  ## Se separa de la etiqueta y se muestra en su propia columna.
  strain <- function(x) sub("_[0-9]+$", "", x)
  snps   <- function(x) sub("^.*_([0-9]+)$", "\\1", x)
  og <- tr$tip.label[strain(tr$tip.label) == outgroup]
  tr$node.label <- NULL
  tr <- root(tr, outgroup = og, resolve.root = TRUE)
  ## Toda la longitud entre la raíz y el outgroup se pone en la rama del
  ## outgroup (en algunos árboles venía en la rama del grupo interno).
  raiz <- Ntip(tr) + 1
  e_raiz <- which(tr$edge[, 1] == raiz)
  e_og <- e_raiz[tr$edge[e_raiz, 2] == which(tr$tip.label == og)]
  tr$edge.length[e_og] <- sum(tr$edge.length[e_raiz])
  tr$edge.length[setdiff(e_raiz, e_og)] <- 0
  tr <- ladderize(tr, right = FALSE)

  falta <- setdiff(setdiff(strain(tr$tip.label), outgroup), meta$label)
  if (length(falta)) warning("Sin metadatos: ", paste(falta, collapse = ", "))

  ## Ramas excesivamente largas (el outgroup) se truncan: se dibujan con una
  ## longitud fija, marcadas con "//" y con su longitud real anotada.
  bl <- tr$edge.length
  thr <- trunc_factor * median(bl[bl > 0])
  larga <- which(bl > thr)
  real <- bl[larga]
  corta_max <- max(node.depth.edgelength(drop_edges(tr, larga)))
  tr$edge.length[larga] <- trunc_frac * corta_max
  nodos_trunc <- tr$edge[larga, 2]

  p <- ggtree(tr, linewidth = 0.4) +
    geom_tiplab(aes(label = strain(label),
                    fontface = ifelse(label == og, "italic", "plain")),
                size = 3, align = TRUE, linesize = 0.2, linetype = "dotted") +
    geom_text(data = function(x) x[x$node %in% nodos_trunc, ],
              aes(x = branch, y = y), label = "//", size = 4, fontface = "bold") +
    geom_text(data = function(x) {
                x <- x[x$node %in% nodos_trunc, ]
                x$lab <- signif(real[match(x$node, nodos_trunc)], 3)
                x },
              aes(x = branch, y = y, label = lab), vjust = -1, size = 2.6) +
    geom_treescale(x = 0, y = -1.1, width = signif(corta_max / 4, 1),
                   fontsize = 2.6, linesize = 0.4, offset = 0.4) +
    hexpand(0.2) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 12))

  ## Una columna de teselas por variable, cada una con su leyenda.
  ## Se alinean con el árbol usando la coordenada y de cada punta.
  ypos <- p$data %>% filter(isTip) %>% transmute(snp = snps(label), label = strain(label), y)
  ylim <- c(-1.6, max(p$data$y) + 0.5)
  d <- meta %>% inner_join(ypos, by = "label")
  tiles <- lapply(unname(cols), function(v) {
    lv <- sort(unique(d[[v]]))
    pal <- setNames(rep_len(palettes[[v]], length(lv)), lv)
    dv <- data.frame(y = d$y, var = v, value = d[[v]])
    ggplot(dv, aes(x = var, y = y, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      scale_fill_manual(values = pal, name = v) +
      scale_x_discrete(position = "top") +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      theme_void() +
      theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9),
            legend.title = element_text(face = "bold", size = 9),
            legend.text  = element_text(size = 8,
                                        face = if (v %in% italic_cols) "italic" else "plain"),
            legend.key.size = unit(0.35, "cm"))
  })
  ## Columna con el número de SNPs de cada hoja (sufijo _n de la etiqueta)
  snp_panel <- ggplot(transmute(ypos, y, snp, var = "Tip SNPs"),
                      aes(x = var, y = y, label = snp)) +
    geom_text(size = 3, fontface = "bold", colour = snp_col) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_void() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9,
                                         colour = snp_col))
  p <- p + scale_y_continuous(limits = ylim, expand = c(0, 0))
  wrap_plots(c(list(p, snp_panel), tiles), nrow = 1,
             widths = c(1, 0.09, rep(0.06, length(tiles)))) +
    plot_layout(guides = "collect")
}

trees <- c("core_SNPs.parsimony" = "Cluster 3: core SNPs (parsimony)",
           "SNPs_all.parsimony"  = "Cluster 3: all SNPs (parsimony)",
           "core_SNPs.ML"        = "Cluster 3: core SNPs (ML)",
           "SNPs_all.ML"         = "Cluster 3: all SNPs (ML)")

for (k in names(trees)) {
  f <- file.path(tree_dir, sprintf("tree_tipAlleleCounts.%s.tre", k))
  if (!file.exists(f)) { message("No encontrado: ", f); next }
  fig <- plot_tree(f, trees[[k]])
  for (ext in c("pdf", "svg", "png")) {
    ggsave(file.path(out_dir, sprintf("arbol_C3_%s.%s", sub("\\.", "_", k), ext)), fig,
           width = 16, height = 12, dpi = 300, bg = "white",
           device = if (ext == "pdf") cairo_pdf else NULL)
  }
  message("OK: ", k)
}
