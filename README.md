## Github Repository for
# Transmission bottlenecks drive seed specialization and outbreak-like population structure in plant-associated bacteria
**by Irene Sanz-Puente1†, Santiago Redondo-Salvo1, Natalia I. García- Tomsig1, , Arancha Peñil-Celis1, Jorge Rodríguez Grande2,3, Gijs Selten4‡, Susana Fernandes 5,6, Andreas Börner7, Alain Ocampo-Sosa2,3, Esther Menendez5,6,7, Óscar Lorenzo5,6, Ronnie de Jonge4,9, Fernando de la Cruz 1*, and Marta Robledo1**

Repository associated with the analysis used in Sanz-Puente et al., 2026, focus on seed transmission as a powerful driver of ecological specialization of plant-associated bacteria. 
<i>This work is under review.</i> 

## Repository structure
- `Data/`: Contains all input data used in the project. 
- `Script/`: Contains all scripts (Bash and R) used for data processing and analysis. 
- `Results/`: Contains all output generated from the analysis, including figures and tables. 
 
## Data analysis

WGS of this study have been deposited in the NCBI Sequence Read Archive (RSA) under the BioProject accession number PRJNA1282304. Other data used in the analysis can be found in the Supplementary material of the paper or as data frame in scripts.

#### R 

R (>= 4.2) with the following packages:  
```r
adegenet         # Population genetics
ape              # Phylogenetic analyses
argparser        # Command-line argument parsing
biomformat       # BIOM format data handling
Biostrings       # Biological sequence analysis
broom            # Tidy statistical model outputs
BSDA             # Basic statistical analyses
circlize         # Circular visualizations
clusterProfiler  # Functional enrichment analysis
ComplexHeatmap   # Heatmaps
ComplexUpset     # UpSet plots
dplyr            # Data manipulation
emmeans          # Estimated marginal means and contrasts
factoextra       # Multivariate analysis visualization
FactoMineR       # Multivariate statistical analysis
forcats          # Categorical data handling
FSA              # Statistical analyses
geodata          # Geographic data
geosphere        # Geographic calculations
ggraph           # Graph visualization
ggnewscale       # Multiple fill and colour scales
ggplot2          # Data visualization
ggrepel          # Non-overlapping plot labels
ggsci            # Scientific colour palettes
ggspatial        # Spatial visualization and map scales
ggtext           # Text formatting in plots
ggtree           # Phylogenetic tree visualization
grid             # Graphics infrastructure
hierfstat        # Population genetics statistics
igraph           # Network analysis
KEGGREST         # KEGG database access
lme4             # Linear mixed-effects models
lmerTest         # Statistical tests for mixed-effects models
mapdata          # Map datasets
maps             # Geographic maps
pairwiseAdonis   # Pairwise PERMANOVA
patchwork        # Figure assembly
pegas            # Population and evolutionary genetics
pheatmap         # Heatmaps
phyloseq         # Microbiome analyses
purrr            # Functional programming
RColorBrewer     # Colour palettes
readr            # Data import
readxl           # Excel data import
reshape2         # Data reshaping
rstatix          # Statistical tests
scales           # Plot scaling utilities
scatterpie       # Pie charts on maps
sf               # Spatial data handling
stringr          # String manipulation
tibble           # Data frame utilities
tidygraph        # Tidy network analysis
tidyr            # Data tidying
tidyverse        # Data science framework
tools            # R tools and utilities
vegan            # Ecological analyses
writexl          # Excel data export
````

#### Citation
If you use this repository or its contents, please cite:
Sanz-Puente et al. (2026). Transmission bottlenecks drive seed specialization and outbreak-like population structure in plant-associated bacteria.
