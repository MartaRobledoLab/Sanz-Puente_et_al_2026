## Github Repository for
# Seed transmission bottlenecks drive the emergence of a seed-specific bacterial lineage in wheat
**by Irene Sanz-Puente, Santiago Redondo-Salvo, Natalia I. García- Tomsig, Gijs Selten, Arancha Peñil-Celis, Jorge Rodríguez Grande, Alain Ocampo-Sosa, Esther Menendez, Óscar Lorenzo, Ronnie de Jonge, Fernando de la Cruz, and Marta Robledo.**

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
ape             # Phylogenetic analyses
ComplexHeatmap  # Heatmaps
dplyr           # Data manipulation
FSA             # Statistical analyses
ggplot2         # Visualization
ggtext          # Text formatting in plots
patchwork       # Figure assembly
phyloseq        # Microbiome analyses
pheatmap        # Heatmaps
purrr           # Functional programming
RColorBrewer    # Colour palettes
readr           # Data import
reshape2        # Data reshaping
rstatix         # Statistical tests
scales          # Plot scaling utilities
tibble          # Data frame utilities
tidyr           # Data tidying
tidyverse       # Data science framework
vegan           # Ecological analyses
````

#### Citation
If you use this repository or its contents, please cite:
Sanz-Puente et al. (2026). Seed transmission bottlenecks drive the emergence of a seed-specific bacterial lineage in wheat.
