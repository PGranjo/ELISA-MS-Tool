# ELISA-MS-Tool

The ELISA-MS tool was designed with Shiny R using the Alphalyse Shiny App template. Documentation on how to use the app can be found in [Wiki ELISA-MS tool](https://github.com/PGranjo/ELISA-MS-Tool/wiki)


## Installation

### Clone repository

```bash
git clone https://github.com/yourusername/yourrepository


packages <- c("shiny", "shinyjs", "DT", "markdown", "roxygen2",
              "shinyWidgets", "svglite", "readxl", "dplyr",
              "seqinr", "stringr", "writexl", "UpSetR", 
              "ggVennDiagram", "ggplot2", "gridExtra", "grid")

# Packages to install
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install) > 0) install.packages(to_install)

