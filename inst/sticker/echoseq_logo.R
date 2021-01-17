rm(list = ls())

require(dplyr)
require(ggpubr)

inst_dir <- "inst/sticker/"
source(file.path(inst_dir, "utils_logo.R"))

seed <- 123
set.seed(seed)

output_dir <- "man/figures/"
dir.create(output_dir, showWarnings = FALSE)

output_file <- "echoseq_logo.png"
col <- "grey30"

my_sticker(file.path(inst_dir, "echoseq_logo_base.png"),
           package="echoseq",
           p_size=5.8,
           p_color = col,
           s_x=0.8,
           s_y=1.1,
           s_width=1.7,
           s_height=1.5,
           p_x = 1.3,
           p_y = 0.7,
           h_size = 0.5,
           h_fill= grDevices::adjustcolor(col, alpha.f = 0.15),
           h_color= col,
           spotlight = TRUE,
           l_x = 1.1,
           l_y = 1.3,
           l_alpha = 0.4,
           white_around_sticker = FALSE,
           filename=file.path(output_dir, output_file),
           asp = 0.9,
           dpi = 1200)

# transparent edges
system(paste0("convert ",
              output_dir, output_file,
              " -transparent '#ffffff' ",
              output_dir, output_file))
