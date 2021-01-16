require(hexSticker)
require(dplyr)
require(ggpubr)


seed <- 123
set.seed(seed)

base_logo <- "inst/sticker/echoseq_logo_base.png"

dir.create("man/figures/", showWarnings = FALSE)

col <- "grey30"#"black" #"#5A6B5E" #"#266040"#"grey30"

sticker(base_logo,
        package="echoseq",
        p_size=5.8,
        p_color = col,
        s_x=0.8,
        s_y=1.1,
        s_width=1.7,
        s_height=1.5,
        p_x = 1.3,
        p_y = 0.7,
        h_size = 1.3,
        h_fill= col,
        h_color= col,
        spotlight = TRUE,
        l_x = 1.1,
        l_y = 1.3,
        l_alpha = 0.6,
        filename="man/figures/echoseq_logo.png",
        asp = 1,
        dpi = 1200)
