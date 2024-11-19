# Preamble ----------------------------------------------------------------

# Dependencies
if (!require('pacman')) install.packages('pacman')
pacman::p_load(
      tidyverse,
      hexSticker,
      ggplot2
    )

p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point()
p <- p + theme_void() + theme_transparent()
sticker(subplot = "./figures/frog.png",
        s_x	= 1,
        s_y	= 0.7,
        s_width = 117/128 * 0.4,
        s_height = 0.4,
        package = "auto_dag",
        p_size = 24,
        p_color = "cornsilk",
        h_fill = "forestgreen",
        h_color = "cornsilk",
        spotlight = TRUE,
        l_x	= 1,
        l_y = 1.25,
        l_width = 6,
        filename = "./figures/badge.png")
