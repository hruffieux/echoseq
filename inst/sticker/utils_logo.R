my_sticker <- function (subplot, s_x = 0.8, s_y = 0.75, s_width = 0.4, s_height = 0.5,
                        package, p_x = 1, p_y = 1.4, p_color = "#FFFFFF", p_family = "Aller_Rg",
                        p_fontface = "plain", p_size = 8, h_size = 1.2, h_fill = "#1881C2",
                        h_color = "#87B13F", spotlight = FALSE, l_x = 1, l_y = 0.5,
                        l_width = 3, l_height = 3, l_alpha = 0.4, url = "", u_x = 1,
                        u_y = 0.08, u_color = "black", u_family = "Aller_Rg", u_size = 1.5,
                        u_angle = 30, white_around_sticker = FALSE, ..., filename = paste0(package,
                                                                                           ".png"), asp = 1, dpi = 300)
{

  require(hexSticker)

  my_whiteTrans <- function(alpha = 0.4) {
    function(n) {
      rgb(red = rep(1, n), green = rep(1, n), blue = rep(1, n),
          alpha = seq(0, alpha, length.out = n))
    }
  }


  my_spotlight <- function(alpha) {
    ## set.seed(123)
    vals_x <- rnorm(500000, sd = 2, mean = 0)
    vals_y <- rnorm(500000, sd = 2, mean = 0)
    hexbinplot(vals_x ~ vals_y, colramp = my_whiteTrans(alpha), colorkey = FALSE,
               bty = "n", scales = list(draw = FALSE), xlab = "", ylab = "",
               border = NA, par.settings = list(axis.line = list(col = NA)))
  }


  hex <- ggplot() + geom_hexagon(size = h_size, fill = h_fill,
                                 color = NA)
  if (inherits(subplot, "character")) {
    d <- data.frame(x = s_x, y = s_y, image = subplot)
    sticker <- hex + geom_image(aes_(x = ~x, y = ~y, image = ~image),
                                d, size = s_width, asp = asp)
  }
  else {
    sticker <- hex + geom_subview(subview = subplot, x = s_x,
                                  y = s_y, width = s_width, height = s_height)
  }
  sticker <- sticker + geom_hexagon(size = h_size, fill = NA,
                                    color = h_color)
  if (spotlight)
    sticker <- sticker + geom_subview(subview = my_spotlight(l_alpha),
                                      x = l_x, y = l_y, width = l_width, height = l_height)
  sticker <- sticker + geom_pkgname(package, p_x, p_y, color = p_color,
                                    family = p_family, fontface = p_fontface, size = p_size,
                                    ...)
  sticker <- sticker + geom_url(url, x = u_x, y = u_y, color = u_color,
                                family = u_family, size = u_size, angle = u_angle)
  if (white_around_sticker)
    sticker <- sticker + white_around_hex(size = h_size)
  sticker <- sticker + theme_sticker(size = h_size) + geom_hexagon(size = h_size, fill = h_fill,
                                                                   color = NA)
  save_sticker(filename, sticker, dpi = dpi)
  class(sticker) <- c("sticker", class(sticker))
  invisible(sticker)
}
