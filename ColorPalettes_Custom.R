# CUSTOM COLOR PALETTES #
#-----------------------#

# Set colours
custom_colours = list(
  darkRainbow = c("#9F248FFF", "#C6242DFF", "#F9791EFF", "#FFCE4EFF", "#017A4AFF", "#244579FF")
)

# Make palettes
custom_palettes = function(name, n, all_palettes = custom_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,
               continuous = grDevices::colorRampPalette(palette)(n),
               discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

# Make ggplot functions
scale_colour_multi_d = function(name) {
  ggplot2::scale_colour_manual(values = custom_palettes(name,
                                                       type = "discrete"))
}
scale_fill_multi_d = function(name) {
  ggplot2::scale_fill_manual(values = custom_palettes(name,
                                                     type = "discrete"))
}
scale_colour_multi_c = function(name) {
  ggplot2::scale_colour_gradientn(colours = custom_palettes(name = name,
                                                           type = "continuous"))
}
scale_fill_multi_c = function(name) {
  ggplot2::scale_fill_gradientn(colours = custom_palettes(name = name,
                                                         type = "continuous"))
}
scale_color_multi_d = scale_colour_multi_d
scale_color_multi_c = scale_colour_multi_c


