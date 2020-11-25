#!/usr/bin/env R

# Custom ggplot2 themes to source in other scripts

# # libraries for testing ----
# library(tidyverse)
# library(viridis)
# library(RColorBrewer)
# library(extrafont)

# define themes ----
# Dark theme
myggthemedark = theme_classic() +
	theme(
		rect = element_rect(fill = "#111122"),
		panel.background = element_rect(fill = "#111122"),
		axis.line = element_line(color = "white"),
		axis.ticks = element_line(color = "white"),
		panel.grid = element_line(color = "#333333"),
		panel.grid.major = element_line(size = rel(0.5)),
		panel.grid.minor = element_line(size = rel(0.25)),
		text = element_text(family = "CMU Bright SemiBold", color = "white"),
		axis.text = element_text(family = "SF Mono", color = "white"),
		legend.text = element_text(family = "SF Mono"),
		legend.position = "bottom",
		strip.background = element_rect(fill = "#111122", color = "white"),
		strip.text = element_text(family = "CMU Bright SemiBold", color = "white")
	)

# Light theme
myggthemelight = theme_light() +
	theme(
		# make plot background an off-white
		rect = element_rect(fill = "ghostwhite"),
		panel.background = element_rect(fill = "ghostwhite"),
		strip.background = element_rect(fill = "ghostwhite",
										color = "black"),
		# assign portable texts
		strip.text = element_text(family = "Helvetica"),
		text = element_text(family = "Helvetica"),
		axis.text = element_text(family = "Courier"),
		legend.text = element_text(family = "Courier")
	)

# # test themes ----
# # boxplot
# mtcars %>%
# 	ggplot() +
# 	aes(x = factor(cyl), y = mpg, fill = factor(cyl)) +
# 	geom_boxplot() +
# 	# scale_fill_brewer(palette = "Set1") +
# 	scale_fill_brewer(palette = "Pastel1") +
# 	# scale_fill_viridis(discrete = T) +
# 	geom_jitter(width = 0.2, shape = 21, color = "black", fill = "white") +
# 	# theme_light() +
# 	myggthemelight +
# 	# theme(legend.position = "none") +
# 	labs(title = "Fuel Efficiency vs. Engine Cylinders",
# 		 x = "Cylinders",
# 		 y = "MPG",
# 		 fill = "Cylinders")
#
# # scatterplot
# mtcars %>%
# 	ggplot() +
# 	aes(x = hp, y = mpg, color = factor(cyl)) +
# 	geom_point() +
# 	geom_smooth(method = "lm",
# 				size = 0.2,
# 				alpha = 0.2,
# 				se = F) +
# 	# geom_smooth(aes(fill = factor(cyl)),
# 	# 			method = "lm",
# 	# 			size = 0.2,
# 	# 			alpha = 0.2) +
# 	geom_smooth(method = "lm",
# 				size = 0.2,
# 				alpha = 0.2,
# 				color = "black") +
# 	# scale_color_brewer(palette = "Accent") +
# 	# scale_color_brewer(palette = "Dark2") +
# 	# scale_fill_brewer(palette = "Dark2") +
# 	# scale_color_brewer(palette = "Set2") +
# 	# scale_fill_brewer(palette = "Set2") +
# 	scale_color_brewer(palette = "Set1") +
# 	scale_fill_brewer(palette = "Set1") +
# 	theme_light() +
# 	labs(title = "Fuel Efficiency vs. Horsepower",
# 		 x = "HP",
# 		 y = "MPG",
# 		 color = "Cylinders")
