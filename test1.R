library(volesti)
library(rgl)


visualize_sampling_2d <- function(n = 1000, step = 20) {
  p <- gen_cube(dimension = 2, representation = "H")
  walk_params <- list("walk" = "BaW", "walk_length" = step)
  samples <- sample_points(P = p, n = n, random_walk = walk_params)

  quartz(title = "2d plot") # disable if not on macOS
  plot(samples[1, ], samples[2, ], col = "blue", pch = 20,
       cex = 0.6,
       main = "",
       xlab = "X-axis", ylab = "Y-axis",
       xlim = c(-1, 1), ylim = c(-1, 1))
  abline(h = 0, v = 0, col = "black", lwd = 2)
  grid()
}


visualize_sampling_3d <- function(n = 1000, step = 20) {
  p <- gen_cube(dimension = 3, representation = "H")
  walk_params <- list("walk" = "BaW", "walk_length" = step)
  samples <- sample_points(P = p, n = n, random_walk = walk_params)

  options(rgl.printRglWidget = FALSE)
  open3d()
  plot3d(samples[1, ], samples[2, ], samples[3, ],
         col = "purple", size = 1.0, type = "s",
         main = "", xlab = "x", ylab = "y", zlab = "z",
         xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1))
  print(rglwidget())
}


visualize_sampling_2d(n = 1500, step = 20)
visualize_sampling_3d(n = 1500, step = 20)
