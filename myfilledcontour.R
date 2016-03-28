my.filled.contour <-
  function (x = seq(0, 1, length.out = nrow(z)), 
            y = seq(0, 1,length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes, ...)
  {
    # filled.contour gives unnecessary legend
    # Used P Lapointe's solution from here:http://stackoverflow.com/questions/16774928/removing-part-of-a-graphic-in-r
    #   also had to changed .Internal(fillcontour) to .filled.contour
    #   and change layout to layout(matrix(c(1, 1), ncol = 1L), widths = c(1, lcm(w)))
    # Created 3/28/16 by Collin Erickson
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    #mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    #on.exit(par(par.orig))
    # Comment out everything that was for legend
    #w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(1, 1), ncol = 1L), widths = c(1, lcm(w)))
    #par(las = las)
    #mar <- mar.orig
    #mar[4L] <- mar[2L]
    #mar[2L] <- 1
    #par(mar = mar)
    #plot.new()
    #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
    #            yaxs = "i")
    #    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    #    if (missing(key.axes)) {
    #        if (axes)
    #            axis(4)
    #    }
    #    else key.axes
    #    box()
    #if (!missing(key.title))
    #  key.title
    #mar <- mar.orig
    #mar[4L] <- 1
    #par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    #.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels),
    #                        col = col))
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot)
      box()
    if (missing(plot.title))
      title(...)
    else plot.title
    invisible()
}

my.filled.contour.func <- function(fn,n=100) {
  x <- seq(0,1,length.out = n)
  y <- seq(0,1,length.out = n)
  z <- matrix(NA,n,n)
  for(xi in 1:n) for(yi in 1:n) z[xi,yi] <- fn(c(x[xi],y[yi]))
  my.filled.contour(x,y,z)
}