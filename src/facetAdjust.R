library(grid)
# pos - where to add new labels newpage, vp - see ?print.ggplot
facetAdjust <- function(x, pos = c("up", "down")) {
    pos <- match.arg(pos)
    p <- ggplot_build(x)
    gtable <- ggplot_gtable(p)
    dev.off()
    dims <- apply(p$panel$layout[2:3], 2, max)
    nrow <- dims[1]
    ncol <- dims[2]
    panels <- sum(grepl("panel", names(gtable$grobs)))
    space <- ncol * nrow
    n <- space - panels
    if (panels != space) {
        idx <- (space - ncol - n + 1):(space - ncol)
        gtable$grobs[paste0("axis_b", idx)] <- list(gtable$grobs[[paste0("axis_b", 
            panels)]])
        if (pos == "down") {
            rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                gtable$layout$name)
            lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
            gtable$layout[rows, c("t", "b")] <- gtable$layout[lastAxis, 
                c("t")]
        }
    }
    class(gtable) <- c("facetAdjust", "gtable", "ggplot")
    gtable
}

print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
    if (newpage) 
        grid.newpage()
    if (is.null(vp)) {
        grid.draw(x)
    } else {
        if (is.character(vp)) 
            seekViewport(vp) else pushViewport(vp)
        grid.draw(x)
        upViewport()
    }
    invisible(x)
} 
