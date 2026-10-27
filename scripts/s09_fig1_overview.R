suppressPackageStartupMessages(library(grid))
out <- "results/figures/Fig1_overview.png"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

## ---- palette ---------------------------------------------------------------
LAYER <- list(
  L1 = list(fill = "#E3EEFB", border = "#3F6FA5", name = "Discovery"),
  L2 = list(fill = "#DCF1EC", border = "#2E8B82", name = "Trajectory"),
  L3 = list(fill = "#E2F0DD", border = "#3F8F3F", name = "In silico"),
  L4 = list(fill = "#FCE7D6", border = "#C2772E", name = "Consensus & validation")
)

## ---- boxes -----------------------------------------------------------------
boxes <- list(
  list(id="B1", layer="L1", x=3.6, y=10.5, w=4.4, h=1.3, badge=1, code="",
       title="Bone-marrow stromal scRNA-seq",
       sub="normal donors & non-union patients (GSE337700)"),
  list(id="B2", layer="L1", x=3.6, y=8.5,  w=4.8, h=1.3, badge=2, code="S02",
       title="Cell atlas & heterogeneity",
       sub="Seurat v5 + Harmony  |  5 stromal subpopulations"),
  list(id="B3", layer="L2", x=3.6, y=6.5,  w=5.0, h=1.3, badge=3, code="S03",
       title="Osteogenic pseudotime trajectory",
       sub="Monocle2 DDRTree  |  root  ->  osteoblast fate"),
  list(id="B4", layer="L3", x=3.6, y=4.5,  w=5.0, h=1.4, badge=4, code="S04",
       title="In silico virtual knockout",
       sub="scTenifoldKnk  |  49 druggable osteogenic brakes"),
  list(id="B5", layer="L3", x=7.6, y=4.5,  w=4.2, h=1.4, badge=5, code="S05",
       title="Compound-target network",
       sub="Ciwujia + Wuweizi (ChEMBL)  |  196 targets"),
  list(id="B6", layer="L4", x=5.0, y=2.5,  w=5.4, h=1.5, badge=6, code="S06",
       title="Three-way consensus",
       sub="KO x trajectory x compound  ->  Tier-1 triad"),
  list(id="B7", layer="L4", x=5.0, y=0.8,  w=5.4, h=1.4, badge=7, code="S07",
       title="OA external validation",
       sub="Tier-1 brakes up-regulated in OA stroma")
)
bands <- list(
  list(layer="L1", x0=0.2, x1=9.8, y0=7.7,  y1=11.3),
  list(layer="L2", x0=0.2, x1=9.8, y0=5.6,  y1=7.4),
  list(layer="L3", x0=0.2, x1=9.8, y0=3.6,  y1=5.4),
  list(layer="L4", x0=0.2, x1=9.8, y0=0.0,  y1=3.45)
)
## segments: c(x1, y1, x2, y2, arrow?)
segs <- list(
  c(3.6, 9.85,  3.6, 9.15, TRUE),   # B1 -> B2
  c(3.6, 7.85,  3.6, 7.15, TRUE),   # B2 -> B3
  c(3.6, 5.85,  3.6, 5.20, TRUE),   # B3 -> B4
  c(3.6, 3.80,  3.6, 3.52, FALSE),  # B4 down to merge line
  c(7.6, 3.80,  7.6, 3.52, FALSE),  # B5 down to merge line
  c(3.6, 3.52,  7.6, 3.52, FALSE),  # merge horizontal
  c(5.0, 3.52,  5.0, 3.25, TRUE),   # merged arrow into B6 top
  c(5.0, 1.75,  5.0, 1.50, TRUE)    # B6 -> B7
)

draw_all <- function() {
  N <- function(v) unit(v, "native")
  pushViewport(viewport(xscale = c(0, 10), yscale = c(0, 11.6)))

  ## header line
  grid.text("Integrative four-layer design linking single-cell discovery to compound-target therapy",
            x = N(5.0), y = N(11.5), just = "centre",
            gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 11.5, col = "#222222"))

  ## layer bands (drawn first, everything else on top)
  for (b in bands) {
    L <- LAYER[[b$layer]]
    grid.rect(x = N(b$x0), y = N(b$y0), width = N(b$x1 - b$x0), height = N(b$y1 - b$y0),
              just = c("left", "bottom"), gp = gpar(fill = L$fill, col = NA))
    yc <- (b$y0 + b$y1) / 2
    grid.text(b$layer, x = N(0.75), y = N(yc), just = "centre",
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 17, col = L$border))
    grid.text(L$name, x = N(0.75), y = N(yc - 0.62), just = "centre",
              gp = gpar(fontfamily = "sans", fontsize = 8, col = "#555555"))
  }

  ## connector arrows (before boxes so heads sit on box borders)
  for (s in segs) {
    grid.lines(x = unit(c(s[1], s[3]), "native"), y = unit(c(s[2], s[4]), "native"),
               arrow = if (s[5]) arrow(angle = 22, length = unit(0.20, "cm"),
                                       type = "closed", ends = "last") else NULL,
               gp = gpar(col = "#3a3a3a", lwd = 2.0, lineend = "round"))
  }

  ## boxes
  for (bx in boxes) {
    L <- LAYER[[bx$layer]]
    x0 <- bx$x - bx$w / 2
    ## soft shadow
    grid.roundrect(x = N(bx$x + 0.12), y = N(bx$y - 0.12), width = N(bx$w), height = N(bx$h),
                   r = unit(0.08, "in"), just = "centre",
                   gp = gpar(fill = "#000000", col = NA, alpha = 0.07))
    ## main rounded box
    grid.roundrect(x = N(bx$x), y = N(bx$y), width = N(bx$w), height = N(bx$h),
                   r = unit(0.08, "in"), just = "centre",
                   gp = gpar(fill = L$fill, col = L$border, lwd = 2.0))
    ## numbered badge
    cx <- x0 + 0.34; cy <- bx$y + bx$h / 2 - 0.34
    grid.circle(x = N(cx), y = N(cy), r = unit(0.13, "in"),
                gp = gpar(fill = L$border, col = "white", lwd = 1.6))
    grid.text(bx$badge, x = N(cx), y = N(cy), just = "centre",
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 11, col = "white"))
    ## title + subtitle
    tx <- x0 + 0.78
    head <- if (nzchar(bx$code)) paste0(bx$code, "  ") else ""
    grid.text(paste0(head, bx$title), x = N(tx), y = N(bx$y + bx$h / 2 - 0.36),
              just = c("left", "centre"),
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 10.5, col = "#1a1a1a"))
    grid.text(bx$sub, x = N(tx), y = N(bx$y + bx$h / 2 - 0.82),
              just = c("left", "centre"),
              gp = gpar(fontfamily = "sans", fontsize = 8, col = "#4a4a4a"))
  }
  popViewport()
}

## ---- render PNG + PDF ------------------------------------------------------
png(out, width = 7.2, height = 8.3, units = "in", res = 300, bg = "white"); draw_all(); dev.off()
pdf(sub("\\.png$", ".pdf", out), width = 7.2, height = 8.3, bg = "white"); draw_all(); dev.off()
cat("Saved", out, "and PDF\n")
