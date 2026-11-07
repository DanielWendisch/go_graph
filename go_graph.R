# go_graph(): radial GO-enrichment tree as a ggplot.
#
#   library(ggplot2); source("go_graph.R")
#   go |> go_graph()                                  # defaults
#   go |> go_graph(twist = 10, text_size = 3) +       # returns a ggplot, so
#     labs(title = "My enrichment")                   # you can keep adding
#   go |> go_graph() |> ggsave("go.png", plot = _, width = 14, height = 12)
#
# Input: one row per child term, with columns for
#   parent term | subcategory (optional, <= 3 terms each) | term | p | count
#
# Reading the plot
#   * every child term is a dot: colour = -log10(p), size = count
#   * every parent is a twisted background slice, coloured by the aggregate
#     (default: mean) of its children
#   * subcategories are small grey nodes between parent and children
#
# Why no lines cross: each child term owns its own angular slot, ordered parent
# by parent and subcategory by subcategory. A subcategory sits at the mean angle
# of its terms and a parent at the mean angle of its subcategories, so every
# node's wedge lies inside its parent's wedge and sibling wedges never overlap.
# The twist is a smooth function of radius applied on top, so it can't create
# crossings either.

library(ggplot2)

go_palette <- c("#a1dab4", "#5ec4b6", "#2c9fc3", "#2f6db5", "#3b3a9e", "#2a1a63")

#' @param data            data.frame, one row per child term
#' @param parent,subcategory,term,p,count   column names in `data`.
#'                        `subcategory` is optional: if the column is missing,
#'                        terms are split into "Subcategory 1, 2, ..." by p-value.
#' @param twist           degrees each ring is rotated further than the previous
#'                        (0 = straight lines, negative = other direction)
#' @param max_per_subnode largest allowed subcategory (only used for auto-split)
#' @param agg             how a parent's value is aggregated from its children
#' @param palette         colours for the -log10(p) scale, low -> high
#' @param wedge_alpha     opacity of the parent slices (0 = hidden, 1 = solid)
#' @param text_size       size of the term labels; other labels scale from it
#' @param radii           ring radii: parent label end, subcategory, term
#' @param gaps            empty angular slots between parents / subcategories
#' @return a ggplot object
go_graph <- function(data,
                     parent      = "parent_term",
                     subcategory = "subcategory",
                     term        = "Description",
                     p           = "p.adjust",
                     count       = "Count",
                     twist       = 7,
                     max_per_subnode = 3,
                     agg         = mean,
                     palette     = go_palette,
                     wedge_alpha = 0.38,
                     text_size   = 3.4,
                     radii       = c(parent = 3.8, subnode = 5.9, child = 7.1),
                     gaps        = c(parent = 1, sub = 0.4)) {

  d    <- go_prepare(data, parent, subcategory, term, p, count, max_per_subnode)
  tree <- go_layout(d, radii, gaps, twist, agg)
  go_draw(tree, palette, wedge_alpha, text_size)
}

# ---- 1. validate + tidy the input --------------------------------------------
go_prepare <- function(data, parent, subcategory, term, p, count,
                       max_per_subnode) {
  missing_cols <- setdiff(c(parent, term, p, count), names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "),
         "\n  Available: ", paste(names(data), collapse = ", "), call. = FALSE)
  }

  d <- data.frame(parent = as.character(data[[parent]]),
                  term   = as.character(data[[term]]),
                  p      = as.numeric(data[[p]]),
                  n      = as.numeric(data[[count]]),
                  stringsAsFactors = FALSE)
  d$subcat <- if (subcategory %in% names(data)) as.character(data[[subcategory]])
              else NA_character_

  if (anyNA(d[c("parent", "term", "p", "n")])) {
    stop("`", parent, "`, `", term, "`, `", p, "` and `", count,
         "` must not contain NA.", call. = FALSE)
  }
  if (any(d$p <= 0)) {                      # p = 0 would give -log10(p) = Inf
    d$p[d$p <= 0] <- min(d$p[d$p > 0]) / 10
    message("go_graph: p-values of 0 were set to ", format(min(d$p) , digits = 2))
  }
  d$logp   <- -log10(d$p)
  d$parent <- factor(d$parent, levels = unique(d$parent))   # keep input order

  if (all(is.na(d$subcat))) {
    # no names given: split each parent's terms into chunks of similar size
    d <- d[order(d$parent, d$p), ]
    d$subcat <- ave(seq_len(nrow(d)), d$parent, FUN = function(i) {
      k <- length(i)
      paste("Subcategory", sort(rep(seq_len(ceiling(k / max_per_subnode)),
                                    length.out = k)))
    })
  } else {
    d$subcat <- factor(d$subcat, levels = unique(d$subcat))
    d <- d[order(d$parent, d$subcat, d$p), ]
    d$subcat <- as.character(d$subcat)
    size <- table(paste(d$parent, d$subcat, sep = " > "))
    if (any(size > max_per_subnode)) {
      stop("Subcategories with more than ", max_per_subnode, " terms: ",
           paste(names(size)[size > max_per_subnode], collapse = "; "),
           call. = FALSE)
    }
  }
  d$sub_id <- paste(d$parent, d$subcat, sep = "#")
  rownames(d) <- NULL
  d
}

# ---- 2. layout: positions, edges, labels, background slices -------------------
go_layout <- function(d, radii, gaps, twist, agg, wedge_pad = 0.35,
                      label_pos = c(parent = 0.58, sub = 0.6)) {
  n <- nrow(d)

  # one angular slot per term, small gap between subcategories, bigger between
  # parents
  new_parent <- c(FALSE, d$parent[-1] != d$parent[-n])
  new_sub    <- c(FALSE, d$sub_id[-1] != d$sub_id[-n]) & !new_parent
  slot  <- cumsum(1 + gaps[["parent"]] * new_parent + gaps[["sub"]] * new_sub)
  slice <- 2 * pi / (max(slot) + gaps[["parent"]])      # angle of one slot
  d$theta <- (slot - 1) * slice

  # inner nodes sit at the mean angle of what they contain
  subs <- unique(d[c("sub_id", "parent", "subcat")])
  subs$theta <- as.numeric(tapply(d$theta, d$sub_id, mean)[subs$sub_id])
  pars <- data.frame(parent = levels(d$parent))
  pars$theta <- as.numeric(tapply(subs$theta, subs$parent, mean)[pars$parent])
  pars$logp  <- as.numeric(tapply(d$logp, d$parent, agg)[pars$parent])

  # twist: ring i is rotated by i * twist degrees; a continuous function of
  # radius, so lines become smooth spirals
  tw <- approxfun(c(0, radii), (0:3) * twist * pi / 180, rule = 2)
  place <- function(r, theta) {
    phi <- theta + tw(r)
    data.frame(x = r * cos(phi), y = r * sin(phi), phi = phi)
  }
  d    <- cbind(d,    place(radii[["child"]],   d$theta))
  subs <- cbind(subs, place(radii[["subnode"]], subs$theta))
  pars <- cbind(pars, place(radii[["parent"]],  pars$theta))

  # an edge is straight in (radius, angle) space and twisted afterwards;
  # `label` is a point on it (at fraction `frac`) with the local text angle
  edge <- function(r0, th0, r1, th1, frac = 0.5, k = 40) {
    t <- seq(0, 1, length.out = k)
    r <- r0 + t * (r1 - r0)
    a <- th0 + t * (th1 - th0) + tw(r)
    x <- r * cos(a); y <- r * sin(a)
    i <- round(frac * (k - 1)) + 1
    ang <- atan2(y[i + 1] - y[i - 1], x[i + 1] - x[i - 1])
    if (cos(ang) < 0) ang <- ang + pi                     # keep text upright
    list(path  = data.frame(x = x, y = y),
         label = data.frame(lx = x[i], ly = y[i], lab_angle = ang * 180 / pi))
  }
  e_par  <- Map(edge, 0, pars$theta, radii[["parent"]], pars$theta,
                label_pos[["parent"]])
  e_sub  <- Map(edge, radii[["parent"]],
                pars$theta[match(subs$parent, pars$parent)],
                radii[["subnode"]], subs$theta, label_pos[["sub"]])
  e_leaf <- Map(edge, radii[["subnode"]],
                subs$theta[match(d$sub_id, subs$sub_id)],
                radii[["child"]], d$theta)

  all_edges <- c(e_par, e_sub, e_leaf)
  edges <- do.call(rbind, Map(function(e, id) cbind(e$path, id = id),
                              all_edges, seq_along(all_edges)))
  labels_of <- function(es) do.call(rbind, lapply(es, `[[`, "label"))
  pars <- cbind(pars, labels_of(e_par));  pars$label <- pars$parent
  subs <- cbind(subs, labels_of(e_sub));  subs$label <- subs$subcat

  # term labels: radial, flipped on the left half so they stay readable
  left <- cos(d$phi) < 0
  d$lab_angle <- (d$phi * 180 / pi) + ifelse(left, 180, 0)
  d$lab_hjust <- ifelse(left, 1, 0)
  d$lx <- (radii[["child"]] + 0.3) * cos(d$phi)
  d$ly <- (radii[["child"]] + 0.3) * sin(d$phi)
  d$label <- d$term

  # background slice per parent: root -> subcategory ring, spiral sides
  wedge <- function(i, k = 40) {
    th <- d$theta[d$parent == pars$parent[i]]
    lo <- min(th) - wedge_pad * slice
    hi <- max(th) + wedge_pad * slice
    R  <- radii[["subnode"]]
    r  <- seq(0, R, length.out = k)
    a  <- c(lo + tw(r), seq(lo, hi, length.out = k) + tw(R), rev(hi + tw(r)))
    rad <- c(r, rep(R, k), rev(r))
    data.frame(x = rad * cos(a), y = rad * sin(a), id = i, logp = pars$logp[i])
  }
  wedges <- do.call(rbind, lapply(seq_len(nrow(pars)), wedge))

  list(terms = d, subs = subs, pars = pars, edges = edges, wedges = wedges,
       radii = radii)
}

# ---- 3. draw ------------------------------------------------------------------
go_draw <- function(tree, palette, wedge_alpha, text_size) {
  terms <- tree$terms; subs <- tree$subs; pars <- tree$pars
  ink <- "#1f2a37"; sub_ink <- "#475569"; muted <- "#8d99a6"
  lims <- range(c(terms$logp, pars$logp))

  # a label's halo has the colour of the slice underneath it
  tint <- function(v) {
    rgb_val <- col2rgb(scales::col_numeric(palette, lims)(v))
    rgb(t(wedge_alpha * rgb_val + (1 - wedge_alpha) * 255), maxColorValue = 255)
  }
  pars$halo <- tint(pars$logp)
  subs$halo <- pars$halo[match(subs$parent, pars$parent)]

  # text with a halo: copies nudged around a circle sit under the real text, so
  # the line behind it looks interrupted by the label
  halo_text <- function(data, size, colour, fontface = "plain", r = 0.09) {
    ang <- seq(0, 2 * pi, length.out = 13)[-13]
    halos <- lapply(ang, function(a) {
      dx <- r * cos(a); dy <- r * sin(a)
      geom_text(data = data, size = size, fontface = fontface,
                aes(.data$lx + dx, .data$ly + dy, label = .data$label,
                    angle = .data$lab_angle, colour = I(.data$halo)))
    })
    c(halos, list(geom_text(data = data, size = size, fontface = fontface,
                            colour = colour,
                            aes(.data$lx, .data$ly, label = .data$label,
                                angle = .data$lab_angle))))
  }

  lim <- tree$radii[["child"]] + 3.4

  ggplot() +
    geom_polygon(data = tree$wedges, alpha = wedge_alpha, colour = NA,
                 aes(.data$x, .data$y, group = .data$id, fill = .data$logp)) +
    geom_path(data = tree$edges, colour = muted, linewidth = 0.75,
              lineend = "round", aes(.data$x, .data$y, group = .data$id)) +
    halo_text(pars, size = text_size - 0.1, colour = ink, fontface = "bold") +
    halo_text(subs, size = text_size - 0.5, colour = sub_ink) +
    geom_point(data = subs, colour = muted, size = 2.6,
               aes(.data$x, .data$y)) +
    geom_point(data = terms,
               aes(.data$x, .data$y, colour = .data$logp, size = .data$n)) +
    geom_text(data = terms, size = text_size, colour = ink, vjust = 0.5,
              aes(.data$lx, .data$ly, label = .data$label,
                  angle = .data$lab_angle, hjust = .data$lab_hjust)) +
    scale_colour_gradientn(colours = palette, limits = lims,
                           name = expression(-log[10]~"(adj. p)")) +
    scale_fill_gradientn(colours = palette, limits = lims) +
    scale_size_area(max_size = 6, name = "Gene count") +
    guides(fill = "none",
           colour = guide_colourbar(order = 1, barheight = 8),
           size = guide_legend(order = 2)) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.title = element_text(colour = ink, size = 11),
          legend.text  = element_text(colour = ink, size = 10),
          plot.margin  = margin(10, 10, 10, 10),
          plot.background = element_rect(fill = "white", colour = NA))
}
