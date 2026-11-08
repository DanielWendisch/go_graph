# tree_graph(): the same radial tree for categorical data (no p-values needed).
#
#   source("tree_graph.R")
#   papers |> tree_graph(parent = "theme", subcategory = "subcategory",
#                        term = "label", group = "role")
#
# Input: one row per item with a parent, an optional subcategory (<= 3 items
# each), a label, and optionally
#   group : a category that colours the dots (e.g. "Joint first author")
#   size  : a number that sizes the dots
# Layout, twist and the no-crossing guarantee are shared with go_graph().

if (!exists("go_layout")) source("go_graph.R")

# five-colour dark theme (swap the values to re-skin)
tree_colors <- list(
  bg        = "#171B18",   # green-black background
  wedge     = "#344A3D",   # forest green parent slices
  node_fill = "#6B1328",   # deep ruby subcategory nodes ...
  node_ring = "#CE284E",   # ... with a bright ruby ring
  ink       = "#CDC7BC",   # warm chalk text
  ink_soft  = "#A09B93",   # dimmed chalk (subcategory labels)
  line      = "#8A978B"    # sage lines, readable on background and slices
)

#' @param data         data.frame, one row per item
#' @param parent,subcategory,term   column names (`subcategory` optional)
#' @param group        optional column; colours the dots by category
#' @param size         optional numeric column; sizes the dots
#' @param group_colors colours for the groups, named or in level order
#' @param colors       list of colours, see `tree_colors`
#' @param wrap         wrap term labels to this many characters (NULL = no wrap)
#' @param ...          further layout options as in go_graph(): twist,
#'                     max_per_subnode, radii, gaps, wedge_pad, label_pos
#' @return a ggplot object
tree_graph <- function(data,
                       parent      = "parent",
                       subcategory = "subcategory",
                       term        = "term",
                       group       = NULL,
                       size        = NULL,
                       group_colors = c("#CE284E", "#CDC7BC", "#8FB09B",
                                        "#E8A0AE", "#6B8F7B"),
                       colors      = tree_colors,
                       text_size   = 3.6,
                       wrap        = 24,
                       label_room  = 4.8,
                       twist       = 7,
                       max_per_subnode = 3,
                       radii       = c(parent = 4.8, subnode = 7.4, child = 9),
                       gaps        = c(parent = 1.3, sub = 0.5),
                       wedge_pad   = 0.35,
                       label_pos   = c(parent = 0.62, sub = 0.5)) {

  d <- go_prepare(data, parent, subcategory, term, p = NULL, count = size,
                  max_per_subnode, colour_fun = identity, keep = group)
  if (!is.null(wrap)) {
    d$term <- vapply(strwrap(d$term, wrap, simplify = FALSE),
                     paste, character(1), collapse = "\n")
  }
  d$group <- if (is.null(group)) "all" else factor(d[[group]], levels = unique(d[[group]]))
  tree <- go_layout(d, radii, gaps, twist, agg = mean, wedge_pad, label_pos)

  terms <- tree$terms; subs <- tree$subs; pars <- tree$pars
  pars$halo <- colors$wedge                        # labels sit on the slices
  subs$halo <- colors$wedge
  lim <- radii[["child"]] + label_room

  p <- ggplot() +
    geom_polygon(data = tree$wedges, fill = colors$wedge, colour = NA,
                 aes(.data$x, .data$y, group = .data$id)) +
    geom_path(data = tree$edges, colour = colors$line, linewidth = 0.75,
              lineend = "round", aes(.data$x, .data$y, group = .data$id)) +
    halo_text(pars, size = text_size - 0.1, colour = colors$ink,
              fontface = "bold") +
    halo_text(subs, size = text_size - 0.5, colour = colors$ink_soft) +
    geom_point(data = subs, shape = 21, size = 2.8, stroke = 0.6,
               fill = colors$node_fill, colour = colors$node_ring,
               aes(.data$x, .data$y)) +
    geom_text(data = terms, size = text_size, colour = colors$ink, vjust = 0.5,
              lineheight = 0.9,
              aes(.data$lx, .data$ly, label = .data$label,
                  angle = .data$lab_angle, hjust = .data$lab_hjust)) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.title = element_text(colour = colors$ink, size = 11),
          legend.text  = element_text(colour = colors$ink, size = 10),
          plot.margin  = margin(10, 10, 10, 10),
          plot.background = element_rect(fill = colors$bg, colour = NA))

  # dots: coloured by group (if given), sized by `size` (if given)
  dot <- if (is.null(size)) aes(.data$x, .data$y, colour = .data$group)
         else aes(.data$x, .data$y, colour = .data$group, size = .data$n)
  p <- p + geom_point(data = terms, mapping = dot,
                      size = if (is.null(size)) 3.6 else NULL)
  if (is.null(group)) {
    p <- p + scale_colour_manual(values = group_colors[1], guide = "none")
  } else {
    cols <- if (is.null(names(group_colors))) {
      setNames(rep_len(group_colors, nlevels(terms$group)), levels(terms$group))
    } else group_colors
    p <- p + scale_colour_manual(values = cols, name = NULL) +
      guides(colour = guide_legend(order = 1))
  }
  if (!is.null(size)) p <- p + scale_size_area(max_size = 6, name = size)
  p
}
