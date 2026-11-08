# Publication tree for a personal website, built with tree_graph().
# papers.csv: theme | subcategory | title | year | role | journal | doi
# Only published papers; edit the csv to change what is shown.

source("tree_graph.R")

papers <- read.csv("papers.csv")
wrap <- function(x, w) vapply(strwrap(x, w, simplify = FALSE), paste,
                              character(1), collapse = "\n")
papers$label <- paste0(wrap(papers$title, 24), "\n",
                       papers$journal, " ", papers$year)

p <- papers |>
  tree_graph(parent = "theme", subcategory = "subcategory", term = "label",
             group = "role", wrap = NULL,
             group_colors = c("Joint first author" = "#CE284E",
                              "Co-author" = "#CDC7BC"))

ggsave("cv_graph.png", p, width = 12, height = 11, dpi = 200)
