source("mock_data.R")
source("go_graph.R")

# the usual way: pipe a GO table in, get a ggplot out
make_mock_go() |>
  go_graph() |>
  ggsave("go_graph.png", plot = _, width = 14, height = 12, dpi = 150)

# tweak it, then keep adding ggplot layers
p <- make_mock_go() |>
  go_graph(twist = 12, wedge_alpha = 0.25) +
  labs(title = "GO enrichment")
