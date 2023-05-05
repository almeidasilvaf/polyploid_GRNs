
# Load packages
library(ggplot2)
library(patchwork)
library(here)
library(ggpubr)
library(ComplexHeatmap)

library(ggtree)
library(colorspace)
library(tidyverse)

# Main figures ====
## Figure 1 ====

# This figure was created manually using Inkscape

## Figure 2 ====
### Load data
load(here(
    "products", "plots",
    "p_sequence_divergence_duplicates.rda"
))

fig2 <- p_sequence_divergence_duplicates

### Export figure
ggsave(
    fig2,
    file = here("products", "figs", "fig2.pdf"),
    width = 10, height = 6
)

ggsave(
    fig2,
    file = here("products", "figs", "fig2.png"),
    width = 10, height = 6, dpi = 300
)


## Figure 3 ====
### Load data
load(here("products", "plots", "p_sorensendice_ppi.rda"))

fig3 <- p_sorensendice_ppi

### Export figure
ggsave(
    fig3,
    file = here("products", "figs", "fig3.png"),
    width = 10, height = 9, dpi = 300
)

## Figure 4 ====
### Load data
load(here("products", "plots", "p_heatmap_z_scores.rda"))
load(here("products", "plots", "p_mean_z_scores.rda"))


### Plot species tree with WGD on branches
get_wgd_coords <- function(peaks_and_nodes, tree_plot) {

    plot_coord <- ggplot_build(tree_plot)$data
    coord_num <- plot_coord[[1]][, c("y", "node")]
    coord_char <- plot_coord[[3]][, c("y", "label")]
    names(coord_char) <- c("y", "node")

    all_coords <- rbind(coord_num, coord_char)

    nodes <- peaks_and_nodes$node
    wgds <- Reduce(rbind, lapply(nodes, function(x) {
        yaxis <- all_coords[all_coords$node == x, "y"]
        peak_info <- peaks_and_nodes[peaks_and_nodes$node == x, ]

        wgd_df <- data.frame(
            xmax = -peak_info$min,
            xmin = -peak_info$max,
            ymin = yaxis - 0.2,
            ymax = yaxis + 0.2,
            id = 1
        )
        return(wgd_df)
    }))
    return(wgds)
}

species_tree <- read.tree(here("data", "species_tree.nwk"))

periods <- data.frame(
    xmin = c(-66, -145, -201, -252, -299, -541),
    xmax = c(-2, -66, -145, -201, -252, -299),
    period = c("paleogene_quaternary", "cretaceous", "jurassic", "triassic",
               "permian", "cambrian_carboniferous"),
    id = 1:6,
    ymin = -Inf,
    ymax = Inf
)
period_col <- c(
    "#f3f4f0", "#dbdfd3", "#c8cebc", "#c3c9b5", "#b7bfa6", "#abb497"
)


mes <- data.frame(
    xmin = c(-75, -260),
    xmax = c(-55, -240),
    period = c("me1", "me2"),
    id = 1,
    ymin = -Inf,
    ymax = Inf
)

p1 <- ggtree(species_tree, color = "dodgerblue4", size = 1) +
    geom_tiplab(size = 3, offset = 4, fontface = "italic") +
    theme_tree2()

p2 <- revts(p1) +
    scale_x_continuous(
        labels = abs, breaks = seq(-350, 0, by = 50)
    ) +
    coord_cartesian(xlim = c(-350, 100), ylim = c(1, 12)) +
    geom_rect(data = periods, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  group = id, fill = factor(id)),
              alpha = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = period_col) +
    geom_vline(xintercept = -66, color = "white", alpha = 0.5) +
    geom_vline(xintercept = -145, color = "white", alpha = 0.5) +
    geom_vline(xintercept = -201, color = "white", alpha = 0.5) +
    geom_vline(xintercept = -252, color = "white", alpha = 0.5) +
    geom_vline(xintercept = -299, color = "white", alpha = 0.5)

# Add rectangles representing WGD events
peaks_and_nodes <- data.frame(
    node = c(
        "Glycine_max", 21, "Populus_trichocarpa",
        19, 19, 16, 14, "Zea_mays", 22, 22, "Solanum_lycopersicum"
    ),
    min = c(
        11.87, 56.04, 32.60, 49.27, 54.58,
        125, 297, 19.71, 63.08, 110, 57.47
    ),
    max = c(
        13.99, 67, 36.34, 50.99, 69.38,
        135, 319, 20.99, 69.89, 117, 64.84
    )
)

wgds <- get_wgd_coords(peaks_and_nodes, p2)

final_tree <- p2 +
    geom_rect(
        data = wgds, inherit.aes = FALSE,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "tomato3", color = "black", size = 0.2
    ) +
    labs(
        title = "Species tree with WGD events highlighted"
    )

final_tree$data <- final_tree$data |>
    mutate(label = str_replace_all(label, "_", " "))

final_tree

### Create figure 4
fig4 <- wrap_plots(
    wrap_plots(
        final_tree, ggplotify::as.ggplot(p_heatmap_z_scores),
        widths = c(1, 2)
    ),
    p_mean_z_scores +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        scale_y_continuous(limits = c(0, 800)),
    nrow = 2,
    heights = c(1.5, 1)
)

fig4 <- ggarrange(
    ggarrange(
        final_tree, ggplotify::as.ggplot(p_heatmap_z_scores),
        widths = c(1, 2)
    ),
    p_mean_z_scores +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        scale_y_continuous(limits = c(0, 800)),
    nrow = 2,
    heights = c(1.2, 1)
)

fig4


### Export figure - manually modify later to include motif icons
ggsave(
    fig4,
    file = here("products", "figs", "fig4.svg"),
    width = 16, height = 10
)


## Figure 5 ====
### Load data
load(here("products", "plots", "p_heatmap_functional_enrichment.rda"))


fig5 <- ComplexHeatmap::draw(
    p_heatmap_functional_enrichment, merge_legend = TRUE
)


png(
    filename = here("products", "figs", "fig5.png"),
    width = 10, height = 12, units = "in", res = 400
)
fig5
dev.off()



# Supplementary Figures ====
## Sup. Fig. S1 ====

# This figure was created manually using Inkscape



## Sup. Fig. S2 ====
### Load data
load(here("products", "plots", "ks_plots_with_boundaries.rda"))

### Create figure
sup_fig2 <- wrap_plots(
    # Row 1: G. max and A. thaliana
    wrap_plots(
        ks_plots_with_boundaries[[1]][[1]],
        ks_plots_with_boundaries[[1]][[2]] + ylab("")
    ),
    # Row 2: O. sativa and Z. mays
    wrap_plots(
        ks_plots_with_boundaries[[1]][[3]],
        ks_plots_with_boundaries[[1]][[4]] + ylab("")
    ),
    # Row 3: P. trichocarpa, S. lycopersicum, V. vinifera
    wrap_plots(
        ks_plots_with_boundaries[[2]][[1]],
        ks_plots_with_boundaries[[2]][[2]] + ylab(""),
        ks_plots_with_boundaries[[2]][[3]] + ylab(""),
        nrow = 1
    ),
    nrow = 3
)

sup_fig2

### Save figure
ggsave(
    sup_fig2,
    file = here("products", "figs", "sup_figure_s2.png"),
    width = 10, height = 10, dpi = 300
)


## Sup. Fig. S3 ====
### Load data
load(here::here("products", "plots", "p_degree_distros_ppi_network.rda"))

sup_fig3 <- p_degree_distros_ppi_network

### Export figure
ggsave(
    sup_fig3,
    file = here("products", "figs", "sup_figure_s3.pdf"),
    width = 10, height = 7
)

ggsave(
    sup_fig3,
    file = here("products", "figs", "sup_figure_s3.png"),
    width = 10, height = 8, dpi = 300
)

## Sup. Fig. S4 ====
### Load data
load(here("products", "plots", "p_degree_distros_grn.rda"))

sup_fig4 <- p_degree_distros_grn

### Export figure
ggsave(
    sup_fig4,
    file = here("products", "figs", "sup_figure_s4.png"),
    width = 11, height = 8, dpi = 300
)


## Sup. Fig. S5 ====
### Load data
load(here("products", "plots", "p_sorensendice_grn_targets.rda"))

### Create figure
sup_fig5 <- p_sorensendice_grn_targets

### Save figure
ggsave(
    sup_fig5,
    file = here("products", "figs", "sup_figure_s5.png"),
    width = 11, height = 10, dpi = 300
)


## Sup. Fig. S6 ====
### Load data
load(here("products", "plots", "p_sorensendice_grn_tfs.rda"))

### Create figure
sup_fig6 <- p_sorensendice_grn_tfs

### Save figure
ggsave(
    sup_fig6,
    file = here("products", "figs", "sup_figure_s6.png"),
    width = 10, height = 10, dpi = 300
)



