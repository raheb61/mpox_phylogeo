##################################################################
# Summaries of the typed tree of B.1.20 subsample from BDMM-Prime
#################################################################

library(treeio)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)



tree_file <- "MCC_trees/B.1.20_subsampled499_combined_downsampled_BDMM_Prime_single_rates.MCC.tree"   # typed tree from BDMM-Prime (NEXUS)
OUTDIR    <- "grafics"

country_cols <- c("Germany" = "#D55E00", "Portugal" = "#0072B2","USA" = "#009E73", "Colombia" = "#9810FA",
                 "Australia" = "#E69F00")


#read tree

td  <- treeio::read.beast(tree_file)
phy <- td@phylo
ann <- td@data %>% mutate(node = as.integer(node))

#lookup vectors 
node_type <- rep(NA_character_, max(phy$edge))
node_type[ann$node] <- ann[["type"]]

node_prob <- rep(NA_real_, max(phy$edge))
node_prob[ann$node] <- as.numeric(ann[["type.prob"]])

depth     <- ape::node.depth.edgelength(phy)
node_date <- as.Date("2025-12-30") - (max(depth) - depth) * 365.25   
ntip      <- 499


edges <- tibble(
  parent      = phy$edge[, 1],
  child       = phy$edge[, 2],
  parent_type = node_type[phy$edge[, 1]],
  child_type  = node_type[phy$edge[, 2]],
  parent_date = node_date[phy$edge[, 1]],
  child_date  = node_date[phy$edge[, 2]],
  parent_prob = node_prob[phy$edge[, 1]],
  child_prob  = node_prob[phy$edge[, 2]],
  child_is_tip = phy$edge[, 2] <= ntip
)


#count type changes

changes <- edges %>% filter(parent_type != child_type)


changes_table <- changes %>%
  group_by(from = parent_type, to = child_type) %>%
  summarise(n_changes = n(),
            median_parent_prob = median(parent_prob, na.rm = TRUE),
            median_child_prob  = median(child_prob,  na.rm = TRUE),
            min_endpoint_prob  = suppressWarnings(
              min(c(parent_prob, child_prob), na.rm = TRUE)),
            .groups = "drop") %>%
  arrange(desc(n_changes))

write.table(changes_table, file.path(OUTDIR, "type_changes_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)



changes_events <- changes %>%
  transmute(from = parent_type, to = child_type,
            date = child_date,
            interval_start = parent_date, interval_end = child_date,
            from_type_prob = parent_prob, to_type_prob = child_prob) %>%
  arrange(date)

write.table(changes_events, file.path(OUTDIR, "type_change_events.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)



germany <- changes %>%
  filter(parent_type == "Germany" | child_type == "Germany") %>%
  mutate(direction = if_else(child_type == "Germany",
                             paste("Into", "Germany"), paste("Out of", "Germany")),
         other = if_else(child_type == "Germany", parent_type, child_type))

germany_table <- germany %>%
  group_by(direction, other) %>%
  summarise(n_changes = n(),
            median_from_prob = median(parent_prob, na.rm = TRUE),
            median_to_prob   = median(child_prob,  na.rm = TRUE),
            .groups = "drop") %>%
  arrange(direction, desc(n_changes))

write.table(germany_table, file.path(OUTDIR, "type_changes_germany.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


#ltt by type plot

grid <- seq(min(node_date), max(node_date), by = "1 day")

ltt <- lapply(unique(edges$child_type), function(ty) {
  e <- edges[edges$child_type == ty, ]
  tibble(date = grid, type = ty,
         n = vapply(grid, function(t)
           sum(e$parent_date <= t & e$child_date > t), integer(1)))
}) %>% bind_rows()



p_ltt <- ggplot(ltt, aes(date, n, colour = type)) +
  geom_step(linewidth = 0.6) +
  scale_colour_manual(values = country_cols, name = "Deme") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Number of lineages") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
p_ltt
ggsave(file.path(OUTDIR, "ltt_by_type.pdf"), p_ltt,
       width = 17, height = 8, units = "cm")
write.table(ltt, file.path(OUTDIR, "ltt_by_type.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

####quantify german import clusters

children_of <- split(phy$edge[, 2], phy$edge[, 1])


focal_cluster_tips <- function(node) {
  out <- integer(0); stack <- node
  while (length(stack)) {
    x <- stack[length(stack)]; stack <- stack[-length(stack)]
    if (node_type[x] != "Germany") next         
    if (x <= ntip) { out <- c(out, x); next }
    kids <- children_of[[as.character(x)]]
    if (!is.null(kids)) stack <- c(stack, kids)
  }
  out
}

import_edges <- edges %>% filter(parent_type != "Germany", child_type == "Germany")

root_node <- 500
cluster_roots <- tibble(
  node        = import_edges$child,
  origin      = import_edges$parent_type,  
  import_date = import_edges$child_date,
  stem_start  = import_edges$parent_date,   
  root_type_prob   = import_edges$child_prob,
  origin_type_prob = import_edges$parent_prob
)

clusters <- lapply(seq_len(nrow(cluster_roots)), function(i) {
  r    <- cluster_roots[i, ]
  tips <- focal_cluster_tips(r$node)
  tibble(
    cluster_id   = i,
    origin       = r$origin,
    import_date  = r$import_date,
    stem_start   = r$stem_start,
    root_type_prob   = r$root_type_prob,
    origin_type_prob = r$origin_type_prob,
    n_tips       = length(tips),
    first_sample = if (length(tips)) min(node_date[tips]) else as.Date(NA),
    last_sample  = if (length(tips)) max(node_date[tips]) else as.Date(NA),
  )
}) %>% bind_rows() %>%
  mutate(span_days = as.numeric(last_sample - first_sample)) %>%
  arrange(desc(n_tips), import_date)

write.table(clusters, 
            file.path(OUTDIR, "german_import_clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


imports <- clusters 




#distribution import cluster sizes

size_dist <- imports %>%
  count(n_tips, name = "n_clusters") %>%
  arrange(n_tips)


lab <- sprintf(
  "%d import clusters | median size %s tips | %d singletons (%.0f%%) | largest %d",
  nrow(imports), format(median(imports$n_tips)),
  sum(imports$n_tips == 1),
  100 * mean(imports$n_tips == 1), max(imports$n_tips))


p_size <- ggplot(size_dist, aes(x = n_tips, y = n_clusters)) +
  geom_segment(aes(xend = n_tips, y = 0, yend = n_clusters),
               colour = "grey70", linewidth = 0.5) +
  geom_point(aes(fill = n_tips), shape = 21, size = 4.5,
             colour = "white", stroke = 0.6) +
  geom_text(aes(label = n_clusters), vjust = -1.4, size = 2.7,
            colour = "grey25") +
  scale_fill_viridis_c(option = "mako", begin = 0.15, end = 0.85,
                       direction = -1, guide = "none",
                       trans = "sqrt") +
  scale_x_continuous(breaks = size_dist$n_tips,
                     expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq(0, max(size_dist$n_clusters), by = 1),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(x = paste("Cluster size (number of German tips)"),
       y = "Number of clusters",
       subtitle = lab) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.25,
                                          colour = "grey90"),
        axis.ticks.x = element_line(linewidth = 0.25, colour = "grey60"),
        axis.ticks.length.x = unit(2, "pt"),
        plot.subtitle = element_text(size = 7, colour = "grey35",
                                     margin = margin(b = 6)))

p_size
ggsave(file.path(OUTDIR, "import_cluster_size_distribution_lollipop.pdf"), p_size,
       width = 19, height = 10, units = "cm")
#timeline of clusters

tip_dates <- lapply(seq_len(nrow(cluster_roots)), function(i) {
  tips <- focal_cluster_tips(cluster_roots$node[i])
  if (!length(tips)) return(NULL)
  tibble(cluster_id = i, date = node_date[tips])
}) %>% bind_rows()


ord <- clusters %>%
  arrange(import_date, desc(n_tips)) %>%
  mutate(lab = sprintf("C%d  (%d)", cluster_id, n_tips),
         lab = factor(lab, levels = lab))


segs <- bind_rows(
  ord %>% transmute(lab, cluster_id,
                    x = import_date, xend = first_sample,
                    phase = "Import to first sample"),
  ord %>% transmute(lab, cluster_id,
                    x = first_sample, xend = last_sample,
                    phase = "First to last sample")
) %>%
  filter(!is.na(x), !is.na(xend)) %>%
  mutate(phase = factor(phase, levels = c("Import to first sample",
                                          "First to last sample")))

tip_pts <- tip_dates %>% inner_join(ord %>% select(cluster_id, lab),
                                    by = "cluster_id")

phase_cols <- c("Import to first sample" = "grey70",
                "First to last sample"   = "#1b7837")


p_tl <- ggplot(segs, aes(y = lab)) +
  geom_segment(aes(x = x, xend = xend, yend = lab,
                   colour = phase, linewidth = phase)) +
  geom_point(data = tip_pts, aes(x = date), shape = 21, size = 1.6,
             fill = "white", colour = "grey25", stroke = 0.5) +
  geom_point(data = ord, aes(x = import_date), shape = 23, size = 2,
             fill = "white", colour = "grey30", stroke = 0.6) +
  scale_colour_manual(values = phase_cols, name = NULL) +
  scale_linewidth_manual(values = c("Import to first sample" = 0.35,
                                    "First to last sample"   = 1.4),
                         guide = "none") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y",
               expand = expansion(mult = 0.03)) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25,
                                          colour = "grey90"),
        axis.text.y = element_text(family = "mono", size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.key.width = unit(14, "pt"))
p_tl
ggsave(file.path(OUTDIR, "cluster_timeline.pdf"), p_tl,
       width = 17, height = 11, units = "cm")

write.table(ord %>% select(cluster_id, origin, import_date, first_sample,
                           last_sample, span_days, n_tips,
                           root_type_prob, origin_type_prob),
            file.path(OUTDIR, "cluster_timeline.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
