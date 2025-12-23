library(tidyverse)
set.seed(1001)
# --- your fixed inputs (unchanged) ---
sc1 <- c(0.35, 0.30, 0.20, 0.10, 0.05)
sc2 <- c(0.10, 0.15, 0.20, 0.25, 0.30)
sc3 <- c(0.10, 0.20, 0.40, 0.20, 0.10)
sc4 <- c(0.20, 0.20, 0.20, 0.20, 0.20)

b_shape <- as.numeric(collect_results[4, 3])
b_rate  <- as.numeric(collect_results[4, 4])
g_shape <- as.numeric(collect_results[4, 5])
g_rate  <- as.numeric(collect_results[4, 6])

no.d <- 5

# --- compute proposal (tidy style) ---
clin_param <- rgamma(10000, b_shape, b_rate)
pat_param  <- rgamma(10000, g_shape, g_rate)

clin_rec <- map_int(clin_param, ~ which.min(abs(skeletonc^.x - 0.25)))
pat_rec  <- map_int(pat_param,  ~ which.min(abs(skeletonp^.x - 0.35)))

# pmin -> pick the minimum of the two recovery indices, then tabulate to proportions
min_rec <- pmin(clin_rec, pat_rec)
proposal_props <- tabulate(min_rec, nbins = no.d) / length(min_rec)

# Create tidy data frames for plotting
plot_df <- tibble(
  cat = factor(1:no.d, labels = paste0("d", 1:no.d)), 
  proposal = proposal_props
)

sc2_df <- tibble(
  cat = factor(1:no.d, labels = paste0("d", 1:no.d)),
  sc2 = sc2
)

# --- ggplot: white bars with black border, red semi-transparent overlay ---
ggplot(plot_df, aes(x = cat, y = proposal)) +
  geom_col(fill = "#8DA0CB", alpha = 0.45) +               # base bars: white fill, black border
  geom_col(data = sc2_df, aes(x = cat, y = sc2),                         # overlay bars
           color="black", fill="white", alpha = 0, position = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "",
       x = "Category", y = "Proportion") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())
