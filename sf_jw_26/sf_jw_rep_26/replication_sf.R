require(tidyverse)
require(ggplot2)
require(fixest)
require(scales)
require(patchwork)
require(usmap)
require(sf)
require(ggrepel)
require(furrr)

# Load Data
remove(list = ls())
timeseries <- read.csv("replication/timeseries.csv")
paneldata <- read.csv("replication/paneldata.csv")

######################################################## TIME SERIES ###########

m1 <- fenegbin(Total_Events ~ scale(lag_Total_Posts) + scale(lag_BLM) + Election + scale(lag_Unemp) | week_num + year_num , data = timeseries)
m2 <- fenegbin(Total_Posts ~ scale(lag_Total_Events) + scale(lag_BLM) + Election + scale(lag_Unemp) | week_num + year_num, data = timeseries)
etable(m1, m2, cluster =~ week_num + year_num)

######################################################## PANEL DATA ############

# Modeling offline to online
m1 <- fenegbin(Number_of_out_main ~ scale(lag_events_1) | state_name + week_num + year, data = paneldata)
m2 <- fenegbin(Number_of_out_state ~ scale(lag_events_1)  | state_name + week_num + year, data = paneldata)
m3 <- fenegbin(Number_of_out_main ~ scale(lag_events_1)  + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_blm_1) + scale(Number_of_original) + scale(lag_unemp_1) + election | state_name + week_num + year, data = paneldata)
m4 <- fenegbin(Number_of_out_state ~ scale(lag_events_1) + scale(lag_inmain_1) + scale(lag_instate_1)  + scale(lag_blm_1)+ scale(Number_of_original) + scale(lag_unemp_1) + election | state_name + week_num + year, data = paneldata)
etable(m1, m2, m3, m4, cluster =~ state_name + week_num)

# Modeling online to offline
m1 <- fenegbin(Number_of_events ~ scale(lag_outmain_1) | state_name + week_num + year, data = paneldata)
m2 <- fenegbin(Number_of_events ~ scale(lag_outstate_1)  | state_name + week_num + year, data = paneldata)
m3 <- fenegbin(Number_of_events ~ scale(lag_orig_1) + scale(lag_outmain_1) + scale(lag_outstate_1) | state_name + week_num + year, data = paneldata)
m4 <- fenegbin(Number_of_events ~ scale(lag_orig_1) + scale(lag_outmain_1) + scale(lag_outstate_1) + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_blm_1) + scale(lag_unemp_1) + election | state_name + week_num + year, data = paneldata)
etable(m1, m2, m3, m4, cluster =~ state_name + week_num) 

# Permutation checks 

permute_fixest <- function(
    data,
    formula,
    fe = NULL,            
    family = "negbin",    
    y_var,                
    param,                
    param_label = NULL,   
    n_perm = 10000,
    cluster = NULL,
    within = NULL,        
    parallel = TRUE,
    seed = 123
) {
  set.seed(seed)
  
  formula_fe <- as.formula(paste(formula, "|", paste(fe, collapse = "+")))
  
  if (family == "negbin") {
    mod_obs <- fenegbin(formula_fe, data = data)
  } else {
    mod_obs <- feglm(formula_fe, data = data, family = family)
  }
  beta_obs <- coef(mod_obs)[param]
  
  if (parallel) plan(multisession, workers = parallel::detectCores() - 1)
  
  message("Running ", n_perm, " permutations...")
  
  perm_coefs <- future_map_dbl(
    1:n_perm,
    function(i) {
      perm_data <- data
      if (!is.null(within)) {
        # Permute outcome within clusters (e.g., states)
        perm_data <- perm_data %>%
          group_by(across(all_of(within))) %>%
          mutate("{y_var}" := sample(.data[[y_var]])) %>%
          ungroup()
      } else {
        perm_data[[y_var]] <- sample(perm_data[[y_var]])
      }
      
      mod_perm <- tryCatch({
        if (family == "negbin") {
          fenegbin(formula_fe, data = perm_data)
        } else {
          feglm(formula_fe, data = perm_data, family = family)
        }
      }, error = function(e) NULL)
      
      if (is.null(mod_perm)) return(NA_real_)
      coef(mod_perm)[param]
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)  
  )
  
  p_val <- mean(abs(perm_coefs) >= abs(beta_obs), na.rm = TRUE)
  
  df_plot <- tibble(coef_perm = perm_coefs) %>% filter(!is.na(coef_perm))
  var_label <- ifelse(is.null(param_label), param, param_label)
  
  p <- ggplot(df_plot, aes(x = coef_perm)) +
    geom_histogram(bins = 40, fill = "grey", color = "white", alpha = 0.8) +
    geom_vline(aes(xintercept = beta_obs), color = "black", linewidth = 1.3) +
    theme_bw(base_size = 14, base_family = "Helvetica") +
    labs(
      title = "",
      subtitle = paste0(var_label),
      x = paste(""),
      y = "Frequency"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0, color = "#333333"),
      plot.subtitle = element_text(size = 12, hjust = 0, color = "#555555"),
      axis.text = element_text(size = 10, color = "#444444"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.major.y = element_line(color = "gray90"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  list(
    observed = beta_obs,
    p_val = p_val,
    perm_coefs = perm_coefs,
    plot = p
  )
}

# Offline to Online: Horizontal
res_out_horiztonal <- permute_fixest(
  formula = "Number_of_out_state ~ scale(lag_events_1) + scale(lag_blm_1) + scale(Number_of_original) + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_unemp_1) + election", fe = c("state_name", "week_num", "year"),
  data = paneldata, family = "negbin", y_var = "Number_of_out_state", param = "scale(lag_events_1)", n_perm = 10000,   within = "state_name", param_label = "Events to Horizontal Exposure", parallel = TRUE)
res_out_horiztonal$p_val
res_out_horiztonal$plot

# Offline to Online: Vertical
res_out_vertical <- permute_fixest(
  formula = "Number_of_out_main~ scale(lag_events_1) + scale(lag_blm_1) + scale(Number_of_original) + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_unemp_1) + election", fe = c("state_name", "week_num", "year"),
  data = paneldata, family = "negbin", y_var = "Number_of_out_main", param = "scale(lag_events_1)", n_perm = 10000, within = "state_name", param_label = "Events to Vertical Exposure", parallel = TRUE)
res_out_vertical$p_val
res_out_vertical$plot

# Online to Offline: Horizontal
res_event_horizontal <- permute_fixest(
  formula = "Number_of_events ~ scale(lag_orig_1) + scale(lag_outmain_1) + scale(lag_outstate_1) + scale(lag_blm_1) + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_unemp_1) + election", fe = c("state_name", "week_num", "year"),
  data = paneldata, family = "negbin", y_var = "Number_of_events", param = "scale(lag_outstate_1)", n_perm = 10000, within = "state_name", param_label = "Horizontal Exposure to Events", parallel = TRUE)
res_event_horizontal$p_val
res_event_horizontal$plot

# Online to Offline: Vertical
res_event_vertical <- permute_fixest(
  formula = "Number_of_events ~ scale(lag_orig_1) + scale(lag_outmain_1) + scale(lag_outstate_1) + scale(lag_blm_1) + scale(lag_inmain_1) + scale(lag_instate_1) + scale(lag_unemp_1) + election", fe = c("state_name", "week_num", "year"),
  data = paneldata, family = "negbin", y_var = "Number_of_events", param = "scale(lag_outmain_1)", n_perm = 10000, within = "state_name", param_label = "Vertical Exposure to Events", parallel = TRUE)
res_event_vertical$p_val
res_event_vertical$plot

perm_grid <- (
  (res_out_horiztonal$plot + res_out_vertical$plot) /
    (res_event_horizontal$plot + res_event_vertical$plot)
) +
  plot_annotation(
    title = "Permutation-Based Hypothesis Tests for Coefficient Significance (N = 10,000)",
    subtitle = "Empirical Distributions of Permuted Estimates vs. Observed β Values",
    theme = theme_bw(base_size = 14, base_family = "Helvetica") +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "#333333"),
        plot.subtitle = element_text(size = 14, hjust = 0, color = "#555555"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "#444444"),
        axis.text.y = element_text(size = 10, color = "#444444"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "gray90"),
        panel.grid.major.y = element_line(color = "gray90"),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 11),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  )

perm_grid

ggsave(
  filename = "replication/permutation_grid.png",  # output file name
  plot = perm_grid,                   
  width = 10,                         
  height = 8,                        
  dpi = 600                           
)

######################################################## PLOTTING #############

#-------------------------------------------- US Map Plot ---------------------

paneldata <- paneldata %>% mutate(state_name = str_replace_all(tolower(state_name), "_", " "))

# Summary state
state_summary <- paneldata %>%
  group_by(state_name) %>%
  summarise(
    total_posts  = sum(Number_of_posts, na.rm = TRUE),
    total_events = sum(Number_of_events, na.rm = TRUE),
    
    first_week_posts = min(Week[Number_of_posts > 0], na.rm = TRUE),
    last_week_posts  = max(Week[Number_of_posts > 0], na.rm = TRUE),
    
    first_week_events = min(Week[Number_of_events > 0], na.rm = TRUE),
    last_week_events  = max(Week[Number_of_events > 0], na.rm = TRUE)
  )
cor.test(state_summary$total_events, state_summary$total_posts)

# Merge geom
map_data <- usmap::us_map(regions = "states") %>%
  mutate(full = str_to_lower(full)) %>%
  left_join(state_summary, by = c("full" = "state_name"))

map_data_tosave <- as.data.frame(map_data) %>% dplyr::select(-geom)

# Compute centroid# Compute centroidgeom
centroids <- map_data %>%
  st_centroid() %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(full, total_posts, total_events, x, y)
head(map_data)

# Plot
mapplot <- ggplot() +
  # Base fill: total posts
  geom_sf(data = map_data, aes(fill = total_posts), color = NA) +
  
  # Add state boundaries
  geom_sf(data = map_data, fill = NA, color = "gray20", size = 0.3) +
  
  # Color and size scales
  scale_fill_gradient(
    name = "Total Posts",
    low = "white", high = "grey",
    labels = scales::comma,
    limits = c(0, 2900)
  ) +
  scale_size_continuous(
    name = "Total Events",
    range = c(3, 16),
    labels = scales::comma
  ) +
  
  # Circles for all states (including zeros)
  geom_point(
    data = centroids,
    aes(x = x, y = y, size = total_events),
    color = "black", alpha = 0.8
  ) +
  
  # Labels only for states with > 0 events
  geom_text(
    data = centroids %>% filter(total_events > 0),
    aes(x = x, y = y, label = total_events),
    color = "white",
    fontface = "bold",
    size = 3.5,
    show.legend = FALSE
  ) +
  
  # Titles and theme
  labs(
    title = "White Lives Matter: Online Posts and Offline Events by State",
    subtitle = "Fill = Total Posts | Circle = Total Events",
  ) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "#333333"),
    plot.subtitle = element_text(size = 14, hjust = 0, color = "#555555"),
    plot.caption = element_text(size = 10, color = "gray40", hjust = 1),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_line(color = "gray90"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 12)
  )
mapplot

ggsave(
  filename = "replication/mapplot_sum.png",  # output file name
  plot = mapplot,                   
  width = 10,                         
  height = 8,                        
  dpi = 600                          
)

#-------------------------------------------- Bubble Plot ------------------------------------------

# Summarize weekly data 
together_summary <- paneldata %>%
  mutate(Week = as.Date(Week)) %>%
  group_by(Week) %>%
  summarise(
    Vertical_Exposure   = mean(Number_of_out_main, na.rm = TRUE),
    Horizontal_Exposure = mean(Number_of_out_state, na.rm = TRUE)
  )

# Compute scaling constants 
max_vert <- max(together_summary$Vertical_Exposure, na.rm = TRUE)
max_horz <- max(together_summary$Horizontal_Exposure, na.rm = TRUE)

state_summary <- paneldata %>%
  mutate(state_name = str_to_title(str_replace_all(state_name, "_", " "))) %>%
  group_by(state_name) %>%
  summarise(
    total_orig      = sum(Number_of_original, na.rm = TRUE),
    total_events      = sum(Number_of_events, na.rm = TRUE),
    vertical_exposure = sum(Number_of_out_main, na.rm = TRUE),
    horizontal_exposure = sum(Number_of_out_state, na.rm = TRUE)
  )

cor.test(state_summary$vertical_exposure, state_summary$horizontal_exposure)
cor.test(state_summary$vertical_exposure, state_summary$total_events)
cor.test(state_summary$total_events, state_summary$horizontal_exposure)

bubble <- ggplot(state_summary, aes(
  x = horizontal_exposure,
  y = vertical_exposure,
  size = total_events,
  fill = total_events    # fill controls shading
)) +
  geom_point(shape = 21, color = "white", alpha = 0.9) +
  scale_fill_gradient(
    name = "Number of Events",
    low = "white",
    high = "grey50",
    guide = guide_colorbar(barwidth = 1, barheight = 8)
  ) +
  scale_size_continuous(
    range = c(3, 25),
    labels = comma
  ) +
  scale_x_continuous(
    labels = comma,
    name = "Horizontal Exposure",
    limits = c(-50, 600)
  ) +
  scale_y_continuous(
    labels = comma,
    name = "Vertical Exposure",
    limits = c(-50, 400)
  ) +
  geom_text_repel(
    data = state_summary %>% filter(total_events > 1),
    aes(label = state_name),
    size = 5,
    color = "black",
    fontface = "bold",
    max.overlaps = 6
  ) +
  labs(
    title = "White Lives Matter: Exposure Metrics and Offline Mobilization",
    subtitle = "Bubble Size & Shading = Number of Events | Position = Exposure Profile"
  ) +
  guides(size = "none") +   # 🔹 remove size legend
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "#333333"),
    plot.subtitle = element_text(size = 14, hjust = 0, color = "#555555"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "#444444"),
    axis.text.y = element_text(size = 10, color = "#444444"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.justification = "right",
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
bubble

ggsave(
  filename = "replication/bubble.png",  # output file name
  plot = bubble,                   
  width = 10,                         
  height = 8,                       
  dpi = 600                          
)

#-------------------------------------------- Post and Protest Plot ------------------------------------------

timeseries <- timeseries %>% mutate(Week = as.Date(Week))  
timeplot <- ggplot(timeseries, aes(x = as.Date(Week))) +
  
  # --- Events as histogram / density background ---
  geom_area(
    aes(y = Total_Events, fill = "Events"),
    alpha = 0.4,
    color = NA
  ) +
  
  # --- Posts as line ---
  geom_line(
    aes(y = Total_Posts / max(Total_Posts) * max(Total_Events),
        color = "Total Posts"),
    size = 1.2
  ) +
  
  # --- Scales and secondary axis ---
  scale_y_continuous(
    name = "Total Events",
    sec.axis = sec_axis(
      ~ . * max(timeseries$Total_Posts) / max(timeseries$Total_Events),
      name = "Total Posts"
    ),
    labels = comma
  ) +
  
  scale_x_date(
    name = "Week",
    date_labels = "%b %Y",
    date_breaks = "4 months"
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c("Total Posts" = "black")
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c("Events" = "grey40")
  ) +
  
  labs(
    title = "White Lives Matter: Weekly Telegram Posts and Offline Events",
    subtitle = "Offline Events = Shaded Density | Online Posts = Line",
    y = "Total Events"
  ) +
  
  theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "#333333"),
    plot.subtitle = element_text(size = 14, hjust = 0, color = "#555555"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "#444444"),
    axis.text.y = element_text(size = 10, color = "#444444"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
timeplot

ggsave(
  filename = "replication/timeseries.png",  # output file name
  plot = timeplot,                  
  width = 10,                        
  height = 8,                       
  dpi = 600                           
)

