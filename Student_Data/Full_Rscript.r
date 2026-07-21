library(damr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggetho)
library(ggplot2)
library(gridExtra)
library(sleepr)
library(zeitgebr)
library(stringr)
library(lme4)
library(survival)
library(survminer)
library(easystats)
library(ggfortify)
library(data.table)
library(ggforce)

# =============================================================================
# To make saved figures look pretty
# =============================================================================

custom_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size * 0.85),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.50),
      legend.position = c(0.9, 0.9),  # Legend now in the top-right
      panel.grid.major = element_line(color = "grey90", size = 0.3)  # Subtle gridlines
    )
}

# =============================================================================
# Set Data Directory and Load Metadata
# =============================================================================

data_dir <- "/Users/clt54/Documents/Monitor_Date" # change to your directory

list.files(data_dir)
setwd(data_dir)

# Load and link metadata to DAM results
metadata <- fread("metadata.csv") 

# Add 'batch' key based on start_datetime (important for activity if you did it in batches)
metadata[, batch := as.character(start_datetime)]

metadata <- link_dam_metadata(metadata, result_dir = data_dir)



# =============================================================================
# Load and Annotate Monitor Data
# =============================================================================

dt <- load_dam(metadata)

# =============================================================================
# Remove Dead Animals
# =============================================================================

dt_curated <- curate_dead_animals(dt, moving_var = activity)
summary(dt_curated)

# Remove dead samples
setdiff(dt[, id, meta = TRUE], dt_curated[, id, meta = TRUE])


# =============================================================================
# Lifespan Analysis
# =============================================================================

dt_lifespan <- dt_curated
lifespan_dt <- dt_lifespan[, .(lifespan = max(t)), by = id]
lifespan_dt[, treatment := xmv(treatment)]
lifespan_dt[, deathdays := lifespan / (24 * 60 * 60)]
lifespan_dt[, treatment := relevel(as.factor(treatment), ref = "Control")]

# Kaplan-Meier survival analysis
km_fit <- survfit(Surv(deathdays) ~ 1, data = lifespan_dt)
autoplot(km_fit)

km_trt_fit <- survfit(Surv(deathdays) ~ treatment, data = lifespan_dt)
autoplot(km_trt_fit)



surv_plot <- ggsurvplot(
  km_trt_fit,
  data = lifespan_dt,
  palette = c("#1f77b4", "#ff7f0e", "#2ca02c"),  # Added a third color for the new treatment
  linetype = "solid",
  size = 1.2,
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  legend.title = "",
  legend.labs = c("Control","Starvation (15%)", "Calorie restriction (20%)"), # Updated legend labels
  ggtheme = custom_theme(base_size = 12),
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  surv.median.line = "hv"
)

# Adjust the legend position if needed
surv_plot$plot <- surv_plot$plot + theme(legend.position = c(0.15, 0.25))
print(surv_plot)

ggsave("surviourship.png")

# Cox proportional hazards model
cox <- coxph(Surv(deathdays) ~ treatment, data = lifespan_dt)
summary(cox)

# =============================================================================
# Activity Analysis
# =============================================================================
dt[, batch := xmv(batch)]

# Get the earliest t per batch
t_starts <- dt[, .(t_start = min(t)), by = batch]
dt <- merge(dt, t_starts, by = "batch", all.x = TRUE)


# Compute POSIXct timestamps per batch
dt[, t_clock := as.POSIXct(batch, format = "%Y-%m-%d %H:%M:%S", tz = "UTC") + seconds(t - t_start)]




# Find first 8:00 AM on or after start time per batch
first_8am_dt <- dt[
  t_clock >= as.POSIXct(batch, format = "%Y-%m-%d %H:%M:%S", tz = "UTC") & hour(t_clock) == 8 & minute(t_clock) == 0,
  .(first_8am_t = min(t)),
  by = batch
]

# Get new zt0 per batch
new_zt0_dt <- dt[first_8am_dt, on = .(batch, t = first_8am_t), .(batch, new_zt0 = t_clock)]




# Update metadata
metadata <- fread("metadata.csv") 

# Ensure both are the same type for joining
metadata[, start_datetime := as.POSIXct(start_datetime)]
new_zt0_dt[, batch := as.POSIXct(batch)]

# Create a named vector: names = old start times, values = new zt0
lookup <- setNames(new_zt0_dt$new_zt0, new_zt0_dt$batch)

# Replace each start_datetime in metadata with corresponding new_zt0
metadata[, start_datetime := lookup[as.character(start_datetime)]]



# Link metadata to DAM data using corrected start_datetime (no file written!)
metadata <- link_dam_metadata(metadata, result_dir = data_dir)


dt <- load_dam(metadata)

dt_curated <- curate_dead_animals(dt, moving_var = activity)

# Limit data to 30 days from that point
dt_curated <- dt_curated[
  t %between% c(0, 30 * 24 * 60 * 60)
]

# Add binary movement indicator
dt_curated[, moving := activity > 0]

# Actogram with light/dark annotation
ggetho(dt_curated, aes(z = activity)) +
  stat_tile_etho() +
  stat_ld_annotations(ld_colours = c("grey", "black"))


# Sample 10% of unique ids (just to help you see the general patterns)
set.seed(123)  # for reproducibility
sample_ids <- sample(unique(dt_curated$id), size = 0.1 * length(unique(dt_curated$id)))

# Subset data
dt_subsampled <- dt_curated[id %in% sample_ids]

# Plot actogram with light/dark annotations
ggetho(dt_subsampled, aes(z = activity)) +
  stat_tile_etho() +
  stat_ld_annotations(ld_colours = c("grey", "black"))

# Double-plotted actograms
dt_curated[, uid := 1:.N, meta = TRUE]

# The below is probably unreadable (go to paginated plots)
ggetho(dt_curated, aes(z = moving), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ treatment + uid, ncol = 8, labeller = label_wrap_gen(multi_line = FALSE))

# =============================================================================
# Paginated Plots (uncomment if required. It usually isn't unless you want all the actograms)
# =============================================================================

  # p <- ggetho(dt_curated, aes(z = moving), multiplot = 2) +
  #   stat_bar_tile_etho() +
  #   facet_wrap_paginate(
  #     ~ treatment + uid,
  #     ncol = 2,
  #     nrow = 4,
  #     labeller = label_wrap_gen(multi_line = FALSE)
  #   )
  # 
  # for (i in seq_len(n_pages(p))) {
  #   p_save <- p +
  #     facet_wrap_paginate(
  #       ~ treatment + uid,
  #       ncol = 2,
  #       nrow = 4,
  #       page = i
  #     )
  #   
  #   ggsave(
  #     plot = p_save,
  #     filename = paste0("actogram_", i, ".pdf"),
  #     width = 210,
  #     height = 297,
  #     units = "mm"
  #   )
  # }

# =============================================================================
# Average activity profile
# =============================================================================

# Average activity profile
ggetho(dt_curated, aes(y = activity, colour = treatment)) +
  stat_pop_etho() +
  stat_ld_annotations() +
  theme_classic()

#Plot each treatment separately 
ggetho(dt_curated, aes(y=activity, colour=treatment)) +
  stat_pop_etho() +
  stat_ld_annotations() + facet_grid(treatment ~ .) +
  coord_cartesian(ylim = c(0, 1))  + theme_pubr()

ggsave("actvity_separate.pdf")




#### averages for weeks
#######################################


# Create a new variable 'week' based on time in seconds, ensuring t = 0 is week = 1
weekdt <- copy(dt_curated)
weekdt[, week := ceiling((t + 1) / (60 * 60 * 24 * 7))] # something weird here

# Plot using ggetho
# Get unique weeks
unique_weeks <- unique(weekdt$week)

# Function to generate plot for a specific week
plot_for_week <- function(week_number) {
  # Filter data for the specific week
  week_data <- weekdt[week == week_number]
  
  # Create the plot
  p <- ggetho(week_data, aes(y = activity, colour = treatment), time_wrap = mins(1440)) +
    stat_pop_etho() +
    stat_ld_annotations() +
    scale_y_continuous(name= "Fraction of time active", limits = c(0,1), labels = scales::percent) +
    theme_pubr() +
    ggtitle(paste("Week", week_number))
  
  # Print the plot
  print(p)
}


# Generate plots for each week using lapply
plots <- lapply(unique_weeks, plot_for_week) 

# Arrange plots into a grid
#grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol to change the number of columns

# Arrange plots into a grid and save as a grob
g <- grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol to change the number of columns

# Save the combined plot as a PNG
ggsave("activity_week.pdf", plot = g, width = 12, height = 8, units = "in", dpi = 300)


# =============================================================================
# Activity Statistics
# =============================================================================

dt_binned <- bin_apply_all(dt_curated, activity)
dt_binned[, treatment := xmv(treatment)]
dt_binned[, phase := ifelse(hour(t) < 12, "L", "D")]
dt_binned[, ttwo := scale(t)]

light <- subset(dt_binned, phase == "L")
dark <- subset(dt_binned, phase == "D")


### Repeated measures with mixed effect models
#light
lm1 <- lmer(activity ~ as.factor(treatment) * ttwo  + (ttwo | id), light)
lm2 <- lmer(activity ~ as.factor(treatment) + ttwo  + (ttwo | id), light)
lm3 <- lmer(activity ~ ttwo  + (ttwo | id), light)
lm4 <- lmer(activity ~ as.factor(treatment)  + (ttwo | id), light)
anova(lm1,lm2) #testing interaction
anova(lm2, lm3) # testing main effect treatment
anova(lm2, lm4) # testing main effect ttwo (time)


library(emmeans)
emmip(lm1, treatment ~ ttwo, cov.reduce = range) +theme_pubr() + xlab("Standarised time") + ylab ("Day time activity (Linear prediction)") # might be too complicated to use in thesis

ggsave("activity_emmip_day.pdf") # not useful for thesis or publication, just to help you understand the stats


#dark
lm1 <- lmer(activity ~ as.factor(treatment) * ttwo  + (ttwo | id), dark)
lm2 <- lmer(activity ~ as.factor(treatment) + ttwo  + (ttwo | id), dark)
lm3 <- lmer(activity ~ ttwo  + (ttwo | id), dark)
lm4 <- lmer(activity ~ as.factor(treatment)  + (ttwo | id), dark)
anova(lm1,lm2) #testing interaction
anova(lm2, lm3) # testing main effect treatment
anova(lm2, lm4) # testing main effect ttwo (time)

emmip(lm1, treatment ~ ttwo, cov.reduce = range) +theme_pubr() + xlab("Standarised time") + ylab ("Night time activity (Linear prediction)") # might be too complicated to use in thesis
ggsave("activity_emmip_night.pdf")# not useful for thesis or publication, just to help you understand the stats

# =============================================================================
# Rest bouts analysis
# =============================================================================

#Upload monitor data with metadata
dt <- load_dam(metadata, FUN = sleepr::sleep_dam_annotation)

# Preliminary Sleep graph - all animals, overall 
ggetho(dt, aes(z=asleep)) +
  stat_tile_etho() 

#Curate to remove dead
dt_curated <- curate_dead_animals(dt)
summary(dt_curated, meta=T)

# Limit data to 30 days from that point
#dt_curated <- dt_curated[
#  t %between% c(0, 30 * 24 * 60 * 60)
#]

bout_dt <- bout_analysis(asleep, dt_curated)
bout_dt <- bout_dt[asleep == TRUE, -"asleep"]


# =============================================================================
# Average rest profile
# =============================================================================

# Average activity profile
ggetho(bout_dt, aes(y=duration / 60, colour = treatment)) +
  stat_pop_etho() +
  stat_ld_annotations() +
  theme_classic()

#Plot each treatment separately 
ggetho(bout_dt, aes(y=duration / 60, colour=treatment)) +
  stat_pop_etho() +
  stat_ld_annotations() + facet_grid(treatment ~ .) +
  #coord_cartesian(ylim = c(0, 1))  + 
  theme_pubr()

ggsave("rest_separate.pdf")

#### Rest averages for weeks
#######################################

weekdt <- copy(bout_dt)

seconds_per_week <- 7 * 24 * 60 * 60

# t already begins at zero after resetting the experiment to the first 08:00
weekdt[, week := floor(t / seconds_per_week) + 1L]

# Remove any accidental week beyond the intended 30-day window
weekdt <- weekdt[week <= 5]

plot_for_week <- function(week_number) {
  
  week_data <- copy(weekdt[week == week_number])
  
  # Reset time so that every weekly plot begins at t = 0
  week_data[, t := t - (week_number - 1) * seconds_per_week]
  
  ggetho(
    week_data,
    aes(y = duration / 60, colour = treatment),
    time_wrap = hours(24)
  ) +
    stat_pop_etho() +
    stat_ld_annotations() +
    labs(
      title = paste("Week", week_number),
      x = "Time of day",
      y = "Rest-bout duration (minutes)",
      colour = "Treatment"
    ) +
    theme_classic()
}

unique_weeks <- sort(unique(weekdt$week))
plots <- lapply(unique_weeks, plot_for_week)

g <- gridExtra::arrangeGrob(
  grobs = plots,
  ncol = 2
)

grid::grid.newpage()
grid::grid.draw(g)

ggsave(
  filename = "rest_week.pdf",
  plot = g,
  width = 12,
  height = 8,
  units = "in"
)
