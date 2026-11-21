# Assignment 4 - LMs/GLMs on focal dataset
# Mac Grosscup
#
# Question: how do plant traits change across a salinity gradient?
#
# Data: my field data (Masters_Data_Modified.xlsx, in the same GitHub repo).
# 5 coastal sites, 31 plots, soil salinity 0.03-4 mS, 38 species.
#
# Traits:
#   leaf thickness - expect thicker/more succulent leaves in salty soil
#   plant height - expect shorter plants from salt stress
#
# Using a Gamma GLM w/ log link for both - traits are positive, right
# skewed, and spread goes up with the mean. exp(coef) = multiplier per 1 mS.
#
# Models (same for both traits):
#   M1  ~ salinity                 one trend for everything
#   M2  ~ species                  species differ, no salinity effect
#   M3  ~ species + salinity       same % change per mS for all species (best from Assignment 3)
#   M4  ~ species + log(salinity)  effect levels off at high salinity (new)
#   M5  ~ species * salinity       each species gets its own slope (new)

install.packages(c("readxl", "ggplot2", "bbmle"))

library(readxl)
library(ggplot2)
library(bbmle)

plants <- read_excel("Masters_Data_Modified.xlsx")

# fix species name typos, claude helped me with this part
plants$species <- trimws(plants$species)
plants$species <- gsub("Sporobulus", "Sporobolus", plants$species)
plants$species <- gsub("Langularia", "Laguncularia", plants$species)

# same rows for all models - drop NAs and the 2 plants w/ height 0
# (Gamma can't handle 0s)
plants <- plants[!is.na(plants$leaf_thickness_mm) &
                   !is.na(plants$plant_height_cm) &
                   plants$plant_height_cm > 0, ]
nrow(plants)   # 1536 plants

# Sesuvium as baseline species (most common, found across whole gradient), Claude helped me with this
plants$species <- relevel(factor(plants$species), ref = "Sesuvium_portulacastrum")

ggplot(plants, aes(x = soil_salinity_ms, y = leaf_thickness_mm)) +
  geom_point(alpha = 0.2) +
  scale_y_log10() +
  labs(x = "soil salinity (mS)", y = "leaf thickness (mm)")

ggplot(plants, aes(x = soil_salinity_ms, y = plant_height_cm)) +
  geom_point(alpha = 0.2) +
  scale_y_log10() +
  labs(x = "soil salinity (mS)", y = "plant height (cm)")


# Trait 1: leaf thickness -------------------------------------------------------

thick.m1 <- glm(leaf_thickness_mm ~ soil_salinity_ms,
                family = Gamma(link = "log"), data = plants)
thick.m2 <- glm(leaf_thickness_mm ~ species,
                family = Gamma(link = "log"), data = plants)
thick.m3 <- glm(leaf_thickness_mm ~ species + soil_salinity_ms,
                family = Gamma(link = "log"), data = plants)
thick.m4 <- glm(leaf_thickness_mm ~ species + log(soil_salinity_ms),
                family = Gamma(link = "log"), data = plants)
thick.m5 <- glm(leaf_thickness_mm ~ species * soil_salinity_ms,
                family = Gamma(link = "log"), data = plants)

AICtab(thick.m1, thick.m2, thick.m3, thick.m4, thick.m5,
       base = TRUE, weights = TRUE)

# multiplier per 1 mS
exp(coef(thick.m1)["soil_salinity_ms"])   # ~1.20
exp(coef(thick.m3)["soil_salinity_ms"])   # ~1.06

# M5 slopes for a few species (Sesuvium = baseline, others add interaction)
b <- coef(thick.m5)
exp(b["soil_salinity_ms"])                                                    # Sesuvium ~1.15
exp(b["soil_salinity_ms"] + b["speciesLaguncularia_racemosa:soil_salinity_ms"])  # ~1.06
exp(b["soil_salinity_ms"] + b["speciesConocarpus_erectus:soil_salinity_ms"])     # ~1.05
exp(b["soil_salinity_ms"] + b["speciesAvicennia_germinans:soil_salinity_ms"])    # ~1.02
exp(b["soil_salinity_ms"] + b["speciesSporobolus_virginicus:soil_salinity_ms"])  # ~0.75

# M1 - mostly just which species grow where
#      (succulents in the salty plots)
# M2 - worse than anything with salinity in it
# M3 - all species ~6% thicker per mS
# M4 - about the same as M3, no real leveling off
# M5 - best by ~490 AIC
#
# Keeping M5. Most species get thicker leaves w/ salinity but at different
# rates. Sesuvium ~15%/mS, mangroves barely change (Avicennia ~2%), which
# makes sense since they exclude/excrete salt. S. virginicus actually gets
# thinner.
# 8 species only found at one salinity -> NA slopes. species only in the
# freshest plots (0.03-0.07 mS) get weird slopes, ignoring those.


# Trait 2: plant height ---------------------------------------------------------

height.m1 <- glm(plant_height_cm ~ soil_salinity_ms,
                 family = Gamma(link = "log"), data = plants)
height.m2 <- glm(plant_height_cm ~ species,
                 family = Gamma(link = "log"), data = plants)
height.m3 <- glm(plant_height_cm ~ species + soil_salinity_ms,
                 family = Gamma(link = "log"), data = plants)
height.m4 <- glm(plant_height_cm ~ species + log(soil_salinity_ms),
                 family = Gamma(link = "log"), data = plants)
height.m5 <- glm(plant_height_cm ~ species * soil_salinity_ms,
                 family = Gamma(link = "log"), data = plants)

AICtab(height.m1, height.m2, height.m3, height.m4, height.m5,
       base = TRUE, weights = TRUE)

exp(coef(height.m1)["soil_salinity_ms"])  # ~0.88
exp(coef(height.m3)["soil_salinity_ms"])  # ~1.01

b <- coef(height.m5)
exp(b["soil_salinity_ms"])                                                    # Sesuvium ~0.99
exp(b["soil_salinity_ms"] + b["speciesLaguncularia_racemosa:soil_salinity_ms"])  # ~1.02
exp(b["soil_salinity_ms"] + b["speciesConocarpus_erectus:soil_salinity_ms"])     # ~0.94
exp(b["soil_salinity_ms"] + b["speciesCoccoloba_uvifera:soil_salinity_ms"])      # ~0.90
exp(b["soil_salinity_ms"] + b["speciesUniola_paniculata:soil_salinity_ms"])      # ~1.40

# M1 - ~12% shorter per mS overall but ignores species, worst model
# M2/M3 - basically tied, shared slope shows almost no effect (~1%/mS)
# M4 - a little better than M3 but still small
# M5 - best by ~300 AIC
#
# Keeping M5. One shared slope hides species going opposite directions.
# Woody species mostly get shorter (Coccoloba ~10%/mS, Conocarpus ~6%),
# some grasses like Uniola get taller - salt tolerant dune grass.


# Model checks for the kept models (M5) -----------------------------------------

# predictions + 95% CI for 3 species found across most of the gradient
thick.pred <- predict(thick.m5, se.fit = TRUE)
plants$thick_fit <- exp(thick.pred$fit)
plants$thick_low <- exp(thick.pred$fit - 1.96 * thick.pred$se.fit)
plants$thick_high <- exp(thick.pred$fit + 1.96 * thick.pred$se.fit)

height.pred <- predict(height.m5, se.fit = TRUE)
plants$height_fit <- exp(height.pred$fit)
plants$height_low <- exp(height.pred$fit - 1.96 * height.pred$se.fit)
plants$height_high <- exp(height.pred$fit + 1.96 * height.pred$se.fit)

focal <- plants[plants$species %in% c("Sesuvium_portulacastrum",
                                      "Laguncularia_racemosa",
                                      "Conocarpus_erectus"), ]

ggplot(focal, aes(x = soil_salinity_ms)) +
  geom_ribbon(aes(ymin = thick_low, ymax = thick_high), fill = "red", alpha = 0.2) +
  geom_point(aes(y = leaf_thickness_mm), alpha = 0.4) +
  geom_line(aes(y = thick_fit), color = "red") +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "soil salinity (mS)", y = "leaf thickness (mm)")

ggplot(focal, aes(x = soil_salinity_ms)) +
  geom_ribbon(aes(ymin = height_low, ymax = height_high), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = plant_height_cm), alpha = 0.4) +
  geom_line(aes(y = height_fit), color = "blue") +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "soil salinity (mS)", y = "plant height (cm)")

# fits follow the data pretty well. leaves get thicker in all 3,
# Conocarpus gets shorter, Sesuvium/Laguncularia height pretty flat.
# Sesuvium at 4 mS is really spread out, model only gets the average.


# Overall -------------------------------------------------------------------------
#
# Saltier soil -> thicker leaves and mostly shorter plants, even within
# species (not just species turnover). Species respond differently though -
# Sesuvium changes the most, mangroves barely, some grasses go the other way.
# M5 (species-specific slopes) best for both traits.
#