##### Fat content model###########
######Estonian winters 

rm(list=ls(all=TRUE)); # graphics.off()

#packages
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)

#data 
Fatt <-read.csv("FatContent.csv", header = T, stringsAsFactors= T)

# check the levels
levels(Fatt$Treatment)

#rescale the response variable, because it needs to be 0 to 1 
Fatt$Fatpercentage <- Fatt$Fatpercentage/ 100 

#GLMM 
model_ <- glmmTMB(Fatpercentage ~ Treatment * Opening + (1|ID),
                  data = Fatt,
                  tweedie())

# Display the summary of the model
summary(model_)

#check assumtions 
#DHARMA test
obj <- simulateResiduals(model_, plot = F)
plot(obj, quantreg = F)
Anova(model_, type="III")

#post hoc posthoc with Tukey correction for the chosen model 
emm <- emmeans(model_, specs = pairwise ~ Opening | Treatment)
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")

emm <- emmeans(model_, specs = pairwise ~ Treatment | Opening)
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")

###plot the model ###############################################
#22.04.24

# Extract unique treatment levels
treatment_levels <- unique(Fatt$Treatment)

# Create a grid of factors for prediction
pred_data <- expand.grid(
  Opening = factor(c("Dec", "Feb"), levels = c("Dec", "Feb")),
  Treatment = factor(c("Cold", "Moderate", "Warm"), levels = c("Cold", "Moderate", "Warm"))
)

# Generate predictions with confidence intervals
predictions <- predict(model_, newdata = pred_data, type = "response")
conf_int <- predict(model_, newdata = pred_data, type = "response", se.fit = TRUE)

# Include the predictions and confidence intervals in pred_data
pred_data$Fatpercentage <- predictions
pred_data$lower <- predictions - 1.96 * conf_int$se.fit  # 95% CI lower bound
pred_data$upper <- predictions + 1.96 * conf_int$se.fit  # 95% CI upper bound
pred_data$DataType <- 'Predicted'

# Transform 'Opening' to numeric in pred_data for consistency
pred_data$Opening_numeric <- as.numeric(pred_data$Opening == "Feb") + 1

# Prepare observed data with necessary columns
observed_data <- Fatt %>%
  mutate(Opening_numeric = case_when(Opening == "Dec" ~ 1, Opening == "Feb" ~ 2),
         DataType = 'Observed',
         lower = NA,  # Not applicable for observed data
         upper = NA)  # Not applicable for observed data

# Ensure both datasets have the same columns, in the same order
if (!"Treatment" %in% names(observed_data)) {
  observed_data$Treatment <- NA  # Assign NA or an appropriate value
}

# Ensure both datasets have identical columns for combining
columns_to_use <- c("Opening_numeric", "Fatpercentage", "DataType", "lower", "upper", "Treatment")

# Add missing columns with NA values where necessary
all_columns <- unique(c(names(observed_data), columns_to_use))
for (col in all_columns) {
  if (!col %in% names(observed_data)) {
    observed_data[[col]] <- NA
  }
  if (!col %in% names(pred_data)) {
    pred_data[[col]] <- NA
  }
}

# Now use rbind with the columns in order
combined_data <- rbind(
  observed_data[columns_to_use],
  pred_data[columns_to_use]
)


# Custom colors
custom_colors <- c("#0072b2", "#009e73", "#d55e00") 

# Plotting the combined data with confidence interval ribbons
ggplot(combined_data, aes(x = Opening_numeric, y = Fatpercentage, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = subset(combined_data, DataType == 'Predicted'),
              aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = subset(combined_data, DataType == 'Observed'),
              aes(x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), # Shift points
              width = 0.05, height = 0, shape = 16, alpha = 0.8) +
  geom_line(data = subset(combined_data, DataType == 'Predicted'), linetype = "solid", size = 1) +
  stat_summary(data = subset(combined_data, DataType == 'Observed'), 
               aes(group = Treatment, x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), 
               fun = median, 
               geom = "point", 
               shape = 22,  # Square shape
               size = 4, 
               position = position_dodge(width = 0.1)) +  # Adjust dodge width if needed
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Mid-Winter", "Post-Winter"))+
  labs(x = NULL, y = "Fat content", color = "Treatment", fill = "Treatment") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black")) +
  ylim(0.03, 0.16)  # Set y-axis limits from 0 to 0.17
ggsave("Fat_content.tiff", width=4, height=6)


