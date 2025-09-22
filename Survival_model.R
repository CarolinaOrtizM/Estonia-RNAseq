#######Survival proportion ########
##Estonian winter
#

rm(list=ls(all=TRUE)); # graphics.off()

#packages
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)

####data
survall <-read.csv("Survival_Alleggsacs.csv", header = T,stringsAsFactors=T) 
survall <- survall[-115,] #take out the row that has no hatching success 

###subsetting 

Egg1 = subset(survall, survall$Eggsac == "1") 
Egg2 = subset(survall, survall$Eggsac == "2")

#survival/100

survall$survivalrate <- survall$survivalrate/100

# Adjusting survival rates slightly to avoid exact 0 or 1
survall$survivalrate_adj <- pmin(pmax(survall$survivalrate, 0.001), 0.999)

# Applying a logit transformation to the adjusted survival rates
survall$logit_survivalrate <- log(survall$survivalrate_adj / (1 - survall$survivalrate_adj))

shapiro.test(survall$logit_survivalrate)
#not normal
leveneTest(survivalrate ~ Treatment, data = survall)

#GLMM model 

TMBmodelS<- glmmTMB(survivalrate ~ Opening + Treatment + scale(Age) + scale(clutchsize) + (1|ID), 
                    family = ordbeta(), 
                    data = survall)


#### testing the model ####
emmeans(TMBmodelS, pairwise ~ Opening)


Tweedie_model<-glmmTMB(logit_survivalrate ~ Opening + Treatment + scale(Age) + scale(clutchsize) + (1|ID), 
                       family = tweedie(), 
                       data = survall)

AIC(TMBmodelS, Tweedie_model, Eggs, Poly, Age)  # ran a comparison


Eggs<- glmmTMB(survivalrate ~ Opening + Treatment + scale(clutchsize), 
                    family = ordbeta(), 
                    data = Egg2)


Poly<- glmmTMB(logit_survivalrate ~ Opening * Treatment  + poly(scale(Age), 2),
        family = tweedie(), 
        data = survall)

Age <- glmmTMB(logit_survivalrate ~ Opening * Treatment * scale(Age) + (1|ID), 
               family = tweedie(), 
               data = survall)

summary(TMBmodelS)

#####DHARMA test####
obj <- simulateResiduals(TMBmodelS, plot = F)
plot(obj, quantreg = F)
Anova(TMBmodelS, type="III")

#post hoc posthoc with Tukey correction for the chosen model 
emm <- emmeans(TMBmodelS, specs = pairwise ~ Opening )
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")

emm <- emmeans(TMBmodelS, specs = pairwise ~ Age)
posthoc_results <- summary(emm$contrasts)
summary(emm$contrasts, adjust = "tukey")


######plot model

model_surv <- glmmTMB(survivalrate ~ Opening + Treatment, 
                      family = ordbeta(), 
                      data = survall)

# Prepare prediction data
treatment_levels <- unique(survall$Treatment)
pred_data <- expand.grid(
  Opening = c("mid-winter", "post-winter"),
  Treatment = treatment_levels
)

# Encode Opening as numeric
pred_data$Opening_numeric <- ifelse(pred_data$Opening == "mid-winter", 1, 2)

# Generate logit predictions with confidence intervals
predictions <- predict(model_surv, newdata = pred_data, type = "link")
conf_int <- predict(model_surv, newdata = pred_data, type = "link", se.fit = TRUE)

# Convert logit to survival rate using the logistic function
pred_data$survivalrate <- exp(predictions) / (1 + exp(predictions))
pred_data$lower <- exp(predictions - 1.96 * conf_int$se.fit) / (1 + exp(predictions - 1.96 * conf_int$se.fit))  # 95% CI lower bound
pred_data$upper <- exp(predictions + 1.96 * conf_int$se.fit) / (1 + exp(predictions + 1.96 * conf_int$se.fit))  # 95% CI upper bound

# Add a DataType column for plotting distinction
pred_data$DataType <- 'Predicted'

# Prepare observed data with the necessary columns
observed_data <- data.frame(
  Opening_numeric = ifelse(survall$Opening == "mid-winter", 1, 2),
  survivalrate = survall$survivalrate,
  lower = NA,  # Not applicable for observed data
  upper = NA,  # Not applicable for observed data
  DataType = 'Observed',
  Treatment = survall$Treatment
)

# Combine observed and predicted data for plotting
combined_data <- rbind(
  observed_data,
  pred_data[, c("Opening_numeric", "survivalrate", "DataType", "lower", "upper", "Treatment")]
)


# Custom colors
custom_colors <- c("#0072b2", "#009e73", "#d55e00")


# Plotting the combined data with confidence interval ribbons
Survival <-ggplot(combined_data, aes(x = Opening_numeric, y = survivalrate, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = subset(combined_data, DataType == 'Predicted'),
              aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_jitter(data = subset(combined_data, DataType == 'Observed'), width=0.05, shape = 19, alpha = 0.8) +
  geom_line(data = subset(combined_data, DataType == 'Predicted'), linetype = "solid", size = 1) +
  stat_summary(data = subset(combined_data, DataType == 'Observed'), 
               aes(group = Treatment, x = ifelse(Opening_numeric == 1, Opening_numeric - 0.05, Opening_numeric + 0.05)), 
               fun = median, 
               geom = "point", 
               shape = 22,  # Square shape
               size = 4, 
               position = position_dodge(width = 0.1)) + 
  scale_color_manual(values = custom_colors) +  # Set custom colors for lines and points
  scale_fill_manual(values = custom_colors) +  # Set custom colors for ribbon fills
  scale_x_continuous(breaks = c(1, 2), labels = c("1 month", "3 months"))+
  labs(x = NULL, y = "Survival Proportion", color = "Treatment") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))
Survival
ggsave("survival_last.tiff", width=4, height=6) ##for saving it 
