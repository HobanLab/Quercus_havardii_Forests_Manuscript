#box plots of morphological variation

library(tibble);library(ggplot2); library(ggpubr)

gg_lobes <- ggplot(QH_total, aes(x = Pop_ID, y=Lobes)) + geom_boxplot() + xlab ("Population") + ylab ("Lobe Number") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_length <- ggplot(QH_total, aes(x = Pop_ID, y=Length)) + geom_boxplot() + xlab ("Population") + ylab ("Leaf Length (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_BLBW <- ggplot(QH_total, aes(x = Pop_ID, y=BLBW)) + geom_boxplot() + xlab ("Population") + ylab ("Basal Lobe Blade Width (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_MLBW <- ggplot(QH_total, aes(x = Pop_ID, y=MLBW)) + geom_boxplot() + xlab ("Population") + ylab ("Middle Lobe Blade Width (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_ALBW <- ggplot(QH_total, aes(x = Pop_ID, y=ALBW)) + geom_boxplot() + xlab ("Population") + ylab ("Apical Lobe Blade Width (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LVL <- ggplot(QH_total, aes(x = Pop_ID, y=LVL)) + geom_boxplot() + xlab ("Population") + ylab ("Lower Vein Length (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_UVL <- ggplot(QH_total, aes(x = Pop_ID, y=UVL)) + geom_boxplot() + xlab ("Population") + ylab ("Upper Vein Length (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_UMLS <- ggplot(QH_total, aes(x = Pop_ID, y=UMLS)) + geom_boxplot() + xlab ("Population") + ylab ("Upper Middle Lobe Sinus Width (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LMLS <- ggplot(QH_total, aes(x = Pop_ID, y=LMLS)) + geom_boxplot() + xlab ("Population") + ylab ("Lower Middle Lobe Sinus Width (mm)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LLA <- ggplot(QH_total, aes(x = Pop_ID, y=LLA)) + geom_boxplot() + xlab ("Population") + ylab ("Angle of Lower Lobe (°)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_ULA <- ggplot(QH_total, aes(x = Pop_ID, y=ULA)) + geom_boxplot() + xlab ("Population") + ylab ("Angle of Upper Lobe (°)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_MLA <- ggplot(QH_total, aes(x = Pop_ID, y=MLA)) + geom_boxplot() + xlab ("Population") + ylab ("Angle of Middle Lobe (°)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))

ggarrange (gg_lobes, gg_length, gg_BLBW, gg_MLBW, gg_ALBW, gg_LVL, gg_UVL, gg_UMLS, gg_LMLS, gg_LLA, gg_ULA, gg_MLA, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I","J", "K", "L"), ncol = 3, nrow = 4)

gg_LWR_1 <- ggplot(QH_total, aes(x = Pop_ID, y=LWR_1)) + geom_boxplot() + xlab ("Population") + ylab ("Length/width ratio 1 (BLBW/LENGTH × 100)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LWR_2 <- ggplot(QH_total, aes(x = Pop_ID, y=LWR_2)) + geom_boxplot() + xlab ("Population") + ylab ("Length/width ratio 2 (MLBW/LENGTH × 100)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LWR_3 <- ggplot(QH_total, aes(x = Pop_ID, y=LWR_3)) + geom_boxplot() + xlab ("Population") + ylab ("Length/width ratio 3 (ALBW/LENGTH × 100)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LODR_1 <- ggplot(QH_total, aes(x = Pop_ID, y=LODR_1)) + geom_boxplot() + xlab ("Population") + ylab ("Lobe depth ratio 1 (MLBW - UMLS/MLBW × 100)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))
gg_LODR_2 <- ggplot(QH_total, aes(x = Pop_ID, y=LODR_2)) + geom_boxplot() + xlab ("Population") + ylab ("Lobe depth ratio 2 (MLBW - LMLS/MLBW × 100)") + coord_flip() + theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8))

ggarrange (gg_LWR_1, gg_LWR_2, gg_LWR_3, gg_LODR_1, gg_LODR_2, labels = c("A", "B", "C", "D", "E"), ncol = 2, nrow = 3)


# for the full PCA - Bethany ran a PCA with averages more recently to make match the environmental and genetic graphs. 

total.pca <- prcomp(QH_total[,c(4:20)], center = TRUE,scale. = TRUE)

two_groups_full <- c(rep("West",627), rep("East",301))

ggbiplot(total.pca, ellipse = FALSE, obs.scale = 1.25, var.scale = 1.5, varname.adjust = 1.5, 
         var.axes = TRUE, labels= NULL, groups=two_groups_full)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=two_groups_full, color=two_groups_full, size=3))+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  xlab(expression("Axis 1 (22.9%)")) + 
  ylab(expression("Axis 2 (17.9%)"))+
  ggtitle("Morphological") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none") +
  scale_x_reverse()


#pearson correlation

cor(QH_total[c(4:15)], use="complete.obs", method="pearson")


#For nested ANOVA


library(nlme)
library(MuMIn)

length.lme <- lme(Length ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(length.lme))
qqline(residuals(length.lme))
summary(length.lme)
r.squaredGLMM(length.lme)

lobes.lme <- lme(Lobes ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(lobes.lme))
qqline(residuals(lobes.lme))
summary(lobes.lme)
r.squaredGLMM(lobes.lme)

BLBW.lme <- lme(BLBW ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(BLBW.lme))
qqline(residuals(BLBW.lme))
summary(BLBW.lme)
r.squaredGLMM(BLBW.lme)

MLBW.lme <- lme(MLBW ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(MLBW.lme))
qqline(residuals(MLBW.lme))
summary(MLBW.lme)
r.squaredGLMM(MLBW.lme)

ALBW.lme <- lme(ALBW ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(ALBW.lme))
qqline(residuals(ALBW.lme))
summary(ALBW.lme)
r.squaredGLMM(ALBW.lme)

LVL.lme <- lme(LVL ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LVL.lme))
qqline(residuals(LVL.lme))
summary(LVL.lme)
r.squaredGLMM(LVL.lme)

UVL.lme <- lme(UVL ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(UVL.lme))
qqline(residuals(UVL.lme))
summary(UVL.lme)
r.squaredGLMM(UVL.lme)

UMLS.lme <- lme(UMLS ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(UMLS.lme))
qqline(residuals(UMLS.lme))
summary(UMLS.lme)
r.squaredGLMM(UMLS.lme)

LMLS.lme <- lme(LMLS ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LMLS.lme))
qqline(residuals(LMLS.lme))
summary(LMLS.lme)
r.squaredGLMM(LMLS.lme)

LLA.lme <- lme(LLA ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LLA.lme))
qqline(residuals(LLA.lme))
summary(LLA.lme)
r.squaredGLMM(LLA.lme)

ULA.lme <- lme(ULA ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(ULA.lme))
qqline(residuals(ULA.lme))
summary(ULA.lme)
r.squaredGLMM(ULA.lme)

MLA.lme <- lme(MLA ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(MLA.lme))
qqline(residuals(MLA.lme))
summary(MLA.lme)
r.squaredGLMM(MLA.lme)

LWR_1.lme <- lme(LWR_1 ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LWR_1.lme))
qqline(residuals(LWR_1.lme))
summary(LWR_1.lme)
r.squaredGLMM(LWR_1.lme)

LWR_2.lme <- lme(LWR_2 ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LWR_2.lme))
qqline(residuals(LWR_2.lme))
summary(LWR_2.lme)
r.squaredGLMM(LWR_2.lme)

LWR_3.lme <- lme(LWR_3 ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LWR_3.lme))
qqline(residuals(LWR_3.lme))
summary(LWR_3.lme)
r.squaredGLMM(LWR_3.lme)

LODR_1.lme <- lme(LODR_1 ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LODR_1.lme))
qqline(residuals(LODR_1.lme))
summary(LODR_1.lme)
r.squaredGLMM(LODR_1.lme)

LODR_2.lme <- lme(LODR_2 ~ Region, random=list(~1|Pop_ID, ~1|Ind_ID),data=QH_total_ANOVA)
qqnorm(residuals(LODR_2.lme))
qqline(residuals(LODR_2.lme))
summary(LODR_2.lme)
r.squaredGLMM(LODR_2.lme)
