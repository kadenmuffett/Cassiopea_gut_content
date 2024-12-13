library(vegan)
library(car)
library(stats)
library(ggfortify)
library(rcompanion)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggstatsplot)
library(FSA)


#Load in##########
cassiopeastomachcontentswmeta<-read.csv("C:/Users/kmmuf/Desktop/Supplementary_Cass_Diet/Combined_diet_07.13.2024.csv")
######dataprep for vegan########
stomvegan<-cassiopeastomachcontentswmeta[,-c(2:11)]
rownames(stomvegan)<-(stomvegan$Specimen.)
stomvegan$Specimen.<-as.factor(stomvegan$Specimen.)
stomvegmin1<-stomvegan[,-1]
stomvegmin1<-stomvegmin1[-c(102:130),-c(23:27)]
#####Prepping version with no diet NAs########
cassiopeastomachcontentswmeta_nona<-cassiopeastomachcontentswmeta[-c(102:130),]
cassiopeastomachcontentswmeta_nona$Simp<-vegan::diversity(stomvegmin1, index = "simpson")
cassiopeastomachcontentswmeta_nona$Shan<-vegan::diversity(stomvegmin1, index = "shannon")
cassiopeastomachcontentswmeta_nona$Sex<-as.factor(cassiopeastomachcontentswmeta_nona$Sex)
cassiopeastomachcontentswmeta_nona$Site_details<-as.factor(cassiopeastomachcontentswmeta_nona$Site_details)
cassiopeastomachcontentswmeta_nona$Location<-as.factor(cassiopeastomachcontentswmeta_nona$Location)
####Linear regression explaining diversity########
globalmodel <- lm(Simp ~ Location + Time_collected + Diameter_cm + Sex + Site_details , data = cassiopeastomachcontentswmeta_nona, family=binomial(logit),
                  na.action = "na.fail")

combinations <- dredge(globalmodel, na.omit)

print(combinations)
hist(cassiopeastomachcontentswmeta_nona$Simp)
simpsonlocaov<- aov(Simp ~ Location,data= cassiopeastomachcontentswmeta_nona)
qqPlot(simpsonlocaov)
shapiro.test(simpsonlocaov$residuals)
boxplot(Simp ~ Location,
        data = cassiopeastomachcontentswmeta_nona)
kruskal.test(Simp ~ Location,
             data = cassiopeastomachcontentswmeta_nona)


dunnTest(Simp ~ Location,
         data = cassiopeastomachcontentswmeta_nona,
         method = "holm"
)

#Shannon diversity plot by location##########
hist(cassiopeastomachcontentswmeta_nona$Shan)
Shansonlocaov<- aov(Shan ~ Location,data= cassiopeastomachcontentswmeta_nona)
qqPlot(Shansonlocaov)
plot(Shansonlocaov)
shapiro.test(Shansonlocaov$residuals)
boxplot(Shan ~ Location,
        data = cassiopeastomachcontentswmeta_nona)
kruskal.test(Shan ~ Location,
             data = cassiopeastomachcontentswmeta_nona)
dunnTest(Shan ~ Location,
         data = cassiopeastomachcontentswmeta_nona,
         method = "holm"
)

### GGStatPlot of Shannon diversity by Location##########
om<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = Shan,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
) + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Shannon Diversity Index")
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/om.png",plot= om,height = 5,width = 6,units = "in")


print(combinations)
hist(cassiopeastomachcontentswmeta_nona$Shan)
shan2locaov<- aov(Shan ~ Site_details,data= cassiopeastomachcontentswmeta_nona)
qqPlot(shan2locaov)
shapiro.test(shan2locaov$residuals)

#Shannon by Site type model testing####
leveneTest(Shan ~ Site_details,
           data = cassiopeastomachcontentswmeta_nona
)
boxplot(Shan ~ Site_details,
        data = cassiopeastomachcontentswmeta_nona)
kruskal.test(Shan ~ Site_details,
             data = cassiopeastomachcontentswmeta_nona)
cassiopeastomachcontentswmeta_nona$Site_details<-revalue(cassiopeastomachcontentswmeta_nona$Site_details, c("Sandy/Silty" = "Sandy", "Sandy (near-dock)" = "Sandy", "Sandy (offshore)" = "Sandy"))

ob<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Site_details,
  y = Shan,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
) 
ob+ scale_colour_manual(values=c("brown4", "#E1BE6A","lightgoldenrod4"))+ ylab("Shannon Diversity Index")
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ob.png",plot= ob,height = 5,width = 6,units = "in")

### Analysis of variance of Shan by Site_details * Location#########
aov2<-aov(Shan ~ Site_details * Location , data= cassiopeastomachcontentswmeta_nona)
par(mfrow=c(2,2))
plot(aov2)
leveneTest(aov2)
summary(aov2)
Anova(aov2, type = "II")

#Crustacean of interest AOV#######
final<-aov(ARTH_Harpacticoid_Copepod/Total_prey_items~Location, data=cassiopeastomachcontentswmeta_nona)
leveneTest(final)
boxplot(ARTH_Harpacticoid_Copepod/Total_prey_items ~ Location,
        data = cassiopeastomachcontentswmeta_nona)
kruskal.test(ARTH_Harpacticoid_Copepod/Total_prey_items ~ Location,
             data = cassiopeastomachcontentswmeta_nona)
cassiopeastomachcontentswmeta_nona$prop_harp<-cassiopeastomachcontentswmeta_nona$ARTH_Harpacticoid_Copepod/cassiopeastomachcontentswmeta_nona$Total_prey_items

### Proportion harpacticoid crustacean by location#####

ombi<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = prop_harp,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE
) + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Harpacticoid proportoion by medusa")
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombi_harpact.png",plot= ombi,height = 5,width = 6,units = "in")

### Amount of microplastic by location########

ombi2<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = Microplastic,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
) + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Microplastic load by medusa")
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombi2_microplastics.png",plot= ombi2,height = 5,width = 6,units = "in")

### Proportion Crustacean by location###########

cassiopeastomachcontentswmeta_nona$arthprop<-rowSums(cassiopeastomachcontentswmeta_nona[,12:27])
cassiopeastomachcontentswmeta_nona$artpct<-cassiopeastomachcontentswmeta_nona$arthprop/cassiopeastomachcontentswmeta_nona$Total_prey_items
cassiopeastomachcontentswmeta_nona$copprop<-(rowSums(cassiopeastomachcontentswmeta_nona[,c(12:14,17,18)]))/cassiopeastomachcontentswmeta_nona$Total_prey_items
ombi3<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = artpct,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE
) + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Proportion arthropod by medusa")
ombi3
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombi3_arthropods.png",plot= ombi3,height = 5,width = 6,units = "in")

ombi4<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = copprop,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE
) + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Proportion copepod by medusa")
ombi4

ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombi4.png",plot= ombi4,height = 5,width = 6,units = "in")


plots <- align_plots(p3, p1, align = 'v', axis = 'l')
bottom_row <- plot_grid(ombi, ombi4, labels = c('B', 'C'), label_size = 12)

mega<-plot_grid(ombi3, bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/mega.png",plot= mega,height = 10,width = 6,units = "in")

#Diameter + Location + Site############
#Diameter + Location
partstom<-lm(Shan~cassiopeastomachcontentswmeta_nona$Diameter_cm+
               cassiopeastomachcontentswmeta_nona$Site_details+cassiopeastomachcontentswmeta_nona$Location,data=cassiopeastomachcontentswmeta_nona)
diaandsite<-lm(Shan~cassiopeastomachcontentswmeta_nona$Diameter_cm
               +cassiopeastomachcontentswmeta_nona$Location,data=cassiopeastomachcontentswmeta_nona)
#Distribution of alpha diversity measures##########
ggplot(cassiopeastomachcontentswmeta_nona, aes(x=Shan)) + geom_histogram()
ggplot(cassiopeastomachcontentswmeta_nona, aes(x=(log(Shan)))) + geom_histogram()

autoplot(diaandsite)
autoplot(partstom)
##Comparison
compareLM(diaandsite,partstom)
##Selected model
anova(partstom)
summary(partstom)
####Site importance
print(summary(aov(Shan~Location, cassiopeastomachcontentswmeta_nona)))
print(summary(aov(Shan~Location*Site_details, cassiopeastomachcontentswmeta_nona)))

################# Total prey by location#############
cass_positive_prey<-subset(cassiopeastomachcontentswmeta_nona,cassiopeastomachcontentswmeta_nona$Total_prey_items>0)
cass_positive_prey$logcons<-log(cass_positive_prey$Total_prey_items)
log_tot_prey<-lm(logcons~Diameter_cm, data=cass_positive_prey)
tot_prey<-lm(Total_prey_items~Diameter_cm, data=cass_positive_prey)
autoplot(log_tot_prey)
autoplot(tot_prey)
#Logged to improve normality
#Data summary
anova(lm(logcons~Diameter_cm, data=cass_positive_prey))
summary(lm(logcons~Diameter_cm, data=cass_positive_prey))
#Graphing location and prey
xdata4<-c(5:25)
ydata4<-c(.08347*xdata4+1.16835)
tf <- ggplot(data = cass_positive_prey, aes(x = Diameter_cm, y = logcons))+
  geom_point(aes(color = Location, shape = Location), size=2)

tf + scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4")) + 
  xlab("Diameter (cm)") + ylab("Total prey items (log10)") + theme_minimal() +
  geom_line(data = data.frame(x = xdata4, y = ydata4),
            aes(x = x, y = y), col = "black", lwd = 1) 

#####
# ######## Aurelia and Cassiopea diet ##############
# diet_o_scyph<-read.csv("C:/Users/kmmuf/Desktop/Supplementary_Cass_Diet/Aurelia_Cassiopea_diet.csv")
# diet_o_scyph$logcons<-log(diet_o_scyph$Total_prey)
# autoplot(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Cassiopea")))
# autoplot(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Aurelia")))
# summary(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Cassiopea")))
# summary(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Aurelia")))
# Pan_diet_o_scyph<-diet_o_scyph[c(44:95),]
# summary(lm(Pan_diet_o_scyph$Tax~Pan_diet_o_scyph$Diam, data= Pan_diet_o_scyph))
# 
# xdata2<-c(1:35)
# ydata2<-c((.33498*xdata2-.7485))
# 
# ydata3<-c((.1701*xdata2+1.2052))    
# tax <- ggplot(data = diet_o_scyph, aes(x = Diam, y = Tax))+
#        geom_point(aes(color = Group, shape = Group), size=2)
#      tax + scale_colour_manual(name="Group", values=c("coral", "black")) + 
#        xlab("Medusa Diameter (cm)") + ylab("Taxa in gut") + theme_minimal() +
#        geom_line(data = data.frame(x = xdata2, y = ydata2),
#                  aes(x = x, y = y), col = "coral", lwd = 1) +
#        geom_line(data = data.frame(x = xdata2, y = ydata3),
#                  aes(x = x, y = y), col = "black", lwd = 1)
#      
############
     #Sex identification#############
     p<-ggplot(cassiopeastomachcontentswmeta_nona) +
       aes(x = Location, fill = factor(Sex)) +
       geom_bar(position = "fill") + theme_bw() + scale_fill_grey(start = 0, end = .9)
     p
     
#########
     #Stat bar for violin plots######
     data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
     }
########## Violin plots###############
     #Shannon diversity by location
     f<-ggplot(cassiopeastomachcontentswmeta_nona, aes(x=Location, y=Shan, fill=Location)) +
       geom_violin(position=position_dodge(.5),trim=FALSE)
     f + geom_point(position=position_jitterdodge(dodge.width=0.5),color="black") + 
       stat_summary(fun.data=data_summary,position=position_dodge(.5),color="white") + 
       xlab("Location") + ylab("Shannon Diversity")  + theme_minimal() +scale_fill_manual(values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))
     #####
     #Shannon diversity by site type
     i<-ggplot(cassiopeastomachcontentswmeta_nona, aes(x=Site_details, y=Shan, fill=Site_details)) +
       geom_violin(position=position_dodge(.5),trim=FALSE)
     i + geom_jitter(height = 0, width = 0.1, aes(colour = factor(Location))) +
       scale_colour_manual(name="Colour", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4")) + 
       stat_summary(fun.data=data_summary,position=position_dodge(.5),color="black") + 
       xlab("Location Type") + ylab("Shannon Diversity") + theme_minimal() +scale_fill_manual(values=c( "lightyellow","bisque3","darkgreen","pink","orange"))
     ####
     #
     cassiopeastomachcontentswmeta$Sexf<-as.factor(cassiopeastomachcontentswmeta$Sex)
     #remove individuals lacking weight
     cassiopeastomachcontentswmeta_noAustralia<-cassiopeastomachcontentswmeta[-c(102:107,114:130),]
     ##Mass by diameter comparisons
     #Identified best fit
     xdata3<-c(1:25)
     y_3<-c(.0557*(xdata3^3.1))
     kl <- ggplot(data = cassiopeastomachcontentswmeta_noAustralia, aes(x = Diameter_cm, y = Weight_g))+
       geom_point(aes(color = Location, shape = Sex), size=2)
     kl + scale_colour_manual(name="Colour", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4")) + 
       xlab("Medusa Diameter (cm)") + ylab("Mass (g)") + theme_minimal() +
       geom_line(data = data.frame(x = xdata3, y = y_3), aes(x = x, y = y), col = 1, lwd = 1)
 
     ############
     #Diameter and Latitude
     ggplot(cassiopeastomachcontentswmeta, aes(x=Lat, y=Diameter_cm,group=Location))+ theme_minimal() +
       geom_boxplot(notch=FALSE, outlier.shape=NA, fill=("grey"), color = c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4","purple"), alpha=0.2)
     
     ##Diameter of medusae by collection locale
     boxy<- ggboxplot(cassiopeastomachcontentswmeta, "Location", "Diameter_cm",
               color = "Location", palette =c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4","purple"),
               add = "jitter", shape = "Location") + theme_minimal() + xlab("Location") + ylab("Medusa Diameter (cm)") 
     boxy
     
     ############
     #Rhopalia###############
     hasrhopalia<-subset(cassiopeastomachcontentswmeta,cassiopeastomachcontentswmeta_nona$Rhopalia.>1)
     likelyCassxam<-subset(hasrhopalia,hasrhopalia$Location!="Philippines")
     rhopalia<-likelyCassxam[,9]
     
     rhopalia<- as.data.frame(rhopalia
     )

     ggplot(rhopalia, aes(x=rhopalia)) + 
       geom_histogram(aes(y=..density..), binwidth = 1, colour="black", fill="white")+
       geom_density(alpha=.2, fill="#FF6666") + theme_minimal() + xlab("Number of Rhopalia") + ylab("Density")
     
     correlation.one <- cor.test(log10(cassiopeastomachcontentswmeta_nona$Total_prey_items+0.001), cassiopeastomachcontentswmeta_nona$Diameter_cm, method = 'pearson')
     correlation.one
####
str(stomvegmin1)
plot(ARTH_Harpacticoid_Copepod/Total_prey_items~Diameter_cm, data=cassiopeastomachcontentswmeta_nona)

#### Plankton group based comparisons ##########
cassiopeastomachcontentswmeta_nona$benthos<-(rowSums(cassiopeastomachcontentswmeta_nona[,c(12,21,23,26,28,32)])/cassiopeastomachcontentswmeta_nona$Total_prey_items)
          
ombiben<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = benthos,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE
) + ylim(0,1.25)+ scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Proportion benthic diet")
ombiben
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombiben.png",plot= ombiben,height = 5,width = 6,units = "in")

cassiopeastomachcontentswmeta_nona$zoo<-(rowSums(cassiopeastomachcontentswmeta_nona[,c(13,14,24,29)])/cassiopeastomachcontentswmeta_nona$Total_prey_items)
ombizoo<-ggbetweenstats(
  data = cassiopeastomachcontentswmeta_nona,
  x = Location,
  y = zoo,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE
) +ylim(0,1.25)+ scale_colour_manual(name="Location", values=c("lightblue", "coral3", "cadetblue", "#E1BE6A","lightgoldenrod4"))+ ylab("Proportion zooplankton diet")
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/ombizoo.png",plot= ombizoo,height = 5,width = 6,units = "in")
zooben<-plot_grid(ombizoo,ombiben,labels = c('A', 'B'), label_size = 12)
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/benzoo.png",plot= zooben,height = 5,width = 11,units = "in")
summary(cassiopeastomachcontentswmeta_nona$benthos)
summary(cassiopeastomachcontentswmeta_nona$zoo)
mega3<-plot_grid(ombizoo,ombiben,labels=c('D','E'),label_size = 12)
mega2<-plot_grid(mega, mega3, ncol = 1, rel_heights = c(2,1))
ggsave("C:/Users/kmmuf/Documents/Panama Gut Content manuscript/mega2_means.png",plot= mega2,height = 14,width = 11,units = "in")
