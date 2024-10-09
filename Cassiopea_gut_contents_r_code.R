library(vegan)
library(car)
library(stats)
library(ggfortify)
library(rcompanion)
library(MuMIn)
library(ggplot2)
library(ggpubr)


#Load in
cassiopeastomachcontentswmeta<-read.csv("C:/Users/kmmuf/Desktop/Supplementary_Cass_Diet/Combined_diet_07.13.2024.csv")
######dataprep for vegan
stomvegan<-cassiopeastomachcontentswmeta[,-c(2:11)]
stomvegan$Specimen.<-as.factor(stomvegan$Specimen.)
stomvegmin1<-stomvegan[,-1]
stomvegmin1<-stomvegmin1[-c(102:130),-c(23:27)]
#####Prepping version with no diet NAs
cassiopeastomachcontentswmeta_nona<-cassiopeastomachcontentswmeta[-c(102:130),]
cassiopeastomachcontentswmeta_nona$Simp<-diversity(stomvegmin1,index = "simpson")
cassiopeastomachcontentswmeta_nona$Shan<-diversity(stomvegmin1,index = "shannon")
cassiopeastomachcontentswmeta_nona$Sex<-as.factor(cassiopeastomachcontentswmeta_nona$Sex)
cassiopeastomachcontentswmeta_nona$Site_details<-as.factor(cassiopeastomachcontentswmeta_nona$Site_details)
cassiopeastomachcontentswmeta_nona$Location<-as.factor(cassiopeastomachcontentswmeta_nona$Location)
####Linear regression explaining diversity
globalmodel <- lm(Simp ~ Location + Time_collected + Diameter_cm + Sex + Site_details , data = cassiopeastomachcontentswmeta_nona, family=binomial(logit),
                  na.action = "na.fail")

combinations <- dredge(globalmodel, na.omit)

print(combinations)
#Diameter + Location + Site
#Diameter + Location
partstom<-lm(Shan~cassiopeastomachcontentswmeta_nona$Diameter_cm+
               cassiopeastomachcontentswmeta_nona$Site_details+cassiopeastomachcontentswmeta_nona$Location,data=cassiopeastomachcontentswmeta_nona)
diaandsite<-lm(Shan~cassiopeastomachcontentswmeta_nona$Diameter_cm
               +cassiopeastomachcontentswmeta_nona$Location,data=cassiopeastomachcontentswmeta_nona)
#Distribution of alpha diversity measures
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

################# Total prey by location
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
######## Aurelia and Cassiopea
diet_o_scyph<-read.csv("C:/Users/kmmuf/Desktop/Supplementary_Cass_Diet/Aurelia_Cassiopea_diet.csv")
diet_o_scyph$logcons<-log(diet_o_scyph$Total_prey)
autoplot(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Cassiopea")))
autoplot(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Aurelia")))
summary(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Cassiopea")))
summary(lm(diet_o_scyph$Tax~diet_o_scyph$Diam, data= diet_o_scyph, subset=(Group=="Aurelia")))
Pan_diet_o_scyph<-diet_o_scyph[c(44:95),]
summary(lm(Pan_diet_o_scyph$Tax~Pan_diet_o_scyph$Diam, data= Pan_diet_o_scyph))

xdata2<-c(1:35)
ydata2<-c((.33498*xdata2-.7485))

ydata3<-c((.1701*xdata2+1.2052))    
tax <- ggplot(data = diet_o_scyph, aes(x = Diam, y = Tax))+
       geom_point(aes(color = Group, shape = Group), size=2)
     tax + scale_colour_manual(name="Group", values=c("coral", "black")) + 
       xlab("Medusa Diameter (cm)") + ylab("Taxa in gut") + theme_minimal() +
       geom_line(data = data.frame(x = xdata2, y = ydata2),
                 aes(x = x, y = y), col = "coral", lwd = 1) +
       geom_line(data = data.frame(x = xdata2, y = ydata3),
                 aes(x = x, y = y), col = "black", lwd = 1)
     
############
     #Sex identification
     p<-ggplot(cassiopeastomachcontentswmeta_nona) +
       aes(x = Location, fill = factor(Sex)) +
       geom_bar(position = "fill") + theme_bw() + scale_fill_grey(start = 0, end = .9)
     p
     
#########
     #Stat bar for violin plots
     data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
     }
##########
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
     ##Mass by diameter
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
     #Rhopalia
     hasrhopalia<-subset(cassiopeastomachcontentswmeta,cassiopeastomachcontentswmeta_nona$Rhopalia.>1)
     likelyCassxam<-subset(hasrhopalia,hasrhopalia$Location!="Philippines")
     rhopalia<-likelyCassxam[,9]
     
     rhopalia<- as.data.frame(rhopalia
     )

     ggplot(rhopalia, aes(x=rhopalia)) + 
       geom_histogram(aes(y=..density..), binwidth = 1, colour="black", fill="white")+
       geom_density(alpha=.2, fill="#FF6666") + theme_minimal() + xlab("Number of Rhopalia") + ylab("Density")
     
     