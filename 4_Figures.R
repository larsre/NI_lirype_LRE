
##################################################

load("data/BES2019_BaseModel_varS_2007_2019.RData")


##################################################
## Plotting density; 

means <- out_real2$mean$D[2:13]
lower <- out_real2$q2.5$D[2:13]
upper <- out_real2$q97.5$D[2:13]
Year <- seq(2008,2019)

dat1 <- as_tibble(data.frame(means=means, lower=lower, upper=upper, Year=Year))

p1 <- ggplot(dat1, aes())

p1 <- ggplot(dat1, aes(x=Year, y=means)) + 
  ylim(0,50) +
  xlab("Year") +
  ylab("Density (95% CI)")+
  geom_point(size=2, color="chartreuse4")+
  geom_line(color="chartreuse4", size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, lwd=1, color="chartreuse4")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12))+
  scale_x_continuous(breaks=c(2007,2010,2013,2016, 2019))

p1 

ggsave("figurer/Figur1.png", p1, ,  bg = "transparent")


##################################################
## Plotting density; 

means <- out_real2$mean$S[2:13]
lower <- out_real2$q2.5$S[2:13]
upper <- out_real2$q97.5$S[2:13]
Year <- seq(2008,2019)

dat1 <- as_tibble(data.frame(means=means, lower=lower, upper=upper, Year=Year))

p1 <- ggplot(dat1, aes())

p1 <- ggplot(dat1, aes(x=Year, y=means)) + 
  ylim(0,1) +
  xlab("Year") +
  ylab("Survival (95% CI)")+
  geom_point(size=2, color="chartreuse4")+
  geom_line(color="chartreuse4", size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, lwd=1, color="chartreuse4")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12))+
  scale_x_continuous(breaks=c(2007,2010,2013,2016, 2019))

p1 

ggsave("figurer/Figur2.png", p1, ,  bg = "transparent")


##################################################
## Plotting R; 

means <- out_real2$mean$R_year[2:13]
lower <- out_real2$q2.5$R_year[2:13]
upper <- out_real2$q97.5$R_year[2:13]
Year <- seq(2008,2019)

dat1 <- as_tibble(data.frame(means=means, lower=lower, upper=upper, Year=Year))



p1 <- ggplot(dat1, aes(x=Year, y=means)) + 
  ylim(0,10) +
  xlab("Year") +
  ylab("Recruitment (95% CI)")+
  geom_point(size=2, color="chartreuse4")+
  geom_line(color="chartreuse4", size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, lwd=1, color="chartreuse4")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12))+
  scale_x_continuous(breaks=c(2007,2010,2013,2016, 2019))

p1 

ggsave("figurer/Figur3.png", p1,  bg = "transparent")

#####################################################################
## 

## Prior on S based on radiotelemetry study;
shape_from_stats <- function(mu , sigma ){
  a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
  b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
  shape_ps <- c(a,b)
  return(shape_ps)
}

S_priors <- shape_from_stats(0.425, 0.035)

## Plott - to see if this works; 

x=seq(0,1,.001)
w=dbeta(x,S_priors[1],S_priors[2])

df <- data.frame(x,w)

p10 <- ggplot(data=df, aes(x=x, y=w))+
      geom_line(colour="chartreuse4", lwd=1.5)+
      xlab("Survival")+
      ylab("betadensity")+
      theme_minimal()

p10

ggsave("figurer/Figur_prior_survival.png", p10,  bg = "transparent")


################################################################################
###

obs_map2 <- obs_map %>% filter(LinjeAvstand>=0)

p11 <- ggplot(data=obs_map2, aes(x=LinjeAvstand)) +
geom_histogram(aes(y=..density..), colour="grey", fill="chartreuse4", binwidth=2)+
  #geom_density(alpha=.5, fill="chartreuse4")+
  xlab("Distance to transect line")+
  ylab("number of obs")+
  theme_minimal()

p11

ggsave("figurer/Figur_Detection.png", p11,  bg = "transparent")

################################################################################

## Plotting Breeding success;
load("data/BES2019_BaseModel_commonS_2015_2019_Nest_Prior.RData")



means <- out_real7$mean$R.nest_year
lower <- out_real7$q2.5$R.nest_year
upper <- out_real7$q97.5$R.nest_year
Year <- seq(2015,2019)

dat1 <- as_tibble(data.frame(means=means, lower=lower, upper=upper, Year=Year))



p7 <- ggplot(dat1, aes(x=Year, y=means)) + 
  ylim(0,10) +
  xlab("Year") +
  ylab("Breeding success")+
  geom_point(size=2, color="chartreuse4")+
  geom_line(color="chartreuse4", size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, lwd=1, color="chartreuse4")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12))

p7 

ggsave("figurer/Figur7.png", p7,  bg = "transparent")

########################################################################
#### juvenile survival

means <- out_real7$mean$S_juv
lower <- out_real7$q2.5$S_juv
upper <- out_real7$q97.5$S_juv
Year <- seq(2015,2019)

dat1 <- as_tibble(data.frame(means=means, lower=lower, upper=upper, Year=Year))



p8 <- ggplot(dat1, aes(x=Year, y=means)) + 
  ylim(0,1) +
  xlab("Year") +
  ylab("Juvenile summer survival")+
  geom_point(size=2, color="chartreuse4")+
  geom_line(color="chartreuse4", size=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, lwd=1, color="chartreuse4")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12))

p8 

ggsave("figurer/Figur8.png", p8,  bg = "transparent")









