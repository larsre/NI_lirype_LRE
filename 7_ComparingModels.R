library(coda)
library(reshape2)
library(ggplot2)
library(viridis)

## Load model fits
modA <- readRDS('mod3b_versionA_realData.rds')
modB <- readRDS('mod3b_versionB_realData.rds')
modC <- readRDS('mod3b_versionC_realData.rds')

## Convert posterior samples to matrix format
out.modA <- as.matrix(modA)
out.modB <- as.matrix(modB)
out.modC <- as.matrix(modC)

## Set numbers of sites and years
N_sites <- 58
N_years <- 6

## Model A: Rename parameters
for(j in 1:N_sites){
  for(t in 1:N_years){
    
    dimnames(out.modA)[[2]][which(dimnames(out.modA)[[2]]==paste0('D1[', j, ', ', t, ']'))] <- paste0('Density[1, ', j, ', ', t, ']')
    dimnames(out.modA)[[2]][which(dimnames(out.modA)[[2]]==paste0('D2[', j, ', ', t, ']'))] <- paste0('Density[2, ', j, ', ', t, ']')
    
  }
}

## Models B and C: Sum N_exp and Density across age classes
for(j in 1:N_sites){
  for(t in 1:N_years){
    N_exp.B <- cbind(out.modB[,paste0('N_exp[1, ', j, ', ', t, ']')] + out.modB[,paste0('N_exp[2, ', j, ', ', t, ']')]) 
    dimnames(N_exp.B)[[2]] <- paste0('N_exp[', j, ', ', t, ']')
    Density.B <- cbind(out.modB[,paste0('Density[1, ', j, ', ', t, ']')] + out.modB[,paste0('Density[2, ', j, ', ', t, ']')]) 
    dimnames(Density.B)[[2]] <- paste0('Density[', j, ', ', t, ']')
    out.modB <- cbind(out.modB, N_exp.B, Density.B)
    
    N_exp.C <- cbind(out.modC[,paste0('N_exp[1, ', j, ', ', t, ']')] + out.modC[,paste0('N_exp[2, ', j, ', ', t, ']')]) 
    dimnames(N_exp.C)[[2]] <- paste0('N_exp[', j, ', ', t, ']')
    Density.C <- cbind(out.modC[,paste0('Density[1, ', j, ', ', t, ']')] + out.modC[,paste0('Density[2, ', j, ', ', t, ']')]) 
    dimnames(Density.C)[[2]] <- paste0('Density[', j, ', ', t, ']')
    out.modC <- cbind(out.modC, N_exp.C, Density.C)
  }
}

## Reformat posterior samples
dataA <- melt(out.modA)
dataA$Model <- 'Original'

dataB <- melt(out.modB)
dataB$Model <- 'Alt process'

dataC <- melt(out.modC)
dataC$Model <- 'Alt process + data'

data.all <- rbind(dataA, dataB, dataC)
colnames(data.all) <- c('Sample', 'Parameter', 'Value', 'Model')

## Identify 0 adundances/densities
#comb <- expand.grid(1:N_sites, 1:N_years)
#Ns <- paste0('N_exp[', comb[,1], ', ', comb[,2], ']')
#N0s <- Ns[which(colSums(out.modA[,Ns]) == 0)]
#jt0 <- matrix(unlist(stringr::str_extract_all(N0s, "\\d+")), ncol = 2, byrow = T)

## Set groups of parameters to plot together
mains <- c('mu.D1', 'sigma.D', 'mu.R', 'sigma.R',
           'Mu.S1', 'Mu.S2', 'ratio.JA1')

detects <- c(paste0('ews[', 1:N_years, ']'), paste0('p[', 1:N_years, ']'))

R_year <- paste0('R_year[', 1:N_years, ']')
D <- paste0('D[', 1:N_years, ']') 


## Plot posterior overlaps from different models

# Main parameters
pdf('Plots/ModelComp_Mains.pdf', width = 8, height = 5)
ggplot(subset(data.all, Parameter %in% mains)) + 
  geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# Detection parameters
pdf('Plots/ModelComp_Detects.pdf', width = 8, height = 5)
ggplot(subset(data.all, Parameter %in% detects)) + 
  geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# Recruitment parameters
pdf('Plots/ModelComp_R_year.pdf', width = 8, height = 5)
ggplot(subset(data.all, Parameter %in% R_year)) + 
  geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# Overall density parameters
pdf('Plots/ModelComp_D.pdf', width = 8, height = 5)
ggplot(subset(data.all, Parameter %in% D)) + 
  geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# Site-specific population sizes (summed over age classes)
pdf('Plots/ModelComp_N_exp.pdf', width = 8, height = 5)
for(j in 1:N_sites){
  print(
    ggplot(subset(data.all, Parameter %in% paste0('N_exp[', j, ', ', 1:N_years, ']'))) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      ggtitle(paste0('Site ', j)) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
}
dev.off()

# Site-specific densities (summed over age classes)
pdf('Plots/ModelComp_Density.pdf', width = 8, height = 5)
for(j in 1:N_sites){
  print(
    ggplot(subset(data.all, Parameter %in% paste0('Density[', j, ', ', 1:N_years, ']') &
                            Value < 1)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      ggtitle(paste0('Site ', j)) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
}
dev.off()

# Site- and age-specific population sizes
pdf('Plots/ModelComp_N_a_exp.pdf', width = 8, height = 5)
for(j in 1:N_sites){
  print(
    ggplot(subset(data.all, Parameter %in% c(paste0('N_exp[1, ', j, ', ', 1:N_years, ']'), 
                                              paste0('N_exp[2, ', j, ', ', 1:N_years, ']')) &
                             Value != 0)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_manual(values = viridis(3)[1:2]) + 
      scale_color_manual(values = viridis(3)[1:2]) + 
      ggtitle(paste0('Site ', j)) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
}
dev.off()


# Site- and age-specific densities
pdf('Plots/ModelComp_Density_a.pdf', width = 8, height = 5)
for(j in 1:N_sites){
  print(
    ggplot(subset(data.dens, Parameter %in% c(paste0('Density[1, ', j, ', ', 1:N_years, ']'), 
                                              paste0('Density[2, ', j, ', ', 1:N_years, ']')) &
                    Value != 0)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      ggtitle(paste0('Site ', j)) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
}
dev.off()