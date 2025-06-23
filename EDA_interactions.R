#%%
rm(list = ls())
library(reshape2)
library(dplyr)
library(stringr)
library(tidyr)
library(car)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(latex2exp)
library(ggridges)
library(tibble)
library(plot.matrix)
# 1) choose dataset name
mydataset <- "CARACOLES"
mydir <- paste0("/data/", mydataset, "/")
setwd(paste0(getwd(),mydir))
#%% 2) load matrices
B <- read.csv("matrix_B.csv", header = TRUE, row.names = 1)
B_sd <- read.csv("matrix_B_sd.csv", header = TRUE, row.names = 1)
A <- read.csv("matrix_A.csv", header = TRUE, row.names = 1)
A_sd <- read.csv("matrix_A_sd.csv", header = TRUE, row.names = 1)
env <- as.vector(read.csv("env.csv", header = TRUE, row.names = 1))
env_name <- names(env)
env <- unlist(env)

names_species <- names(B)
n_sp <- length(names_species)

df <- read.csv(paste0(mydataset, "_clean.csv"))
year <- df$year
nYears <- length(year)

#%% ---- 3) create matrices for each year
matrix.total <- list()
for(i in 1:nYears){
  matrix.total[[i]] <- A + B * env[i]
}
#%% Compare values of matrix.total with histograms each year

# create_hist <- function(df, year) {
#    ggplot(df, aes(x = value, after_stat(ndensity))) + 
#    geom_histogram(bins = 20, alpha = 0.8, position="identity") +
#    labs(title = year) +
#   theme_classic()
# }

# plots_list <- list()

# for (i in 1:nYears){
#   df <- data.frame(value = unlist(matrix.total[[i]]))
#   plots_list[[i]] <- create_hist(df, year[i])
# }

# grid.arrange(grobs = plots_list, ncol = 3, legend = "bottom")

interactions.df <- data.frame()
for(current_year in 1:length(year)){ # 
    # interactions different for each year:    
    interactions.df <- rbind(interactions.df,
                            list(year = rep(as.character(year[current_year]),n_sp*n_sp), 
                            values = unlist(matrix.total[[current_year]]),
                            environment = rep(as.numeric(env[current_year]),n_sp*n_sp)))
    }

myplot<- ggplot(interactions.df, aes(x = values, y = year, height = after_stat(density))) +
  geom_density_ridges(rel_min_height = 0.002, scale = 3,alpha = 0.75) +
  # scale_y_discrete(expand = c(0.01, 0)) +
  # scale_x_continuous(expand = c(0.01, 0)) +
  #scale_fill_brewer(palette = 4) +
  #scale_fill_viridis(discrete = TRUE, option = "D") + 
  labs(title = 'Distribution interactions') +
  theme_ridges() + theme(legend.position = "none", text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) 
print(myplot)

#%% ---- 5) tiles of total matrix, and change of sign and variance of interactions
years_list <- list()
inte_prev <- matrix.total[[1]]
inte_same <- 0
varz <- c()
skw <- c()
lims <- c(abs(min(unlist(matrix.total))),abs(max(unlist(matrix.total))))
idx <- which.max(lims)
for (current_year in 1:nYears){
  interactions <- matrix.total[[current_year]]
  inte_same <- as.integer(sign(inte_prev) == sign(interactions)) + inte_same
  inte_prev <- interactions
  varz[current_year] <- var(unlist(interactions))
  skw[current_year] <- (mean(unlist(interactions)) - median(unlist(interactions)))/sd(unlist(interactions))
  
  df <- interactions %>%
    as.data.frame() %>%
    rownames_to_column("Fila") %>%
    pivot_longer(-Fila, names_to = "Columna", values_to = "Valor")
  
  p <- ggplot(df, aes(Columna, Fila, fill = Valor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
    limits   = c(-lims[idx],  lims[idx]))    +      
    coord_fixed() +
    labs(title = paste("Año", current_year)) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title  = element_blank(),
      plot.title  = element_text(hjust = 0.5)
    )
  
  print(p)
  years_list[[current_year]] <- p
}

#%% 6) cummulative change in signs of interactions
inte_changes <- matrix(nYears,n_sp,n_sp) - matrix(inte_same,n_sp,n_sp)

# 7) variance of interactions with precipitation

df <- data.frame(year = as.numeric(year), sigma = sqrt(as.numeric(varz)), precipitation = env)
myplot <- ggplot(df, aes(x = precipitation, y = sigma, color = year)) +
  geom_point(size = 3) +
  geom_line(aes(group = year), linewidth = 1) +
  labs(x = 'env factor', y = 'Standard deviation of interactions')+
  theme_bw() + theme(legend.position = "right", text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
ggsave("variance_ABP.pdf")      
print(myplot)

# 8) skewness of interactions with precipitation
names(skw) <- year
df <- data.frame(year = as.numeric(names(skw)), skewness = as.numeric(skw), precipitation = env)
myplot <- ggplot(df, aes(x = precipitation, y = skewness, color = year)) +
  geom_point(size = 3) + 
  geom_line(aes(group = year), linewidth = 1) +
  labs(x = 'env factor', y = 'Skewness of interactions')+
  theme_bw() + theme(legend.position = "right", text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
ggsave("skewness_ABP.pdf")      
print(myplot)


# 9) plot of the interactions with precipitation
# delete diagonal of each matrix.total
matrix.total.nodiag <- lapply(matrix.total, function(mat) {
  mat_no_diag <- mat
  diag(mat_no_diag) <- NA
  mat_no_diag
})

# Prepare inter.df for plotting (excluding diagonal)
inter.df <- data.frame()
for(current_year in 1:length(year)){
  mat <- matrix.total.nodiag[[current_year]]
  values <- as.vector(mat)
  # Remove NA (diagonal)
  values <- unlist(values)
  names(values) <- NULL
  values <- values[!is.na(values)]
  inter.df <- rbind(inter.df,
                    data.frame(
                      year = rep(as.character(year[current_year]), n_sp * n_sp - n_sp),
                      values = values,
                      precipitation = rep(as.numeric(env[current_year]), n_sp * n_sp - n_sp)
                    )
  )
}

myplot <- ggplot(inter.df, aes(x = year, y = values, color = year)) +
  geom_point(size = 3) + 
 #¢ scale_color_viridis(discrete = TRUE, option = "D")+ 
  labs(x = 'Total precipitation', y = TeX("Interspecific interactions"))+
  theme_bw() + theme(legend.position = "right", text = element_text(size=20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
print(myplot)


# %% 15) ratio mut and comp links per species


# 16) D, difference btw intra and inter specific interactions
