#Brett Mclain
#

#R script to be paired with MPC.py files and its varients to automatically
#create simple charts to visualize its outputted data

#to function, previously run MPC.py must be run with 

args <- commandArgs(trailingOnly = TRUE)
#args <- c("") #testing line
library(ggridges)
library(ggplot2)
library(gridExtra)
library(hash)
library(magrittr)


loc <- paste("/project/meisel/users/bmmclain/Y1W1Y2W2/OutputArchive/",
              "data_",
              args[1],
              ".txt",
              sep="")
#uncommented for testing
#loc <- "C:/Users/brett/OneDrive - University Of Houston/Work Study Tests/data_775765.txt"
df <- read.table(loc, header = TRUE)

#Find which data types have been returned
types <- c(8, 54, 6, 54, 8, 2, 34, 20)
df_types <- colnames(df)
data_present <- c()
i <- 1
j <- 1
if(df_types[i]=="YM_fem_fitness"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="f1_fit"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="YM_fem_dom"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="f1"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="X_freq"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="female"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="f1_equil"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}
j <- j + 1
if(df_types[i]=="m1_equil"){
  data_present <- append(data_present, i)
  i <- i + types[j]
  data_present <- append(data_present, i-1)
}

h <- hash()
h[["YM_fem_fitness"]] <- "Allele Fitness"
h[["f1_fit"]] <- "Genotype Fitness"
h[["YM_fem_dom"]] <- "Dom Results"
h[["f1"]] <- "Final Simulation Frequencies"
h[["X_freq"]] <- "Allele Frequencies"
h[["female"]] <- "Sex Distribution"
h[["f1_equil"]] <- "Female Genotype Equilibrium"
h[["m1_equil"]] <- "Male Genotype Equilibrium"

plots <- list()
pdf(paste("OutputArchive/data_", args[1], "_graphs.pdf", sep=""),
    width=15,
    height=7.5)
for( i in 1:(length(data_present)/2)){
  group <- stack(df[,data_present[(2*i-1)]:data_present[(2*i)]])
  temp <- group %>%
    ggplot(aes(x=values, y=ind, fill=ind)) +
    geom_density_ridges(alpha=.6) +
    theme_ridges() +
    theme(legend.position = "none",
          plot.background = element_rect(fill="white",
                                         color="black",
                                         linewidth=1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 5),
          plot.title=element_text(size=10,
                                  hjust = 0.5)) +
    labs(title=h[[df_types[data_present[(2*i-1)]]]],
         x=NULL,
         y=NULL)
  plots[[i]] = temp
}
do.call("grid.arrange", c(plots, ncol=4))
dev.off()


# ggsave(paste("OutputArchive/data_", args[1], "_graphs.jpeg", sep=""),
#               type="cairo",
#               width=2000,
#               height=1000,
#               units="px")

