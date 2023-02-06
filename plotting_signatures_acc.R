library(ggplot2)
library(dplyr)
library(scales)
setwd("/u/sauves/gene_signatures")
data = read.csv("RES/SIGNATURES/signatures_2022-11-09T17:44:32.745/tst_accs.csv")
title = "Logistic Regression dans Flux sur GTEX : détermination des 30 types de tissus. 
Échantillons: 17382. Input features tailles variables (Gene Expression RNA-Seq log10(tpm+1)). 
Pas de régularization. 5-Fold Cross-val (test n=3400 (20%)). crossentropy loss."
ggplot(data, aes(x = lengths, y = tst_acc * 100)) + geom_point() +
  theme_light() + 
  scale_x_continuous(trans = log10_trans(), breaks = unique(data$lengths)) +
                      # breaks = trans_breaks("log10", function(x) 10^x),
                      # labels = trans_format("log10", math_format(10^.x))) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^ x), 
  #               labels = trans_format("log10", math_format(10^.x))) + 
  #annotation_logticks() + 
  xlab("Number of genes (randomly picked) - log10 scale") + 
  ylab("Accuracy % on test set") + 
  ggtitle(title) + 
  theme(text = element_text(size = 14, color = "black"),
        title = element_text(size = 10, color = "black") ) 
ggsave("./RES/SIGNATURES/signatures_2022-11-09T17:35:59.818/tst_accs.png")
