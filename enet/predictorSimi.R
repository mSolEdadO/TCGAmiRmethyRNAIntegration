library(VennDiagram)
library(data.table)
interacs=fread("slctdPrdis.tsv")

cpgs=interacs[substr(interacs$predictor,1,1)=="c",]
mirnas=interacs[substr(interacs$predictor,1,1)=="h",]
genes=interacs[substr(interacs$predictor,1,1)=="E",]
cpgs=lapply(unique(cpgs$subtype),function(x) 
    paste(cpgs$pam50[cpgs$subtype==x],cpgs$predictor[cpgs$subtype==x],sep='.'))
mirnas=lapply(unique(mirnas$subtype),function(x) 
        paste(mirnas$pam50[mirnas$subtype==x],mirnas$predictor[mirnas$subtype==x],sep='.'))
genes=lapply(unique(genes$subtype),function(x) 
      paste(genes$pam50[genes$subtype==x],genes$predictor[genes$subtype==x],sep='.'))

gg_color_hue <- function(n) {
   hues = seq(15, 375, length = n + 1)
   hcl(h = hues, l = 65, c = 100)[1:n]
}

venn.diagram(x=list(A=cpgs[[1]],B=cpgs[[2]],C=cpgs[[3]],D=cpgs[[4]],E=cpgs[[5]]),
            filename="cpgs.tiff",col = "transparent",alpha = 0.50,cex = 1.5,
            label.col="black",cat.cex = 1.5,margin = 0.1,
            category.names=unique(interacs$subtype),fill=gg_color_hue(5),
            cat.fontfamily="sans",fontfamily="sans")
venn.diagram(x=list(A=mirnas[[1]],B=mirnas[[2]],C=mirnas[[3]],D=mirnas[[4]],E=mirnas[[5]]),
            filename="mirnas.tiff",col = "transparent",alpha = 0.50,cex = 1.5, 
            label.col="black",cat.cex = 1.5,margin = 0.1,
            category.names=unique(interacs$subtype),fill=gg_color_hue(5),
            cat.fontfamily="sans",fontfamily="sans")
venn.diagram(x=list(A=genes[[1]],B=genes[[2]],C=genes[[3]],D=genes[[4]],E=genes[[5]]),
            filename="genes.tiff",col = "transparent",alpha = 0.50,cex = 1.5, 
            label.col="black",cat.cex = 1.5,margin = 0.1,
            category.names=unique(interacs$subtype),fill=gg_color_hue(5),
            cat.fontfamily="sans",fontfamily="sans")
