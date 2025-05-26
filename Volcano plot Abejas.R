
##Cargarmos la libraria de readr
#Para cargar nuestra base de datos que es tsv
library("readr")

#Cargamos nuestra base de datos
data<-read_tsv("/Users/user/OneDrive/Escritorio/Genómica 2025/GSE28235.top.table.tsv")
View(data)#Si salio

#Hacer nuestro volcano plot
with(data, plot(logFC, -log10(P.Value),#Usamos las columnas directamente
                pch = 20, main = "Volcano Plot",
                xlab = "log2 Fold Change", ylab = "-log10(p-value)"))

# Lineas horizontales - y verticales |
abline(h = -log10(0.05), col = "red", lty = 2)   # p < 0.05
#h es para valor en y
#tipo de línea

abline(v = c(-1, 1), col = "blue", lty = 2)      # logFC > 1 o < -1
#V es para valor en x
