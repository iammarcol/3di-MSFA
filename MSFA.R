library(ggplot2)
library(dplyr)
suppressPackageStartupMessages(library("argparse"))  
parser <- ArgumentParser()
parser$add_argument("--dir", help="directory of the job")
parser$add_argument("--jobid", help="Job Id")
args <- parser$parse_args()

# this added part of the script is used to make 

# function to process the input MSA and create the output data frame
process_msa <- function(file_path) {
  # read the FASTA file
  fasta_data <- readLines(file_path)
  
  # extract sequence IDs and sequences
  ids <- fasta_data[grep("^>", fasta_data)]
  sequences <- fasta_data[!grepl("^>", fasta_data)]
  
  # remove ">" from PDBID
  ids <- gsub(">", "", ids)
  
  df <- data.frame(
    PDBID = rep(ids, sapply(strsplit(sequences, ""), length)),
    `3di` = unlist(strsplit(sequences, ""))
  )
  
  return(df)
}

# process the MSA
df <- process_msa(paste(args$dir, 'MSA_3di_reference.fasta', sep="/"))

# write the data frame to CSV
output_file <- paste(args$dir, paste('AuxFiles/DF_Colores_3di', '.csv', sep=""), sep="/")
write.csv(df, file = output_file, row.names = FALSE, quote = FALSE)

t <- read.table(paste(args$dir,'AuxFiles/DF_Colores', sep=""), sep = ',', header = TRUE)
t2 <- read.table(paste(args$dir, paste('AuxFiles/DF_Colores_3di', '.csv', sep=""), sep="/"), header = TRUE)

merged_df <- cbind(t, t2)
col_names <- colnames(merged_df)
new_col_names <- gsub("\\.", ",", col_names)
colnames(merged_df) <- new_col_names

# here DF_Colored_MODIF is generated - it contains information of 3di added to the original DF_Colores
output_file <- paste(args$dir, paste('AuxFiles/DF_Colores_MODIF', '.csv', sep=""), sep="/")
write.csv(merged_df, file = output_file, row.names = FALSE, quote = FALSE)



######################  ORIGINAL SCRIPT ############################

# Cargar los datos desde el archivo
t <- read.table(paste(args$dir,'AuxFiles/DF_Colores', sep=""), sep = ',', header = TRUE)


# Crear el dataframe
df <- data.frame(
  PDBID = t$PDBID,
  pos = t$pos,
  AA = t$AA,
  FstState= t$FstState,
  FstI=t$FstI
)

# Definir los colores
colores <- c("NEU" = "gray", "MAX" = "red", "MIN" = "green", "-" = "white")

# Filtrar las filas con valores NA en la columna FstI
df_filtered <- df[!is.na(df$FstI), ]

filas <- list()

# Obtener los PDBID únicos
pdbids <- unique(df$PDBID)

# Iterar sobre cada PDBID
for (pdbid in pdbids) {
  # Filtrar los datos por PDBID
  datos_pdbid <- df[df$PDBID == pdbid, ]
  
  # Crear una fila con los valores de FstI
  fila <- datos_pdbid$FstI
  
  # Agregar la fila a la lista
  filas[[pdbid]] <- fila
}

# Crear un nuevo data frame con las filas
df_resultado <- as.data.frame(do.call(rbind, filas))
rownames(df_resultado) <- pdbids

df_resultado[df_resultado == "N/A"] <- NA
df_resultado[is.na(df_resultado)] <- 0

# Realizar el clustering utilizando K-means
num_clusters <- as.integer(length(pdbids)/3)  # Número de clusters deseado
set.seed(123)  # Fijar semilla para reproducibilidad
clusters <- kmeans(df_resultado, centers = num_clusters)

# Obtener los resultados del clustering
df_resultado$cluster <- clusters$cluster

# Ordenar la lista por el clustering
df_ordenado <- df_resultado[order(df_resultado$cluster), ]

# Agregar columna de clustering al dataframe
df_filtered$clustering <- df_resultado$cluster
cluster_labels <- as.character(clusters$cluster)
df_filtered$PDBID <- factor(df_filtered$PDBID, levels = rev(unique(df_filtered$PDBID)))
lfilas = length(filas)
laa=length(t$AA)
if ( lfilas < 20){
  laa=laa*lfilas
  lfilas = 20
}
# Ordenar el gráfico por clustering
p <- ggplot(df_filtered, aes(x = pos, y = reorder(PDBID, PDBID), fill = FstState)) +
  geom_tile(colour = NA, linewidth = 1) +
  geom_text(aes(label = AA), size = lfilas*0.03, color = "black") +
  scale_fill_manual(values = colores) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = lfilas*0.1),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 5),  # Ajusta el tamaño del título de la leyenda aquí
        legend.text = element_text(size = 5),
        legend.margin = margin(t = 1, unit = "lines"),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
warnings()
# Ajustar el tamaño de la leyenda en la parte superior
p <- p + theme(legend.box.margin = margin(t = -lfilas*0.5, unit = "pt"))  # Ajusta el tamaño de la leyenda aquí
h=lfilas*0.08
w=laa*0.0025
ggsave(paste(args$dir,'OutPutFiles/MSFA_',args$jobid,'.png', sep=""), p, width = w, height = h, dpi = 300, bg = "white",limitsize = FALSE)




######################  3di PNG SCRIPT ############################


# Cargar los datos desde el archivo
t <- read.table(paste(args$dir,'AuxFiles/DF_Colores_MODIF.csv', sep=""), sep = ',', header = TRUE)

# Crear el dataframe
df <- data.frame(
  PDBID = t$PDBID,
  pos = t$pos,
  AA = t$X3di,
  FstState= t$FstState,
  FstI=t$FstI
)

# Definir los colores
colores <- c("NEU" = "gray", "MAX" = "red", "MIN" = "green", "-" = "white")

# Filtrar las filas con valores NA en la columna FstI
df_filtered <- df[!is.na(df$FstI), ]

filas <- list()

# Obtener los PDBID únicos
pdbids <- unique(df$PDBID)

# Iterar sobre cada PDBID
for (pdbid in pdbids) {
  # Filtrar los datos por PDBID
  datos_pdbid <- df[df$PDBID == pdbid, ]
  
  # Crear una fila con los valores de FstI
  fila <- datos_pdbid$FstI
  
  # Agregar la fila a la lista
  filas[[pdbid]] <- fila
}

# Crear un nuevo data frame con las filas
df_resultado <- as.data.frame(do.call(rbind, filas))
rownames(df_resultado) <- pdbids

df_resultado[df_resultado == "N/A"] <- NA
df_resultado[is.na(df_resultado)] <- 0

# Realizar el clustering utilizando K-means
num_clusters <- as.integer(length(pdbids)/3)  # Número de clusters deseado
set.seed(123)  # Fijar semilla para reproducibilidad
clusters <- kmeans(df_resultado, centers = num_clusters)

# Obtener los resultados del clustering
df_resultado$cluster <- clusters$cluster

# Ordenar la lista por el clustering
df_ordenado <- df_resultado[order(df_resultado$cluster), ]

# Agregar columna de clustering al dataframe
df_filtered$clustering <- df_resultado$cluster
cluster_labels <- as.character(clusters$cluster)
df_filtered$PDBID <- factor(df_filtered$PDBID, levels = rev(unique(df_filtered$PDBID)))
lfilas = length(filas)
laa=length(t$X3di)
if ( lfilas < 20){
  laa=laa*lfilas
  lfilas = 20
}
# Ordenar el gráfico por clustering
p <- ggplot(df_filtered, aes(x = pos, y = reorder(PDBID, PDBID), fill = FstState)) +
  geom_tile(colour = NA, linewidth = 1) +
  geom_text(aes(label = AA), size = lfilas*0.03, color = "black") +
  scale_fill_manual(values = colores) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = lfilas*0.1),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 5),  # Ajusta el tamaño del título de la leyenda aquí
        legend.text = element_text(size = 5),
        legend.margin = margin(t = 1, unit = "lines"),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
warnings()
# Ajustar el tamaño de la leyenda en la parte superior
p <- p + theme(legend.box.margin = margin(t = -lfilas*0.5, unit = "pt"))  # Ajusta el tamaño de la leyenda aquí
h=lfilas*0.08
w=laa*0.0025
ggsave(paste(args$dir,'OutPutFiles/MSFA_3di_',args$jobid,'.png', sep=""), p, width = w, height = h, dpi = 300, bg = "white",limitsize = FALSE)
