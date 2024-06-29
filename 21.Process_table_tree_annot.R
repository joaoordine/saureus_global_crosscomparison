library(dplyr)
library(readr)

# Read the TSV file
datum <- read_tsv("sccmec_staphopia.tsv")

View(datum)


# Define the columns to check for TRUE values
columns_to_check <- c("I", "Ia")

# Add a new column 'I_all' based on the condition
datum <- mutate(datum, I_all = ifelse(rowSums(select(datum,all_of(columns_to_check)) == "TRUE") > 0, "TRUE", "FALSE"))
View(datum)


# Define the columns to check for TRUE values
columns_to_check <- c("II", "IIa", "IIb")

# Add a new column 'II_all' based on the condition
datum <- mutate(datum, II_all = ifelse(rowSums(select(datum,all_of(columns_to_check)) == "TRUE") > 0, "TRUE", "FALSE"))
View(datum)


# Define the columns to check for TRUE values
columns_to_check <- c("III", "IIIa")

# Add a new column 'III_all' based on the condition
datum <- mutate(datum, III_all = ifelse(rowSums(select(datum,all_of(columns_to_check)) == "TRUE") > 0, "TRUE", "FALSE"))
View(datum)

# Define the columns to check for TRUE values
columns_to_check <- c("IV","IVa", "IVb", "IVc", "IVd", "IVg", "IVh")

# Add a new column 'IV_all' based on the condition
datum <- mutate(datum, IV_all = ifelse(rowSums(select(datum,all_of(columns_to_check)) == "TRUE") > 0, "TRUE", "FALSE"))

View(datum)

# Define the columns to check for TRUE values
columns_to_check <- c("I_all", "II_all", "III_all", "IV_all", "V", "VI", "VII", "VIII", "IX", "meca")

# Add a new column 'III_all' based on the condition
datum <- mutate(datum, NEGATIVE = ifelse(rowSums(select(datum,all_of(columns_to_check)) == "TRUE") > 0, "TRUE", "FALSE"))
View(datum)



# Define the columns to select
columns_to_select <- c("sample", "I_all", "II_all", "III_all", "IV_all", "V", "VI", "VII", "VIII", "IX", "meca", "NEGATIVE")

# Select the specified columns
selected_data <- select(datum, all_of(columns_to_select))

# Write the modified data to a new TSV file
write_tsv(datum, "output_ssmec.tsv")


# Create another collumm for color
next_one <- read_tsv("output_ssmec.tsv")

# Definir as colunas de cores e suas correspondências
cores <- c("I_all" = "#FFEF9F", "II_all" = "#7BF1A8", "III_all" = "#90F1EF", "IV_all" = "#4A7B9D", "V"= "#7C616C", "VI" = "#7E007B",
           "VII"= "#B33951", "VIII"= "#91C7B1", "IX"= "#5158BB")

# Função para atribuir cor com base no valor de uma coluna
atribuir_cor <- function(valor) {
  for (coluna in names(cores)) {
    if (valor[coluna] == "TRUE") {
      return(cores[coluna])
    }
  }
  if (valor["meca"] == "TRUE") {
    return("#000300")
  } else {
    return("#FFFEFF")  # Se nenhuma correspondência for encontrada
  }
}



# Função para atribuir formato com base no valor de uma coluna
atribuir_formato <- function(valor) {
  if (valor["meca"] == "TRUE") {
    return("h")
  } else {
    for (coluna in names(cores)) {
      if (valor[coluna] == "TRUE") {
        return("s")
      }
    }
    return("s")  # Substituir NULL por "s"
  }
}
# Adicionar a nova coluna "shape" com base nas condições especificadas
last_table <- mutate(next_one, shape = sapply(1:nrow(next_one), function(i) atribuir_formato(next_one[i, ])))

# Criar as novas colunas "COLOR" e "formato" com base nas condições especificadas
last_table <- mutate(last_table, 
               COLOR = sapply(1:nrow(last_table), function(i) atribuir_cor(last_table[i, ])))
                     
                     


last_table <- mutate(last_table, level = 1)
last_table <- mutate(last_table, colorir = "ring_color")
last_table <- mutate(last_table, shapado = "ring_shape")
last_table <- mutate(last_table, edge = "ring_edge_color")
last_table <- mutate(last_table, edge_color = "black")
View(last_table)

# Write the modified data to a new TSV file
write_tsv(last_table, "output_final_ssmec.tsv")


#Processamento do arquivo mlst que têm os IDs para os sequences types
#Este arquivo vai ser usado para colorir os ramos da árvore
#Processamento de dividir baseado na coluna ST foi baseado em um script 
st <- read_tsv("/home/itachi/Desktop/Projeto_Paralelo_LZ/Joao_ordine/ST_result.tsv")

#Definir cor por ST e colocar coluna TYPE
st <- mutate(st, TYPE = "branch")
View(st)

# Definir as cores para cada valor de ST
cores <- c("8" = "#E71D36", "5" = "#1B2CC1", "105" = "#2BC016", "Unclassified" = "#D3D3D3", 
           "30" = "#FAF33E", "398" = "#ff8000", "45" = "#62259d", "1" = "#ff89e1", 
           "59" = "#ffd2a3", "72" = "#9b0058", "239" = "#00afad", "1292" = "#006fe6", "22" = "#7d6c06",
           "15"= "#00FFFF", "188" = "#8000ff", "9" = "#EF27A6", "6" = "#0b5394", "25" = "#c8a2c8", "228" = "#d48a3f", "7" = "#a3ff5a", 
           "97" = "#800080", "225" = "#00FFC5", "121" = "#ff7b7b") 

# Criar a coluna "CORES" com as cores correspondentes aos valores de ST

st$CORES <- ifelse(as.character(st$ST) %in% names(cores), cores[as.character(st$ST)], "#000000")

View(st)

#Remover linhas indesejadas
#st <- st[st$CORES != "NULL", ]
View(st)

# Adicionar "_genomic" ao final de cada valor na coluna "Genomes"
st$Genomes <- paste0(st$Genomes, "_genomic")
st <- mutate(st, STYLE = "normal")
st <- mutate(st, LABEL = "label_background")
st
write_tsv(st, "color_branches_new_colors.tsv")


#Conditional stational for Sccmec
#Upload arquivo novo
new_staphopia <- read_tsv("Saureus_Tables_UPDATED_sccmec_sum_staphopia.tsv")
View(new_staphopia)


# Definir as colunas de cores e suas correspondências
columns_check_staphopia_top <- c("I" = "I", "II" = "II", "III" = "III", "IV" = "IV", "V"= "V", "VI" = "VI",
           "VII"= "VII", "VIII"= "VIII", "IX"= "IX")


#segunda tentativa
# Função para atribuir cor com base no valor de uma coluna
atribuir_col <- function(nome) {
  if (nome["Unassigned"] == "TRUE") {
    return("UNASSIGNED")
  } else if (nome["SCCmec_neg"] == "TRUE") {
    return("SCCmec_NEG")
  } else {
    for (coluna in names(columns_check_staphopia_top)) {
      if (nome[coluna] == "TRUE") {
        return(columns_check_staphopia_top[coluna])
      }
    }
    return("FALSE")  # Se não houver correspondência
  }
}




new_staphopia$Sccmec_category <- sapply(1:nrow(new_staphopia), function(i) atribuir_col(new_staphopia[i,]))

View(new_staphopia)
write_tsv(new_staphopia, "new_staphopia.tsv")

#Anotação com novo documento
new_staphopia_annot <- read_tsv("/home/itachi/Desktop/Projeto_Paralelo_LZ/Joao_ordine/new_staphopia.tsv")

cor_names <- c("I" = "#8bd2c6", "II" = "#ffffb1", "III" = "#bdb9d9", "IV" = "#fb7f71", "V"= "#fccce4", "VI" = "#d8d8d8",
               "VII"= "#bb7fbc", "VIII"= "#cbeac4", "IX"= "#7fb0d2", "UNASSIGNED" = "#b2dd69", "SCCmec_NEG" = "#fdb262")

new_staphopia_annot$COR_NAMES <- ifelse(as.character(new_staphopia_annot$Sccmec_category) %in% names(cor_names),
                                        cor_names[as.character(new_staphopia_annot$Sccmec_category)], "NULL")
View(new_staphopia_annot)

new_staphopia_annot$COR_NAMES[new_staphopia_annot$COR_NAMES == 'NULL'] <- '#FFFFFF'
View(new_staphopia_annot)

#If everyone is false, but for mecA is true
#new_staphopia_annot <- new_staphopia_annot %>%
 # mutate(SHAPE = if_else(COR_NAMES == '#1B1B1E', "2", "1"))
View(new_staphopia_annot)

#for all mecA, regards if it's true or not

new_staphopia_annot$star <- ifelse(new_staphopia_annot$mecA == "TRUE", "#1B1B1E", "#FFFFFF")
View(new_staphopia_annot)
write_tsv(new_staphopia_annot, "new_staphopia_1.tsv")

#Option 02 for mecA - add shape instead
new_staphopia_annot$SHAPE_mec <- ifelse(new_staphopia_annot$star == "#1B1B1E", "0", "-1")
write_tsv(new_staphopia_annot, "new_staphopia_2.tsv")

#MRSA versus MSSA

MSSA_table <-read_tsv("Saureus_Tables_UPDATED_phylo_annot.tsv")
data_MSSA<- c("MSSA" = "#C42021", "MRSA" = "#3BC14A")
MSSA_table$beta <- ifelse(as.character(MSSA_table$bLactam_Category) %in% names(data_MSSA),
                                        data_MSSA[as.character(MSSA_table$bLactam_Category)], "NULL")
View(MSSA_table)
MSSA_table$Genome <- paste0(MSSA_table$Genome, "_genomic")
write_tsv(MSSA_table, "MSSA_annot.tsv")

#Fluoroquinolone resistance

fluor_table <-read_tsv("Saureus_Tables_UPDATED_phylo_annot.tsv")
data_fluor<- c("0" = "0", "1" = "1")
fluor_table$fluor <- ifelse(as.character(fluor_table$qac_Fluoroquinolone_Resistance) %in% names(data_fluor),
                          data_fluor[as.character(fluor_table$qac_Fluoroquinolone_Resistance)], "NULL")
View(fluor_table)
fluor_table$Genome <- paste0(fluor_table$Genome, "_genomic")
write_tsv(fluor_table, "fluor_annot.tsv")


##Fluoroquinolone resistance (For another config)

fluor_table <-read_tsv("Saureus_Tables_UPDATED_phylo_annot.tsv")
data_fluor<- c("0" = "#C42021", "1" = "#3BC14A")
fluor_table$fluor <- ifelse(as.character(fluor_table$qac_Fluoroquinolone_Resistance) %in% names(data_fluor),
                            data_fluor[as.character(fluor_table$qac_Fluoroquinolone_Resistance)], "NULL")
View(fluor_table)
fluor_table$Genome <- paste0(fluor_table$Genome, "_genomic")
fluor_table$fluor[fluor_table$fluor == 'NULL'] <- '#3BC14A'
write_tsv(fluor_table, "fluor_annot_2.tsv")



#Virulence factor
#import data
vf_data <- read_tsv("C:/Users/edson/Downloads/Saureus_Tables_UPDATED_abricate_VFs.tsv")
  

#luks
vf_data$lukScolor <- ifelse(vf_data$`lukS-PV` == ".", "#FFFFFF", "#000000")
View(vf_data)
vf_data$FILE <- paste0(vf_data$FILE, "_genomic")

#lukF
vf_data$lukFcolor <- ifelse(vf_data$`lukF-PV` == ".", "#FFFFFF", "#000000")

vf_data[2:85] <- list(NULL)

View(vf_data)
write_tsv(vf_data, "vf_ring.tsv")


#Analysis for carbo transporters
ribose_data <- read_tsv("/home/itachi/Downloads/ribose_transporters_df.tsv")

#For ribose
View(ribose_data)
hist(ribose_data$Transporter_Count)
med_rib <- median(ribose_data$Transporter_Count)
ribose_data$color <- ifelse(ribose_data$Transporter_Count > med_rib, "#2c43b0","#64b6ee")
View(ribose_data)
ribose_data$Genome <- paste0(ribose_data$Genome, "_genomic")
write_tsv(ribose_data, "ribose_ring.tsv")

#For lactose
lactose_data <- read_tsv("/home/itachi/Downloads/lactose_transporters_df.tsv")
hist(lactose_data$Transporter_Count)
med_lac <- median(lactose_data$Transporter_Count)
lactose_data$color <- ifelse(lactose_data$Transporter_Count > med_lac, "#b21500", "#ff4b33")
View(lactose_data)
lactose_data$Genome <- paste0(lactose_data$Genome, "_genomic")
write_tsv(lactose_data, "lactose_ring.tsv")

#For manitol
manitol_data <- read_tsv("/home/itachi/Downloads/mannitol_transporters_df.tsv")
hist(manitol_data$Transporter_Count)
med_mani <- median(manitol_data$Transporter_Count)
manitol_data$color <- ifelse(manitol_data$Transporter_Count > med_mani, "#607c3c", "#b5e550")
View(manitol_data)
manitol_data$Genome <- paste0(manitol_data$Genome, "_genomic")
write_tsv(manitol_data, "manitol_ring.tsv")

#For sia transporters
siatrans_data <- read_tsv("/home/itachi/Downloads/sia_transporters_df.tsv")
hist(siatrans_data$Transporter_Count)
med_siatrans <- median(siatrans_data$Transporter_Count)
siatrans_data$color <- ifelse(siatrans_data$Transporter_Count > med_siatrans, "#a98600", "#e9d700")
View(siatrans_data)
siatrans_data$Genome <- paste0(siatrans_data$Genome, "_genomic")
write_tsv(siatrans_data, "sia_trans_ring.tsv")


