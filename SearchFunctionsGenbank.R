#por: José Juliano Amorim Silva (jose.amorim@aluno.ufabc.edu.br)

# FUNÇÕES DE BUSCA PARA O GENBANK:

#1ª: Função para buscar as informações sobre as espécies e baixar em um arquivo xslx:
search_gb <-function(genus, list_of_species){
  library(rentrez)
  library(qpcR)
  library(xlsx)
  
  final_data_search <- data.frame (species = character(),
                                   matK = character(),
                                   atpB_rbcL = character(),
                                   ycf1 = character(),
                                   ITS = character())
  
  for (species in list_of_species){
    
    matK_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND matK[ALL]"), retmax = "25", idtype = "acc")
    if (is.character(matK_search$ids)){
      matK <- matK_search$ids
    } else {
      matK_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND maturase K[ALL]"), retmax = "25", idtype = "acc")
      if (is.character(matK_search$ids)){
        matK <- matK_search$ids
      } else { 
        matK <- "NA"
      }
    }
    
    atpB_rbcL_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND atpB-rbcL[ALL]"),  retmax = "25", idtype = "acc")
    if (is.character(atpB_rbcL_search$ids)){
      atpB_rbcL <- atpB_rbcL_search$ids
    } else {
      atpB_rbcL_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND atpB[GENE]"),  retmax = "25", idtype = "acc")
      if (is.character(atpB_rbcL_search$ids)){
        atpB_rbcL <- atpB_rbcL_search$ids
      } else {
        atpB_rbcL_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND atpB[ALL]"),  retmax = "25", idtype = "acc")
        if (is.character(atpB_rbcL_search$ids)){
          atpB_rbcL <- atpB_rbcL_search$ids
        } else {
          atpB_rbcL <- "NA"
        }
      }
    }
    
    ycf1_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND ycf1[GENE]"),  retmax = "25", idtype = "acc")
    if (is.character(ycf1_search$ids)){
      ycf1 <- ycf1_search$ids
    } else {
      ycf1_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND hypothetical protein RF1 (ycf1)[ALL]"),  retmax = "25", idtype = "acc")
      if (is.character(ycf1_search$ids)){
        ycf1 <- ycf1_search$ids
      } else {
        ycf1_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND hypothetical protein Ycf1[ALL]"),  retmax = "25", idtype = "acc")
        if (is.character(ycf1_search$ids)){
          ycf1 <- ycf1_search$ids
        } else {
          ycf1 <- "NA"
        }
      }
    }
    
    ITS_search <- entrez_search(db = "nucleotide", term = paste(species, "[ORGN] AND internal transcribed spacer[ALL]"), retmax = "25", idtype = "acc")
    if (is.character(ITS_search$ids)){
      ITS <- ITS_search$ids
    } else {
      ITS <- "NA"
    }
    
    data_search <- qpcR:::cbind.na(species,matK,atpB_rbcL,ycf1,ITS)
    
    final_data_search <- rbind.data.frame(final_data_search, data_search)
  }
  
  
  write.xlsx(final_data_search, file= paste(genus, ".xlsx"))
  
  return(final_data_search)
}

#2ª: Funções para baixar as sequências em FASTA (uma função diferente por tipo de sequência):
fasta_atpB <-function(final_data_search){
  library(rentrez)
  library(xlsx)
  
  atpB_fasta <- entrez_fetch(db="nuccore", id = final_data_search$atpB_rbcL, rettype="fasta")
  write.table(atpB_fasta, paste(genus, "_seqs_atpB.txt"))
}

fasta_ITS <-function(final_data_search){
  library(rentrez)
  library(xlsx)
  
  ITS_fasta <- entrez_fetch(db="nuccore", id = final_data_search$ITS, rettype="fasta")
  write.table(ITS_fasta, paste(genus, "_seqs_ITS.txt"))
}

fasta_ycf1 <-function(final_data_search){
  library(rentrez)
  library(xlsx)
  
  ycf1_fasta <- entrez_fetch(db="nuccore", id = final_data_search$ycf1, rettype="fasta")
  write.table(ycf1_fasta, paste(genus, "_seqs_ycf1.txt"))
}

fasta_matK <-function(final_data_search){
  library(rentrez)
  library(xlsx)
  
  matK_fasta <- entrez_fetch(db="nuccore", id = final_data_search$matK, rettype="fasta")
  write.table(matK_fasta, paste(genus, "_seqs_matK.txt"))
}

# COMO RODAR AS FUNÇÕES (DO JEITO MAIS OTIMIZADO POSSIVEL): Armazenar o genero na variavel "genus" as especies na variavel "list_of_species". Entao, é só selecionar as linhas de codigo seguintes e roda-las
#genus <- "cyrtopodium"
#list_of_species <- c("Cyrtopodium punctatum", "Cyrtopodium andersonii", "Cyrtopodium hatschbachii")

#data_frame_fasta <- search_gb(genus, list_of_species)
#fasta_atpB(data_frame_fasta)
#fasta_ITS(data_frame_fasta)
#fasta_ycf1(data_frame_fasta)
#fasta_matK(data_frame_fasta)

# A cada vez que for fazer uma nova pesquisa com outros generos e outras especies é só trocar o genero em "genus" e as especies em "list_of_species" e rodar as linhas de codigo seguintes
# obs.: é possível executar essas funções sem ter que ficar salvando o genero em "genus" e as espécies em "list_of_species" e, ao invés disso, chamar a função search_gb direto dentro da funções fasta, mas ai vai dar um erro na hora de baixar os arquivos de FASTA já com o nome do gênero automaticamente.