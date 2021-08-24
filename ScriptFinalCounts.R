#Importar o arquivo com as contagens do CCDB(CCDB_counts) e com as contagens em branco de Cymbidieae(Cymbidieae_counts) usando Import Dataset-> From text(base)-> Encoding: automatic-> Heading: No-> Row names: Use numbers-> Separator: Whitespace-> Decimal: Comma-> Quote: Double quote-> Commet: None-> na.strings: NA-> String as factors: FALSE(importante!!!)

#Passar as contagens do CCDB para o arquivo em branco:
partial_counts <- Cymbidieae

for(i in 1:458){
  if(i%%2 != 0){
    for(j in 1:964){
      if(j%%2 != 0){
        if(CCDB[i,1]==Cymbidieae[j,1]){
          partial_counts[j+1,1]<-CCDB[i+1,1]
        }
      }
    }
    
  }
}

#Importar o nosso banco de dados de numero cromossomico(Database) usando Import Dataset-> From text(base)-> Encoding: automatic-> Heading: Yes-> Row names: Use numbers-> Separator: Semicolon-> Decimal: Comma-> Quote: Double quote-> Commet: None-> na.strings: NA-> String as factors: FALSE(importante!!!)

#Passar os numeros do banco de dados para o arquivo com as especies:
partial_counts_2 <- partial_counts

for(i in 1:964){
  if(i%%2 == 0 && partial_counts[i,1] == "X"){
    for(j in 1:732){
      if(partial_counts[i-1,1]==Database[j,1]){
        if(!is.na(Database[j,3])){
          partial_counts_2[i,1]<-Database[j,3]
        }
        else if(!is.na(Database[j,2])){
          partial_counts_2[i,1]<- as.integer(Database[j,2])/2
        }
      }
    }
  }
}

library(rio)
export(partial_counts_2, file = "Cymbidieae_counts_parcial_2.txt")

#As espécies que tem mais de um numero cromossomico registrado tem que ser feitas manualmente mesmo.

