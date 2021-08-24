#por: José Juliano Amorim Silva (jose.amorim@aluno.ufabc.edu.br)

search_genus <- function(x){
  library(chromer)
  library(xlsx)
  y <- chrom_counts(taxa = x, rank = "genus", full = TRUE)
  z <- y[c("resolved_name", "sporophytic", "reference")]
  write.xlsx(z,file="cymbidieae.xlsx")
  return(z)}