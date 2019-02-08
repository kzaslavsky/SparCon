##interfacing with GO to examine expression datasets
library(GO.db)


##function to select correct GOTERM for querying
##input: query string
##output: string with GOterm name
selectGO <- function(query)
{
  #query GO
  keyGO.i <- grep(query, keys(GO.db, keytype="TERM"))
  keyGO.char <- keys(GO.db, keytype="TERM")[keyGO.i]
  
  #show user output and ask to pick
  if(length(keyGO.i)>1)
    {
  print(keyGO.char)
  GOindex <- as.integer(readline(prompt="Pick the number of the GO term you are interested in:"))
  return(keyGO.char[GOindex])}
  else
  {print("Term not found; use another query")}
}

##function to generate a query list of genes for 
##grep from a GO category
##input: GOterm name
##output: geneQuery string to be used for grep

geneQuery.grep <- function(query, querydb = org.Hs.eg.db)
{
  
  
  keyGO.char <- selectGO(query)
  #get all the gene names associated with GOterm
  GOterm.id <- select(GO.db, keyGO.char, keytype="TERM", columns="GOID" )[,2]
  GOterm.idoffspring <- get(GOterm.id, GOBPOFFSPRING)
  #GOterm.idoffspring <- get(GOterm.id, GOMFOFFSPRING)    #use for MF
  GOterm.idchildren <- get(GOterm.id, GOBPCHILDREN)
  GOterm.ids <- c(GOterm.id, GOterm.idoffspring)
  #GOterm.ids <- as.character(GOterm.ids[2,])
  GOterm.genes <- na.omit(select(querydb, GOterm.id, keytype="GOALL", columns="SYMBOL"))
  
  geneQuery.args <- c(as.list(GOterm.genes[,4]),sep="$|^")
  geneQuery <- do.call(paste, geneQuery.args)
}