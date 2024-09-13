#Script with functions to modify function implemented in "DropSeq.util" package for reading in data of DropViz 
#Original function did not work well
if (!require("Matrix")){
  install.packages("Matrix")
  library("Matrix")
}

strStartsWith = function (theString, thePrefix) 
{
  return(strtrim(theString, nchar(thePrefix)) == thePrefix)
}


loadSparseDgeNames = function (file) 
{
  conn = file(file, "r")
  genes = c()
  cell_barcodes = c()
  while (TRUE) {
    line = readLines(con = conn, n = 1, ok = FALSE)
    if (!strStartsWith(line, "%")) {
      break
    }
    if (grepl(pattern="%%GENES\t", x=line)) {
      these_genes = strsplit(line, "\t", fixed = TRUE)[[1]]
      genes = append(genes, these_genes[2:length(these_genes)])
    }
    else if (grepl(pattern="%%CELL_BARCODES\t", x=line)) {
      these_cells = strsplit(line, "\t", fixed = TRUE)[[1]]
      cell_barcodes = append(cell_barcodes, these_cells[2:length(these_cells)])
    }
  }
  close(conn)
  return(list(genes = genes, cell_barcodes = cell_barcodes))
}

loadSparseDge = function(file, column_compressed = FALSE){
  genes_and_cell_barcodes = loadSparseDgeNames(file)
  ret = readMM(file)
  rownames(ret) = genes_and_cell_barcodes$genes
  colnames(ret) = genes_and_cell_barcodes$cell_barcodes
  if (column_compressed) {
    return(convertTripletCompressedMatrixToColumnCompressed(ret))
  }
  else {
    return(ret)
  }
}