munge_sce_mat = function(data_obj,  mapping_df, assay_name = "all",multi_mapping = "sum") {
  #check if the assay exists
  if ( assay_name != "all" & !assay_name %in% SummarizedExperiment::assayNames(data_obj)) {
    stop("The assay you are indicating does not exist")
  }
  #check if the feature name is correct
  if( is.null(rownames(data_obj)) | !any( rownames(data_obj) %in% (mapping_df %>% dplyr::pull(1)))) {
    stop("The feature names do not match the first column of the mapping_df")
  }
  if(any(duplicated(rownames(data_obj)))) {
    stop("The feature names isn't unique. Make it unique or change the feature name into unique ones and then use the function to map non-unique features to unique features.")
  }
  if (assay_name=="all") {
    assay_name =  SummarizedExperiment::assayNames(data_obj)
  }
  
  #merge
  new_assay = list()
  mapping_df = mapping_df %>% dplyr::mutate_all(~as.character(.))  %>% dplyr::distinct_at(c(1,2))
  
  for (assay_name_i in assay_name) {
    data_mat = SummarizedExperiment::assay(data_obj, assay_name_i)
    #filter mapping
    all_mapping = mapping_df %>% 
      dplyr::filter(.[[1]] %in% rownames(data_mat)) %>% #filter by feature names
      tidyr::drop_na(2) #drop features without mapping
    #subset rows
    data_mat = data_mat[match(all_mapping[[1]],rownames(data_mat)),]
    #split matrix into two
    all_mapping = all_mapping %>% 
      dplyr::group_by_at(2) %>%
      dplyr::add_count() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(is_multi_mapping = ifelse(n>=2, TRUE, FALSE))
    
    multiple_mapping = all_mapping %>% 
      dplyr::filter(is_multi_mapping)
    single_mapping = all_mapping %>%
      dplyr::filter(!is_multi_mapping)
    #separate the matrix to the one with multiple mapping and the other
    single_mat = data_mat[which(!all_mapping[["is_multi_mapping"]]),]
    multi_mat = data_mat[which(all_mapping[["is_multi_mapping"]]),]
    
    #transform mapping 
    fac_mat = Matrix::fac2sparse(factor(multiple_mapping[[2]], levels = unique(multiple_mapping[[2]]))) 
    if(multi_mapping=="mean") {
      fac_mat = sweep_sparse(fac_mat, margin=1, stats = Matrix::rowSums(fac_mat), fun = "/")
    }
    #merge and set rownames
    multi_mat = ( fac_mat %*% multi_mat) %>% 
      magrittr::set_rownames(unique(multiple_mapping[[2]])) %>%
      magrittr::set_colnames(colnames(data_mat))
    single_mat = single_mat %>% magrittr::set_rownames(single_mapping[[2]])
    data_mat = rbind(single_mat, multi_mat)
    #rearrange based on the new gene names
    data_mat = data_mat[match(unique(all_mapping[[2]]),rownames(data_mat)),]
    #add assay
    new_assay = new_assay %>% append(data_mat)
  }
  new_assay = new_assay %>% purrr::set_names(assay_name)
  
  #create output data object
  out_data = SingleCellExperiment::SingleCellExperiment(
    assays = new_assay,
    colData = SummarizedExperiment::colData(data_obj),
    rowData = data.frame(rownames(new_assay[[1]]), row.names = rownames(new_assay[[1]]),fix.empty.names = FALSE) %>% 
      magrittr::set_colnames(colnames(mapping_df)[2]) %>%
      S4Vectors::DataFrame()
  )
  
  return(out_data)
}