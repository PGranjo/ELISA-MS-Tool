################################################################################
## Function to retrieve information from Uniprot for a vector of protein accessions.
## Uses base R for accessions length < 30 and parallelization otherwise.
## Is is used in function: fill.protein.names.in.sa.tables().
################################################################################
retrieve.from.uniprot.parallel <- function(accessions, info, cores){
  
  #INPUT:
  #accessions: a vector of protein accessions.
  #info: either "header", "name", "sequence" or "full" to retrieve the fasta header, protein name, 
  #      protein sequence or complete fasta entry, respectively.
  #cores: the number of cores to be used in parallelization.
  
  #OUTPUT:
  #report_data: the data type corresponding to "info" and for all the protein accessions.
  
  #Check for missing arguments.
  if(missing(x = accessions) | missing(x = info)){
    
    stop(crayon::bgRed(crayon::black("One or more arguments are missing at retrieve.from.uniprot.parallel()...\n")))
    
  }
  
  if(!(info %in% c("header", "name", "sequence", "full"))){
    
    info <- "full"
    
  }
  
  #Create the uniprot fasta urls per protein accession.
  url <- paste0('https://rest.uniprot.org/uniprotkb/', accessions, '.fasta', collapse = NULL)
  
  #If there are more than 30 proteins use base R + parallel.
  if(length(x = accessions) >= 30){
    
    #If the cores argument is missing then assign number of cores detected minus one.
    if(missing(x = cores)){
      
      cores <- parallel::detectCores() - 1
      
    }
    
    #Similarly check is cores argument is higher than the available cores. If so leave one core available for other processes.
    detected_cores <- parallel::detectCores()
    
    if(detected_cores < cores){
      
      cluster <- parallel::makePSOCKcluster(names = detected_cores - 1)
      
    } else {
      
      cluster <- parallel::makePSOCKcluster(names = cores)
      
    }
    
    #Initiate core cluster.
    parallel::setDefaultCluster(cl = cluster)
    
    #Import the Uniprot urls to the cluster.
    parallel::clusterExport(cl = cluster, varlist = c("url"), envir = environment())
    
    #Retrieve fasta for the accessions in parallel.
    retrieved_data <- parallel::parLapply(cl = cluster, url, function(x){ 
      
      return(suppressWarnings(expr = try(expr = read.csv(file = x, header = F, sep = "\n"), silent = TRUE)))
      
    })
    
    #Stop the cores.
    parallel::stopCluster(cl = cluster)
    
  } else { #Else use base R only.
    
    #Retrieve fasta for the accessions with just read.csv function.
    retrieved_data <- lapply(url, function(x){ 
      
      return(suppressWarnings(expr = try(expr = read.csv(file = x, header = F, sep = "\n"), silent = TRUE)))
      
    })
    
  }
  
  #Assign NA to accessions that their fasta was not able to be retrieved.
  retrieved_data[sapply(retrieved_data, function(x){ return(class(x = x))}) == "try-error"] <- as.character(x = NA)
  
  #Prepare information to report based on info argument.
  if(info == "full"){
    
    report_data <- sapply(retrieved_data, function(x){  return(paste0(unlist(x = x), collapse = "\n"))}) 
    
  } else if(info == "sequence"){
    
    report_data <- sapply(retrieved_data, function(x){  
      
      sequence <- unlist(x = x)
      
      if(!is.na(x = sequence[1])){
        
        sequence <- paste(unlist(x = x)[-1], collapse = "", sep = "")
        
      }
      
      return(sequence)
      
    })
    
  } else {
    
    report_data <- sapply(retrieved_data, function(x){ return(unlist(x = x)[1])})
    
    if(info == "header"){
      
      report_data <- gsub(pattern = "^>", replacement = "", x = report_data)
      
    } else if(info == "name"){
      
      report_data <- gsub(pattern = "^\\S+\\s+", replacement = "", x = report_data)
      
    }  
    
  }
  
  #Close all connections. (Related to RStudio session termination)
  closeAllConnections()
  
  return(report_data)
  
}

################################################################################
## Function to retrieve signal or transit peptide info from Uniprot for a vector of protein accessions.
## Uses base R for accessions length < 30 and parallelization otherwise.
## Is is used in function: fill.protein.names.in.sa.tables().
################################################################################
retrieve.signal.peptide.from.uniprot.parallel <- function(accessions, cores){
  
  #INPUT:
  #accessions: a vector of protein accessions.
  #info: either "header", "name", "sequence" or "full" to retrieve the fasta header, protein name, 
  #      protein sequence or complete fasta entry, respectively.
  #cores: the number of cores to be used in parallelization.
  
  #OUTPUT:
  #report_data: the data type corresponding to "info" and for all the protein accessions.
  
  #Check for missing arguments.
  if(missing(x = accessions)){
    
    stop(crayon::bgRed(crayon::black("Accessions arguments is missing at retrieve.signal.peptide.from.uniprot.parallel()...\n")))
    
  }
  
  
  #Create the uniprot fasta urls per protein accession.
  url <- paste0('https://rest.uniprot.org/uniprotkb/', accessions, ".txt", collapse = NULL)
  
  #If there are more than 30 proteins use base R + parallel.
  if(length(x = accessions) >= 30){
    
    #If the cores argument is missing then assign number of cores detected minus one.
    if(missing(x = cores)){
      
      cores <- parallel::detectCores() - 1
      
    }
    
    #Similarly check is cores argument is higher than the available cores. If so leave one core available for other processes.
    detected_cores <- parallel::detectCores()
    
    if(detected_cores < cores){
      
      cluster <- parallel::makePSOCKcluster(names = detected_cores - 1)
      
    } else {
      
      cluster <- parallel::makePSOCKcluster(names = cores)
      
    }
    
    #Initiate core cluster.
    parallel::setDefaultCluster(cl = cluster)
    
    #Import the Uniprot urls to the cluster.
    parallel::clusterExport(cl = cluster, varlist = c("url"), envir = environment())
    
    #Retrieve fasta for the accessions in parallel.
    retrieved_data <- parallel::parLapply(cl = cluster, url, function(x){ 
      
      return(suppressWarnings(expr = try(expr = read.csv(file = x, header = F, sep = "\n"), silent = TRUE)))
      
    })
    
    #Stop the cores.
    parallel::stopCluster(cl = cluster)
    
  } else { #Else use base R only.
    
    #Retrieve fasta for the accessions with just read.csv function.
    retrieved_data <- lapply(url, function(x){ 
      
      return(suppressWarnings(expr = try(expr = read.csv(file = x, header = F, sep = "\n"), silent = TRUE)))
      
    })
    
  }
  
  #Assign NA to accessions that their fasta was not able to be retrieved.
  retrieved_data[sapply(retrieved_data, function(x){ return(class(x = x))}) == "try-error"] <- as.character(x = NA)
  
  sp <- sapply(retrieved_data, function(x){
    
    if(!all(is.na(x = x))){
      
      sp <- x[str_detect(string = x$V1, pattern = "FT   SIGNAL"),]
      
      if(identical(x = sp, y = character(0))){
        
        sp <- x[str_detect(string = x$V1, pattern = "FT   TRANSIT"),]
        
      }
      
      if(!identical(x = sp, y = character(0))){
        
        sp <- strsplit(x = sp, split = "          ")[[1]][2]
        
        # Splitting the string
        split_sp <- unlist(strsplit(sp, "\\.\\."))
        
        split_sp <- c(as.numeric(x = split_sp[1]), as.numeric(x = split_sp[2]))
        
      } else {
        
        split_sp <- NA
        
      }
      
    } else {
      
      split_sp <- NA
      
    }
    
    return(split_sp)
    
  })
  
  
  
  
  #Close all connections. (Related to RStudio session termination)
  closeAllConnections()
  
  return(sp)
  
}

############################################################################### #
# Function: Calculates the isoelectric point of a peptide sequence using the 
# Bjellqvist pka scale and weighted C and N terminal pkas for non-charged residues
# Algorithm inspired by expasy
############################################################################### #
Isoelectric.point <- function(seq, weighted) {
  # Input: 
  #   seq: An amino acid sequence
  #   weighted: Boolean deciding if C and N terminal should be weighted for non-charged residues
  #
  # Output: The isoelectric point of the provided sequence.
  
  #ref
  #Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
  #Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F. The focusing positions of polypeptides 
  #in immobilized pH gradients can be predicted from their amino acid sequences. Electrophoresis 1993,
  #14, 1023-1031.
  
  
  #Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E. Reference points for
  #comparisons of two-dimensional maps of proteins from different human cell
  #types defined in a pH scale where isoelectric points correlate with polypeptide
  #compositions. Electrophoresis 1994, 15, 529-539.
  
  if (missing(seq) || missing(weighted)) {
    stop(crayon::bgRed("Error: Arguments in Isoelectric.point() are missing!"))
  }
  
  # Ensure all characters are uppercase and no unwanted characters are present
  seq <- aaCheck(seq = seq)[[1]]
  
  # Count respective amino acids
  aa_count <- table(seq)
  aa_to_check <- c("C", "D", "E", "H", "K", "R", "Y")
  
  # Replace missing values with zeros
  aa_count[aa_to_check] <- replace(aa_count[aa_to_check], is.na(x = aa_count[aa_to_check]), 0)
  
  # Define variables and arguments
  cterm_pka <- 3.54
  Nterm_pka <- 7.5
  
  if (isTRUE(x = weighted)) {
    AA1 <- seq[1]
    
    Nterm_pka <- ifelse(AA1 %in% c("A", "M", "P", "S", "T", "V"), 
                        c("A" = 7.59, "M" = 7.0, "P" = 8.36, "S" = 6.93, "T" = 6.82, "V" = 7.44)[AA1],
                        Nterm_pka)
    
  }
  
  ph_inc <- seq(1, 15, 0.01)
  
  for (pH in ph_inc) {
    QN1 <- -1 / (1 + (10^(cterm_pka - pH)))
    QN2 <- -aa_count["D"] / (1 + (10^(4.05 - pH)))
    QN3 <- -aa_count["E"] / (1 + (10^(4.45 - pH)))
    QN4 <- -aa_count["C"] / (1 + (10^(9 - pH)))
    QN5 <- -aa_count["Y"] / (1 + (10^(10 - pH)))
    QP1 <- aa_count["H"] / (1 + (10^(pH - 5.98)))
    QP2 <- 1 / (1 + (10^(pH - Nterm_pka)))
    QP3 <- aa_count["K"] / (1 + (10^(pH - 10)))
    QP4 <- aa_count["R"] / (1 + (10^(pH - 12)))
    
    NQ <- QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4
    
    if (pH > 14) {
      cat(crayon::bgRed("Something is wrong. pH above 14.... \n"))
      break
    }
    
    if (NQ <= 0) {
      #cat(crayon::green(paste0("Isoelectric point found. pI: ", pH, "\n")))
      return(pH)
      break
    }
  }
}

############################################################################### #
# Function: Check as string for unwanted characters to ensure only accepted amino 
# acids are examined.
# Used in Isoelectric.point()
############################################################################### #
aaCheck <- function(seq){
  #Input:
  # - Seq: Amino acid sequence
  
  #Output: The same sequence as input, if the sequence passes checkpoint critera
  if(!any(lengths(x = seq) > 1)){
    
    seq <- toupper(x = seq)
    seq <- gsub(pattern = "[[:space:]]+",replacement = "",x = seq)
    seq <- strsplit(x = seq,split = "")
    
  } else {
    seq <- lapply(seq,toupper)
  }
  check <- unlist(lapply(seq,function(sequence){
    
    !all(sequence%in%c("A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,
                       "N" ,"P" ,"Q" ,"R" ,"S" ,"T" ,"V" ,"W" ,"Y", "-", "X",
                       "B", "Z", "U"))
    
  }))
  
  if(sum(check) > 0){
    
    sapply(which(check == TRUE),function(sequence){
      
      warning(paste0("Sequence ",sequence," has unrecognized amino acid types. Output value might be wrong calculated"),call. = FALSE)})
    
  }
  
  return(seq)
  
}

############################################################################### #
#Function: Calculates the theoretical molecular weight for the amino acid sequence
# provided. The molecular weight may vary based on the settings and what scale is chosen.
# it is also possible to add a amino acid shift. 
#Used in generate.matrix()
############################################################################### #
molecular.weight <- function(seq, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL) {
  # Split sequence by amino acids
  #seq <- check_string(input_string = toupper(seq))
  # Ensure all characters are uppercase and no unwanted characters are present
  #seq <- strsplit(x = check_string(input_string = toupper(seq)), split = "")[[1]]
  seq <- aaCheck(seq)
  # Create the weight scale
  if (monoisotopic == TRUE) {
    weight <-
      c(
        A = 71.03711,
        R = 156.10111,
        N = 114.04293,
        D = 115.02694,
        C = 103.00919,
        E = 129.04259,
        Q = 128.05858,
        G = 57.02146,
        H = 137.05891,
        I = 113.08406,
        L = 113.08406,
        K = 128.09496,
        M = 131.04049,
        F = 147.06841,
        P = 97.05276,
        S = 87.03203,
        T = 101.04768,
        W = 186.07931,
        Y = 163.06333,
        V = 99.06841,
        U = 150.95363,
        O = 237.14772,
        H2O = 18.01056
      )
  } else if (avgScale == "expasy"){
    weight <-
      c(
        A = 71.0788,
        R = 156.1875,
        N = 114.1038,
        D = 115.0886,
        C = 103.1388,
        E = 129.1155,
        Q = 128.1307,
        G = 57.0519,
        H = 137.1411,
        I = 113.1594,
        L = 113.1594,
        K = 128.1741,
        M = 131.1926,
        F = 147.1766,
        P = 97.1167,
        S = 87.0782,
        T = 101.1051,
        W = 186.2132,
        Y = 163.1760,
        V = 99.1326,
        U = 150.0388,
        O = 237.3018,
        H2O = 18.01524
      )
  } else if (avgScale == "mascot"){
    weight <-
      c(
        A = 71.0779,
        R = 156.1857,
        N = 114.1026,
        D = 115.0874,
        C = 103.1429,
        E = 129.114,
        Q = 128.1292,
        G = 57.0513,
        H = 137.1393,
        I = 113.1576,
        L = 113.1576,
        K = 128.1723,
        M = 131.1961,
        F = 147.1739,
        P = 97.1152,
        S = 87.0773,
        T = 101.1039,
        W = 186.2099,
        Y = 163.1733,
        V = 99.1311,
        U = 150.0379,
        O = 237.2982,
        H2O = 18.01528
      )
  }
  
  # Sum the weight of each amino acid and add H2O weight
  mass <- unlist(lapply(seq, function(seq) {
    sum(weight[c(seq, "H2O")], na.rm = TRUE)
  }))
  
  # Add massShift for labeled proteins
  mass <- mass + massShift(seq = seq, label = label, aaShift = aaShift, monoisotopic = monoisotopic)
  
  return(mass)
}

################################################################################ #
# Function: his function calculates the mass difference of peptides introduced
# by chemical modifications or heavy isotope labelling
################################################################################ #
massShift <- function(seq, label = "none", aaShift = NULL, monoisotopic = TRUE){
  #Input: 
  # - seq: The amino acid sequence
  # - label: Set a predefined heavy isotope label. Accepts "none", "silac_13c",
  #   "silac_13c15n" and "15n". Overwrites input in \code{aaShift}.
  # - monoisotopic: A logical value \code{'TRUE'} or \code{'FALSE'} indicating
  #   if monoisotopic weights of amino-acids should be used
  # - aashift: Define the mass difference in Dalton of given amino acids as a named vector. 
  
  #the 20 standard amino acids in upper case.
  aaList <- c("A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,
              "M" ,"N" ,"P" ,"Q" ,"R" ,"S" ,"T" ,"V" ,"W" ,"Y")
  # Check inputs
  label <- tolower(label) # Renders case-insensitive input string.
  if(!(label %in% c("none", "silac_13c", "silac_13c15n", "15n"))){
    stop("Given label type unknown. Please use one of 'none', '15N', 'Silac_13C15N', or 'Silac_13C' (case-insensitive).")
  }
  if(!is.null(aaShift) & is.null(names(aaShift))){
    stop("'aaShift' must be given as a named vector, e.g. 'aaShift = c(K = 6.020129)'.")
  }
  allowed <- c(aaList, "Cterm", "Nterm")
  if(!is.null(aaShift) & !all(names(aaShift) %in% allowed)){
    stop(paste("Unknown amino acids defined in 'aaShift'. Only the following names are allowed:", paste(allowed, collapse = ", ")))
  }
  
  # Predefined labels
  if (label == "silac_13c"){
    aaShift <- c("K" = 6.020129 - 0.064229*!monoisotopic, "R" = 6.020129- 0.064229*!monoisotopic)
  } else if(label == "silac_13c15n"){
    aaShift <- c("K" = 8.014199 -0.071499*!monoisotopic, "R" = 10.008269-0.078669*!monoisotopic)
  } else if(label == "15n"){
    aaShift <- c(
      #      U = 1,
      #      O = 3,
      A = 1,
      R = 4,
      N = 2,
      D = 1,
      C = 1,
      E = 1,
      Q = 2,
      G = 1,
      H = 3,
      I = 1,
      L = 1,
      K = 2,
      M = 1,
      F = 1,
      P = 1,
      S = 1,
      T = 1,
      W = 2,
      Y = 1,
      V = 1
    ) * 0.997035 -0.003635*!monoisotopic # 0.997035 equals the mass shift from 14N to 15N. 0.9934 equals the average mass shift.
  }
  
  # Split sequence by amino acids
  seq <- aaCheck(seq)
  
  # Calculate mass shifts
  unlist(
    lapply(seq, function(x){
      sum(aaShift[c(x)], aaShift["Nterm"], aaShift["Cterm"], na.rm = TRUE)
    }))
}