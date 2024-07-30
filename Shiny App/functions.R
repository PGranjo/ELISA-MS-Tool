

# List of required libraries
libraries <- c(
  "readxl", "dplyr", "seqinr", "stringr", "writexl", 
  "UpSetR", "ggVennDiagram", "ggplot2", "gridExtra", "grid"
)

# Packages to install (those not already installed)
to_install <- setdiff(libraries, rownames(installed.packages()))

# Install missing packages
if (length(to_install) > 0) {
  install.packages(to_install, repos = "http://cran.us.r-project.org")
}

# Load required libraries
lapply(libraries, library, character.only = TRUE)


################################################################################
################################ Input #########################################
################################################################################


#####################
#' @title Function: coverage_single_prot()
#'
#' @description
#' A function that retrieves the names of antibody kits detected for a given protein accession 
#' within a specified organism from a provided database.
#' 
#' @param database A list containing two data frames: `tcoverage` and `tsample`. 
#'                 `tcoverage` should include a column 'Uniprot' and 'id', while 
#'                 `tsample` should include 'id', 'organism', and 'ab_short_name'.
#' @param acession A character string representing the protein accession number (Uniprot ID).
#' @param org A character string specifying the organism name to filter the samples by.
#' 
#' @returns A vector containing the names of the detected antibody kits.
#' 
#' @examples
#' library(dplyr)
#' database <- list(
#'   tcoverage = data.frame(Uniprot = c("P12345", "Q67890"), id = c(1, 2)),
#'   tsample = data.frame(id = c(1, 2), organism = c("homo sapiens", "mus musculus"), ab_short_name = c("Kit A", "Kit B"))
#' )
#' 
#' kits <- coverage_single_prot(database, "P12345", "homo sapiens")
#' print(kits)
#' 
#' @export
#####################
coverage_single_prot <- function(database, acession, org){
  
  #Unlist the dataframes tcoverage and tsample from the database
  tcoverage <- database[["tcoverage"]]
  tsamples <- database[["tsample"]]
  
  #Filter the tcoverage table based on the Uniprot acession number and retrieves the id
  detected_kits <- tcoverage %>% filter(Uniprot == acession) %>% select(id)
  
  #converts all characters to lowercase and remove all whitespace from a string 
  org <- trimws(tolower(org))
  tsamples$organism <- trimws(tolower(tsamples$organism))
  
  # Check if the specified organism (org) exists in the 'organism' column of the 'tsamples' data frame.
  # If it does, filter 'tsamples' to find the rows where the 'id' matches those in 'detected_kits' 
  # and the 'organism' matches 'org'. Extract the unique 'ab_short_name' values (antibody kit names) 
  if (org %in% tsamples$organism) {
    kits_names <- as.vector(unlist(unique(tsamples %>% filter(id %in% unlist(detected_kits)) %>% filter(organism == org) %>% select(ab_short_name))))
  } else {
    #Otherwise sends a warning stating that there is no organism in the tsample table 
    warning(paste("The organism filter", org, "does not exist in the tsamples data."))
    stop(paste("The organism filter", org, "does not exist in the tsamples data."))
  }
  return(kits_names)
}


##########################################################
#' @title Function: read_fasta_file()
#'
#' @description
#' A function that reads a FASTA file and filters the sequences based on specified terminologies.
#' 
#' @param file_path The path to the FASTA file.
#' @param terminologies A character vector of terminologies to filter the sequences by. 
#'                      The function will keep sequences where the third part after '|' and '_' matches the terminologies.
#' 
#' @returns A vector of filtered sequences from the FASTA file.
#' 
#' @examples
#' 
#' 
#' # Define terminologies to filter by
#' terminologies <- c("SPOFR", "ECOBD")
#' 
#' # Path to your FASTA file
#' fasta_file_path <- "path/to/your/fasta_file.fasta"
#' 
#' # Read and filter the FASTA file
#' filtered_fasta_sequences <- read_fasta_file(fasta_file_path, terminologies)
#################################################################################




read_fasta_file <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("File does not exist.")
  }
  
  # Vector of terminologies that we want to filter in the Fast file names
  terminologies <- c("CRIGR", "ECOBD", "SPOFR")
  
  # Read the FASTA file
  fasta_content <- read.fasta(file = file_path, seqtype = "AA", as.string = TRUE)
  
  # Extract the names from the fasta content
  element_names <- names(fasta_content)
  
  # Split each element by '|'
  split_elements <- str_split(element_names, "\\|")
  
  #After the split there are 3 parts each were initially divided by "|"
  # Extract the second part (ID) and the third part (terminology)
  extracted_terminologies <- sapply(split_elements, function(x) {
    if (length(x) >= 3) {
      id <- x[2]
      terminology <- str_split(x[3], "_")[[1]]
      last_part <- terminology[length(terminology)]
      return(c(id, last_part))
    } else {
      return(c(NA, NA))
    }
  })
  
  # Convert to a data frame for easier manipulation
  extracted_df <- as.data.frame(t(extracted_terminologies), stringsAsFactors = FALSE)
  colnames(extracted_df) <- c("ID", "Terminology")
  
  # Filter elements based on the specified terminologies
  filtered_elements <- extracted_df$ID[extracted_df$Terminology %in% terminologies]
  
  return(filtered_elements)
}


#####################
#' @title Function: coverage_prot()
#'
#' @description
#' A function that processes a given input of protein accession numbers or a file containing them, 
#' and retrieves the detected antibody kits for each protein within a specified organism from a provided database.
#' 
#' @param database A list containing two data frames: `tcoverage` and `tsample`. 
#'                 `tcoverage` should include columns 'Uniprot' and 'id', while 
#'                 `tsample` should include 'id', 'organism', and 'ab_short_name'.
#' @param input A character vector, list, or file path. This can be a comma-separated string of accession numbers,
#'              a path to a FASTA file, or an Excel file containing a "Protein" column with accession numbers.
#' @param organism A character string specifying the organism name to filter the samples by.
#' 
#' @returns A list where each element contains the antibody kits detected for each protein, indexed by the protein accession number.
#' 
#' @examples
#' library(dplyr)
#' library(readxl)
#' database <- list(
#'   tcoverage = data.frame(Uniprot = c("P12345", "Q67890"), id = c(1, 2)),
#'   tsample = data.frame(id = c(1, 2), organism = c("homo sapiens", "mus musculus"), ab_short_name = c("Kit A", "Kit B"))
#' )
#' 
#' # Example using a comma-separated string of accession numbers
#' result <- coverage_prot(database, "P12345,Q67890", "homo sapiens")
#' print(result)
#' 
#' # Example using an Excel file (path must be adjusted to your file location)
#' # result <- coverage_prot(database, "path/to/file.xlsx", "homo sapiens")
#' # print(result)
#' 
#' @export
#####################


coverage_prot <- function(database, input, organism) {
  # Initialize accession_numbers
  accession_numbers <- NULL
  
  # Check if the input is a file that exists
  if (is.character(input) && file.exists(input)) {
    #Extracts file extensio
    file_extension <- tools::file_ext(input)
    
    #Checks if the extension is fasta
    if (file_extension == "fasta") {
    #Process the acession numbers from the fasta file
      accession_numbers <- read_fasta_file(input)
    #Check if the extension is xlsx
    } else if (file_extension == "xlsx") {
      #Reads the excel file and looks for the HCP sheet
      excel_data <- readxl::read_xlsx(input, sheet = "HCP", .name_repair = "minimal")
      #Check if Protein is a colname within the excel file
      if (!"Protein" %in% colnames(excel_data)) {
        stop("The Excel file does not contain a 'Protein' column.")
      } else {
        # Extract accession numbers from the "Protein" column
        accession_numbers <- unique(apply(excel_data["Protein"], 1, function(x) strsplit(x, split = "\\|")[[1]][2]))
      }
    } else {
      stop("Unsupported file format. Please provide a FASTA or Excel file.")
    }
  } else if (is.character(input) && !file.exists(input)) {
    # If input is a character string that is not a file path. The acession number are splited based on white space characters and new lines
    accession_numbers <- unlist(strsplit(input, ",\\s*"))
    accession_numbers <- unlist(strsplit(accession_numbers, "\n"))
    accession_numbers <- trimws(accession_numbers)
  } else if (is.vector(input) || is.list(input)) {
    # If input is a list or vector of accession numbers
    accession_numbers <- unlist(input)
  } else {
    stop("Unsupported input type. Please provide a file path or a list of accession numbers.")
  }
  
  # Process each accession number and collect the results
  results_list <- list()
  for (prot in accession_numbers) {
    #Process each Acession Number and collects a vector of kit(s) that each is coverage
    result <- coverage_single_prot(database, prot, organism)
    #Save kit information in a list where each acession number is a name from the list 
    results_list[[prot]] <- result
    }
  return(results_list)
}


#####################
#' @title Function: hcp_filtration_number()
#'
#' @description
#' A function that filters and retrieves unique protein accession numbers based on the 'ppm' values
#' in the provided dataset, either by a specified number of top entries or by a percentage of the dataset.
#' 
#' @param data A data frame containing at least a 'ppm' column and a 'Protein' column.
#' @param number A numeric value indicating the number of top proteins to select, or a percentage (if the else part is runned).
#' @param type A character string indicating the type of filtration: "numeric" for a specific number of top entries.
#' 
#' @returns A character vector containing the unique protein accession numbers based on the specified filtration method.
#' 
#' @examples
#' # Assuming 'data' is a data frame with columns 'ppm' and 'Protein'
#' data <- data.frame(
#'   ppm = c(1.5, 2.3, 0.8, 4.5),
#'   Protein = c("sp|P12345|Protein1", "sp|Q67890|Protein2", "sp|A12345|Protein3", "sp|B67890|Protein4")
#' )
#' 
#' # Example using numeric filtration
#' top_proteins <- hcp_filtration_number(data, 2, "numeric")
#' print(top_proteins)
#' 
#' # Example using percentage filtration
#' top_proteins_percentage <- hcp_filtration_number(data, 50, "percentage")
#' print(top_proteins_percentage)
#' 
#' @export
#####################

hcp_filtration_number <- function(data, number, type) {
  #Checks if the ppm column exist within the table
  if (!"ppm" %in% colnames(data)) {
    stop("The data does not contain a 'ppm' column.")
  }
  
  # Use order to sort the data by 'ppm' column in descending order
  sorted_data <- data[order(-data$ppm), ]
  
  if (type == "numeric") {
    #Pick of the acession name base on the second element from a split string on the |
    accessions <- unique(sapply(sorted_data$Protein[1:number], function(x) strsplit(x, split = "|", fixed = TRUE)[[1]][2]))
    return(accessions)
  } else {
    #Use a percentage instead where it's fraction is first calculated
    fraction <- number / 100
    #Round the value of the number of hcp filtered
    target_number <- round(nrow(sorted_data) * fraction, 0)
    accessions <- unique(sapply(sorted_data$Protein[1:target_number], function(x) strsplit(x, split = "|", fixed = TRUE)[[1]][2]))
    return(accessions)
  }
}


################################################################################
############################## Output ##########################################
################################################################################



#####################
#' @title Function: count_kits()
#'
#' @description
#' A function that calculates the frequency and percentage of different kit types from a list of protein data.
#' 
#' @param protein_list A list containing protein data, where each entry represents a kit type associated with a protein.
#' 
#' @returns A data frame with three columns: 'kit' for kit types, 'freq' for the frequency of each kit type,
#'          and 'perc' for the percentage each kit type represents in the total list.
#' 
#' @examples
#' # Sample protein list with kit types
#' protein_list <- c("Kit A", "Kit B", "Kit A", "Kit C", "Kit A", "Kit B")
#' 
#' # Calculate the frequency and percentage of each kit type
#' kit_summary <- count_kits(protein_list)
#' print(kit_summary)
#' 
#' @export
#####################
count_kits <- function(protein_list) {
  kits_info <- list()
  kit_types <- c()
  
  # Collect all kit types from the protein list
  for (protein in protein_list) {
    kit_types <- c(kit_types, protein)
  }
  # Calculate the frequency of each kit type
  kit_counts <- table(kit_types)
  # Create a data frame with kit types, their frequency, and their percentage in the total list
  kits_info <- data.frame("kit" = names(kit_counts), "freq" = as.numeric(kit_counts), "perc" = round(as.numeric(kit_counts)/length(protein_data)*100,1))
  return(kits_info)
}


#####################
#' @title Function: ic_hcp_table()
#'
#' @description
#' A function that generates a matrix indicating the detection of proteins by different kits for a specified organism. 
#' Each row represents a protein, and each column represents a unique kit that detected the protein. 
#' The number 1 means that the protein is immunocaptured by the kit
#' 
#' @param database A list containing data frames, including `tsample` which must contain 'organism' and 'ab_short_name' columns.
#' @param protein_list A named list where each name is a protein identifier and each entry is a vector of kit names that detected the protein.
#' @param organism A character string specifying the organism name to filter the samples by.
#' 
#' @returns A binary matrix where rows correspond to proteins and columns correspond to unique kits. A value of 1 indicates the protein was detected by the kit.
#' 
#' @examples
#' # Example database setup
#' database <- list(
#'   tsample = data.frame(
#'     organism = c("homo sapiens", "homo sapiens", "mus musculus"),
#'     ab_short_name = c("Kit A", "Kit B", "Kit C")
#'   )
#' )
#' 
#' # Example protein list with kits that detected each protein
#' protein_list <- list(
#'   "Protein1" = c("Kit A", "Kit B"),
#'   "Protein2" = c("Kit A"),
#'   "Protein3" = c()
#' )
#' 
#' # Generate the binary detection matrix for 'homo sapiens'
#' detection_matrix <- ic_hcp_table(database, protein_list, "homo sapiens")
#' print(detection_matrix)
#' 
#' @export
###################

ic_hcp_table <- function(database, protein_list, organism) {
  # Filter the tsample table to include only entries from the specified organism
  tsamples <- database[["tsample"]] %>% 
    filter(organism == organism) %>%
    distinct(ab_short_name)  # Ensure all unique kits for the specified organism are considered
  
  # Initialize the matrix with NA or zero
  n_kits_database <- nrow(tsamples)
  infopedia <- matrix(0, nrow = length(protein_list), ncol = n_kits_database)
  rownames(infopedia) <- names(protein_list) #rownames are all acession numbers collected from the initial user input
  colnames(infopedia) <- tsamples[,"ab_short_name"] #colnames are all available kits 
  
  # Fill the matrix with detection information
  count <- 1
  for (protein in names(protein_list)) {
    kit_names <- protein_list[[protein]]  # Kits that detected this protein
    if(length(kit_names) > 0) {
      col_idx <- which(colnames(infopedia) %in% kit_names) #based on the available kits name we get a col index
      infopedia[count, col_idx] <- 1 #A 1 is added to those col_indices
    }
    count <- count + 1
  }
  infopedia <- infopedia[, colSums(infopedia) != 0]
  return(infopedia)
}


################################################################################
########################## Venn Diagram ########################################
################################################################################

#####################
#' @title Function: create_set_list()
#'
#' @description
#' A function that creates a list of sets from a binary matrix, where each set contains the names of rows
#' corresponding to a value of 1 in each column.
#' 
#' @param table A binary matrix with row names representing items (e.g., proteins) and column names representing categories or conditions (e.g., kits).
#' 
#' @returns A named list where each name corresponds to a column from the input matrix, and each element is a vector of row names where the matrix has a value of 1 for that column.
#' 
#' @examples
#' # Example binary matrix setup
#' binary_matrix <- matrix(c(1, 0, 1, 1, 0, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' rownames(binary_matrix) <- c("Protein1", "Protein2", "Protein3")
#' colnames(binary_matrix) <- c("Kit A", "Kit B", "Kit C")
#' 
#' # Create a list of sets from the binary matrix
#' set_list <- create_set_list(binary_matrix)
#' print(set_list)
#' 
#' @export
#####################

create_set_list <- function(table) {
  # Initialize an empty list to store the sets
  set_list <- lapply(1:ncol(table), function(i) {
    # For each column, retrieve the row names where the value is 1 (where the protein was immunocaptured)
    rownames(table)[table[, i] == 1]
  })
  # Name the elements of the list according to the column names of the table which are the kits 
  names(set_list) <- colnames(table)
  return(set_list)
}


#####################
#' @title Function: create_venn_diagram()
#'
#' @description
#' A function that generates a Venn diagram from a list of sets, highlighting the intersections among the specified comparisons.
#' 
#' @param set_list A named list where each element which corresponds to a vector of items (acession protein numbers), representing different sets of immunocaptured proteins by Kit.
#' @param comparisons A character vector specifying the names of the sets to be compared in the Venn diagram. This would be the kits short names
#' 
#' @returns Returns a ggplot object representing the Venn diagram. 
#' 
#' @examples
#' # Example set list
#' set_list <- list(
#'   "Group A" = c("A", "B", "C"),
#'   "Group B" = c("B", "C", "D"),
#'   "Group C" = c("C", "D", "E")
#' )
#' 
#' # Create and display the Venn diagram for specified comparisons
#' venn_plot <- create_venn_diagram(set_list, comparisons = c("Group A", "Group B", "Group C"))
#' print(venn_plot)
#' 
#' @export
#####################
create_venn_diagram <- function(set_list, comparisons, filename = NULL) {
  # Ensure comparisons are valid
  if (!all(comparisons %in% names(set_list))) {
    stop("Some comparison names are not in the set list")
  }
  
  # Subset the set list based on the comparisons
  subset_list <- set_list[comparisons]
  
  # Generate the Venn diagram
  venn.plot <- ggVennDiagram(
    subset_list,
    set_size = 6,
    label = "both",
    label_percent_digit = 1,
    label_size = 5
  ) + 
    scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = "white", high = "green4") +
  theme(
      plot.margin = margin(20, 20, 20, 20),  # Adjust margins
      plot.title = element_text(size = 14, face = "bold"),  # Adjust title size if needed
      plot.background = element_rect(fill = "transparent", color = NA),
      text = element_text(size = 16), 
      legend.text = element_text(size = 16)
    )    
  return(venn.plot)
}



################################################################################
############################## Upset Plot ######################################
################################################################################

#####################
#' @title Function: calculate_input_percent()
#'
#' @description
#' A function that calculates the percentage of unique items (e.g., proteins) detected by combinations of kits 
#' from a given set list, relative to the total number of items in the dataset. Input preparation for a Upset Plot
#' 
#' @param set_list A named list where each element which corresponds to a kit contains a vector of items, representing different sets of immunocaptured proteins.
#' @param table A binary matrix with row names representing items (e.g., proteins) and column names representing categories or conditions (e.g., kits).
#' 
#' @returns A named vector with the names representing combinations of kit names (concatenated with '&') and 
#'          values representing the percentage of unique items detected by each combination.
#' 
#' @examples
#' # Example set list
#' set_list <- list(
#'   "Kit A" = c("Protein1", "Protein2"),
#'   "Kit B" = c("Protein2", "Protein3"),
#'   "Kit C" = c("Protein1", "Protein3", "Protein4")
#' )
#' 
#' # Example table representing the total number of items (e.g., proteins)
#' table <- data.frame(
#'   Protein = c("Protein1", "Protein2", "Protein3", "Protein4", "Protein5")
#' )
#' 
#' # Calculate the percentage of unique items detected by each kit combination
#' percentages <- calculate_input_percent(set_list, table)
#' print(percentages)
#' 
#' @export
#####################
calculate_input_percent <- function(set_list, table) {
  kit_names <- names(set_list) # Extract the names of the kits
  input <- c() # Initialize an empty vector to store the counts
  
  for (k in 1:length(kit_names)) {
    combs <- combn(kit_names, k, simplify = FALSE) # Generate all k-combinations of kit names
    for (comb in combs) {
      name <- str_c(comb, collapse = "&") # Concatenate the combination names with '&'
      n <- length(unique(unlist(set_list[comb]))) # Count the unique items detected by this combination
      input <- c(input, n)  # Append the count to the input vector
      names(input)[length(input)] <- name  # Assign the combination name to the corresponding count
    }
  }
  # Calculate the percentage of unique items detected by each combination
  input_percent <- round((input / nrow(table)) * 100, 2)
  return(input_percent)
}




#####################
#' @title Function: create_upset_diagram()
#'
#' @description
#' A function that generates an UpSet diagram to visualize the intersections and coverage percentages of sets 
#' from the provided input percentages.
#' 
#' @param input_percent A named vector containing the percentage of unique items detected by each combination of sets.
#'                      The names of the vector should represent the combinations of sets, concatenated with '&'.
#' 
#' @returns A grid graphics object representing the UpSet plot.
#' 
#' @examples
#' # Example input percentage vector
#' input_percent <- c("Kit A&Kit B" = 40, "Kit A&Kit C" = 30, "Kit B&Kit C" = 20, "Kit A&Kit B&Kit C" = 10)
#' 
#' # Create the UpSet diagram
#' upset_plot <- create_upset_diagram(input_percent)
#' grid::grid.draw(upset_plot)  # Display the UpSet diagram
#' 
#' @export
#####################
create_upset_diagram <- function(input_percent) {
  
  # Generate the UpSetR plot
  ups <- upset(
    fromExpression(input_percent),
    keep.order = TRUE,
    mb.ratio = c(0.7, 0.3),  # Adjust the ratio to allocate more space to the main bar plot
    number.angles = 0, 
    text.scale = c(2, 2.5, 1.1, 1.1, 2.1, 2), 
    point.size = 3, 
    line.size = 1.2,
    mainbar.y.label = "Coverage (%)",
    sets.x.label = "",
    set_size.show = FALSE,  # Hide the set size bar chart
    set_size.numbers_size = 0,  # Hide the set size numbers
    sets.bar.color = "white",  # Set color for set bars
  )
  # Function to skip set size plot in the UpSet diagram
  skip_set_size_plot <- function(ups) {
    main <- ups$Main_bar
    # Find the panel grob
    panel_idx <- grep("panel", main$layout$name, fixed = TRUE)
    # Find the text grob
    text_idx <- which(
      vapply(main$grobs[[panel_idx]]$children,
             function(x) inherits(x, "text"),
             logical(1)))
    tG <- main$grobs[[panel_idx]]$children[[text_idx]]
    
    main$grobs[[panel_idx]]$children[[text_idx]] <- tG
    # Adjust the heights ratio here
    return(arrangeGrob(main, ups$Matrix, heights = unit(c(3, 1), "null")))
  }
  upset_plot <- skip_set_size_plot(ups)
  return(upset_plot)
}




#####################################################################################################
#################################### PI and MW section ##############################################
#####################################################################################################


source("PI_MW_functions.R")

#####################
#' @title Function: process_accessions()
#'
#' @description
#' A function that processes protein accession numbers to retrieve signal peptides and sequences from UniProt.
#' The sequences are then processed to remove the signal peptide regions.
#' 
#' @param accessions A vector of protein accession numbers.
#' 
#' @returns A vector containing the processed sequences for each accession, with signal peptides removed.
#' 
#' @examples
#' # Example list of accession numbers
#' accessions <- c("P12345", "Q67890", "A12345")
#' 
#' # Process the accessions to retrieve and clean sequences
#' processed_sequences <- process_accessions(accessions)
#' print(processed_sequences)
#' 
#' @export
#####################
process_accessions <- function(accessions) {
  # Retrieve signal peptides for each accession from UniProt
  sp <- retrieve.signal.peptide.from.uniprot.parallel(accessions = accessions)
  names(sp) <- accessions
  
  # Retrieve sequences for each accession from UniProt
  seq <- retrieve.from.uniprot.parallel(accessions = accessions, info = "sequence")
  names(seq) <- accessions
  
  # Process each accession to remove signal peptides from the sequences
  processed_seq <- unlist(lapply(accessions, function(x) {
    sp_seq <- sp[[x]]  # Get the signal peptide for the current accession
    seq_seq <- seq[[x]] # Get the sequence for the current accession
    
    # Return NA if the sequence is not available
    if (is.na(seq_seq)) {
      return(NA)
    }
    
    # If signal peptide information is available, remove the signal peptide from the sequence
    if (!all(is.na(sp_seq))) {
      seq_seq <- substring(seq_seq, first = (sp_seq[2] + 1), last = nchar(seq_seq))
    }
    
    return(seq_seq)
  }))
  
  return(processed_seq)
}


#####################
#' @title Function: calculate_protein_properties()
#'
#' @description
#' A function that calculates the isoelectric point (pI) and molecular weight (MW) for a list of protein sequences.
#' 
#' @param sequences A vector of protein sequences.
#' 
#' @returns A list with two named elements: `pI` for the isoelectric points and `MW` for the molecular weights (in kDa).
#' 
#' @examples
#' # Example sequences
#' sequences <- c("Protein1" = "MKTAYIAKQRQISFVKSHFSRQDILDLWQRAP", "Protein2" = "MKAKTFFQEALDAAGDNISAAELEKVT")
#' 
#' # Calculate the pI and MW for the given sequences
#' properties <- calculate_protein_properties(sequences)
#' print(properties)
#' 
#' @export
#####################
calculate_protein_properties <- function(sequences) {
  # Calculate isoelectric point (pI) for each sequence
  pI <- unlist(lapply(sequences, function(x) {
    if (is.na(x)) {
      return(NA)
    }
    
    return(round(Isoelectric.point(seq = x, weighted = TRUE), digits = 1))
  }))
  
  # Calculate molecular weight (MW) for each sequence
  MW <- unlist(lapply(sequences, function(x) {
    if (is.na(x)) {
      return(NA)
    }
    
    return(round(molecular.weight(seq = x), digits = 1) / 1000)  # to kDA
  }))
  
  return(list(pI = pI, MW = MW))
}



#####################
#' @title Function: main_function()
#'
#' @description
#' A function that generates a data frame containing protein accession numbers, sequences, isoelectric points (pI), 
#' molecular weights (MW), and protein names.
#' 
#' @param accessions A character vector of protein accession numbers (Uniprot IDs).
#' 
#' @returns A data frame with columns for accession numbers, pI, MW, and protein names.
#' 
#' @examples
#' # Example accession numbers
#' accessions <- c("P12345", "Q67890")
#' 
#' # Generate the main data frame with all relevant protein information
#' result <- main_function(accessions)
#' print(result)
#' 
#' @export
#####################
main_function <- function(accessions) {
  # Initialize a data frame with the accession numbers
  Gen_matrix <- data.frame("Accession" = accessions)
  
  # Process the accessions to retrieve sequences with signal peptides removed
  sequences <- process_accessions(accessions)
  
  # Calculate protein properties (pI and MW) for the sequences
  protein_properties <- calculate_protein_properties(sequences)
  
  # Add the calculated pI and MW to the Gen_matrix data frame
  Gen_matrix$pI <- protein_properties$pI
  Gen_matrix$MW <- protein_properties$MW
  
  # Retrieve protein names from Uniprot
  Gen_matrix$Protein_name <- retrieve.from.uniprot.parallel(accessions = accessions, info = "name")
  Gen_matrix$Protein_name <- unlist(lapply(Gen_matrix$Protein_name, function(x) {
    return(str_split(string = x, pattern = " OS=")[[1]][1])
  }))
  
  return(Gen_matrix)
}



#####################
#' @title Function: retrieve_protein_info()
#'
#' @description
#' A function that retrieves and formats protein information, including pI and MW, for a given list of accession numbers.
#' 
#' @param accessions A character vector of protein accession numbers (Uniprot IDs).
#' 
#' @returns A data frame with columns for accession numbers, pI, MW (in kDa), and protein names.
#' 
#' @examples
#' # Example accession numbers
#' accessions <- c("P12345", "Q67890")
#' 
#' # Retrieve protein information
#' protein_info <- retrieve_protein_info(accessions)
#' print(protein_info)
#' 
#' @export
#####################
retrieve_protein_info <- function(accessions) {

  # Get protein information using main_function
  protein_info <- main_function(accessions)
  
  # Select only required columns and rename them
  protein_info <- protein_info[, c("Accession", "pI", "MW", "Protein_name")]
  colnames(protein_info) <- c("Accession", "pI", "MW_Kda", "Protein_name")

  return(protein_info)
}

#####################
#' @title Function: retrieve_comprehensive_protein_info()
#'
#' @description
#' A function that combines protein detection information from an IC-HCP table with comprehensive protein data,
#' including pI, MW, and protein names.
#' 
#' @param ic_hcp_table A binary matrix indicating protein detection by different kits, with proteins as rows and kits as columns.
#' 
#' @returns A data frame combining the IC-HCP table with additional protein information (pI, MW, Protein name).
#' 
#' @examples
#' # Example IC-HCP table
#' ic_hcp_table <- matrix(c(1, 0, 1, 1, 0, 1), nrow = 3, dimnames = list(c("P12345", "Q67890", "A12345"), c("Kit1", "Kit2")))
#' 
#' # Retrieve comprehensive protein information
#' comprehensive_info <- retrieve_comprehensive_protein_info(ic_hcp_table)
#' print(comprehensive_info)
#' 
#' @export
#####################
retrieve_comprehensive_protein_info <- function(ic_hcp_table) {
 
  # Get protein information using retrieve_protein_info
  protein_info <- retrieve_protein_info(rownames(ic_hcp_table))
  rownames(protein_info) <- protein_info[,"Accession"]
  
  
  # Combine IC-HCP detection information with protein properties
  comprehensive_table <- cbind(ic_hcp_table,protein_info[rownames(ic_hcp_table),-1])

  return(comprehensive_table)
}


#####################
#' @title Function: visualize_data()
#'
#' @description
#' A function that visualizes protein data based on isoelectric point (pI) and molecular weight (MW), highlighting detection by different kits.
#' 
#' @param comprehensive_table A data frame containing comprehensive protein data including pI, MW, and detection by kits.
#' @param kit1_name A string specifying the name of the first kit for comparison.
#' @param kit2_name Optional. A string specifying the name of the second kit for comparison. If not provided, only the first kit's data is plotted.
#' 
#' @returns A ggplot object visualizing the data. The function also prints the plot.
#' 
#' @examples
#' # Example comprehensive table with protein data
#' comprehensive_table <- data.frame(
#'   pI = c(5.6, 6.2, 7.4),
#'   MW_Kda = c(50, 75, 100),
#'   Kit1 = c(1, 0, 1),
#'   Kit2 = c(0, 1, 1),
#'   row.names = c("P12345", "Q67890", "A12345")
#' )
#' 
#' # Visualize data for a single kit
#' visualize_data(comprehensive_table, "Kit1")
#' 
#' # Visualize data for two kits
#' visualize_data(comprehensive_table, "Kit1", "Kit2")
#' 
#' @export
#####################
visualize_data <- function(comprehensive_table, kit1_name, kit2_name = NULL) {
  # Extract data for the first kit
  data1 <- comprehensive_table[which(comprehensive_table[,c(kit1_name)] == 1),c("pI","MW_Kda")]
  data1 <- data.frame(data1, "Accession" = rownames(comprehensive_table)[which(comprehensive_table[,c(kit1_name)] == 1)])
  data1$Kit <- kit1_name
  if (is.null(kit2_name)) {
    # Plot data for a single kit
    p1 <- ggplot(data = data1, aes(x = pI, y = MW_Kda, color = Kit)) +
      geom_point(size = 4) +
      scale_color_manual(values = setNames(c("darkcyan"), kit1_name)) +
      scale_y_continuous(trans = "log2", name = "Theoretical MW [kDa]",
                         breaks = c(1, 10, 100, 1000)) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), name = "Theoretical pI") +
      expand_limits(x = c(2, 14), y = c(1, 1000)) +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = .1, color = "darkgrey"),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white"),legend.position = "none",
            panel.border = element_blank())
    
    print(p1)
    return(p1)
  } else {
    # Extract data for the second kit and combine with the first kit's data
    data2 <- comprehensive_table[which(comprehensive_table[,c(kit2_name)] == 1),c("pI","MW_Kda")]
    data2 <- data.frame(data2, "Accession" = rownames(comprehensive_table)[which(comprehensive_table[,c(kit2_name)] == 1)])
    data2$Kit <- kit2_name
    
    combined_data <- rbind(data1, data2)
    
    #Status column is create to differentiate in colour hcp which are exclusively immunocaptured either by one kit or both
    combined_data$Status <- with(combined_data, ifelse(Accession %in% data1$Accession & Accession %in% data2$Accession, 
                                                       "Common", ifelse(Accession %in% data1$Accession, "Kit1", "Kit2")))
    
    combined_data$Status <- factor(combined_data$Status, levels = c("Kit1", "Kit2", "Common"))
    
    # Plot data for both kits
    p2 <- ggplot(data = combined_data, aes(x = pI, y = MW_Kda, color = Status)) +
      geom_point(size = 4) +
      scale_color_manual(values = c(Kit1 = "darkcyan", Kit2 = "tomato2", Common = "chartreuse"),
                         labels = c("Kit1" = paste(kit1_name, "exc_hcp", sep = "_"), "Kit2" = paste(kit2_name, "exc_hcp", sep = "_"), "Common" = "Common_hcp")) +
      scale_y_continuous(trans = "log2", name = "Theoretical MW [kDa]",
                         breaks = c(1, 10, 100, 1000)) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14), name = "Theoretical pI") +
      expand_limits(x = c(2, 14), y = c(1, 1000)) +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth = .1, color = "darkgrey"),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white"),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            panel.border = element_blank()) +
      labs(color = "Kit Status")
    
    print(p2)
    return(p2)
  }
}



#################################################################################
############################## Excel File #######################################
##################################################################################


#' @title Function: convert_binary_table()
#'
#' @description
#' A function that converts a binary matrix (containing 0s and 1s) into a matrix with '+' for 1 and '-' for 0.
#' 
#' @param table A binary matrix (0s and 1s) with rows representing proteins and columns representing different conditions or kits.
#' 
#' @returns A data frame with the converted values and an additional 'Uniprot' column containing the row names.
#' 
#' @examples
#' # Example binary table
#' binary_table <- matrix(c(1, 0, 1, 0), nrow = 2, dimnames = list(c("P12345", "Q67890"), c("Kit1", "Kit2")))
#' 
#' # Convert the binary table
#' converted_table <- convert_binary_table(binary_table)
#' print(converted_table)
#' 
#' @export
#####################
convert_binary_table <- function(table) {
  # Convert 1 to '+' and 0 to '-'
  converted_table <- ifelse(table == 1, "+", "-")
  converted_table <- data.frame(Uniprot = rownames(table), converted_table)
  return(as.data.frame(converted_table))
}

#####################
#' @title Function: create_info_sheet()
#'
#' @description
#' A function that creates an information sheet summarizing the number and percentage of immunocaptured proteins per kit.
#' 
#' @param table A binary matrix with rows representing proteins and columns representing different kits.
#' 
#' @returns A data frame with columns for 'Kit', 'IC' (immunocaptured proteins), and 'IC_percent' (percentage of total proteins).
#' 
#' @examples
#' # Example binary table
#' binary_table <- matrix(c(1, 0, 1, 1), nrow = 2, dimnames = list(c("P12345", "Q67890"), c("Kit1", "Kit2")))
#' 
#' # Create the information sheet
#' info_sheet <- create_info_sheet(binary_table)
#' print(info_sheet)
#' 
#' @export
#####################
create_info_sheet <- function(table) {
  #total number of immunocaptured proteins per kit
  total_immunocaptured <- colSums(table)
  
  #Coverage percentage
  percentage <- round(total_immunocaptured / nrow(table) * 100, 2)
  
  info_sheet <- data.frame(
    Kit = colnames(table),
    IC = total_immunocaptured,
    IC_percent = percentage
  )
  
  return(info_sheet)
}


#####################
#' @title Function: write_combined_excel()
#'
#' @description
#' A function that prepares data for export to an Excel file, combining a binary table with protein information.
#' 
#' @param table A binary matrix with rows representing proteins and columns representing different kits.
#' 
#' @returns A list containing two data frames: the combined data (Coverage_table) and an information sheet (Information_Sheet).
#' 
#' @examples
#' # Example binary table
#' binary_table <- matrix(c(1, 0, 1, 1), nrow = 2, dimnames = list(c("P12345", "Q67890"), c("Kit1", "Kit2")))
#' 
#' # Prepare data for Excel
#' table_list <- write_combined_excel(binary_table)
#' print(table_list)
#' 
#' @export
#####################
write_combined_excel <- function(table) {
  
  # Convert the binary table into + and -
  converted_table <- convert_binary_table(table)
  
  # Create the information sheet
  info_sheet <- create_info_sheet(table)
  
  # Get unique accessions from the table
  unique_accessions <- rownames(table)
  
  # Retrieve protein information
  protein_info <- retrieve_protein_info(unique_accessions)
  
  converted_table <- as.data.frame(converted_table)
  converted_table$Accession <- rownames(converted_table)
  
  combined_table <- dplyr::full_join(converted_table, protein_info, by = "Accession")

  # Ensure Accession is set as the rownames again
  rownames(combined_table) <- combined_table$Accession
  combined_table$Accession <- NULL
  
  # Create the list of data frames
  table_list <- list(
    Coverage_table = combined_table,
    Information_Sheet = info_sheet
  )

  return(table_list)
}


#####################
#' @title Function: validate.plot.dimensions()
#'
#' @description
#' A function to validate and set appropriate dimensions for different types of plots.
#' 
#' @param type A string specifying the type of plot, either "rectangular" (e.g., violin plots, box plots) or "pie".
#' @param width A numeric value specifying the width of the plot in inches.
#' @param height A numeric value specifying the height of the plot in inches.
#' 
#' @returns A named vector with the validated width and height values.
#' 
#' @examples
#' # Validate dimensions for a rectangular plot
#' dimensions <- validate.plot.dimensions(type = "rectangular", width = 25, height = 15)
#' print(dimensions)
#' 
#' @export
#####################
validate.plot.dimensions <- function(type, width, height){
  
  #INPUT:
  #type: either "rectangular" (violin plots, box plots, line plots etc) or "pie" (pie-charts).
  #width: width of the plot expressed in inches (numeric).
  #height: height of the plot expressed in inches (numeric).
  
  #OUTPUT:
  #width: the validated width value.
  #height: the validated height value.
  
  #Validate type argument.
  if(missingArg(symbol = type)){
    
    type <- "rectangular"
    
  } else {
    
    if(!type %in% c("rectangular", "pie")){
      
      type <- "rectangular"
      
    }
    
  }
  
  #Validate width and height when the plot is of "rectangular" type.
  if(type == "rectangular"){
    
    if(any(is.na(x = c(width, height)), na.rm = TRUE) | width > 30 | width < 5 | height > 20 | height < 5){
      
      warning(crayon::bgRed(crayon::black(paste0("The width or height arguments have invalid values in validate.plot.dimensions(). ",
                                                 "Assigned default, width = 20.5 inches, height = 10.5 inches.", collapse = NULL))))
      
      width <- 14
      height <- 8
      
    }
    
    #Validate width and height when the plot is of "pie" type.
  } else {
    
    if(any(is.na(x = c(width, height)), na.rm = TRUE) | width > 20 | width < 5 | height > 20 | height < 5){
      
      warning(crayon::bgRed(crayon::black(paste0("The width or height arguments have invalid values in validate.plot.dimensions(). ",
                                                 "Assigned default, width = 10.5 inches, height = 10.5 inches.", collapse = NULL))))
      
      width <- 10
      height <- 8
      
    }
    
  }
  
  return(c(width = width, height = height))
  
}

