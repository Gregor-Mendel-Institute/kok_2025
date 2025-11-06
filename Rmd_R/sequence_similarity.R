# Load necessary libraries
library(msa)        # For multiple sequence alignment
library(ggplot2)    # For plotting

# Function to calculate similarity score
calculate_similarity <- function(seq1, seq2) {
  # Calculate the number of matching residues
  matches <- sum(seq1 == seq2)
  # Return the percentage similarity
  return(matches / length(seq1) * 100)
}

# Function to compute sliding window similarity
sliding_window_similarity <- function(aligned_seq1, aligned_seq2, window_size) {
  # Ensure aligned sequences are of the same length
  if (nchar(aligned_seq1) != nchar(aligned_seq2)) {
    stop("Aligned sequences must be of the same length.")
  }
  
  # Convert aligned sequences to vectors of characters
  seq1 <- unlist(strsplit(as.character(aligned_seq1), ""))
  seq2 <- unlist(strsplit(as.character(aligned_seq2), ""))
  
  # Initialize vectors for similarity and residue numbers
  similarities <- numeric()
  residue_numbers <- numeric()
  
  # Loop through the sequence
  for (i in 1:(length(seq1) - window_size + 1)) {
    # Extract the current window
    window_seq1 <- seq1[i:(i + window_size - 1)]
    window_seq2 <- seq2[i:(i + window_size - 1)]
    
    # Skip windows with gaps ('-')
    if ('-' %in% window_seq1 || '-' %in% window_seq2) {
      next
    }
    
    # Calculate similarity for the current window
    similarity <- calculate_similarity(window_seq1, window_seq2)
    
    # Store results
    similarities <- c(similarities, similarity)
    residue_numbers <- c(residue_numbers, i + (window_size - 1) / 2)
  }
  
  # Return results as a data frame
  return(data.frame(Residue = residue_numbers, Similarity = similarities))
}

# Example protein sequences
protein1 <- "MEPEHSNYDLKNHPTSVQESTLNGSAVTDPQFGNTTNSLPSMNGLLNHENSFASQQSLSSSAFDDSEIVTSTANSTVVSSALSTPKIDDAQNSDDVRVDGTRRSSRAKRPVYRDFSYTEIDEHEIPIPKKRKSKPAPKQKKSVASDDEDAYDKRHRFSINSASGTEIRTSLRSSKGKSVNYNEQEFYDDFEDEEEEVEEQVEEEYEPIIDFVLNHRKRADAQDDDPKSSYQYLIKWQEVSHLHNTWEDYSTLSSVRGYKKVDNYIKQNIIYDREIREDPTTTFEDIEALDIERERKNMLFEEYKIVERIVASETNEEGKTEYFVKWRQLPYDNCTWEDADVIYSMAPNEVYQFLQRENSPYLPYKGVFYNTRPPYRKLEKQPSYIKGGEIRDFQLTGINWMAYLWHRNENGILADEMGLGKTVQTVCFLSYLVHSLKQHGPFLIVVPLSTVPAWQETLANWTPDLNSICYTGNTESRANIREYEFYLSTNSRKLKFNILLTTYEYILKDKQELNNIRWQYLAIDEAHRLKNSESSLYETLSQFRTANRLLITGTPLQNNLKELASLVNFLMPGKFYIRDELNFDQPNAEQERDIRDLQERLQPFILRRLKKDVEKSLPSKSERILRVELSDMQTEWYKNILTKNYRALTGHTDGRGQLSLLNIVVELKKVSNHPYLFPGAAEKWMMGRKMTREDTLRGIIMNSGKMVLLDKLLQRLKHDGHRVLIFSQMVRMLNILGEYMSLRGYNYQRLDGTIPASVRRVSIDHFNAPDSPDFVFLLSTRAGGLGINLNTADTVIIFDSDWNPQADLQAMARAHRIGQKNHVNVYRFLSKDTVEEDILERARRKMILEYAIISLGVTEKSKNSKNDKYDAQELSAILKFGASNMFKATENQKKLENMNLDDILSHAEDRDSSNDVGGASMGGEEFLKQFEVTDYKAEDLNWDDIIPEEEMERIEEEERMLAAQRAKEEERERREEEERENDEDHPSRTYKRTTKSITKRQQRREEMVREKEIRLLYRAMIKFGLVDERFDTIVKEAELQATDPKRIYSLSADMVKACDEAVERLGADDTKNKQPRKAILIEFKGVKNINAETVTLRVKDLTHLHRAYKGLDPLKQIIGYPIRSVHSWNCSWGIKEDSMLLAGINKHGFGCWQAIKNDPDLGLHDKIFLDEAKNDKESRYVPSAVHLVRRGEYLLSVVREHPDLFVVKTDQPTKRKYNRKAPTKSSTRQTTLDGSISNTKKSSRTKKKKEEETNRGDETSPEGTVGEDEVEEEPRQAEPPKRALRSNSGKAASNKRTTRNSMKTHSAMDTLTAVAALDAELDNMSNEKAKEEVDHVKSENGESVNEPNTEDLSLETEENTTVSDISPLVKNEA"
protein2 <- "MSTSAIALALSSSKAIEQLDHVQTETPNLKQEMSESPSNSGVASKRKLQSTEWLDPELYGLRRSGRTRSNPGRYVDTDDQEDVFPSKHRKGTRNGSSFSRHRTIRDLDDEAESVTSEESESDDSSYGGTPKKRSRQKKSNTYVQDEIRFSSRNSKGVNYNEDAYFESFEEEEEEEMYEYATEVSEEPEDTRAIDVVLDHRLIEGHDGSTPSEDYEFLIKWVNFSHLHCTWEPYNNISMIRGSKKVDNHIKQVILLDREIREDPTTTREDIEAMDIEKERKRENYEEYKQVDRIVAKHLNSDGSVEYLVKWKQLLYDFCTWEASSIIEPIAATEIQAFQEREESALSPSRGTNYGNSRPKYRKLEQQPSYITGGELRDFQLTGVNWMAYLWHKNENGILADEMGLGKTVQTVAFLSYLAHSLRQHGPFLVVVPLSTVPAWQETLALWASDMNCISYLGNTTSRQVIRDYEFYVDGTQKIKFNLLLTTYEYVLKDRSVLSNIKWQYMAIDEAHRLKNSESSLYEALSQFKNSNRLLITGTPLQNNIRELAALVDFLMPGKFEIREEINLEAPDEEQEAYIRSLQEHLQPYILRRLKKDVEKSLPSKSERILRVELSDLQMYWYKNILTRNYRVLTQSISSGSQISLLNIVVELKKASNHPYLFDGVEESWMQKINSQGRRDEVLKGLIMNSGKMVLLDKLLSRLRRDGHRVLIFSQMVRMLDILGDYLSLRGYPHQRLDGTVPAAVRRTSIDHFNAPNSPDFVFLLSTRAGGLGINLMTADTVIIFDSDWNPQADLQAMARAHRIGQKNHVMVYRLLSKDTIEEDVLERARRKMILEYAIISLGVTDKQKNSKNDKFSAEELSAILKFGASNMFKAENNQKKLEDMNLDEILEHAEDHDTSNDVGGASMGGEEFLKQFEVTDYKADVSWDDIIPLTEREKFEEEDRLREEEEALKQEIELSSRRGNRPYPSSAVESPSYSGTSERKSKKQMLKDEVLLEKEIRLLYRAMIRYGSLEHRYNDIVKYADLTTQDAHVIKKIAADLVTASRKAVSAAEKDLSNDQSNNKSSRKALLITFKGVKNINAETLVQRLNDLDILYDAMPTSGYSNFQIPMHVRSVHGWSCQWGPREDSMLLSGICKHGFGAWLEIRDDPELKMKDKIFLEDTKQTDNSVPKDKENKEKKVPSAVHLVRRGEYLLSALREHHQNFGIKSSPAISTNGKTQPKKQTANRRQSGKPNVKSAQKIESATRTPSPAISESRKKPSSKDTKIETPSREQSRSQTASPVKSEKDDGNVSLNAEQKARCKELMYPVRKHMKRLRKDSSGLGRAELVKLLTECLTTIGKHIEKTVNDTPSEEKATVRKNLWMFACYFWPKEEVKYTSLISMYEKMK"

# Perform multiple sequence alignment
alignment <- msa(c(protein1, protein2), type = "protein")

# Extract aligned sequences
aligned_seq1 <- alignment@unmasked[1]
aligned_seq2 <- alignment@unmasked[2]

# Define the sliding window size
window_size <- 25

# Calculate sliding window similarity
similarity_data <- sliding_window_similarity(aligned_seq1, aligned_seq2, window_size)

# Plot similarity against residue number
ggplot(similarity_data, aes(x = Residue, y = Similarity)) +
  geom_line() +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(similarity_data$Residue), by = 100)) +
  labs(
    title = "Protein Sequence Similarity (After MSA)",
    x = "Residue Number",
    y = "Similarity (%)"
  ) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))





# Function to assess similarity for specific windows
assess_specific_windows <- function(aligned_seq1, aligned_seq2, windows) {
  # Convert aligned sequences to vectors of characters
  seq1 <- unlist(strsplit(as.character(aligned_seq1), ""))
  seq2 <- unlist(strsplit(as.character(aligned_seq2), ""))
  
  # Initialize vectors for results
  window_labels <- character()
  similarities <- numeric()
  
  for (window in windows) {
    start <- window[1]
    end <- window[2]
    
    # Ensure indices are within bounds of the aligned sequences
    if (start < 1 || end > length(seq1)) {
      cat("Window out of bounds:", start, "-", end, "\n")
      similarities <- c(similarities, NA)
      window_labels <- c(window_labels, paste0("Window_", start, "_to_", end))
      next
    }
    
    # Extract the specific window from each sequence
    window_seq1 <- seq1[start:end]
    window_seq2 <- seq2[start:end]
    
    # Skip windows with gaps ('-')
    if ('-' %in% window_seq1 || '-' %in% window_seq2) {
      similarities <- c(similarities, NA)
    } else {
      similarity <- calculate_similarity(window_seq1, window_seq2)
      similarities <- c(similarities, similarity)
    }
    
    # Store window label
    window_labels <- c(window_labels, paste0("Window_", start, "_to_", end))
  }
  
  # Return results as a data frame
  return(data.frame(Window = window_labels, Similarity = similarities))
}

# Example: Define specific windows of residue numbers
windows <- list(c(1, 50), c(101, 150), c(201, 250))  # Define windows as start-end pairs

# Assess similarity for the specific windows
specific_similarity <- assess_specific_windows(aligned_seq1, aligned_seq2, windows)

# Display the results
print(specific_similarity)