#' @name MOUSE_Background_signaling_interactome
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Intracellular signaling interactions of mouse genes
#' @description Intracellular signaling interactions of mouse genes
#' @format A data frame with 106175 rows and 3 variables:
#' \describe{
#'   \item{\code{Source}}{character Gene symbol of the regulating gene}
#'   \item{\code{Target}}{character Gene symbol of the regulated gene}
#'   \item{\code{Effect}}{integer Whether the effect is activating - denoted as '1' - or inhibiting - denoted as '-1'} 
#'}
#' @details This dataset is used for creating a Markov Chain model of intracellular signaling
#' to identify significantly active paths from receptors to downstream TFs.
"MOUSE_Background_signaling_interactome"

#' @name MOUSE_dummy.var
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Dummy variable (Mouse)
#' @description Dummy variable (Mouse)
#' @format A character with the content "Dummy"
#' @details This variable is used to represent the niche in the Markov Chain model of intracellular
#' signaling, such that the chain becomes ergodic.
"MOUSE_dummy.var"

#' @name MOUSE_intermediates
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Intermediate signaling molecules (Mouse)
#' @description Intermediate signaling molecules (Mouse)
#' @format A data frame with 2493 rows and 1 variables:
#' \describe{
#'   \item{\code{Gene}}{character Gene symbols of intermediate signaling molecules} 
#'}
#' @details Intermediate signaling molecules are used within the Markov Chain model of 
#' intracellular signaling and constitute all genes that are not at the interface to transcriptional
#' regulation. Importantly, they also contain receptors.
"MOUSE_intermediates"

#' @name MOUSE_Ligands
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Secreted ligands (Mouse)
#' @description Secreted ligands (Mouse)
#' @format A character vector containing secreted ligands
#' @details Secreted mouse ligands have been mapped from secreted human ligands using the 
#' homologene R-package (c.f. ./data-raw/process.R).
"MOUSE_Ligands"

#' @name MOUSE_LR
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Ligand-receptor interactions (Mouse)
#' @description Ligand-receptor interactions (Mouse)
#' @format A data frame with 1422 rows and 3 variables:
#' \describe{
#'   \item{\code{L_R}}{character Ligand receptor interactions in the form "Ligand_Receptor"}
#'   \item{\code{Ligand}}{character Gene symbol of the secreted ligand}
#'   \item{\code{Receptor}}{character Gene symbol of the receptor} 
#'}
#' @details Ligand-receptor interactions have been obtained by obtaining mouse orthologs of
#' human ligand-receptor interactions using the "homologene" R-package (c.f. "./data-raw/process.R").
"MOUSE_LR"

#' @name MOUSE_non_interface_TFs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Non-interface TFs (Mouse)
#' @description Non-interface TFs (Mouse)
#' @format A data frame with 323 rows and 1 variables:
#' \describe{
#'   \item{\code{Gene}}{character Gene symbols of non-interface TFs} 
#'}
#' @details Non-interface TFs are those transcription factors that are not direct targets of
#' signaling interactions.
"MOUSE_non_interface_TFs"

#' @name MOUSE_Receptors
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Receptors (Mouse)
#' @description Receptors (Mouse)
#' @format A character vector containing receptors
#' @details Mouse receptors have been mapped from human receptors using the 
#' homologene R-package (c.f. ./data-raw/process.R).
"MOUSE_Receptors"

#' @name MOUSE_TF_TF_interactions
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title TF-TF-interactions (Mouse)
#' @description TF-TF-interactions (Mouse)
#' @format A data frame with 10239 rows and 3 variables:
#' \describe{
#'   \item{\code{Source}}{character Gene Symbol of the regulating transcription factor }
#'   \item{\code{Target}}{character Gene Symbol of the regulated transcription factor}
#'   \item{\code{Effect}}{integer Whether the effect is activating - denoted as '1' - or inhibiting - denoted as '-1'} 
#'}
#' @details Transcriptional regulatory interactions between TFs. Mouse TFs have been mapped from 
#' human TFs using the homologene R-package (c.f. ./data-raw/process.R).
"MOUSE_TF_TF_interactions"

#' @name MOUSE_tf.db
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Database of mouse TFs
#' @description Database of mouse TFs
#' @format A data frame with 1636 rows and 6 variables:
#' \describe{
#'   \item{\code{Species}}{character The species this TF belongs to. Only contains one value ("Mus_musculus")}
#'   \item{\code{Symbol}}{character The gene symbol of the TF}
#'   \item{\code{Ensembl}}{character The ensembl gene id of the TF}
#'   \item{\code{Family}}{character The protein family the TF belongs to}
#'   \item{\code{Protein}}{character The ensembl protein id of the TF}
#'   \item{\code{Entrez.ID}}{character The entrez gene id of the TF} 
#'}
#' @details Dataset of all mouse TFs catalogued in the Animal TF DB.
#' @source Animal TF DB
"MOUSE_tf.db"

#' @name MOUSE_tfs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Transcription factors (Mouse)
#' @description Transcription factors (Mouse)
#' @format Character vector of length 1636.
#' @details The vector of mouse TFs corresponds to the values in the "Symbol" column of the
#' MOUSE_tf.db dataset.
"MOUSE_tfs"

#' @name HUMAN_Background_signaling_interactome
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Intracellular signaling interactions of human genes
#' @description Intracellular signaling interactions of human genes
#' @format A data frame with 128093 rows and 3 variables:
#' \describe{
#'   \item{\code{Source}}{character Gene symbol of the regulating gene}
#'   \item{\code{Target}}{character Gene symbol of the regulated gene}
#'   \item{\code{Effect}}{integer Whether the effect is activating - denoted as '1' - or inhibiting - denoted as '-1'} 
#'}
#' @details This dataset is used for creating a Markov Chain model of intracellular signaling
#' to identify significantly active paths from receptors to downstream TFs.
"HUMAN_Background_signaling_interactome"

#' @name HUMAN_dummy.var
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Dummy variable (Human)
#' @description Dummy variable (Human)
#' @format A character with the content "DUMMY"
#' @details This variable is used to represent the niche in the Markov Chain model of intracellular
#' signaling, such that the chain becomes ergodic.
"HUMAN_dummy.var"

#' @name HUMAN_intermediates
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Intermediate signaling molecules (Human)
#' @description Intermediate signaling molecules (Human)
#' @format A data frame with 2973 rows and 1 variables:
#' \describe{
#'   \item{\code{Gene}}{character Gene symbols of intermediate signaling molecules} 
#'}
#' @details Intermediate signaling molecules are used within the Markov Chain model of 
#' intracellular signaling and constitute all genes that are not at the interface to transcriptional
#' regulation. Importantly, they also contain receptors.
"HUMAN_intermediates"

#' @name HUMAN_Ligands
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Secreted ligands (Human)
#' @description Secreted ligands (Human)
#' @format A character vector containing secreted ligands
#' @details Secreted ligands have been identified based on their Uniprot annotation
#' (c.f. ./data-raw/process.R). Only ligands being annotated as "Secreted" were considered.
"HUMAN_Ligands"

#' @name HUMAN_Receptors
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Receptors (Human)
#' @description Receptors (Human)
#' @format A character vector containing human receptors
#' @details Human receptors
"HUMAN_Receptors"

#' @name HUMAN_tfs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Transcription factors (Human)
#' @description Transcription factors (Human)
#' @format Character vector of length 1665.
#' @details The vector of human TFs corresponds to the values in the "Symbol" column of the
#' HUMAN_tf.db dataset
"HUMAN_tfs"

#' @name HUMAN_LR
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Ligand-receptor interactions (Human)
#' @description Ligand-receptor interactions (Human)
#' @format A data frame with 1450 rows and 3 variables:
#' \describe{
#'   \item{\code{L_R}}{character Ligand receptor interactions in the form "Ligand_Receptor"}
#'   \item{\code{Ligand}}{character Gene symbol of the secreted ligand}
#'   \item{\code{Receptor}}{character Gene symbol of the receptor} 
#'}
"HUMAN_LR"

#' @name HUMAN_non_interface_TFs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Non-interface TFs (Human)
#' @description Non-interface TFs (Human)
#' @format A data frame with 443 rows and 1 variables:
#' \describe{
#'   \item{\code{Gene}}{character Gene symbols of non-interface TFs} 
#'}
#' @details Non-interface TFs are those transcription factors that are not direct targets of
#' signaling interactions.
"HUMAN_non_interface_TFs"

#' @name HUMAN_TF_TF_interactions
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title TF-TF-interactions (Human)
#' @description TF-TF-interactions (Human)
#' @format A data frame with 12203 rows and 3 variables:
#' \describe{
#'   \item{\code{Source}}{character Gene Symbol of the regulating transcription factor}
#'   \item{\code{Target}}{character Gene Symbol of the regulated transcription factor}
#'   \item{\code{Effect}}{integer Whether the effect is activating - denoted as '1' - or inhibiting - denoted as '-1'} 
#'}
"HUMAN_TF_TF_interactions"

#' @name HUMAN_tf.db
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
#' @title Database of human TFs
#' @description Database of human TFs
#' @format A data frame with 1665 rows and 6 variables:
#' \describe{
#'   \item{\code{Species}}{character The species this TF belongs to. Only contains one value ("Homo_sapiens")}
#'   \item{\code{Symbol}}{character The gene symbol of the TF}
#'   \item{\code{Ensembl}}{character The ensembl gene id of the TF}
#'   \item{\code{Family}}{character The protein family the TF belongs to}
#'   \item{\code{Protein}}{character The ensembl protein id of the TF}
#'   \item{\code{Entrez.ID}}{character The entrez gene id of the TF} 
#'}
#' @details Dataset of all human TFs catalogued in the Animal TF DB.
#' @source Animal TF DB
"HUMAN_tf.db"

