# upload data to synapse

# install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))

library(synapser)

synLogin(email="attilagabor", apiKey="/YwdwbcFPKfTSazwKChjXSXd/ZZ8BQ0DXOkh/JuGWvmGtnZCfUfRuZ5Ixid5RjGdvH0J8QlUImDtlBpsC4uRIA==",rememberMe=TRUE)

# live site
live_challenge_data_folder = "syn20564743"

# Median Phospho data ----------------------------------------------------------
median_data_folder = Folder("median_phospho", parentId=live_challenge_data_folder)
stored_median_data_folder = synStore(median_data_folder)

median_data_ent = File("challenge_data/median_phospho/median_phospho_data.csv",
					   parentId=stored_median_data_folder$properties$id)
synStore(median_data_ent)


# Template predictions ---------------------------------------------------------
template_folder_ent = Folder("prediction_template", parentId=live_challenge_data_folder)
stored_template_folder = synStore(template_folder_ent)

for (template_files in list.files("./challenge_data/prediction_templates",full.names = T)){
	
	template_data_ent = File(template_files, parentId=stored_template_folder$properties$id)
	synStore(template_data_ent)
}


# single cell data -------------------------------------------------------------

sc_base_folder_ent = Folder("single_cell_phospo", parentId=live_challenge_data_folder)
stored_sc_base_folder_ent = synStore(sc_base_folder_ent)

# complete ---
complete_cell_lines_folder_ent = Folder("complete_cell_lines", parentId=stored_sc_base_folder_ent$properties$id)
stored_complete_cell_lines_folder_ent = synStore(complete_cell_lines_folder_ent)


for (complete_files in list.files("./challenge_data/single_cell_phospho/complete_cell_lines/",full.names = T)){
	
	template_data_ent = File(complete_files, parentId=stored_complete_cell_lines_folder_ent$properties$id)
	synStore(template_data_ent)
}

# subchallenge 1,2,4 ----
for(scfolders in c("subchallenge_1","subchallenge_2","subchallenge_4") ){
	sc_folder_ent = Folder(scfolders, parentId=stored_sc_base_folder_ent$properties$id)
	stored_sc_folder_ent = synStore(sc_folder_ent)
	
	for (sc_files in list.files(paste0("./challenge_data/single_cell_phospho/",scfolders,"/"),full.names = T)){
		
		template_data_ent = File(sc_files, parentId=stored_sc_folder_ent$properties$id)
		synStore(template_data_ent)
	}
	
	
}


# genomics data -------------------------------------------------------------
base_folder_ent = Folder("transcriptomics_genomics", parentId=live_challenge_data_folder)
stored_base_folder_ent = synStore(base_folder_ent)

for (template_files in list.files("./challenge_data/transcriptomics_genomics/",full.names = T)){
	
	template_data_ent = File(template_files, parentId=stored_base_folder_ent$properties$id)
	synStore(template_data_ent)
}


# Supplementary files ----------------------------------------------------------

file1 = File("./challenge_data/Antibody_table.csv", parentId=live_challenge_data_folder)
synStore(file1)
file2 = File("./challenge_data/CellLines.csv", parentId=live_challenge_data_folder)
synStore(file2)
file3 = File("./challenge_data/FileID_table.csv", parentId=live_challenge_data_folder)
synStore(file3)



### Validation data ------------------------------------------------------------

# This part is only run manually, we dont want to expose the validation data by accident
if(FALSE){
	
	# Staging site
	staging_challenge_project = "syn20366916"
	goldstandard_folder_id = "syn20630281"
	
	
	for (template_files in list.files("./challenge_data/validation_data",full.names = T)){
		
		template_data_ent = File(template_files, parentId=goldstandard_folder_id)
		synStore(template_data_ent)
	}
	
	
	
}


# Create folder in staging
folder_ent = Folder("median_phospho", parentId=staging_challenge_project)
stored_folder_ent = synStore(folder_ent)

