# upload some test cases to synapse.

library(synapser)

synLogin(email="attilagabor")

staging_challenge_project = "syn20366916"

base_folder_ent = Folder("test_submission", parentId=staging_challenge_project)
stored_base_folder_ent = synStore(base_folder_ent)


for (template_files in list.files("./challenge_data/test_scoring/test_data/",full.names = T)){
	
	template_data_ent = File(template_files, parentId=stored_base_folder_ent$properties$id)
	synStore(template_data_ent)
}





### Uploading false predictions


synLogin(email="attilagabor", apiKey="/YwdwbcFPKfTSazwKChjXSXd/ZZ8BQ0DXOkh/JuGWvmGtnZCfUfRuZ5Ixid5RjGdvH0J8QlUImDtlBpsC4uRIA==",rememberMe=TRUE)

staging_challenge_project = "syn20366916"

base_folder_ent = Folder("test_submission", parentId=staging_challenge_project)
stored_base_folder_ent = synStore(base_folder_ent)


for (template_files in list.files("./challenge_data/test_scoring/test_data",pattern = "sc4_test_data_v[6-8]+",full.names = T)){
	
	template_data_ent = File(template_files, parentId=stored_base_folder_ent$properties$id)
	synStore(template_data_ent)
}
