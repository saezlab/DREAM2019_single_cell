# prepare data to test the scoring function



# AIM 1.1 ----------------------------------------------------------------------

# prepare random data for testing and save to the test folder
prediction_aim11 = read_csv("./challenge_data/prediction_templates/AIM_11_template_data.csv")
reporters <- c('p.Akt.Ser473.', 'p.ERK','p.S6', 'p.HER2', 'p.PLCg2') 
prediction_aim11 <- prediction_aim11 %>% mutate_at(reporters,rnorm)
prediction_aim11 %>% write_csv("./challenge_data/test_scoring/test_data/aim_11_test_data.csv")


## continue of needed. 

# AIM 1.2.1     ----------------------------------------------------------------

# prepare random data for testing and save to the test folder
prediction_aim121 = read_csv("./challenge_data/prediction_templates/AIM_121_template_data.csv")
reporters  <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
								  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
								  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
								  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
								  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
								  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
								  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
								  'p.STAT3', 'p.STAT5') 
prediction_aim121 <- prediction_aim121 %>% mutate_at(reporters,rnorm)
prediction_aim121 %>% write_csv("./challenge_data/test_scoring/test_data/aim_121_test_data.csv")





# AIM 1.2.2     ----------------------------------------------------------------
prediction_aim122 = read_csv("./challenge_data/prediction_templates/AIM_122_template_data.csv")
reporters  <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
				'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
				'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
				'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
				'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
				'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
				'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
				'p.STAT3', 'p.STAT5') 
prediction_aim122 <- prediction_aim122 %>% mutate_at(reporters,rnorm)
prediction_aim122 %>% write_csv("./challenge_data/test_scoring/test_data/aim_122_test_data.csv")





# AIM 2    ---------------------------------------------------------------------

prediction_aim2 = read_csv("./challenge_data/prediction_templates/AIM_2_template_data.csv")
reporters  <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
				'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
				'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
				'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
				'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
				'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
				'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
				'p.STAT3', 'p.STAT5') 
prediction_aim2 <- prediction_aim2 %>% mutate_at(reporters,rnorm)
prediction_aim2 %>% write_csv("./challenge_data/test_scoring/test_data/aim_2_test_data.csv")



