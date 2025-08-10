# source(ML_RF_model.R)
	rf_model <- randomForest(train_data, train_target, ntree = 50, nodesize = 3)
	predictions <- predict(rf_model, test_data)	 	  
	writeDTable(paste0("/Users/ecomm77/Dropbox/Data/ANN_data/ANN_et_",region,"_RF_",dynamic_option,".dtable"),data.frame(predictions,test_target))
	if(dynamic_option=="Y")
	{
	INDEX = c((365*0+5):(365*1+5))
	write.csv(importance(rf_model)[INDEX], file = paste0(data_path,"ANN_data/",region,"_",ML_name,"_",year,"_",dynamic_option,"_importance_RF_delt.csv"))
	}
	
	if(dynamic_option=="N")
	{
	write.csv(importance(rf_model), file = paste0(data_path,"ANN_data/",region,"_",ML_name,"_",year,"_",dynamic_option,"_importance_RF.csv"))
	}

	if(dynamic_option=="single")
	{
	write.csv(importance(rf_model), file = paste0(data_path,"ANN_data/",region,"_",ML_name,"_",year,"_",dynamic_option,"_importance_RF_single.csv"))
	}

		save(rf_model, file=paste0("/Users/ecomm77/Dropbox/Data/ANN_data/",region,"_",ML_name,"_",year,"_",dynamic_option,alpha,".Rdata"))