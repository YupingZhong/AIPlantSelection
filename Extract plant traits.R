########################
#Extract flower traits
#To Check the plant names in the 'counts' and integrate the unique plant names extracted from 'counts' and 'interaction' into a single plant list.
data_count_scaled<-readRDS("data_count_scaled_published.rds")#全部互作数据
data_interact<-readRDS("data_interact_published.rds")
merge_plantlist<-data.frame(Plant_species =unique(data_merge$Plant_accepted_name))
merge_plantlist#1246
#write.csv(merge_plantlist,"merge_plantlist20240605.csv")


##delete rows with NA 
null_values <- is.na(merge_plantlist)# 使用 is.na() 函数检查缺失值
has_null <- any(as.vector(null_values))# 通过检查是否有任何缺失值来得知数据框是否包含缺失值
has_null

##flw_shape extract
library(TR8)
class(merge_plantlist)
to_be_downloaded<-c("flw_kugler")
tr8_result<-tr8(species_list = merge_plantlist,download_list = to_be_downloaded,catminat_alternatives=TRUE,  allow_persistent= TRUE)

rows_with_na <- traits_dataframe[which(rowSums(is.na(traits_dataframe)) > 0), ,drop = FALSE]
print(rows_with_na)

print(tr8_result)
#write.csv(tr8_result@results, "merge_plantlist_trait.csv", row.names = TRUE)
#then manual check and revise the plant traits.see "test.merge.trait.csv"
