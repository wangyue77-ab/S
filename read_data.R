## https://portal.gdc.cancer.gov/repository
dir.create("data_in_one")

for (dirname in dir("rawdata/")){  
  file <- list.files(paste0(getwd(),"/rawdata/",dirname),pattern = "*.counts") 
  file.copy(paste0(getwd(),"/rawdata/",dirname,"/",file),"data_in_one")  
}


dir.create("rawcounts")


metadata <- jsonlite::fromJSON("metadata.cart.json")
require(dplyr)
metadata_id <- metadata %>% 
  dplyr::select(c(file_name,associated_entities)) 

metadata_id$associated_entities[[1]]$entity_submitter_id


naid_df <- data.frame()
for (i in 1:nrow(metadata)){
  naid_df[i,1] <- substr(metadata_id$file_name[i],1,nchar(metadata_id$file_name[i])-3)
  naid_df[i,2] <- metadata_id$associated_entities[i][[1]]$entity_submitter_id
}

colnames(naid_df) <- c("filename","TCGA_id")


files <- dir("rawcounts")

myfread <- function(files){
  data.table::fread(paste0("rawcounts/",files))[,2]
}


f <- lapply(files,myfread)
f <- do.call(cbind,f)


rownames(naid_df) <- naid_df[,1]
naid_df <- naid_df[files,]


colnames(f) <- naid_df$TCGA_id

test <- f[1:100,1:4]

gene_id <- data.table::fread(paste0("rawcounts/",files[1]))$V1

expr_df <- cbind(gene_id=gene_id,f)
test <- expr_df[1:100,1:4]

tail(expr_df$gene_id,10)

expr_df <- expr_df[1:(length(expr_df$gene_id)-5),]

expr_df_LUAD <- expr_df
naid_df_LUAD <- naid_df
save(expr_df_LUAD,naid_df_LUAD,file = "LUAD_exprdf.Rdata")

