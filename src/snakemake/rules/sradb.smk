"""
library(SRAdb)
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
res <- dbGetQuery(sra_con, "select * from sra_ft where run_accession='SRR097786'")

srr_ids <- c(
'SRR2096390',
'SRR3098550',
'SRR1019701',
'SRR1143214',
'SRR568222',
'SRR3098554',
'SRR1603650',
'SRR1509753',
'SRR5789189',
'SRR5789197',
'SRR1143139',
'SRR4250311',
'SRR3098557',
'SRR2774675',
'SRR5789213',
'SRR1603656',
'SRR3098559',
'SRR3231768',
'SRR1522114',
'SRR3098563',
'SRR3098565')

res <- list()
for (srr_id in srr_ids){
    res[srr_id] <- dbGetQuery(sra_con, paste0("select * from sra_ft where run_accession='", srr_id , "'"))
}

# To test, should be faster than res
res2 <- dbGetQuery(sra_con, paste("select * from sra_ft where run_accession in (",paste(sra_ids,collapse=","),")",sep="")); 


res <- dbGetQuery(sra_con, "select * from sra_ft where run_accession='SRR1603656'")


"""

