ann <- read.table("~/Meta_Analysis/data/original/LUNGevity_allsamples_curated.txt", sep="\t", header=TRUE)
# lnames <- load("~/Meta_Analysis/kyoto/aracne/data/lungevity_ann_exp_82s_14308g_no_dysgr.rda")

names(ann)[names(ann)=="AGE_at_brush"] <- "age"
names(ann)[names(ann)=="Sex"] <- "sex"
names(ann)[names(ann)=="PACK_YEARS"] <- "pack_years"
names(ann)[names(ann)=="Race"] <- "race"
names(ann)[names(ann)=="SUBGROUP"] <- "smoking_status"
names(ann)[names(ann)=="COPDdef1"] <- "copd_status1"
names(ann)[names(ann)=="COPDdef2"] <- "copd_status2"
names(ann)[names(ann)=="DYSGR_at_Brush"] <- "DYSGR"
names(ann)[names(ann)=="Batch"] <- "batch"
names(ann)[names(ann)=="Seq_Name"] <- "filename"

# Fix smoking status
ann$smoking_status <- as.factor(ifelse(ann$smoking_status=="CURRENT-SMOKER", "current", ifelse(ann$smoking_status=="EX-SMOKER", "former", ann$smoking_status)))

# Fix race
ann$race <- as.factor(ifelse(ann$race=="WHITE", "Caucasian", ann$race))

# Fix sex
ann$sex <- as.factor(ifelse(ann$sex=="F", "female", ifelse(ann$sex=="M", "male", ann$race)))

# Fix COPD
ann$copd_status1 <- as.factor(ifelse(ann$copd_status1==1, "COPD", ifelse(ann$copd_status1==0, "no COPD", ann$copd_status1)))
ann$copd_status2 <- as.factor(ifelse(ann$copd_status2==1, "COPD", ifelse(ann$copd_status2==0, "no COPD", ann$copd_status2)))

# calculate COPD status based on the FEV1/FVC ratio and FEV1pp (per Katie's definition)
ann$copd_status3 <- as.factor(sapply(1:nrow(ann), function(i) as.factor(c("no COPD", "COPD")[as.numeric(ann$'FEV1_FVC.'[i]/100 < 0.70 & ann$'FEV1.Pred'[i]/100 < 0.80)+1L])))

# compare definitions 1 and 2 to 3
table(ann$copd_status1==ann$copd_status3)  # FALSE = 8
table(ann$copd_status2==ann$copd_status3)  # all TRUE; using copd_status2

names(ann)[names(ann)=="copd_status3"] <- "copd_status"  # copd_status 2 and 3 are equivalent, so keeping 1 and 2 as ref 
                                                         # and changing 3 to absolute

# Add cancer status
ann$cancer_status <- "no cancer"
ann$cancer_status <- as.factor(ann$cancer_status)

# Fix Dysplasia grade and add category from Lam
ann$DYSGR_CAT <- dict$DYSGR_CAT[ match(ann$DYSGR, dict$DYSGR)]  # use grade-category dictionary from Lam
ann[ is.na(ann$DYSGR_CAT), c("DYSGR", "DYSGR_CAT")]
ann$DYSGR_CAT[ ann$DYSGR==5.45] <- "SevD"
ann$DYSGR_CAT[ ann$DYSGR==4.10] <- "MildD"
ann$DYSGR_CAT[ ann$DYSGR==4.15] <- "MildD"
ann$DYSGR_CAT <- droplevels(ann$DYSGR_CAT)

# add column with "accepted_hits" (filenames)
ann$sample_id <- ann$filename
ann$filename <- paste0(ann$sample_id, "_accepted_hits")
ann$filename <- gsub("-", "\\.", ann$filename)

# keep only variables of interest
ann.lgv <- ann[, c("filename", "batch", "sex", "age", "race", "cancer_status", "copd_status", "smoking_status", "pack_years", "DYSGR", "DYSGR_CAT")]

lgv.ann.orig <- lgv.ann
lgv.ann <- ann.lgv[ann.lgv$filename %in% colnames(lgv.exp), ]
lgv.ann <- lgv.ann[ match(colnames(lgv.exp), lgv.ann$filename), ]

save(lgv.ann, lgv.exp, file="~/Meta_Analysis/kyoto/aracne/data/lungevity_ann_exp_82s_14308g.rda")









