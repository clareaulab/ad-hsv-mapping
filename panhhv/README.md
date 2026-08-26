# Pan Human Herpesvirus + hg38 Reference 

Following code reproducibily creates a combined viral + human reference. 
Additionally has commands to build a cellranger index, kallisto index, and STAR index. 

# Workflow 

There are 3 main steps to reproduce this reference. 
The first step is has the completed files in this repo. 

However due to size limitations of github the hg38 containing files are not present. 
To recreate all indexes you can download the cellranger hg38 reference in the `hg38_reference.json` and contiue from step 06 onward. 

## 1) Compile Viral References
Run steps 01-05 to create the hardmasked human herpes virus seqeunces. 

## 2) Create Combined Reference 
Once the hg38 files are downloaded, modify step 06 to point to them and create the viral + human reference fasta + gtf . 

## 3) Create Indexes 
With the combined genome files you can now run steps 07-09 to create the indexes for downstream quantification. 
