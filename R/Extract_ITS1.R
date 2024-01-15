# ------------------------------------------------------------------------------------------------#
# Extract ITS1 region for a given SRA project
# Apply this script to each SRA Accession directory on FWD reads only
# Author: Geoffrey Zahn
# Requirements: itsxpress (available in miniconda3)
# ------------------------------------------------------------------------------------------------#


# get directories
datapath <- "./data"
project_directories <- file.path(datapath,list.files(datapath))

getwd()

for(i in project_directories){
# Run raw sequences through itsxpress to extract fungal ITS1 regions (This takes several hours for each project!)
system(paste0("cd ",i,"/seqs/;",
              "for i in *.fastq.gz;",
              "do /home/gzahn/miniconda3/bin/itsxpress ", # needs absolute file path!?
              "--fastq $i --single_end --outfile $i.ITS1.gz --region ITS1 --taxa Fungi --cluster_id 1 --threads 16;", 
              "done"),
       wait = TRUE,intern = FALSE)
}

