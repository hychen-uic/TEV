
# R function to run plink2, d:plink2 is where plink2.exe resides

### runplink2=function(plink2options="") system(paste("d:/plink2/plink2",plink2options))

# set working directory where the data file resides

### setwd("d:/1000Genome/chr1")

#old approach

### system("d:/plink2/plink2 --vcf chr1.vcf --recode --out chr1")

# new simplified way

### runplink2("--vcf chr1.vcf --recode --out chr1")

# command in R to create a directory
#dir.create("newdirctory")
