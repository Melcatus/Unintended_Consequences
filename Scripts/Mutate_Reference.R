#	Script to randomly mutate bases in the given reference sequence. Requires
#	the "Biostrings" package in the BioConductor suite.

#	For UniCon, the reference is a random 1MB region from the Glyma 1.1
#	assembly. The region chosen was Gm07:16717783-17717783. This script will
#	mutate according to estimated nucleotide sequence diversity in soy,
#	approximately 10^-4 per base pair. For each base that is chosen to be
#	mutated, it will mutate according to a 2:1 transition:transversion ratio.

#	Uncomment this when you run this for the first time.
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)

#	Define a function to choose new bases for the mutation. Transitions will be
#	twice as likely as transversions, on average. The weights in 'prob' do not
#	need to sum to 1, and just give relative proabilities.
mutate <- function(base) {
	if(base == "A") {
		newbase <- sample(c("C", "T", "G"), 1, prob=c(1, 1, 2))
	}
	else if(base == "G") {
		newbase <- sample(c("C", "T", "A"), 1, prob=c(1, 1, 2))
	}
	else if(base == "T") {
		newbase <- sample(c("A", "G", "C"), 1, prob=c(1, 1, 2))
	}
	else if(base == "C") {
		newbase <- sample(c("A", "G", "T"), 1, prob=c(1, 1, 2))
	}
	#	If we don't have a normal nucleotide, just kick back an N
	else {
		newbase <- "N"
	}
	return(newbase)
}

#	Set the random seed. Do this with a fixed value so that we can reproduce
#	the results.
set.seed(123)
#	Set the parameters for our mutational process. We need both a numeric and a
#	logical vector. Numeric for subsetting to get the original bases, and
#	logical because the replacement function expects a logical vector.
mutated_bases <- rbinom(1000000, 1, 0.0001)
mutated_log <- mutated_bases == 1
mutated_bases <- as.numeric(which(mutated_bases == 1))

#	Take arguments
args <- commandArgs(TRUE)

#	Read in the reference. readDNAStringSet returns a set of DNAString objects,
#	but we only want the first and only one.
wm82_seq <- readDNAStringSet(args[1])

#	Use our mutate() function to choose new letters to replace in the reference
#	sequence
newbases <- sapply(as.vector(wm82_seq[[1]][mutated_bases]), mutate)
#	And then replace them
mutated_seq <- replaceLetterAt(wm82_seq[[1]], mutated_log, newbases)
#	And write the file out
outfile <- gsub(".fasta", "_Mutated.fasta", args[1])
towrite <- DNAStringSet(mutated_seq)
names(towrite) <- names(wm82_seq)
writeXStringSet(towrite, outfile, format="fasta")

#	We should also save the mutations we created.
subs <- data.frame(
	Positions=mutated_bases,
	Wm82_State=as.vector(wm82_seq[[1]][mutated_bases]),
	Mut_State=as.vector(newbases)
	)
write.table(subs, file="Simulated_Mutations.txt", sep="\t", quote=F, row.names=F)
