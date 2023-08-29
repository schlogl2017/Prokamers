################################################################
## Explore the probabilities of read mapping
##
## Author: Jacques van Helden
## Date: 2013-11-18
## Context: ecole de bioinformatique AVIESAN, Roscoff, France
##
## We position a 35bp read over a random position of the genome. 
## For the sake of simplicity, we assume that nucleotides are equiprobable and independent (beware: this assumption is an over-simplification). 
##
## What are the respective probabilities to obtain the following results ?
## a. A perfect match over the whole length of the read.
## b. Not a single matching residue over the whole length of the read.
## c. A match over the first 30 residues (irrespecfive of what follows)
## d. A match over the first 30 residues, followed by a mismatch.
## e. A match over the first 30 residues, followed by 5 mismatches.
## f. Exactly 30 matches and 5 mismatches, irrespective of the matching/mismatching positions.
## g. At least 30 matching residues.

L <- 35 ## Length of the read
p <- 0.25 ## matching probability for a single residue, assuming equiprobable and independent nucleotides
genome.size <- 3e9
n <- 30 ## Number of matches
read.nb <- 50e6

proba <- vector()

## a. A perfect match over the whole length of the read.
proba["fullmatch"] <- p^L

## b. Not a single matching residue over the whole length of the read.
q <- 1 -p ## mismatch probability 
proba["fullmismatch <- q^L

## c. A match over the first n residues
proba["first.n.matches <- p^n

## d. A match over the first n residues, followed by a mismatch.
proba["first.n.then.one.mis"] <-  p^n*q

## e. A match over the first n residues, followed by 5 mismatches.
proba["first.n.then.all.mis"] <-  p^n*q^(L-n)

## f. Exactly n matches and 5 mismatches, irrespective of the matching/mismatching positions.
proba["n.matches.anywhere"] <-  choose(n=L, k=n)*p^n*q^(L-n)

## The same can be computed with the binomial density
proba["dbinom"] <-  dbinom(x=n, size=L, p=p)

## g. At least n matching residues.
i <- n:L
proba["ge.n.matches"] <-  sum(choose(n=L, k=i)*p^i*q^(L-i))

## The same can be computed with the binomial pvalue
##
## Beware: the upper tail is exclusive : it computes the proba of
## having more than q matches. However, we want the proba to have at
## least n matches. This equals the proba to have more than (q=n-1)
## matches.
proba["pbinom"] <-  pbinom(q=n-1, size=L, p=p, lower=FALSE)


## Cast the proba in a data frame
result <- data.frame(proba)


## Compute e-values for the different matching conditions
result$eval.genome <- result$proba * 2 * genome.size
result$eval.genome.allreads <- result$eval.genome * read.nb

  
