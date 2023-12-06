# Probability of Identity

## Brief Description
The probability of identity (PI) is a measure used in forensic genetics to assess the likelihood that two individuals share the same genetic profile at a set of specific genetic markers, typically STRs (short tandem repeats) or SNPs (single nucleotide polymorphisms). The calculation of the probability of identity depends on the number of markers tested and the allele frequencies at each marker in the population.

For independent genetic markers, the probability of identity is calculated as the product of the probabilities of sharing the same alleles at each marker. This assumes that the markers are unlinked, meaning that the inheritance of alleles at one marker is independent of the inheritance of alleles at another marker.

The formula for the probability of identity (PI) for two individuals at a single genetic marker is:

\[ PI = \sum_{i=1}^{k} P_i^2 \]

Where:
- \( PI \) is the probability of identity.
- \( P_i \) is the frequency of the \(i\)-th allele at the marker.
- \( k \) is the number of different alleles at the marker.

For multiple markers, assuming independence, the probabilities are multiplied together:

\[ P_{\text{total}} = \prod_{j=1}^{m} PI_j \]

Where:
- \( P_{\text{total}} \) is the overall probability of identity.
- \( PI_j \) is the probability of identity at the \(j\)-th marker.
- \( m \) is the total number of markers.

The allele frequencies used in these calculations should ideally come from a representative population database. The smaller the probability of identity, the more discriminating power the set of markers provides for distinguishing individuals.

It's important to note that these formulas simplify assumptions and calculations for illustrative purposes, and the actual implementation in forensic genetics can involve more complexities, such as accounting for population substructure and other statistical considerations.

### Probability of Identity with 50 bi-allelic SNPs with a frequency of 50%

To calculate the probability of identity (\(P_{\text{ID}}\)) for 50 independent SNPs with a frequency of 50%:

\[ P_{\text{ID, total}} = \prod_{j=1}^{m} P_{\text{ID}}^j \]

where \( P_{\text{ID}}^j \) is the probability of identity for the \(j\)-th SNP.

Assuming all SNPs are independent, the probability of identity (\( P_{\text{ID}} \)) for a single SNP is given by:

\[ P_{\text{ID}} = \sum_{i=1}^{k} P_i^2 \]

where:
- \( P_i \) is the frequency of the \(i\)-th allele at the SNP.
- \( k \) is the number of different alleles at the SNP.

Given that the frequency of each allele is 50% (0.5), let's assume each SNP has two alleles. Therefore, \( k = 2 \).

\[ P_{\text{ID}} = (0.5)^2 + (0.5)^2 = 0.25 + 0.25 = 0.5 \]

Now, for 50 independent SNPs:

\[ P_{\text{ID, total}} = (0.5)^{50} \]

Calculating this value:

\[ P_{\text{ID, total}} = 1.12589990684 \times 10^{-15} \]

So, the probability of identity for 50 independent SNPs with a frequency of 50% at each SNP is extremely low, indicating a high level of discrimination between individuals in this hypothetical scenario.This is a simplified calculation, and real-world scenarios may involve additional considerations and complexities.

## Definitions from the litterature

### From [Variability of New STR Loci and Kits in US Population Groups](https://france.promega.com/resources/profiles-in-dna/2012/variability-of-new-str-loci-and-kits-in-us-population-groups/)
(John M. Butler, Carolyn R. Hill and Michael D. Coble) 

The probability of identity (PI) is the chance that two unrelated people selected at random will have the same genotype (43) . The PI value of a single locus is determined by summing the square of the observed genotype frequencies. Thus, the lower the PI value, the more variable the genetic marker is in the measured population because there are more genotypes occurring at a lower frequency. Individual locus PI values can be multiplied together with independently inherited loci to create a profile PI (i.e., the product rule). The PI value is a better measure of locus performance than the total number of observed alleles or genotypes because specific alleles may occur in a relatively high frequency and reduce the overall variability, especially if a number of rare alleles occur at this locus. Ideally, we would like to observe a fairly even level of variation across many genotypes at each locus so that there is a greater chance of finding a difference between two unrelated individuals selected at random.

## The math behind Probability of identity

The calculation of the probability of identity (PI) in forensic genetics involves applying mathematical principles related to probability theory and population genetics. Several rules and concepts contribute to these calculations:

1. **Product Rule of Probability:** When events are independent, the probability of two or more independent events occurring together is calculated by multiplying their individual probabilities. In the context of genetic markers (such as SNPs or STRs), if markers are assumed to be independent, the probabilities of sharing specific alleles across multiple markers are multiplied together to determine the overall probability of identity.

2. **Hardy-Weinberg Equilibrium (HWE):** This principle describes the distribution of alleles and genotypes in a population in the absence of evolutionary forces. It provides the expected frequencies of alleles and genotypes based on allele frequencies. HWE helps determine the probabilities of different genotypes at genetic markers in a population.

3. **Allele Frequencies:** The frequencies of different alleles at a genetic marker in a population are fundamental to calculating the probability of identity. These frequencies can be obtained from population databases or estimated within a specific population.

4. **Independence of Genetic Markers:** The assumption of marker independence is crucial in calculating the probability of identity. If markers are unlinked and segregate independently during inheritance, their probabilities can be multiplied together to obtain the overall probability of identity.

5. **Multiplication Rule:** In probability theory, the multiplication rule states that the probability of multiple independent events occurring together is the product of their individual probabilities. In genetic analysis, this rule is applied to multiply the probabilities of sharing specific alleles across different markers.

6. **Binomial Distribution:** When dealing with allele frequencies, particularly in scenarios with two alleles per marker, the binomial distribution might be utilized to calculate the probability of sharing alleles between individuals at a given marker.

By applying these mathematical rules and principles, forensic geneticists and statisticians compute the probability of identity, which serves as a measure of the likelihood that two individuals share the same genetic profile at a set of genetic markers. These calculations are integral in forensic DNA analysis for assessing the uniqueness and discriminatory power of genetic profiles in identifying individuals.