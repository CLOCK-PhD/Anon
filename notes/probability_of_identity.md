### Probability of Identity

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