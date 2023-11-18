# Overview of anonymity metrics and what can be of interest in our study


## Speacial features of genomics data
- Health / Behavior: DNA contains information about health and behavior
- Phenotypes: Genetic variation is associated with various traits (e.g., eye color)
- Static (Traceable): DNA does not change much over time
- Uniqueness:  DNA of two individuals can be dinstingushed from one another
- Value: Importance of information content in DNA, will increase in time with new discoveries
- Kinship: DNA contains information about blood relatives

## Uses of genomic data
- Health care
- Research
- Legal and Forensics
- Direct-to-Consumer

## Relevance of genome privacy

Leakage of information may have serious implications if misused

- Genetic discrimination (insurance, employment, eduction)
- Blackmail

## Known privacy risks

### Re-identification threats

- Recovering the identities of the individual
- Pseudo-anonymisation: removal of explicit and quasi-identifying attributes
- Genomic data cannot be anonymized by just removing the indentifying information
- There is always a risk for an adversary to infer the phenotype of a DNA-material donnor
- Re-identification can be achieved through inspecting the background information that comes with publicized DNA sequences (public genealogy databases; Personal Genome Project; correlation between Y chromosome and last name)

### Phenotype inference
- ~75 independant SNPs are enough to uniquely distinguish one individual from others
- Correlation of genomic data: *Partially available genomic data* can be used to infer the unpublished genome data due to linkage disequilibrium (LD), a correlation between regions of the genome
- Kin privacy breach: disclosure of relatives can threaten the privacy. Reconstruction attacks can be carryied using genomic data of a subset of family members dans publicly known genomic background information (LD, MAF)

### Other threats
- Anonymous paternity breach
- Legal and forensic

## SNP as biomarkers
- 1st source of polymorphism in the genome
- 1/1000nt
- SNP = allel with a frequency >1% in the population
- Bi-allelic SNPs are the most common
- ~90% of genome variations
- Low mutation rate ($10^{-8}$)

## Solutions

Despite the risks associated wtih genomic data, we can find ways to mitigate them to move forward.

### Health care
- Personalized medicine: homomorphic encryption
- Raw aligned genomic data: too big for encryption
- Pseudo-anonymization

### Research
- GWAS (analyzing the statistical correlation between the variants of a case groupe and a control group) : Applying noise in the data
- Sequence comparison: encryption, secure multiparty computation
- Person-level genome sequence record: masking SNPs
- Sequence alignment: differential privacy (adding controlled noise or distortion to genetic data to protect against re-identification attacks)

## Challenges for genome privacy

Both a lack of and excess of privacy have the potential to derail the expected benefits of genomics in health care and research. 

The efficient and secure handling of individual genotype and sequence data will be central to the implementation of genomic medicine.

An excess of privacy-related hurdles could slow down research and interfere with large-scale adoption of genomics in clinical practice.

- Most end-users are not familiar with computer science and are almost exclusively interested in the clinical utility of test results
- Bioinformaticians are trained to maximize the information to be extracted from genomic data
- Their curricula do not address security and privacy concerns
- Computer scientists are, but lack knowledge in biology and genomics
- Strong potential of cross-fertilization between the two disciplines




## Techniques d'anonymisation
- **Data encryption**: Data encryption involves using mathematical algorithms to convert sensitive genetic data into a format that can only be decoded with a decryption key. This technique can help protect data from unauthorized access or use.
- **Data masking**: Data masking involves removing or obscuring sensitive information from genetic data, such as identifying information or other data that could be used to re-identify an individual. This technique can help protect data while still allowing for useful analysis.
- **Data perturbation**: Data perturbation involves adding random noise or variation to genetic data to make it more difficult to identify specific individuals. This technique can help protect against inference attacks and re-identification attacks.
- **Differential privacy**: Differential privacy is a technique that involves adding controlled noise or distortion to genetic data to protect against re-identification attacks. This technique can help protect against both internal and external threats to privacy.
- **Data aggregation**: Data aggregation involves combining genetic data from multiple individuals to create a summary statistic or aggregate measure. This technique can help protect privacy by making it more difficult to identify specific individuals from their genetic data.
- **Access control**: Access control involves limiting access to sensitive genetic data to authorized individuals or organizations. This technique can help ensure that genetic data is only used for approved purposes and that data is not accessed by unauthorized individuals.

## Definition

Anonymity : the state of being not identifiable within a set of subjects (Pfitzmann & Hansen, 2000)

## Anonymity metrics

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

### Information theory

1. **Shannon's Entropy:**
   - **Definition:** Shannon's entropy is a measure of the uncertainty or information content associated with a random variable. In the context of information theory, it quantifies the average amount of surprise or unpredictability in a set of possible outcomes.
   - **Formula:** For a discrete random variable X with probability mass function P(X), the entropy H(X) is calculated as: \(H(X) = - \sum_{i} P(x_i) \log_2(P(x_i))\).
   - **Interpretation:** Higher entropy values indicate higher uncertainty or randomness in the variable, while lower entropy values suggest more predictability.

2. **Effective Anonymity Set Size:**
   - **Definition:** In the context of privacy and anonymity, the effective anonymity set size refers to the number of individuals or entities that an observer or adversary could reasonably assume as the source of a particular piece of information.
   - **Significance:** A larger effective anonymity set size makes it more challenging for an adversary to identify the specific source of information, thus enhancing privacy.

3. **Degree of Anonymity:**
   - **Definition:** The degree of anonymity is a measure of how well an individual in a group can be protected from identification. It is often associated with the concept of anonymity sets.
   - **Anonymity Sets:** Anonymity sets represent the group of individuals or entities that could potentially be the source of certain information, making it difficult to pinpoint a specific source.
   - **Degree Calculation:** The degree of anonymity is influenced by factors such as the size of the anonymity set, the distribution of individuals within the set, and the observer's ability to differentiate between members of the set.
   - **Higher Degree:** A higher degree of anonymity implies greater protection, as it becomes more challenging to distinguish a specific member within the anonymity set.

In summary, Shannon's entropy measures information uncertainty, effective anonymity set size relates to privacy by assessing the potential sources of information, and the degree of anonymity reflects how well individuals within a group are protected from identification. These concepts are fundamental in understanding and analyzing information and privacy in various contexts.

### Metrics from "Genome privacy metrics: a systematic comparison" (Wagner, 2015)
#### Excluded
- **Differential privacy** : offers privacy guarantees for database queries; could be used to prevent the adversary from acquiring a probability distribution in the first place
- **k-anonymity** [Malin 2005] states that an individual cannot be distinguished among at least k − 1 other individuals.
- **Genomic privacy metric** : assumes that the adversary only aims to infer whether an individual’s genome has a specific SNP or not

#### Included

##### Metrics measuring the adversary's error
- **Expected estimation error** : quantifies the adversary’s correctness by computing the expected distance between the adversary’s estimate and the true value for every SNP. [Humbert et al. 2013]
- **mean squared error** : computed as the squared difference between the true value and the adversary’s estimate, averaged over all SNPs [Oya et al. 2014]
- **Percentage incorrectly classified** measures how often the highest probability in the adversary’s estimate does not correspond to true SNP value [Narayanan and Shmatikov 2009]

##### Metrics Measuring the Adversary’s Uncertainty.
- **Entropy** quantifies the amount of information contained in a random variable. Used as a privacy metric, it indicates the adversary’s uncertainty [Serjantov and Danezis 2002]; Entropy can be normalized to a range of [0, 1] by dividing it by Hartley entropy, that is, the logarithm of the number of outcomes [Humbert et al. 2013]
- **Hartley entropy**, or max-entropy : optimistic metric because it only accounts for the number of outcomes, but not for additional information the adversary may have. In the context of genomics, however, the number of outcomes per SNP is known to be 3, and therefore max-entropy is not useful and has been excluded from the evaluation. [Clauß and Schiffner 2006]
- **Min-entropy** : a pessimistic metric because it is based only on the probability of the most likely outcome, regardless of whether this is also the true outcome [Clauß and Schiffner 2006]. Min-entropy is a conservative measure of how certain the adversary is of his estimate.
- **Cumulative** entropy is based on the notion that the adversary’s uncertainty increases when privacy protection is applied at several independent points. Cumulative entropy is computed as the sum of individual entropies [Freudiger et al. 2007]. In the context of genomics, we sum over the entropies computed for each SNP.
- **Conditional entropy**, or the entropy of X conditioned on Y , measures the amount of information needed to fully describe X, provided that Y is known [Diaz et al. 2007]. For genomic privacy, X can be chosen as the true SNP value and Y as the adversary’s estimate. This measures how much more information the adversary needs to find the true value.
- **Inherent privacy** [Agrawal and Aggarwal 2001; Andersson and Lundin 2008] and **conditional privacy** [Andersson and Lundin 2008] are derivations of base metrics (entropy and
conditional entropy, respectively), each computed as 2base metric. While the base metrics are interpreted as bits of information, these metrics can be interpreted as the number of binary questions an adversary has to ask to resolve his uncertainty. Asymmetric entropy can also be used as a per-SNP metric to measure privacy for individual SNPs.

##### Metrics Measuring Information Gain/Loss
- **Amount of leaked information** [Wang et al. 2009; Ayday et al. 2014] counts the number of leaked SNPs. A SNP is considered leaked when the adversary’s estimate for the true outcome is above the threshold α. A threshold of 1 means that a SNP is considered leaked only if the adversary is absolutely certain.
- **Information surprisal**, or self-information, quantifies how much information is contained in a specific outcome of a random variable [Chen et al. 2013]. In the context of genomics, the outcome is the true value of a SNP, and the information content is the probability the adversary assigns to this outcome. Informally, information surprisal quantifies how surprised the adversary would be upon learning the true value of a SNP.
- **Mutual information** measures how much information is shared between two random variables X and Y [Lin et al. 2002]. As before, X can be chosen as the true SNP value and Y as the adversary’s estimate.
- **Conditional privacy loss** [Andersson and Lundin 2008] is derived from mutual information. While mutual information is interpreted as the bits of information shared between the true value and the adversary’s estimate, conditional privacy loss can be interpreted as the number of binary questions an adversary has to ask to arrive at the true value.
- **Relative entropy**, or Kullback-Leibler divergence, between two random variables Y and X measures the information that is lost when X is used to approximate Y [Deng et al. 2007]. In the context of genomics, good choices for Y and X are the true value and the adversary’s estimate, respectively. This measures how many additional bits of information the adversary needs to reconstruct the true value.
- **Variation of information** is derived from mutual information so that it fulfills the conditions for a distance metric in the mathematical sense, especially the triangle inequality [Meil˘a 2007]. It describes the distance between two random variables, chosen as the true value and the adversary’s estimate.

##### Metrics Measuring the Adversary’s Success Probability.
- **adversary’s success rate** captures how likely it is for the adversary to succeed. In the context of genomics, we can define success on a per-SNP basis as the probability of correctly inferring a SNP value, and aggregate to a per-individual metric by computing the average probability for all SNPs [Ayday et al. 2013].
- **User-specified innocence** can be seen as a counterpart to the amount of leaked information, because it counts the number of SNPs that remain private [Chen and Pang 2012]. A SNP is considered private if the adversary’s estimate for the true outcome is below the threshold α. A threshold of 0 means that a SNP is considered private only if the adversary considers it impossible. Many scenarios will therefore adopt a higher threshold.

##### Metrics Measuring Similarity/Diversity
- **Coefficient of determination r2**: describes how well a statistical model approximates data. It is typically used for linear regression where a value of 1 indicates a perfect fit [Kalogridis et al. 2010]. In the context of genomics, the adversary’s estimate can be used as statistical model, and the true SNP values represent the data.

##### Other metrics
- **Health privacy** focuses on those SNPs known to contribute to a specific disease. Health privacy uses a base metric to compute per-SNP values, and then aggregates to a per-individual metric using a weighted and normalized sum [Humbert et al. 2013]. The weights ci should be chosen to reflect how much each SNP contributes to the disease.

##### Other information
"*the absolute value of genomic privacy depends on the number of SNPs and the choice of severities.*"

## Forensic

### Linkage disequilibrium

Linkage disequilibrium (LD) is often discussed in the context of single nucleotide polymorphisms (SNPs), which are variations at a single position in a DNA sequence among individuals. SNPs are widely used as genetic markers in association studies and genomic research, and understanding the patterns of LD among SNPs is crucial for these applications.

In the context of SNPs, LD refers to the non-random association of alleles at different SNP loci. When two SNPs are in strong LD, it means that certain allelic combinations at these two loci are inherited together more frequently than expected by chance. Conversely, if there is little or no LD between two SNPs, the alleles at one locus are inherited independently of the alleles at the other locus.

Here are some key points regarding LD in SNPs:

1. **Haplotype Blocks:**
   - LD can result in the formation of haplotype blocks, which are regions of the genome where a set of SNPs tends to be inherited together. These blocks are bounded by recombination hotspots, where genetic recombination occurs more frequently.

2. **Association Studies:**
   - LD is a critical factor in association studies that aim to identify genetic variants associated with traits or diseases. If a SNP is in strong LD with a causal variant, it can serve as a proxy for genotyping purposes. This is known as tagging or using a tag SNP.

3. **LD Decay:**
   - LD between SNPs typically decays with physical distance along the chromosome. The rate of LD decay depends on factors such as recombination rates, population history, and selection. In some cases, SNPs that are physically close may exhibit strong LD, while more distant SNPs may show weaker or no LD.

4. **Genome-wide Association Studies (GWAS):**
   - In GWAS, researchers examine thousands to millions of SNPs across the genome to identify associations with diseases or traits. Understanding LD patterns is essential for selecting informative SNPs and interpreting the results.

5. **Imputation:**
   - LD information is used in imputation, a computational method to estimate genotypes at untyped SNPs based on the LD patterns observed in a reference panel.

In summary, LD in SNPs has practical implications for the design and interpretation of genetic studies. It facilitates the identification of informative markers, aids in the mapping of disease genes, and contributes to our understanding of the genomic architecture and evolutionary history of populations.

### Hardy-Weinberg Equilibrium

The Hardy-Weinberg equilibrium is a fundamental concept in population genetics that describes the relationship between the frequencies of alleles in a population and predicts how those frequencies will change over time in the absence of certain factors.

In simple terms, it assumes that a population is not evolving, which means that there are no influences such as mutation, selection, migration, or genetic drift. In such a hypothetical scenario, the frequencies of alleles in the population remain constant from generation to generation.

The Hardy-Weinberg equilibrium is expressed by the equation:

\[ p^2 + 2pq + q^2 = 1 \]

Here, \( p \) and \( q \) represent the frequencies of two alleles for a particular gene in a population, and \( p^2 \), \( 2pq \), and \( q^2 \) represent the frequencies of the three possible genotypes. The equation states that the sum of the frequencies of all genotypes must equal 1 under the assumption of a non-evolving population.

This equilibrium is a useful baseline for understanding genetic processes. If a population's observed genotype frequencies deviate from the expected frequencies under Hardy-Weinberg equilibrium, it suggests that some evolutionary forces are at play, and researchers can investigate the underlying factors influencing the population's genetic composition.

### SNP selection criteria used to establish a SNP panel for re-identification
- MAF > 0.2
- HW p-value > 0.01
- LD < 0.01
- Exclusion of SNPs from HLA (SNP calls prone to artifacts)

# Where we are
## Variant k-mer generation
- Generating all the k-mers containing 1 SNP from dbSNP
- Each k-mer is carrying variant which is a biomarker that can be used for identification
- We generate a number of k k-mers in the better case
- To make sure that a k-mer of the index is present in the raw genomic data, we must find a sufficient number of our k-mers corresponding to the sequencing depth.

## Measures for the score

- Number of SNPs
- SNPs with a MAF close to 0.5
- Probability of Identity



### Other measures that can be interresting after
- Linkage Disequilibrium (LD) can be interresting to infer more SNPs.
- What are the SNPs related to