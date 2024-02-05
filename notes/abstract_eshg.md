# Reducing risk of data genomic privacy breach using k-mer approach
# Reducing re-identification risk of genomic data using k-mer approach

## Points à aborder:

1. séquençage courant
2. augmentation de la qualité et quantité de données, augmentation du risque
3. nécessité du partage des données
4. nécessité de réduire le risque associé au partage
5. méthode: matrice de k-mer, détection des SNPs, estimation du risque
6. résultats attendus: donner une idée du risque de partage
7. conclusion/perspectives: aller plus loin avec l'inférence.

**Introduction:**

With the remarkable advancements in sequencing techniques and the significant reduction in their costs, the field of genomics has experienced an unprecedented increase in the volume of genomic data available for research purposes. This surge in genomic data holds immense potential for enhancing our understanding of complex diseases, advancing personalized medicine, and fostering innovations in health care. In this goal, several countries launched and developped genomic data collection programs, and companies offering Direct-to-Consumer (DtC) genetic testing (Anecestry, 23andMe) lead to people sharing their genomic data online (Personal Genome Project, OpenSNP). However, it simultaneously poses substantial privacy risks, as the highly sensitive nature of genomic information could lead to re-identification or potential misuse if not adequately protected. The balance between leveraging genomic data for scientific advancement and safeguarding individual privacy is a critical concern in the bioinformatics community. Despite the importance of data sharing for fostering research and advancements in genomics, there is a pressing need to develop novel approaches that mitigate privacy risks while preserving the utility of the data for research. This study aims to address this crucial issue by exploring innovative methods that minimize privacy vulnerabilities in genomic data sharing without hindering the accessibility and processing of the data for scientific discovery. Our focus diverges from traditional data access and protection strategies, offering a fresh perspective on ensuring genome privacy in the context of open scientific research.


**Methods:**

Here, we propose a novel approach focused on the utilization of k-mer to accurately estimate the re-identification risks associated with genomic data sharing. 
The utilization of a k-mer inherently reduces the risk of re-identification by circumventing the privacy vulnerabilities posed by repetitive regions, such as transposable elements or Short Tandem Repeats (STRs), which are traditionally used to identify an individual.

Despite the advantages of using k-mers for genomic data analysis, the practice of sharing k-mers poses a residual risk of re-identification. Therefore, it is imperative to evaluate this risk, especially considering the widespread use of tools that analyze k-mers of interest for various research purposes. Assessing the potential for re-identification is crucial to ensuring the privacy and security of individuals whose genomic data is shared, as we aim to understand and mitigate any privacy concerns associated with the use of k-mers in genomic studies.

Given that k-mers are short DNA sequences (usually 31 bases in length, may vary), the primary risk of re-identification emanates from Single Nucleotide Polymorphisms (SNPs) present within these k-mers. SNPs, being the most common type of genetic variation among people, serve as identifiers and are commonly used in population genomics, DtC genetic testing, and forensic genetics.. A combination of a small number of them is sufficient to uniquely identify someone, as only 40 di-allelic SNPs with an heterozygosity equal to 0.5 gives a theoretical average match probability of $10^-17$ (Pakstis et al, 2007).

Our analysis begins with the detection of SNPs within the k-mers, utilizing publicly available databases such as dbSNP. This step is crucial for identifying which k-mers contain variations that could contribute to re-identification risk. Following SNP identification, we estimate the probability of individual identity by analyzing the frequency and distribution of these SNPs. Our goal in this process is the establishment of a risk score: a quantitative measure reflecting the likelihood of re-identification from shared k-mer matrices.

**Expected Results:**

Through our k-mer based approach, we anticipate developing a comprehensive risk score that will enable researchers to assess the re-identification risks associated with sharing genomic data in the form of k-mer sets. 

Our work is focused on estimating the re-identification risk inherent in sharing k-mer sets. By the conclusion of our project, we aim to offer a tool that mitigates this risk by facilitating the filtration of subsets of sensitive k-mers. This approach will allow for the reduction of privacy risks associated with genomic data sharing, without diminishing the data's value for predicting genetic diseases. The tool's capability to filter sensitive k-mers will empower researchers to share genomic data more safely, ensuring that privacy considerations are adequately addressed.

The validation of our approach would mark a significant advancement in the field, offering a pragmatic solution to the privacy concerns that hinder the sharing of genomic data. Moving forward, the next phase of our research will involve leveraging the detected SNPs within the k-mer sets to further assess re-identification risks. We plan to explore additional methods, such as inference or imputation techniques using tools like IMPUTE2, to refine our risk assessment.

**Conclusion:**

The endeavor to balance the advancement of genomic research with the imperative of protecting individual privacy is at the forefront of our study. A pivotal challenge in ensuring the reproducibility of experiments involving patient data lies in the complexities of data sharing. Our project addresses this challenge head-on by introducing a metric that quantifies the risk of re-identification, coupled with a practical solution for excising overly sensitive data from datasets. This approach should enhance the privacy and security of genomic data but also represents a contribution to open science and the reproducibility of scientific research. By providing these tools, we aim to foster a safer environment for data sharing, thereby supporting the continuous progress of genomic studies and the broader scientific community's pursuit of knowledge.