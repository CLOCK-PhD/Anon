# SNP Severity

considering the individual characteristics of the real
SNPs (i.e., their severity levels), we study the level of genomic
privacy of a patient against a curious party at the SPU. The
severity of a SNP i can be defined as the privacy-sensitivity of
the SNP when SNP P
i = 1 (i.e., when it exists as a variant at the
patient P). For example, a real SNP revealing the predisposition
of a patient for Alzheimer’s disease can be considered more
severe than another real SNP revealing his predisposition to a
more benign disease. Severity values of the SNPs are determined
as a result of medical studies (depending on their contributions
to various diseases) and tables of disease severities provided by
insurance companies (e.g., percentage of invalidity). We denote
the severity of a real SNP i as V i , and 0 ≤ V i ≤ 1 (1 denotes
the highest severity). Thus, we define the genomic privacy of the
patient P as below:

We do not use the traditional entropy metric [21], [22] to quantify
privacy, as only one state of SNP P
i poses privacy risks (i.e.,
SNP P
i = 1 ), as discussed before.
First, we study the relationship between the storage redun-
dancy and the severity of the real SNPs by focusing on three
types of patients: (i) patient A, carrying mostly low severity real
SNPs (in Υ A ), (ii) patient B, carrying mostly high severity real
SNPs (in Υ B ), and (iii) patient C, carrying mixed severity real
SNPs (in Υ C ). For each patient, the highest level of privacy is
achieved when the storage redundancy is maximum (i.e., when
all potential SNPs of the patient are stored at the SPU). Thus, we
recognize this level as 100% genomic privacy for the patient. For
the evaluation, we take the highest privacy level of patient C as
the base and normalize everything with respect to this value. We
use the following parameters for the simulation. The severities of
patient A’s and patient B’s real SNPs are represented as truncated
Gaussian random variables with (µ A , σ A ) = (0.25, 0.15) and
(µ B , σ B ) = (0.75, 0.15) , respectively. Furthermore, the severity
of patient C’s real SNPs are represented as a uniform distribution
between 0 and 1 . We also set µ(l) = 0.8 , σ(l) = 0.25 , µ(k) = 2 ,
and σ(k) = 0.75 . In Fig. 7, we illustrate the increase in privacy
with increments in the storage redundancy for these three types
of patients (A, B, and C). We observe that by increasing the
storage redundancy, a patient with high severity real SNPs gains
more privacy than a patient with lower severity real SNPs,
hence the storage redundancy can be customized for each patient
differently based on the types of his real SNPs. It can be argued
that the amount of storage redundancy for a patient can leak
information (to the curious party the SPU) about the severities
of his real SNPs. However, the severity of the SNPs is not the
only criteria to determine the storage redundancy for a desired
level of genomic privacy as we discuss next.