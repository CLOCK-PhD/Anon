#!/usr/bin/python3

class Variant :

    def __init__(self, rsid:str, chr:str, snp:str, snpPos:int, relPos:int, kmersCount:int, ambiguousKmersCount:int) -> None:
        """Générateur de l'objet Variant.

        Parameters :
            rsid                (str):  Identifiant du SNP
            chr                 (str):  Chromosome où se situe le SNP
            snp                 (str):  Variation du SNP
            snpPos              (int):  Position du SNP sur le chromosome
            relPos              (int):  Position relative du k-mer sur le chromosome
            kmersCount          (int):  Nombre de k-mers créés à partir du SNP
            ambiguousKmersCount (int):  Nombre de k-mers identiques
        """
        self._rsid = rsid
        self._chr = chr
        self._snp = snp
        self._snpPos = snpPos
        self._relPos = relPos
        self._kmersCount = kmersCount
        self._ambiguousKmersCount = ambiguousKmersCount

    # Getters / Setters
    @property
    def rsid(self):
        return self._rsid

    @rsid.setter
    def rsid(self, id:str):
        self._rsid = id

    @property
    def chr(self):
        return self._chr

    @chr.setter
    def chr(self, chrom:str):
        self._chr = chrom

    @property
    def snp(self):
        return self._snp

    @snp.setter
    def snp(self, s:str):
        self._snp = s

    @property
    def snpPos(self):
        return self._snpPos

    @snpPos.setter
    def snpPos(self, pos:int):
        self._snpPos = pos

    @property
    def relPos(self):
        return self._relPos

    @relPos.setter
    def relPos(self, pos:int):
        self._relPos = pos

    @property
    def kmersCount(self):
        return self._kmersCount

    @kmersCount.setter
    def kmersCount(self, count:int):
        self._kmersCount = count

    @property
    def ambiguousKmersCount(self):
        return self._ambiguousKmersCount

    @ambiguousKmersCount.setter
    def ambiguousKmersCount(self, count:int):
        self._ambiguousKmersCount = count

    @property
    def variantProperties(self):
        """# Afficher les propriétés de l'objet dans une chaine de caractères
        Ancienne version : toutes les caractéristiques
        Nouvelle version : 
            - rs_id
            - Position de la variation dans le kmer (snpPos - relPos)
            - kmers_count
            - ambiguous k-mers count
        """
        #ppty = f"{self._rsid}\t{self._chr}\t{self._snp}\t{self._snpPos}\t{self._relPos}\t{self._kmersCount}\t{self._ambiguousKmersCount}"
        # Pour indiquer la position du variant dans le k-mer:
        #varPosInKmer = self._snpPos - self._relPos
        ppty = f"{self._rsid}\t{self._kmersCount}\t{self._ambiguousKmersCount}"
        return ppty

if __name__ == '__main__':
    print("le chat")

    variant = Variant("rs1", "X", "T", 1, 1, 0, 0)

    # tests getters/setters
    print("RSID")
    print(variant.rsid)
    variant.rsid = "rs123"
    print(variant.rsid)

    print("CHROMOSOME")
    print(variant.chr)
    variant.chr = "Y"
    print(variant.chr)

    print("SNP")
    print(variant.snp)
    variant.snp = "N"
    print(variant.snp)

    print("SNP POS")
    print(variant.snpPos)
    variant.snpPos = 20
    print(variant.snpPos)

    print("RelPos")
    print(variant.relPos)
    variant.relPos = 5
    print(variant.relPos)

    print("kmers count")
    print(variant.kmersCount)
    variant.kmersCount = 100
    print(variant.kmersCount)

    print("ambiguousKmersCount")
    print(variant.ambiguousKmersCount)
    variant.ambiguousKmersCount = 2000
    print(variant.ambiguousKmersCount)

    variant.variantProperties

    print("LES CHATS")

    