#!/usr/bin/python3

class Variant :

    def __init__(self, rsid:str, chr:str, snp:str, snpPos:int, relPos:int, kmersCount:int, ambiguousKmersCount:int) -> None:
        self._rsid = rsid
        self._chr = chr
        self._snp = snp
        self._snpPos = snpPos
        self._relPos = relPos
        self._kmersCount = kmersCount
        self._ambiguousKmersCount = ambiguousKmersCount

    """def __init__(self) -> None:
        self._rsid = "rs0"
        self._chr = "Z"
        self._snp = "N"
        self._snpPos = 0
        self._relPos = 0
        self._kmersCount = 0
        self._ambiguousKmersCount = 0"""

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
        ppty = f"{self._rsid}\t{self._chr}\t{self._snp}\t{self._snpPos}\t{self._relPos}\t{self._kmersCount}\t{self._ambiguousKmersCount}"
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

    