#!/usr/bin/python3

# En fait on va pas utiliser cette classe.

from variant import Variant

class Kmer :

    def __init__(self, sequence:str) -> None:
        self._sequence = sequence
        self._variants = []
        self._in_genome = False

    # Getters/Setters
    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        self._sequence = seq

    @property
    def variants(self):
        return self._variants

    @property
    def in_genome(self):
        return self._in_genome

    @in_genome.setter
    def in_genome(self, ingen:bool):
        self.in_genome = ingen

    # Fonctions
    def addVariant(self, var:Variant):
        self._variants.append(var)

    def show_kmer(self):
        print(self._sequence)
        for v in self._variants:
            print(v.variantProperties)
        print(f"In genome ? : {self._in_genome}")
        

if __name__ == '__main__':
    print("le chat")

    # Test pour l'import de la classe Variant
    variant1 = Variant("rs1", "X", "T", 1, 1, 0, 0)
    print(variant1.rsid)
    print(variant1.chr)

    variant2 = Variant("rs2", "Y", "A", 2, 2, 1, 1)

    # Création de la classe :
    kmer = Kmer("ATGC")
    print(kmer.sequence)
    print(kmer.in_genome)
    kmer.addVariant(variant1)
    print(kmer.variants)
    kmer.addVariant(variant2)
    print(kmer.variants)
    print(kmer.variants[0].rsid)
    print(kmer.variants[1].rsid)
    print(f"longueur séquence : {len(kmer.sequence)}")

    kmer.show_kmer()

