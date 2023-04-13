#!/usr/bin/python3

"""
Encodeur JSON pour la classe Kmer
"""

import json
from json import JSONEncoder

class KmerEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__

