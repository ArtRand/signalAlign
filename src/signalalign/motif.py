"""Classes and functions for handling variant calling
"""
import re

from pandas import read_table


class AbstractSequenceMotif(object):
    def __init__(self, motif, dna_sequence, offset):
        self.dna_sequence         = dna_sequence.upper()
        self.forward_positions    = [m.start() for m in re.finditer(motif, self.dna_sequence)]
        self.complement_positions = [m + offset for m in self.forward_positions]
        #correct_positiosn         = [m.start() for m in re.finditer(motif, self.dna_sequence)]

    def forwardPositions(self):
        return self.forward_positions

    def complementPositions(self):
        return self.complement_positions

    def forwardSubstitutedSequence(self, sub_char):
        raise NotImplementedError

    def complementSubstitutedSequence(self, sub_char):
        raise NotImplementedError

    def _checkSubstitution(self, position):
        raise NotImplementedError

    def _substituteSequence(self, positions, sequence, ambig_nuc, orig_nuc):
        l_seq = list(sequence)
        for pos in positions:
            assert(l_seq[pos] == orig_nuc), "Illegal substitution of %s with %s" % (l_seq[pos], ambig_nuc)
            l_seq[pos] = ambig_nuc
        return "".join(l_seq)

    def substitutionPositionCount(self):
        return len(self.forward_positions)

    @staticmethod
    def _complementBase(base):
        # canonicals
        if base == "A":
            return "T"
        elif base == "C":
            return "G"
        elif base == "G":
            return "C"
        elif base == "T":
            return "A"
        # two redundant
        elif base == "M":
            return "K"
        elif base == "K":
            return "M"
        elif base == "R":
            return "Y"
        elif base == "Y":
            return "R"

        # three redundant
        elif base == "B":
            return "V"
        elif base == "V":
            return "B"
        elif base == "D":
            return "H"
        elif base == "H":
            return "D"

        elif base == "N":
            return "N"
        else:
            raise RuntimeError("Got base I don't understand %s" % base)

    @staticmethod
    def parseVariantCalls(file_):
        raise NotImplementedError


class CpG(AbstractSequenceMotif):
    def __init__(self, dna_sequence):
        super(CpG, self).__init__(motif="CG", dna_sequence=dna_sequence, offset=1)

    def forwardSubstitutedSequence(self, ambig_nuc):
        return self._substituteSequence(self.forward_positions, self.dna_sequence, ambig_nuc, "C")

    def complementSubstitutedSequence(self, ambig_nuc):
        complement_seq = "".join([self._complementBase(b) for b in list(self.dna_sequence)])
        return self._substituteSequence(self.complement_positions, complement_seq, ambig_nuc, "C")


class AbstractDegenerateType(object):
    def __init__(self):
        pass  # just a function container

    @staticmethod
    def enum():
        raise NotImplementedError

    @staticmethod
    def baseOptions():
        raise NotImplementedError

    @staticmethod
    def parseVariantCalls(file_):
        raise NotImplementedError

    def writeVariantCalls(outfile):
        raise NotImplementedError


class CytosineMethylationVariantCall(AbstractDegenerateType):
    def __init__(self):
        pass  # just a function container

    @staticmethod
    def enum():
        return 0

    @staticmethod
    def baseOptions():
        return ("C", "E")

    @staticmethod
    def parseVariantCalls(file_):
        return read_table(file_,
                          usecols=(0, 1, 2, 3, 4),
                          names=["contig", "ref_pos", "variant_call", "pC", "p5mC"],
                          header=None,
                          sep="\t")

    @staticmethod
    def writeVariantCalls(table, outfile):
        for _, row in table.iterrows():
            outfile.write("%s\t%s\t%s\t%s\t%s\n" % (row["contig"],
                                                    row["ref_pos"],
                                                    row["variant_call"],
                                                    row["pC"],
                                                    row["p5mC"]))
        return


class CanonicalNucleotideVariantCall(AbstractDegenerateType):
    def __init__(self):
        pass  # just a function container

    @staticmethod
    def enum():
        return 3

    @staticmethod
    def baseOptions():
        return ("A", "C", "G", "T")

    @staticmethod
    def parseVariantCalls(file_):
        raise NotImplementedError


class DegenerateEnum(object):
    """represents a mapping of string keys to enums for signalAlign variant calling
    """
    def __init__(self):
        self.degenerate_type = {
            "cytosine2": 0,
            "cytosine3": 1,
            "adenosine": 2,
            "variant": 3,
        }

    def get(self, degnerate_key):
        # throws KeyError
        return self.degenerate_type[degnerate_key]

    def check(self, degenerate_key):
        return degenerate_key in self.degenerate_type.keys()


# getters
def getMotif(key, dna_sequence):
    """returns a SequenceMotif object and an error
    """
    if key == "CG":
        return CpG(dna_sequence), True
    else:
        return None, False


def getDegenerateEnum(degenerate_key):
    # throws RuntimeError
    enum = DegenerateEnum()
    ok   = enum.check(degenerate_key)
    if not ok:
        raise RuntimeError("[getDegenerateEnum]Degenerate key %s not allowed, options are"
                           "%s " % (degenerate_key, enum.keys()))
    else:
        return enum.get(degenerate_key)


def getVariantCallFunctions(degenerate_key):
    if degenerate_key == "cytosine2":
        return CytosineMethylationVariantCall()
    elif degenerate_key == "variant":
        return CanonicalNucleotideVariantCall()
    else:
        raise RuntimeError("[getVariantCallFunctions]Degenerate %s not implemented" % degenerate_key)


def checkDegenerate(degenerate_key):
    try:
        getVariantCallFunctions(degenerate_key)
    except RuntimeError:
        return False
    return True
