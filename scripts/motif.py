import re


class AbstractSequenceMotif(object):
    def __init__(self, motif, dna_sequence, offset):
        self.dna_sequence         = dna_sequence.upper()
        self.forward_positions    = [m.start() for m in re.finditer(motif, dna_sequence)]
        self.complement_positions = [m + offset for m in self.forward_positions]

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

    @staticmethod
    def _complementBase(base):
        if base == "A":
            return "T"
        elif base == "C":
            return "G"
        elif base == "G":
            return "C"
        elif base == "T":
            return "A"
        elif base == "N":
            return "N"
        else:
            raise RuntimeError


class CpG(AbstractSequenceMotif):
    def __init__(self, dna_sequence):
        super(CpG, self).__init__(motif="CG", dna_sequence=dna_sequence, offset=1)

    def forwardSubstitutedSequence(self, ambig_nuc):
        return self._substituteSequence(self.forward_positions, self.dna_sequence, ambig_nuc, "C")

    def complementSubstitutedSequence(self, ambig_nuc):
        complement_seq = "".join([self._complementBase(b) for b in list(self.dna_sequence)])
        return self._substituteSequence(self.complement_positions, complement_seq, ambig_nuc, "C")


def getMotif(key, dna_sequence):
    if key == "CG":
        return CpG(dna_sequence)
