import os
import subprocess


class Bwa(object):
    """Wrapper of BWA aligner, requires bwa to be in path.
    Citation:
        Program: bwa (alignment via Burrows-Wheeler transformation)
        Contact: Heng Li <lh3@sanger.ac.uk>
    """
    def __init__(self, target):
        self.target = target
        self.db_handle = ''

    def build_index(self, destination, output=None):
        self.db_handle = destination + '/temp_bwaIndex'
        cmd = "bwa index -p {0} {1}".format(self.db_handle, self.target)
        if output is None:
            output = open(os.devnull, 'w')
        else:
            output = open(output, 'w')
        try:
            subprocess.check_call(cmd.split(), stdout=output, stderr=output)
            output.close()
            return True
        except subprocess.CalledProcessError:
            output.close()
            return False

    @staticmethod
    def suffixes():
        return [".amb", ".ann", ".bwt", ".pac", ".sa"]

    @staticmethod
    def align(bwa_index, query, output_sam_path, outerr=None):
        for suff in Bwa.suffixes():
            assert os.path.exists(bwa_index + suff),\
                "[Bwa:.lign] Didn't find index files {}".format(bwa_index + suff)
        assert os.path.exists(query), "[Bwa::align] Didn't find query file {}".format(query)
        cmd = "bwa mem -x ont2d {idx} {query}".format(idx=bwa_index, query=query)
        if outerr is None:
            outerr = open(os.devnull, 'w')
        else:
            outerr = open(outerr, 'w')
        try:
            with open(output_sam_path, 'w') as fH:
                fH.write(subprocess.check_output(cmd.split(), stderr=outerr))
            outerr.close()
            return True
        except subprocess.CalledProcessError:
            outerr.close()
            return False


def getBwaIndex(reference, dest, output=None):
    bwa = Bwa(reference)
    bwa.build_index(dest, output=output)
    bwa_ref_index = dest + "temp_bwaIndex"
    return bwa_ref_index
