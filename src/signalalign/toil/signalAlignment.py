import os
import string
from itertools import chain

from toil_lib import require
from toil_lib.programs import docker_call

from margin.toil.localFileManager import \
    LocalFile,\
    LocalFileManager,\
    urlDownloadToLocalFile,\
    unzipLocalFile,\
    deliverOutput
from margin.toil.realign import shardSamJobFunction, DOCKER_DIR

from signalalign.motif import getMotif, degenerateEnum


def _reverseComplement(dna, reverse=True, complement=True):
    trans_table = string.maketrans('ATGCatgc', 'TACGtacg')
    # Make complement to DNA
    comp_dna = dna.translate(trans_table)
    # Output all as strings
    if reverse and complement:
        return comp_dna[::-1]
    if reverse and not complement:
        return dna[::-1]
    if complement and not reverse:
        return comp_dna
    if not complement and not reverse:
        return dna


def signalAlignJobFunction(job, config, alignment_shards):
    alignment_shards = chain(*alignment_shards)
    # each shard is a region of the genome/chromosome/contig and can be methylation called
    # independently
    all_methylation_probs = []  # contains the methylation probabilites for all of the shards together
    count = 0
    for aln_shard in alignment_shards:
        # TODO disk requirement
        # TODO memory requirement
        # TODO cores requirement
        methylation_probs = job.addChildJobFn(shardSamJobFunction,
                                              config, aln_shard, None,
                                              calculateMethylationProbabilityJobFunction,
                                              callMethylationJobFunction).rv()
        all_methylation_probs.append(methylation_probs)
        count += 1
    job.fileStore.logToMaster("[signalAlignJobFunction]Issued methylation calling for %s alignment shards"
                              % count)
    # TODO followOn 'consolidate job function'
    return


def processReferenceSequence(ref_seq, workdir, motif_key=None, sub_char="X"):
    # make the forward and backward sequences, substituting the necessary motifs
    if motif_key is not None:
        motif, ok = getMotif(motif_key, ref_seq)
        require(ok, "[processReferenceSequence]Illegal motif_key given %s" % motif_key)
        try:
            fw_refseq = motif.forwardSubstitutedSequence(sub_char)
            bw_refseq = motif.complementSubstitutedSequence(sub_char)
        except AssertionError:
            return None, None, False
    else:
        fw_refseq = ref_seq.upper()
        bw_refseq = _reverseComplement(fw_refseq, reverse=False, complement=True)

    fw_refseqfile  = LocalFile(workdir=workdir)
    bw_refseqfile  = LocalFile(workdir=workdir)
    sequences      = [fw_refseq, bw_refseq]
    sequence_files = [fw_refseqfile, bw_refseqfile]

    for f, s in zip(sequence_files, sequences):
        _h = open(f.fullpathGetter(), "w")
        _h.write(s + "\n")
        _h.close()

    [require(os.path.exists(f.fullpathGetter()), "[processReferenceSequence]Missing %s" % f.filenameGetter())
        for f in sequence_files]

    return fw_refseqfile, bw_refseqfile, True


def calculateMethylationProbabilityJobFunction(job, config, cPecan_config, ignore_hmm, batch_number,
                                               signalMachine_image="506217c84696"):
    def runSignalMachine(read_label, cigar, nanopore_read):
        guide_aln = LocalFile(workdir=workdir)
        _handle   = open(guide_aln.fullpathGetter(), "w")
        _handle.write(cigar)
        _handle.close()
        require(os.path.exists(guide_aln.fullpathGetter()), "NO guide aln file")
        signalMachine_args = [
            "--sm3Hdp",
            "-s", "1",
            "-o", "%s" % degenerateEnum(config["degenerate"]),
            "-L", "%s" % read_label,
            "-T", "%s%s" % (DOCKER_DIR, models.localFileName(hmmfid)),
            "-q", "%s%s" % (DOCKER_DIR, nanopore_read.filenameGetter()),
            "-f", "%s%s" % (DOCKER_DIR, fw_seqfile.filenameGetter()),
            "-b", "%s%s" % (DOCKER_DIR, bw_seqfile.filenameGetter()),
            "-p", "%s%s" % (DOCKER_DIR, guide_aln.filenameGetter()),
            "-u", "%s%s" % (DOCKER_DIR, posteriors.filenameGetter()),
            "-v", "%s%s" % (DOCKER_DIR, models.localFileName(hdpfid)),
        ]
        docker_call(job=job,
                    tool=signalMachine_image,
                    parameters=signalMachine_args,
                    work_dir=(workdir + "/"))

    job.fileStore.logToMaster("[calculateMethylationProbabilityJobFunction]Running on batch %s" % batch_number)
    workdir = job.fileStore.getLocalTempDir()
    fw_seqfile, bw_seqfile, ok = processReferenceSequence(cPecan_config["contig_seq"],
                                                          workdir,
                                                          config["motif_key"],
                                                          config["substitute_char"])
    if not ok:
        raise RuntimeError("[calculateMethylationProbabilityJobFunction]ERROR processing reference sequences")
    # get the models
    hmmfid = config["HMM_fid"]
    hdpfid = config["HDP_fid"]
    try:
        models = LocalFileManager(job=job, fileIds_to_get=[hmmfid, hdpfid], workdir=workdir)
    except AssertionError:
        raise RuntimeError("[calculateMethylationProbabilityJobFunction]ERROR getting models locally")
    # download the npRead files
    ledger    = config["ledger"]
    read_urls = [ledger[label.strip()] for label in cPecan_config["query_labels"]]
    # TODO could make this error check, and continue even when some of the reads aren't downloaded
    npReads    = [unzipLocalFile(f) for f in [urlDownloadToLocalFile(job, workdir, url) for url in read_urls]]
    posteriors = LocalFile(workdir=workdir)
    #map(lambda l, c, n: runSignalMachine(l, c, n), zip(cPecan_config["query_labels"],
    #                                                   cPecan_config["exonerate_cigars"],
    #                                                   npReads))
    for label, cigar, npread in zip(cPecan_config["query_labels"], cPecan_config["exonerate_cigars"], npReads):
        runSignalMachine(label.strip(), cigar, npread)

    deliverOutput(job, posteriors, config["output_dir"])


def callMethylationJobFunction(job, config, alignment_shard, methylation_probs):
    job.fileStore.logToMaster("[callMethylationJobFunction]Calling methylation")
