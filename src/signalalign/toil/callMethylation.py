import os
import uuid
import string
import subprocess
import fileinput

import pandas as pd
import numpy as np

from itertools import chain

from bd2k.util.humanize import human2bytes
from toil_lib import require
from toil_lib.programs import docker_call

from margin.toil.localFileManager import \
    LocalFile,\
    LocalFileManager,\
    urlDownloadToLocalFile,\
    unzipLocalFile,\
    deliverOutput

from margin.toil.shardAlignment import shardSamJobFunction
from margin.toil.realign import DOCKER_DIR

from signalalign import exonerateCigarWithStrandOrientation
from signalalign.motif import getMotif, getVariantCallFunctions


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
        disk       = (2 * config["reference_FileStoreID"].size)
        memory     = (6 * aln_shard.FileStoreID.size)
        batch_disk = human2bytes("250M") + config["reference_FileStoreID"].size
        methylation_probs = job.addChildJobFn(shardSamJobFunction,
                                              config, aln_shard, None,
                                              calculateMethylationProbabilityJobFunction,
                                              callMethylationJobFunction,
                                              exonerateCigarStringFn=exonerateCigarWithStrandOrientation,
                                              batch_disk=batch_disk,
                                              disk=disk, memory=memory).rv()
        all_methylation_probs.append(methylation_probs)
        count += 1
    job.fileStore.logToMaster("[signalAlignJobFunction]Issued methylation calling for %s alignment shards"
                              % count)
    job.addFollowOnJobFn(consolidateVariantCallsJobFunction, config, all_methylation_probs)
    return


def consolidateVariantCallsJobFunction(job, config, posterior_prob_fids):
    variants  = getVariantCallFunctions(config["degenerate"])
    parser    = variants.parseVariantCalls
    file_iter = (job.fileStore.readGlobalFile(fid) for fid in posterior_prob_fids)
    table     = pd.concat([parser(f) for f in file_iter]).sort_values(["contig", "ref_pos"])
    outfile   = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                          filename="%s_%s.tsv" % (config["sample_label"], config["degenerate"]))
    _handle   = open(outfile.fullpathGetter(), "w")
    variants.writeVariantCalls(table, _handle)
    _handle.close()
    deliverOutput(job, outfile, config["output_dir"])


def consolidateMethylationCallsJobFunction(job, config, methylation_prob_fids):
    outfile = LocalFile(workdir=job.fileStore.getLocalTempDir(),
                        filename="%s_%s.tsv" % (config["sample_label"], config["degenerate"]))
    _handle = open(outfile.fullpathGetter(), "w")
    files   = fileinput.input([job.fileStore.readGlobalFile(fid) for fid in methylation_prob_fids])
    map(lambda l: _handle.write(l), files)
    files.close()
    _handle.close()
    deliverOutput(job, outfile, config["output_dir"])
    return


def processReferenceSequence(ref_seq, workdir, motif_key=None, sub_char="X", parent_job=None):
    # make the forward and backward sequences, substituting the necessary motifs
    if motif_key is not None:
        motif, ok = getMotif(motif_key, ref_seq)
        require(ok, "[processReferenceSequence]Illegal motif_key given %s" % motif_key)
        if parent_job is not None:
            parent_job.fileStore.logToMaster("[processReferenceSequence]Made %s substitutions" %
                                             motif.substitutionPositionCount())
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
                                               signalMachine_image="quay.io/artrand/signalmachine"):
    def _get_url(read_label):
        try:
            return ledger[read_label]
        except KeyError:
            return None

    def _SignalMachine(read_label, cigar, nanopore_read):
        guide_aln = LocalFile(workdir=workdir)
        _handle   = open(guide_aln.fullpathGetter(), "w")
        _handle.write(cigar)
        _handle.close()
        require(os.path.exists(guide_aln.fullpathGetter()), "NO guide aln file")
        signalMachine_args = [
            "--sm3Hdp",
            "-s", "1",
            "-o", "%s" % degenerate_enum,
            "-L", "%s" % read_label,
            "-T", "%s%s" % (DOCKER_DIR, models.localFileName(hmmfid)),
            "-q", "%s%s" % (DOCKER_DIR, nanopore_read.filenameGetter()),
            "-f", "%s%s" % (DOCKER_DIR, fw_seqfile.filenameGetter()),
            "-b", "%s%s" % (DOCKER_DIR, bw_seqfile.filenameGetter()),
            "-p", "%s%s" % (DOCKER_DIR, guide_aln.filenameGetter()),
            "-u", "%s%s" % (DOCKER_DIR, posteriors.filenameGetter()),
            "-v", "%s%s" % (DOCKER_DIR, models.localFileName(hdpfid)),
        ]
        try:
            docker_call(job=job,
                        tool=signalMachine_image,
                        parameters=signalMachine_args,
                        work_dir=(workdir + "/"))
        except subprocess.CalledProcessError:
            pass

    def _parse_probabilities():
        return pd.read_table(posteriors.fullpathGetter(),
                             usecols=(1, 2, 3, 6),
                             names=["ref_pos", "base", "posterior", "read_label"],
                             dtype={"ref_pos"    : np.int,
                                    "base"       : np.str,
                                    "posterior"  : np.float64,
                                    "read_label" : np.str})

    def _sumExpectationsOverColumns():
        f  = LocalFile(workdir=workdir)
        _h = open(f.fullpathGetter(), "w")
        for pos, pos_df in aligned_pairs.groupby(["ref_pos"]):
            for base, base_df in pos_df.groupby("base"):
                marginal_prob = base_df["posterior"].sum()
                coverage      = len(base_df["read_label"].unique())
                l = "%s\t%s\t%s\t%s\t%s\n" % (cPecan_config["contig_name"], pos, base, marginal_prob, coverage)
                _h.write(l)
        _h.close()
        return f

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
    url_iter  = (_get_url(l.strip()) for l in cPecan_config["query_labels"])
    read_urls = [u for u in url_iter if u is not None]

    if config["debug"]:
        job.fileStore.logToMaster("[calculateMethylationProbabilityJobFunction]Got %s URLs" % len(read_urls))

    npReads   = [unzipLocalFile(f)
                 for f in [urlDownloadToLocalFile(job, workdir, url) for url in read_urls]
                 if f is not None]
    failed    = len(read_urls) - len(npReads)

    if failed > 0 and config["stop_at_failed_reads"]:
        raise RuntimeError("[calculateMethylationProbabilityJobFunction]Got %s failed npRead"
                           "downloads and stop_at_failed_reads is True" % failed)
    else:
        if config["debug"]:
            job.fileStore.logToMaster("[calculateMethylationProbabilityJobFunction]"
                                      "Failed to download and upzip %s NanoporeReads" % failed)

    # file to collect the posterior probs
    posteriors      = LocalFile(workdir=workdir, filename="%s_%s.dat" % (config["sample_label"], uuid.uuid4()))
    degenerate_enum = getVariantCallFunctions(config["degenerate"]).enum()

    # do the signal alignment, and get the posterior probabilities
    map(lambda (l, c, n): _SignalMachine(l.strip(), c, n), zip(cPecan_config["query_labels"],
                                                               cPecan_config["exonerate_cigars"],
                                                               npReads))

    # the reads may not produce any posteriors, if, for example, they don't align to a region where
    # there are any ambiguity characters the posteriors file will be empty and we just return
    # None, which is the convention
    if not os.path.exists(posteriors.fullpathGetter()) or os.stat(posteriors.fullpathGetter()).st_size == 0:
        return None

    # reminder: the convention is that 'expectations' are un-normalized posterior probabilities
    # so this file is a table of expectatiosn, I also use the convention that the trailing
    # underscore means `file` or `file-path`
    aligned_pairs = _parse_probabilities()
    expectations_ = _sumExpectationsOverColumns()
    if config["probs_output_dir"] is not None:
        deliverOutput(job, posteriors, config["probs_output_dir"])

    return job.fileStore.writeGlobalFile(expectations_.fullpathGetter())


def callMethylationJobFunction(job, config, alignment_shard, prob_fids):
    def parse_expectations(_file):
        return pd.read_table(_file, usecols=(0, 1, 2, 3, 4),
                             names=["contig", "ref_pos", "base", "prob", "coverage"],
                             dtype={"contig"    : np.str,
                                    "ref_pos"   : np.int,
                                    "base"      : np.str,
                                    "prob"      : np.float64,
                                    "coverage"  : np.int})

    def write_down_methylation(contig, position, posterior_probs):
        base_call = max(posterior_probs, key=posterior_probs.get)  # the base with highest marginal prob
        _handle.write("%s\t%s\t%s\t" % (contig, position, base_call))
        for base in base_options:
            try:
                prob, coverage = posterior_probs[base]
                _handle.write("%s,%s,%s\t" % (base, prob, coverage))
            except KeyError:
                _handle.write("%s,%s,%s\t" % (base, 0.0, 0))
        _handle.write("\n")

    job.fileStore.logToMaster("[callMethylationJobFunction]Calling methylation")
    prob_fids = [x[0] for x in prob_fids if x[0] is not None]

    # inner list comp: downloads the proability file outer list comp: parses the files, then 
    # they are concatenated
    expectations = pd.concat([parse_expectations(f)
                              for f in [job.fileStore.readGlobalFile(f)
                                        for f in prob_fids]])

    # to decide on the methylation status to call, we first sum up the expectations for all 
    # bases at a given reference position. Then we divide the total expectations for each 
    # base by the total to give the final (normalized) posterior probability for the base
    # at that position given the reference (and the model... etc)
    base_options = getVariantCallFunctions(config["degenerate"]).baseOptions()
    calls_file   = job.fileStore.getLocalTempFile()
    _handle      = open(calls_file, "w")
    for contig, contig_df in expectations.groupby(["contig"]):
        for ref_pos, pos_df in contig_df.groupby(["ref_pos"]):
            if ref_pos >= alignment_shard.start and ref_pos < alignment_shard.end:
                posterior_probabilities = {}
                total_prob = pos_df["prob"].sum()
                for base, base_df in pos_df.groupby(["base"]):
                    posterior_probability = base_df["prob"].sum() / total_prob
                    coverage = base_df["coverage"].max()
                    posterior_probabilities[base] = (posterior_probability, coverage)
                write_down_methylation(contig, ref_pos, posterior_probabilities)
            else:
                pass
    _handle.close()
    return job.fileStore.writeGlobalFile(calls_file)
