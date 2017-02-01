import os
import uuid
import shutil
import gzip
import tarfile
import cPickle

from toil_lib import require
from toil_lib.files import tarball_files

from margin.toil.localFileManager import LocalFile, urlDownload, deliverOutput

from signalalign.nanoporeRead import NanoporeRead


def makeReadstoreJobFunction(job, config, samples):
    ledger_fids = [job.addChildJobFn(makeNanoporeReadLedgerJobFunction, s,
                                     config["batchsize"],
                                     config["readstore_dir"]).rv()
                   for s in samples]
    job.addFollowOnJobFn(deliverLedgerJobFunction, config, ledger_fids)
    return


def deliverLedgerJobFunction(job, config, ledger_fids):
    fHs    = [open(job.fileStore.readGlobalFile(f), "r") for f in ledger_fids]
    ls     = [cPickle.load(f) for f in fHs]
    ledger = ls[0]
    [ledger.update(d) for d in ls[1:]]

    fn = LocalFile(workdir=job.fileStore.getLocalTempDir(), filename="%s_ledger.pkl" % config["ledger_name"])
    _h = open(fn.fullpathGetter(), "w")
    cPickle.dump(ledger, _h)
    _h.close()
    deliverOutput(job, fn, config["readstore_ledger_dir"])


def makeNanoporeReadLedgerJobFunction(job, sample, batchsize, readstore_dir):
    def check_and_write(filepath):
        require(os.path.exists(filepath), "[makeNanoporeReadLedgerJobFunction]Missing %s" % filepath)
        return job.fileStore.writeGlobalFile(filepath)

    def tar_and_make_send_batch(batch):
        tarname = "%s.tmp" % uuid.uuid4().hex
        tarpath = os.path.join(workdir, tarname)
        tarball_files(tar_name=tarname, file_paths=batch, output_dir=workdir)
        require(os.path.exists(tarpath), "[makeNanoporeReadLedgerJobFunction]Didn't make smaller tar")
        return job.fileStore.writeGlobalFile(tarpath)
    job.fileStore.logToMaster("[makeNanoporeReadLedgerJobFunction]Working on sample %s" % sample.sample_label)
    # download the sample locally
    workdir        = job.fileStore.getLocalTempDir()
    minion_archive = LocalFile(workdir=workdir, filename="%s.tmp" % uuid.uuid4().hex)
    urlDownload(job, sample.URL, minion_archive)
    require(os.path.exists(minion_archive.fullpathGetter()),
            "[makeNanoporeReadLedgerJobFunction]Didn't download tar from {url} to {path}"
            "".format(url=sample.URL, path=minion_archive.fullpathGetter()))

    # make batches and send child jobs to make NanoporeReads out of them
    tar_handle   = tarfile.TarFile(minion_archive.fullpathGetter(), "r")
    members      = tar_handle.getmembers()[1:]
    member_paths = [os.path.join(workdir, m.name) for m in members]
    tar_handle.extractall(path=workdir)
    require(batchsize <= len(member_paths),
            "[makeNanoporeReadLedgerJobFunction]Cannot split %s members into batches of %s"
            % (len(member_paths), batchsize))

    member_iter   = [member_paths[i:i + batchsize]
                     for i in xrange(0, len(member_paths), batchsize)]
    tar_fids      = map(tar_and_make_send_batch, member_iter)
    ledger_shards = [job.addChildJobFn(makeNanoporeReadsJobFunction, fid, readstore_dir, cores=0.5).rv()
                     for fid in tar_fids]

    return job.addFollowOnJobFn(consolidateLedgerShardsJobFunction, ledger_shards).rv()


def makeNanoporeReadsJobFunction(job, tar_fid, readstore_dir):
    def makeNanoporeRead(f5_path):
        # here we load the NanoporeRead and write it to a file
        np = NanoporeRead(fast_five_file=f5_path, twoD=False)
        _l = np.read_label
        tf = job.fileStore.getLocalTempFile()
        fH = open(tf, "w")
        np.write_npRead(fH)
        fH.close()
        # then we gzip it and deliver it to the readstore and return the ledger line
        fn = LocalFile(workdir=workdir, filename="%s.np.gz" % _l)
        fH = open(tf, "rb")
        gz = gzip.open(fn.fullpathGetter(), "wb")
        shutil.copyfileobj(fH, gz)
        fH.close()
        gz.close()
        deliverOutput(job, fn, readstore_dir)
        return (_l, "%s%s\n" % (readstore_dir, fn.filenameGetter()))

    def write_ledger_line(line, fH):
        l = "%s\t%s" % (line[0], line[1])  # read_label, npread URL
        fH.write(l)

    workdir = job.fileStore.getLocalTempDir()
    tar = job.fileStore.readGlobalFile(tar_fid)
    # TODO make this a function, opens the tar, extracts, and returns members
    tar_handle = tarfile.open(tar, "r:gz")
    members = tar_handle.getmembers()
    members = [os.path.join(workdir, m.name) for m in members]
    tar_handle.extractall(path=workdir)
    [require(os.path.exists(m), "[makeNanoporeReadsJobFunction]Missing member %s" % m) for m in members]
    ledger_lines = map(makeNanoporeRead, members)

    ledger_shard = job.fileStore.getLocalTempFile()
    _handle      = open(ledger_shard, "w")
    [write_ledger_line(l, _handle) for l in ledger_lines]
    _handle.close()

    return job.fileStore.writeGlobalFile(ledger_shard)


def consolidateLedgerShardsJobFunction(job, ledger_shards):
    def parse_line(line):
        l = line.strip().split("\t")
        label, url = l
        return label, url

    def parse_ledger_shard(fid):
        flat_file = job.fileStore.readGlobalFile(fid)
        _handle   = open(flat_file, "r")
        pairs     = map(parse_line, [x for x in _handle if not x.isspace()])
        return dict(pairs)

    ledger = parse_ledger_shard(ledger_shards[0])
    [ledger.update(x) for x in [parse_ledger_shard(s) for s in ledger_shards[1:]]]

    fn = job.fileStore.getLocalTempFile()
    _h = open(fn, "w")
    cPickle.dump(ledger, _h, protocol=cPickle.HIGHEST_PROTOCOL)
    _h.close()
    return job.fileStore.writeGlobalFile(fn)
