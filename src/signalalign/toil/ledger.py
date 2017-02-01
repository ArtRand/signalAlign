import os
import uuid
import shutil
import gzip
import tarfile

from toil_lib import require
from toil_lib.files import tarball_files

from margin.toil.localFileManager import LocalFile, urlDownload, deliverOutput

from signalalign.nanoporeRead import NanoporeRead


def makeNanoporeReadLedgerJobFunction(job, config, sample):
    def check_and_write(filepath):
        require(os.path.exists(filepath), "[makeNanoporeReadLedgerJobFunction]Missing %s" % filepath)
        return job.fileStore.writeGlobalFile(filepath)

    def tar_and_make_send_batch(batch):
        tarname = "%s.tmp" % uuid.uuid4().hex
        tarpath = os.path.join(workdir, tarname)
        tarball_files(tar_name=tarname, file_paths=batch, output_dir=workdir)
        require(os.path.exists(tarpath), "[makeNanoporeReadLedgerJobFunction]Didn't make smaller tar")
        return job.fileStore.writeGlobalFile(tarpath)

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

    member_iter   = [member_paths[i:i + config["batchsize"]] for i in xrange(0, len(member_paths), config["batchsize"])]
    tar_fids      = map(tar_and_make_send_batch, member_iter)
    ledger_shards = [job.addChildJobFn(makeNanoporeReadsJobFunction, fid, config["readstore_dir"], cores=0.5).rv()
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
        return (_l, "%s%s" % (readstore_dir, fn.filenameGetter()))

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
    job.fileStore.logToMaster("REDUCE")
    a = []
    aex = a.extend
    return [aex(i) for i in ledger_shards]
