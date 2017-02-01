import os
import uuid
import tarfile

from toil_lib import require

from margin.toil.localFileManager import LocalFile, urlDownload


def makeNanoporeReadLedgerJobFunction(job, batchsize, sample):
    def check_and_write(filepath):
        require(os.path.exists(filepath), "[makeNanoporeReadLedgerJobFunction]Missing %s" % filepath)
        return job.fileStore.writeGlobalFile(filepath)
    # download the sample locally
    workdir        = job.fileStore.getLocalTempDir()
    minion_archive = LocalFile(workdir=workdir, filename="%s.tmp" % uuid.uuid4().hex)
    urlDownload(job, sample.minion_data_url, minion_archive)
    require(os.path.exists(minion_archive.fullpathGetter()),
            "[makeNanoporeReadLedgerJobFunction]Didn't download tar from {url} to {path}"
            "".format(url=sample.minion_data_url, path=minion_archive.fullpathGetter()))

    # make batches and send child jobs to make NanoporeReads out of them
    tar_handle   = tarfile.TarFile(minion_archive.fullpathGetter(), "r")
    members      = tar_handle.getmembers()[1:]
    member_paths = (os.path.join(workdir, m.name) for m in members)
    tar_handle.extractall(path=workdir)

    fids          = [check_and_write(f) for f in member_paths]
    batch_iter    = [fids[i:i + batchsize] for i in xrange(0, len(fids), batchsize)]
    ledger_shards = [job.addChildJobFn(makeNanoporeReadsJobFunction, batch).rv()
                     for batch in batch_iter]

    return job.addFollowOnJobFn(consolidateLedgerShardsJobFunction, ledger_shards).rv()


def makeNanoporeReadsJobFunction(job, fids):
    job.fileStore.logToMaster("got %s fids to work on" % len(fids))
    return range(10)


def consolidateLedgerShardsJobFunction(job, ledger_shards):
    job.fileStore.logToMaster("REDUCE")
    a = []
    aex = a.extend
    return [aex(i) for i in ledger_shards]
