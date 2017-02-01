import os
import uuid
import tarfile

from toil_lib import require

from margin.toil.localFileManager import LocalFile, urlDownload

from signalalign.nanoporeRead import NanoporeRead


def makeNanoporeReadLedgerJobFunction(job, config, sample):
    def check_and_write(filepath):
        require(os.path.exists(filepath), "[makeNanoporeReadLedgerJobFunction]Missing %s" % filepath)
        return job.fileStore.writeGlobalFile(filepath)
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
    member_paths = (os.path.join(workdir, m.name) for m in members)
    tar_handle.extractall(path=workdir)

    fids          = [check_and_write(f) for f in member_paths]
    batch_iter    = [fids[i:i + config["batchsize"]] for i in xrange(0, len(fids), config["batchsize"])]
    ledger_shards = [job.addChildJobFn(makeNanoporeReadsJobFunction, batch, cores=0.5).rv()
                     for batch in batch_iter]

    return job.addFollowOnJobFn(consolidateLedgerShardsJobFunction, ledger_shards).rv()


def makeNanoporeReadsJobFunction(job, fids):
    def makeNanoporeRead(fid):
        f5    = job.fileStore.readGlobalFile(fid)
        fn    = job.fileStore.getLocalTempFileName()
        fH    = open(fn, "w")
        np    = NanoporeRead(fast_five_file=f5, twoD=False)  # todo make this part of config
        label = np.read_label
        np.write_npRead(fH)
        fH.close()
        require(os.path.exists(fn), "[makeNanoporeReadsJobFunction]Didn't write %s to %s" % (label, fn))
        return (label, job.fileStore.writeGlobalFile(fn))
    job.fileStore.logToMaster("got %s fids to work on" % len(fids))
    ledger_entries = map(makeNanoporeRead, fids)
    job.fileStore.logToMaster("ledger entries %s" % ledger_entries)
    return range(10)


def consolidateLedgerShardsJobFunction(job, ledger_shards):
    job.fileStore.logToMaster("REDUCE")
    a = []
    aex = a.extend
    return [aex(i) for i in ledger_shards]
