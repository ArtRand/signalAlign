import os
import uuid
import shutil
import gzip
import tarfile
import cPickle

from itertools import chain

from bd2k.util.humanize import human2bytes
from toil_lib import require
from toil_lib.files import tarball_files

from margin.toil.localFileManager import LocalFile, urlDownload, deliverOutput
from signalalign.nanoporeRead import NanoporeRead


def makeReadstoreJobFunction(job, config, samples):
    cores    = config["download_cores"]
    tar_fids = [job.addChildJobFn(prepareFast5Tarfile,
                                  human2bytes(config["split_tars_bigger_than_this"]),
                                  config["put_this_many_reads_in_a_tar"],  # batchsize
                                  config["max_download_slots"],
                                  config["download_part_size"],
                                  sample, cores=cores,
                                  disk=(3 * sample.size)).rv()
                for sample in samples]
    job.addFollowOnJobFn(makeLedgerJobFunction, config, tar_fids)


def prepareFast5Tarfile(job, split_tars_bigger_than_this, batchsize, download_slots, part_size, rs_sample):
    job.fileStore.logToMaster("[prepareFast5Tarfile]Working on sample %s" % rs_sample.sample_label)
    workdir = job.fileStore.getLocalTempDir()
    archive = LocalFile(workdir=workdir, filename="%s.tar" % uuid.uuid4().hex)
    urlDownload(job, rs_sample.URL, archive, download_slots=str(download_slots), part_size=str(part_size))

    _handle = tarfile.open(archive.fullpathGetter(), "r")
    members = _handle.getmembers()[1:]  # the first member is often just the directory with the fast5s
    paths   = [os.path.join(workdir, m.name) for m in members]
    _handle.extractall(path=workdir)

    if rs_sample.size >= split_tars_bigger_than_this:
        _iter    = [paths[i:i + batchsize] for i in xrange(0, len(paths), batchsize)]
        tar_fids = [archiveBatchAndUploadToFileStore(job, b, workdir) for b in _iter]
        _handle.close()
        job.fileStore.logToMaster("[prepareFast5Tarfile]Split %s into %s smaller tars"
                                  % (rs_sample.sample_label, len(tar_fids)))
        return tar_fids
    else:
        tar_fid = archiveBatchAndUploadToFileStore(job, paths, workdir)
        _handle.close()
        return [tar_fid]


def makeLedgerJobFunction(job, config, tar_fids):
    tar_fid_iter = chain(*tar_fids)
    ledger_fids  = [job.addChildJobFn(makeNanoporeReadLedgerJobFunction,
                                      fid,
                                      config["NanoporeRead_batchsize"],
                                      config["readstore_dir"]).rv()
                    for fid in tar_fid_iter]
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


def archiveBatchAndUploadToFileStore(parent_job, batch, workdir):
    tarname = "%s.tmp" % uuid.uuid4().hex
    tarpath = os.path.join(workdir, tarname)
    tarball_files(tar_name=tarname, file_paths=batch, output_dir=workdir)
    require(os.path.exists(tarpath), "[archiveBatchAndUploadToFileStore]Didn't make smaller tar")
    return parent_job.fileStore.writeGlobalFile(tarpath)


def makeNanoporeReadLedgerJobFunction(job, tar_fid, batchsize, readstore_dir):
    workdir        = job.fileStore.getLocalTempDir()
    minion_archive = job.fileStore.readGlobalFile(tar_fid)
    tar_handle     = tarfile.open(minion_archive, "r:gz")
    members        = tar_handle.getmembers()
    member_paths   = [os.path.join(workdir, m.name) for m in members]
    tar_handle.extractall(path=workdir)
    require(batchsize <= len(member_paths),
            "[makeNanoporeReadLedgerJobFunction]Cannot split %s members into batches of %s"
            % (len(member_paths), batchsize))

    member_iter   = [member_paths[i:i + batchsize]
                     for i in xrange(0, len(member_paths), batchsize)]
    tar_fids      = [archiveBatchAndUploadToFileStore(job, b, workdir) for b in member_iter]
    ledger_shards = [job.addChildJobFn(makeNanoporeReadsJobFunction, fid, readstore_dir, cores=0.5).rv()
                     for fid in tar_fids]
    tar_handle.close()
    return job.addFollowOnJobFn(consolidateLedgerShardsJobFunction, ledger_shards).rv()


def makeNanoporeReadsJobFunction(job, tar_fid, readstore_dir):
    def makeNanoporeRead(f5_path):
        # here we load the NanoporeRead and write it to a file
        np = NanoporeRead(fast_five_file=f5_path, twoD=False)  # make this a config arg
        ok = np.Initialize(job)
        if not ok:
            return None
        _l = np.read_label
        tF = job.fileStore.getLocalTempFile()
        fH = open(tF, "w")
        ok = np.Write(job, fH, initialize=False)
        if not ok:
            fH.close()
            return None
        fH.close()
        # then we gzip it and deliver it to the readstore and return the ledger line
        fn = LocalFile(workdir=workdir, filename="%s.np.gz" % _l)
        fH = open(tF, "rb")
        gz = gzip.open(fn.fullpathGetter(), "wb")
        shutil.copyfileobj(fH, gz)
        fH.close()
        gz.close()
        try:
            deliverOutput(job, fn, readstore_dir)
        except RuntimeError:
            job.fileStore.logToMaster("[makeNanoporeReadsJobFunction]Read %s failed to upload" % _l)
            return None
        return (_l, "%s%s\n" % (readstore_dir, fn.filenameGetter()))

    def write_ledger_line(line, fH):
        l = "%s\t%s" % (line[0], line[1])  # read_label, npread URL
        fH.write(l)

    workdir    = job.fileStore.getLocalTempDir()
    tar        = job.fileStore.readGlobalFile(tar_fid)
    tar_handle = tarfile.open(tar, "r:gz")
    members    = tar_handle.getmembers()
    members    = [os.path.join(workdir, m.name) for m in members]
    tar_handle.extractall(path=workdir)
    [require(os.path.exists(m), "[makeNanoporeReadsJobFunction]Missing member %s" % m) for m in members]
    ledger_lines = map(makeNanoporeRead, members)
    tar_handle.close()

    ledger_shard = job.fileStore.getLocalTempFile()
    _handle      = open(ledger_shard, "w")
    [write_ledger_line(l, _handle) for l in ledger_lines if l is not None]
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
