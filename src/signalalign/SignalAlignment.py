from __future__ import print_function

import sys
import os

from sonLib.bioio import fastaWrite
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.bwaWrapper import generateGuideAlignment
from signalalign.utils.fileHandlers import FolderHandler


class SignalAlignment(object):
    def __init__(self,
                 in_fast5,
                 reference_map,
                 path_to_EC_refs,
                 destination,
                 stateMachineType,
                 bwa_index,
                 in_templateHmm,
                 in_complementHmm,
                 in_templateHdp,
                 in_complementHdp,
                 threshold,
                 diagonal_expansion,
                 constraint_trim,
                 degenerate,
                 twoD_chemistry,
                 target_regions=None,
                 output_format="full"):
        self.in_fast5           = in_fast5            # fast5 file to align
        self.reference_map      = reference_map       # map with paths to reference sequences
        self.path_to_EC_refs    = path_to_EC_refs     # place where the reference sequence with ambiguous characters is
        self.destination        = destination         # place where the alignments go, should already exist
        self.stateMachineType   = stateMachineType    # flag for signalMachine
        self.bwa_index          = bwa_index           # index of reference sequence
        self.threshold          = threshold           # min posterior probability to keep
        self.diagonal_expansion = diagonal_expansion  # alignment algorithm param
        self.constraint_trim    = constraint_trim     # alignment algorithm param
        self.target_regions     = target_regions      # only signal-align reads that map to these positions
        self.output_format      = output_format       # smaller output files
        self.degenerate         = degenerate          # set of nucleotides for degenerate characters
        self.twoD_chemistry     = twoD_chemistry      # flag for 2D sequencing runs
        self.temp_folder        = FolderHandler()     # object for holding temporary files (non-toil)
        self.read_name          = self.in_fast5.split("/")[-1][:-6]  # get the name without the '.fast5'

        if (in_templateHmm is not None) and os.path.isfile(in_templateHmm):
            self.in_templateHmm = in_templateHmm
        else:
            self.in_templateHmm = None
        if (in_complementHmm is not None) and os.path.isfile(in_complementHmm):
            self.in_complementHmm = in_complementHmm
        else:
            self.in_complementHmm = None

        # similarly for HDPs
        if (in_templateHdp is not None) and os.path.isfile(in_templateHdp):
            self.in_templateHdp = in_templateHdp
        else:
            self.in_templateHdp = None
        if (in_complementHdp is not None) and os.path.isfile(in_complementHdp):
            self.in_complementHdp = in_complementHdp
        else:
            self.in_complementHdp = None

    def run(self, get_expectations=False):
        print("[SignalAlignment.run]INFO: Starting on {read}".format(read=self.in_fast5), file=sys.stderr)
        if get_expectations:
            assert self.in_templateHmm is not None and self.in_complementHmm is not None,\
                "Need HMM files for model training"
        # file checks
        if os.path.isfile(self.in_fast5) is False:
            print("[SignalAlignment.run]ERROR: Did not find .fast5 at{file}".format(file=self.in_fast5))
            return False

        self.openTempFolder("tempFiles_%s" % self.read_name)
        npRead_ = self.addTempFilePath("temp_%s.npRead" % self.read_name)
        npRead  = NanoporeRead(fast_five_file=self.in_fast5, twoD=self.twoD_chemistry)
        fH      = open(npRead_, "w")
        ok      = npRead.Write(parent_job=None, out_file=fH, initialize=True)
        fH.close()
        if not ok:
            self.failStop("[SignalAlignment.run]File: %s did not pass initial checks" % self.read_name, npRead)
            return False

        read_label    = npRead.read_label  # use this to identify the read throughout
        read_fasta_   = self.addTempFilePath("temp_seq_%s.fa" % read_label)
        temp_samfile_ = self.addTempFilePath("temp_sam_file_%s.sam" % read_label)
        cigar_file_   = self.addTempFilePath("temp_cigar_%s.txt" % read_label)
        if self.twoD_chemistry:
            ok, version, pop1_complement = self.prepare_twod(nanopore_read=npRead, twod_read_path=read_fasta_)
        else:
            ok, version, _ = self.prepare_oned(nanopore_read=npRead, oned_read_path=read_fasta_)
            pop1_complement = None

        # add an indicator for the model being used
        if self.stateMachineType == "threeState":
            model_label = ".sm"
            stateMachineType_flag = ""
        elif self.stateMachineType == "threeStateHdp":
            model_label = ".sm3Hdp"
            stateMachineType_flag = "--sm3Hdp "
            if self.twoD_chemistry:
                assert (self.in_templateHdp is not None) and (self.in_complementHdp is not None), "Need to provide HDPs"
            else:
                assert self.in_templateHdp is not None, "Need to provide Template HDP"
        else:  # make invalid stateMachine control?
            model_label = ".sm"
            stateMachineType_flag = ""

        cigar_string, strand, mapped_refernce = generateGuideAlignment(bwa_index=self.bwa_index,
                                                                       query=read_fasta_,
                                                                       temp_sam_path=temp_samfile_,
                                                                       target_regions=self.target_regions)
        cig_handle = open(cigar_file, "w")
        cig_handle.write(cigar_string + "\n")
        cig_handle.close()
        if mapped_refernce not in self.reference_map.keys():
            if mapped_refernce is False:
                print("[SignalAlignment::run]Read {read} didn't map"
                      "".format(read=read_label), file=sys.stderr)
            else:
                self.failStop("[SignalAlignment::run]Reference {ref} not found in contigs"
                              "{keys}".format(ref=mapped_refernce, keys=self.reference_map.keys()),
                              npRead)
            return False

        # this gives the format: /directory/for/files/file.model.orientation.tsv
        posteriors_file_path = ''

        # forward strand
        if strand == "+":
            if self.output_format == "full":
                posteriors_file_path = self.destination + read_name + model_label + ".forward.tsv"
            elif self.output_format == "variantCaller":
                posteriors_file_path = self.destination + read_name + model_label + ".tsv"
            else:
                posteriors_file_path = self.destination + read_name + model_label + ".assignments"

        # backward strand
        if strand == "-":
            if self.output_format == "full":
                posteriors_file_path = self.destination + read_name + model_label + ".backward.tsv"
            elif self.output_format == "variantCaller":
                posteriors_file_path = self.destination + read_name + model_label + ".tsv"
            else:
                posteriors_file_path = self.destination + read_name + model_label + ".assignments"

        # didn't map
        elif (strand != "+") and (strand != "-"):
            print("[SignalAlignment::run]- {read} gave unrecognizable strand flag: {flag}".format(read=read_label, flag=strand),
                  file=sys.stderr)
            temp_folder.remove_folder()
            return False

        # Alignment/Expectations routine

        # containers and defaults
        path_to_signalAlign = "./signalMachine"

        # flags

        # input (match) models
        if self.in_templateHmm is None:
            self.in_templateHmm = default_template_model_from_version(version=version)
        if self.twoD_chemistry:
            if self.in_complementHmm is None:
                self.in_complementHmm = default_complement_model_from_version(version=version,
                                                                              pop1_complement=pop1_complement)

        assert self.in_templateHmm is not None
        if self.twoD_chemistry:
            assert self.in_complementHmm is not None

        template_model_flag = "-T {} ".format(self.in_templateHmm)
        if self.twoD_chemistry:
            complement_model_flag = "-C {} ".format(self.in_complementHmm)
        else:
            complement_model_flag = ""

        print("signalAlign - NOTICE: template model {t} complement model {c}"
              "".format(t=self.in_templateHmm, c=self.in_complementHmm), file=sys.stderr)

        # reference sequences
        assert self.reference_map[mapped_refernce]["forward"] is not None
        assert self.reference_map[mapped_refernce]["backward"] is not None
        forward_reference = self.reference_map[mapped_refernce]["forward"]
        backward_reference = self.reference_map[mapped_refernce]["backward"]
        assert os.path.isfile(forward_reference)
        assert os.path.isfile(backward_reference)
        forward_ref_flag = "-f {f_ref} ".format(f_ref=forward_reference)
        backward_ref_flag = "-b {b_ref} ".format(b_ref=backward_reference)

        # input HDPs
        if (self.in_templateHdp is not None) or (self.in_complementHdp is not None):
            hdp_flags = "-v {tHdp_loc} ".format(tHdp_loc=self.in_templateHdp)
            if self.twoD_chemistry and self.in_complementHdp is not None:
                hdp_flags += "-w {cHdp_loc} ".format(cHdp_loc=self.in_complementHdp)
        else:
            hdp_flags = ""

        # threshold
        if self.threshold is not None:
            threshold_flag = "-D {threshold} ".format(threshold=self.threshold)
        else:
            threshold_flag = ""

        # diagonal expansion
        if self.diagonal_expansion is not None:
            diag_expansion_flag = "-x {expansion} ".format(expansion=self.diagonal_expansion)
        else:
            diag_expansion_flag = ""

        # constraint trim
        if self.constraint_trim is not None:
            trim_flag = "-m {trim} ".format(trim=self.constraint_trim)
        else:
            trim_flag = ""

        # NOTE: to turn off banded alignment, uncomment this flag, it just trimms away all of the anchors
        #trim_flag = "-m 9999"

        # output format
        fmts = {"full": 0, "variantCaller": 1, "assignments": 2}
        if self.output_format not in fmts.keys():
            temp_folder.remove_folder()
            return False
        out_fmt = "-s {fmt} ".format(fmt=fmts[self.output_format])

        # degenerate nucleotide information
        if self.degenerate is not None:
            degenerate_flag = "-o {} ".format(self.degenerate)
        else:
            degenerate_flag = ""

        if self.twoD_chemistry:
            twoD_flag = "--twoD"
        else:
            twoD_flag = ""
        # commands
        if get_expectations:
            template_expectations_file_path = self.destination + read_name + ".template.expectations"
            complement_expectations_file_path = self.destination + read_name + ".complement.expectations"

            command = \
                "{vA} {td} {degen}{sparse}{model}{f_ref}{b_ref} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} {hdp}-L {readLabel} -p {cigarFile} " \
                "-t {templateExpectations} -c {complementExpectations}"\
                .format(cigar=cigar_string, vA=path_to_signalAlign, model=stateMachineType_flag,
                        f_ref=forward_ref_flag, b_ref=backward_ref_flag, cigarFile=cigar_file,
                        npRead=temp_npRead, readLabel=read_label, td=twoD_flag,
                        templateExpectations=template_expectations_file_path, hdp=hdp_flags,
                        complementExpectations=complement_expectations_file_path, t_model=template_model_flag,
                        c_model=complement_model_flag, thresh=threshold_flag, expansion=diag_expansion_flag,
                        trim=trim_flag, degen=degenerate_flag, sparse=out_fmt)
        else:
            command = \
                "{vA} {td} {degen}{sparse}{model}{f_ref}{b_ref} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} -p {cigarFile} " \
                "-u {posteriors} {hdp}-L {readLabel}"\
                .format(cigar=cigar_string, vA=path_to_signalAlign, model=stateMachineType_flag, sparse=out_fmt,
                        f_ref=forward_ref_flag, b_ref=backward_ref_flag, cigarFile=cigar_file,
                        readLabel=read_label, npRead=temp_npRead, td=twoD_flag,
                        t_model=template_model_flag, c_model=complement_model_flag,
                        posteriors=posteriors_file_path, thresh=threshold_flag, expansion=diag_expansion_flag,
                        trim=trim_flag, hdp=hdp_flags, degen=degenerate_flag)

        # run
        print("signalAlign - running command: ", command, end="\n", file=sys.stderr)
        os.system(command)
        temp_folder.remove_folder()
        return True

    def prepare_oned(self, nanopore_read, oned_read_path):
        try:
            read_file = open(oned_read_path, "w")
            fastaWrite(fileHandleOrFile=read_file,
                       name=nanopore_read.read_label,
                       seq=nanopore_read.template_read)
            version = nanopore_read.version
            read_file.close()
            nanopore_read.close()
            return True, version, False
        except Exception:
            return False, None, False

    def prepare_twod(self, nanopore_read, twod_read_path):
        # check for table to make 'assembled' 2D alignment table fasta with
        if nanopore_read.has2D_alignment_table is False:
            nanopore_read.close()
            return False, None, False
        fasta_handle = open(twod_read_path, "w")
        fastaWrite(fileHandleOrFile=fasta_handle,
                   name=nanopore_read.read_label,
                   seq=nanopore_read.alignment_table_sequence)
        if nanopore_read.complement_model_id == "complement_median68pA_pop1.model":
            pop1_complement = True
        else:
            pop1_complement = False
        version = nanopore_read.version
        fasta_handle.close()
        nanopore_read.close()
        return True, version, pop1_complement

    def openTempFolder(self, temp_dir):
        self.temp_folder.open_folder("%s%s" % (self.destination, temp_dir))

    def addTempFilePath(self, path_to_add):
        return self.temp_folder.add_file_path(path_to_add)

    def failStop(self, message, nanopore_read=None):
        self.temp_folder.remove_folder()
        if nanopore_read is not None:
            nanopore_read.close()
        print(message, file=sys.stderr)
