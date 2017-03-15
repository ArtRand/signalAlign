from __future__ import print_function

import sys
import h5py

from itertools import islice


TEMPLATE_BASECALL_KEY_0 = "/Analyses/Basecall_1D_000"
TWOD_BASECALL_KEY_0     = "/Analyses/Basecall_2D_000"
VERSION_KEY             = "dragonet version"
SUPPORTED_1D_VERSIONS   = ("1.23.0", "1.22.4")


class NanoporeRead(object):
    def __init__(self, fast_five_file, twoD=False):
        # load the fast5
        self.filename = fast_five_file         # fast5 file path
        self.is_open = self.open()             # bool, is the member .fast5 open?
        self.read_label = ""                   # read label, template read by default
        self.alignment_table_sequence = ""     # the sequence made by assembling the alignment table
        self.template_events = []              # template event sequence
        self.complement_events = []            # complement event sequence
        self.template_read = ""                # template strand base (from fastq) sequence
        self.complement_read = ""              # complement strand base (from fastq) sequence
        self.template_strand_event_map = []    # map of events to kmers in the 1D template read
        self.complement_strand_event_map = []  # map of events to kmers in the 1D complement read
        self.template_event_map = []           # map of template events to kmers in 2D read
        self.complement_event_map = []         # map of complement events to kmers in 2D read
        self.stay_prob = 0                     # TODO, do I need this?
        self.template_model_name = ""          # legacy, for reads base-called by a specific model (template)
        self.complement_model_name = ""        # legacy, for reads base-called by a specific model (complement)
        self.template_scale = 1                # initial values for scaling parameters
        self.template_shift = 1                # 
        self.template_drift = 0                # Template Parameters
        self.template_var = 1                  #
        self.template_scale_sd = 1             #
        self.template_var_sd = 1               # --------------------------------------
        self.complement_scale = 1              #
        self.complement_shift = 1              # Complement Parameters
        self.complement_drift = 0              # 
        self.complement_var = 1                #
        self.complement_scale_sd = 1           #
        self.complement_var_sd = 1             #
        self.twoD = twoD                       # 2D read flag, necessary right now, and the client should know

    def open(self):
        try:
            self.fastFive = h5py.File(self.filename, 'r')
            return True
        except Exception, e:
            self.close()
            print("Error opening file {filename}, {e}".format(filename=self.filename, e=e), file=sys.stderr)
            return False

    def get_latest_basecall_edition(self, address):
        highest = 0
        while(highest < 10):
            if address.format(highest) in self.fastFive:
                highest += 1
                continue
            else:
                return address.format(highest - 1)  # the last base-called version we saw

    def Initialize(self, parent_job):
        if not self.is_open:
            ok = self.open()
            if not ok:
                self.logError("[NanoporeRead:Initialize]ERROR opening %s" % self.filename, parent_job)
                self.close()
                return False
        if self.twoD:
            ok = self._initialize_twoD(parent_job)
        else:
            ok = self._initialize(parent_job)
        return ok

    def _initialize(self, parent_job):
        """Routine setup 1D NanoporeReads, returns false if basecalled with upsupported
        version or is not base-called
        """
        if TEMPLATE_BASECALL_KEY_0 not in self.fastFive:  # not base-called
            self.logError("[NanoporeRead:_initialize]ERROR %s not basecalled" % self.filename, parent_job)
            self.close()
            return False

        oneD_root_address = self.get_latest_basecall_edition("/Analyses/Basecall_1D_00{}")

        if VERSION_KEY not in self.fastFive[oneD_root_address].attrs.keys():
            self.logError("[NanoporeRead:_initialize]ERROR %s missing version" % self.filename, parent_job)
            self.close()
            return False

        self.version = self.fastFive[oneD_root_address].attrs["dragonet version"]

        if self.version not in SUPPORTED_1D_VERSIONS:
            self.logError("[NanoporeRead:_initialize]ERROR %s unsupported version %s " % (self.filename, self.version),
                          parent_job)
            self.close()
            return False

        self.template_event_table_address = "%s/BaseCalled_template/Events" % oneD_root_address
        self.template_model_address       = "%s/BaseCalled_template/Model" % oneD_root_address
        self.template_model_id            = None

        fastq_sequence_address = "%s/BaseCalled_template/Fastq" % oneD_root_address
        if fastq_sequence_address not in self.fastFive:
            self.logError("[NanoporeRead:_initialize]ERROR %s missing fastq" % self.filename, parent_job)
            self.close()
            return False

        self.template_read        = self.fastFive[fastq_sequence_address][()].split()[2]
        self.read_label           = self.fastFive[fastq_sequence_address][()].split()[0][1:]
        self.kmer_length          = len(self.fastFive[self.template_event_table_address][0][4])
        self.template_read_length = len(self.template_read)
        if self.template_read_length <= 0 or not self.read_label or self.kmer_length <= 0:
            self.logError("[NanoporeRead:_initialize]ERROR %s illegal read parameters "
                          "template_read_length: %s, read_label: %s, kmer_length: %s"
                          % (self.template_read_length, self.read_label, self.kmer_length), parent_job)
            self.close()
            return False

        return True

    def _initialize_twoD(self, parent_job=None):
        self.has2D = False
        self.has2D_alignment_table = False

        if TWOD_BASECALL_KEY_0 not in self.fastFive:
            self.close()
            return False

        twoD_address = self.get_latest_basecall_edition("/Analyses/Basecall_2D_00{}")
        if twoD_address not in self.fastFive:
            self.logError("[NanoporeRead::initialize_twoD] Didn't find twoD address, looked here %s " % twoD_address, 
                          parent_job)
            self.close()
            return False

        self.version = self.fastFive[twoD_address].attrs["dragonet version"]

        supported_versions = ["1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4", "1.23.0"]
        if self.version not in supported_versions:
            self.logError("[NanoporeRead::initialize_twoD]Unsupported Version {} (1.15.0, 1.19.0, 1.20.0, "
                          "1.22.2, 1.22.4, 1.23.0 supported)".format(self.version), parent_job)
            self.close()
            return False

        if self.version == "1.15.0":
            oneD_address = self.get_latest_basecall_edition("/Analyses/Basecall_2D_00{}")
        else:
            oneD_address = self.get_latest_basecall_edition("/Analyses/Basecall_1D_00{}")

        twoD_alignment_table_address = twoD_address + "/BaseCalled_2D/Alignment"
        if twoD_alignment_table_address in self.fastFive:
            self.twoD_alignment_table = self.fastFive[twoD_alignment_table_address]
            if len(self.twoD_alignment_table) > 0:
                self.has2D_alignment_table = True
            self.kmer_length = len(self.twoD_alignment_table[0][2])

        twoD_read_sequence_address = twoD_address + "/BaseCalled_2D/Fastq"
        if twoD_read_sequence_address in self.fastFive:
            self.has2D = True
            self.twoD_read_sequence = self.fastFive[twoD_read_sequence_address][()].split()[2]
            self.read_label = self.fastFive[twoD_read_sequence_address][()].split()[0:2][0][1:]

        # initialize version-specific paths
        if self.version == "1.15.0":
            self.template_event_table_address = twoD_address + '/BaseCalled_template/Events'
            self.template_model_address = twoD_address + "/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id(twoD_address + "/Summary/basecall_1d_template")
            self.template_read = self.fastFive[twoD_address + "/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = twoD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = twoD_address + "/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id(twoD_address + "/Summary/basecall_1d_complement")
            self.complement_read = self.fastFive[twoD_address + "/BaseCalled_complement/Fastq"][()].split()[2]
            return True

        elif self.version == "1.19.0" or self.version == "1.20.0":
            self.template_event_table_address = oneD_address + '/BaseCalled_template/Events'
            self.template_model_address = oneD_address + "/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id(oneD_address + "/Summary/basecall_1d_template")
            self.template_read = self.fastFive[oneD_address + "/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = oneD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = oneD_address + "/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id(oneD_address + "/Summary/basecall_1d_complement")
            self.complement_read = self.fastFive[oneD_address + "/BaseCalled_complement/Fastq"][()].split()[2]
            return True

        elif self.version == "1.22.2" or self.version == "1.22.4" or self.version == "1.23.0":
            self.template_event_table_address = oneD_address + '/BaseCalled_template/Events'
            self.template_model_address = ""
            self.template_model_id = None
            self.template_read = self.fastFive[oneD_address + "/BaseCalled_template/Fastq"][()].split()[2]

            self.complement_event_table_address = oneD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = ""
            self.complement_model_id = None
            self.complement_read = self.fastFive[oneD_address + "/BaseCalled_complement/Fastq"][()].split()[2]
            return True
        else:
            self.logError("Unsupported Version (1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4 supported)", parent_job)
            return False

    def assemble_2d_sequence_from_table(self):
        """The 2D read sequence contains kmers that may not map to a template or complement event, which can make
        mapping difficult downstream. This function makes a sequence from the 2D alignment table, which is usually
        pretty similar to the 2D read, except it is guaranteed to have an event map to every position.

        returns: sequence made from alignment table
        """
        def find_kmer_overlap(k_i, k_j):
            """ finds the overlap between two non-identical kmers.
            k_i: one kmer
            k_j: another kmer
            returns: The number of positions not matching
            """
            for i in xrange(1, len(k_i)):
                sk_i = k_i[i:]
                sk_j = k_j[:-i]
                if sk_i == sk_j:
                    return i
            return len(k_i)

        self.alignment_table_sequence = ''
        self.alignment_table_sequence = self.twoD_alignment_table[0][2]
        p_kmer = self.twoD_alignment_table[0][2]

        # iterate through the k-mers in the alignment table
        for t, c, kmer in self.twoD_alignment_table:
            # if we're at a new 6-mer
            if kmer != p_kmer:
                # find overlap, could move up to len(k-mer) - 1 bases
                i = find_kmer_overlap(p_kmer, kmer)
                # append the suffix of the new 6-mer to the sequence
                self.alignment_table_sequence += kmer[-i:]
                # update
                p_kmer = kmer
            else:
                continue
        return

    def init_1d_event_maps(self):
        """Maps the events from the template and complement strands to their base called kmers the map
        generated by this function is called the "strand_event_map" because it only works for mapping the
        strand read (1D read) to to it's events. Uses the same fields as 'get_twoD_event_map' below.
        """
        def make_map(events):
            event_map = [0]
            previous_prob = 0
            for i, line in islice(enumerate(events), 1, None):
                move = line['move']
                this_prob = line['p_model_state']
                if move == 1:
                    event_map.append(i)
                if move > 1:
                    for skip in xrange(move - 1):
                        event_map.append(i - 1)
                    event_map.append(i)
                if move == 0:
                    if this_prob > previous_prob:
                        event_map[-1] = i
                previous_prob = this_prob
            final_event_index = [event_map[-1]]
            padding = final_event_index * (self.kmer_length - 1)
            event_map = event_map + padding
            return event_map

        self.template_strand_event_map = make_map(self.template_events)
        assert len(self.template_strand_event_map) == len(self.template_read)

        if self.twoD:
            self.complement_strand_event_map = make_map(self.complement_events)
            assert len(self.complement_strand_event_map) == len(self.complement_read)

        return True

    def get_twoD_event_map(self):
        """Maps the kmers in the alignment table sequence read to events in the template and complement strand reads
        """
        def kmer_iterator(dna, k):
            for i in xrange(len(dna)):
                kmer = dna[i:(i + k)]
                if len(kmer) == k:
                    yield kmer
        # initialize
        alignment_row = 0
        prev_alignment_kmer = ''
        nb_template_gaps = 0
        previous_complement_event = None
        previous_template_event = None

        #twoD_init = self.initialize_twoD()
        #if twoD_init is False:
        #    return False

        if not self.has2D_alignment_table:
            print("{file} doesn't have 2D alignment table".format(file=self.filename))
            return False

        self.assemble_2d_sequence_from_table()

        # go thought the kmers in the read sequence and match up the events
        for i, seq_kmer in enumerate(kmer_iterator(self.alignment_table_sequence, self.kmer_length)):
            # assign the current row's kmer
            current_alignment_kmer = self.twoD_alignment_table[alignment_row][2]

            # in the situation where there is a repeat kmer in the alignment then
            # we want to pick the best event to kmer alignment, TODO implement this
            # right now we just use the first alignment
            while current_alignment_kmer == prev_alignment_kmer:
                alignment_row += 1
                current_alignment_kmer = self.twoD_alignment_table[alignment_row][2]

            # a match
            if seq_kmer == current_alignment_kmer:
                template_event = self.twoD_alignment_table[alignment_row][0]
                complement_event = self.twoD_alignment_table[alignment_row][1]

                # handle template event
                # if there is a gap, count it and don't add anything to the map
                if template_event == -1:
                    nb_template_gaps += 1

                # if there is an aligned event
                if template_event != -1:
                    # if it is an aligned event and there are no gaps, add it to the map
                    if nb_template_gaps == 0:
                        self.template_event_map.append(template_event)
                        # update
                        previous_template_event = template_event
                    # if there were gaps in the alignment we have to add 'best guess'
                    # event alignments to the map which is the current aligned event
                    if nb_template_gaps > 0:
                        self.template_event_map += [template_event] * (nb_template_gaps + 1)
                        # reset template gaps
                        nb_template_gaps = 0
                        # update
                        previous_template_event = template_event

                # handle complement event
                # if there is a gap, add the last aligned complement event to the map
                if complement_event == -1:
                    self.complement_event_map.append(previous_complement_event)

                # if there is an aligned complement event add it to the map
                if complement_event != -1:
                    self.complement_event_map.append(complement_event)
                    # update the most recent aligned complement event
                    previous_complement_event = complement_event

                # update previous alignment kmer and increment alignment row
                prev_alignment_kmer = current_alignment_kmer
                alignment_row += 1
                continue

            # not a match, meaning that this kmer in the read sequence is not
            # in the event alignment but we need to assign an event to it so
            # we use the heuristic that we use the alignment of the most
            # recent aligned events to this base
            if seq_kmer != current_alignment_kmer:
                self.template_event_map.append(previous_template_event)
                self.complement_event_map.append(previous_complement_event)
                continue

        # fill in the final events for the partial last kmer
        for _ in xrange(self.kmer_length - 1):
            self.template_event_map += [previous_template_event] * (nb_template_gaps + 1)
            self.complement_event_map.append(previous_complement_event)
            nb_template_gaps = 0

        # check that we have mapped all of the bases in the 2D read
        assert(len(self.template_event_map) == len(self.alignment_table_sequence))
        assert(len(self.complement_event_map) == len(self.alignment_table_sequence))
        return True

    def get_template_events(self):
        if self.template_event_table_address in self.fastFive:
            self.template_events = self.fastFive[self.template_event_table_address]
            return True

        if self.template_event_table_address not in self.fastFive:
            return False

    def get_complement_events(self):
        if self.complement_event_table_address in self.fastFive:
            self.complement_events = self.fastFive[self.complement_event_table_address]
            return True

        if self.complement_event_table_address not in self.fastFive:
            return False

    def get_template_model_adjustments(self):
        if self.template_model_address in self.fastFive:
            self.has_template_model = True
            self.template_scale = self.fastFive[self.template_model_address].attrs["scale"]
            self.template_shift = self.fastFive[self.template_model_address].attrs["shift"]
            self.template_drift = self.fastFive[self.template_model_address].attrs["drift"]
            self.template_var = self.fastFive[self.template_model_address].attrs["var"]
            self.template_scale_sd = self.fastFive[self.template_model_address].attrs["scale_sd"]
            self.template_var_sd = self.fastFive[self.template_model_address].attrs["var_sd"]

        if self.template_model_address not in self.fastFive:
            self.has_template_model = False
        return

    def get_complement_model_adjustments(self):
        if self.complement_model_address in self.fastFive:
            self.has_complement_model = True
            self.complement_scale = self.fastFive[self.complement_model_address].attrs["scale"]
            self.complement_shift = self.fastFive[self.complement_model_address].attrs["shift"]
            self.complement_drift = self.fastFive[self.complement_model_address].attrs["drift"]
            self.complement_var = self.fastFive[self.complement_model_address].attrs["var"]
            self.complement_scale_sd = self.fastFive[self.complement_model_address].attrs["scale_sd"]
            self.complement_var_sd = self.fastFive[self.complement_model_address].attrs["var_sd"]

        if self.complement_model_address not in self.fastFive:
            self.has_complement_model = False
        return

    def get_model_id(self, address):
        if address in self.fastFive:
            model_name = self.fastFive[address].attrs["model_file"]
            model_name = model_name.split('/')[-1]
            return model_name
        else:
            return None

    def Write(self, parent_job, out_file, initialize=True):
        if initialize:
            ok = self.Initialize(parent_job)
            if not ok:
                self.close()
                return False

        if self.twoD:
            twoD_map_check          = self.get_twoD_event_map()
            complement_events_check = self.get_complement_events()
        else:
            twoD_map_check          = True
            complement_events_check = True

        template_events_check = self.get_template_events()
        oneD_event_map_check  = self.init_1d_event_maps()

        ok = False not in [twoD_map_check, template_events_check, complement_events_check, oneD_event_map_check]
        if not ok:
            self.close()
            return False

        # get model params
        self.get_template_model_adjustments()
        if self.twoD:
            self.get_complement_model_adjustments()

        # Make the npRead
        # line 1 parameters
        print(len(self.alignment_table_sequence), end=' ', file=out_file)  # 0alignment read length
        print(len(self.template_events), end=' ', file=out_file)           # 1nb of template events
        print(len(self.complement_events), end=' ', file=out_file)         # 2nb of complement events
        print(len(self.template_read), end=' ', file=out_file)             # 3length of template read
        print(len(self.complement_read), end=' ', file=out_file)           # 4length of complement read
        print(self.template_scale, end=' ', file=out_file)                 # 5template scale
        print(self.template_shift, end=' ', file=out_file)                 # 6template shift
        print(self.template_var, end=' ', file=out_file)                   # 7template var
        print(self.template_scale_sd, end=' ', file=out_file)              # 8template scale_sd
        print(self.template_var_sd, end=' ', file=out_file)                # 9template var_sd
        print(self.template_drift, end=' ', file=out_file)                 # 0template_drift
        print(self.complement_scale, end=' ', file=out_file)               # 1complement scale
        print(self.complement_shift, end=' ', file=out_file)               # 2complement shift
        print(self.complement_var, end=' ', file=out_file)                 # 3complement var
        print(self.complement_scale_sd, end=' ', file=out_file)            # 4complement scale_sd
        print(self.complement_var_sd, end=' ', file=out_file)              # 5complement var_sd
        print(self.complement_drift, end=' ', file=out_file)               # 6complement_drift
        print((1 if self.twoD else 0), end='\n', file=out_file)            # has 2D

        # line 2 alignment table sequence
        print(self.alignment_table_sequence, end='\n', file=out_file)

        # line 3 template read
        print(self.template_read, end='\n', file=out_file)

        # line 4 template strand map
        for _ in self.template_strand_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 5 complement read
        print(self.complement_read, end='\n', file=out_file)

        # line 6 complement strand map
        for _ in self.complement_strand_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 7 template 2D event map
        for _ in self.template_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 8 template events
        template_start_time = self.template_events[0]['start']
        for mean, stdev, length, start in self.template_events['mean', 'stdv', 'length', 'start']:
            print(mean, stdev, length, (start - template_start_time), sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 9 complement 2D event map
        for _ in self.complement_event_map[::-1]:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 10 complement events
        if self.twoD:
            complement_start_time = self.complement_events[0]['start']
            for mean, stdev, length, start in self.complement_events['mean', 'stdv', 'length', 'start']:
                print(mean, stdev, length, (start - complement_start_time), sep=' ', end=' ', file=out_file)
        else:
            pass
        print("", end="\n", file=out_file)

        # line 11 model_state (template)
        for _ in self.template_events['model_state']:
            print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 12 p(model) (template)
        for _ in self.template_events['p_model_state']:
            print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 13 model_state (complement)
        if self.twoD:
            for _ in self.complement_events['model_state']:
                print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 14 p(model) (complement)
        if self.twoD:
            for _ in self.complement_events['p_model_state']:
                print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        return True

    def close(self):
        self.fastFive.close()

    @staticmethod
    def logError(message, parent_job=None):
        if parent_job is None:
            print(message, file=sys.stderr)
        else:
            parent_job.fileStore.logToMaster(message)
