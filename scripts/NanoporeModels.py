from __future__ import print_function
import h5py
import sys
import numpy as np


class NanoporeModel(object):
    def __init__(self, fast5File):
        self.fastFive = h5py.File(fast5File, "r")
        self.stay_prob = 0
        self.skip_prob_bins = []
        self.model_name = ''
        self.model = None

    def export_model(self, destination_path):
        """Exports the model to a file. Format:
        line 1: [correlation coefficient] [level_mean] [level_sd] [noise_mean]
                    [noise_sd] [noise_lambda ] (.../kmer) \n
        line 2: skip bins \n
        line 3: [correlation coefficient] [level_mean] [level_sd, scaled]
                    [noise_mean] [noise_sd] [noise_lambda ] (.../kmer) \n
        """
        def calculate_lambda(noise_mean, noise_stdev):
            return (np.power(noise_mean, 3)) / (np.power(noise_stdev, 2))

        if self.model is None:
            print("This method is meant to be used as part of the child class TemplateModel or ComplementModel",
                  file=sys.stderr)
            # output the model for cPecan to a file
            model_path = destination_path + self.model_name
            out_file = open(model_path, 'w')

            # line 1
            print("0", end=' ', file=out_file) # placeholder for correlation parameter
            for kmer, level_mean, level_stdev, sd_mean, sd_stdev, weight in self.model:
                lam = calculate_lambda(sd_mean, sd_stdev)
                print(level_mean, level_stdev, sd_mean, sd_stdev, lam, end=' ', file=out_file)
            print("", end="\n", file=out_file)
            # line 2
            for _ in self.skip_prob_bins:
                print(_, end=' ', file=out_file)
                print("", end="\n", file=out_file)
            # line 3
            print("0", end=' ', file=out_file) # placeholder for correlation parameter
            for kmer, level_mean, level_stdev, sd_mean, sd_stdev, weight in self.model:
                lam = calculate_lambda(sd_mean, sd_stdev)
                print(level_mean, (level_stdev*1.75), sd_mean, sd_stdev, lam, end=' ', file=out_file)
                print("", end="\n", file=out_file)
        return

    def get_model_dict(self):
        # check
        if self.model is None:
            print("This method is meant to be used as part of the child class TemplateModel or ComplementModel",
                  file=sys.stderr)
            # go through the model and build a lookup table
            model_dict = {}
            for kmer, level_mean, level_stdev, sd_mean, sd_stdev, weight in self.model:
                model_dict[kmer] = [level_mean, level_stdev, sd_mean, sd_stdev]
            return model_dict

    def close(self):
        self.fastFive.close()


class TemplateModel(NanoporeModel):
    def __init__(self, fast5File):
        super(TemplateModel, self).__init__(fast5File=fast5File)
        self.model = self.fastFive['/Analyses/Basecall_2D_000/BaseCalled_template/Model']
        self.stay_prob = np.log2(
            self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_template/Model"].attrs["stay_prob"])
        self.skip_prob_bins = [0.487, 0.412, 0.311, 0.229, 0.174, 0.134, 0.115, 0.103, 0.096, 0.092,
                               0.088, 0.087, 0.084, 0.085, 0.083, 0.082, 0.085, 0.083, 0.084, 0.082,
                               0.080, 0.085, 0.088, 0.086, 0.087, 0.089, 0.085, 0.090, 0.087, 0.096]
        self.parse_model_name()

    def parse_model_name(self):
        model_name = self.fastFive["/Analyses/Basecall_2D_000/Summary/basecall_1d_template"].attrs["model_file"]
        model_name = model_name.split('/')[-1]
        self.model_name = model_name
        return


class ComplementModel(NanoporeModel):
    def __init__(self, fast5File):
        super(ComplementModel, self).__init__(fast5File=fast5File)
        self.model = self.fastFive['/Analyses/Basecall_2D_000/BaseCalled_complement/Model']
        self.stay_prob = np.log2(
            self.fastFive["/Analyses/Basecall_2D_000/BaseCalled_complement/Model"].attrs["stay_prob"])
        self.skip_prob_bins = [0.531, 0.478, 0.405, 0.327, 0.257, 0.207, 0.172, 0.154, 0.138, 0.132,
                               0.127, 0.123, 0.117, 0.115, 0.113, 0.113, 0.115, 0.109, 0.109, 0.107,
                               0.104, 0.105, 0.108, 0.106, 0.111, 0.114, 0.118, 0.119, 0.110, 0.119]
        self.parse_model_name()

    def parse_model_name(self):
        model_name = self.fastFive["/Analyses/Basecall_2D_000/Summary/basecall_1d_complement"].attrs["model_file"]
        model_name = model_name.split('/')[-1]
        self.model_name = model_name
        return
