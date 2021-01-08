# David R Thompson

import os, sys, argparse
import json 
import numpy as np
import pathlib

emit_l1b_datapath = pathlib.Path(__file__).parent.absolute()+'/../data/'

class EMITL1Config():

    def __init__(self, filepath):
        with open(filepath,'r') as fin:
            for name,val in json.load(fin).items():
                setattr(self,name,val)
        self.validate()

        self.rcc, self.rcc_uncert, q = \
            np.loadtxt(self.radiometric_coefficient_file).T

        q, self.wl, self.fwhm = \
            np.loadtxt(self.spectral_calibration_file).T

        self.crf = np.fromfile(self.crf_correction_file, dtype=np.float32, 
              count=(self.columns_raw**2)).reshape((self.columns_raw, 
                                                    self.columns_raw))

        self.srf = np.fromfile(self.srf_correction_file, dtype=np.float32, 
              count=(self.channels_raw**2)).reshape((self.channels_raw, 
                                                     self.channels_raw))

        self.bad = np.fromfile(self.bad_element_file, dtype=np.int16, 
              count=(self.channels_raw**2)).reshape((self.channels_raw, 
                                                     self.columns_raw))

        self.linearity = np.loadtxt(self.linearity_file)

        self.dark = np.fromfile(self.dark_frame_file, dtype=np.float32, 
              count=(self.channels_raw**2)).reshape((self.channels_raw, 
                                                     self.columns_raw))


    def validate(self):
        assert(self.channels_raw == self.channels + \
            (self.channels_masked[0] + self.channels_masked[1]))
        assert(self.channels_total == self.channels_raw + \
            self.channels_header)
        assert(len(self.rcc) == self.channels_raw)
        assert(len(self.wl) == self.channels_raw)
        assert(len(self.linearity) == self.max_DN)



def main():

    description = "Validate a configuration file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_path')
    args = parser.parse_args()
    config = EMITConfig(args.input_path)
    config.validate()
    print('OKAY!')


if __name__ == '__main__':
    main()
