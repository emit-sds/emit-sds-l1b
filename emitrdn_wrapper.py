#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A wrapper script to execute radiance calibration (emitrdn.py) followed by destriping (utils/buildflat.py,
utils/medianflat.py, and utils/applyflat.py).

Author: Winston Olson-Duvall, winston.olson-duvall@jpl.nasa.gov
"""

import json
import logging
import os
import shutil
import subprocess
import sys


def set_up_logging(level, log_path):
    # Set up console logging using root logger
    logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", level=level)
    logger = logging.getLogger("emit-sds-l1b")
    # Set up file handler logging
    handler = logging.FileHandler(log_path)
    handler.setLevel(level)
    formatter = logging.Formatter("%(asctime)s %(levelname)s [%(module)s]: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def main():
    """
        This script takes a runconfig.json file containing all required inputs to the scripts below and it executes the
        scripts in this order:
          * emitrdn.py
          * utils/buildflat.py
          * utils/medianflat.py
          * utils/applyflat.py
    """

    in_file = sys.argv[1]

    # Read in runconfig
    print("Reading in runconfig")
    with open(in_file, "r") as f:
        runconfig = json.load(f)

    logger = set_up_logging(runconfig["level"], runconfig["tmp_log_path"])
    logger.info(f"Starting emitrdn_wrapper.py with runconfig {in_file}")

    # Create emitrdn.py command
    cmd = ["python",
           runconfig["emitrdn_exe"],
           f"--mode {runconfig['instrument_mode']}"
           f"--level {runconfig['level']}"
           f"--log_file {runconfig['tmp_log_path']}",
           runconfig["raw_img_path"],
           runconfig["dark_img_path"],
           runconfig["tmp_l1b_config_path"],
           runconfig["tmp_rdn_img_path"],
           runconfig["tmp_rdn_bandmask_img_path"]]
    print("Running emitrdn.py command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

    # Create buildflat.py command
    cmd = ["python",
           runconfig["buildflat_exe"],
           f"--config {runconfig['tmp_l1b_config_path']}",
           runconfig["tmp_rdn_img_path"],
           runconfig["tmp_ffupdate_img_path"]]
    print("Running utils/buildflat.py command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

    # Create medianflat.py command
    cmd = ["python",
           runconfig["medianflat_exe"],
           runconfig["tmp_ff_list_path"],
           runconfig["tmp_ffmedian_img_path"]]
    print("Running utils/medianflat.py command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

    # Create applyflat.py command
    cmd = ["python",
           runconfig["applyflat_exe"],
           runconfig["tmp_rdn_img_path"],
           runconfig["tmp_l1b_config_path"],
           runconfig["tmp_ffmedian_img_path"],
           runconfig["tmp_rdn_destripe_img_path"]]
    print("Running utils/applyflat.py command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)


if __name__ == "__main__":
    main()

