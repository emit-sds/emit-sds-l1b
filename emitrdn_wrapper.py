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
    with open(in_file, "r") as f:
        runconfig = json.load(f)

    # output dir should be created by workflow manager, but if not, create it here
    output_dir = f"{runconfig['tmp_dir']}/output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create path names
    rdn_basename = os.path.basename(runconfig["raw_img_path"]).replace("l1a", "l1b").replace("raw", "rdn")
    log_path = f"{output_dir}/{rdn_basename.replace('.img', '_pge.log')}"
    rdn_img_path = f"{output_dir}/{rdn_basename}"
    bandmask_img_path = f"{output_dir}/{rdn_basename.replace('rdn', 'bandmask')}"
    ffupdate_img_path = f"{output_dir}/{rdn_basename.replace('rdn', 'ffupdate')}"
    ffmedian_img_path = f"{output_dir}/{rdn_basename.replace('rdn', 'ffmedian')}"
    rdn_destripe_img_path = f"{output_dir}/{rdn_basename.replace('rdn', 'rdndestripe')}"
    utils_path = f"{runconfig['repository_dir']}/utils"
    emitrdn_exe = f"{runconfig['repository_dir']}/emitrdn.py"
    buildflat_exe = f"{runconfig['repository_dir']}/utils/buildflat.py"
    medianflat_exe = f"{runconfig['repository_dir']}/utils/medianflat.py"
    applyflat_exe = f"{runconfig['repository_dir']}/utils/applyflat.py"

    logger = set_up_logging(runconfig["level"], log_path)
    logger.info(f"Starting emitrdn_wrapper.py with runconfig {in_file}")

    # Write out tmp_l1b_config and tmp_ff_list
    l1b_config_path = f"{runconfig['tmp_dir']}/l1b_config.json"
    with open(l1b_config_path, "w") as f:
        json.dump(runconfig["l1b_config"], f, indent=4)

    # Set environment variables
    env = os.environ.copy()
    env["PYTHONPATH"] = f"$PYTHONPATH:{utils_path}"
    env["RAY_worker_register_timeout_seconds"] = "600"

    # Create emitrdn.py command
    cmd = ["python",
           emitrdn_exe,
           f"--mode {runconfig['instrument_mode']}",
           f"--level {runconfig['level']}",
           f"--log_file {log_path}",
           runconfig["raw_img_path"],
           runconfig["dark_img_path"],
           l1b_config_path,
           rdn_img_path,
           bandmask_img_path]
    logger.info("Running emitrdn.py command: " + " ".join(cmd))
    output = subprocess.run(" ".join(cmd), shell=True, env=env)
    if output.returncode == 0:
        logger.info("emitrdn.py command SUCCEEDED")
    else:
        logger.error("emitrdn.py command FAILED, exiting...")
        raise RuntimeError(output.stderr.decode("utf-8"))

    # Create buildflat.py command
    cmd = ["python",
           buildflat_exe,
           f"--config {l1b_config_path}",
           rdn_img_path,
           ffupdate_img_path]
    logger.info("Running utils/buildflat.py command: " + " ".join(cmd))
    output = subprocess.run(" ".join(cmd), shell=True, env=env)
    if output.returncode == 0:
        logger.info("buildflat.py command SUCCEEDED")
    else:
        logger.error("buildflat.py command FAILED, but processing will continue and there will be no flat field "
                     "update for this scene.  Removing flat field if it exists.")
        if os.path.exists(ffupdate_img_path):
            os.remove(ffupdate_img_path)

    # Only proceed if we have at least one flat field update path.
    # NOTE: The SDS workflow code enforces a minimum of 100 paths
    if len(runconfig["flat_field_update_paths"]) > 0:
        # Create ff_list
        ff_list_path = f"{runconfig['tmp_dir']}/ff_list.txt"
        with open(ff_list_path, "w") as f:
            for p in runconfig["flat_field_update_paths"]:
                f.write(f"{p}\n")

        # Create medianflat.py command
        cmd = ["python",
               medianflat_exe,
               ff_list_path,
               ffmedian_img_path]
        logger.info("Running utils/medianflat.py command: " + " ".join(cmd))
        output = subprocess.run(" ".join(cmd), shell=True, env=env)
        if output.returncode == 0:
            logger.info("medianflat.py command SUCCEEDED")
        else:
            logger.error("medianflat.py command FAILED, exiting...")
            raise RuntimeError(output.stderr.decode("utf-8"))

        # Create applyflat.py command
        cmd = ["python",
               applyflat_exe,
               rdn_img_path,
               l1b_config_path,
               ffmedian_img_path,
               rdn_destripe_img_path]
        logger.info("Running utils/applyflat.py command: " + " ".join(cmd))
        output = subprocess.run(" ".join(cmd), shell=True, env=env)
        if output.returncode == 0:
            logger.info("applyflat.py command SUCCEEDED")
        else:
            logger.error("applyflat.py command FAILED, exiting...")
            raise RuntimeError(output.stderr.decode("utf-8"))


if __name__ == "__main__":
    main()
