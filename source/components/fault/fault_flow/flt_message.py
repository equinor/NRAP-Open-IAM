#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains functions to write miscellaneous information.

Author: Ernest N. Lindner
Date: 08/23/2022

Module Name
    flt_message

Contents (5)
    write_fault_list(fault_controls, fault_list)
    echo_status(alive, phase):
    echo_time_step(desired, alive, current)
    echo_simulation_step(alive, sim_step)
    closeout_time(alive, fault_controls)

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import time                     # For calculating time
from datetime import timedelta  # For calculating runtime

import flt_file as fileop       # For paths and error reports
import flt_config as scfg       # IO directory and file names

# Constants for file handling & writing to Summary File
LEAD_IN = "\n  --> "            # Message start
ZIPX_IN = "  --> "              # Parameter start
INFO_IN = "    > "              # Other messages start


def write_fault_list(fault_controls, fault_list):
    """Printout list of fault existence of simulations to file.

    Parameters
    ----------
    fault_controls = (dict) dictionary of fault parameters
    fault_list = (list) fault existence record at end of simulations

    Returns
    -------
    N/A

    Notes
    -----
    FAULT_EXIST_FILE = Defined destination file name for file record.
    """
    # Construct full path to results directory and summary file.
    #  -- Check for extension "txt".
    sub_path, destination = fileop.get_path_name(scfg.OUTPUT_DIR,
                                                 scfg.FAULT_EXIST_FILE,
                                                 scfg.EXTENSION_TXT)

    # Write  information to file on run.
    try:
        with open(destination, 'w', encoding="utf8") as zum:
            print('\n  Fault_Flo Fault Record', file=zum)

            print(ZIPX_IN +
                  f'Analysis Title = {fault_controls.get("title"):<}.',
                  file=zum)

            # Write time stamp for start of the run.
            now = fault_controls['record_date']
            clump = ZIPX_IN + 'Simulation Run Date: %a %d %b %Y - Time: %H:%M'
            stamp = now.strftime(clump)
            print(stamp, file=zum)

            # Write existence records.
            print('\n  > Existence of Fault in Simulations', file=zum)
            sim = 0
            for list_item in fault_list:
                sim += 1
                zum.write(INFO_IN + f"{sim}. {list_item}\n")

    except OSError as err:
        fileop.io_snag(err, sub_path, scfg.SUMMARY_FILE)

    # return None


def echo_status(alive, phrase):
    """Send progress to console.

    Parameters
    ----------
    alive = (bool) flag to indicate stand-alone operation
    phase = (str) text string indicating position in program

    Returns
    -------
    N/A
    """
    # Echo status to user.
    if alive:
        if phrase == 'BEGIN':
            print(LEAD_IN + 'READING CONTROL FILE.', flush=True)
        elif phrase == 'CLEAN':
            print(LEAD_IN + 'DELETING OLD FILES & PERFORMING INPUT CHECK.',
                  flush=True)
        elif phrase == 'SURFACE':
            print(LEAD_IN + 'DEFINING SURFACE PROFILES.', flush=True)
        elif phrase == 'RUNNING':
            print(LEAD_IN + 'STARTING FAULT_FLO COMPUTATIONS.', flush=True)
        elif phrase == 'SUMMARY':
            print(LEAD_IN + 'CREATED SUMMARY FILE. ', flush=True)
        else:
            print('\n Code Error in Echo!')

    # return None


def echo_time_step(desired, alive, current):
    """Echo time step progress to console.

    Parameters
    ----------
    desired = (bool) flag to indicate header is to be printed
    alive = (bool) stand-alone operation
    current = (float) current time step

    Returns
    -------
    N/A
    """
    # Echo progress line to user on console.
    if desired and alive:
        print(f'    -->> Time Step = {current}', flush=True)

    # return None


def echo_simulation_step(alive, sim_step):
    """Echo simulation progress to console.

    Parameters
    ----------
    alive = (bool) stand-alone operation
    sim_step = (int) current simulation step

    Returns
    -------
    N/A
    """
    # Echo progress line to user on console.
    if alive:
        print(ZIPX_IN + f'Realization Loop = {sim_step + 1}',
              flush=True)

    # return None


def closeout_time(alive, fault_controls):
    """Write analysis time information to console.

    Parameters
    ----------
    alive = (bool) flag to indicate stand-alone operation
    fault_controls = (dict) dictionary of fault parameters

    Returns
    -------
    real_time = simulation time
    """
    # Stop clock on runtime.
    end_code = time.monotonic()

    # Echo end.
    print(LEAD_IN + 'COMPLETED ANALYSIS COMPUTATIONS.', flush=True)

    # Echo execution time.
    start_code = fault_controls['code_start']
    elapsed_time = timedelta(seconds=(end_code - start_code))
    fault_controls['elapsed'] = elapsed_time

    if alive:
        print(INFO_IN + f'Analysis Execution Time = {elapsed_time}')

    return fault_controls


#
# -----------------------------------------------------------------------------
# End of module
