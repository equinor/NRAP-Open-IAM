#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains warranty and acknowledgement for fault ROM.

Author: Ernest N. Lindner
Date: 07/29/2022

Module Name
    flt_warranty

Contents (1)
    show_warranty()

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import textwrap as tr                    # For warranty header statement

# Startup Text - 3 groups of text.
INFORM_0 = """
The software is provided "As Is", without warranty of any kind, express
or implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall
the authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other dealings
in the software.
"""
INFORM_1 = """
This software development was funded by the United States Department
of Energy, National Energy Technology Laboratory, in part, through a site
support contract. Neither the United States Government nor any agency
thereof, nor any of their employees, nor the support contractor, nor any of
their employees, makes any warranty, express or implied, or assumes any
legal liability or responsibility for the accuracy, completeness, or
usefulness of any information, apparatus, product, or process disclosed,
or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service
by trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the
United States Government or any agency thereof. The views and opinions of
authors expressed herein do not necessarily state or reflect those of the
United States Government or any agency thereof.
"""
HEADER_0 = (" "*26) + "Warranty"
HEADER_1 = (" "*25) + "Disclaimer"
HEADER_2 = "    Copyright(c) 2022 by Ernest N. Lindner"
GAP = "  "                      # Gap from left margin
LINER = ("-"*60)                # Divider for printout
TEXT_WIDTH = 60                 # width of text on console


def show_warranty():
    """Write warranty, acknowledgment and disclaimer to console.

    Parameters
    ----------
    N/A

    Returns
    -------
    None
    """
    # Define lists of text to print and wrapper.
    headers = [HEADER_0, HEADER_1]
    informers = [INFORM_0, INFORM_1]
    wrapper = tr.TextWrapper(width=TEXT_WIDTH, fix_sentence_endings=True)

    # Print separator at start.
    print(GAP + LINER)

    # Print three groups of text.
    for indx in range(2):
        title = headers[indx]
        body = informers[indx]

        # Print header of text.
        print("\n", title, "\n")

        # Revise statement to remove leading space and newlines.
        revised = body.replace("\n", " ")
        revised = revised.strip()

        # Setup text at specified width and print each line.
        word_list = wrapper.wrap(text=revised)
        for element in word_list:
            print(GAP + element)

    # Add separator at end, show copyright, and then add last separator.
    print("\n" + GAP + LINER)
    print(GAP + HEADER_2)
    print(GAP + LINER)

    # Return None


#
# -----------------------------------------------------------------------------
# - End of module
