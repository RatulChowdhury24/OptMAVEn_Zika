#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Solvation Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This module contains functions for the use of solvation in the IPRO Suite of
Programs."""

# Include some standard PYTHON modules
import os
import sys
# Also include IPRO Suite Modules
from STANDARDS import *

def SolvationError(IPRO_Error):
    """An error for problems in the Solvation Module."""
    def __init__(self, error = ''):
        """The initialization function of the Solvation Error class."""
        IPRO_Error.__init__(error)

# How to use solvation in CHARMM
#CHARMM_Solvation = """
#GBMV BETA -8 FAST 1 EPSILON 80 DN 1.0 WATR 1.4 GEOM -
# TOL 1e-8 BUFR 0.5 Mem 10 CUTA 20 HSX1 -0.125 HSX2 0.25 -
# ALFRQ 5 EMP 1.5 P4 -4.5 P6 4.0 P3 0.35 ONX 1.9 OFFX 2.1 -
# WTYP 2 NPHI 38 SHIFT -0.102 SLOPE 0.9085 CORR 1 -
# SA 0.00542 SB 0.92
#"""

CHARMM_Solvation = """
GBMV BETA -20 EPSILON 80 DN 1.0 watr 1.4 GEOM -
    TOL 1e-8 BUFR 0.5 Mem 10 CUTA 20 HSX1 -0.125 HSX2 0.25 -
    ALFRQ 1 EMP 0.25 P4 0.0 P6 8.0 P3 0.70 ONX 1.9 OFFX 2.1 -
    WTYP 2 NPHI 38 SHIFT -0.102 SLOPE 0.9085 CORR 1
"""

def get_string(experiment, procedure):
    """Retrieve the string describing how to use solvation in a force field"""
    # First, determine if solvation should be used
    try:
        # If there are specific directions for this procedure
        if procedure.capitalize() + " Solvation" in experiment:
            use = experiment[procedure.capitalize() + " solvation"]
        # Or if there are general directions about how to use solvation
        elif "Use Solvation" in experiment:
            use = experiment["Use Solvation"]
        # Otherwise use the default setting
        else:
            use = defaultUseSolvation
    # If that didn't work, use the default
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        use = defaultUseSolvation
    # Get the force field
    try:
        if "Force Field" in experiment:
            forceField = experiment["Force Field"]
        else:
            forceField = defaultField
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        forceField = defaultField
    # Return the proper string based on the force field
    if not use:
        return ''
    elif forceField == "CHARMM":
        return CHARMM_Solvation
    else:
        text = "The get_string function does not support the " + str(forceField)
        text += " force field."
        raise SolvationError(text)
