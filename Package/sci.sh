#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by:
#
#    GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiologia e Historia Forestal)
#    Dpto. Sistemas y Recursos Naturales
#    ETSI Montes, Forestal y del Medio Natural
#    Universidad Politecnica de Madrid
#    https://github.com/ggfhf/
#
# Licence: GNU General Public Licence Version 3.

#-------------------------------------------------------------------------------

# This script saves the StartCluster information in a file.

#-------------------------------------------------------------------------------

# -- set -o nounset
# -- set -o verbose
# -- set -o xtrace

#--------------------------------------------------------------------------------

PYTHON2PATH=python2

#-------------------------------------------------------------------------------

function end
{

    exit 0

}

#-------------------------------------------------------------------------------

function manage_error
{

    echo "*** ERROR: sci.py has ended with return code $1."
    exit 1

}

#--------------------------------------------------------------------------------

$PYTHON2PATH ./sci.py
RC=$?; [ $RC -ne 0 ] && manage_error $RC

end

#-------------------------------------------------------------------------------
