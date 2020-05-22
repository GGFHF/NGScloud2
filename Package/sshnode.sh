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

#--------------------------------------------------------------------------------

# This script executes remote commands on a cluster node without logging in interactively.

#-------------------------------------------------------------------------------

# -- set -o nounset
# -- set -o verbose
# -- set -o xtrace

#--------------------------------------------------------------------------------

export STARCLUSTER_CONFIG="./config/"$1"-ngscloud2-config.txt"

#-------------------------------------------------------------------------------

function end
{

    exit 0

}

#-------------------------------------------------------------------------------

function manage_error
{

    echo "*** ERROR: StarCluster has ended with return code $1."
    exit 1

}

#--------------------------------------------------------------------------------

starcluster sshnode $2 $3
RC=$?; [ $RC -ne 0 ] && manage_error $RC

end

#-------------------------------------------------------------------------------
