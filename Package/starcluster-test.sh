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

# This script executes actions on a cluster in the test environment.

#-------------------------------------------------------------------------------

# -- set -o nounset
# -- set -o verbose
# -- set -o xtrace

#--------------------------------------------------------------------------------

export STARCLUSTER_CONFIG=./config/test-ngscloud2-config.txt

#-------------------------------------------------------------------------------

function end
{

    exit 0

}

#-------------------------------------------------------------------------------

function manage_error
{

    if [ $1 == 'params' ]; then
        echo 'This scripts needs params.'
    else
        echo "*** ERROR: starcluster has ended with return code $1."
    fi
    exit 1

}

#--------------------------------------------------------------------------------

if [ -z "$*" ]; then
    manage_error params
elif [ $1 == 'sshmaster' ] && [ -n "$3" ]; then
    starcluster $1 $2 "$3"
    RC=$?; [ $RC -ne 0 ] && manage_error $RC
elif [ $1 == 'sshnode' ] && [ -n "$4" ]; then
    starcluster $1 $2 $3 "$4"
    RC=$?; [ $RC -ne 0 ] && manage_error $RC
else
    starcluster $*
    RC=$?; [ $RC -ne 0 ] && manage_error $RC
fi

#-------------------------------------------------------------------------------
