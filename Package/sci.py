#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''
This software has been developed by:

    GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiologia e Historia Forestal)
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politecnica de Madrid
    https://github.com/ggfhf/

Licence: GNU General Public Licence Version 3.
'''

#-------------------------------------------------------------------------------

'''
This source saves the StartCluster information in a file.
'''

#-------------------------------------------------------------------------------

import os
import StringIO
import sys

#-------------------------------------------------------------------------------

def main(argv):
    '''
    Main line of the program.
    '''

    # save the sdtout device
    old_stdout = sys.stdout

    # set the sdtout device as a string buffer
    output = StringIO.StringIO()
    sys.stdout = output

    # show StarCluster package information
    help('starcluster')

    # restore the saved stdout device
    sys.stdout = old_stdout 

    # retrieve the entire contents of the output
    output_string = output.getvalue()

    try:
        if not os.path.exists(os.path.dirname(argv[0])):
            os.makedirs(os.path.dirname(argv[0]))
        with open(argv[0], mode='w') as file_id:
            for line in output_string.splitlines():
                file_id.write( '{0}\n'.format(line.encode('utf8')))
    except Exception as e:
        for line in output_string.splitlines():
            print line

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
