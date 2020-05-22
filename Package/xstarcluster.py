#!/usr/bin/env python3
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
This file contains the functions related to StarCluster used in both console mode and gui mode.
'''
#-------------------------------------------------------------------------------

import re
import sys

import xlib

#-------------------------------------------------------------------------------

def get_starcluster_dir():
    '''
    Get the directory where StarCluster is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the StarCluster directory
    starcluster_dir = ''

    # set the path of the file used to store the StarCluster information
    sci_file_path = '{0}/sci.txt'.format(xlib.get_temp_dir())

    # store the StarCluster information
    if OK:
        command = '{0} {1}'.format(xlib.get_sci(), sci_file_path)
        rc = xlib.run_command(command, sys.stdout)
        if rc != 0:
            error_list.append('*** ERROR: Return code {0} in command -> {1}\n'.format(rc, command))
            OK = False

    # find the StarCluster path
    if OK:
        try:
            file_id = open(sci_file_path, mode='r', encoding='iso-8859-1', newline='\n')
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be opened.\n'.format(sci_file_path))
            OK = False
        else:
            OK = False
            record = file_id.readline()
            while record != '':
                if record.startswith('FILE'):
                    record = file_id.readline().strip()
                    starcluster_dir = record[:record.find('__init__.py')]
                    OK = True
                    break
                record = file_id.readline()

    # return the control variable, error list and StarCluster directory
    return (OK, error_list, starcluster_dir)

#-------------------------------------------------------------------------------

def get_static_py_path():
    '''
    Get the path of the file static.py.
    '''

    # initialize the static.py path
    static_py_path = ''

    # get the StarCluster directory
    (OK, error_list, starcluster_dir) = get_starcluster_dir()

    # set the static.py path
    if OK and starcluster_dir != '':
        static_py_path = starcluster_dir + 'static.py'

    # return the control variable and static.py path
    return (OK, static_py_path)

#-------------------------------------------------------------------------------

def get_instance_type_list_in_static_py():
    '''
    Get the instance type list define in the file static.py.
    '''

    # initialize the instance type list
    instance_type_list = []

    # get the static.py path
    (OK, static_py_path) = get_static_py_path()

    # find the instance types in static.py
    if OK and static_py_path != '':
        try:
            file_id = open(static_py_path, mode='r', encoding='iso-8859-1', newline='\n')
        except Exception as e:
            pass
        else:
            pattern = r'^[ ]+\'(.*)\': \[.*$'
            record = file_id.readline()
            while record != '':
                if record.startswith('INSTANCE_TYPES'):
                    record = file_id.readline()
                    while record != ''  and record.find('}') < 0:
                        try:
                            mo = re.search(pattern, record)
                            instance_type = mo.group(1).strip()
                            instance_type_list.append(instance_type)
                        except Exception as e:
                            pass
                        record = file_id.readline()
                    break
                record = file_id.readline()

    # return the instance type list
    return instance_type_list

#-------------------------------------------------------------------------------

def is_instance_type_in_static_py_path(instance_type):
    '''
    Check if a instance type is defined in the file static.py.
    '''

    # get the instance type list define in the file static.py
    instance_type_list = get_instance_type_list_in_static_py()

    # find the instance type in static.py
    instance_type_found = instance_type in instance_type_list

    # return the control variable
    return instance_type_found

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the functions related to the SSH used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
