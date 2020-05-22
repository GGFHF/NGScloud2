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
This source contains the general functions and classes used in NGScloud2
software package used in console mode.
'''

#-------------------------------------------------------------------------------

import os
import sys

import xconfiguration
import xlib

#-------------------------------------------------------------------------------

def view_file(file, text):
    '''
    View the contents of a file.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clear_screen()
    print_headers_with_environment(text)

    # read and print the record of file
    print('*' * 20 + '   ' + file + '   ' + '*' * 20)
    try:
        with open(file, mode='r', encoding='iso-8859-1', newline='\n') as file_id:
            print(file_id.read())
    except Exception as e:
        print('*** ERROR: The file {0} can not be read.'.format(file))
        OK = False
    else:
        print('*' * 20 + '*' * (len(file) + 6) + '*' * 20)
        print()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def clear_screen():
    '''
    Clear the screen depending on the Operating System.
    '''

    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        os.system('clear')
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        os.system('cls')

#-------------------------------------------------------------------------------

def print_headers_without_environment(process_name):
    '''
    Print the headers of a screen withtout the environment information.
    '''

    # print the project name, version and the process name
    title = '{0} v {1} - {2}'.format(xlib.get_project_name(), xlib.get_project_version(), process_name)
    line = '-' * len(title)
    print('+-{0}-+'.format(line))
    print('| {0} |'.format(title))
    print('+-{0}-+'.format(line))
    print()

#-------------------------------------------------------------------------------

def print_headers_with_environment(process_name):
    '''
    Print the headers of a screen with environmen information.
    '''

    # print the project name, version and the process name
    title = '{0} v {1} - {2}'.format(xlib.get_project_name(), xlib.get_project_version(), process_name)
    line = '-' * len(title)
    print('+-{0}-+'.format(line))
    print('| {0} |'.format(title))
    print('+-{0}-+'.format(line))
    print()

    # get current region and zone names
    region_name = xconfiguration.get_current_region_name()
    zone_name = xconfiguration.get_current_zone_name()

    # print the environment and the current region and zone names
    print('Environment: {0} - Region: {1} - Zone: {2}'.format(xconfiguration.environment, region_name, zone_name))
    print()

#-------------------------------------------------------------------------------

def confirm_action(action):
    '''
    Ask for the confirmation to do an action.
    '''

    OK = False

    print(action)
    sure = ''
    while sure not in ['y', 'n']:
        sure = input('Are you sure to continue? (y or n): ').lower()
    if sure == 'y':
        OK = True

    return OK

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This source contains the general functions and classes used in {0} software package used in console mode.'.format(xlib.get_project_name()))
    sys.exit(0)

#-------------------------------------------------------------------------------
