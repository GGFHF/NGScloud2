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
This file contains functions related to BioInfo applications used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import subprocess
import sys

import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def is_installed_miniconda3(cluster_name, passed_connection, ssh_client):
    '''
    Check if Miniconda3 is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the installation control variable
    is_installed = False

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if not OK:
            for error in error_list:
                error_list.append('{0}\n'.format(error))
                OK = False

    # check the Miniconda3 directory is created
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
            is_installed = True
        else:
            OK = True
            is_installed = False

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and installation control variable
    return (OK, error_list, is_installed)

#-------------------------------------------------------------------------------

def install_miniconda3(cluster_name, log, function=None):
    '''
    Install the Miniconda3 in the cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # warn that the requirements are being verified 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking process requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_miniconda3_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Miniconda3 installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_miniconda3_installation_script()))
        (OK, error_list) = build_miniconda3_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Miniconda3 installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_miniconda3_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_miniconda3_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_miniconda3_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Miniconda3 installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_miniconda3_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_miniconda3_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Miniconda3 installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_miniconda3_installation_starter()))
        (OK, error_list) = build_miniconda3_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Miniconda3 installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_miniconda3_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_miniconda3_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_miniconda3_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Miniconda3 installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_miniconda3_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_miniconda3_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Miniconda3 installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_miniconda3_installation_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_miniconda3_installation_starter()), log)

    # close the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH transport connection ...\n')
        xssh.close_ssh_transport_connection(ssh_transport)
        log.write('The connection is closed.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_miniconda3_installation_script(cluster_name, current_run_dir ):
    '''
    Build the Miniconda3 installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the name, version and download URL of Miniconda3
    (miniconda3_version, miniconda3_url) = xconfiguration.get_bioinfo_app_data(xlib.get_miniconda3_name())

    # write the Miniconda3 installation script
    try:
        if not os.path.exists(os.path.dirname(get_miniconda3_installation_script())):
            os.makedirs(os.path.dirname(get_miniconda3_installation_script()))
        with open(get_miniconda3_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('STATUS_DIR={0}'.format(xlib.get_status_dir(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_OK={0}'.format(xlib.get_status_ok(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_WRONG={0}'.format(xlib.get_status_wrong(current_run_dir))))
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATUS_DIR'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi'))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function remove_miniconda3_directory'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing {0} directory ..."'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    if [ -d "{0}" ]; then'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('        rm -rf {0}'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('        echo "The directory is removed."'))
            script_file_id.write( '{0}\n'.format('    else'))
            script_file_id.write( '{0}\n'.format('        echo "The directory is not found."'))
            script_file_id.write( '{0}\n'.format('    fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function download_miniconda3_package'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Downloading the {0} installation package ..."'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    wget --quiet --output-document {0}.sh {1}'.format(xlib.get_miniconda3_name(), miniconda3_url)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo'))
            script_file_id.write( '{0}\n'.format('    echo "The package is downloaded."'))
            script_file_id.write( '{0}\n'.format('    chmod u+x {0}.sh'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    echo "The run permision is set on."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_miniconda3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing {0} to create Python 3 environment ..."'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    ./{0}.sh -b -p {0}'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error {0} $RC; fi'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    echo "Python 3 environment is created."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function remove_miniconda3_package'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing the {0} installation package ..."'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    rm -f {0}.sh'.format(xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    echo "The software is removed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_gffutils_python3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package gffutils in Python 3 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./pip install --quiet gffutils'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_joblib_python3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package joblib in Python 3 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes joblib'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_matplotlib_python3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package matplotlib in Python 3 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes matplotlib'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_biopython_python3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package biopython in Python 3 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes biopython'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_requests_python3'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package requests in Python 3 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes requests'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function create_python2_environment'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Creating the Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda create --yes --quiet --name py27 python=2.7'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The environment is created."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_gffutils_python2'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package gffutils in Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    source activate py27'))
            script_file_id.write( '{0}\n'.format('    pip install --quiet gffutils'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    conda deactivate'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_joblib_python2'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package joblib in Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    source activate py27'))
            script_file_id.write( '{0}\n'.format('    conda install --quiet --yes joblib'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    conda deactivate'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_matplotlib_python2'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package matplotlib in Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    source activate py27'))
            script_file_id.write( '{0}\n'.format('    conda install --quiet --yes matplotlib'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    conda deactivate'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_biopython_python2'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package biopython in Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    source activate py27'))
            script_file_id.write( '{0}\n'.format('    conda install --quiet --yes biopython'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    conda deactivate'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_requests_python2'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing package requests in Python 2 environment ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    source activate py27'))
            script_file_id.write( '{0}\n'.format('    conda install --quiet --yes requests'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pip $RC; fi'))
            script_file_id.write( '{0}\n'.format('    conda deactivate'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function end'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} installation'.format(xlib.get_miniconda3_name()), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_OK'))
            script_file_id.write( '{0}\n'.format('    exit 0'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function manage_error'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "ERROR: $1 returned error $2"'))
            script_file_id.write( '{0}\n'.format('    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} installation'.format(xlib.get_miniconda3_name()), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_WRONG'))
            script_file_id.write( '{0}\n'.format('    exit 3'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_duration\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
            script_file_id.write( '    HH=`expr $DURATION / 3600`\n')
            script_file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
            script_file_id.write( '    SS=`expr $DURATION % 60`\n')
            script_file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'init\n')
            script_file_id.write( '{0}\n'.format('remove_miniconda3_directory'))
            script_file_id.write( '{0}\n'.format('download_miniconda3_package'))
            script_file_id.write( '{0}\n'.format('install_miniconda3'))
            script_file_id.write( '{0}\n'.format('remove_miniconda3_package'))
            script_file_id.write( '{0}\n'.format('install_gffutils_python3'))
            script_file_id.write( '{0}\n'.format('install_joblib_python3'))
            script_file_id.write( '{0}\n'.format('install_matplotlib_python3'))
            script_file_id.write( '{0}\n'.format('install_biopython_python3'))
            script_file_id.write( '{0}\n'.format('install_requests_python3'))
            script_file_id.write( '{0}\n'.format('create_python2_environment'))
            script_file_id.write( '{0}\n'.format('install_gffutils_python2'))
            script_file_id.write( '{0}\n'.format('install_joblib_python2'))
            script_file_id.write( '{0}\n'.format('install_matplotlib_python2'))
            script_file_id.write( '{0}\n'.format('install_biopython_python2'))
            script_file_id.write( '{0}\n'.format('install_requests_python2'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_miniconda3_installation_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_miniconda3_installation_starter(current_run_dir):
    '''
    Build the starter of the Miniconda3 installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Miniconda3 installation starter
    try:
        if not os.path.exists(os.path.dirname(get_miniconda3_installation_starter())):
            os.makedirs(os.path.dirname(get_miniconda3_installation_starter()))
        with open(get_miniconda3_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_miniconda3_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_miniconda3_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_miniconda3_installation_script():
    '''
    Get the Miniconda3 installation path in the local computer.
    '''

    # assign the Miniconda3 installation path
    miniconda3_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_miniconda3_name())

    # return the Miniconda3 installation path
    return miniconda3_installation_script

#-------------------------------------------------------------------------------

def get_miniconda3_installation_starter():
    '''
    Get the Miniconda3 installation starter path in the local computer.
    '''

    # assign the Miniconda3 installation starter path
    miniconda3_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_miniconda3_name())

    # return the Miniconda3 installation starter path
    return miniconda3_installation_starter

#-------------------------------------------------------------------------------

def is_installed_conda_package(python_version, channel_code, package_code, cluster_name, passed_connection, ssh_client):
    '''
    Check if a Conda package is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the installation control variable
    is_installed = False

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if not OK:
            for error in error_list:
                error_list.append('{0}\n'.format(error))
                OK = False

    # check if the Conda package is installed
    if OK:
        if python_version in [2, 3]:
            if python_version == 3:
                command = 'cd {0}/{1}/bin; ./conda list {2}'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), package_code)
            elif python_version == 2:
                command = 'cd {0}/{1}/bin; source activate py27; ./conda list {2}'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), package_code)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                for line in stdout:
                    if package_code in line:
                        is_installed = True
            else:
                error_list.append('*** ERROR: Wrong command ---> {0}\n'.format(command))
        else:
            error_list.append('Invalid Python version {0}\n'.format(python_version))
            OK = False

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and installation control variable
    return (OK, error_list, is_installed)

#-------------------------------------------------------------------------------

def install_conda_package_list(app_code, app_name, python_version, channel_code, package_code_list, cluster_name, log, function=None):
    '''
    Install a Conda package list.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does no have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # warn that Conda package installation requirements are being verified
    if OK: 
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking the Conda package list ({0}) installation requirements ...\n'.format(str(package_code_list).strip('[]').replace('\'','')))

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the app directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
        else:
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(xlib.get_cluster_app_dir(), cluster_name))
            OK = False

    # check the Miniconda3 installation
    if OK:
        (OK, error_list, is_installed) = is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** error: {0} is not installation\n'.format(xlib.get_miniconda3_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', app_code)
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Conda package installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_conda_package_installation_script()))
        (OK, error_list) = build_conda_package_installation_script(app_name, python_version, channel_code, package_code_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Conda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_conda_package_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_conda_package_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_conda_package_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Conda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_conda_package_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_conda_package_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Conda package installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_conda_package_installation_starter()))
        (OK, error_list) = build_conda_package_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Conda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_conda_package_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_conda_package_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_conda_package_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Conda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_conda_package_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_conda_package_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Conda package installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_conda_package_installation_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_conda_package_installation_starter()), log)

    # close the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH transport connection ...\n')
        xssh.close_ssh_transport_connection(ssh_transport)
        log.write('The connection is closed.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_conda_package_installation_script(app_name, python_version, channel_code, package_code_list, cluster_name, current_run_dir):
    '''
    Build the Conda package installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the name, version and download URL of Miniconda3
    (miniconda3_version, miniconda3_url) = xconfiguration.get_bioinfo_app_data(xlib.get_miniconda3_name())

    # write the Conda package installation script
    try:
        if not os.path.exists(os.path.dirname(get_conda_package_installation_script())):
            os.makedirs(os.path.dirname(get_conda_package_installation_script()))
        with open(get_conda_package_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('STATUS_DIR={0}'.format(xlib.get_status_dir(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_OK={0}'.format(xlib.get_status_ok(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_WRONG={0}'.format(xlib.get_status_wrong(current_run_dir))))
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATUS_DIR'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi'))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            for package_code in package_code_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function install_conda_package_{0}'.format(package_code)))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Installing {0} package {1} ..."'.format(xlib.get_conda_name(), package_code)))
                script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                if python_version == 2:
                    script_file_id.write( '{0}\n'.format('    source activate py27'))
                if channel_code is None:
                    script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes {0}'.format(package_code)))
                else:
                    script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes --channel {0} {1}'.format(channel_code, package_code)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
                if python_version == 2:
                    script_file_id.write( '{0}\n'.format('    conda deactivate'))
                script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function end'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), app_name)))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} installation'.format(app_name), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_OK'))
            script_file_id.write( '{0}\n'.format('    exit 0'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function manage_error'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "ERROR: $1 returned error $2"'))
            script_file_id.write( '{0}\n'.format('    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), app_name)))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} installation'.format(app_name), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_WRONG'))
            script_file_id.write( '{0}\n'.format('    exit 3'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_duration\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
            script_file_id.write( '    HH=`expr $DURATION / 3600`\n')
            script_file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
            script_file_id.write( '    SS=`expr $DURATION % 60`\n')
            script_file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'init\n')
            for package_code in package_code_list:
                script_file_id.write( '{0}\n'.format('install_conda_package_{0}'.format(package_code)))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_conda_package_installation_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_conda_package_installation_starter(current_run_dir):
    '''
    Build the starter of the Conda package installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Conda package installation starter
    try:
        if not os.path.exists(os.path.dirname(get_conda_package_installation_starter())):
            os.makedirs(os.path.dirname(get_conda_package_installation_starter()))
        with open(get_conda_package_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_conda_package_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_conda_package_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_conda_package_installation_script():
    '''
    Get the Conda package installation path in the local computer.
    '''

    # assign the Conda package installation path
    conda_package_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_conda_name())

    # return the Conda package sinstallation path
    return conda_package_installation_script

#-------------------------------------------------------------------------------

def get_conda_package_installation_starter():
    '''
    Get the Conda package installation starter path in the local computer.
    '''

    # assign the Conda package installation starter path
    conda_package_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_conda_name())

    # return the Conda package installation starter path
    return conda_package_installation_starter

#-------------------------------------------------------------------------------

def is_installed_bioconda_package(package_code, cluster_name, passed_connection, ssh_client):
    '''
    Check if a Bioconda package is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the installation control variable
    is_installed = False

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if not OK:
            for error in error_list:
                error_list.append('{0}\n'.format(error))
                OK = False

    # check the Bioconda package directory is created
    if OK:
        command = '[ -d {0}/{1}/envs/{2} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), package_code)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
            is_installed = True
        else:
            OK = True
            is_installed = False

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and installation control variable
    return (OK, error_list, is_installed)

#-------------------------------------------------------------------------------

def install_bioconda_package_list(app_code, app_name, package_list, cluster_name, log, function=None):
    '''
    Install the Bioconda package list in the cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # warn that Bioconda package installation requirements are being verified
    if OK: 
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking the Bioconda package list ({0}) installation requirements ...\n'.format(str(package_list).strip('[]').replace('\'','')))

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the app directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
        else:
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(xlib.get_cluster_app_dir(), cluster_name))
            OK = False

    # check the Miniconda3 installation
    if OK:
        (OK, error_list, is_installed) = is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** error: {0} is not installed.\n'.format(xlib.get_miniconda3_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', app_code)
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Bioconda package installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_bioconda_package_installation_script()))
        (OK, error_list) = build_bioconda_package_installation_script(app_name, package_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Bioconda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_bioconda_package_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_bioconda_package_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Bioconda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Bioconda package installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_bioconda_package_installation_starter()))
        (OK, error_list) = build_bioconda_package_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Bioconda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_bioconda_package_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_bioconda_package_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Bioconda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Bioconda package installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_bioconda_package_installation_starter()), log)

    # close the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH transport connection ...\n')
        xssh.close_ssh_transport_connection(ssh_transport)
        log.write('The connection is closed.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_bioconda_package_installation_script(app_name, package_list, cluster_name, current_run_dir):
    '''
    Build the Bioconda package installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the name, version and download URL of Miniconda3
    (miniconda3_version, miniconda3_url) = xconfiguration.get_bioinfo_app_data(xlib.get_miniconda3_name())

    # write the Bioconda package installation script
    try:
        if not os.path.exists(os.path.dirname(get_bioconda_package_installation_script())):
            os.makedirs(os.path.dirname(get_bioconda_package_installation_script()))
        with open(get_bioconda_package_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('STATUS_DIR={0}'.format(xlib.get_status_dir(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_OK={0}'.format(xlib.get_status_ok(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_WRONG={0}'.format(xlib.get_status_wrong(current_run_dir))))
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATUS_DIR'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi'))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_defaults'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel defaults ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels defaults'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_conda_forge'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel conda-forge ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels conda-forge'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_r'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel r ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels r'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_bioconda'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel bioconda ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels bioconda'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            for package in package_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function remove_bioconda_package_{0}'.format(package[0])))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Removing {0} package {1} ..."'.format(xlib.get_bioconda_name(), package[0])))
                script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('    ./conda env remove --yes --quiet --name {0}'.format(package[0])))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -eq 0 ]; then'))
                script_file_id.write( '{0}\n'.format('      echo "The old package is removed."'))
                script_file_id.write( '{0}\n'.format('    else'))
                script_file_id.write( '{0}\n'.format('      echo "The old package is not found."'))
                script_file_id.write( '{0}\n'.format('    fi'))
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function install_bioconda_package_{0}'.format(package[0])))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                if package[1] == 'last':
                    script_file_id.write( '{0}\n'.format('    echo "Installing {0} package {1} - last version ..."'.format(xlib.get_bioconda_name(), package[0])))
                    script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                    script_file_id.write( '{0}\n'.format('    ./conda create --yes --quiet --name {0} {0}'.format(package[0])))
                    script_file_id.write( '{0}\n'.format('    RC=$?'))
                    script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
                    script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
                else:
                    script_file_id.write( '{0}\n'.format('    echo "Installing {0} package {1} - version {2} ..."'.format(xlib.get_bioconda_name(), package[0], package[1])))
                    script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                    script_file_id.write( '{0}\n'.format('    ./conda create --yes --quiet --name {0} {0}={1}'.format(package[0], package[1])))
                    script_file_id.write( '{0}\n'.format('    RC=$?'))
                    script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then '))
                    script_file_id.write( '{0}\n'.format('        echo "Installing {0} package {1} - last version ..."'.format(xlib.get_bioconda_name(), package[0])))
                    script_file_id.write( '{0}\n'.format('        cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                    script_file_id.write( '{0}\n'.format('        ./conda create --yes --quiet --name {0} {0}'.format(package[0])))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
                    script_file_id.write( '{0}\n'.format('    fi'))
                    script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function end'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), app_name)))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} installation'.format(app_name), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_OK'))
            script_file_id.write( '{0}\n'.format('    exit 0'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function manage_error'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "ERROR: $1 returned error $2"'))
            script_file_id.write( '{0}\n'.format('    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), app_name)))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} installation'.format(app_name), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_WRONG'))
            script_file_id.write( '{0}\n'.format('    exit 3'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_duration\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
            script_file_id.write( '    HH=`expr $DURATION / 3600`\n')
            script_file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
            script_file_id.write( '    SS=`expr $DURATION % 60`\n')
            script_file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'init\n')
            script_file_id.write( '{0}\n'.format('add_channel_defaults'))
            script_file_id.write( '{0}\n'.format('add_channel_conda_forge'))
            script_file_id.write( '{0}\n'.format('add_channel_r'))
            script_file_id.write( '{0}\n'.format('add_channel_bioconda'))
            for package in package_list:
                script_file_id.write( '{0}\n'.format('remove_bioconda_package_{0}'.format(package[0])))
                script_file_id.write( '{0}\n'.format('install_bioconda_package_{0}'.format(package[0])))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_bioconda_package_installation_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_bioconda_package_installation_starter(current_run_dir):
    '''
    Build the starter of the Bioconda package installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Bioconda package installation starter
    try:
        if not os.path.exists(os.path.dirname(get_bioconda_package_installation_starter())):
            os.makedirs(os.path.dirname(get_bioconda_package_installation_starter()))
        with open(get_bioconda_package_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_bioconda_package_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_bioconda_package_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_bioconda_package_installation_script():
    '''
    Get the Bioconda package installation path in the local computer.
    '''

    # assign the Bioconda package installation path
    bioconda_package_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_bioconda_name())

    # return the Bioconda package installation path
    return bioconda_package_installation_script

#-------------------------------------------------------------------------------

def get_bioconda_package_installation_starter():
    '''
    Get the Bioconda package installation starter path in the local computer.
    '''

    # assign the Bioconda package installation starter path
    bioconda_package_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_bioconda_name())

    # return the Bioconda package installation starter path
    return bioconda_package_installation_starter

#-------------------------------------------------------------------------------

def is_installed_r(cluster_name, passed_connection, ssh_client):
    '''
    Check if R is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the installation control variable
    is_installed = False

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if not OK:
            for error in error_list:
                error_list.append('{0}\n'.format(error))
                OK = False

    # check the Bioconda package directory is created
    if OK:
        command = '[ -d {0}/{1}/envs/{2} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_r_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
            is_installed = True
        else:
            OK = True
            is_installed = False

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and installation control variable
    return (OK, error_list, is_installed)

#-------------------------------------------------------------------------------

def install_r(cluster_name, log, function=None):
    '''
    Install the Bioconda package list in the cluster.
    '''

    # initialize the control variable
    OK = True

    # set the addicional R package code list
    package_code_list = []

    # warn that the log window does not have to not be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # warn that R installation requirements are being verified
    if OK: 
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking the R and analysis packages ({0}) installation requirements ...\n'.format(str(package_code_list).strip('[]').replace('\'','')))

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the app directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
        else:
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(xlib.get_cluster_app_dir(), cluster_name))
            OK = False

    # check the Miniconda3 installation
    if OK:
        (OK, error_list, is_installed) = is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** error: {0} is not installed.\n'.format(xlib.get_miniconda3_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_r_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the R installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_r_installation_script()} ...\n')
        (OK, error_list) = build_r_installation_script(package_code_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the R installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_r_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_r_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_r_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the R installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_r_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_r_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the R installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_r_installation_starter()))
        (OK, error_list) = build_r_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the R installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_r_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_r_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_r_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the R installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_r_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_r_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the R installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_r_installation_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_r_installation_starter()), log)

    # close the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH transport connection ...\n')
        xssh.close_ssh_transport_connection(ssh_transport)
        log.write('The connection is closed.\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # execute final function
    if function is not None:
        function()

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_r_installation_script(package_code_list, cluster_name, current_run_dir):
    '''
    Build the R installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the name, version and download URL of Miniconda3
    (miniconda3_version, miniconda3_url) = xconfiguration.get_bioinfo_app_data(xlib.get_miniconda3_name())

    # write the R installation script
    try:
        if not os.path.exists(os.path.dirname(get_r_installation_script())):
            os.makedirs(os.path.dirname(get_r_installation_script()))
        with open(get_r_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('STATUS_DIR={0}'.format(xlib.get_status_dir(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_OK={0}'.format(xlib.get_status_ok(current_run_dir))))
            script_file_id.write( '{0}\n'.format('SCRIPT_STATUS_WRONG={0}'.format(xlib.get_status_wrong(current_run_dir))))
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATUS_DIR'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi'))
            script_file_id.write( '{0}\n'.format('if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi'))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_defaults'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel defaults ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels defaults'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_conda_forge'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel conda-forge ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels conda-forge'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function add_channel_r'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Adding channel r ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda config --add channels r'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The channel is added."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function remove_r'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing {0} ..."'.format(xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda env remove --yes --quiet --name {0}'.format(xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -eq 0 ]; then'))
            script_file_id.write( '{0}\n'.format('      echo "The old package is removed."'))
            script_file_id.write( '{0}\n'.format('    else'))
            script_file_id.write( '{0}\n'.format('      echo "The old package is not found."'))
            script_file_id.write( '{0}\n'.format('    fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function install_r'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Installing {0} ..."'.format(xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('    ./conda create --yes --quiet --name {0} r-essentials'.format(xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
            script_file_id.write( '}\n')
            for package_code in package_code_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function install_r_package_{0}'.format(package_code)))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Installing {0} package {1} ..."'.format(xlib.get_conda_name(), package_code)))
                script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('    source activate {0}'.format(xlib.get_r_name())))
                script_file_id.write( '{0}\n'.format('    ./conda install --quiet --yes {0}'.format(package_code)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
                script_file_id.write( '{0}\n'.format('    conda deactivate'))
                script_file_id.write( '{0}\n'.format('    echo "The package is installed."'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function end'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} and Analysis Packages installation"'.format(xlib.get_project_name(), xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} And Analysis Packages ({1}) installation'.format(xlib.get_r_name(), str(package_code_list).strip('[]').replace('\'','')), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_OK'))
            script_file_id.write( '{0}\n'.format('    exit 0'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function manage_error'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    END_DATETIME=`date --utc +%s`'))
            script_file_id.write( '{0}\n'.format('    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`'))
            script_file_id.write( '{0}\n'.format('    calculate_duration'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "ERROR: $1 returned error $2"'))
            script_file_id.write( '{0}\n'.format('    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} and Analysis Packages installation"'.format(xlib.get_project_name(), xlib.get_r_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} And Analysis Packages ({1}) installation'.format(xlib.get_r_name(), str(package_code_list).strip('[]').replace('\'','')), cluster_name))))
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '{0}\n'.format('    touch $SCRIPT_STATUS_WRONG'))
            script_file_id.write( '{0}\n'.format('    exit 3'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_duration\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    DURATION=`expr $END_DATETIME - $INIT_DATETIME`\n')
            script_file_id.write( '    HH=`expr $DURATION / 3600`\n')
            script_file_id.write( '    MM=`expr $DURATION % 3600 / 60`\n')
            script_file_id.write( '    SS=`expr $DURATION % 60`\n')
            script_file_id.write( '    FORMATTED_DURATION=`printf "%03d:%02d:%02d\\n" $HH $MM $SS`\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'init\n')
            script_file_id.write( '{0}\n'.format('add_channel_defaults'))
            script_file_id.write( '{0}\n'.format('add_channel_conda_forge'))
            script_file_id.write( '{0}\n'.format('add_channel_r'))
            script_file_id.write( '{0}\n'.format('remove_r'))
            script_file_id.write( '{0}\n'.format('install_r'))
            for package_code in package_code_list:
                script_file_id.write( '{0}\n'.format('install_r_package_{0}'.format(package_code)))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_r_installation_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_r_installation_starter(current_run_dir):
    '''
    Build the starter of the R installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the R installation starter
    try:
        if not os.path.exists(os.path.dirname(get_r_installation_starter())):
            os.makedirs(os.path.dirname(get_r_installation_starter()))
        with open(get_r_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_r_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_r_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_r_installation_script():
    '''
    Get the R installation path in the local computer.
    '''

    # assign the R installation path
    r_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_r_name())

    # return the R installation path
    return r_installation_script

#-------------------------------------------------------------------------------

def get_r_installation_starter():
    '''
    Get the R installation starter path in the local computer.
    '''

    # assign the R installation starter path
    r_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_r_name())

    # return the R installation starter path
    return r_installation_starter

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to BioInfo applications used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
