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
This file contains functions related to the NGShelper process used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys
import urllib

import xbioinfoapp
import xcluster
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def is_installed_ngshelper(cluster_name, passed_connection, ssh_client):
    '''
    Check if NGShelper is installed.
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

    # check the NGShelper directory is created
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ngshelper_name())
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

def install_ngshelper(cluster_name, log, function=None):
    '''
    Install the NGShelper software in the cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('Do not close this window, please wait!\n')

    # create the SSH client connection
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
        log.write('Checking installation requirements ...\n')

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
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(xlib.get_cluster_app_dir(), cluster_name))
            OK = False

    # check the Miniconda3 installation
    if OK:
        miniconda3_name = xlib.get_miniconda3_name()
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** error: {0} is not installed.\n'.format(miniconda3_name))
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # initialize the Bioconda package list
    package_code_list = []

    # check the BLAST+ installation
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_blastplus_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('{0} is not installed.\n'.format(xlib.get_blastplus_name()))
                package_code_list.append(xlib.get_blastplus_bioconda_code())
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_ngshelper_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the NGShelper installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_ngshelper_installation_script()))
        (OK, error_list) = build_ngshelper_installation_script(package_code_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the NGShelper installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_ngshelper_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ngshelper_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ngshelper_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the NGShelper installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ngshelper_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ngshelper_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the NGShelper installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ngshelper_installation_starter()))
        (OK, error_list) = build_ngshelper_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the NGShelper installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_ngshelper_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ngshelper_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ngshelper_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the NGShelper installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ngshelper_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ngshelper_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the NGShelper installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ngshelper_installation_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ngshelper_installation_starter()), log)

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # execute final function
    if function is not None:
        function()

    # warn that the log window can be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write(f'{xlib.get_separator()}\n')
        log.write('You can close this window now.\n')

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_ngshelper_installation_script(package_code_list, cluster_name, current_run_dir):
    '''
    Build the NGShelper installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of NGShelper
    (ngshelper_version, ngshelper_url) = xconfiguration.get_bioinfo_app_data(xlib.get_ngshelper_name())

    # write the NGShelper installation script
    try:
        if not os.path.exists(os.path.dirname(get_ngshelper_installation_script())):
            os.makedirs(os.path.dirname(get_ngshelper_installation_script()))
        with open(get_ngshelper_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write(f'PYTHON3_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
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
            script_file_id.write( '{0}\n'.format('function remove_ngshelper_directory'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing {0} directory ..."'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    if [ -d "{0}" ]; then'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('        rm -rf {0}'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('        echo "The directory is removed."'))
            script_file_id.write( '{0}\n'.format('    else'))
            script_file_id.write( '{0}\n'.format('        echo "The directory is not found."'))
            script_file_id.write( '{0}\n'.format('    fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function download_ngshelper_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Downloading the {0} installation file ..."'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            download_script = f'import requests; r = requests.get(\'{ngshelper_url}\') ; open(\'{xlib.get_ngshelper_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $PYTHON3_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is downloaded."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function decompress_ngshelper_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Decompressing the {0} installation file ..."'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write(f'    unzip -u {xlib.get_ngshelper_name()}.zip\n')
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error tar $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is decompressed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function rename_ngshelper_directory'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Renaming the {0} directory ..."'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write(f'    mv {xlib.get_ngshelper_name()}-master {xlib.get_ngshelper_name()}\n')
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error mv $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The directory is renamed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function set_execution_permissions'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Setting execution permissions ..."'))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    chmod u+x {0}/Package/*.py'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The directory is renamed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function remove_ngshelper_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing the {0} installation file ..."'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    rm -f {0}.zip'.format(xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is removed."'))
            script_file_id.write( '}\n')
            if len(package_code_list) > 0:
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
            for package_code in package_code_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function remove_bioconda_package_{0}'.format(package_code)))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Removing {0} package {1} ..."'.format(xlib.get_bioconda_name(), package_code)))
                script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('    ./conda env remove --yes --quiet --name {0}'.format(package_code)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -eq 0 ]; then'))
                script_file_id.write( '{0}\n'.format('      echo "The old package is removed."'))
                script_file_id.write( '{0}\n'.format('    else'))
                script_file_id.write( '{0}\n'.format('      echo "The old package is not found."'))
                script_file_id.write( '{0}\n'.format('    fi'))
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function install_bioconda_package_{0}'.format(package_code)))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Installing {0} package {1} ..."'.format(xlib.get_bioconda_name(), package_code)))
                script_file_id.write( '{0}\n'.format('    cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('    ./conda create --yes --quiet --name {0} {0}'.format(package_code)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error conda $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} installation'.format(xlib.get_ngshelper_name()), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_ngshelper_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} installation'.format(xlib.get_ngshelper_name()), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('remove_ngshelper_directory'))
            script_file_id.write( '{0}\n'.format('download_ngshelper_installation_file'))
            script_file_id.write( '{0}\n'.format('decompress_ngshelper_installation_file'))
            script_file_id.write( '{0}\n'.format('rename_ngshelper_directory'))
            script_file_id.write( '{0}\n'.format('set_execution_permissions'))
            script_file_id.write( '{0}\n'.format('remove_ngshelper_installation_file'))
            if len(package_code_list) > 0:
                script_file_id.write( '{0}\n'.format('add_channel_defaults'))
                script_file_id.write( '{0}\n'.format('add_channel_conda_forge'))
                script_file_id.write( '{0}\n'.format('add_channel_r'))
                script_file_id.write( '{0}\n'.format('add_channel_bioconda'))
            for package_code in package_code_list:
                script_file_id.write( '{0}\n'.format('remove_bioconda_package_{0}'.format(package_code)))
                script_file_id.write( '{0}\n'.format('install_bioconda_package_{0}'.format(package_code)))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ngshelper_installation_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ngshelper_installation_starter(current_run_dir):
    '''
    Build the starter of the NGShelper installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the NGShelper installation starter
    try:
        if not os.path.exists(os.path.dirname(get_ngshelper_installation_starter())):
            os.makedirs(os.path.dirname(get_ngshelper_installation_starter()))
        with open(get_ngshelper_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_ngshelper_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ngshelper_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ngshelper_installation_script():
    '''
    Get the NGShelper installation path in the local computer.
    '''

    # assign the NGShelper installation path
    ngshelper_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_ngshelper_name())

    # return the NGShelper installation path
    return ngshelper_installation_script

#-------------------------------------------------------------------------------

def get_ngshelper_installation_starter():
    '''
    Get the NGShelper installation starter path in the local computer.
    '''

    # assign the NGShelper installation starter path
    ngshelper_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ngshelper_name())

    # return the NGShelper installation starter path
    return ngshelper_installation_starter

#-------------------------------------------------------------------------------

def create_transcript_filter_config_file(experiment_id='exp001', rsem_eval_dataset_id='rsemeval-170101-235959'):
    '''
    Create transcript-filter config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the transcript-filter config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_transcript_filter_config_file())):
            os.makedirs(os.path.dirname(get_transcript_filter_config_file()))
        with open(get_transcript_filter_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The RSEM-EVAL files have to be located in the cluster directory {0}/experiment_id/rsem_eval_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and rsem_eval_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of transcript-filter (NGShelper package) and their meaning in https://github.com/GGFHF/.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('rsem_eval_dataset_id = {0}'.format(rsem_eval_dataset_id), '# rsem_eval dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the transcript-filter parameters'))
            file_id.write( '{0}\n'.format('[transcript-filter parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('minlen = 200', '# transcript with length values less than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format('maxlen = 10000', '# transcript with length values greater than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format('fpkm = 1.0', '# transcript with FPKM values less than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format('tpm = 1.0', '# transcript with TPM values less than this value will be filtered'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_transcript_filter_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_transcript_filter_process(cluster_name, log, function=None):
    '''
    Run a transcript-filter process.
    '''

    # initialize the control variable
    OK = True

    # get the transcript-filter option dictionary
    transcript_filter_option_dict = xlib.get_option_dict(get_transcript_filter_config_file())

    # get the experiment identification
    experiment_id = transcript_filter_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the transcript-filter config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_transcript_filter_name()))
    (OK, error_list) = check_transcript_filter_config_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file is not valid.\n')
        log.write('Please correct this file or recreate the config files.\n')

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

    # check the NGShelper is installed
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ngshelper_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ngshelper_name()))
            OK = False

    # check BLAST+ is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_blastplus_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_blastplus_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_blastplus_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_transcript_filter_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcript-filter process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_transcript_filter_process_script()))
        (OK, error_list) = build_transcript_filter_process_script(cluster_name, current_run_dir, sftp_client)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the transcript-filter process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_transcript_filter_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_transcript_filter_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_transcript_filter_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcript-filter process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcript_filter_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_transcript_filter_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcript-filter process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_transcript_filter_process_starter()))
        (OK, error_list) = build_transcript_filter_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the transcript-filter process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_transcript_filter_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_transcript_filter_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_transcript_filter_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcript-filter process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcript_filter_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_transcript_filter_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the transcript-filter process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcript_filter_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_transcript_filter_process_starter()), log)

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

def check_transcript_filter_config_file(strict):
    '''
    Check the transcript-filter config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        transcript_filter_option_dict = xlib.get_option_dict(get_transcript_filter_config_file())
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in transcript_filter_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = transcript_filter_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "rsem_eval_dataset_id"
            rsem_eval_dataset_id = transcript_filter_option_dict.get('identification', {}).get('rsem_eval_dataset_id', not_found)
            if rsem_eval_dataset_id == not_found:
                error_list.append('*** ERROR: the key "rsem_eval_dataset_id" is not found in the section "identification".')
                OK = False
            elif not rsem_eval_dataset_id.startswith(xlib.get_rsem_eval_code()):
                error_list.append('*** ERROR: the key "rsem_eval_dataset_id" value is not a {0} assembly assessment.'.format(xlib.get_rsem_eval_name()))
                OK = False

        # check section "transcript-filter parameters"
        if 'transcript-filter parameters' not in sections_list:
            error_list.append('*** ERROR: the section "transcript-filter parameters" is not found.')
            OK = False
        else:

            # check section "transcript-filter parameters" - key "minlen"
            minlen = transcript_filter_option_dict.get('transcript-filter parameters', {}).get('minlen', not_found)
            is_ok_minlen = False
            if minlen == not_found:
                error_list.append('*** ERROR: the key "minlen" is not found in the section "transcript-filter parameters".')
                OK = False
            elif not xlib.check_int(minlen, minimum=1):
                error_list.append('*** ERROR: the key "minlen" has to be an integer number greater than or equal to 1.')
                OK = False
            else:
                is_ok_minlen = True

            # check section "transcript-filter parameters" - key "maxlen"
            maxlen = transcript_filter_option_dict.get('transcript-filter parameters', {}).get('maxlen', not_found)
            is_ok_maxlen = False
            if maxlen == not_found:
                error_list.append('*** ERROR: the key "maxlen" is not found in the section "transcript-filter parameters".')
                OK = False
            elif not xlib.check_int(maxlen, minimum=2):
                error_list.append('*** ERROR: the key "maxlen" has to be an integer number greater than or equal to 2.')
                OK = False
            else:
                is_ok_maxlen = True

            # check if maxlen value is greater than or equal than minlen value
            if is_ok_minlen and is_ok_maxlen and int(maxlen) < int(minlen):
                error_list.append('*** ERROR: The value maxlen value ({0}) is less than the minlen value ({1}).'.format(maxlen, minlen))
                OK = False

            # check section "transcript-filter parameters" - key "fpkm"
            fpkm = transcript_filter_option_dict.get('transcript-filter parameters', {}).get('fpkm', not_found)
            if fpkm == not_found:
                error_list.append('*** ERROR: the key "fpkm" is not found in the section "transcript-filter parameters".')
                OK = False
            elif not xlib.check_float(fpkm, minimum=0.):
                error_list.append('*** ERROR: the key "fpkm" has to be a float number greater than or equal to 0.0.')
                OK = False

            # check section "transcript-filter parameters" - key "tpm"
            tpm = transcript_filter_option_dict.get('transcript-filter parameters', {}).get('tpm', not_found)
            if tpm == not_found:
                error_list.append('*** ERROR: the key "tpm" is not found in the section "transcript-filter parameters".')
                OK = False
            elif not xlib.check_float(tpm, minimum=0.):
                error_list.append('*** ERROR: the key "tpm" has to be a float number greater than or equal to 0.0.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_transcript_filter_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transcript_filter_process_script(cluster_name, current_run_dir, sftp_client):
    '''
    Build the current transcript-filter process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the transcript-filter option dictionary
    transcript_filter_option_dict = xlib.get_option_dict(get_transcript_filter_config_file())

    # get the options
    experiment_id = transcript_filter_option_dict['identification']['experiment_id']
    rsem_eval_dataset_id = transcript_filter_option_dict['identification']['rsem_eval_dataset_id']
    minlen = transcript_filter_option_dict['transcript-filter parameters']['minlen']
    maxlen = transcript_filter_option_dict['transcript-filter parameters']['maxlen']
    fpkm = transcript_filter_option_dict['transcript-filter parameters']['fpkm']
    tpm = transcript_filter_option_dict['transcript-filter parameters']['tpm']

    # get the log file name and build local and cluster paths
    local_log_file = '{0}/{1}'.format(xlib.get_temp_dir(), xlib.get_cluster_log_file())
    cluster_log_file = '{0}/{1}/{2}'.format(xlib.get_cluster_experiment_result_dir(experiment_id), rsem_eval_dataset_id, xlib.get_cluster_log_file())

    # download the RSEM-EVAL log file from the cluster
    OK = xssh.get_file(sftp_client, cluster_log_file, local_log_file)
    if not OK:
        error_list.append('*** ERROR: The file {0} does not have been downloaded.'.format(cluster_log_file))

    # get the assembly software, result_data_set_id and assembly type from the RSEM-EVAL log file
    if OK:
        assembly_software = ''
        assembly_dataset_id = ''
        assembly_type = ''
        filtering_data_id = 'FILTERING_DATA'
        pattern = '{0} - (.+): (.+)'.format(filtering_data_id)
        with open(local_log_file, mode='r', encoding='iso-8859-1', newline='\n') as script_file_id:
            for record in script_file_id:
                if record.startswith(filtering_data_id):
                    mo = re.match(pattern, record)
                    data = mo.group(1)
                    value = mo.group(2)
                    if data == 'ASSEMBLY_SOFTWARE':
                        assembly_software = value
                    elif data == 'ASSEMBLY_DATASET_ID':
                        assembly_dataset_id = value
                    elif data == 'ASSEMBLY_TYPE':
                        assembly_type = value
        if assembly_software == '' or assembly_dataset_id == '' or assembly_type == '':
            error_list.append('*** ERROR: Some filtering are not in the file {0}.'.format(cluster_log_file))
            OK = False

    # set the transcriptome file path
    if assembly_software == xlib.get_soapdenovotrans_code():
        if assembly_type == 'CONTIGS':
            transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.contig'
        elif  assembly_type == 'SCAFFOLDS':
            transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.scafSeq'
    elif assembly_software == xlib.get_transabyss_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/transabyss-final.fa'
    elif assembly_software == xlib.get_trinity_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity.fasta'
    elif assembly_software == xlib.get_ggtrinity_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity-GG.fasta'
    elif assembly_software == xlib.get_cd_hit_est_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/clustered-transcriptome.fasta'
    elif assembly_software == xlib.get_transcript_filter_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/filtered-transcriptome.fasta'

    # set the score file path
    if OK:
        score_file = '{0}/{1}.genes.results'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, rsem_eval_dataset_id), rsem_eval_dataset_id)

    # set the output file path
    if OK:
        output_file = '{0}/filtered-transcriptome.fasta'.format(current_run_dir)

    # write the transcript-filter process script
    if OK:
        try:
            if not os.path.exists(os.path.dirname(get_transcript_filter_process_script())):
                os.makedirs(os.path.dirname(get_transcript_filter_process_script()))
            with open(get_transcript_filter_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
                script_file_id.write( '#!/bin/bash\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( '{0}\n'.format('PYTHON3_PATH={0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('NGSHELPER_PATH={0}/{1}/Package'.format(xlib.get_cluster_app_dir(), xlib.get_ngshelper_name())))
                script_file_id.write( '{0}\n'.format('export PATH=$PYTHON3_PATH:$NGSHELPER_PATH:$PATH'))
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
                script_file_id.write( '{0}\n'.format('function run_transcript_filter_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Filtering the transcripts ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        transcript-filter.py \\'))
                script_file_id.write( '{0}\n'.format('            --assembler={0} \\'.format('ngscloud')))
                script_file_id.write( '{0}\n'.format('            --transcriptome={0} \\'.format(transcriptome_file)))
                script_file_id.write( '{0}\n'.format('            --score={0} \\'.format(score_file)))
                script_file_id.write( '{0}\n'.format('            --output={0} \\'.format(output_file)))
                script_file_id.write( '{0}\n'.format('            --minlen={0} \\'.format(minlen)))
                script_file_id.write( '{0}\n'.format('            --maxlen={0} \\'.format(maxlen)))
                script_file_id.write( '{0}\n'.format('            --FPKM={0} \\'.format(fpkm)))
                script_file_id.write( '{0}\n'.format('            --TPM={0} \\'.format(tpm)))
                script_file_id.write( '{0}\n'.format('            --verbose=n'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error transcript-filter.py $RC; fi'))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_transcript_filter_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_transcript_filter_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_transcript_filter_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_transcript_filter_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('run_transcript_filter_process'))
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be created.'.format(get_transcript_filter_process_script()))
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transcript_filter_process_starter(current_run_dir):
    '''
    Build the starter of the current transcript_filter process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the transcript-filter process starter
    try:
        if not os.path.exists(os.path.dirname(get_transcript_filter_process_starter())):
            os.makedirs(os.path.dirname(get_transcript_filter_process_starter()))
        with open(get_transcript_filter_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_transcript_filter_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_transcript_filter_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transcript_filter_config_file():
    '''
    Get the transcript-filter config file path.
    '''

    # assign the transcript-filter config file path
    transcript_filter_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_transcript_filter_code())

    # return the transcript-filter config file path
    return transcript_filter_config_file

#-------------------------------------------------------------------------------

def get_transcript_filter_process_script():
    '''
    Get the transcript-filter process script path in the local computer.
    '''

    # assign the transcript-filter script path
    transcript_filter_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_transcript_filter_code())

    # return the transcript-filter script path
    return transcript_filter_process_script

#-------------------------------------------------------------------------------

def get_transcript_filter_process_starter():
    '''
    Get the transcript-filter process starter path in the local computer.
    '''

    # assign the transcript-filter process starter path
    transcript_filter_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_transcript_filter_code())

    # return the transcript-filter starter path
    return transcript_filter_process_starter

#-------------------------------------------------------------------------------

def create_transcriptome_blastx_config_file(database_dataset_id='RefSeq_Plant_Protein', protein_database_name='RefSeq_Plant_Protein', experiment_id='exp001', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create transcriptome-blastx config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the assembly software
    if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()):
        assembly_software = xlib.get_soapdenovotrans_code()
    elif assembly_dataset_id.startswith(xlib.get_transabyss_code()):
        assembly_software = xlib.get_transabyss_code()
    elif assembly_dataset_id.startswith(xlib.get_trinity_code()):
        assembly_software = xlib.get_trinity_code()
    elif assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()):
        assembly_software = xlib.get_cd_hit_est_code()
    elif assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
        assembly_software = xlib.get_transcript_filter_code()

    # create the transcript-filter config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_transcriptome_blastx_config_file())):
            os.makedirs(os.path.dirname(get_transcriptome_blastx_config_file()))
        with open(get_transcriptome_blastx_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The database files have to be located in the cluster directory {0}/database_dataset_id'.format(xlib.get_cluster_database_dir())))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, database_dataset_id and assembly_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of transcriptome-blastx (NGShelper package) and their meaning in https://github.com/GGFHF/.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('database_dataset_id = {0}'.format(database_dataset_id), '# database dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('protein_database_name = {0}'.format(protein_database_name), '# protein database name'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_type = {0}'.format(assembly_type), '# assembly type: CONTIGS or SCAFFOLDS in {0}; NONE in any other case'.format(xlib.get_soapdenovotrans_name())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the transcriptome-blastx parameters'))
            file_id.write( '{0}\n'.format('[transcriptome-blastx parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('node_number = 1', '# node number (previously they have to be started)'))
            file_id.write( '{0:<50} {1}\n'.format('blastx_thread_number = 1', '# threads number using by blastx in every node'))
            file_id.write( '{0:<50} {1}\n'.format('e_value = 1E-6', '# expectation value (E-value) threshold for saving hits'))
            file_id.write( '{0:<50} {1}\n'.format('max_target_seqs = 10', '# maximum number of aligned sequences to keep'))
            file_id.write( '{0:<50} {1}\n'.format('max_hsps = 999999', '# maximum number of HSPs per subject sequence to save for each query'))
            file_id.write( '{0:<50} {1}\n'.format('qcov_hsp_perc = 0.0', '# alignments below the specified query coverage per HSPs are removed'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_transcriptome_blastx_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_transcriptome_blastx_process(cluster_name, log, function=None):
    '''
    Run a transcriptome-blastx process.
    '''

    # initialize the control variable
    OK = True

    # get the transcriptome-blastx option dictionary
    transcriptome_blastx_option_dict = xlib.get_option_dict(get_transcriptome_blastx_config_file())

    # get the experiment identification
    experiment_id = transcriptome_blastx_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the transcriptome-blastx config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_transcriptome_blastx_name()))
    (OK, error_list) = check_transcriptome_blastx_config_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file is not valid.\n')
        log.write('Please correct this file or recreate the config files.\n')

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

    # check the node number
    if OK:

        # get the running node number
        running_node_number = len(xec2.get_cluster_node_list(cluster_name)) - 1

        # get the transcriptome-blastx option dictionary
        transcriptome_blastx_option_dict = xlib.get_option_dict(get_transcriptome_blastx_config_file())

        # get the requested nodes number
        requested_node_number = int(transcriptome_blastx_option_dict['transcriptome-blastx parameters']['node_number'])

        # check the requested nodes number is less than or equal than the running node number
        if requested_node_number > running_node_number:
            log.write('*** ERROR: The requested node number ({0}) is greater than running node number ({1}).\n'.format(requested_node_number,running_node_number))
            OK = False

    # check the NGShelper is installed
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ngshelper_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ngshelper_name()))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_transcriptome_blastx_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcriptome-blastx process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_transcriptome_blastx_process_script()))
        (OK, error_list) = build_transcriptome_blastx_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the transcriptome-blastx process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_transcriptome_blastx_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_transcriptome_blastx_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcriptome-blastx process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcriptome-blastx process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_transcriptome_blastx_process_starter()))
        (OK, error_list) = build_transcriptome_blastx_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the transcriptome-blastx process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_transcriptome_blastx_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_transcriptome_blastx_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcriptome-blastx process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the transcriptome-blastx process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_transcriptome_blastx_process_starter()), log)

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

def check_transcriptome_blastx_config_file(strict):
    '''
    Check the transcriptome-blastx config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        transcriptome_blastx_option_dict = xlib.get_option_dict(get_transcriptome_blastx_config_file())
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in transcriptome_blastx_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = transcriptome_blastx_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "database_dataset_id"
            database_dataset_id = transcriptome_blastx_option_dict.get('identification', {}).get('database_dataset_id', not_found)
            if database_dataset_id == not_found:
                error_list.append('*** ERROR: the key "database_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "protein_database_name"
            protein_database_name = transcriptome_blastx_option_dict.get('identification', {}).get('protein_database_name', not_found)
            if protein_database_name == not_found:
                error_list.append('*** ERROR: the key "protein_database_name" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = transcriptome_blastx_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "assembly_software" has to be {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = transcriptome_blastx_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "assembly_dataset_id" has to start with {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = transcriptome_blastx_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {0} or NONE in any other case.'.format(xlib.get_soapdenovotrans_name()))
                    OK = False

        # check section "transcriptome-blastx parameters"
        if 'transcriptome-blastx parameters' not in sections_list:
            error_list.append('*** ERROR: the section "transcriptome-blastx parameters" is not found.')
            OK = False
        else:

            # check section "transcriptome-blastx parameters" - key "node_number"
            node_number = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('node_number', not_found)
            if node_number == not_found:
                error_list.append('*** ERROR: the key "node_number" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_int(node_number, minimum=1):
                error_list.append('*** ERROR: the key "node_number" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "transcriptome-blastx parameters" - key "blastx_thread_number"
            blastx_thread_number = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('blastx_thread_number', not_found)
            if blastx_thread_number == not_found:
                error_list.append('*** ERROR: the key "blastx_thread_number" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_int(blastx_thread_number, minimum=1):
                error_list.append('*** ERROR: the key "blastx_thread_number" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "transcriptome-blastx parameters" - key "e_value"
            e_value = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('e_value', not_found)
            if e_value == not_found:
                error_list.append('*** ERROR: the key "e_value" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_float(e_value, minimum=0., mne=1E-12):
                error_list.append('*** ERROR: the key "e_value" has to be a float number greater than to 0.0.')
                OK = False

            # check section "transcriptome-blastx parameters" - key "max_target_seqs"
            max_target_seqs = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('max_target_seqs', not_found)
            if max_target_seqs == not_found:
                error_list.append('*** ERROR: the key "max_target_seqs" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_int(max_target_seqs, minimum=1):
                error_list.append('*** ERROR: the key "max_target_seqsr" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "transcriptome-blastx parameters" - key "max_hsps"
            max_hsps = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('max_hsps', not_found)
            if max_hsps == not_found:
                error_list.append('*** ERROR: the key "max_hsps" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_int(max_hsps, minimum=1):
                error_list.append('*** ERROR: the key "max_hsps" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "transcriptome-blastx parameters" - key "qcov_hsp_perc"
            qcov_hsp_perc = transcriptome_blastx_option_dict.get('transcriptome-blastx parameters', {}).get('qcov_hsp_perc', not_found)
            if qcov_hsp_perc == not_found:
                error_list.append('*** ERROR: the key "qcov_hsp_perc" is not found in the section "transcriptome-blastx parameters".')
                OK = False
            elif not xlib.check_float(qcov_hsp_perc, minimum=0., maximum=100., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "qcov_hsp_perc" has to be a float number greater than or equal to 0.0 and less than 100.0.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_transcriptome_blastx_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transcriptome_blastx_process_script(cluster_name, current_run_dir):
    '''
    Build the current transcriptome-blastx process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the transcriptome-blastx option dictionary
    transcriptome_blastx_option_dict = xlib.get_option_dict(get_transcriptome_blastx_config_file())

    # get the options
    database_dataset_id = transcriptome_blastx_option_dict['identification']['database_dataset_id']
    protein_database_name = transcriptome_blastx_option_dict['identification']['protein_database_name']
    experiment_id = transcriptome_blastx_option_dict['identification']['experiment_id']
    assembly_software = transcriptome_blastx_option_dict['identification']['assembly_software']
    assembly_dataset_id = transcriptome_blastx_option_dict['identification']['assembly_dataset_id']
    assembly_type = transcriptome_blastx_option_dict['identification']['assembly_type']
    node_number = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['node_number']
    blastx_thread_number = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['blastx_thread_number']
    e_value = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['e_value']
    max_target_seqs = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['max_target_seqs']
    max_hsps = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['max_hsps']
    qcov_hsp_perc = transcriptome_blastx_option_dict['transcriptome-blastx parameters']['qcov_hsp_perc']

    # set the transcriptome file path
    if assembly_software == xlib.get_soapdenovotrans_code():
        if assembly_type == 'CONTIGS':
            transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.contig'
        elif  assembly_type == 'SCAFFOLDS':
            transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.scafSeq'
    elif assembly_software == xlib.get_transabyss_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/transabyss-final.fa'
    elif assembly_software == xlib.get_trinity_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity.fasta'
    elif assembly_software == xlib.get_ggtrinity_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity-GG.fasta'
    elif assembly_software == xlib.get_cd_hit_est_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/clustered-transcriptome.fasta'
    elif assembly_software == xlib.get_transcript_filter_code():
        transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/filtered-transcriptome.fasta'

    # write the transcriptome-blastx process script
    if OK:
        try:
            if not os.path.exists(os.path.dirname(get_transcriptome_blastx_process_script())):
                os.makedirs(os.path.dirname(get_transcriptome_blastx_process_script()))
            with open(get_transcriptome_blastx_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
                script_file_id.write( '#!/bin/bash\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( '{0}\n'.format('PYTHON3_PATH={0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('NGSHELPER_PATH={0}/{1}/Package'.format(xlib.get_cluster_app_dir(), xlib.get_ngshelper_name())))
                script_file_id.write( '{0}\n'.format('BLASTPLUS_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_blastplus_bioconda_code())))
                script_file_id.write( '{0}\n'.format('export PATH=$PYTHON3_PATH:$NGSHELPER_PATH:$BLASTPLUS_PATH:$PATH'))
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
                script_file_id.write( '{0}\n'.format('function run_transcriptome_blastx_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Running the transcriptome blastx process ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        transcriptome-blastx.py \\'))
                script_file_id.write( '{0}\n'.format('            --machine_type="ngscloud" \\'))
                script_file_id.write( '{0}\n'.format('            --node_number={0} \\'.format(node_number)))
                script_file_id.write( '{0}\n'.format('            --blastx_thread_number={0} \\'.format(blastx_thread_number)))
                script_file_id.write( '{0}\n'.format('            --blast_db={0} \\'.format(xlib.get_cluster_database_dataset_dir(database_dataset_id))))
                script_file_id.write( '{0}\n'.format('            --protein_database_name={0} \\'.format(protein_database_name)))
                script_file_id.write( '{0}\n'.format('            --transcriptome={0} \\'.format(transcriptome_file)))
                script_file_id.write( '{0}\n'.format('            --e_value={0} \\'.format(e_value)))
                script_file_id.write( '{0}\n'.format('            --max_target_seqs={0} \\'.format(max_target_seqs)))
                script_file_id.write( '{0}\n'.format('            --max_hsps={0} \\'.format(max_hsps)))
                script_file_id.write( '{0}\n'.format('            --qcov_hsp_perc={0} \\'.format(qcov_hsp_perc)))
                script_file_id.write( '{0}\n'.format('            --output={0} \\'.format(current_run_dir)))
                script_file_id.write( '{0}\n'.format('            --email={0} \\'.format(xconfiguration.get_contact_data())))
                script_file_id.write( '{0}\n'.format('            --verbose="n"'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error transcriptome-blastx.py $RC; fi'))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_transcriptome_blastx_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_transcript_filter_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_transcriptome_blastx_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_transcript_filter_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('run_transcriptome_blastx_process'))
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be created.'.format(get_transcriptome_blastx_process_script()))
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transcriptome_blastx_process_starter(current_run_dir):
    '''
    Build the starter of the current transcriptome-blastx process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the transcriptome-blastx process starter
    try:
        if not os.path.exists(os.path.dirname(get_transcriptome_blastx_process_starter())):
            os.makedirs(os.path.dirname(get_transcriptome_blastx_process_starter()))
        with open(get_transcriptome_blastx_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_transcriptome_blastx_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_transcriptome_blastx_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_config_file():
    '''
    Get the transcriptome-blastx config file path.
    '''

    # assign the transcriptome-blastx config file path
    transcriptome_blastx_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_transcriptome_blastx_code())

    # return the transcriptome-blastx config file path
    return transcriptome_blastx_config_file

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_process_script():
    '''
    Get the transcriptome-blastx process script path in the local computer.
    '''

    # assign the transcriptome-blastx script path
    transcriptome_blastx_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_transcriptome_blastx_code())

    # return the transcriptome-blastx script path
    return transcriptome_blastx_process_script

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_process_starter():
    '''
    Get the transcriptome-blastx process starter path in the local computer.
    '''

    # assign the transcriptome-blastx process starter path
    transcriptome_blastx_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_transcriptome_blastx_code())

    # return the transcriptome-blastx starter path
    return transcriptome_blastx_process_starter

#-------------------------------------------------------------------------------
    
def get_assembly_software_code_list():
    '''
    Get the code list of "assembly_software".
    '''

    return [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(),  xlib.get_transcript_filter_code()]

#-------------------------------------------------------------------------------
    
def get_assembly_software_code_list_text():
    '''
    Get the code list of "assembly_software" as text.
    '''

    return '{0} ({1}) or {2} ({3}) or {4} ({5}) or {6} ({7}) or {8} ({9}) or {10} ({11})'.format(xlib.get_soapdenovotrans_code(), xlib.get_soapdenovotrans_name(), xlib.get_transabyss_code(), xlib.get_transabyss_name(), xlib.get_trinity_code(), xlib.get_trinity_name(), xlib.get_ggtrinity_code(), xlib.get_ggtrinity_name(), xlib.get_cd_hit_est_code(), xlib.get_cd_hit_est_name(), xlib.get_transcript_filter_code(), xlib.get_transcript_filter_name())

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the NGShelper process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
