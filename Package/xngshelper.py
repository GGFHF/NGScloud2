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
                error_list.append(f'{error}\n')
                OK = False

    # check the NGShelper directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ngshelper_name()} ] && echo RC=0 || echo RC=1'
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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the app directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write(f'You have to link a volume in the mounting point {xlib.get_cluster_app_dir()} for the cluster {cluster_name}.\n')
            OK = False

    # check the Miniconda3 installation
    if OK:
        miniconda3_name = xlib.get_miniconda3_name()
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** error: {miniconda3_name} is not installed.\n')
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # initialize the Anaconda package list
    package_list = []

    # check the BLAST+ installation
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_blastplus_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'{xlib.get_blastplus_name()} is not installed.\n')
                (bioinfoapp_version, bioinfoapp__url, bioinfoapp_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_blastplus_anaconda_code())
                package_list.append([xlib.get_blastplus_anaconda_code(), 'last', bioinfoapp_channel])
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), xlib.get_ngshelper_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the NGShelper installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_ngshelper_installation_script()} ...\n')
        (OK, error_list) = build_ngshelper_installation_script(package_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the NGShelper installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_ngshelper_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ngshelper_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ngshelper_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the NGShelper installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ngshelper_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ngshelper_installation_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the NGShelper installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_ngshelper_installation_starter()} ...\n')
        (OK, error_list) = build_ngshelper_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the NGShelper installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_ngshelper_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ngshelper_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ngshelper_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the NGShelper installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ngshelper_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ngshelper_installation_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the NGShelper installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ngshelper_installation_starter())} ...\n')
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

def build_ngshelper_installation_script(package_list, cluster_name, current_run_dir):
    '''
    Build the NGShelper installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of NGShelper
    (ngshelper_version, ngshelper_url, ngshelper_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_ngshelper_name())

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
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write( '    echo "HOST NAME: $HOSTNAME"\n')
            script_file_id.write( '    echo "HOST IP: $HOST_IP"\n')
            script_file_id.write( '    echo "HOST ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{function remove_ngshelper_directory}\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_ngshelper_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_ngshelper_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_ngshelper_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_ngshelper_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_ngshelper_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            download_script = f'import requests; r = requests.get(\'{ngshelper_url}\') ; open(\'{xlib.get_ngshelper_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $MINICONDA3_BIN_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write( '    echo "The file is downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function decompress_ngshelper_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Decompressing the {xlib.get_ngshelper_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    unzip -u {xlib.get_ngshelper_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write( '    echo "The file is decompressed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_ngshelper_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Renaming the {xlib.get_ngshelper_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    mv {xlib.get_ngshelper_name()}-master {xlib.get_ngshelper_name()}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '    echo "The directory is renamed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function set_execution_permissions\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Setting execution permissions ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    chmod u+x {xlib.get_ngshelper_name()}/Package/*.py\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi\n')
            script_file_id.write( '    echo "The directory is renamed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_ngshelper_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_ngshelper_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_ngshelper_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "The file is removed."\n')
            script_file_id.write( '}\n')
            if len(package_list) > 0:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function add_channel_defaults\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Adding channel defaults ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write( '    ./conda config --add channels defaults\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    echo "The channel is added."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function add_channel_conda_forge\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Adding channel conda-forge ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write( '    ./conda config --add channels conda-forge\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    echo "The channel is added."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function add_channel_r\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Adding channel r ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write( '    ./conda config --add channels r\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    echo "The channel is added."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function add_channel_bioconda\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Adding channel bioconda ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write( '    ./conda config --add channels bioconda\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    echo "The channel is added."\n')
                script_file_id.write( '}\n')
            for package in package_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function install_anaconda_package_{package[0]}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "Installing {xlib.get_anaconda_name()} package {package[0]} ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write(f'    ./conda create --yes --quiet --channel {package[2]} --name {package[0]} {package[0]}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    echo "The package is installed."\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function end\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail ok\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_OK\n')
            script_file_id.write( '    exit 0\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function manage_error\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "ERROR: $1 returned error $2"\n')
            script_file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail wrong\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
            script_file_id.write( '    exit 3\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            process_name = f'{xlib.get_ngshelper_name()} installation'
            mail_message_ok = xlib.get_mail_message_ok(process_name, cluster_name)
            mail_message_wrong = xlib.get_mail_message_wrong(process_name, cluster_name)
            script_file_id.write( 'function send_mail\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {process_name}"\n')
            script_file_id.write( '    if [ "$1" == "ok" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_ok}"\n')
            script_file_id.write( '    elif [ "$1" == "wrong" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_wrong}"\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '         MESSAGE=""\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '    DESTINATION_FILE=mail-destination.json\n')
            script_file_id.write( '    echo "{" > $DESTINATION_FILE\n')
            script_file_id.write(f'    echo "    \\\"ToAddresses\\\":  [\\\"{xconfiguration.get_contact_data()}\\\"]," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"CcAddresses\\\":  []," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"BccAddresses\\\":  []" >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "}" >> $DESTINATION_FILE\n')
            script_file_id.write( '    MESSAGE_FILE=mail-message.json\n')
            script_file_id.write( '    echo "{" > $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Subject\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Data\\\":  \\\"$SUBJECT\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Body\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Html\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Data\\\":  \\\"$MESSAGE\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "}" >> $MESSAGE_FILE\n')
            script_file_id.write(f'    aws ses send-email --from {xconfiguration.get_contact_data()} --destination file://$DESTINATION_FILE --message file://$MESSAGE_FILE\n')
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
            script_file_id.write( 'remove_ngshelper_directory\n')
            script_file_id.write( 'download_ngshelper_installation_file\n')
            script_file_id.write( 'decompress_ngshelper_installation_file\n')
            script_file_id.write( 'rename_ngshelper_directory\n')
            script_file_id.write( 'set_execution_permissions\n')
            script_file_id.write( 'remove_ngshelper_installation_file\n')
            if len(package_list) > 0:
                script_file_id.write( 'add_channel_defaults\n')
                script_file_id.write( 'add_channel_conda_forge\n')
                script_file_id.write( 'add_channel_r\n')
                script_file_id.write( 'add_channel_bioconda\n')
            for package in package_list:
                script_file_id.write(f'install_anaconda_package_{package[0]}\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ngshelper_installation_script()} can not be created')
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
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_ngshelper_installation_script())} &>{current_run_dir}/{xlib.get_cluster_log_file()}')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ngshelper_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ngshelper_installation_script():
    '''
    Get the NGShelper installation path in the local computer.
    '''

    # assign the NGShelper installation path
    ngshelper_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_ngshelper_name()}-installation.sh'

    # return the NGShelper installation path
    return ngshelper_installation_script

#-------------------------------------------------------------------------------

def get_ngshelper_installation_starter():
    '''
    Get the NGShelper installation starter path in the local computer.
    '''

    # assign the NGShelper installation starter path
    ngshelper_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_ngshelper_name()}-installation-starter.sh'

    # return the NGShelper installation starter path
    return ngshelper_installation_starter

#-------------------------------------------------------------------------------

def create_vcf_sample_file():
    '''
    Create the file of VCF samples with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of restriction sites and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_vcf_sample_file())):
            os.makedirs(os.path.dirname(get_vcf_sample_file()))
        with open(get_vcf_sample_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains the sample data of a VCF file.\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: sample_id;species_id;mother_id\n')
            file_id.write( '\n')
            file_id.write( '# WARNINGS when this file is used by RADdesigner:\n')
            file_id.write( '#     - Each sample has to have a replica.\n')
            file_id.write( '#     - The replica identification have to be composed by the sample identification ended by "-dupl".\n')
            file_id.write( '\n')
            file_id.write( 'A09;AL;NONE\n')
            file_id.write( 'A09-dupl;AL;NONE\n')
            file_id.write( 'E96;AL;NONE\n')
            file_id.write( 'E96-dupl;AL;NONE\n')
            file_id.write( 'FS16;AL;NONE\n')
            file_id.write( 'FS16-dupl;AL;NONE\n')
            file_id.write( 'FS22;AL;NONE\n')
            file_id.write( 'FS22-dupl;AL;NONE\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_vcf_sample_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_vcf_sample_file(strict):
    '''
    Check the file of VCF samples.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the data records
    # format: sample_id;species_id;mother_id
    record_pattern = re.compile(r'^(.+);(.+);(.+)$')

    # open the file of VCF samples
    try:
        vcf_sample_file_id = open(get_vcf_sample_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_vcf_sample_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = vcf_sample_file_id.readline()

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try:
                    mo = record_pattern.match(record)
                    sample_id = mo.group(1).strip()
                    species_id = mo.group(2).strip()
                    mother_id = mo.group(3).strip()
                except Exception as e:
                    error_list.append(f'*** EXCEPTION: "{e}".')
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: There is a format error in the record: {record}')
                    OK = False
                    break

                # check x
                pass

            # read the next record
            record = vcf_sample_file_id.readline()

    # close the file of VCF samples
    if OK:
        vcf_sample_file_id.close()

    # warn that the file of VCF samples is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_vcf_sample_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_vcf_sample_file():
    '''
    Get the vcf_sample file path.
    '''

    # assign the VCF sample file path
    vcf_sample_file = f'{xlib.get_config_dir()}/vcf-samples.txt'

    # return the VCF sample file path
    return vcf_sample_file

#-------------------------------------------------------------------------------

def get_vcf_sample_dict():
    '''
    Get a dictionary of VCF sample enzymes with their rectriction sites.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of VCF samples
    vcf_sample_dict = {}

    # set the pattern of the data records
    # format: sample_id;species_id;mother_id
    record_pattern = re.compile(r'^(.+);(.+);(.+)$')

    # check the VCF sample file
    (OK, error_list) = check_vcf_sample_file(strict=True)

    # read every record of VCF sample file
    if OK:
        with open(get_vcf_sample_file(), mode='r', encoding='iso-8859-1', newline='\n') as file_id:

            for record in file_id:

                # if the record is not a comment nor a line with blank characters
                if not record.lstrip().startswith('#') and record.strip() != '':

                    # extract the data and add VCF sample data to the dictionary of VCF samples
                    mo = record_pattern.match(record)
                    sample_id = mo.group(1).strip()
                    species_id = mo.group(2).strip()
                    mother_id = mo.group(3).strip()

                    # add data to the dictionary of VCF samples
                    vcf_sample_dict[sample_id] = {'sample_id': sample_id, 'species_id': species_id, 'mother_id': mother_id}

    # return the control variable, the error list and the dictionary of VCF samples
    return (OK, error_list, vcf_sample_dict)

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
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The RSEM-EVAL files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/rsem_eval_dataset_id\n')
            file_id.write( '# The experiment_id and rsem_eval_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of transcript-filter (NGShelper package) and their meaning in "https://github.com/GGFHF/".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'rsem_eval_dataset_id = {rsem_eval_dataset_id}', '# rsem_eval dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the transcript-filter parameters\n')
            file_id.write( '[transcript-filter parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'minlen = 200', '# transcript with length values less than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format( 'maxlen = 10000', '# transcript with length values greater than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format( 'fpkm = 1.0', '# transcript with FPKM values less than this value will be filtered'))
            file_id.write( '{0:<50} {1}\n'.format( 'tpm = 1.0', '# transcript with TPM values less than this value will be filtered'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_transcript_filter_config_file()} can not be recreated')
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
    log.write(f'Checking the {xlib.get_transcript_filter_name()} config file ...\n')
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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the NGShelper is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ngshelper_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_ngshelper_name()} is not installed.\n')
            OK = False

    # check BLAST+ is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_blastplus_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_blastplus_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_blastplus_name()} installation could not be performed.\n')

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
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcript-filter process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_transcript_filter_process_script()} ...\n')
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
        log.write(f'Uploading the process script {get_transcript_filter_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transcript_filter_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transcript_filter_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcript-filter process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transcript_filter_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transcript_filter_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcript-filter process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_transcript_filter_process_starter()} ...\n')
        (OK, error_list) = build_transcript_filter_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the transcript-filter process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_transcript_filter_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transcript_filter_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transcript_filter_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcript-filter process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transcript_filter_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transcript_filter_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the transcript-filter process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_transcript_filter_process_starter())} ...\n')
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
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
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
                error_list.append(f'*** ERROR: the key "rsem_eval_dataset_id" value is not a {xlib.get_rsem_eval_name()} assembly assessment.')
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
                error_list.append(f'*** ERROR: The value maxlen value ({maxlen}) is less than the minlen value ({minlen}).')
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
        error_list.append(f'\nThe {xlib.get_transcript_filter_name()} config file is not valid. Please, correct this file or recreate it.')

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
    local_log_file = f'{xlib.get_temp_dir()}/{xlib.get_cluster_log_file()}'
    cluster_log_file = f'{xlib.get_cluster_experiment_result_dir(experiment_id)}/{rsem_eval_dataset_id}/{xlib.get_cluster_log_file()}'

    # download the RSEM-EVAL log file from the cluster
    OK = xssh.get_file(sftp_client, cluster_log_file, local_log_file)
    if not OK:
        error_list.append(f'*** ERROR: The file {cluster_log_file} does not have been downloaded.')

    # get the assembly software, result_data_set_id and assembly type from the RSEM-EVAL log file
    if OK:
        assembly_software = ''
        assembly_dataset_id = ''
        assembly_type = ''
        filtering_data_id = 'FILTERING_DATA'
        pattern = f'{filtering_data_id} - (.+): (.+)'
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
            error_list.append(f'*** ERROR: Some filtering are not in the file {cluster_log_file}.')
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
        score_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, rsem_eval_dataset_id)}/{rsem_eval_dataset_id}.genes.results'

    # set the output file path
    if OK:
        output_file = f'{current_run_dir}/filtered-transcriptome.fasta'

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
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write(f'NGSHELPER_PATH={xlib.get_cluster_app_dir()}/{xlib.get_ngshelper_name()}/Package\n')
                script_file_id.write(f'export PATH=$MINICONDA3_BIN_PATH:$NGSHELPER_PATH:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function init\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
                script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
                script_file_id.write( '    echo "HOST NAME: $HOSTNAME"\n')
                script_file_id.write( '    echo "HOST IP: $HOST_IP"\n')
                script_file_id.write( '    echo "HOST ADDRESS: $HOST_ADDRESS"\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function run_transcript_filter_process\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Filtering the transcripts ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '        transcript-filter.py \\\n')
                script_file_id.write( '            --assembler=ngscloud \\\n')
                script_file_id.write(f'            --transcriptome={transcriptome_file} \\\n')
                script_file_id.write(f'            --score={score_file} \\\n')
                script_file_id.write(f'            --output={output_file} \\\n')
                script_file_id.write(f'            --minlen={minlen} \\\n')
                script_file_id.write(f'            --maxlen={maxlen} \\\n')
                script_file_id.write(f'            --FPKM={fpkm} \\\n')
                script_file_id.write(f'            --TPM={tpm} \\\n')
                script_file_id.write( '            --verbose=n\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error transcript-filter.py $RC; fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function end\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
                script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
                script_file_id.write( '    calculate_duration\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    send_mail ok\n')
                script_file_id.write( '    touch $SCRIPT_STATUS_OK\n')
                script_file_id.write( '    exit 0\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function manage_error\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
                script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
                script_file_id.write( '    calculate_duration\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "ERROR: $1 returned error $2"\n')
                script_file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    send_mail wrong\n')
                script_file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
                script_file_id.write( '    exit 3\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                process_name = f'{xlib.get_transcript_filter_name()} process'
                mail_message_ok = xlib.get_mail_message_ok(process_name, cluster_name)
                mail_message_wrong = xlib.get_mail_message_wrong(process_name, cluster_name)
                script_file_id.write( 'function send_mail\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {process_name}"\n')
                script_file_id.write( '    if [ "$1" == "ok" ]; then\n')
                script_file_id.write(f'        MESSAGE="{mail_message_ok}"\n')
                script_file_id.write( '    elif [ "$1" == "wrong" ]; then\n')
                script_file_id.write(f'        MESSAGE="{mail_message_wrong}"\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '         MESSAGE=""\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '    DESTINATION_FILE=mail-destination.json\n')
                script_file_id.write( '    echo "{" > $DESTINATION_FILE\n')
                script_file_id.write(f'    echo "    \\\"ToAddresses\\\":  [\\\"{xconfiguration.get_contact_data()}\\\"]," >> $DESTINATION_FILE\n')
                script_file_id.write( '    echo "    \\\"CcAddresses\\\":  []," >> $DESTINATION_FILE\n')
                script_file_id.write( '    echo "    \\\"BccAddresses\\\":  []" >> $DESTINATION_FILE\n')
                script_file_id.write( '    echo "}" >> $DESTINATION_FILE\n')
                script_file_id.write( '    MESSAGE_FILE=mail-message.json\n')
                script_file_id.write( '    echo "{" > $MESSAGE_FILE\n')
                script_file_id.write( '    echo "    \\\"Subject\\\": {" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "        \\\"Data\\\":  \\\"$SUBJECT\\\"," >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "        \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "    }," >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "    \\\"Body\\\": {" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "        \\\"Html\\\": {" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "            \\\"Data\\\":  \\\"$MESSAGE\\\"," >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "            \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "        }" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "    }" >> $MESSAGE_FILE\n')
                script_file_id.write( '    echo "}" >> $MESSAGE_FILE\n')
                script_file_id.write(f'    aws ses send-email --from {xconfiguration.get_contact_data()} --destination file://$DESTINATION_FILE --message file://$MESSAGE_FILE\n')
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
                script_file_id.write(f'run_transcript_filter_process\n')
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {get_transcript_filter_process_script()} can not be created.')
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
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_transcript_filter_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_transcript_filter_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transcript_filter_config_file():
    '''
    Get the transcript-filter config file path.
    '''

    # assign the transcript-filter config file path
    transcript_filter_config_file = f'{xlib.get_config_dir()}/{xlib.get_transcript_filter_code()}-config.txt'

    # return the transcript-filter config file path
    return transcript_filter_config_file

#-------------------------------------------------------------------------------

def get_transcript_filter_process_script():
    '''
    Get the transcript-filter process script path in the local computer.
    '''

    # assign the transcript-filter script path
    transcript_filter_process_script = f'{xlib.get_temp_dir()}/{xlib.get_transcript_filter_code()}-process.sh'

    # return the transcript-filter script path
    return transcript_filter_process_script

#-------------------------------------------------------------------------------

def get_transcript_filter_process_starter():
    '''
    Get the transcript-filter process starter path in the local computer.
    '''

    # assign the transcript-filter process starter path
    transcript_filter_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_transcript_filter_code()}-process-starter.sh'

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
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The database files have to be located in the cluster directory {xlib.get_cluster_database_dir()}/database_dataset_id\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write( '# The experiment_id, database_dataset_id and assembly_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of transcriptome-blastx (NGShelper package) and their meaning in "https://github.com/GGFHF/NGShelper".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'database_dataset_id = {database_dataset_id}', '# database dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'protein_database_name = {protein_database_name}', '# protein database name'))
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the transcriptome-blastx parameters\n')
            file_id.write( '[transcriptome-blastx parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'node_number = 1', '# node number (previously they have to be started)'))
            file_id.write( '{0:<50} {1}\n'.format( 'blastx_thread_number = 1', '# threads number using by blastx in every node'))
            file_id.write( '{0:<50} {1}\n'.format( 'e_value = 1E-6', '# expectation value (E-value) threshold for saving hits'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_target_seqs = 10', '# maximum number of aligned sequences to keep'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_hsps = 999999', '# maximum number of HSPs per subject sequence to save for each query'))
            file_id.write( '{0:<50} {1}\n'.format( 'qcov_hsp_perc = 0.0', '# alignments below the specified query coverage per HSPs are removed'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_transcriptome_blastx_config_file()} can not be recreated')
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
    log.write(f'Checking the {xlib.get_transcriptome_blastx_name()} config file ...\n')
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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
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
            log.write(f'*** ERROR: The requested node number ({requested_node_number}) is greater than running node number ({running_node_number}).\n')
            OK = False

    # check the NGShelper is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ngshelper_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_ngshelper_name()} is not installed.\n')
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
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcriptome-blastx process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_transcriptome_blastx_process_script()} ...\n')
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
        log.write(f'Uploading the process script {get_transcriptome_blastx_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transcriptome_blastx_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcriptome-blastx process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the transcriptome-blastx process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_transcriptome_blastx_process_starter()} ...\n')
        (OK, error_list) = build_transcriptome_blastx_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the transcriptome-blastx process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_transcriptome_blastx_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transcriptome_blastx_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the transcriptome-blastx process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the transcriptome-blastx process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_starter())} ...\n')
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
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
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
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = transcriptome_blastx_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" has to start with {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = transcriptome_blastx_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
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
        error_list.append(f'\nThe {xlib.get_transcriptome_blastx_name()} config file is not valid. Please, correct this file or recreate it.')

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
    try:
        if not os.path.exists(os.path.dirname(get_transcriptome_blastx_process_script())):
            os.makedirs(os.path.dirname(get_transcriptome_blastx_process_script()))
        with open(get_transcriptome_blastx_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'NGSHELPER_PATH={xlib.get_cluster_app_dir()}/{xlib.get_ngshelper_name()}/Package\n')
            script_file_id.write(f'BLASTPLUS_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_blastplus_anaconda_code()}/bin\n')
            script_file_id.write( 'export PATH=$MINICONDA3_BIN_PATH:$NGSHELPER_PATH:$BLASTPLUS_PATH:$PATH\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function init\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    INIT_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_INIT_DATETIME=`date --date="@$INIT_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME+00:00."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "CLUSTER: {cluster_name}"\n')
            script_file_id.write( '    echo "HOST NAME: $HOSTNAME"\n')
            script_file_id.write( '    echo "HOST IP: $HOST_IP"\n')
            script_file_id.write( '    echo "HOST ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_transcriptome_blastx_process\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running the transcriptome blastx process ..."\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '        transcriptome-blastx.py \\\n')
            script_file_id.write( '            --machine_type=ngscloud \\\n')
            script_file_id.write(f'            --node_number={node_number} \\\n')
            script_file_id.write(f'            --blastx_thread_number={blastx_thread_number} \\\n')
            script_file_id.write(f'            --blast_db={xlib.get_cluster_database_dataset_dir(database_dataset_id)} \\\n')
            script_file_id.write(f'            --protein_database_name={protein_database_name} \\\n')
            script_file_id.write(f'            --transcriptome={transcriptome_file} \\\n')
            script_file_id.write(f'            --e_value={e_value} \\\n')
            script_file_id.write(f'            --max_target_seqs={max_target_seqs} \\\n')
            script_file_id.write(f'            --max_hsps={max_hsps} \\\n')
            script_file_id.write(f'            --qcov_hsp_perc={qcov_hsp_perc} \\\n')
            script_file_id.write(f'            --output={current_run_dir} \\\n')
            script_file_id.write(f'            --email={xconfiguration.get_contact_data()} \\\n')
            script_file_id.write( '            --verbose=n\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error transcriptome-blastx.py $RC; fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function end\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Script ended OK at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail ok\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_OK\n')
            script_file_id.write( '    exit 0\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function manage_error\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    END_DATETIME=`date --utc +%s`\n')
            script_file_id.write( '    FORMATTED_END_DATETIME=`date --date="@$END_DATETIME" "+%Y-%m-%d %H:%M:%S"`\n')
            script_file_id.write( '    calculate_duration\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "ERROR: $1 returned error $2"\n')
            script_file_id.write( '    echo "Script ended WRONG at $FORMATTED_END_DATETIME+00:00 with a run duration of $DURATION s ($FORMATTED_DURATION)."\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    send_mail wrong\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
            script_file_id.write( '    exit 3\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            process_name = f'{xlib.get_transcriptome_blastx_name()} process'
            mail_message_ok = xlib.get_mail_message_ok(process_name, cluster_name)
            mail_message_wrong = xlib.get_mail_message_wrong(process_name, cluster_name)
            script_file_id.write( 'function send_mail\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {process_name}"\n')
            script_file_id.write( '    if [ "$1" == "ok" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_ok}"\n')
            script_file_id.write( '    elif [ "$1" == "wrong" ]; then\n')
            script_file_id.write(f'        MESSAGE="{mail_message_wrong}"\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '         MESSAGE=""\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '    DESTINATION_FILE=mail-destination.json\n')
            script_file_id.write( '    echo "{" > $DESTINATION_FILE\n')
            script_file_id.write(f'    echo "    \\\"ToAddresses\\\":  [\\\"{xconfiguration.get_contact_data()}\\\"]," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"CcAddresses\\\":  []," >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "    \\\"BccAddresses\\\":  []" >> $DESTINATION_FILE\n')
            script_file_id.write( '    echo "}" >> $DESTINATION_FILE\n')
            script_file_id.write( '    MESSAGE_FILE=mail-message.json\n')
            script_file_id.write( '    echo "{" > $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Subject\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Data\\\":  \\\"$SUBJECT\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    \\\"Body\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        \\\"Html\\\": {" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Data\\\":  \\\"$MESSAGE\\\"," >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "            \\\"Charset\\\":  \\\"UTF-8\\\"" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "        }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "    }" >> $MESSAGE_FILE\n')
            script_file_id.write( '    echo "}" >> $MESSAGE_FILE\n')
            script_file_id.write(f'    aws ses send-email --from {xconfiguration.get_contact_data()} --destination file://$DESTINATION_FILE --message file://$MESSAGE_FILE\n')
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
            script_file_id.write(f'run_transcriptome_blastx_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_transcriptome_blastx_process_script()} can not be created.')
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
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_transcriptome_blastx_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_transcriptome_blastx_process_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_config_file():
    '''
    Get the transcriptome-blastx config file path.
    '''

    # assign the transcriptome-blastx config file path
    transcriptome_blastx_config_file = f'{xlib.get_config_dir()}/{xlib.get_transcriptome_blastx_code()}-config.txt'

    # return the transcriptome-blastx config file path
    return transcriptome_blastx_config_file

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_process_script():
    '''
    Get the transcriptome-blastx process script path in the local computer.
    '''

    # assign the transcriptome-blastx script path
    transcriptome_blastx_process_script = f'{xlib.get_temp_dir()}/{xlib.get_transcriptome_blastx_code()}-process.sh'

    # return the transcriptome-blastx script path
    return transcriptome_blastx_process_script

#-------------------------------------------------------------------------------

def get_transcriptome_blastx_process_starter():
    '''
    Get the transcriptome-blastx process starter path in the local computer.
    '''

    # assign the transcriptome-blastx process starter path
    transcriptome_blastx_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_transcriptome_blastx_code()}-process-starter.sh'

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

    return f'{xlib.get_soapdenovotrans_code()} ({xlib.get_soapdenovotrans_name()}) or {xlib.get_transabyss_code()} ({xlib.get_transabyss_name()}) or {xlib.get_trinity_code()} ({xlib.get_trinity_name()}) or {xlib.get_ggtrinity_code()} ({xlib.get_ggtrinity_name()}) or {xlib.get_cd_hit_est_code()} ({xlib.get_cd_hit_est_name()}) or {xlib.get_transcript_filter_code()} ({xlib.get_transcript_filter_name()})'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the NGShelper process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
