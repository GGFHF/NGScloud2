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
This file contains functions related to the RADdesigner process used in both console
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
import xngshelper
import xssh

#-------------------------------------------------------------------------------

def is_installed_raddesigner(cluster_name, passed_connection, ssh_client):
    '''
    Check if RADdesigner is installed.
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

    # check the RADdesigner directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_raddesigner_name()} ] && echo RC=0 || echo RC=1'
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

def install_raddesigner(cluster_name, log, function=None):
    '''
    Install the RADdesigner software in the cluster.
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

    # check the Miniconda3 setup
    if OK:
        miniconda3_name = xlib.get_miniconda3_name()
        (OK, error_list, is_setup) = xbioinfoapp.is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_setup:
                log.write(f'*** error: {miniconda3_name} is not setup. Previously it has to be set up.\n')
                OK = False
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # initialize the Anaconda package list
    package_list = []

    # check the VCFtools Perl libraries setup
    if OK:
        (OK, error_list, is_setup) = xbioinfoapp.is_installed_anaconda_package(xlib.get_vcftools_perl_libraries_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_setup:
                log.write(f'{xlib.get_vcftools_perl_libraries_name()} is not set up. It has to be set up.\n')
                (bioinfoapp_version, bioinfoapp__url, bioinfoapp_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_vcftools_perl_libraries_name())
                package_list.append([xlib.get_vcftools_perl_libraries_anaconda_code(), 'last', bioinfoapp_channel])
        else:
            log.write('*** ERROR: The verification can not run.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_raddesigner_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the RADdesigner installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_raddesigner_installation_script()} ...\n')
        (OK, error_list) = build_raddesigner_installation_script(package_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the RADdesigner installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_raddesigner_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_raddesigner_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_raddesigner_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RADdesigner installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_raddesigner_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_raddesigner_installation_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the RADdesigner installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_raddesigner_installation_starter()} ...\n')
        (OK, error_list) = build_raddesigner_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the RADdesigner installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_raddesigner_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_raddesigner_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_raddesigner_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RADdesigner installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_raddesigner_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_raddesigner_installation_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the RADdesigner installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_raddesigner_installation_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_raddesigner_installation_starter()), log)

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

def build_raddesigner_installation_script(package_list, cluster_name, current_run_dir):
    '''
    Build the RADdesigner installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of RADdesigner
    (raddesigner_version, raddesigner_url, raddesigner_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_raddesigner_name())

    # write the RADdesigner installation script
    try:
        if not os.path.exists(os.path.dirname(get_raddesigner_installation_script())):
            os.makedirs(os.path.dirname(get_raddesigner_installation_script()))
        with open(get_raddesigner_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function remove_raddesigner_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_raddesigner_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_raddesigner_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_raddesigner_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_raddesigner_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_raddesigner_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            download_script = f'import requests; r = requests.get(\'{raddesigner_url}\') ; open(\'{xlib.get_raddesigner_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $MINICONDA3_BIN_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write( '    echo "The file is downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function decompress_raddesigner_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Decompressing the {xlib.get_raddesigner_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    unzip -u {xlib.get_raddesigner_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write( '    echo "The file is decompressed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_raddesigner_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Renaming the {xlib.get_raddesigner_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    mv {xlib.get_raddesigner_name()}-master {xlib.get_raddesigner_name()}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '    echo "The directory is renamed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_raddesigner_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_raddesigner_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_raddesigner_name()}.zip\n')
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
            process_name = f'{xlib.get_raddesigner_name()} installation'
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
            script_file_id.write( 'remove_raddesigner_directory\n')
            script_file_id.write( 'download_raddesigner_installation_file\n')
            script_file_id.write( 'decompress_raddesigner_installation_file\n')
            script_file_id.write( 'rename_raddesigner_directory\n')
            script_file_id.write( 'remove_raddesigner_installation_file\n')
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
        error_list.append(f'*** ERROR: The file {get_raddesigner_installation_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_raddesigner_installation_starter(current_run_dir):
    '''
    Build the starter of the RADdesigner installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the RADdesigner installation starter
    try:
        if not os.path.exists(os.path.dirname(get_raddesigner_installation_starter())):
            os.makedirs(os.path.dirname(get_raddesigner_installation_starter()))
        with open(get_raddesigner_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_raddesigner_installation_script())} &>{current_run_dir}/{xlib.get_cluster_log_file()}')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_raddesigner_installation_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_raddesigner_installation_script():
    '''
    Get the RADdesigner installation path in the local computer.
    '''

    # assign the RADdesigner installation path
    raddesigner_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_raddesigner_name()}-installation.sh'

    # return the RADdesigner installation path
    return raddesigner_installation_script

#-------------------------------------------------------------------------------

def get_raddesigner_installation_starter():
    '''
    Get the RADdesigner installation starter path in the local computer.
    '''

    # assign the RADdesigner installation starter path
    raddesigner_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_raddesigner_name()}-installation-starter.sh'

    # return the RADdesigner installation starter path
    return raddesigner_installation_starter

#-------------------------------------------------------------------------------

def create_condition_file():
    '''
    Create the file of conditions with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of restriction sites and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_condition_file())):
            os.makedirs(os.path.dirname(get_condition_file()))
        with open(get_condition_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( 'Read_type,Enzyme,Methodology,Read_number,parameters\n')
            file_id.write( 'PE,MslI,GBS,1_5M,095_10_10\n')
            file_id.write( 'PE,MslI,GBS,3M,095_10_10\n')
            file_id.write( 'PE,PstI_MspI,ddRAD,1_5M,095_10_10\n')
            file_id.write( 'PE,PstI_MspI,ddRAD,3M,095_10_10\n')
            file_id.write( 'SE,MslI,GBS,1_5M,095_10_10\n')
            file_id.write( 'SE,MslI,GBS,3M,095_10_10\n')
            file_id.write( 'SE,PstI_MspI,ddRAD,1_5M,095_10_10\n')
            file_id.write( 'SE,PstI_MspI,ddRAD,3M,095_10_10\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_condition_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_condition_file(strict):
    '''
    Check the file of conditions.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the data records
    # format: read_type,enzyme,methodology,read_number,parameters
    record_pattern = re.compile(r'^(.+),(.+),(.+),(.+),(.*)$')

    # set the first record
    first_record = 'Read_type,Enzyme,Methodology,Read_number,parameters'

    # open the file of conditions
    try:
        condition_file_id = open(get_condition_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_condition_file()} can not be opened.')
        OK = False

    # intitialize the record counter
    record_counter = 0

    # check that all records are OK
    if OK:

        # read the first record
        record = condition_file_id.readline()

        # while there are records
        while record != '':

            # check the first record
            if record_counter == 0:
                if record.strip() != first_record:
                    error_list.append(f'*** ERROR: The first record has to be: {first_record}')

            # check remainder records
            else:

                # extract the data
                try:
                    mo = record_pattern.match(record)
                    read_type = mo.group(1).strip()
                    enzyme = mo.group(2).strip()
                    methodology = mo.group(3).strip()
                    read_number = mo.group(4).strip()
                    parameters = mo.group(5).strip()
                except Exception as e:
                    error_list.append(f'*** EXCEPTION: "{e}".')
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: There is a format error in the record: {record}')
                    OK = False
                    break

            # add 1 to the record counter
            record_counter += 1

            # read the next record
            record = condition_file_id.readline()

    # close the file of conditions
    condition_file_id.close()

    # warn that the file of conditions is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_condition_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_condition_file():
    '''
    Get the condition file path.
    '''

    # assign the condition file path
    condition_file = f'{xlib.get_config_dir()}/{xlib.get_raddesigner_code()}-conditions.txt'

    # return the condition file path
    return condition_file

#-------------------------------------------------------------------------------

def create_raddesigner_config_file(experiment_id='exp001', vcf_location_dataset_id='ipyrad-170101-235959'):
    '''
    Create RADdesigner config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the RADdesigner config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_raddesigner_config_file())):
            os.makedirs(os.path.dirname(get_raddesigner_config_file()))
        with open(get_raddesigner_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The VCF files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/vcf_location_dataset_id\n')
            file_id.write( '# The experiment_id AND vcf_location_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of RADdesigner and their meaning in "https://github.com/GGFHF/RADdesigner".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'vcf_location_dataset_id = {vcf_location_dataset_id}', '# dataset identification where the VCF files are located.'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the RADdesigner parameters\n')
            file_id.write( '[RADdesigner parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use (greater than or equal to 4)'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_raddesigner_config_file()} can not be recreated.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_raddesigner_process(cluster_name, log, function=None):
    '''
    Run a RADdesigner process.
    '''

    # initialize the control variable
    OK = True

    # get the RADdesigner option dictionary
    raddesigner_option_dict = xlib.get_option_dict(get_raddesigner_config_file())

    # get the experiment identification
    experiment_id = raddesigner_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the RADdesigner config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_raddesigner_name()} config file ...\n')
    (OK, error_list) = check_raddesigner_config_file(strict=True)
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

    # check the RADdesigner is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_raddesigner_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_raddesigner_name()} is not installed.\n')
            OK = False

    # check R is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_r_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_r_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_r_name()} installation could not be performed.\n')

    # check VCFtools Perl libraries is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_vcftools_perl_libraries_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_vcftools_perl_libraries_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_vcftools_perl_libraries_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_raddesigner_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # check the RADdesigner condition file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the file {get_condition_file()} ...\n')
    (OK, error_list) = check_condition_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The condition file is not valid.\n')
        log.write('Please correct this file or recreate it.\n')

    # upload the RADdesigner condition file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {get_condition_file()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_condition_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_condition_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the the file with sample identifications
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the file {get_sample_file()} ...\n')
        (OK, error_list) = build_sample_file()
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the file with sample identifications to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {get_sample_file()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_sample_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_sample_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the RADdesigner process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_raddesigner_process_script()} ...\n')
        (OK, error_list) = build_raddesigner_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the RADdesigner process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_raddesigner_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_raddesigner_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_raddesigner_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RADdesigner process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_raddesigner_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_raddesigner_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the RADdesigner process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_raddesigner_process_starter()} ...\n')
        (OK, error_list) = build_raddesigner_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the RADdesigner process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_raddesigner_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_raddesigner_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_raddesigner_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RADdesigner process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_raddesigner_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_raddesigner_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the RADdesigner process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_raddesigner_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_raddesigner_process_starter()), log)

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

def check_raddesigner_config_file(strict):
    '''
    Check the RADdesigner config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        raddesigner_option_dict = xlib.get_option_dict(get_raddesigner_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in raddesigner_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = raddesigner_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "vcf_location_dataset_id"
            vcf_location_dataset_id = raddesigner_option_dict.get('identification', {}).get('vcf_location_dataset_id', not_found)
            if vcf_location_dataset_id == not_found:
                error_list.append('*** ERROR: the key "vcf_location_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "RADdesigner parameters"
        if 'RADdesigner parameters' not in sections_list:
            error_list.append('*** ERROR: the section "RADdesigner parameters" is not found.')
            OK = False
        else:

            # check section "RADdesigner parameters" - key "threads"
            threads = raddesigner_option_dict.get('RADdesigner parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "RADdesigner parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=4):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 4.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_raddesigner_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_sample_file():
    '''
    Build the sample file from VCF sample file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize sample and replica lists
    sample_list = []
    replica_list = []

    # get 
    (OK, error_list, vcf_sample_dict) = xngshelper.get_vcf_sample_dict()

    # get the sample identifications and load the sample and replica lists
    if OK:
        for sample_id in vcf_sample_dict.keys():
            if not sample_id.endswith('-dupl'):
                sample_list.append(sample_id)
            else:
                replica_list.append(sample_id)

    # check the item number of the sample and replica lists
    if OK:
        if len(sample_list) == 0:
            error_list.append(f'*** ERROR: There are not samples.')
            OK = False
        elif len(sample_list) != len(replica_list):
            error_list.append(f'*** ERROR: The samples number and the replica number are different.')
            OK = False

    # sort sample and replica lists
    if OK:
        sample_list.sort()
        replica_list.sort()

    # check items of the sample and replica lists
    if OK:
        for i in range(len(sample_list)):
            if replica_list[i] != f'{sample_list[i]}-dupl':
                error_list.append(f'*** ERROR: Samples and replicas are mismatched.')
                OK = False
                break

    # open the file of samples
    if OK:
        try:
            sample_file_id = open(get_sample_file(), mode='w', encoding='iso-8859-1', newline='\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {get_sample_file()} can not be opened.')
            OK = False

    # write samples and replicas in the file of samples
    if OK:
        for i in range(len(sample_list)):
            sample_file_id.write(f'{sample_list[i]}\n')
            sample_file_id.write(f'{replica_list[i]}\n')

    # close the file of samples
    if OK:
        sample_file_id.close()

    # warn that the file of individuals is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {xngshelper.get_vcf_sample_file()} is not valid for {xlib.get_raddesigner_name()}. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_raddesigner_process_script(cluster_name, current_run_dir):
    '''
    Build the current RADdesigner process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the RADdesigner option dictionary
    raddesigner_option_dict = xlib.get_option_dict(get_raddesigner_config_file())

    # get the options
    experiment_id = raddesigner_option_dict['identification']['experiment_id']
    vcf_location_dataset_id = raddesigner_option_dict['identification']['vcf_location_dataset_id']
    threads = raddesigner_option_dict['RADdesigner parameters']['threads']

    # set the input VCF directory path
    if vcf_location_dataset_id.startswith(xlib.get_ipyrad_code()):
        input_vcf_dir = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, vcf_location_dataset_id)}/{experiment_id.replace("-","_")}_outfiles'
    else:
        input_vcf_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, vcf_location_dataset_id)

    # write the RADdesigner process script
    try:
        if not os.path.exists(os.path.dirname(get_raddesigner_process_script())):
            os.makedirs(os.path.dirname(get_raddesigner_process_script()))
        with open(get_raddesigner_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'R_DIR={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_r_anaconda_code()}/bin\n')
            script_file_id.write(f'VCFTOOLS_PERL_LIBRARIES_DIR={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_vcftools_perl_libraries_anaconda_code()}/bin\n')
            script_file_id.write(f'RADDESIGNER_PATH={xlib.get_cluster_app_dir()}/{xlib.get_raddesigner_name()}/bin\n')
            script_file_id.write(f'export PATH=$MINICONDA3_BIN_PATH:$R_DIR:$VCFTOOLS_PERL_LIBRARIES_DIR:$PATH\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STEP_1_DIR={current_run_dir}/1.Filtering_multiallelicSNPs\n')
            script_file_id.write( 'if [ ! -d "$STEP_1_DIR" ]; then mkdir --parents $STEP_1_DIR; fi\n')
            script_file_id.write( 'STEP_1_FILTERED_DIR=$STEP_1_DIR/FILTERED_OUTPUT\n')
            script_file_id.write( 'if [ ! -d "$STEP_1_FILTERED_DIR" ]; then mkdir --parents $STEP_1_FILTERED_DIR; fi\n')
            script_file_id.write( 'STEP_1_BIALLELIC_DIR=$STEP_1_DIR/Bi_allelic\n')
            script_file_id.write( 'if [ ! -d "$STEP_1_BIALLELIC_DIR" ]; then mkdir --parents $STEP_1_BIALLELIC_DIR; fi\n')
            script_file_id.write( 'STEP_1_MULTIALLELIC_DIR=$STEP_1_DIR/Multi_allelic\n')
            script_file_id.write( 'if [ ! -d "$STEP_1_MULTIALLELIC_DIR" ]; then mkdir --parents $STEP_1_MULTIALLELIC_DIR; fi\n')
            script_file_id.write(f'STEP_2_DIR={current_run_dir}/2.Depth_locusSNPerror\n')
            script_file_id.write( 'if [ ! -d "$STEP_2_DIR" ]; then mkdir --parents $STEP_2_DIR; fi\n')
            script_file_id.write(f'STEP_3_DIR={current_run_dir}/3.Prepare_vcf\n')
            script_file_id.write( 'if [ ! -d "$STEP_3_DIR" ]; then mkdir --parents $STEP_3_DIR; fi\n')
            script_file_id.write( 'STEP_3_TXT_DIR=$STEP_3_DIR/OUTPUT_txt\n')
            script_file_id.write( 'if [ ! -d "$STEP_3_TXT_DIR" ]; then mkdir --parents $STEP_3_TXT_DIR; fi\n')
            script_file_id.write(f'STEP_4_DIR={current_run_dir}/4.Filtering_locusSNPerror\n')
            script_file_id.write( 'if [ ! -d "$STEP_4_DIR" ]; then mkdir --parents $STEP_4_DIR; fi\n')
            script_file_id.write(f'SAMPLE_FILE={current_run_dir}/{os.path.basename(get_sample_file())}\n')
            script_file_id.write(f'CONDITION_FILE={current_run_dir}/{os.path.basename(get_condition_file())}\n')
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
            script_file_id.write( 'function copy_r_scripts\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/copy_r_scripts.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Copying the RADdesigner R scripts to the current directory ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        cp $RADDESIGNER_PATH/1.Filtering_multiallelicSNPs/*.R .\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error cp $RC; fi\n')
            script_file_id.write( '        cp $RADDESIGNER_PATH/2.Depth_locusSNPerror/*.R .\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error cp $RC; fi\n')
            script_file_id.write( '        cp $RADDESIGNER_PATH/2.Depth_locusSNPerror/Mastretta_functions/*.R .\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error cp $RC; fi\n')
            script_file_id.write( '        cp $RADDESIGNER_PATH/4.Filtering_locusSNPerror/*.R .\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error cp $RC; fi\n')
            script_file_id.write( '        echo "The file are copied."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function fix_r_scripts\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/fix_r_scripts.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Fixing the RADdesigner R scripts ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        sed -i "s|library(githubinstall)|#library(githubinstall)|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dir<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs\\")|dir<-(\\"{current_run_dir}/1.Filtering_multiallelicSNPs\\")|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dirInput<-(\\"/PATH_TO_WORKING_DIRECTORY/Input_files\\") ##Same for all the scripts|dirInput<-(\\"{current_run_dir}\\")|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|Ch_names <- read.csv(paste(dirInput, \\"Ch_combinations.csv\\", sep = \\"/\\"))|Ch_names <- read.csv(\\"$CONDITION_FILE\\")|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|S_names <- read.table(paste(dirInput,\\"Samples_names.txt\\", sep = \\"/\\"))|S_names <- read.table(\\"$SAMPLE_FILE\\")|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|filelist <- list.files((paste(dirInput,\\"INPUT_VCF_DATA\\", sep = \\"/\\"))|filelist <- list.files(\\"{input_vcf_dir}\\"|g" multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '        sed -i "s|library(githubinstall)|#library(githubinstall)|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dir<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs\\")|dir<-(\\"{current_run_dir}/1.Filtering_multiallelicSNPs\\")|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dirInput<-(\\"/PATH_TO_WORKING_DIRECTORY/Input_files\\") ##Same for all the scripts|dirInput<-(\\"{current_run_dir}\\")|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|Ch_names <- read.csv(paste(dirInput, \\"Ch_combinations.csv\\", sep = \\"/\\"))|Ch_names <- read.csv(\\"$CONDITION_FILE\\")|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|S_names <- read.table(paste(dirInput,\\"Samples_names.txt\\", sep = \\"/\\"))|S_names <- read.table(\\"$SAMPLE_FILE\\")|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|filelist <- list.files((paste(dirInput,\\"INPUT_VCF_DATA\\", sep = \\"/\\"))|filelist <- list.files(\\"{input_vcf_dir}\\"|g" plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '        sed -i "s|dir<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/2.Depth_locusSNPerror\\")|dir<-(\\"$STEP_2_DIR\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dirInput<-(\\"/PATH_TO_WORKING_DIRECTORY/Input_files\\") ##Same for all the scripts|dirInput<-(\\"{current_run_dir}\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dirInput_vcf<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs/FILTERED_OUTPUT\\")|dirInput_vcf<-(\\"{current_run_dir}/1.Filtering_multiallelicSNPs/FILTERED_OUTPUT\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|Ch_names <- read.csv(paste(dirInput, \\"Ch_combinations.csv\\", sep = \\"/\\"))|Ch_names <- read.csv(\\"$CONDITION_FILE\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|S_names <- read.table(paste(dirInput,\\"Samples_names.txt\\", sep = \\"/\\"))|S_names <- read.table(\\"$SAMPLE_FILE\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|  myBoots[[i]] <- boot.phylo(dendro_all[[i]], genlight_b[[i]], dendrofunction, mc.cores = 4)|  myBoots[[i]] <- boot.phylo(dendro_all[[i]], genlight_b[[i]], dendrofunction, mc.cores = {threads})|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|source(\\"Mastretta_functions/LociAllele_error_pyrad.R\\")|source(\\"{current_run_dir}/LociAllele_error_pyrad.R\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|source(\\"Mastretta_functions/SNPs_error.R\\")|source(\\"{current_run_dir}/SNPs_error.R\\")|g" 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '        sed -i "s|dir<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/4.Filtering_locusSNPerror\\")|dir<-(\\"$STEP_4_DIR\\")|g" 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|dirInput<-(\\"/PATH_TO_WORKING_DIRECTORY/Input_files\\") ##Same for all the scripts|dirInput<-(\\"{current_run_dir}\\")|g" 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|Ch_names <- read.csv(paste(dirInput, \\"Ch_combinations.csv\\", sep = \\"/\\"))|Ch_names <- read.csv(\\"$CONDITION_FILE\\")|g" 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write(f'        sed -i "s|S_names <- read.table(paste(dirInput,\\"Samples_names.txt\\", sep = \\"/\\"))|S_names <- read.table(\\"$SAMPLE_FILE\\")|g" 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '        sed -i "s|dirInput_txt<-(\\"/PATH_TO_WORKING_DIRECTORY/bin/3.Prepare_vcf/OUTPUT_txt\\")|dirInput_txt<-(\\"$STEP_3_TXT_DIR\\")|g" 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '        echo "R scripts are fixed."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_1_1_save_headers\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_1_1_save_headers.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 1.1: Saving the header (first ten rows) of files *.vcf in files header*.txt ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        ls {input_vcf_dir}/*.vcf > IpyRAD_output_vcf.txt\n')
            script_file_id.write( '        while read FILE; do\n')
            script_file_id.write( '            head $FILE > $STEP_1_DIR/header`basename "$FILE"`.txt\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error head $RC; fi\n')
            script_file_id.write( '        done < IpyRAD_output_vcf.txt\n')
            script_file_id.write( '        echo "The headers are save."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_1_2_separate_snps\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_1_2_separate_snps.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 1.2: Separating the bi- and multi-allelic SNPs ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_r_anaconda_code()}\n')
            script_file_id.write( '        Rscript --vanilla --verbose multi_bi_allelic_filtered.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error Rscript $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "SNPs are separated."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_1_3_paste_headers\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_1_3_paste_headers.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 1.3: Pasting the header to the bi_allelic*.vcf files ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        ls $STEP_1_DIR/*bi_allelic.vcf > bi_allelic_files.txt\n')
            script_file_id.write( '        ls $STEP_1_DIR/header*.txt > header_files.txt\n')
            script_file_id.write( '        while read FILE1; do\n')
            script_file_id.write( '            while read FILE2; do\n')
            script_file_id.write( '                cat $FILE2 $FILE1 > $STEP_1_DIR/filtered_`basename "$FILE1"`\n')
            script_file_id.write( '                RC=$?\n')
            script_file_id.write( '                if [ $RC -ne 0 ]; then manage_error cat $RC; fi\n')
            script_file_id.write( '            done < header_files.txt\n')
            script_file_id.write( '        done < bi_allelic_files.txt\n')
            script_file_id.write( '        echo "Headers are pasted."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_1_4_move_files\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_1_4_move_files.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 1.4: Moving the bi-allelic, multi-allelic and filtered VCF files to output directories ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        mv $STEP_1_DIR/filtered*.vcf $STEP_1_FILTERED_DIR\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '        mv $STEP_1_DIR/*bi_allelic.vcf $STEP_1_BIALLELIC_DIR\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '        mv $STEP_1_DIR/*multi_allelic.vcf $STEP_1_MULTIALLELIC_DIR\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '        echo "Files are moved."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_1_5_plot_snp_data\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_1_5_plot_snp_data.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 1.5: Plotting the percentage of multi-allelic SNPs and the number of loci/SNPs post multi-allelic filtered ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_r_anaconda_code()}\n')
            script_file_id.write( '        Rscript --vanilla --verbose plot_percentage_multi_allelicSNPs.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error Rscript $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Plots are done."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_2_quantify_read_depth_and_error_rates\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_2_quantify_read_depth_and_error_rates.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 2: Quantifying and plotting read-depth and locus/SNP error rates ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_r_anaconda_code()}\n')
            script_file_id.write( '        Rscript --vanilla --verbose 2.Depth_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error Rscript $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        cho "Quantification and plots are done."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_3_reorganize_files\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_3_reorganize_files.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 3: Reorganizing files ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        INDIVIDUALS=""\n')
            script_file_id.write( '        while read LINE; do\n')
            script_file_id.write( '            if [ "$INDIVIDUALS" == "" ]; then \n')
            script_file_id.write( '                INDIVIDUALS=$LINE\n')
            script_file_id.write( '            else\n')
            script_file_id.write( '                INDIVIDUALS=$INDIVIDUALS","$LINE\n')
            script_file_id.write( '            fi\n')
            script_file_id.write( '        done < $SAMPLE_FILE\n')
            script_file_id.write(f'        source activate {xlib.get_vcftools_perl_libraries_anaconda_code()}\n')
            script_file_id.write( '        ls $STEP_1_FILTERED_DIR/*.vcf > vcf_files.txt\n')
            script_file_id.write( '        while read FILE; do\n')
            script_file_id.write( '            vcf-contrast +"$INDIVIDUALS" -"$INDIVIDUALS" $FILE | vcf-query -f "%CHROM\\t%POS[\\t%GTR]\\n" > $STEP_3_TXT_DIR/vcf_contrast_`basename "$FILE"`.txt\n')
            script_file_id.write( '            if [ $? -ne 0 ]; then manage_error perl-vcftools-vcf $RC; fi\n')
            script_file_id.write( '        done < vcf_files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Reorganization is done."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function step_4_filter_rates\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/step_4_filter_rates.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Step 4: Filtering locus-error-rates and SNP-error-rates ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_r_anaconda_code()}\n')
            script_file_id.write( '        Rscript --vanilla --verbose 4.Filtering_locusSNPerror.R\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error Rscript $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Filter is done."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_tmp_files\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/remove_tmp_files.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Removing temporal files ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        rm $STEP_1_DIR/*header*.txt\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        rm header_files.txt\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        rm bi_allelic_files.txt\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        rm IpyRAD_output_vcf.txt\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        rm vcf_files.txt\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        echo "Files are removed."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
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
            process_name = f'{xlib.get_raddesigner_name()} process'
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
            script_file_id.write(f'copy_r_scripts\n')
            script_file_id.write(f'fix_r_scripts\n')
            script_file_id.write(f'step_1_1_save_headers\n')
            script_file_id.write(f'step_1_2_separate_snps\n')
            script_file_id.write(f'step_1_3_paste_headers\n')
            script_file_id.write(f'step_1_4_move_files\n')
            script_file_id.write(f'step_1_5_plot_snp_data\n')
            script_file_id.write(f'step_2_quantify_read_depth_and_error_rates\n')
            script_file_id.write(f'step_3_reorganize_files\n')
            script_file_id.write(f'step_4_filter_rates\n')
            script_file_id.write(f'remove_tmp_files\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_raddesigner_process_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_raddesigner_process_starter(current_run_dir):
    '''
    Build the starter of the current RADdesigner process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the RADdesigner process starter
    try:
        if not os.path.exists(os.path.dirname(get_raddesigner_process_starter())):
            os.makedirs(os.path.dirname(get_raddesigner_process_starter()))
        with open(get_raddesigner_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_raddesigner_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_raddesigner_process_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_raddesigner_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a RADdesigner process from the last step ended OK.
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

    # get the current run directory
    if OK:
        current_run_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, result_dataset_id)

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_raddesigner_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_raddesigner_process_starter()), log)

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

def get_raddesigner_config_file():
    '''
    Get the RADdesigner config file path.
    '''

    # assign the RADdesigner config file path
    raddesigner_config_file = f'{xlib.get_config_dir()}/{xlib.get_raddesigner_code()}-config.txt'

    # return the RADdesigner config file path
    return raddesigner_config_file

#-------------------------------------------------------------------------------

def get_sample_file():
    '''
    Get the RADdesigner sample file path.
    '''

    # assign the RADdesigner sample file path
    sample_file = f'{xlib.get_temp_dir()}/{xlib.get_raddesigner_code()}-samples.txt'

    # return the RADdesigner sample file path
    return sample_file

#-------------------------------------------------------------------------------

def get_raddesigner_process_script():
    '''
    Get the RADdesigner process script path in the local computer.
    '''

    # assign the RADdesigner script path
    raddesigner_process_script = f'{xlib.get_temp_dir()}/{xlib.get_raddesigner_code()}-process.sh'

    # return the RADdesigner script path
    return raddesigner_process_script

#-------------------------------------------------------------------------------

def get_raddesigner_process_starter():
    '''
    Get the RADdesigner process starter path in the local computer.
    '''

    # assign the RADdesigner process starter path
    raddesigner_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_raddesigner_code()}-process-starter.sh'

    # return the RADdesigner starter path
    return raddesigner_process_starter

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the RADdesigner process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
