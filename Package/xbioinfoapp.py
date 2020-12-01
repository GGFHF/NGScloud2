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
                error_list.append(f'{error}\n')
                OK = False

    # check the Miniconda3 directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), xlib.get_miniconda3_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Miniconda3 installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_miniconda3_installation_script()} ...\n')
        (OK, error_list) = build_miniconda3_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Miniconda3 installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_miniconda3_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_miniconda3_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_miniconda3_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Miniconda3 installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_miniconda3_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_miniconda3_installation_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Miniconda3 installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_miniconda3_installation_starter()} ...\n')
        (OK, error_list) = build_miniconda3_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Miniconda3 installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_miniconda3_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_miniconda3_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_miniconda3_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Miniconda3 installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_miniconda3_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_miniconda3_installation_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Miniconda3 installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_miniconda3_installation_starter())} ...\n')
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
    (_, miniconda3_url, _) = xconfiguration.get_bioinfo_app_data(xlib.get_miniconda3_name())

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
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
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
            script_file_id.write( 'function remove_miniconda3_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_miniconda3_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_miniconda3_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_miniconda3_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_miniconda3_package\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_miniconda3_name()} installation package ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    wget --quiet --output-document {xlib.get_miniconda3_name()}.sh {miniconda3_url}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
            script_file_id.write( '    echo "The package is downloaded."\n')
            script_file_id.write(f'    chmod u+x {xlib.get_miniconda3_name()}.sh\n')
            script_file_id.write( '    echo "The run permision is set on."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_miniconda3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Installing {xlib.get_miniconda3_name()} to create Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    ./{xlib.get_miniconda3_name()}.sh -b -p {xlib.get_miniconda3_name()}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write(f'    if [ $RC -ne 0 ]; then manage_error {xlib.get_miniconda3_name()} $RC; fi\n'.format(''.format()))
            script_file_id.write( '    echo "Python 3 environment is created."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_miniconda3_package\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_miniconda3_name()} installation package ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_miniconda3_name()}.sh\n')
            script_file_id.write( '    echo "The package is removed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_gffutils_python3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing package gffutils in Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '    ./pip install --quiet gffutils\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error pip $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_joblib_python3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing package joblib in Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '    ./conda install --quiet --yes joblib\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error pip $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_matplotlib_python3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing package matplotlib in Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '    ./conda install --quiet --yes matplotlib\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error pip $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_biopython_python3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing package biopython in Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '    ./conda install --quiet --yes biopython\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error pip $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_requests_python3\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing package requests in Python 3 environment ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( '    ./conda install --quiet --yes requests\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error pip $RC; fi\n')
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
            process_name = f'{xlib.get_miniconda3_name()} installation'
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
            script_file_id.write( 'remove_miniconda3_directory\n')
            script_file_id.write( 'download_miniconda3_package\n')
            script_file_id.write( 'install_miniconda3\n')
            script_file_id.write( 'remove_miniconda3_package\n')
            script_file_id.write( 'install_gffutils_python3\n')
            script_file_id.write( 'install_joblib_python3\n')
            script_file_id.write( 'install_matplotlib_python3\n')
            script_file_id.write( 'install_biopython_python3\n')
            script_file_id.write( 'install_requests_python3\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_miniconda3_installation_script()} can not be created')
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
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_miniconda3_installation_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_miniconda3_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_miniconda3_installation_script():
    '''
    Get the Miniconda3 installation path in the local computer.
    '''

    # assign the Miniconda3 installation path
    miniconda3_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_miniconda3_name()}-installation.sh'

    # return the Miniconda3 installation path
    return miniconda3_installation_script

#-------------------------------------------------------------------------------

def get_miniconda3_installation_starter():
    '''
    Get the Miniconda3 installation starter path in the local computer.
    '''

    # assign the Miniconda3 installation starter path
    miniconda3_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_miniconda3_name()}-installation-starter.sh'

    # return the Miniconda3 installation starter path
    return miniconda3_installation_starter

#-------------------------------------------------------------------------------

def is_installed_anaconda_package(package_code, cluster_name, passed_connection, ssh_client):
    '''
    Check if a Anaconda package is installed.
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

    # check the Anaconda package directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{package_code} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
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

def install_anaconda_package_list(app_code, app_name, package_list, cluster_name, log, function=None):
    '''
    Install the Anaconda package list in the cluster.
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

    # warn that Anaconda package installation requirements are being verified
    if OK: 
        log.write(f'{xlib.get_separator()}\n')
        package_list_text = str(package_list).strip('[]').replace('\'','')
        log.write(f'Checking the Anaconda package list ({package_list_text}) installation requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the app directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
        else:
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write(f'You have to link a volume in the mounting point {xlib.get_cluster_app_dir()} for the cluster {cluster_name}.\n')
            OK = False

    # check the Miniconda3 installation
    if OK:
        (OK, error_list, is_installed) = is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_miniconda3_name()} is not installed.\n')
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
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), app_code)
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Anaconda package installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_anaconda_package_installation_script()} ...\n')
        (OK, error_list) = build_anaconda_package_installation_script(app_name, package_list, cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Anaconda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_anaconda_package_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_anaconda_package_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_anaconda_package_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Anaconda package installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_anaconda_package_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_anaconda_package_installation_script())}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Anaconda package installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_anaconda_package_installation_starter()} ...\n')
        (OK, error_list) = build_anaconda_package_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Anaconda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_anaconda_package_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_anaconda_package_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_anaconda_package_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Anaconda package installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_anaconda_package_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_anaconda_package_installation_starter())}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Anaconda package installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_anaconda_package_installation_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_anaconda_package_installation_starter()), log)

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

def build_anaconda_package_installation_script(app_name, package_list, cluster_name, current_run_dir):
    '''
    Build the Anaconda package installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Anaconda package installation script
    try:
        if not os.path.exists(os.path.dirname(get_anaconda_package_installation_script())):
            os.makedirs(os.path.dirname(get_anaconda_package_installation_script()))
        with open(get_anaconda_package_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
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
            for package in package_list:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function remove_anaconda_package_{package[0]}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "Removing {xlib.get_anaconda_name()} package {package[0]} ..."\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write(f'    ./conda env remove --yes --quiet --name {package[0]}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -eq 0 ]; then\n')
                script_file_id.write( '        echo "The old package is removed."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "The old package is not found."\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function install_anaconda_package_{package[0]}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                channel_list = package[2].split(';')
                channel_param_list_text = ''
                for channel in channel_list:
                    channel_param_list_text = f'{channel_param_list_text} --channel {channel}'
                if package[1] == 'last':
                    script_file_id.write(f'    echo "Installing {xlib.get_anaconda_name()} package {package[0]} - last version ..."\n')
                    script_file_id.write(f'    ./conda create --yes --quiet {channel_param_list_text} --name {package[0]} {package[0]}\n')
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                else:
                    script_file_id.write(f'    echo "Installing {xlib.get_anaconda_name()} package {package[0]} - version {package[1]} ..."\n')
                    script_file_id.write(f'    ./conda create --yes --quiet {channel_param_list_text} --name {package[0]} {package[0]}={package[1]}\n')
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '    if [ $RC -ne 0 ]; then\n')
                    script_file_id.write(f'        echo "Installing {xlib.get_anaconda_name()} package {package[0]} - last version ..."\n')
                    script_file_id.write(f'        ./conda create --yes --quiet {channel_param_list_text} --name {package[0]} {package[0]}\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                    script_file_id.write( '    fi\n')
                script_file_id.write(f'    echo "The package is installed."\n')
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
            process_name = f'{app_name} installation'
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
            for package in package_list:
                script_file_id.write(f'remove_anaconda_package_{package[0]}\n')
                script_file_id.write(f'install_anaconda_package_{package[0]}\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_anaconda_package_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_anaconda_package_installation_starter(current_run_dir):
    '''
    Build the starter of the Anaconda package installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Anaonda package installation starter
    try:
        if not os.path.exists(os.path.dirname(get_anaconda_package_installation_starter())):
            os.makedirs(os.path.dirname(get_anaconda_package_installation_starter()))
        with open(get_anaconda_package_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_anaconda_package_installation_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_anaconda_package_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_anaconda_package_installation_script():
    '''
    Get the Anaconda package installation path in the local computer.
    '''

    # assign the Anaconda package installation path
    anaconda_package_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_anaconda_name()}-installation.sh'

    # return the Anaconda package sinstallation path
    return anaconda_package_installation_script

#-------------------------------------------------------------------------------

def get_anaconda_package_installation_starter():
    '''
    Get the Anaconda package installation starter path in the local computer.
    '''

    # assign the Anaconda package installation starter path
    anaconda_package_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_anaconda_name()}-installation-starter.sh'

    # return the Anaconda package installation starter path
    return anaconda_package_installation_starter

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
                error_list.append(f'{error}\n')
                OK = False

    # check the Anaconda package directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_r_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
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
    Install the Anaconda package list in the cluster.
    '''

    # initialize the control variable
    OK = True

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
        log.write(f'Checking the R and analysis packages installation requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the app directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] == 'RC=0':
            OK = True
        else:
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write(f'You have to link a volume in the mounting point {xlib.get_cluster_app_dir()} for the cluster {cluster_name}.\n')
            OK = False

    # check the Miniconda3 installation
    if OK:
        (OK, error_list, is_installed) = is_installed_miniconda3(cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_miniconda3_name()} is not installed.\n')
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
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), xlib.get_r_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the R installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_r_installation_script()} ...\n')
        (OK, error_list) = build_r_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the R installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_r_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_r_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_r_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the R installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_r_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_r_installation_script())}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the R installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_r_installation_starter()} ...\n')
        (OK, error_list) = build_r_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the R installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_r_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_r_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_r_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the R installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_r_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_r_installation_starter())}'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the R installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_r_installation_starter())} ...\n')
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

def build_r_installation_script(cluster_name, current_run_dir):
    '''
    Build the R installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the R software dictionary
    r_software_dict = xconfiguration.get_r_software_dict()

    # get R software dictionary list in Anaconda and GitHub
    anaconda_r_software_list = []
    github_r_software_list = []
    for package_code in r_software_dict.keys():
        if package_code != 'r-essentials' and r_software_dict[package_code]['url'] == xlib.get_anaconda_url():
            anaconda_r_software_list.append(package_code)
        elif package_code != 'r-essentials' and r_software_dict[package_code]['url'] == xlib.get_github_url():
            github_r_software_list.append(package_code)

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
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'GITHUB_INSTALLATION_SCRIPT={xlib.get_cluster_app_dir()}/install_github_packages.R\n')
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
            script_file_id.write( 'function remove_r\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_r_name()} ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'    ./conda env remove --yes --quiet --name {xlib.get_r_name()}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -eq 0 ]; then\n')
            script_file_id.write( '        echo "The old package is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '      echo "The old package is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            (r_software_version, _, r_software_channels) = xconfiguration.get_r_software_data('r-essentials')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_package_r-essentials\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            channel_list = r_software_channels.split(';')
            channel_param_list_text = ''
            for channel in channel_list:
                channel_param_list_text = f'{channel_param_list_text} --channel {channel}'
            if r_software_version == 'last':
                script_file_id.write(f'    echo "Installing package r-essential - last version - channel {r_software_channels} ..."\n')
                script_file_id.write(f'    ./conda create --yes --quiet {channel_param_list_text} --name {xlib.get_r_name()} r-essentials\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
            else:
                script_file_id.write(f'    echo "Installing package r-essential - version {r_software_version} - channel {r_software_channels} ..."\n')
                script_file_id.write(f'    ./conda create --yes --quiet {channel_param_list_text} --name {xlib.get_r_name()} r-essentials={r_software_version}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then\n')
                script_file_id.write(f'        echo "Installing package r-essential - last version - channel {r_software_channels} ..."\n')
                script_file_id.write(f'        ./conda create --yes --quiet {channel_param_list_text} --name {xlib.get_r_name()} r-essentials\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                script_file_id.write( '    fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            for package_code in anaconda_r_software_list:
                (r_software_version, _, r_software_channels) = xconfiguration.get_r_software_data(package_code)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function install_package_{package_code}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
                script_file_id.write(f'    source activate {xlib.get_r_name()}\n')
                channel_list = r_software_channels.split(';')
                channel_param_list_text = ''
                for channel in channel_list:
                    channel_param_list_text = f'{channel_param_list_text} --channel {channel}'
                if r_software_version == 'last':
                    script_file_id.write(f'    echo "Installing package {package_code} - last version - channel {r_software_channels} ..."\n')
                    script_file_id.write(f'    ./conda install --yes --quiet {channel_param_list_text} {package_code}\n')
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                else:
                    script_file_id.write(f'    echo "Installing package {package_code} - version {r_software_version} - channel {r_software_channels} ..."\n')
                    script_file_id.write(f'    ./conda install --yes --quiet {channel_param_list_text} {package_code}={r_software_version}\n')
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '    if [ $RC -ne 0 ]; then\n')
                    script_file_id.write(f'        echo "Installing package {package_code} - last version - channel {r_software_channels} ..."\n')
                    script_file_id.write(f'        ./conda create --yes --quiet {channel_param_list_text} --name {xlib.get_r_name()} {package_code}\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error conda $RC; fi\n')
                    script_file_id.write( '    fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "The package is installed."\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function build_github_installation_script\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Creating the R script to install GitHuB packages..."\n')
            script_file_id.write(f'    echo "library(\\"devtools\\")" >$GITHUB_INSTALLATION_SCRIPT \n')
            for package_code in github_r_software_list:
                (r_software_version, _, r_software_channel) = xconfiguration.get_r_software_data(package_code)
                script_file_id.write(f'    echo "print(\\"Installing package {package_code} - GitHub: {r_software_channel} ...\\")" >>$GITHUB_INSTALLATION_SCRIPT \n')
                script_file_id.write(f'    echo "install_github(\\"{r_software_channel}\\")" >>$GITHUB_INSTALLATION_SCRIPT \n')
            script_file_id.write( '    echo "The R script is buit."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_github_installation_script\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running the R script to install GitHuB packages..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'    source activate {xlib.get_r_name()}\n')
            script_file_id.write( '    Rscript --vanilla --verbose $GITHUB_INSTALLATION_SCRIPT\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error Rscript $RC; fi\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '    echo "The R script is run."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_github_installation_script\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Removing the R script to install GitHuB packages..."\n')
            script_file_id.write( '    rm $GITHUB_INSTALLATION_SCRIPT\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "The R script is removed."\n')
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
            process_name = f'{xlib.get_r_name()} and Analysis Packages installation'
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
            script_file_id.write( 'remove_r\n')
            script_file_id.write( 'install_package_r-essentials\n')
            for package_code in anaconda_r_software_list:
                script_file_id.write(f'install_package_{package_code}\n')
            script_file_id.write( 'build_github_installation_script\n')
            script_file_id.write( 'run_github_installation_script\n')
            script_file_id.write( 'remove_github_installation_script\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_r_installation_script()} can not be created')
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
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_r_installation_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_r_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_r_installation_script():
    '''
    Get the R installation path in the local computer.
    '''

    # assign the R installation path
    r_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_r_name()}-installation.sh'

    # return the R installation path
    return r_installation_script

#-------------------------------------------------------------------------------

def get_r_installation_starter():
    '''
    Get the R installation starter path in the local computer.
    '''

    # assign the R installation starter path
    r_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_r_name()}-installation-starter.sh'

    # return the R installation starter path
    return r_installation_starter

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to BioInfo applications used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
