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
This file contains functions related to the ddRADseqTools process used in both
console mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys
import urllib

import xbioinfoapp
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def is_installed_ddradseqtools(cluster_name, passed_connection, ssh_client):
    '''
    Check if ddRADseqTools is installed.
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

    # check the ddRADseqTools directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()} ] && echo RC=0 || echo RC=1'
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

def install_ddradseqtools(cluster_name, log, function=None):
    '''
    Install the ddRADseqTools software in the cluster.
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

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), xlib.get_ddradseqtools_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseqTools installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_ddradseqtools_installation_script()} ...\n')
        (OK, error_list) = build_ddradseqtools_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the ddRADseqTools installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_ddradseqtools_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ddradseqtools_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseqtools_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseqTools installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ddradseqtools_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ddradseqtools_installation_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseqTools installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_ddradseqtools_installation_starter()} ...\n')
        (OK, error_list) = build_ddradseqtools_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the ddRADseqTools installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_ddradseqtools_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ddradseqtools_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseqtools_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseqTools installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ddradseqtools_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ddradseqtools_installation_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the ddRADseqTools installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ddradseqtools_installation_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ddradseqtools_installation_starter()), log)

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

def build_ddradseqtools_installation_script(cluster_name, current_run_dir):
    '''
    Build the ddRADseqTools installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of ddRADseqTools
    (ddradseqtools_version, ddradseqtools_url, ddradseqtools_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_ddradseqtools_name())

    # write the ddRADseqTools installation script
    try:
        if not os.path.exists(os.path.dirname(get_ddradseqtools_installation_script())):
            os.makedirs(os.path.dirname(get_ddradseqtools_installation_script()))
        with open(get_ddradseqtools_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function remove_ddradseqtools_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_ddradseqtools_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_ddradseqtools_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_ddradseqtools_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_ddradseqtools_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_ddradseqtools_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            download_script = f'import requests; r = requests.get(\'{ddradseqtools_url}\') ; open(\'{xlib.get_ddradseqtools_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $MINICONDA3_BIN_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write( '    echo "The file is downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function decompress_ddradseqtools_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Decompressing the {xlib.get_ddradseqtools_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    unzip -u {xlib.get_ddradseqtools_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write( '    echo "The file is decompressed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_ddradseqtools_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Renaming the {xlib.get_ddradseqtools_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    mv {xlib.get_ddradseqtools_name()}-master {xlib.get_ddradseqtools_name()}\n')
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
            script_file_id.write(f'    chmod u+x {xlib.get_ddradseqtools_name()}/Package/*.py\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi\n')
            script_file_id.write( '    echo "Permissions are set."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_ddradseqtools_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_ddradseqtools_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_ddradseqtools_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "The file is removed."\n')
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
            process_name = f'{xlib.get_ddradseqtools_name()} installation'
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
            script_file_id.write( 'remove_ddradseqtools_directory\n')
            script_file_id.write( 'download_ddradseqtools_installation_file\n')
            script_file_id.write( 'decompress_ddradseqtools_installation_file\n')
            script_file_id.write( 'rename_ddradseqtools_directory\n')
            script_file_id.write( 'set_execution_permissions\n')
            script_file_id.write( 'remove_ddradseqtools_installation_file\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ddradseqtools_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ddradseqtools_installation_starter(current_run_dir):
    '''
    Build the starter of the ddRADseqTools installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the ddRADseqTools installation starter
    try:
        if not os.path.exists(os.path.dirname(get_ddradseqtools_installation_starter())):
            os.makedirs(os.path.dirname(get_ddradseqtools_installation_starter()))
        with open(get_ddradseqtools_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_ddradseqtools_installation_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ddradseqtools_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ddradseqtools_installation_script():
    '''
    Get the ddRADseqTools installation path in the local computer.
    '''

    # assign the ddRADseqTools installation path
    ddradseqtools_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_ddradseqtools_name()}-installation.sh'

    # return the ddRADseqTools installation path
    return ddradseqtools_installation_script

#-------------------------------------------------------------------------------

def get_ddradseqtools_installation_starter():
    '''
    Get the ddRADseqTools installation starter path in the local computer.
    '''

    # assign the ddRADseqTools installation starter path
    ddradseqtools_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_ddradseqtools_name()}-installation-starter.sh'

    # return the ddRADseqTools installation starter path
    return ddradseqtools_installation_starter

#-------------------------------------------------------------------------------

def create_restriction_site_file():
    '''
    Create the file of restriction sites with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of restriction sites and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_restriction_site_file())):
            os.makedirs(os.path.dirname(get_restriction_site_file()))
        with open(get_restriction_site_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains the restriction sites recognition motifs and their cut sites\n')
            file_id.write( '# (a cut site is represented by * in the sequence).\n')
            file_id.write( '# New enzymes can be included at the end of the file\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: enzyme_id;restriction_site_seq(5\'->3\')\n')
            file_id.write( '\n')
            file_id.write( 'AatII;GACGT*C\n')
            file_id.write( 'Acc65I;G*GTACC\n')
            file_id.write( 'AclI;AA*CGTT\n')
            file_id.write( 'AatII;GACGT*C\n')
            file_id.write( 'Acc65I;G*GTACC\n')
            file_id.write( 'AclI;AA*CGTT\n')
            file_id.write( 'AfeI;AGC*GCT\n')
            file_id.write( 'AflII;C*TTAAG\n')
            file_id.write( 'AgeI;A*CCGGT\n')
            file_id.write( 'ApaI;GGGCC*C\n')
            file_id.write( 'ApaLI;G*TGCAC\n')
            file_id.write( 'AscI;GG*CGCGCC\n')
            file_id.write( 'AseI;AT*TAAT\n')
            file_id.write( 'AsiSI;GCGAT*CGC\n')
            file_id.write( 'AvrII;C*CTAGG\n')
            file_id.write( 'BamHI;G*GATCC\n')
            file_id.write( 'BclI;T*GATCA\n')
            file_id.write( 'BglII;A*GATCT\n')
            file_id.write( 'BmtI;GCTAG*C\n')
            file_id.write( 'BsiWI;C*GTACG\n')
            file_id.write( 'BspEI;T*CCGGA\n')
            file_id.write( 'BspHI;T*CATGA\n')
            file_id.write( 'BsrGI;T*GTACA\n')
            file_id.write( 'BssHII;G*CGCGC\n')
            file_id.write( 'BstBI;TT*CGAA\n')
            file_id.write( 'BstZ17I;GTA*TAC\n')
            file_id.write( 'ClaI;AT*CGAT\n')
            file_id.write( 'Csp6I;G*TAC\n')
            file_id.write( 'DraI;TTT*AAA\n')
            file_id.write( 'EagI;C*GGCCG\n')
            file_id.write( 'EcoRI;G*AATTC\n')
            file_id.write( 'EcoRV;GAT*ATC\n')
            file_id.write( 'FseI;GGCCGG*CC\n')
            file_id.write( 'FspI;TGC*GCA\n')
            file_id.write( 'HindIII;A*AGCTT\n')
            file_id.write( 'HpaI;GTT*AAC\n')
            file_id.write( 'KpnI;GGTAC*C\n')
            file_id.write( 'MfeI;C*AATTG\n')
            file_id.write( 'MluI;A*CGCGT\n')
            file_id.write( 'MscI;TGG*CCA\n')
            file_id.write( 'MseI;T*TAA\n')
            file_id.write( 'NaeI;GCC*GGC\n')
            file_id.write( 'NarI;GG*CGCC\n')
            file_id.write( 'NcoI;C*CATGG\n')
            file_id.write( 'NdeI;CA*TATG\n')
            file_id.write( 'NgoMIV;G*CCGGC\n')
            file_id.write( 'NheI;G*CTAGC\n')
            file_id.write( 'NlaIII;CATG*\n')
            file_id.write( 'NotI;GC*GGCCGC\n')
            file_id.write( 'NruI;TCG*CGA\n')
            file_id.write( 'NsiI;ATGCA*T\n')
            file_id.write( 'PacI;TTAAT*TAA\n')
            file_id.write( 'PciI;A*CATGT\n')
            file_id.write( 'PmeI;GTTT*AAAC\n')
            file_id.write( 'PmlI;CAC*GTG\n')
            file_id.write( 'PsiI;TTA*TAA\n')
            file_id.write( 'PspOMI;G*GGCCC\n')
            file_id.write( 'PstI;CTGCA*G\n')
            file_id.write( 'PvuI;CGAT*CG\n')
            file_id.write( 'PvuII;CAG*CTG\n')
            file_id.write( 'SacI;GAGCT*C\n')
            file_id.write( 'SacII;CCGC*GG\n')
            file_id.write( 'SalI;G*TCGAC\n')
            file_id.write( 'SbfI;CCTGCA*GG\n')
            file_id.write( 'ScaI;AGT*ACT\n')
            file_id.write( 'SfoI;GGC*GCC\n')
            file_id.write( 'SmaI;CCC*GGG\n')
            file_id.write( 'SnaBI;TAC*GTA\n')
            file_id.write( 'SpeI;A*CTAGT\n')
            file_id.write( 'SphI;GCATG*C\n')
            file_id.write( 'SspI;AAT*ATT\n')
            file_id.write( 'StuI;AGG*CCT\n')
            file_id.write( 'SwaI;ATTT*AAAT\n')
            file_id.write( 'XbaI;T*CTAGA\n')
            file_id.write( 'XhoI;C*TCGAG\n')
            file_id.write( 'XmaI;C*CCGGG\n')
            file_id.write( '\n')
            file_id.write( 'AccI;GT*MKAC\n')
            file_id.write( 'AccB1I;G*GYRCC\n')
            file_id.write( 'AccB2I;RGCGC*Y\n')
            file_id.write( 'AceI;G*CWGC\n')
            file_id.write( 'AdeI;CACNNN*GTG\n')
            file_id.write( 'MslI;CAYNN*NNRTG\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_restriction_site_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_restriction_site_file(strict):
    '''
    Check the file of restriction sites.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the record of restriction site file (record format: enzyme_id;restriction_site_seq)
    pattern = r'^(.*);(.*)$'

    # open the file of restriction sites
    try:
        restriction_site_file_id = open(get_restriction_site_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_restriction_site_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = restriction_site_file_id.readline()

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try: 
                    mo = re.search(pattern, record)
                    enzyme_id = mo.group(1).strip()
                    restriction_site_seq = mo.group(2).strip()
                except Exception as e:
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: There is a format error in the record: {record}')
                    OK = False
                    break

                # check that the restriction site sequence is correct
                if not xlib.is_valid_sequence(restriction_site_seq, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The sequence is invalid in the record: {record}')
                    OK = False
                    break

                # check that the enzyme identification is formed by alphanumeric characters
                if not xlib.is_name_valid(enzyme_id):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The enzyme id has characters non-alphanumeric in the record: {record}')
                    OK = False
                    break

            # read the next record
            record = restriction_site_file_id.readline()

    # close the file of restriction sites
    restriction_site_file_id.close()

    # warn that the file of restriction sites is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_restriction_site_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_restriction_site_file():
    '''
    Get the restriction site file path.
    '''

    # assign the restriction site file path
    restriction_site_file = f'{xlib.get_config_dir()}/restrictionsites.txt'

    # return the restriction site file path
    return restriction_site_file

#-------------------------------------------------------------------------------

def get_restriction_enzyme_dict():
    '''
    Get a dictionary of restriction enzymes with their rectriction sites.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of restriction enzymes
    restriction_enzyme_dict = {}

    # set the pattern of the restriction site file record (record format: enzyme_id;restriction_site_seq)
    pattern = r'^(.*);(.*)$'

    # check the restriction site file
    (OK, error_list) = check_restriction_site_file(strict=True)

    # read every record of restriction site file
    if OK:
        with open(get_restriction_site_file(), mode='r', encoding='iso-8859-1', newline='\n') as file_id:
            for record in file_id:

                # if the record is not a comment nor a line with blank characters
                if not record.lstrip().startswith('#') and record.strip() != '':

                    # extract the data and add enzyme data to the dictionary of restriction enzymes
                    mo = re.search(pattern, record)
                    enzyme_id = mo.group(1).strip()
                    restriction_site_seq = mo.group(2).strip()

                    # add data to the dictionary of restriction enzymes
                    enzyme_id_seq = f'{enzyme_id} ({restriction_site_seq})'
                    restriction_enzyme_dict[enzyme_id] = {'enzyme_id': enzyme_id, 'restriction_site_seq': restriction_site_seq, 'enzyme_id_seq': enzyme_id_seq}

    # return the control variable, the error list and the dictionary of restriction enzymes
    return (OK, error_list, restriction_enzyme_dict)

#-------------------------------------------------------------------------------

def get_enzyme_id_seq_list(restriction_enzyme_dict):
    '''
    Get a list of enzyme identification and restriction site seq.
    '''

    # initialize the list of the enzyme identification and restriction site seq
    enzyme_id_seq_list = []

    # build the list of the enzyme identification and restriction site seq
    for enzyme_id in restriction_enzyme_dict.keys():
        enzyme_id_seq_list.append(restriction_enzyme_dict[enzyme_id]['enzyme_id_seq'])

    # sort the list of the enzyme identification and restriction site seq
    if enzyme_id_seq_list != []:
        enzyme_id_seq_list.sort()

    # return the list of the enzyme identification and restriction site seq
    return enzyme_id_seq_list

#-------------------------------------------------------------------------------

def get_enzyme_id(enzyme_id_seq, restriction_enzyme_dict):
    '''
    Get the enzyme identification from the enzyme identification and restriction site seq.
    '''

    # initialize the control variable
    enzyme_id_found = None

    # search the result dataset identification
    for enzyme_id in restriction_enzyme_dict.keys():
        if restriction_enzyme_dict[enzyme_id]['enzyme_id_seq'] == enzyme_id_seq:
            enzyme_id_found = enzyme_id
            break

    # return the enzyme identification
    return enzyme_id_found

#-------------------------------------------------------------------------------

def is_restriction_enzyme_code(enzyme_id, restriction_enzyme_dict):
    '''
    Check if an enzyme identification is in the dictionary of restriction enzymes.
    '''

    # initialize the control variable
    OK = True

    try:
        restriction_site_seq = restriction_enzyme_dict[enzyme_id]
    except:
        OK = False

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def create_end_file():
    '''
    Create the file of ends with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of ends and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_end_file())):
            os.makedirs(os.path.dirname(get_end_file()))
        with open(get_end_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains end sequences integrated by indexes, degenerate nucleotides to indentify\n')
            file_id.write( '# the PCR duplicates (DBR), adapter and primer. A read has two ends: one where the adapter 1\n')
            file_id.write( '#  is and another where adapter 2 is.\n')
            file_id.write( '#\n')
            file_id.write( '# IND1_IND2_DBR technique: The sequence of the end corresponding to the adapter 1 includes\n')
            file_id.write( '# a index1 sequence (111..., each digit 1 represents a nucleotide of the index1), the sequence\n')
            file_id.write( '# of the end corresponing to the adapter 2 include a index2 sequence (222..., each digit 2\n')
            file_id.write( '# represents a nucleotide of the index2). A DBR sequence (333..., each digit 3 represents a\n')
            file_id.write( '# nucleotide of the DBR) has to be included in the end sequence of adapter 1 or adapter2.\n')
            file_id.write( '#\n')
            file_id.write( '# IND1_IND2 technique: The sequence of the end corresponding to the adapter 1 includes\n')
            file_id.write( '# a index1 sequence (111...), the sequence of the end corresponing to the adapter 2 include a\n')
            file_id.write( '# index2 sequence (222...). The DBR sequence (333...) is not considered.\n')
            file_id.write( '#\n')
            file_id.write( '# IND1_DBR technique: The sequence of the end corresponding to the adapter 1 includes a index1\n')
            file_id.write( '# sequence (111...) and a DBR sequence (333...).  The index2 sequence (222...) is not\n')
            file_id.write( '# considered.\n')
            file_id.write( '#\n')
            file_id.write( '# IND1 technique: The sequence of the end corresponding to the adapter 1 includes a index1\n')
            file_id.write( '# sequence (111...). The index2 sequence (222...) and the DBR sequence (333...) are not\n')
            file_id.write( '# considered.\n')
            file_id.write( '#\n')
            file_id.write( '# The length of index1 (111...), index2 (222...) and DBR (333...) of ends used by a program\n')
            file_id.write( '# have to be equal to the value of the options index1len, index2len and dbrlen received by\n')
            file_id.write( '# that program.\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: end_id;end_seq(5\'->3\')\n')
            file_id.write( '\n')
            file_id.write( 'end01;AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT111111C\n')
            file_id.write( 'end02;CAAGCAGAAGACGGCATACGAGAT3333222222GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC\n')
            file_id.write( 'end03;CAAGCAGAAGACGGCATACGAGAT111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC\n')
            file_id.write( 'end04;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC\n')
            file_id.write( 'end05;CAAGCAGAAGACGGCATACGAGAT3333111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC\n')
            file_id.write( 'end06;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC\n')
            file_id.write( 'end07;CAAGCAGAAGACGGCATACGAGAT3333111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC\n')
            file_id.write( 'end08;AATGATACGGCGACCACCGAGATCTACAC222222ACACTCTTTCCCTACACGACGCTCTTCCGATC\n')
            file_id.write( 'end09;CAAGCAGAAGACGGCATACGAGAT111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC\n')
            file_id.write( 'end10;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC\n')
            file_id.write( '\n')
            file_id.write( 'end51;111111\n')
            file_id.write( 'end52;\n')
            file_id.write( '\n')
            file_id.write( 'end61;3333111111\n')
            file_id.write( 'end62;\n')
            file_id.write( 'end63;33331111111\n')
            file_id.write( 'end64;\n')
            file_id.write( 'end65;333331111111\n')
            file_id.write( 'end66;\n')
            file_id.write( '\n')
            file_id.write( 'end71;111111\n')
            file_id.write( 'end72;222222\n')
            file_id.write( 'end73;11111\n')
            file_id.write( 'end74;22222\n')
            file_id.write( '\n')
            file_id.write( 'end81;111111\n')
            file_id.write( 'end82;3333222222\n')
            file_id.write( '\n')
            file_id.write( 'end91;3333111111\n')
            file_id.write( 'end92;222222\n')
            file_id.write( 'end93;33331111111\n')
            file_id.write( 'end94;2222222\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_end_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_end_file(strict):
    '''
    Check the file of ends.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the record of end file (record format: end_id;end_seq)
    pattern = r'^(.*);(.*)$'

    # open the file of ends
    try:
        end_file_id = open(get_end_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_end_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = end_file_id.readline()

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try: 
                    mo = re.search(pattern, record)
                    end_id = mo.group(1).strip()
                    end_seq = mo.group(2).strip()
                except Exception as e:
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: There is a format error in the record: {record}')
                    OK = False
                    break

                # check that the end sequence is correct
                if not xlib.is_valid_sequence(end_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=['1', '2', '3'], cut_tag_check=False):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The sequence is invalid in the record: {record}')
                    OK = False
                    break

                # check that the end identification is formed by alphanumeric characters
                if not xlib.is_name_valid(end_id):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The end id has characters non-alphanumeric in the record: {record}')
                    OK = False
                    break

            # read the next record
            record = end_file_id.readline()

    # close the file of ends
    end_file_id.close()

    # warn that the file of ends is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_end_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_end_file():
    '''
    Get the end file path.
    '''

    # assign the end file path
    end_file = f'{xlib.get_config_dir()}/ends.txt'

    # return the end file path
    return end_file

#-------------------------------------------------------------------------------

def get_end_dict():
    '''
    Get a dictionary of ends with their sequences.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of ends
    end_dict = {}

    # set the pattern of the end file record (record format: end_id;end_seq)
    pattern = r'^(.*);(.*)$'

    # check the end file
    (OK, error_list) = check_end_file(strict=True)

    # read every record of end file
    if OK:
        with open(get_end_file(), mode='r', encoding='iso-8859-1', newline='\n') as file_id:
            for record in file_id:

                # if the record is not a comment nor a line with blank characters
                if not record.lstrip().startswith('#') and record.strip() != '':

                    # extract the data and add end data to the dictionary of ends
                    mo = re.search(pattern, record)
                    end_id = mo.group(1).strip()
                    end_seq = mo.group(2).strip()

                    # add data to the dictionary of end
                    end_dict[end_id] = {'end_id': end_id, 'end_seq': end_seq}

    # return the control variable, the error list and the dictionary of ends
    return (OK, error_list, end_dict)

#-------------------------------------------------------------------------------

def create_individual_file():
    '''
    Create the file of individuals with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of restriction sites and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_individual_file())):
            os.makedirs(os.path.dirname(get_individual_file()))
        with open(get_individual_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains the sequences of each invidual: index1 sequences that are attached to the\n')
            file_id.write( '# 5\' end of the Watson strand and, optionally, the index2 sequences that are attached to the\n')
            file_id.write( '# 5\' end of the Crick strand.\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: individual_id;replicated_individual_id or NONE;population_id;index1_seq(5\'->3\' in Watson strand);[index2_seq(5\'->3\' in Crick strand)]\n')
            file_id.write( '\n')
            file_id.write( '# These data corresponding to an example with 8 index1 (6bp) and 6 index2 (6bp) by index1 -> up 48 individuals.\n')
            file_id.write( '\n')
            file_id.write( 'ind0101;NONE;pop01;GGTCTT;ATCACG\n')
            file_id.write( 'ind0102;NONE;pop01;GGTCTT;CGATGT\n')
            file_id.write( 'ind0103;NONE;pop01;GGTCTT;TTAGGC\n')
            file_id.write( 'ind0104;NONE;pop01;GGTCTT;TGACCA\n')
            file_id.write( 'ind0105;NONE;pop01;GGTCTT;ACAGTG\n')
            file_id.write( 'ind0103r;ind0103;pop01;GGTCTT;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0201;NONE;pop01;CTGGTT;ATCACG\n')
            file_id.write( 'ind0202;NONE;pop01;CTGGTT;CGATGT\n')
            file_id.write( 'ind0203;NONE;pop01;CTGGTT;TTAGGC\n')
            file_id.write( 'ind0204;NONE;pop01;CTGGTT;TGACCA\n')
            file_id.write( 'ind0205;NONE;pop01;CTGGTT;ACAGTG\n')
            file_id.write( 'ind0204r;ind0204;pop01;CTGGTT;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0301;NONE;pop01;AAGATA;ATCACG\n')
            file_id.write( 'ind0302;NONE;pop01;AAGATA;CGATGT\n')
            file_id.write( 'ind0303;NONE;pop01;AAGATA;TTAGGC\n')
            file_id.write( 'ind0304;NONE;pop01;AAGATA;TGACCA\n')
            file_id.write( 'ind0305;NONE;pop01;AAGATA;ACAGTG\n')
            file_id.write( 'ind0306;NONE;pop01;AAGATA;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0401;NONE;pop01;ACTTCC;ATCACG\n')
            file_id.write( 'ind0402;NONE;pop01;ACTTCC;CGATGT\n')
            file_id.write( 'ind0403;NONE;pop01;ACTTCC;TTAGGC\n')
            file_id.write( 'ind0404;NONE;pop01;ACTTCC;TGACCA\n')
            file_id.write( 'ind0405;NONE;pop01;ACTTCC;ACAGTG\n')
            file_id.write( 'ind0406;NONE;pop01;ACTTCC;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0501;NONE;pop01;TTACGG;ATCACG\n')
            file_id.write( 'ind0502;NONE;pop01;TTACGG;CGATGT\n')
            file_id.write( 'ind0503;NONE;pop01;TTACGG;TTAGGC\n')
            file_id.write( 'ind0504;NONE;pop01;TTACGG;TGACCA\n')
            file_id.write( 'ind0505;NONE;pop01;TTACGG;ACAGTG\n')
            file_id.write( 'ind0506;NONE;pop01;TTACGG;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0601;NONE;pop01;AACGAA;ATCACG\n')
            file_id.write( 'ind0602;NONE;pop01;AACGAA;CGATGT\n')
            file_id.write( 'ind0603;NONE;pop01;AACGAA;TTAGGC\n')
            file_id.write( 'ind0604;NONE;pop01;AACGAA;TGACCA\n')
            file_id.write( 'ind0605;NONE;pop01;AACGAA;ACAGTG\n')
            file_id.write( 'ind0606;NONE;pop01;AACGAA;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0701;NONE;pop01;ATTCAT;ATCACG\n')
            file_id.write( 'ind0702;NONE;pop01;ATTCAT;CGATGT\n')
            file_id.write( 'ind0703;NONE;pop01;ATTCAT;TTAGGC\n')
            file_id.write( 'ind0704;NONE;pop01;ATTCAT;TGACCA\n')
            file_id.write( 'ind0705;NONE;pop01;ATTCAT;ACAGTG\n')
            file_id.write( 'ind0706;NONE;pop01;ATTCAT;GCCAAT\n')
            file_id.write( '\n')
            file_id.write( 'ind0801;NONE;pop01;CCGACC;ATCACG\n')
            file_id.write( 'ind0802;NONE;pop01;CCGACC;CGATGT\n')
            file_id.write( 'ind0803;NONE;pop01;CCGACC;TTAGGC\n')
            file_id.write( 'ind0804;NONE;pop01;CCGACC;TGACCA\n')
            file_id.write( 'ind0805;NONE;pop01;CCGACC;ACAGTG\n')
            file_id.write( 'ind0806;NONE;pop01;CCGACC;GCCAAT\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_individual_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_individual_file(strict):
    '''
    Check the file of individuals.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the record of individual file (record format: individual_id;replicated_individual_id;population_id;index1_seq;[index2_seq])
    pattern = r'^(.+);(.+);(.+);(.+);(.*)$'

    # open the file of individuals
    try:
        individual_file_id = open(get_individual_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_individual_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = individual_file_id.readline()
        is_there_index2 = None

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try:
                    mo = re.search(pattern, record)
                    individual_id = mo.group(1).strip()
                    replicated_individual_id = mo.group(2).strip()
                    population_id = mo.group(3).strip()
                    index1_seq = mo.group(4).strip()
                    index2_seq = mo.group(5).strip()
                except Exception as e:
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: There is a format error in the record: {record}')
                    OK = False
                    break

                # check that the index1 sequence is correct
                if not xlib.is_valid_sequence(index1_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The index 1 sequence is invalid in the record: {record}')
                    OK = False
                    break

                # check that the index2 sequence is correct
                if index2_seq != '':
                    if not xlib.is_valid_sequence(index2_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                        record = record.replace('\n', '')
                        error_list.append(f'*** ERROR: The index 2 sequence is invalid in the record: {record}')
                        OK = False
                        break
                    if is_there_index2 == None:
                        is_there_index2 = True
                    else:
                        if not is_there_index2:
                            error_list.append('*** ERROR: All the records have to have the same format.')
                            OK = False
                            break
                else:
                    if is_there_index2 == None:
                        is_there_index2 = False
                    else:
                        if is_there_index2:
                            error_list.append('*** ERROR: All the records have to have the same format.')
                            OK = False
                            break

                # check that the individual identification is formed by alphanumeric characters
                if not xlib.is_name_valid(individual_id):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The individual id has characters non-alphanumeric in the record: {record}')
                    OK = False
                    break

                # check that the replicated individual identification is formed by alphanumeric characters
                if not xlib.is_name_valid(replicated_individual_id):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The replicated individual id has characters non-alphanumeric in the record: {record}')
                    OK = False
                    break

                # check that the population identification is formed by alphanumeric characters
                if not xlib.is_name_valid(population_id):
                    record = record.replace('\n', '')
                    error_list.append(f'*** ERROR: The population id has characters non-alphanumeric in the record: {record}')
                    OK = False
                    break

            # read the next record
            record = individual_file_id.readline()

    # close the file of individuals
    individual_file_id.close()

    # warn that the file of individuals is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_individual_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_individual_file():
    '''
    Get the individual file path.
    '''

    # assign the individual file path
    individual_file = f'{xlib.get_config_dir()}/individuals.txt'

    # return the individual file path
    return individual_file

#-------------------------------------------------------------------------------

def get_individual_dict():
    '''
    Get a dictionary of individual enzymes with their rectriction sites.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of individuals
    individual_dict = {}

    # set the pattern of the record of individual file (record format: individual_id;replicated_individual_id;population_id;index1_seq;[index2_seq])
    pattern = r'^(.+);(.+);(.+);(.+);(.*)$'

    # check the individual file
    (OK, error_list) = check_individual_file(strict=True)

    # read every record of individual file
    if OK:
        with open(get_individual_file(), mode='r', encoding='iso-8859-1', newline='\n') as file_id:
            for record in file_id:

                # if the record is not a comment nor a line with blank characters
                if not record.lstrip().startswith('#') and record.strip() != '':

                    # extract the data and add individual data to the dictionary of individuals
                    mo = re.search(pattern, record)
                    individual_id = mo.group(1).strip()
                    replicated_individual_id = mo.group(2).strip()
                    population_id = mo.group(3).strip()
                    index1_seq = mo.group(4).strip()
                    index2_seq = mo.group(5).strip()

                    # add data to the dictionary of individuals
                    individual_dict[individual_id] = {'individual_id': individual_id, 'replicated_individual_id': replicated_individual_id, 'population_id': population_id, 'index1_seq': index1_seq, 'index2_seq': index2_seq}

    # return the control variable, the error list and the dictionary of individuals
    return (OK, error_list, individual_dict)

#-------------------------------------------------------------------------------

def create_rsitesearch_config_file(reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', enzyme1='EcoRI', enzyme2='MseI'):
    '''
    Create rsitesearch config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the rsitesearch config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_rsitesearch_config_file())):
            os.makedirs(os.path.dirname(get_rsitesearch_config_file()))
        with open(get_rsitesearch_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/reference_dataset_id\n')
            file_id.write( '# The reference_dataset_id and reference_file_name are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of rsitesearch (ddRADseqTools package) and their meaning in "https://github.com/GGFHF/ddRADseqTools".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the rsitesearch parameters\n')
            file_id.write( '[rsitesearch parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'enzyme1 = {enzyme1}', '# id of 1st restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format(f'enzyme2 = {enzyme2}', '# id of 2nd restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format( 'minfragsize = 101', '# lower loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format( 'maxfragsize = 300', '# upper loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format( 'fragstinterval = 25', '# interval length of fragment size'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_rsitesearch_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_rsitesearch_process(cluster_name, log, function=None):
    '''
    Run a rsitesearch process.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the file of restriction sites
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the file {get_restriction_site_file()} ...\n')
    (OK, error_list) = check_restriction_site_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        for error in error_list:
            log.write(f'{error}\n')
        OK = False

    # check the rsitesearch config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Checking the {xlib.get_rsitesearch_name()} config file ...\n')
        (OK, error_list) = check_rsitesearch_config_file(strict=True)
        if OK:
            log.write('The file is OK.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            OK = False

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

    # check the ddRADseqTools is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_ddradseqtools_name()} is not installed.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_design_dataset_name(), xlib.get_rsitesearch_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the file of restriction sites to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {get_restriction_site_file()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_restriction_site_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_restriction_site_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the rsitesearch process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_rsitesearch_process_script()} ...\n')
        (OK, error_list) = build_rsitesearch_process_script(cluster_name, current_run_dir, sftp_client)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the rsitesearch process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_rsitesearch_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_rsitesearch_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_rsitesearch_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the rsitesearch process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_rsitesearch_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_rsitesearch_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the rsitesearch process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_rsitesearch_process_starter()} ...\n')
        (OK, error_list) = build_rsitesearch_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the rsitesearch process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_rsitesearch_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_rsitesearch_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_rsitesearch_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the rsitesearch process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_rsitesearch_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_rsitesearch_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the rsitesearch process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_rsitesearch_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_rsitesearch_process_starter()), log)

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

def check_rsitesearch_config_file(strict):
    '''
    Check the rsitesearch config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the dictionary of restriction enzymes
    (is_OK_restriction_enzyme_dict, enzyme_dict_error_list, restriction_enzyme_dict) = get_restriction_enzyme_dict()
    if enzyme_dict_error_list != []:
        error_list = error_list + enzyme_dict_error_list + ['\n']

    # get the enzime identification list
    if is_OK_restriction_enzyme_dict:
        enzyme_id_list = list(restriction_enzyme_dict.keys())

    # get the option dictionary
    try:
        rsitesearch_option_dict = xlib.get_option_dict(get_rsitesearch_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in rsitesearch_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = rsitesearch_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = rsitesearch_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

        # check section "rsitesearch parameters"
        if 'rsitesearch parameters' not in sections_list:
            error_list.append('*** ERROR: the section "rsitesearch parameters" is not found.')
            OK = False
        else:

            # check section "rsitesearch parameters" - key "enzyme1"
            enzyme1 = rsitesearch_option_dict.get('rsitesearch parameters', {}).get('enzyme1', not_found)
            is_ok_enzyme1 = False
            if enzyme1 == not_found:
                error_list.append('*** ERROR: the key "enzyme1" is not found in the section "rsitesearch parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if enzyme1 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme1, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: the key "enzyme1" has to be a restriction enzyme id or a restriction site seq.')
                    OK = False
                else:
                    if enzyme1 in enzyme_id_list:
                        enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
                    else:
                        enzyme1_seq = enzyme1
                    is_ok_enzyme1 = True

            # check section "rsitesearch parameters" - key "enzyme2"
            enzyme2 = rsitesearch_option_dict.get('rsitesearch parameters', {}).get('enzyme2', not_found)
            is_ok_enzyme2 = False
            if enzyme2 == not_found:
                error_list.append('*** ERROR: the key "enzyme2" is not found in the section "rsitesearch parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if enzyme2 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme2, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: the key "enzyme2" has to be a restriction enzyme id or a restriction site seq.')
                    OK = False
                else:
                    if enzyme2 in enzyme_id_list:
                        enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
                    else:
                        enzyme2_seq = enzyme2
                    is_ok_enzyme2 = True

            # check that enzyme1 has to be different from enzyme2
            if is_ok_enzyme1 and is_ok_enzyme2 and enzyme1_seq == enzyme2_seq:
                error_list.append('*** ERROR: Both enzymes have the same sequence.')
                OK = False

            # check section "rsitesearch parameters" - key "minfragsize"
            minfragsize = rsitesearch_option_dict.get('rsitesearch parameters', {}).get('minfragsize', not_found)
            is_ok_minfragsize = False
            if minfragsize == not_found:
                error_list.append('*** ERROR: the key "minfragsize" is not found in the section "rsitesearch parameters".')
                OK = False
            elif not xlib.check_int(minfragsize, minimum=1):
                error_list.append('*** ERROR: the key "minfragsize" has to be an integer number greater than or equal to 1.')
                OK = False
            else:
                is_ok_minfragsize = True

            # check section "rsitesearch parameters" - key "maxfragsize"
            maxfragsize = rsitesearch_option_dict.get('rsitesearch parameters', {}).get('maxfragsize', not_found)
            is_ok_maxfragsize = False
            if maxfragsize == not_found:
                error_list.append('*** ERROR: the key "maxfragsize" is not found in the section "rsitesearch parameters".')
                OK = False
            elif not xlib.check_int(maxfragsize, minimum=1):
                error_list.append('*** ERROR: the key "maxfragsize" has to be an integer number greater than or equal to 1.')
                OK = False
            else:
                is_ok_maxfragsize = True

            # check if maxfragsize value is greater than or equal than minfragsize value
            if is_ok_minfragsize and is_ok_maxfragsize and int(maxfragsize) < int(minfragsize):
                error_list.append(f'*** ERROR: The value maxfragsize value ({maxfragsize}) is less than the minfragsize value ({minfragsize}).')
                OK = False

            # check section "rsitesearch parameters" - key "fragstinterval"
            fragstinterval = rsitesearch_option_dict.get('rsitesearch parameters', {}).get('fragstinterval', not_found)
            if fragstinterval == not_found:
                error_list.append('*** ERROR: the key "fragstinterval" is not found in the section "rsitesearch parameters".')
                OK = False
            elif not xlib.check_int(fragstinterval, minimum=1):
                error_list.append('*** ERROR: the key "fragstinterval" has to be an integer number greater than or equal to 1.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'The {xlib.get_rsitesearch_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the esult of the control variables and the error list
    return (OK and is_OK_restriction_enzyme_dict, error_list)

#-------------------------------------------------------------------------------

def build_rsitesearch_process_script(cluster_name, current_run_dir, sftp_client):
    '''
    Build the current rsitesearch process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the rsitesearch option dictionary
    rsitesearch_option_dict = xlib.get_option_dict(get_rsitesearch_config_file())

    # get the options
    reference_dataset_id = rsitesearch_option_dict['identification']['reference_dataset_id']
    reference_file = rsitesearch_option_dict['identification']['reference_file']
    enzyme1 = rsitesearch_option_dict['rsitesearch parameters']['enzyme1']
    enzyme2 = rsitesearch_option_dict['rsitesearch parameters']['enzyme2']
    minfragsize = rsitesearch_option_dict['rsitesearch parameters']['minfragsize']
    maxfragsize = rsitesearch_option_dict['rsitesearch parameters']['maxfragsize']
    fragstinterval = rsitesearch_option_dict['rsitesearch parameters']['fragstinterval']

    # set file paths
    genfile = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    rsfile = f'{current_run_dir}/{os.path.basename(get_restriction_site_file())}'

    fragsfile = f'{current_run_dir}/genome-fragments.txt'
    fragstfile = f'{current_run_dir}/genome-fragment-stats.txt'

    # write the rsitesearch process script
    try:
        if not os.path.exists(os.path.dirname(get_rsitesearch_process_script())):
            os.makedirs(os.path.dirname(get_rsitesearch_process_script()))
        with open(get_rsitesearch_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'PYTHON3_DIR={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'DDRADSEQTOOLS_DIR={xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()}/Package\n')
            script_file_id.write( 'export PATH=$PYTHON3_DIR:$DDRADSEQTOOLS_DIR:$PATH\n')
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
            script_file_id.write( 'function run_rsitesearch_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Analyzing enzymes of a ddRADseq experiment ..."\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '        rsitesearch.py \\\n')
            script_file_id.write(f'            --genfile={genfile} \\\n')
            script_file_id.write(f'            --fragsfile={fragsfile} \\\n')
            script_file_id.write(f'            --rsfile={rsfile} \\\n')
            script_file_id.write(f'            --enzyme1={enzyme1} \\\n')
            script_file_id.write(f'            --enzyme2={enzyme2} \\\n')
            script_file_id.write(f'            --minfragsize={minfragsize} \\\n')
            script_file_id.write(f'            --maxfragsize={maxfragsize} \\\n')
            script_file_id.write(f'            --fragstfile={fragstfile} \\\n')
            script_file_id.write(f'            --fragstinterval={fragstinterval} \\\n')
            script_file_id.write( '            --plot=YES \\\n')
            script_file_id.write( '            --verbose=NO \\\n')
            script_file_id.write( '            --trace=NO\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rsitesearch.py $RC; fi\n')
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
            process_name = f'{xlib.get_rsitesearch_name()} process'
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
            script_file_id.write( 'run_rsitesearch_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_rsitesearch_process_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_rsitesearch_process_starter(current_run_dir):
    '''
    Build the starter of the current rsitesearch process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the rsitesearch process starter
    try:
        if not os.path.exists(os.path.dirname(get_rsitesearch_process_starter())):
            os.makedirs(os.path.dirname(get_rsitesearch_process_starter()))
        with open(get_rsitesearch_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_rsitesearch_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_rsitesearch_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_rsitesearch_config_file():
    '''
    Get the rsitesearch config file path.
    '''

    # assign the rsitesearch config file path
    rsitesearch_config_file = f'{xlib.get_config_dir()}/{xlib.get_rsitesearch_code()}-config.txt'

    # return the rsitesearch config file path
    return rsitesearch_config_file

#-------------------------------------------------------------------------------

def get_rsitesearch_process_script():
    '''
    Get the rsitesearch process script path in the local computer.
    '''

    # assign the rsitesearch script path
    rsitesearch_process_script = f'{xlib.get_temp_dir()}/{xlib.get_rsitesearch_code()}-process.sh'

    # return the rsitesearch script path
    return rsitesearch_process_script

#-------------------------------------------------------------------------------

def get_rsitesearch_process_starter():
    '''
    Get the rsitesearch process starter path in the local computer.
    '''

    # assign the rsitesearch process starter path
    rsitesearch_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_rsitesearch_code()}-process-starter.sh'

    # return the rsitesearch starter path
    return rsitesearch_process_starter

#-------------------------------------------------------------------------------

def create_ddradseq_simulation_config_file(reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', enzyme1='EcoRI', enzyme2='MseI'):
    '''
    Create ddRADseq simulation config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the ddRADseq simulation config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_ddradseq_simulation_config_file())):
            os.makedirs(os.path.dirname(get_ddradseq_simulation_config_file()))
        with open(get_ddradseq_simulation_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/reference_dataset_id\n')
            file_id.write( '# The reference_dataset_id and reference_file_name are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of ddRADseqTools programs and their meaning in "https://github.com/GGFHF/ddRADseqTools".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format( 'rsfile = ./config/restrictionsites.txt', '# local path of the restriction sites file'))
            file_id.write( '{0:<50} {1}\n'.format( 'endsfile = ./config/ends.txt', '# local path oh the end sequences file'))
            file_id.write( '{0:<50} {1}\n'.format( 'individualsfile = ./config/individuals.txt', '# local path oh the end sequences file'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the ddRADseq simulation parameters\n')
            file_id.write( '[ddRADseq simulation parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'enzyme1 = {enzyme1}', '# id of 1st restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format(f'enzyme2 = {enzyme2}', '# id of 2nd restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format( 'minfragsize = 101', '# lower loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format( 'maxfragsize = 300', '# upper loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format( 'fragstinterval = 25', '# interval length of fragment size'))
            file_id.write( '{0:<50} {1}\n'.format( 'technique = IND1_IND2_DBR', f'# technique: {get_technique_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'readtype = PE', f'# read type: {get_read_type_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'index1len = 6', '# index sequence length in the adapter 1'))
            file_id.write( '{0:<50} {1}\n'.format( 'index2len = 6', '# index sequence length in the adapter 2 (it has to be 0 when technique is IND1)'))
            file_id.write( '{0:<50} {1}\n'.format( 'dbrlen = 4', '# DBR sequence length (it has to be 0 when technique is IND1 or IND1_IND2)'))
            file_id.write( '{0:<50} {1}\n'.format( 'wend = end91', '# code used in endsfile corresponding to the end where the adapter 1 is'))
            file_id.write( '{0:<50} {1}\n'.format( 'cend = end92', '# code used in endsfile corresponding to the end where the adapter 2 is'))
            file_id.write( '{0:<50} {1}\n'.format( 'locinum = 3000', '# loci number to sample'))
            file_id.write( '{0:<50} {1}\n'.format( 'readsnum = 300000', '# reads number'))
            file_id.write( '{0:<50} {1}\n'.format( 'minreadvar = 0.8', '# lower variation on reads number per locus (0.5 <= minreadvar <= 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'maxreadvar = 1.2', '# upper variation on reads number per locus (1.0 <= maxreadvar <= 1.5)'))
            file_id.write( '{0:<50} {1}\n'.format( 'insertlen = 100', '# read length, i. e. genome sequence length inserted in reads'))
            file_id.write( '{0:<50} {1}\n'.format( 'mutprob = 0.2', '# mutation probability (0.0 <= mutprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'locusmaxmut = 1', '# maximum mutations number by locus (1 <= locusmaxmut <= 5)'))
            file_id.write( '{0:<50} {1}\n'.format( 'indelprob = 0.1', '# insertion/deletion probability (0.0 <= indelprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'maxindelsize = 10', '# upper insertion/deletion size (1 <= maxindelsize < 30)'))
            file_id.write( '{0:<50} {1}\n'.format( 'dropout = 0.0', '# mutation probability in the enzyme recognition sites (0.0 <= dropout < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'pcrdupprob = 0.2', '# probability of loci bearing PCR duplicates (0.0 <= pcrdupprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'pcrdistribution = MULTINOMIAL', f'# distribution type to calculate the PCR duplicates number: {get_pcrdistribution_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'multiparam = 0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015', '# probability values to multinomial distribution with format prob1,prob2,...,probn (they have to sum 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format( 'poissonparam = 1.0', '# lambda value of the Poisson distribution'))
            file_id.write( '{0:<50} {1}\n'.format( 'gcfactor = 0.2', '# weight factor of GC ratio in a locus with PCR duplicates (0.0 <= gcfactor < 1.0)'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ddradseq_simulation_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_ddradseq_simulation_process(cluster_name, log, function=None):
    '''
    Run a ddRADseq simulation process.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the ddRADseq simulation config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Checking the {xlib.get_ddradseq_simulation_name()} config file ...\n')
        (OK, error_list) = check_ddradseq_simulation_config_file(strict=True)
        if OK:
            log.write('The file is OK.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            OK = False

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

    # check the ddRADseq simulation is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_ddradseqtools_name()} is not installed.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_design_dataset_name(), xlib.get_ddradseq_simulation_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the file of restriction sites to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {get_restriction_site_file()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_restriction_site_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_restriction_site_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the file of ends to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {current_run_dir} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_end_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_end_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the file of individuals to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the file {get_individual_file()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_individual_file())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_individual_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the ddRADseq simulation process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_ddradseq_simulation_process_script()} ...\n')
        (OK, error_list) = build_ddradseq_simulation_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the ddRADseq simulation process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_ddradseq_simulation_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseq_simulation_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseq simulation process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseq simulation process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_ddradseq_simulation_process_starter()} ...\n')
        (OK, error_list) = build_ddradseq_simulation_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the ddRADseq simulation process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_ddradseq_simulation_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseq_simulation_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseq simulation process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the ddRADseq simulation process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter()), log)

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

def check_ddradseq_simulation_config_file(strict):
    '''
    Check the ddRADseq simulation config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the dictionary of restriction enzymes
    (is_OK_restriction_enzyme_dict, enzyme_dict_error_list, restriction_enzyme_dict) = get_restriction_enzyme_dict()
    if enzyme_dict_error_list != []:
        error_list = error_list + enzyme_dict_error_list + ['\n']

    # get the enzime identification list
    if is_OK_restriction_enzyme_dict:
        enzyme_id_list = list(restriction_enzyme_dict.keys())

    # get the dictionary of end
    (is_OK_end_dict, end_dict_error_list, end_dict) = get_end_dict()
    if end_dict_error_list != []:
        error_list = error_list + end_dict_error_list + ['\n']

    # get the end identification list
    if is_OK_end_dict:
        end_id_list = list(end_dict.keys())

    # get the dictionary of individuals
    (is_OK_individual_dict, individual_dict_error_list, individual_dict) = get_individual_dict()
    if individual_dict_error_list != []:
        error_list = error_list + individual_dict_error_list + ['\n']

    # get the individual identification list
    if is_OK_individual_dict:
        individual_id_list = list(individual_dict.keys())

    # get the option dictionary
    try:
        ddradseq_simulation_option_dict = xlib.get_option_dict(get_ddradseq_simulation_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in ddradseq_simulation_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = ddradseq_simulation_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = ddradseq_simulation_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

        # check section "ddRADseq simulation parameters"
        if 'ddRADseq simulation parameters' not in sections_list:
            error_list.append('*** ERROR: the section "ddRADseq simulation parameters" is not found.')
            OK = False
        else:

            # check section "ddRADseq simulation parameters" - key "enzyme1"
            enzyme1 = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('enzyme1', not_found)
            is_ok_enzyme1 = False
            if enzyme1 == not_found:
                error_list.append('*** ERROR: the key "enzyme1" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if enzyme1 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme1, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: the key "enzyme1" has to be a restriction enzyme id or a restriction site seq.')
                    OK = False
                else:
                    if enzyme1 in enzyme_id_list:
                        enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
                    else:
                        enzyme1_seq = enzyme1
                    is_ok_enzyme1 = True

            # check section "ddRADseq simulation parameters" - key "enzyme2"
            enzyme2 = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('enzyme2', not_found)
            is_ok_enzyme2 = False
            if enzyme2 == not_found:
                error_list.append('*** ERROR: the key "enzyme2" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if enzyme2 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme2, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: the key "enzyme2" has to be a restriction enzyme id or a restriction site seq.')
                    OK = False
                else:
                    if enzyme2 in enzyme_id_list:
                        enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
                    else:
                        enzyme2_seq = enzyme2
                    is_ok_enzyme2 = True

            # check that enzyme1 has to be different from enzyme2
            if is_ok_enzyme1 and is_ok_enzyme2 and enzyme1_seq == enzyme2_seq:
                error_list.append('*** ERROR: Both enzymes have the same sequence.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "minfragsize"
            minfragsize = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('minfragsize', not_found)
            is_ok_minfragsize = False
            if minfragsize == not_found:
                error_list.append('*** ERROR: the key "minfragsize" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(minfragsize, minimum=1):
                error_list.append('*** ERROR: the key "minfragsize" has to be an integer number greater than or equal to 1.')
                OK = False
            else:
                is_ok_minfragsize = True

            # check section "ddRADseq simulation parameters" - key "maxfragsize"
            maxfragsize = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('maxfragsize', not_found)
            is_ok_maxfragsize = False
            if maxfragsize == not_found:
                error_list.append('*** ERROR: the key "maxfragsize" is not found in the section "ddRADseq simulation parameters".')
                OK = False
                maxfragsize = 99999999999999
            elif not xlib.check_int(maxfragsize, minimum=2):
                error_list.append('*** ERROR: the key "maxfragsize" has to be an integer number greater than or equal to 2.')
                OK = False
            else:
                is_ok_maxfragsize = True

            # check if maxfragsize value is greater than or equal than minfragsize value
            if is_ok_minfragsize and is_ok_maxfragsize and int(maxfragsize) < int(minfragsize):
                error_list.append(f'*** ERROR: The value maxfragsize value ({maxfragsize}) is less than the minfragsize value ({minfragsize}).')
                OK = False

            # check section "ddRADseq simulation parameters" - key "fragstinterval"
            fragstinterval = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('fragstinterval', not_found)
            if fragstinterval == not_found:
                error_list.append('*** ERROR: the key "fragstinterval" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(fragstinterval, minimum=1):
                error_list.append('*** ERROR: the key "fragstinterval" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "technique"
            technique = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('technique', not_found)
            is_ok_technique = False
            if technique == not_found:
                error_list.append('*** ERROR: the key "technique" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(technique, get_technique_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "technique" has to be {get_technique_code_list_text()}.')
                OK = False
            else:
                is_ok_technique = True

            # check section "ddRADseq simulation parameters" - key "format"
            format = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "readtype"
            readtype = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('readtype', not_found)
            if readtype == not_found:
                error_list.append('*** ERROR: the key "readtype" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(readtype, get_read_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "readtype" has to be {get_read_type_code_list_text()}.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "index1len"
            index1len = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('index1len', not_found)
            is_ok_index1len = False
            if index1len == not_found:
                error_list.append('*** ERROR: the key "index1len" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(index1len, minimum=1):
                error_list.append('*** ERROR: the key "index1len" has to be an integer number greater than or equal to 1.')
                OK = False
            else:
                is_ok_index1len = True

            # check that index1len is equal to the lenght of the sequence of index1 in the individual file
            if is_OK_individual_dict and is_ok_index1len and int(index1len) != len(individual_dict[individual_id_list[0]]['index1_seq']):
                error_list.append('*** ERROR: the index1 sequence in individual file does not have a lenght equal to index1len.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "index2len"
            index2len = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('index2len', not_found)
            is_ok_index2len = False
            if index2len == not_found:
                error_list.append('*** ERROR: the key "index2len" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(index2len, minimum=1):
                error_list.append('*** ERROR: the key "index2len" has to be an integer number greater than or equal to 1.')
                OK = False
            elif technique in ['IND1', 'IND1_DBR'] and int(index2len) != 0:
                error_list.append('*** ERROR: index2len has to be equal to 0 if techinque is IND1 or IND1_DBR.')
                OK = False
            elif technique in ['IND1_IND2', 'IND1_IND2_DBR'] and int(index2len) == 0:
                error_list.append('*** ERROR: index2len has to be greater than or equal to 1 if techinque is IND1_IND2 or IND1_IND2_DBR.')
                OK = False
            else:
                is_ok_index2len = True

            # check that index2len is equal to the lenght of the sequence of index2 in the individual file
            if is_OK_individual_dict and is_ok_index1len and int(index2len) != len(individual_dict[individual_id_list[0]]['index2_seq']):
                error_list.append('*** ERROR: the index2 sequence in individual file does not have a lenght equal to index2len.')
                OK = False

            # check index2len depending on technique
            if is_ok_technique and is_ok_index2len and technique in ['IND1', 'IND1_DBR'] and int(index2len) != 0:
                error_list.append('*** ERROR: index2len has to be equal to 0 if techinque is IND1 or IND1_DBR.')
                OK = False
            elif is_ok_technique and is_ok_index2len and technique in ['IND1_IND2', 'IND1_IND2_DBR'] and int(index2len) == 0:
                error_list.append('*** ERROR: index2len has to be greater than or equal to 1 if techinque is IND1_IND2 or IND1_IND2_DBR.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "dbrlen"
            dbrlen = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('dbrlen', not_found)
            is_ok_dbrlen = False
            if dbrlen == not_found:
                error_list.append('*** ERROR: the key "dbrlen" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(dbrlen, minimum=0):
                error_list.append('*** ERROR: the key "dbrlen" has to be an integer number greater than or equal to 0.')
                OK = False
            else:
                is_ok_dbrlen = True
                    
            # check dbrlen depending on technique
            if is_ok_technique and is_ok_dbrlen and technique in ['IND1', 'IND1_IND2'] and int(dbrlen) != 0:
                error_list.append('*** ERROR: dbrlen has to be equal to 0 if techinque is IND1 or IND1_IND2.')
                OK = False
            elif is_ok_technique and is_ok_dbrlen and technique in ['IND1_DBR', 'IND1_IND2_DBR'] and int(dbrlen) == 0:
                error_list.append('*** ERROR: dbrlen has to be greater than or equal to 1 if techinque is IND1_DBR or IND1_IND2_DBR.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "wend"
            wend = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('wend', not_found)
            is_ok_wend = False
            if wend == not_found:
                error_list.append('*** ERROR: the key "wend" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif is_OK_end_dict and wend not in end_id_list:
                error_list.append('*** ERROR: the key "wend" has to be an end id.')
                OK = False
            else:
                is_ok_wend = True

            # check section "ddRADseq simulation parameters" - key "cend"
            cend = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('cend', not_found)
            is_ok_cend = False
            if cend == not_found:
                error_list.append('*** ERROR: the key "cend" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif is_OK_end_dict and cend not in end_id_list:
                error_list.append('*** ERROR: the key "cend" has to be an end id.')
                OK = False
            else:
                is_ok_cend = True

            # check that wend has to be different from cend
            if is_ok_wend and is_ok_cend:
                if wend == cend:
                    error_list.append('*** ERROR: Both ends have the same id.')
                    OK = False

            # check section "ddRADseq simulation parameters" - key "locinum"
            locinum = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('locinum', not_found)
            if locinum == not_found:
                error_list.append('*** ERROR: the key "locinum" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(locinum, minimum=0):
                error_list.append('*** ERROR: the key "locinum" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "readsnum"
            readsnum = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('readsnum', not_found)
            if readsnum == not_found:
                error_list.append('*** ERROR: the key "readsnum" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(readsnum, minimum=0):
                error_list.append('*** ERROR: the key "readsnum" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "minreadvar"
            minreadvar = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('minreadvar', not_found)
            is_ok_minreadvar = False
            if minreadvar == not_found:
                error_list.append('*** ERROR: the key "minreadvar" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(minreadvar, minimum=0.5, maximum=1.):
                error_list.append('*** ERROR: the key "minreadvar" has to be a float number between 0.5 and 1.0.')
                OK = False
            else:
                is_ok_minreadvar = True

            # check section "ddRADseq simulation parameters" - key "maxreadvar"
            maxreadvar = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('maxreadvar', not_found)
            is_ok_maxreadvar = False
            if maxreadvar == not_found:
                error_list.append('*** ERROR: the key "maxreadvar" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(maxreadvar, minimum=0.5, maximum=1.5):
                error_list.append('*** ERROR: the key "maxreadvar" has to be a float number between 1.0 and 1.5.')
                OK = False
            else:
                is_ok_maxreadvar = True

            # check that maxreadvar has to be greater than or equal to minreadvar
            if is_ok_minreadvar and is_ok_maxreadvar and maxreadvar < minreadvar:
                error_list.append('*** ERROR: maxreadvar has to be greater than or equal to minreadvar.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "insertlen"
            insertlen = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('insertlen', not_found)
            if insertlen == not_found:
                error_list.append('*** ERROR: the key "insertlen" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(insertlen, minimum=1):
                error_list.append('*** ERROR: the key "insertlen" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "mutprob"
            mutprob = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('mutprob', not_found)
            if mutprob == not_found:
                error_list.append('*** ERROR: the key "mutprob" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(mutprob, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "mutprob" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "locusmaxmut"
            locusmaxmut = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('locusmaxmut', not_found)
            if locusmaxmut == not_found:
                error_list.append('*** ERROR: the key "locusmaxmut" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(locusmaxmut, minimum=1, maximum=5):
                error_list.append('*** ERROR: the key "locusmaxmut" has to be an integer number between 1 and 5.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "indelprob"
            indelprob = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('indelprob', not_found)
            if indelprob == not_found:
                error_list.append('*** ERROR: the key "indelprob" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(indelprob, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "indelprob" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "maxindelsize"
            maxindelsize = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('maxindelsize', not_found)
            if maxindelsize == not_found:
                error_list.append('*** ERROR: the key "maxindelsize" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_int(maxindelsize, minimum=1, maximum=30):
                error_list.append('*** ERROR: the key "maxindelsize" has to be an integer number between 1 and 30.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "dropout"
            dropout = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('dropout', not_found)
            if dropout == not_found:
                error_list.append('*** ERROR: the key "dropout" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(dropout, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "dropout" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "pcrdupprob"
            pcrdupprob = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('pcrdupprob', not_found)
            if pcrdupprob == not_found:
                error_list.append('*** ERROR: the key "pcrdupprob" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(pcrdupprob, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "pcrdupprob" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "pcrdistribution"
            pcrdistribution = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('pcrdistribution', not_found)
            if pcrdistribution == not_found:
                error_list.append('*** ERROR: the key "pcrdistribution" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(pcrdistribution, get_pcrdistribution_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "pcrdistribution" has to be {get_pcrdistribution_code_list_text()}.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "multiparam"
            multiparam = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('multiparam', not_found)
            if multiparam == not_found:
                error_list.append('*** ERROR: the key "multiparam" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            else:
                multiparam_list = xlib.split_literal_to_float_list(multiparam)
                if multiparam_list == []:
                    error_list.append('*** ERROR: the key "multiparam" has to be a float number or a float values list.')
                    OK = False

            # check section "ddRADseq simulation parameters" - key "poissonparam"
            poissonparam = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('poissonparam', not_found)
            if poissonparam == not_found:
                error_list.append('*** ERROR: the key "poissonparam" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(poissonparam):
                error_list.append('*** ERROR: the key "poissonparam" has to be a float number.')
                OK = False

            # check section "ddRADseq simulation parameters" - key "gcfactor"
            gcfactor = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('gcfactor', not_found)
            if gcfactor == not_found:
                error_list.append('*** ERROR: the key "gcfactor" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_float(gcfactor, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "gcfactor" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_ddradseq_simulation_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the result of the control variables and the error list
    return (OK and is_OK_restriction_enzyme_dict and is_OK_end_dict and is_OK_individual_dict, error_list)

#-------------------------------------------------------------------------------

def build_ddradseq_simulation_process_script(cluster_name, current_run_dir):
    '''
    Build the current ddRADseq simulation process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the ddRADseq simulation option dictionary
    ddradseq_simulation_option_dict = xlib.get_option_dict(get_ddradseq_simulation_config_file())

    # get the options
    reference_dataset_id = ddradseq_simulation_option_dict['identification']['reference_dataset_id']
    reference_file = ddradseq_simulation_option_dict['identification']['reference_file']
    enzyme1 = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['enzyme1']
    enzyme2 = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['enzyme2']
    minfragsize = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['minfragsize']
    maxfragsize = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['maxfragsize']
    fragstinterval = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['fragstinterval']
    technique = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['technique'].upper()
    format = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['format'].upper()
    readtype = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['readtype'].upper()
    index1len = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['index1len']
    index2len = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['index2len']
    dbrlen = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['dbrlen']
    wend = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['wend']
    cend = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['cend']
    locinum = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['locinum']
    readsnum = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['readsnum']
    minreadvar = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['minreadvar']
    maxreadvar = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['maxreadvar']
    insertlen = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['insertlen']
    mutprob = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['mutprob']
    locusmaxmut = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['locusmaxmut']
    indelprob = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['indelprob']
    maxindelsize = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['maxindelsize']
    dropout = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['dropout']
    pcrdupprob = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['pcrdupprob']
    pcrdistribution = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['pcrdistribution'].upper()
    multiparam = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['multiparam']
    poissonparam = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['poissonparam']
    gcfactor = ddradseq_simulation_option_dict['ddRADseq simulation parameters']['gcfactor']

    # set file paths
    genfile = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    rsfile = f'{current_run_dir}/{os.path.basename(get_restriction_site_file())}'
    endsfile = f'{current_run_dir}/{os.path.basename(get_end_file())}'
    individualsfile = f'{current_run_dir}/{os.path.basename(get_individual_file())}'
    fragsfile = f'{current_run_dir}/genome-fragments.txt'
    fragstfile = f'{current_run_dir}/genome-fragment-stats.txt'
    readsfile = f'{current_run_dir}/reads'
    readsfile1 = f'{current_run_dir}/reads-1'
    readsfile2 = f'{current_run_dir}/reads-2'
    clearfile = f'{current_run_dir}/reads-cleared'
    clearfile1 = f'{current_run_dir}/reads-cleared-1'
    clearfile2 = f'{current_run_dir}/reads-cleared-2'
    dupstfile = f'{current_run_dir}/pcrduplicates-stats'
    read_file_list = f'{current_run_dir}/reads-files.txt'
    demultiplexed_file = f'{current_run_dir}/demultiplexed-ind'

    # write the ddRADseq simulation process script
    try:
        if not os.path.exists(os.path.dirname(get_ddradseq_simulation_process_script())):
            os.makedirs(os.path.dirname(get_ddradseq_simulation_process_script()))
        with open(get_ddradseq_simulation_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'PYTHON3_DIR={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'DDRADSEQTOOLS_DIR={xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()}/Package\n')
            script_file_id.write( 'export PATH=$PYTHON3_DIR:$DDRADSEQTOOLS_DIR:$PATH\n')
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
            script_file_id.write( 'function run_rsitesearch_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_rsitesearch_process.ok\n')
            script_file_id.write( '    echo "Generating fragments ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            rsitesearch.py \\\n')
            script_file_id.write(f'                --genfile={genfile} \\\n')
            script_file_id.write(f'                --fragsfile={fragsfile} \\\n')
            script_file_id.write(f'                --rsfile={rsfile} \\\n')
            script_file_id.write(f'                --enzyme1={enzyme1} \\\n')
            script_file_id.write(f'                --enzyme2={enzyme2} \\\n')
            script_file_id.write(f'                --minfragsize={minfragsize} \\\n')
            script_file_id.write(f'                --maxfragsize={maxfragsize} \\\n')
            script_file_id.write(f'                --fragstfile={fragstfile} \\\n')
            script_file_id.write(f'                --fragstinterval={fragstinterval} \\\n')
            script_file_id.write( '                --plot=YES \\\n')
            script_file_id.write( '                --verbose=NO \\\n')
            script_file_id.write( '                --trace=NO\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rsitesearch.py $RC; fi\n')
            script_file_id.write( '        echo "Fragments are generated."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_simddradseq_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_simddradseq_process.ok\n')
            script_file_id.write( '    echo "Generating reads ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            simddradseq.py \\\n')
            script_file_id.write(f'                --fragsfile={fragsfile} \\\n')
            script_file_id.write(f'                --technique={technique} \\\n')
            script_file_id.write(f'                --format={format} \\\n')
            script_file_id.write(f'                --readsfile={readsfile} \\\n')
            script_file_id.write(f'                --readtype={readtype} \\\n')
            script_file_id.write(f'                --rsfile={rsfile} \\\n')
            script_file_id.write(f'                --enzyme1={enzyme1} \\\n')
            script_file_id.write(f'                --enzyme2={enzyme2} \\\n')
            script_file_id.write(f'                --endsfile={endsfile} \\\n')
            script_file_id.write(f'                --index1len={index1len} \\\n')
            script_file_id.write(f'                --index2len={index2len} \\\n')
            script_file_id.write(f'                --dbrlen={dbrlen} \\\n')
            script_file_id.write(f'                --wend={wend} \\\n')
            script_file_id.write(f'                --cend={cend} \\\n')
            script_file_id.write(f'                --individualsfile={individualsfile} \\\n')
            script_file_id.write(f'                --locinum={locinum} \\\n')
            script_file_id.write(f'                --readsnum={readsnum} \\\n')
            script_file_id.write(f'                --minreadvar={minreadvar} \\\n')
            script_file_id.write(f'                --maxreadvar={maxreadvar} \\\n')
            script_file_id.write(f'                --insertlen={insertlen} \\\n')
            script_file_id.write(f'                --mutprob={mutprob} \\\n')
            script_file_id.write(f'                --locusmaxmut={locusmaxmut} \\\n')
            script_file_id.write(f'                --indelprob={indelprob} \\\n')
            script_file_id.write(f'                --maxindelsize={maxindelsize} \\\n')
            script_file_id.write(f'                --dropout={dropout} \\\n')
            script_file_id.write(f'                --pcrdupprob={pcrdupprob} \\\n')
            script_file_id.write(f'                --pcrdistribution={pcrdistribution} \\\n')
            script_file_id.write(f'                --multiparam={multiparam} \\\n')
            script_file_id.write(f'                --poissonparam={poissonparam} \\\n')
            script_file_id.write(f'                --gcfactor={gcfactor} \\\n')
            script_file_id.write( '                --verbose=NO \\\n')
            script_file_id.write( '                --trace=NO\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error simddradseq.py $RC; fi\n')
            script_file_id.write( '        echo "Reads are generated."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            if int(dbrlen) > 0:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function run_pcrdupremoval_process\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_pcrdupremoval_process.ok\n')
                script_file_id.write( '    echo "Removing PCR duplicates ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            pcrdupremoval.py \\\n')
                script_file_id.write(f'                --format={format} \\\n')
                script_file_id.write(f'                --readtype={readtype} \\\n')
                if format == 'FASTA':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={readsfile}.fasta \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={readsfile1}.fasta \\\n')
                        script_file_id.write(f'                --readsfile2={readsfile2}.fasta \\\n')
                elif format == 'FASTQ':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={readsfile}.fastq \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={readsfile1}.fastq \\\n')
                        script_file_id.write(f'                --readsfile2={readsfile2}.fastq \\\n')
                script_file_id.write(f'                --clearfile={clearfile} \\\n')
                script_file_id.write(f'                --dupstfile={dupstfile} \\\n')
                script_file_id.write( '                --plot=YES \\\n')
                script_file_id.write( '                --verbose=NO \\\n')
                script_file_id.write( '                --trace=NO\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error pcrdupremoval.py $RC; fi\n')
                script_file_id.write( '        echo "PCR duplicates are removed."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_indsdemultiplexing_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_indsdemultiplexing_process.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Demultiplexing individuals ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            indsdemultiplexing.py \\\n')
            script_file_id.write(f'                --technique={technique} \\\n')
            script_file_id.write(f'                --format={format} \\\n')
            script_file_id.write(f'                --readtype={readtype} \\\n')
            script_file_id.write(f'                --endsfile={endsfile} \\\n')
            script_file_id.write(f'                --index1len={index1len} \\\n')
            script_file_id.write(f'                --index2len={index2len} \\\n')
            script_file_id.write(f'                --dbrlen={dbrlen} \\\n')
            script_file_id.write(f'                --wend={wend} \\\n')
            script_file_id.write(f'                --cend={cend} \\\n')
            script_file_id.write(f'                --individualsfile={individualsfile} \\\n')
            if int(dbrlen) > 0:
                if format == 'FASTA':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={clearfile}.fasta \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={clearfile1}.fasta \\\n')
                        script_file_id.write(f'                --readsfile2={clearfile2}.fasta \\\n')
                elif format == 'FASTQ':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={clearfile}.fastq \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={clearfile1}.fastq \\\n')
                        script_file_id.write(f'                --readsfile2={clearfile2}.fastq \\\n')
            else:
                if format == 'FASTA':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={readsfile}.fasta \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={readsfile1}.fasta \\\n')
                        script_file_id.write(f'                --readsfile2={readsfile2}.fasta \\\n')
                elif format == 'FASTQ':
                    if readtype == 'SE':
                        script_file_id.write(f'                --readsfile1={readsfile}.fastq \\\n')
                        script_file_id.write( '                --readsfile2=NONE \\\n')
                    elif readtype == 'PE':
                        script_file_id.write(f'                --readsfile1={readsfile1}.fastq \\\n')
                        script_file_id.write(f'                --readsfile2={readsfile2}.fastq \\\n')
            script_file_id.write( '                --verbose=NO \\\n')
            script_file_id.write( '                --trace=NO\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error indsdemultiplexing.py $RC; fi\n')
            script_file_id.write( '        echo "Individuals are demultiplexed."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_readstrim_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_readstrim_process.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Trimming reads ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            if format == 'FASTA':
                if readtype == 'SE':
                    script_file_id.write(f'        ls {demultiplexed_file}*.fasta > {read_file_list}\n')
                elif readtype == 'PE':
                    script_file_id.write(f'        ls {demultiplexed_file}*-1.fasta > {read_file_list}\n')
            elif format == 'FASTQ':
                if readtype == 'SE':
                    script_file_id.write(f'        ls {demultiplexed_file}*.fastq > {read_file_list}\n')
                elif readtype == 'PE':
                    script_file_id.write(f'        ls {demultiplexed_file}*-1.fastq > {read_file_list}\n')
            script_file_id.write( '        while read FILE_1; do\n')
            script_file_id.write( '            if [[ $FILE_1 =~ .*errors.* ]]; then continue; fi\n')
            script_file_id.write( '            echo "$SEP"\n')
            script_file_id.write( '            echo "... file $FILE_1 ..."\n')
            if readtype == 'SE':
                script_file_id.write( '            FILE_2=NONE\n')
            elif readtype == 'PE':
                if format == 'FASTA':
                    script_file_id.write( '            FILE_2=`echo $FILE_1 | sed "s/-1.fasta/-2.fasta/g"`\n')
                elif format == 'FASTQ':
                    script_file_id.write( '            FILE_2=`echo $FILE_1 | sed "s/-1.fastq/-2.fastq/g"`\n')
            if format == 'FASTA':
                script_file_id.write( '            FILE_TRIMMED=`echo $FILE_1 | sed "s/-1.fasta/-trimmed/g"`\n')
            elif format == 'FASTQ':
                script_file_id.write( '            FILE_TRIMMED=`echo $FILE_1 | sed "s/-1.fastq/-trimmed/g"`\n')
            script_file_id.write( '            /usr/bin/time \\\n')
            script_file_id.write(f'                --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write(f'                readstrim.py \\\n')
            script_file_id.write(f'                    --technique={technique} \\\n')
            script_file_id.write(f'                    --format={format} \\\n')
            script_file_id.write(f'                    --readtype={readtype} \\\n')
            script_file_id.write(f'                    --endsfile={endsfile} \\\n')
            script_file_id.write(f'                    --index1len={index1len} \\\n')
            script_file_id.write(f'                    --index2len={index2len} \\\n')
            script_file_id.write(f'                    --dbrlen={dbrlen} \\\n')
            script_file_id.write(f'                    --wend={wend} \\\n')
            script_file_id.write(f'                    --cend={cend} \\\n')
            script_file_id.write( '                    --readsfile1=$FILE_1 \\\n')
            if readtype == 'SE':
                script_file_id.write( '                    --readsfile2=NONE \\\n')
            elif readtype == 'PE':
                script_file_id.write( '                    --readsfile2=$FILE_2 \\\n')
            script_file_id.write( '                    --trimfile=$FILE_TRIMMED \\\n')
            script_file_id.write( '                    --verbose=NO \\\n')
            script_file_id.write( '                    --trace=NO\n')
            script_file_id.write( '                RC=$?\n')
            script_file_id.write( '                if [ $RC -ne 0 ]; then manage_error readstrim.py $RC; fi\n')
            script_file_id.write(f'        done < {read_file_list}\n')
            script_file_id.write( '        echo "Reads are trimmed."\n')
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
            process_name = f'{xlib.get_ddradseq_simulation_name()} process'
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
            script_file_id.write( 'run_rsitesearch_process\n')
            script_file_id.write( 'run_simddradseq_process\n')
            if int(dbrlen) > 0:
                script_file_id.write( 'run_pcrdupremoval_process\n')
            script_file_id.write( 'run_indsdemultiplexing_process\n')
            script_file_id.write( 'run_readstrim_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ddradseq_simulation_process_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ddradseq_simulation_process_starter(current_run_dir):
    '''
    Build the starter of the current ddRADseq simulation process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the ddRADseq simulation process starter
    try:
        if not os.path.exists(os.path.dirname(get_ddradseq_simulation_process_starter())):
            os.makedirs(os.path.dirname(get_ddradseq_simulation_process_starter()))
        with open(get_ddradseq_simulation_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ddradseq_simulation_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_ddradseq_simulation_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a ddRADseq simulation process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ddradseq_simulation_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter()), log)

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

def get_ddradseq_simulation_config_file():
    '''
    Get the ddRADseq simulation config file path.
    '''

    # assign the ddRADseq simulation config file path
    ddradseq_simulation_config_file = f'{xlib.get_config_dir()}/{xlib.get_ddradseq_simulation_code()}-config.txt'

    # return the ddRADseq simulation config file path
    return ddradseq_simulation_config_file

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_process_script():
    '''
    Get the ddRADseq simulation process script path in the local computer.
    '''

    # assign the ddRADseq simulation script path
    ddradseq_simulation_process_script = f'{xlib.get_temp_dir()}/{xlib.get_ddradseq_simulation_code()}-process.sh'

    # return the ddRADseq simulation script path
    return ddradseq_simulation_process_script

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_process_starter():
    '''
    Get the ddRADseq simulation process starter path in the local computer.
    '''

    # assign the ddRADseq simulation process starter path
    ddradseq_simulation_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_ddradseq_simulation_code()}-process-starter.sh'

    # return the ddRADseq simulation starter path
    return ddradseq_simulation_process_starter

#-------------------------------------------------------------------------------

def create_variant_calling_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', assembly_dataset_id='NONE', assembly_type='NONE', alignment_dataset_id='hisat2-170101-235959'):
    '''
    Create Variant calling config file with the default options. It is necessary
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
    elif assembly_dataset_id.startswith(xlib.get_ggtrinity_code()):
        assembly_software = xlib.get_ggtrinity_code()
    elif assembly_dataset_id.startswith(xlib.get_cd_hit_est_code()):
        assembly_software = xlib.get_cd_hit_est_code()
    elif assembly_dataset_id.startswith(xlib.get_transcript_filter_code()):
        assembly_software = xlib.get_transcript_filter_code()
    elif assembly_dataset_id.startswith(xlib.get_soapdenovo2_code()):
        assembly_software = xlib.get_soapdenovo2_code()
    elif assembly_dataset_id.startswith(xlib.get_starcode_code()):
        assembly_software = xlib.get_starcode_code()
    elif assembly_dataset_id.upper() == 'NONE':
        assembly_software = 'NONE'

    # set the alignment
    if alignment_dataset_id.startswith(xlib.get_bowtie2_code()):
        alignment_software = xlib.get_bowtie2_code()
    elif alignment_dataset_id.startswith(xlib.get_gsnap_code()):
        alignment_software = xlib.get_gsnap_code()
    elif alignment_dataset_id.startswith(xlib.get_hisat2_code()):
        alignment_software = xlib.get_hisat2_code()
    elif alignment_dataset_id.startswith(xlib.get_star_code()):
        alignment_software = xlib.get_star_code()
    elif alignment_dataset_id.startswith(xlib.get_tophat_code()):
        alignment_software = xlib.get_tophat_code()

    # create the Variant calling config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_variant_calling_config_file())):
            os.makedirs(os.path.dirname(get_variant_calling_config_file()))
        with open(get_variant_calling_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/reference_dataset_id\n')
            file_id.write(f'# The alignment file has to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/alignment_dataset_id\n')
            file_id.write( '# The experiment_id, reference_dataset_id reference_file_name and alignment_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of Variant calling (ddRADseqTools package) and their meaning in "https://github.com/GGFHF/ddRADseqTools"\n')
            file_id.write( '# and the ones of BCFtools and their meaning in "https://samtools.github.io/bcftools/bcftools.html".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "bcftools call parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of bcftools call and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --keep-alts; --pval-threshold=0.5\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_extended_assembly_software_code_list_text()}; or NONE if a reference is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification or NONE if a reference is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_software = {alignment_software}', f'# alignment software: {get_alignment_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_dataset_id = {alignment_dataset_id}', '# alignment dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the Variant calling parameters\n')
            file_id.write( '[Variant calling parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'vcf-merger = NO', f'# merger of the VCF files: {get_vcf_merger_code_list_text()}'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the bcftools call parameters\n')
            file_id.write( '[bcftools call parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'variants-only = YES', f'# output variant sites only: {get_variants_only_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'consensus-caller = YES', f'# the old samtools calling model (conflicts with multiallelic-caller): {get_consensus_caller_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'multiallelic-caller = NO', f'#  the alternative model for multiallelic and rare-variant calling  (conflicts with consensus-caller): {get_multiallelic_caller_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_variant_calling_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_variant_calling_process(cluster_name, log, function=None):
    '''
    Run a Variant calling process.
    '''

    # initialize the control variable
    OK = True

    # get the Variant calling option dictionary
    variant_calling_option_dict = xlib.get_option_dict(get_variant_calling_config_file())

    # get the experiment identification
    experiment_id = variant_calling_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the Variant calling config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Checking the {xlib.get_variant_calling_name()} config file ...\n')
        (OK, error_list) = check_variant_calling_config_file(strict=True)
        if OK:
            log.write('The file is OK.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            OK = False

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

    # check the ddRADseqTools is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_ddradseqtools_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_ddradseqtools_name()} is not installed.\n')
            OK = False

    # check SAMtools is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_samtools_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_samtools_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_samtools_name()} installation could not be performed.\n')

    # check BCFtools is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_bcftools_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_bcftools_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_bcftools_name()} installation could not be performed.\n')

    # check BEDtools is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_bedtools_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_bedtools_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_bedtools_name()} installation could not be performed.\n')

    # check VCFtools is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_vcftools_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_vcftools_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_vcftools_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_variant_calling_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Variant calling process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_variant_calling_process_script()} ...\n')
        (OK, error_list) = build_variant_calling_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Variant calling process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_variant_calling_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_variant_calling_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_variant_calling_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Variant calling process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_variant_calling_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_variant_calling_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Variant calling process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_variant_calling_process_starter()} ...\n')
        (OK, error_list) = build_variant_calling_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the Variant calling process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_variant_calling_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_variant_calling_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_variant_calling_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Variant calling process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_variant_calling_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_variant_calling_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Variant calling process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_variant_calling_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_variant_calling_process_starter()), log)

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

def check_variant_calling_config_file(strict):
    '''
    Check the Variant calling config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        variant_calling_option_dict = xlib.get_option_dict(get_variant_calling_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in variant_calling_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = variant_calling_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = variant_calling_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            is_ok_reference_dataset_id = False
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False
            else:
                is_ok_reference_dataset_id = True

            # check section "identification" - key "reference_file"
            reference_file = variant_calling_option_dict.get('identification', {}).get('reference_file', not_found)
            is_ok_reference_file = False
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False
            else:
                is_ok_reference_file = True

            # check that "reference_file" has to be NONE if "reference_dataset_id" is NONE
            if is_ok_reference_dataset_id and is_ok_reference_file and reference_dataset_id.upper() == 'NONE' and reference_file.upper() != 'NONE':
                error_list.append('*** ERROR: "reference_file" has to be NONE if "reference_dataset_id" is NONE.')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = variant_calling_option_dict.get('identification', {}).get('assembly_software', not_found)
            is_ok_assembly_software = False
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif assembly_software.upper() != 'NONE' and not xlib.check_code(assembly_software, get_extended_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_extended_assembly_software_code_list_text()}; or NONE if a reference is used.')
                OK = False
            else:
                is_ok_assembly_software = True

            # check that "assembly_software" has to be NONE if "reference_dataset_id" is not NONE, and vice versa
            if is_ok_reference_dataset_id and is_ok_assembly_software and (reference_dataset_id.upper() == 'NONE' and assembly_software.upper() == 'NONE' or reference_dataset_id.upper() != 'NONE' and assembly_software.upper() != 'NONE'):
                error_list.append('*** ERROR: "assembly_software" has to be NONE if "reference_dataset_id" is not NONE, and vice versa.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = variant_calling_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            is_ok_assembly_dataset_id = False
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.upper() != 'NONE' and not xlib.check_startswith(assembly_dataset_id, get_extended_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" does not have to start with {get_extended_assembly_software_code_list_text()}.')
                OK = False
            else:
                is_ok_assembly_dataset_id = True

            # check that "assembly_dataset_id" has to be NONE if "assembly_software" is NONE
            if is_ok_assembly_software and is_ok_assembly_dataset_id and assembly_software.upper() == 'NONE' and assembly_dataset_id.upper() != 'NONE':
                error_list.append('*** ERROR: "assembly_dataset_id" has to be NONE if "assembly_software" is NONE.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = variant_calling_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif (assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or assembly_dataset_id.startswith(xlib.get_soapdenovo2_code())) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not (assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or assembly_dataset_id.startswith(xlib.get_soapdenovo2_code())) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} and {xlib.get_soapdenovo2_name()}; or NONE in any other case.')
                    OK = False

            # check section "identification" - key "alignment_software"
            alignment_software = variant_calling_option_dict.get('identification', {}).get('alignment_software', not_found)
            if alignment_software == not_found:
                error_list.append(f'*** ERROR: the key "alignment_software" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_code(alignment_software, get_alignment_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "alignment_software" has to be {get_alignment_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "alignment_dataset_id"
            alignment_dataset_id = variant_calling_option_dict.get('identification', {}).get('alignment_dataset_id', not_found)
            if alignment_dataset_id == not_found:
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_startswith(alignment_dataset_id, get_alignment_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" has to start with {get_alignment_software_code_list_text()}.')
                OK = False

        # check section "Variant calling parameters"
        if 'Variant calling parameters' not in sections_list:
            error_list.append('*** ERROR: the section "Variant calling parameters" is not found.')
            OK = False
        else:

            # check section "Variant calling parameters" - key "vcf-merger"
            vcf_merger = variant_calling_option_dict.get('Variant calling parameters', {}).get('vcf-merger', not_found)
            if vcf_merger == not_found:
                error_list.append('*** ERROR: the key "vcf-merger" is not found in the section "Variant calling parameters".')
                OK = False
            elif not xlib.check_code(vcf_merger, get_vcf_merger_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "vcf-merger" has to be {get_vcf_merger_code_list_text()}.')
                OK = False

        # check section "bcftools call parameters"
        if 'bcftools call parameters' not in sections_list:
            error_list.append('*** ERROR: the section "bcftools call parameters" is not found.')
            OK = False
        else:

            # check section "bcftools call parameters" - key "variants-only"
            variants_only = variant_calling_option_dict.get('bcftools call parameters', {}).get('variants-only', not_found)
            if variants_only == not_found:
                error_list.append('*** ERROR: the key "variants-only" is not found in the section "bcftools call parameters".')
                OK = False
            elif not xlib.check_code(variants_only, get_variants_only_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "variants-only" has to be {get_variants_only_code_list_text()}.')
                OK = False

            # check section "bcftools call parameters" - key "consensus-caller"
            consensus_caller = variant_calling_option_dict.get('bcftools call parameters', {}).get('consensus-caller', not_found)
            is_ok_consensus_caller = False
            if consensus_caller == not_found:
                error_list.append('*** ERROR: the key "consensus-caller" is not found in the section "bcftools call parameters".')
                OK = False
            elif not xlib.check_code(consensus_caller, get_consensus_caller_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "consensus-caller" has to be {get_consensus_caller_code_list_text()}.')
                OK = False
            else:
                is_ok_consensus_caller = True

            # check section "bcftools call parameters" - key "multiallelic-caller"
            multiallelic_caller = variant_calling_option_dict.get('bcftools call parameters', {}).get('multiallelic-caller', not_found)
            is_ok_multiallelic_caller = False
            if multiallelic_caller == not_found:
                error_list.append('*** ERROR: the key "multiallelic-caller" is not found in the section "bcftools call parameters".')
                OK = False
            elif not xlib.check_code(multiallelic_caller, get_multiallelic_caller_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "multiallelic-caller" has to be {get_multiallelic_caller_code_list_text()}.')
                OK = False
            else:
                is_ok_multiallelic_caller = True

            # check if consensus-caller value is not equal to multiallelic-caller value
            if is_ok_consensus_caller and is_ok_multiallelic_caller and consensus_caller == multiallelic_caller:
                error_list.append(f'*** ERROR: The value consensus-caller value ({consensus_caller}) and multiallelic-caller ({multiallelic_caller} have to be different).')
                OK = False

            # check section "bcftools call parameters" - key "other_parameters"
            not_allowed_parameters_list = ['threads', '', '', '']
            other_parameters = variant_calling_option_dict.get('bcftools call parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "bcftools call parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'The {xlib.get_variant_calling_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the esult of the control variables and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_variant_calling_process_script(cluster_name, current_run_dir):
    '''
    Build the current Variant calling process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the transcript-filter option dictionary
    variant_calling_option_dict = xlib.get_option_dict(get_variant_calling_config_file())

    # get the options
    experiment_id = variant_calling_option_dict['identification']['experiment_id']
    reference_dataset_id = variant_calling_option_dict['identification']['reference_dataset_id']
    reference_file = variant_calling_option_dict['identification']['reference_file']
    assembly_software = variant_calling_option_dict['identification']['assembly_software']
    assembly_dataset_id = variant_calling_option_dict['identification']['assembly_dataset_id']
    assembly_type = variant_calling_option_dict['identification']['assembly_type']
    alignment_software = variant_calling_option_dict['identification']['alignment_software']
    alignment_dataset_id = variant_calling_option_dict['identification']['alignment_dataset_id']
    vcf_merger = variant_calling_option_dict['Variant calling parameters']['vcf-merger']
    variants_only = variant_calling_option_dict['bcftools call parameters']['variants-only']
    consensus_caller = variant_calling_option_dict['bcftools call parameters']['consensus-caller']
    multiallelic_caller = variant_calling_option_dict['bcftools call parameters']['multiallelic-caller']
    other_parameters = variant_calling_option_dict['bcftools call parameters']['other_parameters']

    # set the cluster reference file
    if reference_dataset_id.upper() != 'NONE':
        cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    else:
        if assembly_software in [xlib.get_soapdenovotrans_code(), xlib.get_soapdenovo2_code()]:
            if assembly_type == 'CONTIGS':
                cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.contig'
            elif  assembly_type == 'SCAFFOLDS':
                cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.scafSeq'
        elif assembly_software == xlib.get_transabyss_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/transabyss-final.fa'
        elif assembly_software == xlib.get_trinity_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity.fasta'
        elif assembly_software == xlib.get_ggtrinity_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity-GG.fasta'
        elif assembly_software == xlib.get_cd_hit_est_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/clustered-transcriptome.fasta'
        elif assembly_software == xlib.get_transcript_filter_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/filtered-transcriptome.fasta'
        elif assembly_software == xlib.get_starcode_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/starcode.fasta'

    # set the alignment file paths
    if alignment_software == xlib.get_bowtie2_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/alignment.sam'
        bam_files = '$BAM_DIR/alignment.bam'
    elif alignment_software == xlib.get_gsnap_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/*-split.concordant_uniq'
        sam_files_2 = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/*-split.uniq'
        bam_files = '$BAM_DIR/*.bam'
    elif alignment_software == xlib.get_hisat2_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/alignment.sam'
        bam_files = '$BAM_DIR/alignment.bam'
    elif alignment_software == xlib.get_star_code():
        bam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/*-Aligned.sortedByCoord.out.bam'
    elif alignment_software == xlib.get_tophat_code():
        bam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/accepted_hits.bam'

    # set parameters of bcftools call
    bcftools_call_parameter = ''
    if variants_only.upper() == 'YES':
        bcftools_call_parameter = f'{bcftools_call_parameter} --variants-only'
    if consensus_caller.upper() == 'YES':
        bcftools_call_parameter = f'{bcftools_call_parameter} --consensus-caller'
    if multiallelic_caller.upper() == 'YES':
        bcftools_call_parameter = f'{bcftools_call_parameter} --multiallelic-caller'
    if other_parameters.upper() != 'NONE':
        parameter_list = [x.strip() for x in other_parameters.split(';')]
        for i in range(len(parameter_list)):
            if parameter_list[i].find('=') > 0:
                pattern = r'^--(.+)=(.+)$'
                mo = re.search(pattern, parameter_list[i])
                parameter_name = mo.group(1).strip()
                parameter_value = mo.group(2).strip()
                bcftools_call_parameter = f'{bcftools_call_parameter} --{parameter_name}={parameter_value}'
            else:
                pattern = r'^--(.+)$'
                mo = re.search(pattern, parameter_list[i])
                parameter_name = mo.group(1).strip()
                bcftools_call_parameter = f'{bcftools_call_parameter} --{parameter_name}'

    # write the transcript-filter process script
    try:
        if not os.path.exists(os.path.dirname(get_variant_calling_process_script())):
            os.makedirs(os.path.dirname(get_variant_calling_process_script()))
        with open(get_variant_calling_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'export PATH=$MINICONDA3_BIN_PATH:$PATH\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'CURRENT_DIR={current_run_dir}\n')
            script_file_id.write( 'BAM_DIR=$CURRENT_DIR/BAM\n')
            script_file_id.write( 'if [ ! -d "$BAM_DIR" ]; then mkdir --parents $BAM_DIR; fi\n')
            script_file_id.write( 'BED_DIR=$CURRENT_DIR/BED\n')
            script_file_id.write( 'if [ ! -d "$BED_DIR" ]; then mkdir --parents $BED_DIR; fi\n')
            script_file_id.write( 'VCF_DIR=$CURRENT_DIR/VCF\n')
            script_file_id.write( 'if [ ! -d "$VCF_DIR" ]; then mkdir --parents $VCF_DIR; fi\n')
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
            if alignment_software in [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code()]:
                script_file_id.write( 'function convert_sam2bam\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/convert_sam2bam.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Converting SAM files to BAM format ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write(f'        source activate {xlib.get_samtools_anaconda_code()}\n')
                script_file_id.write(f'        ls {sam_files} > sam-files.txt\n')
                if alignment_software == xlib.get_gsnap_code():
                    script_file_id.write(f'        if ! [ -s sam-files.txt ]; then\n')
                    script_file_id.write(f'            ls {sam_files_2} > sam-files.txt\n')
                    script_file_id.write(f'        fi\n')
                script_file_id.write( '        while read FILE_SAM; do\n')
                if alignment_software == xlib.get_gsnap_code():
                    script_file_id.write( '            FILE_BAM=$BAM_DIR/`basename $FILE_SAM`.bam\n')
                else:
                    script_file_id.write( '            FILE_BAM=$BAM_DIR/`basename $FILE_SAM | sed "s|.sam|.bam|g"`\n')
                script_file_id.write( '            samtools view -b -S -o $FILE_BAM $FILE_SAM\n')
                script_file_id.write( '            RC=$?\n')
                script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error samtools-view $RC; fi\n')
                script_file_id.write( '        done < sam-files.txt\n')
                script_file_id.write( '        conda deactivate\n')
                script_file_id.write( '        echo "SAM files are converted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function get_alignment_stats\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/get_alignment_stats.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Getting alignment statistic ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_samtools_anaconda_code()}\n')
            script_file_id.write(f'        ls {bam_files} > bam-files.txt\n')
            script_file_id.write( '        while read FILE_SAM; do\n')
            script_file_id.write( '            FILE_BAM_STATS=$BAM_DIR/`basename $FILE_BAM | sed "s|.bam|.bam-stats.txt|g"`\n')
            script_file_id.write( '            samtools flagstat $FILE_BAM >$FILE_BAM_STATS\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error samtools-flagstat $RC; fi\n')
            script_file_id.write( '        done < sam-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Statistics are got."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function sort_and_index_bam_files\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/sort_and_index_bam_files.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Sorting and indexing BAM files ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_samtools_anaconda_code()}\n')
            script_file_id.write(f'        ls {bam_files} > bam-files.txt\n')
            script_file_id.write( '        while read FILE_BAM; do\n')
            script_file_id.write( '            FILE_SORTED_BAM=$BAM_DIR/`basename $FILE_BAM | sed "s|.bam|.sorted.bam|g"`\n')
            script_file_id.write( '            samtools sort $FILE_BAM -o $FILE_SORTED_BAM\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error samtools-sort $RC; fi\n')
            script_file_id.write( '            samtools index $FILE_SORTED_BAM\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error samtools-index $RC; fi\n')
            script_file_id.write( '        done < bam-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "BAM files are sorted and indexed."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function convert_bam2bed\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/convert_bam2bed.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Converting BAM files to BED format ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_bedtools_anaconda_code()}\n')
            script_file_id.write(f'        ls {bam_files} > bam-files.txt\n')
            script_file_id.write( '        while read FILE_BAM; do\n')
            script_file_id.write( '            FILE_BED=$BED_DIR/`basename $FILE_BAM | sed "s|.bam|.bed|g"`\n')
            script_file_id.write( '            bedtools bamtobed -i $FILE_BAM > $FILE_BED\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error bedtools-bamtobed $RC; fi\n')
            script_file_id.write( '        done < bam-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "BAM files are converted."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function convert_bam2vcf\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/convert_bam2vcf.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Converting sorted BAM files to VCF format ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_bcftools_anaconda_code()}\n')
            script_file_id.write( '        ls $BAM_DIR/*.sorted.bam > sorted-bam-files.txt\n')
            script_file_id.write( '        while read FILE_SORTED_BAM; do\n')
            script_file_id.write( '            FILE_VCF=$VCF_DIR/`basename $FILE_SORTED_BAM | sed "s|.sorted.bam|.vcf|g"`\n')
            script_file_id.write(f'            bcftools mpileup --output-type u --fasta-ref {cluster_reference_file} $FILE_SORTED_BAM | bcftools call {bcftools_call_parameter} - > $FILE_VCF\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error bcftools-mpileup_bcftools-call $RC; fi\n')
            script_file_id.write( '        done < sorted-bam-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "BAM files are converted."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function compress_vcf_files\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/compress_vcf_files.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Compressing VCF files ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        ls $VCF_DIR/*.vcf > vcf-files.txt\n')
            script_file_id.write(f'        source activate {xlib.get_tabix_anaconda_code()}\n')
            script_file_id.write( '        while read FILE_VCF; do\n')
            script_file_id.write( '            bgzip $FILE_VCF\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error bgzip $RC; fi\n')
            script_file_id.write( '        done < vcf-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "VCF files are compressed."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            if vcf_merger.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function index_vcf_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/index_vcf_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Indexing VCF files for the merger ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        FILENUM=`ls -l $VCF_DIR/*.vcf.gz | grep -v ^d | wc -l`\n')
                script_file_id.write( '        echo "File number: $FILENUM"\n')
                script_file_id.write( '        if [ "$FILENUM" -gt 1 ]; then\n')
                script_file_id.write(f'            source activate {xlib.get_tabix_anaconda_code()}\n')
                script_file_id.write( '            parallel tabix -p vcf ::: $VCF_DIR/*.vcf.gz\n')
                script_file_id.write( '            RC=$?\n')
                script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error parallel_tabix $RC; fi\n')
                script_file_id.write( '            conda deactivate\n')
                script_file_id.write( '            echo "VCF files are indexed."\n')
                script_file_id.write( '        else\n')
                script_file_id.write( '           echo "The merger needs because two or more VCF files are necessary. The indexing is not done."\n')
                script_file_id.write( '        fi\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_vcf_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_vcf_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Merging VCF files ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        FILENUM=`ls -l $VCF_DIR/*.vcf.gz | grep -v ^d | wc -l`\n')
                script_file_id.write( '        echo "File number: $FILENUM"\n')
                script_file_id.write( '        if [ "$FILENUM" -gt 1 ]; then\n')
                script_file_id.write(f'            source activate {xlib.get_bcftools_anaconda_code()}\n')
                script_file_id.write( '            bcftools merge --merge all --output-type z $VCF_DIR/*.vcf.gz > $VCF_DIR/merged_samples.vcf.gz\n')
                script_file_id.write( '            RC=$?\n')
                script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error vcfutils.pl $RC; fi\n')
                script_file_id.write( '            conda deactivate\n')
                script_file_id.write( '            echo "VCF files are merged."\n')
                script_file_id.write( '            touch $STEP_STATUS\n')
                script_file_id.write( '        else\n')
                script_file_id.write( '           echo "The merger is not done because two or more VCF files are necessary."\n')
                script_file_id.write( '        fi\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function get_variant_stats\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/get_variant_stats.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Getting variant statistics ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_bcftools_anaconda_code()}\n')
            script_file_id.write( '        ls $VCF_DIR/*.vcf.gz > compressed-vcf-files.txt\n')
            script_file_id.write( '        while read COMPRESSED_FILE_VCF; do\n')
            script_file_id.write( '            STATS_FILE=`echo $COMPRESSED_FILE_VCF | sed "s|.vcf.gz|.vcf-stats.txt|g"`\n')
            script_file_id.write( '            PLOT_DIR=`echo $COMPRESSED_FILE_VCF | sed "s|.vcf.gz|-stats-plots|g"`\n')
            script_file_id.write( '            bcftools stats $COMPRESSED_FILE_VCF > $STATS_FILE\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error bcftools-stats $RC; fi\n')
            script_file_id.write( '            plot-vcfstats --prefix $PLOT_DIR $STATS_FILE\n')
            script_file_id.write( '            RC=$?\n')
            script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error plot-vcfstats $RC; fi\n')
            script_file_id.write( '        done < compressed-vcf-files.txt\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Statistics are got."\n')
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
            process_name = f'{xlib.get_variant_calling_name()} process'
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
            if alignment_software in [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code()]:
                script_file_id.write( 'convert_sam2bam\n')
            script_file_id.write( 'get_alignment_stats\n')
            script_file_id.write( 'sort_and_index_bam_files\n')
            script_file_id.write( 'convert_bam2bed\n')
            script_file_id.write( 'convert_bam2vcf\n')
            script_file_id.write( 'compress_vcf_files\n')
            if vcf_merger.upper() == 'YES':
                script_file_id.write( 'index_vcf_files\n')
                script_file_id.write( 'merge_vcf_files\n')
            script_file_id.write( 'get_variant_stats\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_variant_calling_process_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_variant_calling_process_starter(current_run_dir):
    '''
    Build the starter of the current Variant calling process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Variant calling process starter
    try:
        if not os.path.exists(os.path.dirname(get_variant_calling_process_starter())):
            os.makedirs(os.path.dirname(get_variant_calling_process_starter()))
        with open(get_variant_calling_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_variant_calling_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_variant_calling_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_variant_calling_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a Variant calling process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_variant_calling_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_variant_calling_process_starter()), log)

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

def get_variant_calling_config_file():
    '''
    Get the Variant calling config file path.
    '''

    # assign the Variant calling config file path
    variant_calling_config_file = f'{xlib.get_config_dir()}/{xlib.get_variant_calling_code()}-config.txt'

    # return the Variant calling config file path
    return variant_calling_config_file

#-------------------------------------------------------------------------------

def get_variant_calling_process_script():
    '''
    Get the Variant calling process script path in the local computer.
    '''

    # assign the Variant calling script path
    variant_calling_process_script = f'{xlib.get_temp_dir()}/{xlib.get_variant_calling_code()}-process.sh'

    # return the Variant calling script path
    return variant_calling_process_script

#-------------------------------------------------------------------------------

def get_variant_calling_process_starter():
    '''
    Get the Variant calling process starter path in the local computer.
    '''

    # assign the Variant calling process starter path
    variant_calling_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_variant_calling_code()}-process-starter.sh'

    # return the Variant calling starter path
    return variant_calling_process_starter

#-------------------------------------------------------------------------------
    
def get_technique_code_list():
    '''
    Get the code list of "technique".
    '''

    return ['IND1', 'IND1_DBR', 'IND1_IND2', 'IND1_IND2_DBR']

#-------------------------------------------------------------------------------
    
def get_technique_code_list_text():
    '''
    Get the code list of "technique" as text.
    '''

    return 'IND1 (only index1) or IND1_DBR (index1 + DBR) or IND1_IND2 (index1 + index2) or IND1_IND2_DBR (index1 + index2 + DBR)'

#-------------------------------------------------------------------------------
    
def get_format_code_list():
    '''
    Get the code list of "format".
    '''

    return ['FASTA', 'FASTQ']

#-------------------------------------------------------------------------------
    
def get_format_code_list_text():
    '''
    Get the code list of "format" as text.
    '''

    return str(get_format_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_read_type_code_list():
    '''
    Get the code list of "readtype".
    '''

    return ['SE', 'PE']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "readtype" as text.
    '''

    return 'SE (single-end) or PE (pair-end)'

#-------------------------------------------------------------------------------
    
def get_pcrdistribution_code_list():
    '''
    Get the code list of "pcrdistribution".
    '''

    return ['MULTINOMIAL', 'POISSON']

#-------------------------------------------------------------------------------
    
def get_pcrdistribution_code_list_text():
    '''
    Get the code list of "pcrdistribution" as text.
    '''

    return str(get_pcrdistribution_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_extended_assembly_software_code_list():
    '''
    Get the code list of "assembly_software".
    '''

    return [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(),  xlib.get_transcript_filter_code(), xlib.get_soapdenovo2_code(), xlib.get_starcode_name()]

#-------------------------------------------------------------------------------

def get_extended_assembly_software_code_list_text():
    '''
    Get the code list of "assembly_software" as text.
    '''

    return f'{xlib.get_soapdenovotrans_code()} ({xlib.get_soapdenovotrans_name()}) or {xlib.get_transabyss_code()} ({xlib.get_transabyss_name()}) or {xlib.get_trinity_code()} ({xlib.get_trinity_name()}) or {xlib.get_ggtrinity_code()} ({xlib.get_ggtrinity_name()}) or {xlib.get_cd_hit_est_code()} ({xlib.get_cd_hit_est_name()}) or {xlib.get_transcript_filter_code()} ({xlib.get_transcript_filter_name()}) or {xlib.get_soapdenovo2_code()} ({xlib.get_soapdenovo2_name()}) or {xlib.get_starcode_code()} ({xlib.get_starcode_name()})'

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list():
    '''
    Get the code list of "alignment_software".
    '''

    return [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code(), xlib.get_star_code(), xlib.get_tophat_code()]

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list_text():
    '''
    Get the code list of "alignment_software" as text.
    '''

    return f'{xlib.get_bowtie2_code()} ({xlib.get_bowtie2_name()}) or {xlib.get_gsnap_code()} ({xlib.get_gsnap_name()}) or {xlib.get_hisat2_code()} ({xlib.get_hisat2_name()}) or {xlib.get_star_code()} ({xlib.get_star_name()}) or {xlib.get_tophat_code()} ({xlib.get_tophat_name()})'

#-------------------------------------------------------------------------------
    
def get_vcf_merger_code_list():
    '''
    Get the code list of "vcf-merger".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_vcf_merger_code_list_text():
    '''
    Get the code list of "variants-only" as text.
    '''

    return str(get_vcf_merger_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_variants_only_code_list():
    '''
    Get the code list of "variants-only".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_variants_only_code_list_text():
    '''
    Get the code list of "variants-only" as text.
    '''

    return str(get_variants_only_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_consensus_caller_code_list():
    '''
    Get the code list of "consensus-caller".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_consensus_caller_code_list_text():
    '''
    Get the code list of "consensus-caller" as text.
    '''

    return str(get_consensus_caller_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_multiallelic_caller_code_list():
    '''
    Get the code list of "multiallelic-caller".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_multiallelic_caller_code_list_text():
    '''
    Get the code list of "multiallelic_caller" as text.
    '''

    return str(get_multiallelic_caller_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the ddRADseqTools process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
