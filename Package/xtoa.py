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
This file contains functions related to the TOA (Tree-oriented Annotation) process used in both
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

def is_installed_toa(cluster_name, passed_connection, ssh_client):
    '''
    Check if TOA is installed.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the installation control variable
    is_installed = False

    # get the aplication directory in the cluster
    cluster_app_dir = xlib.get_cluster_app_dir()

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        if not OK:
            for error in error_list:
                error_list.append(f'{error}\n')
                OK = False

    # check the TOA directory is created
    if OK:
        command = f'[ -d {cluster_app_dir}/{xlib.get_toa_name()} ] && echo RC=0 || echo RC=1'
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

def install_toa(cluster_name, log, function=None):
    '''
    Install the TOA software in the cluster.
    '''

    # initialize the control variable
    OK = True

    # get the aplication directory in the cluster
    cluster_app_dir = xlib.get_cluster_app_dir()

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
        command = f'[ -d {cluster_app_dir} ] && echo RC=0 || echo RC=1'.format()
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write(f'You have to link a volume in the mounting point {cluster_app_dir} for the cluster {cluster_name}.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_installation_dir(), xlib.get_toa_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TOA installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_toa_installation_script()} ...\n')
        (OK, error_list) = build_toa_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the TOA installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_toa_installation_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_toa_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_toa_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TOA installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_toa_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_toa_installation_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TOA installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_toa_installation_starter()} ...\n')
        (OK, error_list) = build_toa_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the TOA installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_toa_installation_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_toa_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_toa_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TOA installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_toa_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_toa_installation_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the TOA installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_toa_installation_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_toa_installation_starter()), log)

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

def build_toa_installation_script(cluster_name, current_run_dir):
    '''
    Build the TOA installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of TOA
    (_, toa_url, _) = xconfiguration.get_bioinfo_app_data(xlib.get_toa_name())

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_toa_installation_script())):
            os.makedirs(os.path.dirname(get_toa_installation_script()))
        with open(get_toa_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function remove_toa_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_toa_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_toa_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_toa_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_toa_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_toa_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            download_script = f'import requests; r = requests.get(\'{toa_url}\') ; open(\'{xlib.get_toa_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $MINICONDA3_BIN_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write( '    echo "The file is downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function decompress_toa_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Decompressing the {xlib.get_toa_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    unzip {xlib.get_toa_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write( '    echo "The file is decompressed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_toa_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Renaming the {xlib.get_toa_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    mv {xlib.get_toa_name()}-master {xlib.get_toa_name()}\n')
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
            script_file_id.write(f'    chmod u+x {xlib.get_toa_name()}/Package/*.py\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi\n')
            script_file_id.write( '    echo "The directory is renamed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_toa_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_toa_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_toa_name()}.zip\n')
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
            process_name = f'{xlib.get_toa_name()} installation'
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
            script_file_id.write( 'remove_toa_directory\n')
            script_file_id.write( 'download_toa_installation_file\n')
            script_file_id.write( 'decompress_toa_installation_file\n')
            script_file_id.write( 'rename_toa_directory\n')
            script_file_id.write( 'set_execution_permissions\n')
            script_file_id.write( 'remove_toa_installation_file\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_toa_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_toa_installation_starter(current_run_dir):
    '''
    Build the starter of the TOA installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_toa_installation_starter())):
            os.makedirs(os.path.dirname(get_toa_installation_starter()))
        with open(get_toa_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_toa_installation_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_toa_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_toa_installation_script():
    '''
    Get the TOA installation script path in the local computer.
    '''

    # assign the script path
    toa_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_name()}-installation.sh'

    # return the script path
    return toa_installation_script

#-------------------------------------------------------------------------------

def get_toa_installation_starter():
    '''
    Get the TOA installation starter path in the local computer.
    '''

    # assign the starter path
    toa_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_name()}-installation-starter.sh'

    # return the starter path
    return toa_installation_starter

#-------------------------------------------------------------------------------

def create_toa_config_file(toa_dir=f'{xlib.get_cluster_app_dir()}/TOA/Package', miniconda3_dir=f'{xlib.get_cluster_app_dir()}/Miniconda3', db_dir=xlib.get_cluster_database_dir(), result_dir=xlib.get_cluster_result_dir()):
    '''
    Create the TOA config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the Miniconda3 bin and envs directories
    miniconda3_bin_dir = f'{miniconda3_dir}/bin'
    miniconda3_envs_dir = f'{miniconda3_dir}/envs'

    # get the NGScloud config file
    toa_config_file = get_toa_config_file()

    # create the TOA config file
    if OK:
        try:
            if not os.path.exists(os.path.dirname(toa_config_file)):
                os.makedirs(os.path.dirname(toa_config_file))
            with open(toa_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
                file_id.write( '# environment\n')
                file_id.write(f'TOA_DIR={toa_dir}\n')
                file_id.write(f'MINICONDA3_DIR={miniconda3_dir}\n')
                file_id.write(f'MINICONDA3_BIN_DIR={miniconda3_bin_dir}\n')
                file_id.write(f'MINICONDA3_ENVS_DIR={miniconda3_envs_dir}\n')
                file_id.write(f'RESULT_DIR={result_dir}\n')
                file_id.write(f'DB_DIR={db_dir}\n')
                file_id.write(f'TOA_DB_DIR={db_dir}/TOA\n')
                file_id.write(f'DATA_DIR={db_dir}/data\n')
                file_id.write(f'PLAZA_DIR={db_dir}/PLAZA\n')
                file_id.write(f'NCBI_DIR={db_dir}/NCBI\n')
                file_id.write(f'INTERPRO_DIR={db_dir}/InterPro\n')
                file_id.write(f'GO_DIR={db_dir}/GO\n')
                file_id.write(f'EC_DIR={db_dir}/EC\n')
                file_id.write(f'KEGG_DIR={db_dir}/KEGG\n')
                file_id.write( '\n')
                file_id.write( '# TOA database\n')
                file_id.write(f'TOA_DB={db_dir}/TOA/toa.db\n')
                file_id.write( '\n')
                file_id.write( '# basic data\n')
                file_id.write(f'DATASET_FILE={db_dir}/data/datasets.txt\n')
                file_id.write(f'SPECIES_FILE={db_dir}/data/species.txt\n')
                file_id.write( '\n')
                file_id.write( '# other basic data\n')
                file_id.write( '#EC_IDS_FTP=https://www.enzyme-database.org/downloads/enzyme-data.sql.gz\n')
                file_id.write(f'#EC_IDS_FILE={db_dir}/EC/enzyme-data.sql.gz\n')
                file_id.write( 'EC_IDS_FTP=ftp://ftp.expasy.org/databases/enzyme/enzyme.dat\n')
                file_id.write(f'EC_IDS_FILE={db_dir}/EC/enzyme.dat\n')
                file_id.write( 'KEGG_IDS_FTP=ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz\n')
                file_id.write(f'KEGG_IDS_FILE={db_dir}/KEGG/ko_list.gz\n')
                file_id.write( '\n')
                file_id.write( '# Gymno PLAZA 1.0\n')
                file_id.write( 'GYMNO_01_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Fasta/cds.csv.gz\n')
                file_id.write(f'GYMNO_01_CDS_FILE={db_dir}/PLAZA/gymno_01-cds.fasta.gz\n')
                file_id.write( 'GYMNO_01_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Fasta/proteome.csv.gz\n')
                file_id.write(f'GYMNO_01_PROTEOME_FILE={db_dir}/PLAZA/gymno_01-proteome.fasta.gz\n')
                file_id.write( 'GYMNO_01_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Descriptions\n')
                file_id.write(f'GYMNO_01_GENEDESC_DIR={db_dir}/PLAZA/gymno_01-descriptions\n')
                file_id.write( 'GYMNO_01_GENEDESC_FILE_PATTERN="gene_description.*.csv.gz"\n')
                file_id.write( 'GYMNO_01_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/InterPro/interpro.csv.gz\n')
                file_id.write(f'GYMNO_01_INTERPRO_FILE={db_dir}/PLAZA/gymno_01-interpro.csv.gz\n')
                file_id.write( 'GYMNO_01_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/GO/go.csv.gz\n')
                file_id.write(f'GYMNO_01_GO_FILE={db_dir}/PLAZA/gymno_01-go.csv.gz\n')
                file_id.write( 'GYMNO_01_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/MapMan/mapman.csv.gz\n')
                file_id.write(f'GYMNO_01_MAPMAN_FILE={db_dir}/PLAZA/gymno_01-mapman.csv.gz\n')
                file_id.write( 'GYMNO_01_BLASTPLUS_DB_NAME=gymno_01\n')
                file_id.write(f'GYMNO_01_BLASTPLUS_DB_DIR={db_dir}/PLAZA/gymno_01-blastplus-db\n')
                file_id.write(f'GYMNO_01_BLASTPLUS_DB_FILE={db_dir}/PLAZA/gymno_01-blastplus-db/gymno_01\n')
                file_id.write(f'GYMNO_01_DIAMOND_DB_DIR={db_dir}/PLAZA/gymno_01-diamond-db\n')
                file_id.write(f'GYMNO_01_DIAMOND_DB_FILE={db_dir}/PLAZA/gymno_01-diamond-db/gymno_01\n')
                file_id.write(f'GYMNO_01_BLAST_XML=$OUTPUT_DIR/gymno_01-alignment.xml\n')
                file_id.write( 'GYMNO_01_ANNOTATION_FILE=$OUTPUT_DIR/gymno_01-annotation.csv\n')
                file_id.write( 'GYMNO_01_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/gymno_01-nonann-transcripts.fasta\n')
                file_id.write( 'GYMNO_01_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/gymno_01-nonann-peptides.fasta\n')
                file_id.write( '\n')
                file_id.write( '# Dicots PLAZA 4.0\n')
                file_id.write( 'DICOTS_04_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/cds.all_transcripts.fasta.gz\n')
                file_id.write(f'DICOTS_04_CDS_FILE={db_dir}/PLAZA/dicots_04-cds.fasta.gz\n')
                file_id.write( 'DICOTS_04_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/proteome.all_transcripts.fasta.gz\n')
                file_id.write(f'DICOTS_04_PROTEOME_FILE={db_dir}/PLAZA/dicots_04-proteome.fasta.gz\n')
                file_id.write( 'DICOTS_04_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Descriptions\n')
                file_id.write(f'DICOTS_04_GENEDESC_DIR={db_dir}/PLAZA/dicots_04-descriptions\n')
                file_id.write( 'DICOTS_04_GENEDESC_FILE_PATTERN="gene_description.*.csv.gz"\n')
                file_id.write( 'DICOTS_04_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/InterPro/interpro.csv.gz\n')
                file_id.write(f'DICOTS_04_INTERPRO_FILE={db_dir}/PLAZA/dicots_04-interpro.csv.gz\n')
                file_id.write( 'DICOTS_04_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/GO/go.csv.gz\n')
                file_id.write(f'DICOTS_04_GO_FILE={db_dir}/PLAZA/dicots_04-go.csv.gz\n')
                file_id.write( 'DICOTS_04_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/MapMan/mapman.csv.gz\n')
                file_id.write(f'DICOTS_04_MAPMAN_FILE={db_dir}/PLAZA/dicots_04-mapman.csv.gz\n')
                file_id.write( 'DICOTS_04_BLASTPLUS_DB_NAME=dicots_04\n')
                file_id.write(f'DICOTS_04_BLASTPLUS_DB_DIR={db_dir}/PLAZA/dicots_04-blastplus-db\n')
                file_id.write(f'DICOTS_04_BLASTPLUS_DB_FILE={db_dir}/PLAZA/dicots_04-blastplus-db/dicots_04\n')
                file_id.write(f'DICOTS_04_DIAMOND_DB_DIR={db_dir}/PLAZA/dicots_04-diamond-db\n')
                file_id.write(f'DICOTS_04_DIAMOND_DB_FILE={db_dir}/PLAZA/dicots_04-diamond-db/dicots_04\n')
                file_id.write( 'DICOTS_04_BLAST_XML=$OUTPUT_DIR/dicots_04-alignment.xml\n')
                file_id.write( 'DICOTS_04_ANNOTATION_FILE=$OUTPUT_DIR/dicots_04-annotation.csv\n')
                file_id.write( 'DICOTS_04_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/dicots_04-nonann-transcripts.fasta\n')
                file_id.write( 'DICOTS_04_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/dicots_04-nonann-peptides.fasta\n')
                file_id.write( '\n')
                file_id.write( '# Monocots PLAZA 4.0\n')
                file_id.write( 'MONOCOTS_04_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Fasta/cds.all_transcripts.fasta.gz\n')
                file_id.write(f'MONOCOTS_04_CDS_FILE={db_dir}/PLAZA/monocots_04-cds.fasta.gz\n')
                file_id.write( 'MONOCOTS_04_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Fasta/proteome.all_transcripts.fasta.gz\n')
                file_id.write(f'MONOCOTS_04_PROTEOME_FILE={db_dir}/PLAZA/monocots_04-proteome.fasta.gz\n')
                file_id.write( 'MONOCOTS_04_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Descriptions\n')
                file_id.write(f'MONOCOTS_04_GENEDESC_DIR={db_dir}/PLAZA/monocots_04-descriptions\n')
                file_id.write( 'MONOCOTS_04_GENEDESC_FILE_PATTERN="gene_description.*.csv.gz"\n')
                file_id.write( 'MONOCOTS_04_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/InterPro/interpro.csv.gz\n')
                file_id.write(f'MONOCOTS_04_INTERPRO_FILE={db_dir}/PLAZA/monocots_04-interpro.csv.gz\n')
                file_id.write( 'MONOCOTS_04_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/GO/go.csv.gz\n')
                file_id.write(f'MONOCOTS_04_GO_FILE={db_dir}/PLAZA/monocots_04-go.csv.gz\n')
                file_id.write( 'MONOCOTS_04_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/MapMan/mapman.csv.gz\n')
                file_id.write(f'MONOCOTS_04_MAPMAN_FILE={db_dir}/PLAZA/monocots_04-mapman.csv.gz\n')
                file_id.write( 'MONOCOTS_04_BLASTPLUS_DB_NAME=monocots_04\n')
                file_id.write(f'MONOCOTS_04_BLASTPLUS_DB_DIR={db_dir}/PLAZA/monocots_04-blasplus-db\n')
                file_id.write(f'MONOCOTS_04_BLASTPLUS_DB_FILE={db_dir}/PLAZA/monocots_04-blasplus-db/monocots_04\n')
                file_id.write(f'MONOCOTS_04_DIAMOND_DB_DIR={db_dir}/PLAZA/monocots_04-diamond-db\n')
                file_id.write(f'MONOCOTS_04_DIAMOND_DB_FILE={db_dir}/PLAZA/monocots_04-diamond-db/monocots_04\n')
                file_id.write( 'MONOCOTS_04_BLAST_XML=$OUTPUT_DIR/monocots_04-alignment.xml\n')
                file_id.write( 'MONOCOTS_04_ANNOTATION_FILE=$OUTPUT_DIR/monocots_04-annotation.csv\n')
                file_id.write( 'MONOCOTS_04_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/monocots_04-nonann-transcripts.fasta\n')
                file_id.write( 'MONOCOTS_04_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/monocots_04-nonann-peptides.fasta\n')
                file_id.write( '\n')
                file_id.write( '# NCBI RefSeq Plant\n')
                file_id.write( 'REFSEQ_PLANT_FTP=ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/\n')
                file_id.write(f'REFSEQ_PLANT_LOCAL={db_dir}/NCBI/ftp.ncbi.nlm.nih.gov/refseq/release/plant/\n')
                file_id.write( 'REFSEQ_PROTEIN_FILE_PATTERN=".protein.faa.gz"\n')
                file_id.write(f'REFSEQ_PLANT_PROTEOME_FILE={db_dir}/NCBI/refseq_plant-proteome.fasta\n')
                file_id.write( 'REFSEQ_PLANT_BLASTPLUS_DB_NAME=refseq_plant\n')
                file_id.write(f'REFSEQ_PLANT_BLASTPLUS_DB_DIR={db_dir}/NCBI/refseq_plant-blastplus-db\n')
                file_id.write(f'REFSEQ_PLANT_BLASTPLUS_DB_FILE={db_dir}/NCBI/refseq_plant-blastplus-db/refseq_plant\n')
                file_id.write(f'REFSEQ_PLANT_DIAMOND_DB_DIR={db_dir}/NCBI/refseq_plant-diamond-db\n')
                file_id.write(f'REFSEQ_PLANT_DIAMOND_DB_FILE={db_dir}/NCBI/refseq_plant-diamond-db/refseq_plant\n')
                file_id.write(f'REFSEQ_PLANT_FILE_LIST={db_dir}/NCBI/refseq_plant-file-list.txt\n')
                file_id.write( 'REFSEQ_PLANT_BLAST_XML=$OUTPUT_DIR/refseq_plant-alignment.xml\n')
                file_id.write( 'REFSEQ_PLANT_ANNOTATION_FILE=$OUTPUT_DIR/refseq_plant-annotation.csv\n')
                file_id.write( 'REFSEQ_PLANT_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/refseq_plant-nonann-transcripts.fasta\n')
                file_id.write( 'REFSEQ_PLANT_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/refseq_plant-nonann-peptides.fasta\n')
                file_id.write( '\n')
                file_id.write( '# NCBI Taxonomy\n')
                file_id.write( 'TAXONOMY_TAXDMP_FTP=ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip\n')
                file_id.write(f'TAXONOMY_TAXDMP_FILE={db_dir}/NCBI/taxdmp.zip\n')
                file_id.write(f'TAXONOMY_TAXONNODES_FILE={db_dir}/NCBI/nodes.dmp\n')
                file_id.write(f'TAXONOMY_TAXONNAMES_FILE={db_dir}/NCBI/names.dmp\n')
                file_id.write( 'TAXONOMY_PROTACCESSION_2_TAXID_FTP=ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz\n')
                file_id.write(f'TAXONOMY_PROTACCESSION_2_TAXID_FILE={db_dir}/NCBI/prot.accession2taxid.gz\n')
                file_id.write( 'VIRIDIPLANTAE_TAXID=33090\n')
                file_id.write(f'VIRIDIPLANTAE_TAXID_LIST_FILE={db_dir}/NCBI/taxids_$VIRIDIPLANTAE_TAXID.txt\n')
                file_id.write( '\n')
                file_id.write( '# NCBI BLAST database\n')
                file_id.write( 'BLAST_DATABASES_FTP=ftp://ftp.ncbi.nlm.nih.gov/blast/db/\n')
                file_id.write(f'BLAST_DATABASES_LOCAL={db_dir}/NCBI/ftp.ncbi.nlm.nih.gov/blast/db/\n')
                file_id.write( '\n')
                file_id.write( '# NCBI BLAST database NT\n')
                file_id.write( 'NT_FILE_PATTERN="nt.*.tar.gz"\n')
                file_id.write(f'NT_FILE_LIST={db_dir}/NCBI/nt-file-list.txt\n')
                file_id.write( 'NT_BLASTPLUS_DB_NAME=nt\n')
                file_id.write(f'NT_BLASTPLUS_DB_DIR={db_dir}/NCBI/nt-blastplus-db\n')
                file_id.write(f'NT_BLASTPLUS_DB_FILE={db_dir}/NCBI/nt-blastplus-db/nt\n')
                file_id.write( 'NT_BLAST_XML=$OUTPUT_DIR/nt-alignment.xml\n')
                file_id.write( 'NT_VIRIDIPLANTAE_ANNOTATION_FILE=$OUTPUT_DIR/nt-viridiplantae-annotation.csv\n')
                file_id.write( 'NT_CONTAMINATION_ANNOTATION_FILE=$OUTPUT_DIR/nt-contamination-annotation.csv\n')
                file_id.write( 'NT_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/nt-nonann-transcripts.fasta\n')
                file_id.write( '\n')
                file_id.write( '# NCBI Nucleotide GenInfo identifier lists\n')
                file_id.write(f'NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST={db_dir}/NCBI/nucleotide_viridiplantae.gi\n')
                file_id.write( '\n')
                file_id.write( '# NCBI BLAST database NR\n')
                file_id.write( 'NR_PROTEOME_FTP=ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz\n')
                file_id.write(f'NR_PROTEOME_FILE={db_dir}/NCBI/nr.gz\n')
                file_id.write( 'NR_FILE_PATTERN="nr.*.tar.gz"\n')
                file_id.write(f'NR_FILE_LIST={db_dir}/NCBI/nr-file-list.txt\n')
                file_id.write( 'NR_BLASTPLUS_DB_NAME=nr\n')
                file_id.write(f'NR_BLASTPLUS_DB_DIR={db_dir}/NCBI/nr-blastplus-db\n')
                file_id.write(f'NR_BLASTPLUS_DB_FILE={db_dir}/NCBI/nr-blastplus-db/nr\n')
                file_id.write(f'NR_DIAMOND_DB_DIR={db_dir}/NCBI/nr-diamond-db\n')
                file_id.write(f'NR_DIAMOND_DB_FILE={db_dir}/NCBI/nr-diamond-db/nr\n')
                file_id.write( 'NR_BLAST_XML=$OUTPUT_DIR/nr-alignment.xml\n')
                file_id.write( 'NR_VIRIDIPLANTAE_ANNOTATION_FILE=$OUTPUT_DIR/nr-viridiplantae-annotation.csv\n')
                file_id.write( 'NR_CONTAMINATION_ANNOTATION_FILE=$OUTPUT_DIR/nr-contamination-annotation.csv\n')
                file_id.write( 'NR_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/nr-nonann-peptides.fasta\n')
                file_id.write( '\n')
                file_id.write( '# NCBI Protein GenInfo identifier lists\n')
                file_id.write(f'PROTEIN_VIRIDIPLANTAE_GI_LIST={db_dir}/NCBI/protein_viridiplantae.gi\n')
                file_id.write( '\n')
                file_id.write( '# NCBI Gene\n')
                file_id.write( 'GENE_GENE2GO_FTP=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz\n')
                file_id.write(f'GENE_GENE2GO_FILE={db_dir}/NCBI/gene-gene2go.gz\n')
                file_id.write( 'GENE_GENE2REFSEQ_FTP=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz\n')
                file_id.write(f'GENE_GENE2REFSEQ_FILE={db_dir}/NCBI/gene-gene2refseq.gz\n')
                file_id.write( '\n')
                file_id.write( '# InterPro\n')
                file_id.write( 'INTERPRO_INTERPRO2GO_FTP=ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go\n')
                file_id.write(f'INTERPRO_INTERPRO2GO_FILE={db_dir}/InterPro/interpro2go\n')
                file_id.write( '\n')
                file_id.write( '# Gene Onlogolgy\n')
                file_id.write( 'GO_ONTOLOGY_FTP=http://purl.obolibrary.org/obo/go.obo\n')
                file_id.write(f'GO_ONTOLOGY_FILE={db_dir}/GO/go.obo\n')
                file_id.write( 'GO_EC2GO_FTP=http://geneontology.org/external2go/ec2go\n')
                file_id.write( '#GO_EC2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/ec2go\n')
                file_id.write(f'GO_EC2GO_FILE={db_dir}/GO/ec2go.txt\n')
                file_id.write( 'GO_KEGG2GO_FTP=http://geneontology.org/external2go/kegg2go\n')
                file_id.write( '#GO_KEGG2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/kegg2go\n')
                file_id.write(f'GO_KEGG2GO_FILE={db_dir}/GO/kegg2go.txt\n')
                file_id.write( 'GO_METACYC2GO_FTP=http://geneontology.org/external2go/metacyc2go\n')
                file_id.write( '#GO_METACYC2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/metacyc2go\n')
                file_id.write(f'GO_METACYC2GO_FILE={db_dir}/GO/metacyc2go.txt\n')
                file_id.write( 'GO_INTERPRO2GO_FTP=http://geneontology.org/external2go/interpro2go\n')
                file_id.write( '#GO_INTERPRO2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/interpro2go\n')
                file_id.write(f'GO_INTERPRO2GO_FILE={db_dir}/GO/interpro2go.txt\n')
                file_id.write( '\n')
                file_id.write( '# other transcriptome files\n')
                file_id.write( 'REIDENTIFIED_TRANSCRIPTOME_FILE=$OUTPUT_DIR/reidentified-transcriptome.fasta\n')
                file_id.write( 'TOA_TRANSCRIPTOME_RELATIONSHIP_FILE=$OUTPUT_DIR/toa_transcriptome_relationships.csv\n')
                file_id.write( 'PURGED_TRANSCRIPTOME_FILE=$OUTPUT_DIR/purged-transcriptome.fasta\n')
                file_id.write( '\n')
                file_id.write( '# peptide sequence files\n')
                file_id.write( 'TRANSDECODER_OUTPUT_DIR=$OUTPUT_DIR/transdecoder\n')
                file_id.write( 'ORF_FILE=$TRANSDECODER_OUTPUT_DIR/longest_orfs.pep\n')
                file_id.write( 'PEPTIDE_FILE=$OUTPUT_DIR/`basename "$REIDENTIFIED_TRANSCRIPTOME_FILE"`.transdecoder.pep\n')
                file_id.write( 'REIDENTIFIED_PEPTIDE_FILE=$OUTPUT_DIR/reidentified-peptides.fasta\n')
                file_id.write( 'TOA_TRANSDECODER_RELATIONSHIP_FILE=$OUTPUT_DIR/toa_transdecoder_relationships.csv\n')
                file_id.write( '\n')
                file_id.write( '# merger files\n')
                file_id.write( 'MERGED_ANNOTATION_FILE=$OUTPUT_DIR/merged-annotation.csv\n')
                file_id.write( 'PLANT_ANNOTATION_FILE=$OUTPUT_DIR/plant-annotation.csv\n')
                file_id.write( 'MERGED_BLAST_XML=$OUTPUT_DIR/merged-alignment.xml\n')
                # -- file_id.write( 'RESTOREDIDS_MERGED_BLAST_XML=$OUTPUT_DIR/restoredids-merged-alignment.xml\n')
                file_id.write( '\n')
                file_id.write( '# statistics\n')
                file_id.write( 'STATS_SUBDIR_NAME=stats\n')
                file_id.write( 'STATS_DIR=$OUTPUT_DIR/stats\n')
                file_id.write( 'STATS_BASE_NAME=stats\n')
                file_id.write( 'ANNOTATION_STATS_FILE=$OUTPUT_DIR/stats/stats.csv\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {toa_config_file} can not be created.')
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_toa_config_file():
    '''
    Get the path of the TOA config file corresponding to the environment.
    '''

    # assign the config file
    toa_config_file = f'{xlib.get_config_dir()}/{xlib.get_toa_code()}-config.txt'

    # return the config file
    return toa_config_file

#-------------------------------------------------------------------------------

def get_toa_config_dict():
    '''
    Get the dictionary of TOA configuration.
    '''

    # initialize the dictionary of TOA configuration
    toa_config_dict = {}

    # open the TOA config file
    try:
        toa_config_file_id = open(get_toa_config_file(), mode='r', encoding='iso-8859-1')
    except Exception:
        raise xlib.ProgramException('F001', get_toa_config_file())

    # read the first record
    record = toa_config_file_id.readline()

    # while there are records
    while record != '':

        # process data records
        if not record.lstrip().startswith('#') and record.strip() != '':
            equal_position = record.find('=')
            if equal_position != -1:
                key = record[:equal_position].strip()
                value = record[equal_position + 1:].strip()
            else:
                pass
            toa_config_dict[key] = value

        # read the next record
        record = toa_config_file_id.readline()

    # return the dictionary of TOA configuration
    return toa_config_dict

#-------------------------------------------------------------------------------

def create_dataset_file():
    '''
    Create the file of datasets with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of datasets and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_dataset_file())):
            os.makedirs(os.path.dirname(get_dataset_file()))
        with open(get_dataset_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains the data of the genomic datasets used by Taxonomy-oriented Annotation (TOA) software package.\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: "dataset_id";"dataset_name";"repository_id";"ftp_adress"\n')
            file_id.write( '\n')
            file_id.write( '"dicots_04";"Dicots PLAZA 4.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/"\n')
            file_id.write( '"gene";"Gene";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"\n')
            file_id.write( '"gymno_01";"Gymno PLAZA 1.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/"\n')
            file_id.write( '"monocots_04";"Monocots PLAZA 4.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/"\n')
            file_id.write( '"nr";"nr";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"\n')
            file_id.write( '"nt";"nt";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"\n')
            file_id.write( '"refseq_plant";"RefSeq Plant";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/"\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {get_dataset_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_dataset_file(strict):
    '''
    Check the file of datasets.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the data records
    # format: "dataset_id";"dataset_name";"repository_id";"ftp_adress"
    record_pattern = re.compile(r'^"(.*)";"(.*)";"(.*)";"(.*)"$')

    # open the file of datasets
    try:
        dataset_file_id = open(get_dataset_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dataset_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = dataset_file_id.readline()

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try: 
                    mo = record_pattern.match(record)
                    # -- dataset_id = mo.group(1).strip()
                    # -- dataset_name = mo.group(2).strip()
                    # -- repository_id = mo.group(3).strip()
                    # -- ftp_adress = mo.group(4).strip()
                except Exception as e:
                    record_text = record.replace("\n", "")
                    error_list.append(f'*** ERROR: There is a format error in the record: {record_text}.')
                    OK = False
                    break

            # read the next record
            record = dataset_file_id.readline()

    # close the file of datasets
    dataset_file_id.close()

    # warn that the file of datasets is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_dataset_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dataset_file():
    '''
    Get the dataset file path.
    '''

    # assign the dataset file path
    dataset_file = f'{xlib.get_config_dir()}/datasets.txt'

    # return the dataset file path
    return dataset_file

#-------------------------------------------------------------------------------

def create_species_file():
    '''
    Create the file of species with the default data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the file of species and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_species_file())):
            os.makedirs(os.path.dirname(get_species_file()))
        with open(get_species_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# This file contains the data of species used by Taxonomy-oriented Annotation (TOA) software package.\n')
            file_id.write( '\n')
            file_id.write( '# RECORD FORMAT: "species_name";"plaza_id"\n')
            file_id.write( '\n')
            file_id.write( '"Actinidia chinensis";"ach"\n')
            file_id.write( '"Amborella trichopoda";"atr"\n')
            file_id.write( '"Ananas comosus";"aco"\n')
            file_id.write( '"Arabidopsis lyrata";"aly"\n')
            file_id.write( '"Arabidopsis thaliana";"ath"\n')
            file_id.write( '"Arachis ipaensis";"aip"\n')
            file_id.write( '"Beta vulgaris";"bvu"\n')
            file_id.write( '"Brachypodium distachyon";"bdi"\n')
            file_id.write( '"Brassica oleracea";"bol"\n')
            file_id.write( '"Brassica rapa";"bra"\n')
            file_id.write( '"Cajanus cajan";"ccaj"\n')
            file_id.write( '"Capsella rubella";"cru"\n')
            file_id.write( '"Capsicum annuum";"can"\n')
            file_id.write( '"Carica papaya";"cpa"\n')
            file_id.write( '"Chenopodium quinoa";"cqu"\n')
            file_id.write( '"Chlamydomonas reinhardtii";"cre"\n')
            file_id.write( '"Cicer arietinum";"car"\n')
            file_id.write( '"Citrullus lanatus";"cla"\n')
            file_id.write( '"Citrus clementina";"ccl"\n')
            file_id.write( '"Coffea canephora";"ccan"\n')
            file_id.write( '"Corchorus olitorius";"col"\n')
            file_id.write( '"Cucumis melo";"cme"\n')
            file_id.write( '"Cucumis sativus L.";"csa"\n')
            file_id.write( '"Cycas micholitzii";"cmi"\n')
            file_id.write( '"Daucus carota";"dca"\n')
            file_id.write( '"Elaeis guineensis";"egu"\n')
            file_id.write( '"Erythranthe guttata";"egut"\n')
            file_id.write( '"Eucalyptus grandis";"egr"\n')
            file_id.write( '"Fragaria vesca";"fve"\n')
            file_id.write( '"Ginkgo biloba";"gbi"\n')
            file_id.write( '"Glycine max";"gma"\n')
            file_id.write( '"Gnetum montanum";"gmo"\n')
            file_id.write( '"Gossypium raimondii";"gra"\n')
            file_id.write( '"Hevea brasiliensis";"hbr"\n')
            file_id.write( '"Hordeum vulgare";"hvu"\n')
            file_id.write( '"Malus domestica";"mdo"\n')
            file_id.write( '"Manihot esculenta";"mes"\n')
            file_id.write( '"Marchantia polymorpha";"mpo"\n')
            file_id.write( '"Medicago truncatula";"mtr"\n')
            file_id.write( '"Micromonas commoda";"mco"\n')
            file_id.write( '"Musa acuminata";"mac"\n')
            file_id.write( '"Nelumbo nucifera";"nnu"\n')
            file_id.write( '"Oropetium thomaeum";"oth"\n')
            file_id.write( '"Oryza brachyantha";"obr"\n')
            file_id.write( '"Oryza sativa ssp. indica";"osaindica"\n')
            file_id.write( '"Oryza sativa ssp. japonica";"osa"\n')
            file_id.write( '"Petunia axillaris";"pax"\n')
            file_id.write( '"Phalaenopsis equestris";"peq"\n')
            file_id.write( '"Phyllostachys edulis";"ped"\n')
            file_id.write( '"Physcomitrella patens";"ppa"\n')
            file_id.write( '"Picea abies";"pab"\n')
            file_id.write( '"Picea glauca";"pgl"\n')
            file_id.write( '"Picea sitchensis";"psi"\n')
            file_id.write( '"Pinus pinaster";"ppi"\n')
            file_id.write( '"Pinus sylvestris";"psy"\n')
            file_id.write( '"Pinus taeda";"pta"\n')
            file_id.write( '"Populus trichocarpa";"ptr"\n')
            file_id.write( '"Prunus persica";"ppe"\n')
            file_id.write( '"Pseudotsuga menziesii";"pme"\n')
            file_id.write( '"Pyrus bretschneideri";"pbr"\n')
            file_id.write( '"Ricinus communis";"rco"\n')
            file_id.write( '"Schrenkiella parvula";"spa"\n')
            file_id.write( '"Selaginella moellendorffii";"smo"\n')
            file_id.write( '"Setaria italica";"sit"\n')
            file_id.write( '"Solanum lycopersicum";"sly"\n')
            file_id.write( '"Solanum tuberosum";"stu"\n')
            file_id.write( '"Sorghum bicolor";"sbi"\n')
            file_id.write( '"Spirodela polyrhiza";"spo"\n')
            file_id.write( '"Tarenaya hassleriana";"tha"\n')
            file_id.write( '"Taxus baccata";"tba"\n')
            file_id.write( '"Theobroma cacao";"tca"\n')
            file_id.write( '"Trifolium pratense";"tpr"\n')
            file_id.write( '"Triticum aestivum";"tae"\n')
            file_id.write( '"Utricularia gibba";"ugi"\n')
            file_id.write( '"Vigna radiata var. radiata";"vra"\n')
            file_id.write( '"Vitis vinifera";"vvi"\n')
            file_id.write( '"Zea mays";"zma"\n')
            file_id.write( '"Ziziphus jujuba";"zju"\n')
            file_id.write( '"Zostera marina";"zosmarina"\n')
            file_id.write( '"Zoysia japonica ssp. nagirizaki";"zjn"\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_species_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_species_file(strict):
    '''
    Check the file of species.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the pattern of the data records
    # format: "species_name";"plaza_id"
    record_pattern = re.compile(r'^"(.*)";"(.*)"$')

    # open the file of species
    try:
        species_file_id = open(get_species_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_species_file()} can not be opened.')
        OK = False

    # check that all records are OK
    if OK:

        # read the first record
        record = species_file_id.readline()

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try:
                    mo = record_pattern.match(record)
                    # -- species_name = mo.group(1).strip()
                    # -- plaza_species_id = mo.group(2).strip()
                except Exception as e:
                    record_text = record.replace("\n", "")
                    error_list.append(f'*** ERROR: There is a format error in the record: {record_text}.')
                    OK = False
                    break

            # read the next record
            record = species_file_id.readline()

    # close the file of species
    species_file_id.close()

    # warn that the file of species is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe file {get_species_file()} is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_species_file():
    '''
    Get the species file path.
    '''

    # assign the species file path
    species_file = f'{xlib.get_config_dir()}/species.txt'

    # return the species file path
    return species_file

#-------------------------------------------------------------------------------

def manage_toa_database(cluster_name, process_type, log, function=None):
    '''
   Manage processes of the TOA database.
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

    # check the TOA config file
    if OK:
        if not os.path.isfile(get_toa_config_file()):
            log.write('*** ERROR: The TOA config file does not exist. Please, recreate it.\n')
            OK = False

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the TOA is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_toa_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_toa_name()} is not installed.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        if process_type == xlib.get_toa_type_recreate(): 
            current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_recreate_toa_database_code())
        elif process_type == xlib.get_toa_type_rebuild(): 
            current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_rebuild_toa_database_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if process_type == xlib.get_toa_type_recreate():
            script = get_recreate_toa_database_script()
            log.write(f'Building the process script {script} ...\n')
            (OK, error_list) = build_recreate_toa_database_script(cluster_name, current_run_dir)
        elif process_type == xlib.get_toa_type_rebuild(): 
            script = get_rebuild_toa_database_script()
            log.write(f'Building the process script {script} ...\n')
            (OK, error_list) = build_rebuild_toa_database_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {script} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(script)}'
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(script)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(script)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if process_type == xlib.get_toa_type_recreate():
            starter = get_recreate_toa_database_starter()
            log.write(f'Building the process starter {starter} ...\n')
            (OK, error_list) = build_recreate_toa_database_starter(current_run_dir)
        elif process_type == xlib.get_toa_type_rebuild():
            starter = get_rebuild_toa_database_starter()
            log.write(f'Building the process starter {starter} ...\n')
            (OK, error_list) = build_rebuild_toa_database_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_recreate_toa_database_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_recreate_toa_database_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_recreate_toa_database_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_recreate_toa_database_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_recreate_toa_database_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_recreate_toa_database_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_recreate_toa_database_starter()), log)

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

def build_recreate_toa_database_script(cluster_name, current_run_dir):
    '''
    Build the script to recreate a TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_recreate_toa_database_script())):
            os.makedirs(os.path.dirname(get_recreate_toa_database_script()))
        with open(get_recreate_toa_database_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function create_toa_database_dir\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Creating the database directory ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'        mkdir --parents {toa_config_dict["TOA_DB_DIR"]}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mkdir $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function recreate_toa_database\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        recreate-database.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error recreate-database.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_get_toa_process_recreate_toa_database_name()}'
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
                script_file_id.write( 'create_toa_database_dir\n')
                script_file_id.write( 'recreate_toa_database\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_recreate_toa_database_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_recreate_toa_database_starter(current_run_dir):
    '''
    Build the starter of the script to recreate a TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_recreate_toa_database_starter())):
            os.makedirs(os.path.dirname(get_recreate_toa_database_starter()))
        with open(get_recreate_toa_database_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_recreate_toa_database_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_recreate_toa_database_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_recreate_toa_database_script():
    '''
    Get the script path in the local computer to recreate a TOA database.
    '''

    # assign the script path
    recreate_toa_database_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_recreate_toa_database_code()}-process.sh'

    # return the script path
    return recreate_toa_database_script

#-------------------------------------------------------------------------------

def get_recreate_toa_database_starter():
    '''
    Get the starter path in the local computer to recreate a TOA database.
    '''

    # assign the starter path
    recreate_toa_database_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_recreate_toa_database_code()}-process-starter.sh'

    # return the starter path
    return recreate_toa_database_starter

#-------------------------------------------------------------------------------

def build_rebuild_toa_database_script(cluster_name, current_run_dir):
    '''
    Build the script to rebuild a TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_rebuild_toa_database_script())):
            os.makedirs(os.path.dirname(get_rebuild_toa_database_script()))
        with open(get_rebuild_toa_database_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function rebuild_toa_database\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        rebuild-database.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rebuild-database.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_get_toa_process_rebuild_toa_database_name()}'
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
                script_file_id.write( 'rebuild_toa_database\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_rebuild_toa_database_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_rebuild_toa_database_starter(current_run_dir):
    '''
    Build the starter of the script to rebuild a TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_rebuild_toa_database_starter())):
            os.makedirs(os.path.dirname(get_rebuild_toa_database_starter()))
        with open(get_rebuild_toa_database_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{current_run_dir}/{os.path.basename(get_rebuild_toa_database_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_rebuild_toa_database_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_rebuild_toa_database_script():
    '''
    Get the script path in the local computer to rebuild a TOA database.
    '''

    # assign the script path
    rebuild_toa_database_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_rebuild_toa_database_code()}-process.sh'

    # return the script path
    return rebuild_toa_database_script

#-------------------------------------------------------------------------------

def get_rebuild_toa_database_starter():
    '''
    Get the starter path in the local computer to rebuild a TOA database.
    '''

    # assign the starter path
    rebuild_toa_database_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_rebuild_toa_database_code()}-process-starter.sh'

    # return the starter path
    return rebuild_toa_database_starter

#-------------------------------------------------------------------------------

def manage_genomic_database(cluster_name, process_type, genomic_database, log, function=None):
    '''
    Manage processes of genomic database.
    '''

    # initialize the control variable
    OK = True

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # get the data directory
    data_dir = toa_config_dict['DATA_DIR']

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

    # check the TOA config file
    if OK:
        if not os.path.isfile(get_toa_config_file()):
            log.write('*** ERROR: The TOA config file does not exist. Please, recreate it.\n')
            OK = False

    # check the genomic dataset and species file
    if OK:
        if process_type == xlib.get_toa_type_load_data() and genomic_database == xlib.get_toa_data_basic_data_code():
            if  not os.path.isfile(get_dataset_file()):
                log.write('*** ERROR: The genomic dataset file does not exist. Please, recreate it.\n')
                OK = False
            if  not os.path.isfile(get_species_file()):
                log.write('*** ERROR: The species file does not exist. Please, recreate it.\n')
                OK = False

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name)
        if master_state_code != 16:
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
            OK = False

    # check the TOA is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_toa_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_toa_name()} is not installed.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')

        # processes to build proteomes
        if process_type == xlib.get_toa_type_build_proteome():
            if genomic_database == xlib.get_toa_data_gymno_01_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_proteome_gymno_01_code())
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_proteome_dicots_04_code())
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_proteome_monocots_04_code())
            elif genomic_database == xlib.get_toa_data_refseq_plant_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_proteome_refseq_plant_code())

        # processes to download functional annotations from a genomic database server
        elif process_type == xlib.get_toa_type_download_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_basic_data_code())
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_gymno_01_code())
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_dicots_04_code())
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_monocots_04_code())
            elif genomic_database == xlib.get_toa_data_taxonomy_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_taxonomy_code())
            elif genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code())
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_protein_gi_code())
            elif genomic_database == xlib.get_toa_data_gene_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_gene_code())
            elif genomic_database == xlib.get_toa_data_interpro_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_interpro_code())
            elif genomic_database == xlib.get_toa_data_go_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_download_go_code())

        # processes to load data of a genomic database into TOA database
        elif process_type == xlib.get_toa_type_load_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_basic_data_code())
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_gymno_01_code())
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_dicots_04_code())
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_monocots_04_code())
            elif genomic_database == xlib.get_toa_data_gene_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_gene_code())
            elif genomic_database == xlib.get_toa_data_interpro_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_interpro_code())
            elif genomic_database == xlib.get_toa_data_go_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_load_go_code())

        # processes to build BLAST databases for BLAST+
        elif process_type == xlib.get_toa_type_build_blastplus_db():
            if genomic_database == xlib.get_toa_data_nt_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_nt_blastplus_db_code())
            elif genomic_database == xlib.get_toa_data_nr_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_nr_blastplus_db_code())

        # processes to build BLAST databases for DIAMOND
        elif process_type == xlib.get_toa_type_build_diamond_db():
            if genomic_database == xlib.get_toa_data_nr_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_nr_diamond_db_code())

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code())
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_protein_gi_code())

        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # create the data subdirectory from the database directory
    if OK:
        if process_type == xlib.get_toa_type_load_data() and genomic_database == xlib.get_toa_data_basic_data_code():
            log.write(f'{xlib.get_separator()}\n')
            log.write('Creating the TOA data directory ...\n')
            command = f'mkdir --parents {data_dir}'
            (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write(f'The directory path {data_dir} is created.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the file of datasets to the cluster
    if OK:
        if process_type == xlib.get_toa_type_load_data() and genomic_database == xlib.get_toa_data_basic_data_code():
            if  not os.path.isfile(get_dataset_file()):
                log.write('*** ERROR: The genomic dataset file does not exist. Please, recreate it.\n')
                OK = False
            else:
                log.write(f'{xlib.get_separator()}\n')
                log.write(f'Uploading the file {get_dataset_file()} to the directory {data_dir} ...\n')
                cluster_path = f'{data_dir}/{os.path.basename(get_dataset_file())}'
                (OK, error_list) = xssh.put_file(sftp_client, get_dataset_file(), cluster_path)
                if OK:
                    log.write('The file is uploaded.\n')
                else:
                    for error in error_list:
                        log.write(f'{error}\n')

    # upload the file of species to the cluster
    if OK:
        if process_type == xlib.get_toa_type_load_data() and genomic_database == xlib.get_toa_data_basic_data_code():
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the file {get_species_file()} to the directory {data_dir} ...\n')
            cluster_path = f'{data_dir}/{os.path.basename(get_species_file())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_species_file(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')

    # build the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')

        # processes to build proteomes
        if process_type == xlib.get_toa_type_build_proteome():
            if genomic_database == xlib.get_toa_data_gymno_01_code():
                script = get_gymno_01_proteome_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_gymno_01_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_proteome_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_dicots_04_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_proteome_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_monocots_04_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_refseq_plant_code():
                script = get_refseq_plant_proteome_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_refseq_plant_proteome_script(cluster_name, current_run_dir)

        # processes to download functional annotations from a genomic database server
        elif process_type == xlib.get_toa_type_download_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                script = get_basic_data_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_basic_data_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                script = get_gymno_01_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_gymno_01_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_dicots_04_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_monocots_04_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_taxonomy_code():
                script = get_taxonomy_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_taxonomy_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                script = get_gene_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_gene_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                script = get_interpro_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_interpro_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                script = get_go_download_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_go_download_script(cluster_name, current_run_dir)

        # processes to load data of a genomic database into TOA database
        elif process_type == xlib.get_toa_type_load_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                script = get_basic_data_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_basic_data_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                script = get_gymno_01_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_gymno_01_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_dicots_04_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_monocots_04_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                script = get_gene_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_gene_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                script = get_interpro_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_interpro_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                script = get_go_load_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_go_load_script(cluster_name, current_run_dir)

        # processes to build BLAST databases for BLAST+
        elif process_type == xlib.get_toa_type_build_blastplus_db():
            if genomic_database == xlib.get_toa_data_nt_code():
                script = get_nt_blastplus_db_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_nt_blastplus_db_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_nr_code():
                script = get_nr_blastplus_db_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_nr_blastplus_db_script(cluster_name, current_run_dir)

        # processes to build BLAST databases for DIAMOND
        elif process_type == xlib.get_toa_type_build_diamond_db():
            if genomic_database == xlib.get_toa_data_nr_code():
                script = get_nr_diamond_db_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_nr_diamond_db_script(cluster_name, current_run_dir)

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                script = get_viridiplantae_nucleotide_gi_gilist_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_viridiplantae_nucleotide_gi_gilist_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                script = get_viridiplantae_protein_gi_gilist_script()
                log.write(f'Building the process script {script} ...\n')
                (OK, error_list) = build_viridiplantae_protein_gi_gilist_script(cluster_name, current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {script} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(script)}'
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(script)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(script)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')

        # processes to build proteomes
        if process_type == xlib.get_toa_type_build_proteome():
            if genomic_database == xlib.get_toa_data_gymno_01_code():
                starter = get_gymno_01_proteome_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_gymno_01_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_proteome_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_dicots_04_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_proteome_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_monocots_04_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_refseq_plant_code():
                starter = get_refseq_plant_proteome_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_refseq_plant_proteome_starter(current_run_dir)

        # processes to download functional annotations from a genomic database server
        elif process_type == xlib.get_toa_type_download_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                starter = get_basic_data_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_basic_data_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                starter = get_gymno_01_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_gymno_01_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_dicots_04_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_monocots_04_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_taxonomy_code():
                starter = get_taxonomy_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_taxonomy_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                starter = get_gene_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_gene_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                starter = get_interpro_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_interpro_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                starter = get_go_download_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_go_download_starter(current_run_dir)

        # processes to load data of a genomic database into TOA database
        elif process_type == xlib.get_toa_type_load_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                starter = get_basic_data_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_basic_data_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                starter = get_gymno_01_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_gymno_01_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_dicots_04_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_monocots_04_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                starter = get_gene_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_gene_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                starter = get_interpro_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_interpro_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                starter = get_go_load_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_go_load_starter(current_run_dir)

        # processes to build BLAST databases for BLAST+
        elif process_type == xlib.get_toa_type_build_blastplus_db():
            if genomic_database == xlib.get_toa_data_nt_code():
                starter = get_nt_blastplus_db_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_nt_blastplus_db_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_nr_code():
                starter = get_nr_blastplus_db_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_nr_blastplus_db_starter(current_run_dir)

        # processes to build BLAST databases for DIAMOND
        elif process_type == xlib.get_toa_type_build_diamond_db():
            if genomic_database == xlib.get_toa_data_nr_code():
                starter = get_nr_diamond_db_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_nr_diamond_db_starter(current_run_dir)

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                starter = get_viridiplantae_nucleotide_gi_gilist_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_viridiplantae_nucleotide_gi_gilist_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                starter = get_viridiplantae_protein_gi_gilist_starter()
                log.write(f'Building the process starter {starter} ...\n')
                (OK, error_list) = build_viridiplantae_protein_gi_gilist_starter(current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {starter} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(starter)}'
        (OK, error_list) = xssh.put_file(sftp_client, starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(starter)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(starter)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(starter)} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(starter), log)

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

def build_basic_data_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download other basic data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_basic_data_download_script())):
            os.makedirs(os.path.dirname(get_basic_data_download_script()))
        with open(get_basic_data_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$EC_DIR" ]; then mkdir --parents $EC_DIR; fi\n')
                script_file_id.write( 'if [ ! -d "$KEGG_DIR" ]; then mkdir --parents $KEGG_DIR; fi\n')
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
                script_file_id.write( 'function download_basic_data\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Enzyme Commission (EC) ids ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document  $EC_IDS_FILE \\\n')
                script_file_id.write( '            $EC_IDS_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading KEGG ids ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $KEGG_IDS_FILE \\\n')
                script_file_id.write( '            $KEGG_IDS_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_basic_data_name()}'
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
                script_file_id.write( 'download_basic_data\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_basic_data_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_basic_data_download_starter(current_run_dir):
    '''
    Build the starter of the script to download other basic data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_basic_data_download_starter())):
            os.makedirs(os.path.dirname(get_basic_data_download_starter()))
        with open(get_basic_data_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_basic_data_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_basic_data_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_basic_data_download_script():
    '''
    Get the script path in the local computer to download other basic data.
    '''

    # assign the script path
    basic_data_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_basic_data_code()}-process.sh'

    # return the script path
    return basic_data_download_script

#-------------------------------------------------------------------------------

def get_basic_data_download_starter():
    '''
    Get the script path in the local computer to download other basic data.
    '''

    # assign the starter path
    basic_data_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_basic_data_code()}-process-starter.sh'

    # return the starter path
    return basic_data_download_starter

#-------------------------------------------------------------------------------

def build_basic_data_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load basic data into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_basic_data_load_script())):
            os.makedirs(os.path.dirname(get_basic_data_load_script()))
        with open(get_basic_data_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$TOA_DB_DIR" ]; then mkdir --parents $TOA_DB_DIR; fi\n')
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
                script_file_id.write( 'function load_basic_data\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading basic data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-basic-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --datasets=$DATASET_FILE \\\n')
                script_file_id.write( '            --species=$SPECIES_FILE \\\n')
                script_file_id.write( '            --ecids=$EC_IDS_FILE \\\n')
                script_file_id.write( '            --keggids=$KEGG_IDS_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_basic_data_name()}'
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
                script_file_id.write( 'load_basic_data\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_basic_data_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_basic_data_load_starter(current_run_dir):
    '''
    Build the starter of the script to load basic data into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_basic_data_load_starter())):
            os.makedirs(os.path.dirname(get_basic_data_load_starter()))
        with open(get_basic_data_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_basic_data_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_basic_data_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_basic_data_load_script():
    '''
    Get the script path in the local computer to load basic data into TOA database.
    '''

    # assign the script path
    basic_data_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_basic_data_code()}-process.sh'

    # return the script path
    return basic_data_load_script

#-------------------------------------------------------------------------------

def get_basic_data_load_starter():
    '''
    Get the starter path in the local computer to load basic data into TOA database.
    '''

    # assign the starter path
    basic_data_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_basic_data_code()}-process-starter.sh'

    # return the starter path
    return basic_data_load_starter

#-------------------------------------------------------------------------------

def build_gymno_01_proteome_script(cluster_name, current_run_dir):
    '''
    Build the script to build the Gymno PLAZA 1.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_proteome_script())):
            os.makedirs(os.path.dirname(get_gymno_01_proteome_script()))
        with open(get_gymno_01_proteome_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$GYMNO_01_BLASTPLUS_DB_DIR" ]; then rm -rf $GYMNO_01_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $GYMNO_01_BLASTPLUS_DB_DIR\n')
                script_file_id.write( 'if [ -d "$GYMNO_01_DIAMOND_DB_DIR" ]; then rm -rf $GYMNO_01_DIAMOND_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $GYMNO_01_DIAMOND_DB_DIR\n')
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
                script_file_id.write( 'function build_gymno01_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading proteome file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GYMNO_01_PROTEOME_FILE \\\n')
                script_file_id.write( '            $GYMNO_01_PROTEOME_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating BLAST+ database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        gunzip -c $GYMNO_01_PROTEOME_FILE | \\\n')
                script_file_id.write( '        makeblastdb \\\n')
                script_file_id.write( '            -title $GYMNO_01_BLASTPLUS_DB_NAME \\\n')
                script_file_id.write( '            -dbtype prot \\\n')
                script_file_id.write( '            -input_type fasta \\\n')
                script_file_id.write( '            -hash_index \\\n')
                script_file_id.write( '            -in - \\\n')
                script_file_id.write( '            -out $GYMNO_01_BLASTPLUS_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error makeblastdb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "BLAST+ database is generated."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating DIAMOND database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_diamond_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        diamond makedb \\\n')
                script_file_id.write( '            --threads 4 \\\n')
                script_file_id.write( '            --in $GYMNO_01_PROTEOME_FILE \\\n')
                script_file_id.write( '            --db $GYMNO_01_DIAMOND_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error diamond-makedb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "DIAMOND database is generated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_proteome_gymno_01_name()}'
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
                script_file_id.write( 'build_gymno01_proteome\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_proteome_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gymno_01_proteome_starter(current_run_dir):
    '''
    Build the starter of script to build the Gymno PLAZA 1.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_proteome_starter())):
            os.makedirs(os.path.dirname(get_gymno_01_proteome_starter()))
        with open(get_gymno_01_proteome_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gymno_01_proteome_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_proteome_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_proteome_script():
    '''
    Get the script path in the local computer to build the Gymno PLAZA 1.0 proteome.
    '''

    # assign the script path
    gymno_01_proteome_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_gymno_01_code()}-process.sh'

    # return the script path
    return gymno_01_proteome_script

#-------------------------------------------------------------------------------

def get_gymno_01_proteome_starter():
    '''
    Get the starter path in the local computer to build the Gymno PLAZA 1.0 proteome.
    '''

    # assign the starter path
    gymno_01_proteome_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_gymno_01_code()}-process-starter.sh'

    # return the starter path
    return gymno_01_proteome_starter

#-------------------------------------------------------------------------------

def build_gymno_01_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_download_script())):
            os.makedirs(os.path.dirname(get_gymno_01_download_script()))
        with open(get_gymno_01_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$GYMNO_01_GENEDESC_DIR" ]; then rm -rf $GYMNO_01_GENEDESC_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $GYMNO_01_GENEDESC_DIR\n')
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
                script_file_id.write( 'function download_gymno01_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading gene description files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=0 \\\n')
                script_file_id.write( '            --no-host-directories \\\n')
                script_file_id.write( '            --cut-dirs=4 \\\n')
                script_file_id.write( '            --accept=$GYMNO_01_GENEDESC_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$GYMNO_01_GENEDESC_DIR \\\n')
                script_file_id.write( '            $GYMNO_01_GENEDESC_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading InterPro file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GYMNO_01_INTERPRO_FILE \\\n')
                script_file_id.write( '            $GYMNO_01_INTERPRO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GYMNO_01_GO_FILE \\\n')
                script_file_id.write( '            $GYMNO_01_GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GYMNO_01_MAPMAN_FILE \\\n')
                script_file_id.write( '            $GYMNO_01_MAPMAN_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_gymno_01_name()}'
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
                script_file_id.write( 'download_gymno01_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gymno_01_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_download_starter())):
            os.makedirs(os.path.dirname(get_gymno_01_download_starter()))
        with open(get_gymno_01_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gymno_01_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_download_script():
    '''
    Get the script path in the local computer to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # assign the script path
    gymno_01_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_gymno_01_code()}-process.sh'

    # return the script path
    return gymno_01_download_script

#-------------------------------------------------------------------------------

def get_gymno_01_download_starter():
    '''
    Get the script path in the local computer to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # assign the starter path
    gymno_01_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_gymno_01_code()}-process-starter.sh'

    # return the starter path
    return gymno_01_download_starter

#-------------------------------------------------------------------------------

def build_gymno_01_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_load_script())):
            os.makedirs(os.path.dirname(get_gymno_01_load_script()))
        with open(get_gymno_01_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_gymno01_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-plaza-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --dataset=gymno_01 \\\n')
                script_file_id.write( '            --species=all \\\n')
                script_file_id.write( '            --genedesc=$GYMNO_01_GENEDESC_DIR \\\n')
                script_file_id.write( '            --interpro=$GYMNO_01_INTERPRO_FILE \\\n')
                script_file_id.write( '            --go=$GYMNO_01_GO_FILE \\\n')
                script_file_id.write( '            --mapman=$GYMNO_01_MAPMAN_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_gymno_01_name()}'
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
                script_file_id.write( 'load_gymno01_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gymno_01_load_starter(current_run_dir):
    '''
    Build the starter of the script to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_gymno_01_load_starter())):
            os.makedirs(os.path.dirname(get_gymno_01_load_starter()))
        with open(get_gymno_01_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gymno_01_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gymno_01_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_load_script():
    '''
    Get the script path in the local computer to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # assign the script path
    gymno_01_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_gymno_01_code()}-process.sh'

    # return the script path
    return gymno_01_load_script

#-------------------------------------------------------------------------------

def get_gymno_01_load_starter():
    '''
    Get the starter path in the local computer to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # assign the starter path
    gymno_01_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_gymno_01_code()}-process-starter.sh'

    # return the starter path
    return gymno_01_load_starter

#-------------------------------------------------------------------------------

def build_dicots_04_proteome_script(cluster_name, current_run_dir):
    '''
    Build the script to build the Dicots PLAZA 4.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_proteome_script())):
            os.makedirs(os.path.dirname(get_dicots_04_proteome_script()))
        with open(get_dicots_04_proteome_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$DICOTS_04_BLASTPLUS_DB_DIR" ]; then rm -rf $DICOTS_04_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $DICOTS_04_BLASTPLUS_DB_DIR\n')
                script_file_id.write( 'if [ -d "$DICOTS_04_DIAMOND_DB_DIR" ]; then rm -rf $DICOTS_04_DIAMOND_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $DICOTS_04_DIAMOND_DB_DIR\n')
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
                script_file_id.write( 'function build_dicots04_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading proteome file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $DICOTS_04_PROTEOME_FILE \\\n')
                script_file_id.write( '            $DICOTS_04_PROTEOME_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating BLAST+ database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        gunzip -c $DICOTS_04_PROTEOME_FILE | \\\n')
                script_file_id.write( '        makeblastdb \\\n')
                script_file_id.write( '            -title $DICOTS_04_BLASTPLUS_DB_NAME \\\n')
                script_file_id.write( '            -dbtype prot \\\n')
                script_file_id.write( '            -input_type fasta \\\n')
                script_file_id.write( '            -hash_index \\\n')
                script_file_id.write( '            -in - \\\n')
                script_file_id.write( '            -out $DICOTS_04_BLASTPLUS_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error makeblastdb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "BLAST+ database is generated."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating DIAMOND database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_diamond_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        diamond makedb \\\n')
                script_file_id.write( '            --threads 4 \\\n')
                script_file_id.write( '            --in $DICOTS_04_PROTEOME_FILE \\\n')
                script_file_id.write( '            --db $DICOTS_04_DIAMOND_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error diamond-makedb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "DIAMOND database is generated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_proteome_dicots_04_name()}'
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
                script_file_id.write( 'build_dicots04_proteome\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_proteome_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_dicots_04_proteome_starter(current_run_dir):
    '''
    Build the starter of script to build the Dicots PLAZA 4.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_proteome_starter())):
            os.makedirs(os.path.dirname(get_dicots_04_proteome_starter()))
        with open(get_dicots_04_proteome_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_dicots_04_proteome_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_proteome_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_proteome_script():
    '''
    Get the script path in the local computer to build the Dicots PLAZA 4.0 proteome.
    '''

    # assign the script path
    dicots_04_proteome_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_dicots_04_code()}-process.sh'

    # return the script path
    return dicots_04_proteome_script

#-------------------------------------------------------------------------------

def get_dicots_04_proteome_starter():
    '''
    Get the starter path in the local computer to build the Dicots PLAZA 4.0 proteome.
    '''

    # assign the starter path
    dicots_04_proteome_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_dicots_04_code()}-process-starter.sh'

    # return the starter path
    return dicots_04_proteome_starter

#-------------------------------------------------------------------------------

def build_dicots_04_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_download_script())):
            os.makedirs(os.path.dirname(get_dicots_04_download_script()))
        with open(get_dicots_04_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$DICOTS_04_GENEDESC_DIR" ]; then rm -rf $DICOTS_04_GENEDESC_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $DICOTS_04_GENEDESC_DIR\n')
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
                script_file_id.write( 'function download_dicots04_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading gene description files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=0 \\\n')
                script_file_id.write( '            --no-host-directories \\\n')
                script_file_id.write( '            --cut-dirs=4 \\\n')
                script_file_id.write( '            --accept=$DICOTS_04_GENEDESC_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$DICOTS_04_GENEDESC_DIR \\\n')
                script_file_id.write( '            $DICOTS_04_GENEDESC_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading InterPro file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $DICOTS_04_INTERPRO_FILE \\\n')
                script_file_id.write( '            $DICOTS_04_INTERPRO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $DICOTS_04_GO_FILE \\\n')
                script_file_id.write( '            $DICOTS_04_GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $DICOTS_04_MAPMAN_FILE \\\n')
                script_file_id.write( '            $DICOTS_04_MAPMAN_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_dicots_04_name()}'
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
                script_file_id.write( 'download_dicots04_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_dicots_04_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_download_starter())):
            os.makedirs(os.path.dirname(get_dicots_04_download_starter()))
        with open(get_dicots_04_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_dicots_04_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_download_script():
    '''
    Get the script path in the local computer to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # assign the script path
    dicots_04_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_dicots_04_code()}-process.sh'

    # return the script path
    return dicots_04_download_script

#-------------------------------------------------------------------------------

def get_dicots_04_download_starter():
    '''
    Get the script path in the local computer to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # assign the starter path
    dicots_04_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_dicots_04_code()}-process-starter.sh'

    # return the starter path
    return dicots_04_download_starter

#-------------------------------------------------------------------------------

def build_dicots_04_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_load_script())):
            os.makedirs(os.path.dirname(get_dicots_04_load_script()))
        with open(get_dicots_04_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_dicots04_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-plaza-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --dataset=dicots_04 \\\n')
                script_file_id.write( '            --species=all \\\n')
                script_file_id.write( '            --genedesc=$DICOTS_04_GENEDESC_DIR \\\n')
                script_file_id.write( '            --interpro=$DICOTS_04_INTERPRO_FILE \\\n')
                script_file_id.write( '            --go=$DICOTS_04_GO_FILE \\\n')
                script_file_id.write( '            --mapman=$DICOTS_04_MAPMAN_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_dicots_04_name()}'
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
                script_file_id.write( 'load_dicots04_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_dicots_04_load_starter(current_run_dir):
    '''
    Build the starter of the script to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_dicots_04_load_starter())):
            os.makedirs(os.path.dirname(get_dicots_04_load_starter()))
        with open(get_dicots_04_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_dicots_04_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_dicots_04_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_load_script():
    '''
    Get the script path in the local computer to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the script path
    dicots_04_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_dicots_04_code()}-process.sh'

    # return the script path
    return dicots_04_load_script

#-------------------------------------------------------------------------------

def get_dicots_04_load_starter():
    '''
    Get the starter path in the local computer to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the starter path
    dicots_04_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_dicots_04_code()}-process-starter.sh'

    # return the starter path
    return dicots_04_load_starter

#-------------------------------------------------------------------------------

def build_monocots_04_proteome_script(cluster_name, current_run_dir):
    '''
    Build the script to build the Monocots PLAZA 4.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_proteome_script())):
            os.makedirs(os.path.dirname(get_monocots_04_proteome_script()))
        with open(get_monocots_04_proteome_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$MONOCOTS_04_BLASTPLUS_DB_DIR" ]; then rm -rf $MONOCOTS_04_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $MONOCOTS_04_BLASTPLUS_DB_DIR\n')
                script_file_id.write( 'if [ -d "$MONOCOTS_04_DIAMOND_DB_DIR" ]; then rm -rf $MONOCOTS_04_DIAMOND_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $MONOCOTS_04_DIAMOND_DB_DIR\n')
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
                script_file_id.write( 'function build_monocots04_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading proteome file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $MONOCOTS_04_PROTEOME_FILE \\\n')
                script_file_id.write( '            $MONOCOTS_04_PROTEOME_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating BLAST+ database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        gunzip -c $MONOCOTS_04_PROTEOME_FILE | \\\n')
                script_file_id.write( '        makeblastdb \\\n')
                script_file_id.write( '            -title $MONOCOTS_04_BLASTPLUS_DB_NAME \\\n')
                script_file_id.write( '            -dbtype prot \\\n')
                script_file_id.write( '            -input_type fasta \\\n')
                script_file_id.write( '            -hash_index \\\n')
                script_file_id.write( '            -in - \\\n')
                script_file_id.write( '            -out $MONOCOTS_04_BLASTPLUS_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error makeblastdb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "BLAST+ database is generated."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating DIAMOND database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_diamond_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        diamond makedb \\\n')
                script_file_id.write( '            --threads 4 \\\n')
                script_file_id.write( '            --in $MONOCOTS_04_PROTEOME_FILE \\\n')
                script_file_id.write( '            --db $MONOCOTS_04_DIAMOND_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error diamond-makedb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "DIAMOND database is generated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_proteome_monocots_04_name()}'
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
                script_file_id.write( 'build_monocots04_proteome\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_proteome_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_monocots_04_proteome_starter(current_run_dir):
    '''
    Build the starter of script to build the Monocots PLAZA 4.0 proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_proteome_starter())):
            os.makedirs(os.path.dirname(get_monocots_04_proteome_starter()))
        with open(get_monocots_04_proteome_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_monocots_04_proteome_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_proteome_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_proteome_script():
    '''
    Get the script path in the local computer to build the Monocots PLAZA 4.0 proteome.
    '''

    # assign the script path
    monocots_04_proteome_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_monocots_04_code()}-process.sh'

    # return the script path
    return monocots_04_proteome_script

#-------------------------------------------------------------------------------

def get_monocots_04_proteome_starter():
    '''
    Get the starter path in the local computer to build the Monocots PLAZA 4.0 proteome.
    '''

    # assign the starter path
    monocots_04_proteome_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_monocots_04_code()}-process-starter.sh'

    # return the starter path
    return monocots_04_proteome_starter

#-------------------------------------------------------------------------------

def build_monocots_04_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_download_script())):
            os.makedirs(os.path.dirname(get_monocots_04_download_script()))
        with open(get_monocots_04_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$MONOCOTS_04_GENEDESC_DIR" ]; then rm -rf $MONOCOTS_04_GENEDESC_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $MONOCOTS_04_GENEDESC_DIR\n')
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
                script_file_id.write( 'function download_monocots04_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading gene description files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=0 \\\n')
                script_file_id.write( '            --no-host-directories \\\n')
                script_file_id.write( '            --cut-dirs=4 \\\n')
                script_file_id.write( '            --accept=$MONOCOTS_04_GENEDESC_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$MONOCOTS_04_GENEDESC_DIR \\\n')
                script_file_id.write( '            $MONOCOTS_04_GENEDESC_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading InterPro file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $MONOCOTS_04_INTERPRO_FILE \\\n')
                script_file_id.write( '            $MONOCOTS_04_INTERPRO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $MONOCOTS_04_GO_FILE \\\n')
                script_file_id.write( '            $MONOCOTS_04_GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $MONOCOTS_04_MAPMAN_FILE \\\n')
                script_file_id.write( '            $MONOCOTS_04_MAPMAN_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_monocots_04_name()}'
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
                script_file_id.write( 'download_monocots04_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_monocots_04_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_download_starter())):
            os.makedirs(os.path.dirname(get_monocots_04_download_starter()))
        with open(get_monocots_04_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_monocots_04_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_download_script():
    '''
    Get the script path in the local computer to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # assign the script path
    monocots_04_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_monocots_04_code()}-process.sh'

    # return the script path
    return monocots_04_download_script

#-------------------------------------------------------------------------------

def get_monocots_04_download_starter():
    '''
    Get the script path in the local computer to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # assign the starter path
    monocots_04_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_monocots_04_code()}-process-starter.sh'

    # return the starter path
    return monocots_04_download_starter

#-------------------------------------------------------------------------------

def build_monocots_04_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_load_script())):
            os.makedirs(os.path.dirname(get_monocots_04_load_script()))
        with open(get_monocots_04_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_monocots04_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-plaza-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --dataset=monocots_04 \\\n')
                script_file_id.write( '            --species=all \\\n')
                script_file_id.write( '            --genedesc=$MONOCOTS_04_GENEDESC_DIR \\\n')
                script_file_id.write( '            --interpro=$MONOCOTS_04_INTERPRO_FILE \\\n')
                script_file_id.write( '            --go=$MONOCOTS_04_GO_FILE \\\n')
                script_file_id.write( '            --mapman=$MONOCOTS_04_MAPMAN_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_monocots_04_name()}'
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
                script_file_id.write( 'load_monocots04_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_monocots_04_load_starter(current_run_dir):
    '''
    Build the starter of the script to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_monocots_04_load_starter())):
            os.makedirs(os.path.dirname(get_monocots_04_load_starter()))
        with open(get_monocots_04_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_monocots_04_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_monocots_04_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_load_script():
    '''
    Get the script path in the local computer to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the script path
    monocots_04_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_monocots_04_code()}-process.sh'

    # return the script path
    return monocots_04_load_script

#-------------------------------------------------------------------------------

def get_monocots_04_load_starter():
    '''
    Get the starter path in the local computer to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the starter path
    monocots_04_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_monocots_04_code()}-process-starter.sh'

    # return the starter path
    return monocots_04_load_starter

#-------------------------------------------------------------------------------

def build_refseq_plant_proteome_script(cluster_name, current_run_dir):
    '''
    Build the script to build the NCBI RefSeq Plant proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_refseq_plant_proteome_script())):
            os.makedirs(os.path.dirname(get_refseq_plant_proteome_script()))
        with open(get_refseq_plant_proteome_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$REFSEQ_PLANT_LOCAL" ]; then mkdir --parents $REFSEQ_PLANT_LOCAL; fi\n')
                script_file_id.write( 'if [ -d "$REFSEQ_PLANT_BLASTPLUS_DB_DIR" ]; then rm -rf $REFSEQ_PLANT_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $REFSEQ_PLANT_BLASTPLUS_DB_DIR\n')
                script_file_id.write( 'if [ -d "$REFSEQ_PLANT_DIAMOND_DB_DIR" ]; then rm -rf $REFSEQ_PLANT_DIAMOND_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $REFSEQ_PLANT_DIAMOND_DB_DIR\n')
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
                script_file_id.write( 'function build_refseqplant_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading protein FASTA files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=1 \\\n')
                script_file_id.write( '            --accept=$REFSEQ_PROTEIN_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$NCBI_DIR \\\n')
                script_file_id.write( '            $REFSEQ_PLANT_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building proteome file ..."\n')
                script_file_id.write( '    > $REFSEQ_PLANT_PROTEOME_FILE\n')
                script_file_id.write( '    ls `echo $REFSEQ_PLANT_LOCAL/"*"$REFSEQ_PROTEIN_FILE_PATTERN` > $REFSEQ_PLANT_FILE_LIST\n')
                script_file_id.write( '    while read FILE_GZ; do\n')
                script_file_id.write( '        FILE_FASTA=`echo $FILE_GZ | sed "s|.gz||g"`\n')
                script_file_id.write( '        gzip --decompress --force  $FILE_GZ\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error gzip $RC; fi\n')
                script_file_id.write( '        cat $FILE_FASTA >> $REFSEQ_PLANT_PROTEOME_FILE\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error cat $RC; fi\n')
                script_file_id.write( '        rm -f $FILE_FASTA\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '    done < $REFSEQ_PLANT_FILE_LIST\n')
                script_file_id.write( '    echo "Proteome is buit."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating BLAST+ database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        makeblastdb \\\n')
                script_file_id.write( '            -title $REFSEQ_PLANT_BLASTPLUS_DB_NAME \\\n')
                script_file_id.write( '            -dbtype prot \\\n')
                script_file_id.write( '            -input_type fasta \\\n')
                # -- script_file_id.write( '            -hash_index \\\n')
                script_file_id.write( '            -in $REFSEQ_PLANT_PROTEOME_FILE \\\n')
                script_file_id.write( '            -out $REFSEQ_PLANT_BLASTPLUS_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error makeblastdb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "BLAST+ database is generated."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating DIAMOND database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_diamond_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        diamond makedb \\\n')
                script_file_id.write( '            --threads 4 \\\n')
                script_file_id.write( '            --in $REFSEQ_PLANT_PROTEOME_FILE \\\n')
                script_file_id.write( '            --db $REFSEQ_PLANT_DIAMOND_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error diamond-makedb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "DIAMOND database is generated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_proteome_refseq_plant_name()}'
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
                script_file_id.write( 'build_refseqplant_proteome\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_refseq_plant_proteome_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_refseq_plant_proteome_starter(current_run_dir):
    '''
    Build the starter of script to build the NCBI RefSeq Plant proteome.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_refseq_plant_proteome_starter())):
            os.makedirs(os.path.dirname(get_refseq_plant_proteome_starter()))
        with open(get_refseq_plant_proteome_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_refseq_plant_proteome_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_refseq_plant_proteome_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_refseq_plant_proteome_script():
    '''
    Get the script path in the local computer to build the NCBI RefSeq Plant proteome.
    '''

    # assign the script path
    refseq_plant_proteome_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_refseq_plant_code()}-process.sh'

    # return the script path
    return refseq_plant_proteome_script

#-------------------------------------------------------------------------------

def get_refseq_plant_proteome_starter():
    '''
    Get the starter path in the local computer to build the NCBI RefSeq Plant proteome.
    '''

    # assign the starter path
    refseq_plant_proteome_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_proteome_refseq_plant_code()}-process-starter.sh'

    # return the starter path
    return refseq_plant_proteome_starter

#-------------------------------------------------------------------------------

def build_taxonomy_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the NCBI Taxonomy data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_taxonomy_download_script())):
            os.makedirs(os.path.dirname(get_taxonomy_download_script()))
        with open(get_taxonomy_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi\n')
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
                script_file_id.write( 'function download_taxonomy_data\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading compressed NCBI Taxonomy database dump files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $TAXONOMY_TAXDMP_FILE \\\n')
                script_file_id.write( '            $TAXONOMY_TAXDMP_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "Decompressing NCBI Taxonomy database dump files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        unzip -o -d $NCBI_DIR $TAXONOMY_TAXDMP_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is decompressed."\n')
                script_file_id.write( '    echo "Downloading NCBI TaxID mapping for live protein sequence record. ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $TAXONOMY_PROTACCESSION_2_TAXID_FILE \\\n')
                script_file_id.write( '            $TAXONOMY_PROTACCESSION_2_TAXID_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "Downloading NCBI taxonomy identifications of Viridiplantae ..."\n')
                script_file_id.write(f'    cp $MINICONDA3_ENVS_DIR/{xlib.get_blastplus_anaconda_code()}/bin/get_species_taxids.sh .\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error cp $RC; fi\n')
                script_file_id.write(f'    sed -i "s|export PATH=|#export PATH=|g" ./get_species_taxids.sh\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
                script_file_id.write(f'    source activate {xlib.get_entrez_direct_anaconda_code()}\n')
                script_file_id.write( '    ./get_species_taxids.sh -n Viridiplantae\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error get_species_taxids.sh $RC; fi\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        get_species_taxids.sh -t 33090 >$VIRIDIPLANTAE_TAXID_LIST_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error get_species_taxids.sh $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "Taxids ared downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_proteome_refseq_plant_name()}'
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
                script_file_id.write( 'download_taxonomy_data\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_taxonomy_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_taxonomy_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the NCBI Taxonomy data.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_taxonomy_download_starter())):
            os.makedirs(os.path.dirname(get_taxonomy_download_starter()))
        with open(get_taxonomy_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_taxonomy_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_taxonomy_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_taxonomy_download_script():
    '''
    Get the script path to download the NCBI Taxonomy data.
    '''

    # assign the script path
    taxonomy_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_taxonomy_code()}-process.sh'

    # return the script path
    return taxonomy_download_script

#-------------------------------------------------------------------------------

def get_taxonomy_download_starter():
    '''
    Get the script path to download the NCBI Taxonomy data.
    '''

    # assign the starter path
    taxonomy_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_taxonomy_code()}-process-starter.sh'

    # return the starter path
    return taxonomy_download_starter

#-------------------------------------------------------------------------------

def build_nt_blastplus_db_script(cluster_name, current_run_dir):
    '''
    Build the script to build BLAST database NT for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nt_blastplus_db_script())):
            os.makedirs(os.path.dirname(get_nt_blastplus_db_script()))
        with open(get_nt_blastplus_db_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$NT_BLASTPLUS_DB_DIR" ]; then rm -rf $NT_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $NT_BLASTPLUS_DB_DIR\n')
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
                script_file_id.write( 'function build_database_nt\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading nt database files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=1 \\\n')
                script_file_id.write( '            --accept=$NT_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$NCBI_DIR \\\n')
                script_file_id.write( '            $BLAST_DATABASES_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Decompressing nt database files ..."\n')
                script_file_id.write( '    ls `echo $BLAST_DATABASES_LOCAL/$NT_FILE_PATTERN` > $NT_FILE_LIST\n')
                script_file_id.write( '    while read NT_FILE; do\n')
                script_file_id.write( '        tar --extract --gzip --file=$NT_FILE --directory=$NT_BLASTPLUS_DB_DIR\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
                script_file_id.write( '        rm -f $NT_FILE\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '    done < $NT_FILE_LIST\n')
                script_file_id.write( '    echo "Files are decompressed."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_nt_blastplus_db_name()}'
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
                script_file_id.write( 'build_database_nt\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nt_blastplus_db_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nt_blastplus_db_starter(current_run_dir):
    '''
    Build the starter of the script to build BLAST database NT for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nt_blastplus_db_starter())):
            os.makedirs(os.path.dirname(get_nt_blastplus_db_starter()))
        with open(get_nt_blastplus_db_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_nt_blastplus_db_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nt_blastplus_db_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nt_blastplus_db_script():
    '''
    Get the script path to build BLAST database NT for BLAST+.
    '''

    # assign the script path
    nt_blastplus_db_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nt_blastplus_db_code()}-process.sh'

    # return the script path
    return nt_blastplus_db_script

#-------------------------------------------------------------------------------

def get_nt_blastplus_db_starter():
    '''
    Get the starter path to build BLAST database NT for BLAST+.
    '''

    # assign the starter path
    nt_blastplus_db_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nt_blastplus_db_code()}-process-starter.sh'

    # return the starter path
    return nt_blastplus_db_starter

#-------------------------------------------------------------------------------

def build_viridiplantae_nucleotide_gi_gilist_script(cluster_name, current_run_dir):
    '''
    Build the script to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_viridiplantae_nucleotide_gi_gilist_script())):
            os.makedirs(os.path.dirname(get_viridiplantae_nucleotide_gi_gilist_script()))
        with open(get_viridiplantae_nucleotide_gi_gilist_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi\n')
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
                script_file_id.write( 'function build_viridiplantae_nucleotide_gilist\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building Viridiplantae nucleotide GI list ..."\n')
                script_file_id.write(f'    source activate {xlib.get_entrez_direct_code()}\n')
                script_file_id.write( '    esearch -db nucleotide -query "Viridiplantae[Organism]" | efetch -format uid > $NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error entrez-direct $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "List is built."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_name()}'
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
                script_file_id.write( 'build_viridiplantae_nucleotide_gilist\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_viridiplantae_nucleotide_gi_gilist_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_viridiplantae_nucleotide_gi_gilist_starter(current_run_dir):
    '''
    Build the starter of the script to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_viridiplantae_nucleotide_gi_gilist_starter())):
            os.makedirs(os.path.dirname(get_viridiplantae_nucleotide_gi_gilist_starter()))
        with open(get_viridiplantae_nucleotide_gi_gilist_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_viridiplantae_nucleotide_gi_gilist_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_viridiplantae_nucleotide_gi_gilist_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_viridiplantae_nucleotide_gi_gilist_script():
    '''
    Get the script path in the local computer to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # assign the script path
    viridiplantae_nucleotide_gi_gilist_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code()}-process.sh'

    # return the script path
    return viridiplantae_nucleotide_gi_gilist_script

#-------------------------------------------------------------------------------

def get_viridiplantae_nucleotide_gi_gilist_starter():
    '''
    Get the starter path in the local computer to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # assign the starter path
    viridiplantae_nucleotide_gi_gilist_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code()}-process-starter.sh'

    # return the starter path
    return viridiplantae_nucleotide_gi_gilist_starter

#-------------------------------------------------------------------------------

def build_nr_blastplus_db_script(cluster_name, current_run_dir):
    '''
    Build the script to build BLAST database NR for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nr_blastplus_db_script())):
            os.makedirs(os.path.dirname(get_nr_blastplus_db_script()))
        with open(get_nr_blastplus_db_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$NR_BLASTPLUS_DB_DIR" ]; then rm -rf $NR_BLASTPLUS_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $NR_BLASTPLUS_DB_DIR\n')
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
                script_file_id.write( 'function build_database_nr\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading nr database files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --recursive \\\n')
                script_file_id.write( '            --level=1 \\\n')
                script_file_id.write( '            --accept=$NR_FILE_PATTERN \\\n')
                script_file_id.write( '            --directory-prefix=$NCBI_DIR \\\n')
                script_file_id.write( '            $BLAST_DATABASES_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "Files are downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Decompressing nr database files ..."\n')
                script_file_id.write( '    ls `echo $BLAST_DATABASES_LOCAL/$NR_FILE_PATTERN` > $NR_FILE_LIST\n')
                script_file_id.write( '    while read NR_FILE; do\n')
                script_file_id.write( '        tar --extract --gzip --file=$NR_FILE --directory=$NR_BLASTPLUS_DB_DIR\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
                script_file_id.write( '        rm -f $NR_FILE\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '    done < $NR_FILE_LIST\n')
                script_file_id.write( '    echo "Files are decompressed."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_nr_blastplus_db_name()}'
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
                script_file_id.write( 'build_database_nr\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nr_blastplus_db_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nr_blastplus_db_starter(current_run_dir):
    '''
    Build the starter of the script to build BLAST database NR for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nr_blastplus_db_starter())):
            os.makedirs(os.path.dirname(get_nr_blastplus_db_starter()))
        with open(get_nr_blastplus_db_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_nr_blastplus_db_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nr_blastplus_db_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nr_blastplus_db_script():
    '''
    Get the script path to build BLAST database NR for BLAST+.
    '''

    # assign the script path
    nr_blastplus_db_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nr_blastplus_db_code()}-process.sh'

    # return the script path
    return nr_blastplus_db_script

#-------------------------------------------------------------------------------

def get_nr_blastplus_db_starter():
    '''
    Get the starter path to build BLAST database NR for BLAST+.
    '''

    # assign the starter path
    nr_blastplus_db_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nr_blastplus_db_code()}-process-starter.sh'

    # return the starter path
    return nr_blastplus_db_starter

#-------------------------------------------------------------------------------

def build_nr_diamond_db_script(cluster_name, current_run_dir):
    '''
    Build the script to build BLAST database NR for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nr_diamond_db_script())):
            os.makedirs(os.path.dirname(get_nr_diamond_db_script()))
        with open(get_nr_diamond_db_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ -d "$NR_DIAMOND_DB_DIR" ]; then rm -rf $NR_DIAMOND_DB_DIR; fi\n')
                script_file_id.write( 'mkdir --parents $NR_DIAMOND_DB_DIR\n')
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
                script_file_id.write( 'function build_database_nr\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading proteome file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $NR_PROTEOME_FILE \\\n')
                script_file_id.write( '            $NR_PROTEOME_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Generating DIAMOND database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_diamond_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        diamond makedb \\\n')
                script_file_id.write( '            --in $NR_PROTEOME_FILE \\\n')
                script_file_id.write( '            --db $NR_DIAMOND_DB_FILE\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error diamond-makedb $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "DIAMOND database is generated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_nr_diamond_db_name()}'
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
                script_file_id.write( 'build_database_nr\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nr_diamond_db_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nr_diamond_db_starter(current_run_dir):
    '''
    Build the starter of the script to build BLAST database NR for BLAST+.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nr_diamond_db_starter())):
            os.makedirs(os.path.dirname(get_nr_diamond_db_starter()))
        with open(get_nr_diamond_db_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_nr_diamond_db_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nr_diamond_db_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nr_diamond_db_script():
    '''
    Get the script path to build BLAST database NR for BLAST+.
    '''

    # assign the script path
    nr_diamond_db_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nr_diamond_db_code()}-process.sh'

    # return the script path
    return nr_diamond_db_script

#-------------------------------------------------------------------------------

def get_nr_diamond_db_starter():
    '''
    Get the starter path to build BLAST database NR for BLAST+.
    '''

    # assign the starter path
    nr_diamond_db_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_nr_diamond_db_code()}-process-starter.sh'

    # return the starter path
    return nr_diamond_db_starter

#-------------------------------------------------------------------------------

def build_viridiplantae_protein_gi_gilist_script(cluster_name, current_run_dir):
    '''
    Build the script to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_viridiplantae_protein_gi_gilist_script())):
            os.makedirs(os.path.dirname(get_viridiplantae_protein_gi_gilist_script()))
        with open(get_viridiplantae_protein_gi_gilist_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$NCBI_DIR" ]; then mkdir  --parents $NCBI_DIR; fi\n')
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
                script_file_id.write( 'function build_viridiplantae_protein_gilist\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building Viridiplantae protein GI list ..."\n')
                script_file_id.write(f'        source activate {xlib.get_entrez_direct_anaconda_code()}\n')
                script_file_id.write( '    esearch -db protein -query "Viridiplantae[Organism]" | efetch -format uid > $PROTEIN_VIRIDIPLANTAE_GI_LIST\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error entrez-direct $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "List is built."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_name()}'
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
                script_file_id.write( 'build_viridiplantae_protein_gilist\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_viridiplantae_protein_gi_gilist_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_viridiplantae_protein_gi_gilist_starter(current_run_dir):
    '''
    Build the starter of the script to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_viridiplantae_protein_gi_gilist_starter())):
            os.makedirs(os.path.dirname(get_viridiplantae_protein_gi_gilist_starter()))
        with open(get_viridiplantae_protein_gi_gilist_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_viridiplantae_protein_gi_gilist_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_viridiplantae_protein_gi_gilist_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_viridiplantae_protein_gi_gilist_script():
    '''
    Get the script path in the local computer to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # assign the script path
    viridiplantae_protein_gi_gilist_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_gilist_viridiplantae_protein_gi_code()}-process.sh'

    # return the script path
    return viridiplantae_protein_gi_gilist_script

#-------------------------------------------------------------------------------

def get_viridiplantae_protein_gi_gilist_starter():
    '''
    Get the starter path in the local computer to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # assign the starter path
    viridiplantae_protein_gi_gilist_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_gilist_viridiplantae_protein_gi_code()}-process-starter.sh'

    # return the starter path
    return viridiplantae_protein_gi_gilist_starter

#-------------------------------------------------------------------------------

def build_gene_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the NCBI Gene functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_gene_download_script())):
            os.makedirs(os.path.dirname(get_gene_download_script()))
        with open(get_gene_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi\n')
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
                script_file_id.write( 'function download_gene_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene to RefSeq file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GENE_GENE2REFSEQ_FILE \\\n')
                script_file_id.write( '            $GENE_GENE2REFSEQ_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading Gene Ontology file ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GENE_GENE2GO_FILE \\\n')
                script_file_id.write( '            $GENE_GENE2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_gene_name()}'
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
                script_file_id.write( 'download_gene_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gene_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gene_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the NCBI Gene functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_gene_download_starter())):
            os.makedirs(os.path.dirname(get_gene_download_starter()))
        with open(get_gene_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gene_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gene_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gene_download_script():
    '''
    Get the script path in the local computer to download the NCBI Gene functional annotation.
    '''

    # assign the script path
    gene_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_gene_code()}-process.sh'

    # return the script path
    return gene_download_script

#-------------------------------------------------------------------------------

def get_gene_download_starter():
    '''
    Get the script path in the local computer to download the NCBI Gene functional annotation.
    '''

    # assign the starter path
    gene_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_gene_code()}-process-starter.sh'

    # return the starter path
    return gene_download_starter

#-------------------------------------------------------------------------------

def build_gene_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load NCBI Gene functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_gene_load_script())):
            os.makedirs(os.path.dirname(get_gene_load_script()))
        with open(get_gene_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_gene_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-ncbi-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --dataset=gene \\\n')
                script_file_id.write( '            --gene2refseq=$GENE_GENE2REFSEQ_FILE \\\n')
                script_file_id.write( '            --gene2go=$GENE_GENE2GO_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-ncbi-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_gene_name()}'
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
                script_file_id.write( 'load_gene_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gene_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gene_load_starter(current_run_dir):
    '''
    Build the starter of the script to load NCBI Gene functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_gene_load_starter())):
            os.makedirs(os.path.dirname(get_gene_load_starter()))
        with open(get_gene_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gene_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gene_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gene_load_script():
    '''
    Get the script path to load NCBI Gene functional annotation into TOA database.
    '''

    # assign the script path
    gene_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_gene_code()}-process.sh'

    # return the script path
    return gene_load_script

#-------------------------------------------------------------------------------

def get_gene_load_starter():
    '''
    Get the starter path to load NCBI Gene functional annotation into TOA database.
    '''

    # assign the starter path
    gene_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_gene_code()}-process-starter.sh'

    # return the starter path
    return gene_load_starter

#-------------------------------------------------------------------------------

def build_interpro_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the InterPro functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_interpro_download_script())):
            os.makedirs(os.path.dirname(get_interpro_download_script()))
        with open(get_interpro_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$INTERPRO_DIR" ]; then mkdir --parents $INTERPRO_DIR; fi\n')
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
                script_file_id.write( 'function download_interpro_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading file of mappings of InterPro entries to Gene Ontology terms ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $INTERPRO_INTERPRO2GO_FILE \\\n')
                script_file_id.write( '            $INTERPRO_INTERPRO2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_interpro_name()}'
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
                script_file_id.write( 'download_interpro_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_interpro_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_interpro_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the InterPro functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_interpro_download_starter())):
            os.makedirs(os.path.dirname(get_interpro_download_starter()))
        with open(get_interpro_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_interpro_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_interpro_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_interpro_download_script():
    '''
    Get the script path to download the InterPro functional annotation.
    '''

    # assign the script path
    interpro_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_interpro_code()}-process.sh'

    # return the script path
    return interpro_download_script

#-------------------------------------------------------------------------------

def get_interpro_download_starter():
    '''
    Get the script path to download the InterPro functional annotation.
    '''

    # assign the starter path
    interpro_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_interpro_code()}-process-starter.sh'

    # return the starter path
    return interpro_download_starter

#-------------------------------------------------------------------------------

def build_interpro_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load InterPro functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_interpro_load_script())):
            os.makedirs(os.path.dirname(get_interpro_load_script()))
        with open(get_interpro_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_interpro_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-interpro-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --interpro2go=$INTERPRO_INTERPRO2GO_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-interpro-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_interpro_name()}'
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
                script_file_id.write( 'load_interpro_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_interpro_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_interpro_load_starter(current_run_dir):
    '''
    Build the starter of the script to load InterPro functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_interpro_load_starter())):
            os.makedirs(os.path.dirname(get_interpro_load_starter()))
        with open(get_interpro_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_interpro_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_interpro_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_interpro_load_script():
    '''
    Get the script path to load InterPro functional annotation into TOA database.
    '''

    # assign the script path
    interpro_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_interpro_code()}-process.sh'

    # return the script path
    return interpro_load_script

#-------------------------------------------------------------------------------

def get_interpro_load_starter():
    '''
    Get the starter path to load InterPro functional annotation into TOA database.
    '''

    # assign the starter path
    interpro_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_interpro_code()}-process-starter.sh'

    # return the starter path
    return interpro_load_starter

#-------------------------------------------------------------------------------

def build_go_download_script(cluster_name, current_run_dir):
    '''
    Build the script to download the Gene Ontology functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_go_download_script())):
            os.makedirs(os.path.dirname(get_go_download_script()))
        with open(get_go_download_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'if [ ! -d "$GO_DIR" ]; then mkdir --parents $GO_DIR; fi\n')
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
                script_file_id.write( 'function download_go_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading ontology ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GO_ONTOLOGY_FILE \\\n')
                script_file_id.write( '            $GO_ONTOLOGY_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading GO to Enzyme Commission (EC) ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document  $GO_EC2GO_FILE \\\n')
                script_file_id.write( '            $GO_EC2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading GO to KEGG ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GO_KEGG2GO_FILE \\\n')
                script_file_id.write( '            $GO_KEGG2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading GO to MetaCyc ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GO_METACYC2GO_FILE \\\n')
                script_file_id.write( '            $GO_METACYC2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Downloading GO to InterPro ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        wget \\\n')
                script_file_id.write( '            --quiet \\\n')
                script_file_id.write( '            --output-document $GO_INTERPRO2GO_FILE \\\n')
                script_file_id.write( '            $GO_INTERPRO2GO_FTP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error wget $RC; fi\n')
                script_file_id.write( '    echo "File is downloaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_download_go_name()}'
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
                script_file_id.write( 'download_go_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_go_download_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_go_download_starter(current_run_dir):
    '''
    Build the starter of the script to download the Gene Ontology functional annotation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_go_download_starter())):
            os.makedirs(os.path.dirname(get_go_download_starter()))
        with open(get_go_download_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_go_download_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_go_download_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_go_download_script():
    '''
    Get the script path to download the Gene Ontology functional annotation.
    '''

    # assign the script path
    go_download_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_go_code()}-process.sh'

    # return the script path
    return go_download_script

#-------------------------------------------------------------------------------

def get_go_download_starter():
    '''
    Get the script path to download the Gene Ontology functional annotation.
    '''

    # assign the starter path
    go_download_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_download_go_code()}-process-starter.sh'

    # return the starter path
    return go_download_starter

#-------------------------------------------------------------------------------

def build_go_load_script(cluster_name, current_run_dir):
    '''
    Build the script to load Gene Ontology functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_go_load_script())):
            os.makedirs(os.path.dirname(get_go_load_script()))
        with open(get_go_load_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
                script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
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
                script_file_id.write( 'function load_go_functional_annotation\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Loading functional annotation data into TOA database ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        load-go-data.py \\\n')
                script_file_id.write( '            --db=$TOA_DB \\\n')
                script_file_id.write( '            --ontology=$GO_ONTOLOGY_FILE \\\n')
                script_file_id.write( '            --ec2go=$GO_EC2GO_FILE \\\n')
                script_file_id.write( '            --kegg2go=$GO_KEGG2GO_FILE \\\n')
                script_file_id.write( '            --metacyc2go=$GO_METACYC2GO_FILE \\\n')
                script_file_id.write( '            --interpro2go=$GO_INTERPRO2GO_FILE \\\n')
                script_file_id.write( '            --verbose=N \\\n')
                script_file_id.write( '            --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error load-go-data.py $RC; fi\n')
                script_file_id.write( '    echo "Data are loaded."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_load_go_name()}'
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
                script_file_id.write( 'load_go_functional_annotation\n')
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_go_load_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_go_load_starter(current_run_dir):
    '''
    Build the starter of the script to load Gene Ontology functional annotation into TOA database.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_go_load_starter())):
            os.makedirs(os.path.dirname(get_go_load_starter()))
        with open(get_go_load_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_go_load_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_go_load_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_go_load_script():
    '''
    Get the script path to load Gene Ontology functional annotation into TOA database.
    '''

    # assign the script path
    go_load_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_go_code()}-process.sh'

    # return the script path
    return go_load_script

#-------------------------------------------------------------------------------

def get_go_load_starter():
    '''
    Get the starter path to load Gene Ontology functional annotation into TOA database.
    '''

    # assign the starter path
    go_load_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_load_go_code()}-process-starter.sh'

    # return the starter path
    return go_load_starter

#-------------------------------------------------------------------------------

def create_pipeline_config_file(pipeline_type, assembly_origin='NGSCLOUD', experiment_id='exp001', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS', reference_dataset_id='NONE', transcriptome_file='NONE', database_list=['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']):
    '''
    Create nucleotide pipeline config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the assembly software
    if assembly_origin == 'NGSCLOUD':
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
    elif assembly_origin == 'EXTERNAL':
        assembly_software = 'NONE'

    # get the order of the database in the list
    try:
        dicots_04_order = database_list.index('dicots_04') + 1
    except Exception as e:
        dicots_04_order = 0
    try:
        gymno_01_order = database_list.index('gymno_01') + 1
    except Exception as e:
        gymno_01_order = 0
    try:
        monocots_04_order = database_list.index('monocots_04') + 1
    except Exception as e:
        monocots_04_order = 0
    try:
        refseq_plant_order = database_list.index('refseq_plant') + 1
    except Exception as e:
        refseq_plant_order = 0
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        try:
            nt_order = database_list.index('nt') + 1
        except Exception as e:
            nt_order = 0
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        try:
            nr_order = database_list.index('nr') + 1
        except Exception as e:
            nr_order = 0

    # get the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = get_nucleotide_pipeline_config_file()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = get_aminoacid_pipeline_config_file()

    # create the transcript-filter config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(config_file)):
            os.makedirs(os.path.dirname(config_file))
        with open(config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write( '# The assembly files have to be located in:\n')
            file_id.write(f'#    {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id  (when the assembly origin is NGSCLOUD)\n')
            file_id.write(f'#    {xlib.get_cluster_reference_dir()}/reference_dataset_id  (when the assembly origin is EXTERNAL)\n')
            file_id.write( '# The experiment_id, database_dataset_id, assembly_dataset_id, reference_dataset_id and trasncriptome_file names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of NCBI BLAST+ and their meaning in "https://blast.ncbi.nlm.nih.gov/"\n')
            file_id.write( '# and the ones of DIAMOND in "http://www.diamondsearch.org/".\n')
            file_id.write( '#\n')
            file_id.write( '# In sections "BLAST+ parameters" and "DIAMOND parameters", the key "other_parameters_blast?" allows you to input additional\n')
            file_id.write( '# parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of BLAST+/DIAMOND and value-i a valid value of parameter-i, e.g. in BLAST+:\n')
            file_id.write( '#\n')
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '#    other_parameters_blastx = --ungapped; --best_hit_score_edge=0.1\n')
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '#    other_parameters_blastp = --lcase_masking; --matrix=BLOSUM62\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_origin = {assembly_origin}', f'# assembly origin: {get_assembly_origin_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification; NONE when the assembly origin is EXTERNAL'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}; NONE when the assembly origin is EXTERNAL'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification; NONE when the assembly origin is EXTERNAL'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification; NONE when the assembly origin is NGSCLOUD'))
            file_id.write( '{0:<50} {1}\n'.format(f'transcriptome_file = {transcriptome_file}', '# transcriptome file name; NONE when the assembly origin is NGSCLOUD'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the database parameters\n')
            file_id.write( '[database parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'gymno_01 = {gymno_01_order}', '# order of Gymno PLAZA 1.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format(f'dicots_04 = {dicots_04_order}', '# order of Dicots PLAZA 4.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format(f'monocots_04 = {monocots_04_order}', '# order of Monocots PLAZA 4.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format(f'refseq_plant = {refseq_plant_order}', '# order of NCBI RefSeq Plant in the annotation; 0 if it is not used'))
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '{0:<50} {1}\n'.format(f'nt = {nt_order}', '# order of NCBI BLAST database NT in the annotation; 0 if it is not used'))
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '{0:<50} {1}\n'.format(f'nr = {nr_order}', '# order of NCBI BLAST database NR in the annotation; 0 if it is not used'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the pipeline building parameters\n')
            file_id.write( '[pipeline parameters]\n')
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '{0:<50} {1}\n'.format(f'alignment_tool = {xlib.get_blastplus_name()}', f'# tool used in blastx alignments: {xlib.get_alignment_tool_code_list_text()}; blastn alignments will use BLAST+'))
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '{0:<50} {1}\n'.format(f'alignment_tool = {xlib.get_blastplus_name()}', f'# tool used in blastp alignments: {xlib.get_alignment_tool_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the NCBI BLAST+ parameters\n')
            file_id.write( '[BLAST+ parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'evalue = 1E-6', '# expectation value threshold for saving hits'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_target_seqs = 20', '# maximum number of aligned sequences to keep'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_hsps = 999999', '# maximum number of HSPs per subject sequence to save for each query'))
            file_id.write( '{0:<50} {1}\n'.format( 'qcov_hsp_perc = 0.0', '# alignments below the specified query coverage per HSPs are removed'))
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '{0:<50} {1}\n'.format( 'other_parameters_blastx = NONE', '# blastx additional parameters or NONE'))
                file_id.write( '{0:<50} {1}\n'.format( 'other_parameters_blastn = NONE', '# blastn additional parameters or NONE'))
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '{0:<50} {1}\n'.format( 'other_parameters_blastp = NONE', '# blastp additional parameters or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the DIAMOND parameters\n')
            file_id.write( '[DIAMOND parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'evalue = 1E-6', '# expectation value threshold for saving hits'))
            file_id.write( '{0:<50} {1}\n'.format( 'max-target-seqs = 20', '# maximum number of aligned sequences to keep'))
            file_id.write( '{0:<50} {1}\n'.format( 'max-hsps = 999999', '# maximum number of HSPs per subject sequence to save for each query'))
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '{0:<50} {1}\n'.format( 'other_parameters_blastx = NONE', '# blastx additional parameters or NONE'))
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '{0:<50} {1}\n'.format( 'other_parameters_blastp = NONE', '# blastp additional parameters or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {config_file} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_pipeline_config_file(pipeline_type, strict):
    '''
    Check a TOA pipeline config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = get_nucleotide_pipeline_config_file()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = get_aminoacid_pipeline_config_file()

    # get the option dictionary
    try:
        pipeline_option_dict = xlib.get_option_dict(config_file)
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in pipeline_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "assembly_origin"
            assembly_origin = pipeline_option_dict.get('identification', {}).get('assembly_origin', not_found)
            if assembly_origin == not_found:
                error_list.append('*** ERROR: the key "assembly_origin" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_origin, get_assembly_origin_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_origin" has to be {get_assembly_origin_code_list_text()}.')
                OK = False

            # check section "identification" - key "experiment_id"
            experiment_id = pipeline_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if experiment_id.upper() == 'NONE':
                    error_list.append('*** ERROR: the key "experiment_id" has to be a valid experiment identification when the assembly origin is NGSCLOUD.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if experiment_id.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "experiment_id" has to be NONE when the assembly origin is EXTERNAL.')
                    OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = pipeline_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                    error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()} when the assembly origin is NGSCLOUD.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if assembly_software.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_dataset_id" has to be NONE when the assembly origin is EXTERNAL.')
                    OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = pipeline_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                    error_list.append(f'*** ERROR: the key "assembly_software" has to start with {get_assembly_software_code_list_text()} when the assembly origin is NGSCLOUD.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if assembly_dataset_id.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_dataset_id" has to be NONE when the assembly origin is EXTERNAL.')
                    OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = pipeline_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be NONE when the assembly origin is EXTERNAL.')
                    OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = pipeline_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if reference_dataset_id.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "reference_dataset_id" has to be NONE when the assembly origin is NGSCLOUD.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if reference_dataset_id.upper() == 'NONE':
                    error_list.append('*** ERROR: the key "reference_dataset_id" has to be a valid reference dataset identification when the assembly origin is EXTERNAL.')
                    OK = False

            # check section "identification" - key "transcriptome_file"
            transcriptome_file = pipeline_option_dict.get('identification', {}).get('transcriptome_file', not_found)
            if transcriptome_file == not_found:
                error_list.append('*** ERROR: the key "transcriptome_file" is not found in the section "identification".')
                OK = False
            elif assembly_origin == 'NGSCLOUD':
                if transcriptome_file.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "transcriptome_file" has to be NONE when the assembly origin is NGSCLOUD.')
                    OK = False
            elif assembly_origin == 'EXTERNAL':
                if transcriptome_file.upper() == 'NONE':
                    error_list.append('*** ERROR: the key "transcriptome_file" has to be a valid reference file when the assembly origin is EXTERNAL.')
                    OK = False

        # check section "database parameters"
        if 'database parameters' not in sections_list:
            error_list.append('*** ERROR: the section "database parameters" is not found.')
            OK = False
        else:

            # check section "database parameters" - key "gymno_01"
            is_ok_gymno_01 = False
            gymno_01 = pipeline_option_dict.get('database parameters', {}).get('gymno_01', not_found)
            if gymno_01 == not_found:
                error_list.append('*** ERROR: the key "gymno_01" is not found in the section "database parameters".')
                OK = False
            elif not xlib.check_int(gymno_01, minimum=0, maximum=6):
                error_list.append('*** ERROR: the key "gymno_01" has to be an integer number between 0 and 6.')
                OK = False
            else:
                is_ok_gymno_01 = True

            # check section "database parameters" - key "dicots_04"
            is_ok_dicots_04 = False
            dicots_04 = pipeline_option_dict.get('database parameters', {}).get('dicots_04', not_found)
            if dicots_04 == not_found:
                error_list.append('*** ERROR: the key "dicots_04" is not found in the section "database parameters".')
                OK = False
            elif not xlib.check_int(dicots_04, minimum=0, maximum=6):
                error_list.append('*** ERROR: the key "dicots_04" has to be an integer number between 0 and 6.')
                OK = False
            else:
                is_ok_dicots_04 = True

            # check section "database parameters" - key "monocots_04"
            is_ok_monocots_04 = False
            monocots_04 = pipeline_option_dict.get('database parameters', {}).get('monocots_04', not_found)
            if monocots_04 == not_found:
                error_list.append('*** ERROR: the key "monocots_04" is not found in the section "database parameters".')
                OK = False
            elif not xlib.check_int(monocots_04, minimum=0, maximum=6):
                error_list.append('*** ERROR: the key "monocots_04" has to be an integer number between 0 and 6.')
                OK = False
            else:
                is_ok_monocots_04 = True

            # check section "database parameters" - key "refseq_plant"
            is_ok_refseq_plant = False
            refseq_plant = pipeline_option_dict.get('database parameters', {}).get('refseq_plant', not_found)
            if refseq_plant == not_found:
                error_list.append('*** ERROR: the key "refseq_plant" is not found in the section "database parameters".')
                OK = False
            elif not xlib.check_int(refseq_plant, minimum=0, maximum=6):
                error_list.append('*** ERROR: the key "refseq_plant" has to be an integer number between 0 and 6.')
                OK = False
            else:
                is_ok_refseq_plant = True

            # check section "database parameters" - key "nt"/"nr"
            is_ok_nx = False
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                nt = pipeline_option_dict.get('database parameters', {}).get('nt', not_found)
                if nt == not_found:
                    error_list.append('*** ERROR: the key "nt" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nt, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nt" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx = True
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                nr = pipeline_option_dict.get('database parameters', {}).get('nr', not_found)
                if nr == not_found:
                    error_list.append('*** ERROR: the key "nr" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nr, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nr" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx = True

            # check the order of databases
            if is_ok_gymno_01 and is_ok_dicots_04 and is_ok_monocots_04 and is_ok_refseq_plant and is_ok_nx:
                if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                    database_code_list = get_nucleotide_annotation_database_code_list()
                    database_order_list = [int(dicots_04), int(gymno_01), int(monocots_04), int(refseq_plant), int(nt)]
                    last_database_code = 'nt'
                elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                    database_code_list = get_aminoacid_annotation_database_code_list()
                    database_order_list = [int(dicots_04), int(gymno_01), int(monocots_04), int(refseq_plant), int(nr)]
                    last_database_code = 'nr'
                (OK2, error_list2) = check_database_order(database_code_list, database_order_list, last_database_code)
                if not OK2:
                    OK = False
                    error_list = error_list + error_list2

        # check section "pipeline parameters"
        if 'pipeline parameters' not in sections_list:
            error_list.append('*** ERROR: the section "pipeline parameters" is not found.')
            OK = False
        else:

            # check section "pipeline parameters" - key "alignment_tool"
            alignment_tool = pipeline_option_dict.get('pipeline parameters', {}).get('alignment_tool', not_found)
            if alignment_tool == not_found:
                error_list.append('*** ERROR: the key "alignment_tool" is not found in the section "pipeline parameters".')
                OK = False
            elif not xlib.check_code(alignment_tool, xlib.get_alignment_tool_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "alignment_tool" has to be {xlib.get_alignment_tool_code_list_text()}.')
                OK = False

            # check section "pipeline parameters" - key "threads"
            threads = pipeline_option_dict.get('pipeline parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "pipeline parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

        # check section "BLAST+ parameters"
        if 'BLAST+ parameters' not in sections_list:
            error_list.append('*** ERROR: the section "BLAST+ parameters" is not found.')
            OK = False
        else:

            # check section "BLAST+ parameters" - key "evalue"
            blastplus_evalue = pipeline_option_dict.get('BLAST+ parameters', {}).get('evalue', not_found)
            if blastplus_evalue == not_found:
                error_list.append('*** ERROR: the key "evalue" is not found in the section "BLAST+ parameters".')
                OK = False
            elif not xlib.check_float(blastplus_evalue, minimum=0., mne=1E-12):
                error_list.append('*** ERROR: the key "evalue" has to be a float number greater than to 0.0.')
                OK = False

            # check section "BLAST+ parameters" - key "max_target_seqs"
            blastplus_max_target_seqs = pipeline_option_dict.get('BLAST+ parameters', {}).get('max_target_seqs', not_found)
            if blastplus_max_target_seqs == not_found:
                error_list.append('*** ERROR: the key "max_target_seqs" is not found in the section "BLAST+ parameters".')
                OK = False
            elif not xlib.check_int(blastplus_max_target_seqs, minimum=1):
                error_list.append('*** ERROR: the key "max_target_seqsr" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BLAST+ parameters" - key "max_hsps"
            blastplus_max_hsps = pipeline_option_dict.get('BLAST+ parameters', {}).get('max_hsps', not_found)
            if blastplus_max_hsps == not_found:
                error_list.append('*** ERROR: the key "max_hsps" is not found in the section "BLAST+ parameters".')
                OK = False
            elif not xlib.check_int(blastplus_max_hsps, minimum=1):
                error_list.append('*** ERROR: the key "max_hsps" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BLAST+ parameters" - key "qcov_hsp_perc"
            blastplus_qcov_hsp_perc = pipeline_option_dict.get('BLAST+ parameters', {}).get('qcov_hsp_perc', not_found)
            if blastplus_qcov_hsp_perc == not_found:
                error_list.append('*** ERROR: the key "qcov_hsp_perc" is not found in the section "BLAST+ parameters".')
                OK = False
            elif not xlib.check_float(blastplus_qcov_hsp_perc, minimum=0., maximum=100., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "qcov_hsp_perc" has to be a float number greater than or equal to 0.0 and less than 100.0.')
                OK = False

            # check section "BLAST+ parameters" - key "other_parameters_blastx"
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                not_allowed_parameters_list = ['num_threads', 'db', 'query', 'evalue', 'max_target_seqs', 'max_hsps', 'qcov_hsp_perc', 'outfmt', 'out']
                blastplus_other_parameters_blastx = pipeline_option_dict.get('BLAST+ parameters', {}).get('other_parameters_blastx', not_found)
                if blastplus_other_parameters_blastx == not_found:
                    error_list.append('*** ERROR: the key "other_parameters_blastx" is not found in the section "BLAST+ parameters".')
                    OK = False
                elif blastplus_other_parameters_blastx.upper() != 'NONE':
                    (OK, error_list2) = xlib.check_parameter_list(blastplus_other_parameters_blastx, "other_parameters_blastx", not_allowed_parameters_list)
                    error_list = error_list + error_list2

            # check section "BLAST+ parameters" - key "other_parameters_blastn"
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                not_allowed_parameters_list = ['num_threads', 'db', 'query', 'evalue', 'max_target_seqs', 'max_hsps', 'qcov_hsp_perc', 'outfmt', 'out']
                blastplus_other_parameters_blastn = pipeline_option_dict.get('BLAST+ parameters', {}).get('other_parameters_blastn', not_found)
                if blastplus_other_parameters_blastn == not_found:
                    error_list.append('*** ERROR: the key "other_parameters_blastn" is not found in the section "BLAST+ parameters".')
                    OK = False
                elif blastplus_other_parameters_blastn.upper() != 'NONE':
                    (OK, error_list2) = xlib.check_parameter_list(blastplus_other_parameters_blastn, "other_parameters_blastn", not_allowed_parameters_list)
                    error_list = error_list + error_list2

            # check section "BLAST+ parameters" - key "other_parameters_blastp"
            if pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                not_allowed_parameters_list = ['num_threads', 'db', 'query', 'evalue', 'max_target_seqs', 'max_hsps', 'qcov_hsp_perc', 'outfmt', 'out']
                blastplus_other_parameters_blastp = pipeline_option_dict.get('BLAST+ parameters', {}).get('other_parameters_blastp', not_found)
                if blastplus_other_parameters_blastp == not_found:
                    error_list.append('*** ERROR: the key "other_parameters_blastp" is not found in the section "BLAST+ parameters".')
                    OK = False
                elif blastplus_other_parameters_blastp.upper() != 'NONE':
                    (OK, error_list2) = xlib.check_parameter_list(blastplus_other_parameters_blastp, "other_parameters_blastp", not_allowed_parameters_list)
                    error_list = error_list + error_list2

        # check section "DIAMOND parameters"
        if 'DIAMOND parameters' not in sections_list:
            error_list.append('*** ERROR: the section "DIAMOND parameters" is not found.')
            OK = False
        else:

            # check section "DIAMOND parameters" - key "evalue"
            diamond_evalue = pipeline_option_dict.get('DIAMOND parameters', {}).get('evalue', not_found)
            if diamond_evalue == not_found:
                error_list.append('*** ERROR: the key "evalue" is not found in the section "DIAMOND parameters".')
                OK = False
            elif not xlib.check_float(diamond_evalue, minimum=0., mne=1E-12):
                error_list.append('*** ERROR: the key "evalue" has to be a float number greater than to 0.0.')
                OK = False

            # check section "DIAMOND parameters" - key "max-target-seqs"
            diamond_max_target_seqs = pipeline_option_dict.get('DIAMOND parameters', {}).get('max-target-seqs', not_found)
            if diamond_max_target_seqs == not_found:
                error_list.append('*** ERROR: the key "max-target-seqs" is not found in the section "DIAMOND parameters".')
                OK = False
            elif not xlib.check_int(diamond_max_target_seqs, minimum=1):
                error_list.append('*** ERROR: the key "max-target-seqsr" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "DIAMOND parameters" - key "max-hsps"
            diamond_max_hsps = pipeline_option_dict.get('DIAMOND parameters', {}).get('max-hsps', not_found)
            if diamond_max_hsps == not_found:
                error_list.append('*** ERROR: the key "max-hsps" is not found in the section "DIAMOND parameters".')
                OK = False
            elif not xlib.check_int(diamond_max_hsps, minimum=1):
                error_list.append('*** ERROR: the key "max-hsps" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "DIAMOND parameters" - key "other_parameters_blastx"
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                not_allowed_parameters_list = ['threads', 'db', 'query', 'evalue', 'max-target-seqs', 'max-hsps', 'outfmt', 'out']
                diamond_other_parameters_blastx = pipeline_option_dict.get('DIAMOND parameters', {}).get('other_parameters_blastx', not_found)
                if diamond_other_parameters_blastx == not_found:
                    error_list.append('*** ERROR: the key "other_parameters_blastx" is not found in the section "DIAMOND parameters".')
                    OK = False
                elif diamond_other_parameters_blastx.upper() != 'NONE':
                    (OK, error_list2) = xlib.check_parameter_list(diamond_other_parameters_blastx, "other_parameters_blastx", not_allowed_parameters_list)
                    error_list = error_list + error_list2

            # check section "DIAMOND parameters" - key "other_parameters_blastp"
            if pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                not_allowed_parameters_list = ['threads', 'db', 'query', 'evalue', 'max-target-seqs', 'max-hsps', 'outfmt', 'out']
                diamond_other_parameters_blastp = pipeline_option_dict.get('DIAMOND parameters', {}).get('other_parameters_blastp', not_found)
                if diamond_other_parameters_blastp == not_found:
                    error_list.append('*** ERROR: the key "other_parameters_blastp" is not found in the section "DIAMOND parameters".')
                    OK = False
                elif diamond_other_parameters_blastp.upper() != 'NONE':
                    (OK, error_list2) = xlib.check_parameter_list(diamond_other_parameters_blastp, "other_parameters_blastp", not_allowed_parameters_list)
                    error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_toa_process_pipeline_nucleotide_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def check_database_order(database_code_list, database_order_list, last_database_code):
    '''
    Check database order in a selelected database list.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the database dictionary
    database_dict = {}

    # set the last database order
    last_database_order = database_order_list[database_code_list.index(last_database_code)]

    # define the function to add database to the dictionary
    def add_database_to_dict(dict, name, value):
        is_db_added = True
        if value > 0:
            old_name = database_dict.get(value, '')
            if old_name == '':
                database_dict[value] = name
            else:
                is_db_added = False
        return (is_db_added, dict)

    # add databases to the dictionary
    for i in range(len(database_code_list)):
        (is_db_added, database_dict) = add_database_to_dict(database_dict, database_code_list[i], database_order_list[i])
        if not is_db_added:
            error_list.append('*** ERROR: there are databases with the same order number.')
            OK = False
            break

    # if all databases are added to the dictionary, check the order
    if is_db_added:
        order_list = sorted(database_dict.keys())
        if database_dict == {}:
            error_list.append('*** ERROR: there are not databases selected to annotate.')
            OK = False
        elif order_list[0] != 1:
            error_list.append('*** ERROR: the first database has to have 1 in the order to annotate.')
            OK = False
        elif order_list[len(order_list) -1] != len(order_list):
            error_list.append('*** ERROR: the database orders has to have sequencial numbers.')
            OK = False
        elif last_database_order > 0 and last_database_order != len(order_list):
            error_list.append(f'*** ERROR: the {last_database_code} order has to be the last one.')
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_selected_database_list(pipeline_type):
    '''
    Get the selected database list order by the database order to annotate of a pipeline config file.
    '''

    # initialize the selected database list
    selected_database_list = []

    # set the candidate database list 
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        candidate_database_list = get_nucleotide_annotation_database_code_list()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        candidate_database_list = get_aminoacid_annotation_database_code_list()

    # get the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = get_nucleotide_pipeline_config_file()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = get_aminoacid_pipeline_config_file()

    # get the option dictionary
    pipeline_option_dict = xlib.get_option_dict(config_file)

    # get the select database list order by the database order to annotate
    temp_dict = {}
    for database_code in candidate_database_list:
        if int(pipeline_option_dict['database parameters'][database_code]) > 0:
            temp_dict[int(pipeline_option_dict['database parameters'][database_code])] = database_code
    for database_order in sorted(temp_dict.keys()):
        selected_database_list.append(temp_dict[database_order])

    # return the selected database list
    return selected_database_list

#-------------------------------------------------------------------------------

def get_nucleotide_pipeline_config_file():
    '''
    Get the nucleotide pipeline config file path.
    '''

    # assign the nucleotide pipeline config file path
    nucleotide_pipeline_config_file = f'{xlib.get_config_dir()}/{xlib.get_toa_process_pipeline_nucleotide_code()}-config.txt'

    # return the nucleotide pipeline config file path
    return nucleotide_pipeline_config_file

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_config_file():
    '''
    Get the amino acid pipeline config file path.
    '''

    # assign the amino acid pipeline config file path
    aminoacid_pipeline_config_file = f'{xlib.get_config_dir()}/{xlib.get_toa_process_pipeline_aminoacid_code()}-config.txt'

    # return the amino acid pipeline config file path
    return aminoacid_pipeline_config_file

#-------------------------------------------------------------------------------

def run_pipeline_process(cluster_name, pipeline_type, log, function=None):
    '''
    Run a TOA pipeline process.
    '''

    # initialize the control variable
    OK = True

    # get the dictionary of TOA configuration
    toa_config_dict = get_toa_config_dict()

    # get the config file
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        config_file = get_nucleotide_pipeline_config_file()
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        config_file = get_aminoacid_pipeline_config_file()

    # get the option dictionary
    pipeline_option_dict = xlib.get_option_dict(config_file)

    # build the sentence to set the cluster PATH to check requirements
    path_sentence = f'export PATH={toa_config_dict["MINICONDA3_BIN_DIR"]}:{toa_config_dict["TOA_DIR"]}:$PATH'

    # get the selected database list
    selected_database_list = get_selected_database_list(pipeline_type)

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the TOA config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_toa_name()} config file ...\n')
    OK = os.path.isfile(get_toa_config_file())
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file does not exist.\n')
        log.write('Please recreate this file.\n')
        OK = False

    # check the pipeline config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            log.write(f'Checking the {xlib.get_toa_process_pipeline_nucleotide_name()} config file ...\n')
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            log.write(f'Checking the {xlib.get_toa_process_pipeline_aminoacid_name()} config file ...\n')
        (OK, error_list) = check_pipeline_config_file(pipeline_type, strict=True)
        if OK:
            log.write('The file is OK.\n')
        else:
            log.write('*** ERROR: The config file is not valid.\n')
            log.write('Please correct this file or recreate it.\n')

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

    # check the TOA is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_toa_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_toa_name()} is not installed.\n')
            OK = False

    # check the basic data load in TOA databasep
    if OK:
        check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=basic'
        command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: The basic data load in TOA database is wrong.\n')
            OK = False
        else:
            log.write('... basic data load in TOA database is OK ...\n')

    # check the Gymno PLAZA 1.0 proteome
    if OK:
        if 'gymno_01' in selected_database_list:
            command = f'[ -d {toa_config_dict["GYMNO_01_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Gymno PLAZA 1.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Gymno PLAZA 1.0 proteome is OK ...\n')

    # check the Gymno PLAZA 1.0 data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=gymno_01'
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Gymno PLAZA 1.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Gymno PLAZA 1.0 load data in TOA database is OK ...\n')

    # check the Dicots PLAZA 4.0 proteome
    if OK:
        if 'dicots_04' in selected_database_list:
            command = f'[ -d {toa_config_dict["DICOTS_04_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Dicots PLAZA 4.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Dicots PLAZA 4.0 proteome is OK ...\n')

    # check the Dicots PLAZA 4.0 data load in TOA database
    if OK:
        if 'dicots_04' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=dicots_04'
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Dicots PLAZA 4.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Dicots PLAZA 4.0 data load in TOA database is OK ...\n')

    # check the Monocots PLAZA 4.0 proteome
    if OK:
        if 'monocots_04' in selected_database_list:
            command = f'[ -d {toa_config_dict["MONOCOTS_04_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Monocots PLAZA 4.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Monocots PLAZA 4.0 proteome is OK ...\n')

    # check the Monocots PLAZA 4.0 data load in TOA database
    if OK:
        if 'monocots_04' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=monocots_04'
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Monocots PLAZA 4.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Monocots PLAZA 4.0 data load in TOA database is OK ...\n')

    # check the NCBI RefSeq Plant proteome
    if OK:
        if 'refseq_plant' in selected_database_list:
            command = f'[ -d {toa_config_dict["REFSEQ_PLANT_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI RefSeq Plant proteome is not built.\n')
                OK = False
            else:
                log.write('... NCBI RefSeq Plant proteome is OK ...\n')

    # check the NCBI BLAST database NT
    if OK:
        if 'nt' in selected_database_list and pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            command = f'[ -d {toa_config_dict["NT_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI BLAST database NT for BLAST+ is not built.\n')
                OK = False
            else:
                log.write('... NCBI BLAST database NTfor BLAST+ is OK ...\n')

    # check the NCBI BLAST database NR
    if OK:
        if 'nr' in selected_database_list and pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            if pipeline_option_dict['pipeline parameters']['alignment_tool'] == 'BLAST+':
                command = f'[ -d {toa_config_dict["NR_BLASTPLUS_DB_DIR"]} ] && echo RC=0 || echo RC=1'
                (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
                if stdout[len(stdout) - 1] != 'RC=0':
                    log.write('*** ERROR: NCBI BLAST database NR for BLAST+ is not built.\n')
                    OK = False
                else:
                    log.write('... NCBI BLAST database NR for BLAST+ is OK ...\n')
            elif pipeline_option_dict['pipeline parameters']['alignment_tool'] == 'DIAMOND':
                command = f'[ -d {toa_config_dict["NR_DIAMOND_DB_DIR"]} ] && echo RC=0 || echo RC=1'
                (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
                if stdout[len(stdout) - 1] != 'RC=0':
                    log.write('*** ERROR: NCBI BLAST database NR for DIAMOND is not built.\n')
                    OK = False
                else:
                    log.write('... NCBI BLAST database NR for DIAMOND is OK ...\n')

    # check the NCBI Gene data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=gene'.format()
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI Gene load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... NCBI Gene data load in TOA database is OK ...\n')

    # check the InterPro data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=interpro'
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: InterPro load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... InterPro data load in TOA database is OK ...\n')

    # check the Gene Ontology data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = f'check-data-load.py --db={toa_config_dict["TOA_DB"]} --group=go'
            command = f'{path_sentence}; {check_sentence} && echo RC=0 || echo RC=1'
            (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Gene Ontology data load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Gene Ontology load in TOA database is OK ...\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')

        # nucleotide pipelines
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_pipeline_dir(), xlib.get_toa_process_pipeline_nucleotide_code())

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_pipeline_dir(), xlib.get_toa_process_pipeline_aminoacid_code())

        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')

        # nucleotide pipelines
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            script = get_nucleotide_pipeline_script()
            log.write(f'Building the process script {script} ...\n')
            (OK, error_list) = build_nucleotide_pipeline_script(cluster_name, current_run_dir)

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            script = get_aminoacid_pipeline_script()
            log.write(f'Building the process script {script} ...\n')
            (OK, error_list) = build_aminoacid_pipeline_script(cluster_name, current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {script} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(script)}'
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(script)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(script)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')

        # nucleotide pipelines
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            starter = get_nucleotide_pipeline_starter()
            log.write(f'Building the process starter {starter} ...\n')
            (OK, error_list) = build_nucleotide_pipeline_starter(current_run_dir)

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            starter = get_aminoacid_pipeline_starter()
            log.write(f'Building the process starter {starter} ...\n')
            (OK, error_list) = build_aminoacid_pipeline_starter(current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {starter} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(starter)}'
        (OK, error_list) = xssh.put_file(sftp_client, starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(starter)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(starter)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(starter)} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(starter), log)

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

def build_nucleotide_pipeline_script(cluster_name, current_run_dir):
    '''
    Build the script to process a nucleotide pipeline.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # get the pipeline option dictionary
    pipeline_option_dict = xlib.get_option_dict(get_nucleotide_pipeline_config_file())

    # get the options
    assembly_origin = pipeline_option_dict['identification']['assembly_origin']
    experiment_id = pipeline_option_dict['identification']['experiment_id']
    assembly_software = pipeline_option_dict['identification']['assembly_software']
    assembly_dataset_id = pipeline_option_dict['identification']['assembly_dataset_id']
    assembly_type = pipeline_option_dict['identification']['assembly_type']
    reference_dataset_id = pipeline_option_dict['identification']['reference_dataset_id']
    transcriptome_file = pipeline_option_dict['identification']['transcriptome_file']
    alignment_tool = pipeline_option_dict['pipeline parameters']['alignment_tool']
    threads = pipeline_option_dict['pipeline parameters']['threads']
    blastplus_evalue = pipeline_option_dict['BLAST+ parameters']['evalue']
    blastplus_max_target_seqs = pipeline_option_dict['BLAST+ parameters']['max_target_seqs']
    blastplus_max_hsps = pipeline_option_dict['BLAST+ parameters']['max_hsps']
    blastplus_qcov_hsp_perc = pipeline_option_dict['BLAST+ parameters']['qcov_hsp_perc']
    blastplus_other_parameters_blastx = pipeline_option_dict['BLAST+ parameters']['other_parameters_blastx']
    blastplus_other_parameters_blastn = pipeline_option_dict['BLAST+ parameters']['other_parameters_blastn']
    diamond_evalue = pipeline_option_dict['DIAMOND parameters']['evalue']
    diamond_max_target_seqs = pipeline_option_dict['DIAMOND parameters']['max-target-seqs']
    diamond_max_hsps = pipeline_option_dict['DIAMOND parameters']['max-hsps']
    diamond_other_parameters_blastx = pipeline_option_dict['DIAMOND parameters']['other_parameters_blastx']

    # get the all selected database list
    database_list = get_selected_database_list(xlib.get_toa_process_pipeline_nucleotide_code())

    # get the database type dictionary
    database_type_dict = {}
    for i in range(len(database_list)):
        if database_list[i] in ['gymno_01', 'dicots_04', 'monocots_04']:
            database_type_dict[database_list[i]] = 'PLAZA'
        elif database_list[i] == 'refseq_plant':
            database_type_dict[database_list[i]] = 'REFSEQ'
        elif database_list[i] in ['nt']:
            database_type_dict[database_list[i]] = 'NT'

    # get the database list for the merger of plant annotation files
    database_list2 = database_list.copy()
    for i in range(len(database_list2)):
        if database_list2[i] == 'nt':
            database_list2[i] = 'nt_viridiplantae'

    # set the transcriptome file path
    if OK:
        if assembly_origin == 'NGSCLOUD':
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
        elif assembly_origin == 'EXTERNAL':
            transcriptome_file = xlib.get_cluster_reference_file(reference_dataset_id, transcriptome_file)

    # get the non annotation file list
    non_annotation_file_list = []
    for database in database_list:
        non_annotation_file_list.append(f'${database.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE')

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nucleotide_pipeline_script())):
            os.makedirs(os.path.dirname(get_nucleotide_pipeline_script()))
        with open(get_nucleotide_pipeline_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '# transcriptome file\n')
            script_file_id.write(f'TRANSCRIPTOME_FILE={transcriptome_file}\n')
            script_file_id.write( '\n')
            script_file_id.write( '# pipeline parameters\n')
            script_file_id.write(f'THREADS={threads}\n')
            script_file_id.write( '\n')
            script_file_id.write( '# BLAST+ parameters\n')
            script_file_id.write(f'BLASTPLUS_EVALUE={blastplus_evalue}\n')
            script_file_id.write(f'BLASTPLUS_MAX_TARGET_SEQS={blastplus_max_target_seqs}\n')
            script_file_id.write(f'BLASTPLUS_MAX_HSPS={blastplus_max_hsps}\n')
            script_file_id.write(f'BLASTPLUS_QCOV_HSP_PERC={blastplus_qcov_hsp_perc}\n')
            script_file_id.write(f'BLASTPLUS_OTHER_PARAMETERS_BLASTX="{blastplus_other_parameters_blastx}"\n')
            script_file_id.write(f'BLASTPLUS_OTHER_PARAMETERS_BLASTN="{blastplus_other_parameters_blastn}"\n')
            script_file_id.write( '\n')
            script_file_id.write( '# DIAMOND parameters\n')
            script_file_id.write(f'DIAMOND_EVALUE={diamond_evalue}\n')
            script_file_id.write(f'DIAMOND_MAX_TARGET_SEQS={diamond_max_target_seqs}\n')
            script_file_id.write(f'DIAMOND_MAX_HSPS={diamond_max_hsps}\n')
            script_file_id.write(f'DIAMOND_OTHER_PARAMETERS_BLASTX="{diamond_other_parameters_blastx}"\n')
            script_file_id.write( '\n')
            script_file_id.write( '# output directory\n')
            script_file_id.write(f'OUTPUT_DIR={current_run_dir}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
            script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
            script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'mkdir --parents $STATS_DIR\n')
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
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "TRANSCRIPTOME FILE: $TRANSCRIPTOME_FILE"\n')
            database_list_text = ','.join(database_list)
            script_file_id.write(f'    echo "ALIGNMENT DATASETS: {database_list_text}"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "ALIGNMENT TOOLS (blastp and blastx): {alignment_tool}"\n')
            script_file_id.write(f'    echo "ALIGNMENT TOOLS (blastn): {xlib.get_blastplus_name()}"\n')
            script_file_id.write( '    echo "THREADS: $THREADS"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "BLAST+ EVALUE: $BLASTPLUS_EVALUE"\n')
            script_file_id.write( '    echo "BLAST+ MAX_TARGET_SEQS: $BLASTPLUS_MAX_TARGET_SEQS"\n')
            script_file_id.write( '    echo "BLAST+ MAX_HSPS: $BLASTPLUS_MAX_HSPS"\n')
            script_file_id.write( '    echo "BLAST+ QCOV_HSP_PERC: $BLASTPLUS_QCOV_HSP_PERC"\n')
            script_file_id.write( '    echo "BLAST+ BLASTX OTHER PARAMETERS: $BLASTPLUS_OTHER_PARAMETERS_BLASTX"\n')
            script_file_id.write( '    echo "BLAST+ BLASTN OTHER PARAMETERS: $BLASTPLUS_OTHER_PARAMETERS_BLASTN"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "DIAMOND EVALUE: $DIAMOND_EVALUE"\n')
            script_file_id.write( '    echo "DIAMOND MAX_TARGET_SEQS: $DIAMOND_MAX_TARGET_SEQS"\n')
            script_file_id.write( '    echo "DIAMOND MAX_HSPS: $DIAMOND_MAX_HSPS"\n')
            script_file_id.write( '    echo "DIAMOND BLASTX OTHER PARAMETERS: $DIAMOND_OTHER_PARAMETERS_BLASTX"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function reidentify_transcript_sequences\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/reidentify_transcript_sequences.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "RE-INDENTIFY SEQUENCES OF THE TRANSCRIPTOME FILE"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Re-identifing sequences ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            reid-fasta-file.py \\\n')
            script_file_id.write( '                --fasta=$TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --type=NT \\\n')
            script_file_id.write( '                --out=$REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error reid-fasta-file.py $RC; fi\n')
            script_file_id.write( '        echo "Sequences are re-identified."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            for i in range(len(database_list)):
                current_code = database_list[i]
                previus_code = database_list[i - 1] if i > 0 else ''
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function align_transcripts_{current_code}_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/align_transcripts_{current_code}_proteome.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "ALIGNMENT OF TRANSCRIPTS TO {current_code.upper()} PROTEOME"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Aligning transcripts ..."\n')
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    if alignment_tool == xlib.get_blastplus_name():
                        script_file_id.write( '        source activate blast\n')
                        script_file_id.write(f'        export BLASTDB=${current_code.upper()}_BLASTPLUS_DB_DIR\n')
                        script_file_id.write( '        /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write( '            blastx \\\n')
                        script_file_id.write( '                -num_threads $THREADS \\\n')
                        script_file_id.write(f'                -db ${current_code.upper()}_BLASTPLUS_DB_NAME \\\n')
                        if i == 0:
                            script_file_id.write( '                -query $REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
                        else:
                            script_file_id.write(f'                -query ${previus_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
                        script_file_id.write( '                -evalue $BLASTPLUS_EVALUE \\\n')
                        script_file_id.write( '                -max_target_seqs $BLASTPLUS_MAX_TARGET_SEQS \\\n')
                        script_file_id.write( '                -max_hsps $BLASTPLUS_MAX_HSPS \\\n')
                        script_file_id.write( '                -qcov_hsp_perc $BLASTPLUS_QCOV_HSP_PERC \\\n')
                        if blastplus_other_parameters_blastx.upper() != 'NONE':
                            parameter_list = [x.strip() for x in blastplus_other_parameters_blastx.split(';')]
                            for parameter in parameter_list:
                                if parameter.find('=') > 0:
                                    pattern = r'^--(.+)=(.+)$'
                                    mo = re.search(pattern, parameter)
                                    parameter_name = mo.group(1).strip()
                                    parameter_value = mo.group(2).strip()
                                    script_file_id.write(f'                -{parameter_name} {parameter_value} \\\n')
                                else:
                                    pattern = r'^--(.+)$'
                                    mo = re.search(pattern, parameter)
                                    parameter_name = mo.group(1).strip()
                                    script_file_id.write(f'                -{parameter_name} \\\n')
                        script_file_id.write( '                -outfmt 5 \\\n')
                        script_file_id.write(f'                -out ${current_code.upper()}_BLAST_XML\n')
                        script_file_id.write( '        RC=$?\n')
                        script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastx $RC; fi\n')
                        script_file_id.write( '        echo "Alignment is done."\n')
                        script_file_id.write( '        conda deactivate\n')
                    elif alignment_tool == xlib.get_diamond_name():
                        if i == 0:
                            script_file_id.write( '        if [[ -s $REIDENTIFIED_TRANSCRIPTOME_FILE ]]; then\n')
                        else:
                            script_file_id.write(f'        if [[ -s ${previus_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE ]]; then\n')
                        script_file_id.write( '            source activate diamond\n')
                        script_file_id.write( '            /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write( '                diamond blastx \\\n')
                        script_file_id.write( '                    --threads $THREADS \\\n')
                        script_file_id.write(f'                    --db ${current_code.upper()}_DIAMOND_DB_FILE \\\n')
                        if i == 0:
                            script_file_id.write( '                    --query $REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
                        else:
                            script_file_id.write(f'                    --query ${previus_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
                        script_file_id.write( '                    --evalue $DIAMOND_EVALUE \\\n')
                        script_file_id.write( '                    --max-target-seqs $DIAMOND_MAX_TARGET_SEQS \\\n')
                        script_file_id.write( '                    --max-hsps $DIAMOND_MAX_HSPS \\\n')
                        if diamond_other_parameters_blastx.upper() != 'NONE':
                            parameter_list = [x.strip() for x in diamond_other_parameters_blastx.split(';')]
                            for parameter in parameter_list:
                                if parameter.find('=') > 0:
                                    pattern = r'^--(.+)=(.+)$'
                                    mo = re.search(pattern, parameter)
                                    parameter_name = mo.group(1).strip()
                                    parameter_value = mo.group(2).strip()
                                    script_file_id.write(f'                    --{parameter_name} {parameter_value} \\\n')
                                else:
                                    pattern = r'^--(.+)$'
                                    mo = re.search(pattern, parameter)
                                    parameter_name = mo.group(1).strip()
                                    script_file_id.write(f'                    --{parameter_name} \\\n')
                        script_file_id.write( '                    --outfmt 5 \\\n')
                        script_file_id.write(f'                    --out ${current_code.upper()}_BLAST_XML\n')
                        script_file_id.write( '            RC=$?\n')
                        script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error diamond-blastx $RC; fi\n')
                        script_file_id.write( '        else\n')
                        script_file_id.write(f'            touch ${current_code.upper()}_BLAST_XML\n')
                        script_file_id.write( '        fi\n')
                        script_file_id.write( '        echo "Alignment is done."\n')
                        script_file_id.write( '        conda deactivate\n')
                elif current_code == 'nt':
                    script_file_id.write( '        source activate blast\n')
                    script_file_id.write(f'        export BLASTDB=$NT_BLASTPLUS_DB_DIR\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            blastn \\\n')
                    script_file_id.write( '                -num_threads $THREADS \\\n')
                    script_file_id.write(f'                -db $NT_BLASTPLUS_DB_NAME \\\n')
                    if i == 0:
                        script_file_id.write( '                -query $REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
                    else:
                        script_file_id.write(f'                -query ${previus_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
                    script_file_id.write( '                -evalue $BLASTPLUS_EVALUE \\\n')
                    script_file_id.write( '                -max_target_seqs $BLASTPLUS_MAX_TARGET_SEQS \\\n')
                    script_file_id.write( '                -max_hsps $BLASTPLUS_MAX_HSPS \\\n')
                    script_file_id.write( '                -qcov_hsp_perc $BLASTPLUS_QCOV_HSP_PERC \\\n')
                    if blastplus_other_parameters_blastn.upper() != 'NONE':
                        parameter_list = [x.strip() for x in blastplus_other_parameters_blastn.split(';')]
                        for parameter in parameter_list:
                            if parameter.find('=') > 0:
                                pattern = r'^--(.+)=(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                parameter_value = mo.group(2).strip()
                                script_file_id.write(f'                -{parameter_name} {parameter_value} \\\n')
                            else:
                                pattern = r'^--(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                script_file_id.write(f'                -{parameter_name} \\\n')
                    script_file_id.write( '                -outfmt 5 \\\n')
                    script_file_id.write(f'                -out ${current_code.upper()}_BLAST_XML\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastn $RC; fi\n')
                    script_file_id.write( '        echo "Alignment is done."\n')
                    script_file_id.write( '        conda deactivate\n')
                if len(database_list) == 1:
                    script_file_id.write( '        echo "Restoring sequence identifications in alignment file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            restore-ids.py \\\n')
                    script_file_id.write(f'                --in=${current_code.upper()}_BLAST_XML \\\n')
                    script_file_id.write( '                --format=XML \\\n')
                    script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                    script_file_id.write( '                --relationships2=NONE \\\n')
                    script_file_id.write( '                --out=$MERGED_BLAST_XML \\\n')
                    script_file_id.write( '                --verbose=N \\\n')
                    script_file_id.write( '                --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi\n')
                    script_file_id.write( '        echo "Identifications are restored."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function load_alignment_{current_code}_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/load_alignment_{current_code}_proteome.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "LOAD OF TRANSCRIPT ALIGNMENT TO {current_code.upper()} PROTEOME INTO TOA DATABASE"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Loading alignmnet data ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            load-blast-data.py \\\n')
                script_file_id.write( '                --db=$TOA_DB \\\n')
                script_file_id.write(f'                --dataset={current_code} \\\n')
                script_file_id.write( '                --format=5 \\\n')
                script_file_id.write(f'                --blast=${current_code.upper()}_BLAST_XML \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error load-blast-data.py $RC; fi\n')
                script_file_id.write( '        echo "Data are loaded."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function annotate_transcripts_{current_code}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/annotate_transcripts_{current_code}.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "ANNOTATION OF TRANSCRIPTS WITH {current_code.upper()}"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Annotating transcripts ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            annotate-sequences.py \\\n')
                script_file_id.write( '                --db=$TOA_DB \\\n')
                script_file_id.write(f'                --dataset={current_code} \\\n')
                if current_code == 'nt':
                    script_file_id.write(f'                --aligner=BLAST+ \\\n')
                else:
                    script_file_id.write(f'                --aligner={alignment_tool} \\\n')
                if i == 0:
                    script_file_id.write( '                --seqs=$REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
                else:
                    script_file_id.write(f'                --seqs=${previus_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
                script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                script_file_id.write( '                --relationships2=NONE \\\n')
                if current_code == 'nt':
                    script_file_id.write(f'                --annotation=$NT_VIRIDIPLANTAE_ANNOTATION_FILE \\\n')
                    script_file_id.write(f'                --annotation2=$NT_CONTAMINATION_ANNOTATION_FILE \\\n')
                else:
                    script_file_id.write(f'                --annotation=${current_code.upper()}_ANNOTATION_FILE \\\n')
                    script_file_id.write(f'                --annotation2=NONE \\\n')
                script_file_id.write(f'                --nonann=${current_code.upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error annotate-sequences.py $RC; fi\n')
                script_file_id.write( '        echo "Annotation is done."\n')
                if len(database_list) == 1:
                    if current_code == 'nt':
                        script_file_id.write( '        ANNOTATION_FILE_TMP=$NT_VIRIDIPLANTAE_ANNOTATION_FILE.tmp\n')
                        script_file_id.write( '        echo "Deleting the header record of $NT_VIRIDIPLANTAE_ANNOTATION_FILE ..."\n')
                        script_file_id.write( '        /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write( '            tail -n +2 $NT_VIRIDIPLANTAE_ANNOTATION_FILE > $ANNOTATION_FILE_TMP\n')
                    else:
                        script_file_id.write(f'        ANNOTATION_FILE_TMP=${current_code.upper()}_ANNOTATION_FILE.tmp\n')
                        script_file_id.write(f'        echo "Deleting the header record of ${current_code.upper()}_ANNOTATION_FILE ..."\n')
                        script_file_id.write( '        /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write(f'            tail -n +2 ${current_code.upper()}_ANNOTATION_FILE > $ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                    script_file_id.write( '        echo "Record is deleted."\n')
                    script_file_id.write( '        echo "Creating plant annotation file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            merge-annotation-files.py \\\n')
                    script_file_id.write( '            --file1=$ANNOTATION_FILE_TMP \\\n')
                    script_file_id.write(f'            --type1={database_type_dict[database_list[0]]} \\\n')
                    script_file_id.write( '            --file2=NONE \\\n')
                    script_file_id.write( '            --type2=NONE \\\n')
                    script_file_id.write( '            --operation=SAVE1 \\\n')
                    script_file_id.write( '            --mfile=$PLANT_ANNOTATION_FILE \\\n')
                    script_file_id.write( '            --header=Y \\\n')
                    script_file_id.write( '            --verbose=N \\\n')
                    script_file_id.write( '            --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                    script_file_id.write( '        echo "File is created."\n')
                    script_file_id.write( '        echo "Deleting temporal file  ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            rm $ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                    script_file_id.write( '        echo "File is deleted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            if len(database_list) > 1:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_alignment_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_alignment_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "MERGER OF ALIGNMENT FILES"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Merging alignment files ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            merge-xml-files.py \\\n')
                blast_xml_list = []
                for database_code in database_list:
                    blast_xml_list.append(f'${database_code.upper()}_BLAST_XML')
                script_file_id.write(f'                --list={",".join(blast_xml_list)} \\\n')
                script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                script_file_id.write( '                --relationships2=NONE \\\n')
                script_file_id.write( '                --mfile=$MERGED_BLAST_XML \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-xml-files.py $RC; fi\n')
                script_file_id.write( '        echo "Files are merged."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_annotation_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_annotation_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "MERGER OF PLANT ANNOTATION FILES"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                for database_code in database_list2:
                    script_file_id.write(f'        {database_code.upper()}_ANNOTATION_FILE_TMP=${database_code.upper()}_ANNOTATION_FILE".tmp"\n')
                    script_file_id.write(f'        {database_code.upper()}_ANNOTATION_FILE_SORTED=${database_code.upper()}_ANNOTATION_FILE".sorted"\n')
                tmp_file_list = []
                for i in range(len(database_list2) - 2):
                    script_file_id.write(f'        MERGED_ANNOTATION_FILE_TMP{i + 1}=$MERGED_ANNOTATION_FILE".tmp{i + 1}"\n')
                    tmp_file_list.append(f'$MERGED_ANNOTATION_FILE_TMP{i + 1}')
                script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            tail -n +2 ${database_list2[0].upper()}_ANNOTATION_FILE > ${database_list2[0].upper()}_ANNOTATION_FILE_TMP\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write( '        echo "Record is deleted."\n')
                script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[0].upper()}_ANNOTATION_FILE_TMP > ${database_list2[0].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '        echo "Records are sorted."\n')
                script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            tail -n +2 ${database_list2[1].upper()}_ANNOTATION_FILE > ${database_list2[1].upper()}_ANNOTATION_FILE_TMP\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        echo "Record is deleted."\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[1].upper()}_ANNOTATION_FILE_TMP > ${database_list2[1].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '        echo "Records are sorted."\n')
                script_file_id.write(f'        echo "Merging `basename ${database_list2[0].upper()}_ANNOTATION_FILE` and `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            merge-annotation-files.py \\\n')
                script_file_id.write(f'                --file1=${database_list2[0].upper()}_ANNOTATION_FILE_SORTED \\\n')
                script_file_id.write(f'                --type1={database_type_dict[database_list[0]]} \\\n')
                script_file_id.write(f'                --file2=${database_list2[1].upper()}_ANNOTATION_FILE_SORTED \\\n')
                script_file_id.write(f'                --type2={database_type_dict[database_list[1]]} \\\n')
                script_file_id.write( '                --operation=1AND2 \\\n')
                if len(database_list2) > 2:
                    script_file_id.write( '                --mfile=$MERGED_ANNOTATION_FILE_TMP1 \\\n')
                    script_file_id.write( '                --header=N \\\n')
                else:
                    script_file_id.write( '                --mfile=$PLANT_ANNOTATION_FILE \\\n')
                    script_file_id.write( '                --header=Y \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                script_file_id.write( '        echo "Files are merged."\n')
                script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm ${database_list2[0].upper()}_ANNOTATION_FILE_TMP ${database_list2[0].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm ${database_list2[1].upper()}_ANNOTATION_FILE_TMP ${database_list2[1].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                for i in range(2, len(database_list2)):
                    script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            tail -n +2 ${database_list2[i].upper()}_ANNOTATION_FILE > ${database_list2[i].upper()}_ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        echo "Record is deleted."\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                    script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[i].upper()}_ANNOTATION_FILE_TMP > ${database_list2[i].upper()}_ANNOTATION_FILE_SORTED\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                    script_file_id.write( '        echo "Records are sorted."\n')
                    script_file_id.write(f'        echo "Adding annotation of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` to the merged annotation file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            merge-annotation-files.py \\\n')
                    script_file_id.write(f'                --file1=$MERGED_ANNOTATION_FILE_TMP{i - 1} \\\n')
                    script_file_id.write( '                --type1=MERGER \\\n')
                    script_file_id.write(f'                --file2=${database_list2[i].upper()}_ANNOTATION_FILE_SORTED \\\n')
                    script_file_id.write(f'                --type2={database_type_dict[database_list[i]]} \\\n')
                    script_file_id.write( '                --operation=1AND2 \\\n')
                    if i < len(database_list2) - 1:
                        script_file_id.write(f'                --mfile=$MERGED_ANNOTATION_FILE_TMP{i} \\\n')
                        script_file_id.write( '                --header=N \\\n')
                    else:
                        script_file_id.write( '                --mfile=$PLANT_ANNOTATION_FILE \\\n')
                        script_file_id.write( '                --header=Y \\\n')
                    script_file_id.write( '                --verbose=N \\\n')
                    script_file_id.write( '                --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                    script_file_id.write( '        echo "Files are merged."\n')
                    script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            rm ${database_list2[i].upper()}_ANNOTATION_FILE_TMP ${database_list2[i].upper()}_ANNOTATION_FILE_SORTED\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                    script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write( '        echo "Deleting temporal annotation files ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm {" ".join(tmp_file_list)}\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function split_merged_plant_annotation_file\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/split_merged_plant_annotation_file.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "SPLIT OF MERGED PLANT ANNOTATION FILE"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Splitting merged file `basename $PLANT_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            split-annotation-file.py \\\n')
                script_file_id.write( '                --annotation=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '                --type=MERGER \\\n')
                script_file_id.write( '                --header=Y \\\n')
                script_file_id.write( '                --rnum=250000 \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error split-annotation-file.py $RC; fi\n')
                script_file_id.write( '        echo "File is splitted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function purge_transcriptome\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/purge_transcriptome.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "PURGE TRANSCRIPTOME REMOVING NON-ANNOTATED TRANSCRIPTS"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Purging transcriptome files ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            merge-fasta-files.py \\\n')
            script_file_id.write( '                --file1=$REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write(f'                --file2=${database_list[len(database_list) - 1].upper()}_NON_ANNOTATED_TRANSCRIPT_FILE \\\n')
            script_file_id.write( '                --mfile=$PURGED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --operation=1LESS2 \\\n')
            script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-fasta-files.py $RC; fi\n')
            script_file_id.write( '        echo "Transcriptome is purged."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_annotation_stats\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/calculate_annotation_stats.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "CALCULATE ANNOTATION STATISTICS"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Calculating stats ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            calculate-annotation-stats.py \\\n')
            script_file_id.write( '                --db=$TOA_DB \\\n')
            script_file_id.write( '                --transcriptome=$TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --peptides=NONE \\\n')
            script_file_id.write(f'                --dslist={",".join(database_list)} \\\n')
            script_file_id.write(f'                --nonannlist={",".join(non_annotation_file_list)} \\\n')
            if len(database_list) > 1:
                script_file_id.write( '                --annotation=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '                --type=MERGER \\\n')
            else:
                script_file_id.write(f'                --annotation=${database_list[0].upper()}_ANNOTATION_FILE \\\n')
                script_file_id.write(f'                --type={database_type_dict[database_list[0]]} \\\n')
            script_file_id.write( '                --stats=$ANNOTATION_STATS_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error calculate-annotation-stats.py $RC; fi\n')
            script_file_id.write( '        echo "Stats are calculated."\n')
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
            process_name = f'TOA - {xlib.get_toa_process_pipeline_nucleotide_name()}'
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
            script_file_id.write( '\n')
            script_file_id.write( '# re-identify sequences of the transcriptome file\n')
            script_file_id.write( 'reidentify_transcript_sequences\n')
            script_file_id.write( '\n')
            for i in range(len(database_list)):
                current_code = database_list[i]
                previus_code = database_list[i - 1] if i > 0 else ''
                if i == 0:
                    script_file_id.write(f'# complete transcriptome -> {current_code}\n')
                else:
                    script_file_id.write(f'# transcripts not annotated with {previus_code} -> {current_code}\n')
                script_file_id.write(f'align_transcripts_{current_code}_proteome\n')
                script_file_id.write(f'load_alignment_{current_code}_proteome\n')
                script_file_id.write(f'annotate_transcripts_{current_code}\n')
                script_file_id.write( '\n')
            if len(database_list) > 1:
                script_file_id.write( '# merged files\n')
                script_file_id.write( 'merge_alignment_files\n')
                script_file_id.write( 'merge_annotation_files\n')
                script_file_id.write( '# -- split_merged_plant_annotation_file\n')
                script_file_id.write( '\n')
            script_file_id.write( '# transcriptome with plant transcripts\n')
            script_file_id.write( 'purge_transcriptome\n')
            script_file_id.write( '\n')
            script_file_id.write( '# annotation statistics\n')
            script_file_id.write( 'calculate_annotation_stats\n')
            script_file_id.write( '\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nucleotide_pipeline_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nucleotide_pipeline_starter(current_run_dir):
    '''
    Build the starter of the script to process a nucleotide pipeline.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nucleotide_pipeline_starter())):
            os.makedirs(os.path.dirname(get_nucleotide_pipeline_starter()))
        with open(get_nucleotide_pipeline_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_nucleotide_pipeline_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_nucleotide_pipeline_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nucleotide_pipeline_script():
    '''
    Get the script path to process a nucleotide pipeline.
    '''

    # assign the script path
    nucleotide_pipeline_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_pipeline_nucleotide_code()}-process.sh'

    # return the script path
    return nucleotide_pipeline_script

#-------------------------------------------------------------------------------

def get_nucleotide_pipeline_starter():
    '''
    Get the starter path to process a nucleotide pipeline.
    '''

    # assign the starter path
    nucleotide_pipeline_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_pipeline_nucleotide_code()}-process-starter.sh'

    # return the starter path
    return nucleotide_pipeline_starter

#-------------------------------------------------------------------------------

def build_aminoacid_pipeline_script(cluster_name, current_run_dir):
    '''
    Build the script to process a amino acid pipeline.
    '''


    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # get the pipeline option dictionary
    pipeline_option_dict = xlib.get_option_dict(get_aminoacid_pipeline_config_file())

    # get the options
    assembly_origin = pipeline_option_dict['identification']['assembly_origin']
    experiment_id = pipeline_option_dict['identification']['experiment_id']
    assembly_software = pipeline_option_dict['identification']['assembly_software']
    assembly_dataset_id = pipeline_option_dict['identification']['assembly_dataset_id']
    assembly_type = pipeline_option_dict['identification']['assembly_type']
    reference_dataset_id = pipeline_option_dict['identification']['reference_dataset_id']
    transcriptome_file = pipeline_option_dict['identification']['transcriptome_file']
    alignment_tool = pipeline_option_dict['pipeline parameters']['alignment_tool']
    threads = pipeline_option_dict['pipeline parameters']['threads']
    blastplus_evalue = pipeline_option_dict['BLAST+ parameters']['evalue']
    blastplus_max_target_seqs = pipeline_option_dict['BLAST+ parameters']['max_target_seqs']
    blastplus_max_hsps = pipeline_option_dict['BLAST+ parameters']['max_hsps']
    blastplus_qcov_hsp_perc = pipeline_option_dict['BLAST+ parameters']['qcov_hsp_perc']
    blastplus_other_parameters_blastp = pipeline_option_dict['BLAST+ parameters']['other_parameters_blastp']
    diamond_evalue = pipeline_option_dict['DIAMOND parameters']['evalue']
    diamond_max_target_seqs = pipeline_option_dict['DIAMOND parameters']['max-target-seqs']
    diamond_max_hsps = pipeline_option_dict['DIAMOND parameters']['max-hsps']
    diamond_other_parameters_blastp = pipeline_option_dict['DIAMOND parameters']['other_parameters_blastp']

    # get the all selected database list
    database_list = get_selected_database_list(xlib.get_toa_process_pipeline_aminoacid_code())

    # get the database type dictionary
    database_type_dict = {}
    for i in range(len(database_list)):
        if database_list[i] in ['gymno_01', 'dicots_04', 'monocots_04']:
            database_type_dict[database_list[i]] = 'PLAZA'
        elif database_list[i] == 'refseq_plant':
            database_type_dict[database_list[i]] = 'REFSEQ'
        elif database_list[i] in ['nr']:
            database_type_dict[database_list[i]] = 'NR'

    # get the database list for the merger of plant annotation files
    database_list2 = database_list.copy()
    for i in range(len(database_list2)):
        if database_list2[i] == 'nr':
            database_list2[i] = 'nr_viridiplantae'

    # set the transcriptome file path
    if OK:
        if assembly_origin == 'NGSCLOUD':
            if assembly_software == xlib.get_soapdenovotrans_code():
                if assembly_type == 'CONTIGS':
                    transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.contig'
                elif  assembly_type == 'SCAFFOLDS':
                    transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.scafSeq'
            elif assembly_software == xlib.get_transabyss_code():
                transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/transabyss-final.fa'
            elif assembly_software == xlib.get_trinity_code():
                transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity.fasta'
            elif assembly_software == xlib.get_cd_hit_est_code():
                transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/clustered-transcriptome.fasta'
            elif assembly_software == xlib.get_transcript_filter_code():
                transcriptome_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/filtered-transcriptome.fasta'
        elif assembly_origin == 'EXTERNAL':
            transcriptome_file = xlib.get_cluster_reference_file(reference_dataset_id, transcriptome_file)

    # get the non annotation file list
    non_annotation_file_list = []
    for database in database_list:
        non_annotation_file_list.append(f'${database.upper()}_NON_ANNOTATED_PEPTIDE_FILE')

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_aminoacid_pipeline_script())):
            os.makedirs(os.path.dirname(get_aminoacid_pipeline_script()))
        with open(get_aminoacid_pipeline_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '# transcriptome file\n')
            script_file_id.write(f'TRANSCRIPTOME_FILE={transcriptome_file}\n')
            script_file_id.write( '\n')
            script_file_id.write( '# pipeline parameters\n')
            script_file_id.write(f'THREADS={threads}\n')
            script_file_id.write( '\n')
            script_file_id.write( '# BLAST+ parameters\n')
            script_file_id.write(f'BLASTPLUS_EVALUE={blastplus_evalue}\n')
            script_file_id.write(f'BLASTPLUS_MAX_TARGET_SEQS={blastplus_max_target_seqs}\n')
            script_file_id.write(f'BLASTPLUS_MAX_HSPS={blastplus_max_hsps}\n')
            script_file_id.write(f'BLASTPLUS_QCOV_HSP_PERC={blastplus_qcov_hsp_perc}\n')
            script_file_id.write(f'BLASTPLUS_OTHER_PARAMETERS_BLASTP="{blastplus_other_parameters_blastp}"\n')
            script_file_id.write( '\n')
            script_file_id.write( '# DIAMOND parameters\n')
            script_file_id.write(f'DIAMOND_EVALUE={diamond_evalue}\n')
            script_file_id.write(f'DIAMOND_MAX_TARGET_SEQS={diamond_max_target_seqs}\n')
            script_file_id.write(f'DIAMOND_MAX_HSPS={diamond_max_hsps}\n')
            script_file_id.write(f'DIAMOND_OTHER_PARAMETERS_BLASTP="{diamond_other_parameters_blastp}"\n')
            script_file_id.write( '\n')
            script_file_id.write( '# output directory\n')
            script_file_id.write(f'OUTPUT_DIR={current_run_dir}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                records = toa_config_file_id.readlines()
                for record in records:
                    script_file_id.write(record)
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
            script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
            script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'mkdir --parents $STATS_DIR\n')
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
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "TRANSCRIPTOME FILE: $TRANSCRIPTOME_FILE"\n')
            script_file_id.write( '    echo "$SEP"\n')
            database_list_text = ','.join(database_list)
            script_file_id.write(f'    echo "ALIGNMENT DATASETS: {database_list_text}"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "ALIGNMENT TOOLS (blastp): {alignment_tool}"\n')
            script_file_id.write( '    echo "THREADS: $THREADS"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "BLAST+ EVALUE: $BLASTPLUS_EVALUE"\n')
            script_file_id.write( '    echo "BLAST+ MAX_TARGET_SEQS: $BLASTPLUS_MAX_TARGET_SEQS"\n')
            script_file_id.write( '    echo "BLAST+ MAX_HSPS: $BLASTPLUS_MAX_HSPS"\n')
            script_file_id.write( '    echo "BLAST+ QCOV_HSP_PERC: $BLASTPLUS_QCOV_HSP_PERC"\n')
            script_file_id.write( '    echo "BLAST+ BLASTP OTHER PARAMETERS: $BLASTPLUS_OTHER_PARAMETERS_BLASTP"\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "DIAMOND EVALUE: $DIAMOND_EVALUE"\n')
            script_file_id.write( '    echo "DIAMOND MAX_TARGET_SEQS: $DIAMOND_MAX_TARGET_SEQS"\n')
            script_file_id.write( '    echo "DIAMOND MAX_HSPS: $DIAMOND_MAX_HSPS"\n')
            script_file_id.write( '    echo "DIAMOND BLASTP OTHER PARAMETERS: $DIAMOND_OTHER_PARAMETERS_BLASTP"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function reidentify_transcript_sequences\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/reidentify_transcript_sequences.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "RE-INDENTIFY SEQUENCES OF THE TRANSCRIPTOME FILE"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Re-identifing sequences ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            reid-fasta-file.py \\\n')
            script_file_id.write( '                --fasta=$TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --type=NT \\\n')
            script_file_id.write( '                --out=$REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error reid-fasta-file.py $RC; fi\n')
            script_file_id.write( '        echo "Sequences are re-identified."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function extract_orfs\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/extract_orfs.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "EXTRACT THE LONG OPEN READING FRAMES"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        source activate transdecoder\n')
            script_file_id.write( '        echo "Extracting ORFs ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            TransDecoder.LongOrfs \\\n')
            script_file_id.write( '                -t $REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                -m 100 \\\n')
            script_file_id.write( '                --output_dir $TRANSDECODER_OUTPUT_DIR\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error TransDecoder.LongOrfs $RC; fi\n')
            script_file_id.write( '        echo "ORFs are extracted."\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function predict_coding_regions\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/predict_coding_regions.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "PREDICT THE LIKELY CODING REGIONS"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        source activate transdecoder\n')
            script_file_id.write( '        echo "Predicting codign regions ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            TransDecoder.Predict \\\n')
            script_file_id.write( '                -t $REIDENTIFIED_TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --output_dir $TRANSDECODER_OUTPUT_DIR\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error TransDecoder.Predict $RC; fi\n')
            script_file_id.write( '        echo "Coding regions are predicted."\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function reidentify_peptide_sequences\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/reidentify_peptide_sequences.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "RE-INDENTIFY SEQUENCES OF THE PEPTIDE FILE"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Re-identifing sequences ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            reid-fasta-file.py \\\n')
            script_file_id.write( '                --fasta=$PEPTIDE_FILE \\\n')
            script_file_id.write( '                --type=AA \\\n')
            script_file_id.write( '                --out=$REIDENTIFIED_PEPTIDE_FILE \\\n')
            script_file_id.write( '                --relationships=$TOA_TRANSDECODER_RELATIONSHIP_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error reid-fasta-file.py $RC; fi\n')
            script_file_id.write( '        echo "Sequences are re-identified."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            for i in range(len(database_list)):
                current_code = database_list[i]
                previus_code = database_list[i - 1] if i > 0 else ''
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function align_peptides_{current_code}_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/align_peptides_{current_code}_proteome.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "ALIGNMENT OF PEPTIDES TO {current_code.upper()} PROTEOME"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        source activate blast\n')
                script_file_id.write( '        echo "Aligning peptides ..."\n')
                if alignment_tool == xlib.get_blastplus_name():
                    script_file_id.write(f'        export BLASTDB=${current_code.upper()}_BLASTPLUS_DB_DIR\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            blastp \\\n')
                    script_file_id.write( '                -num_threads $THREADS \\\n')
                    script_file_id.write(f'                -db ${current_code.upper()}_BLASTPLUS_DB_NAME \\\n')
                    if i == 0:
                        script_file_id.write( '                -query $REIDENTIFIED_PEPTIDE_FILE \\\n')
                    else:
                        script_file_id.write(f'                -query ${previus_code.upper()}_NON_ANNOTATED_PEPTIDE_FILE \\\n')
                    script_file_id.write( '                -evalue $BLASTPLUS_EVALUE \\\n')
                    script_file_id.write( '                -max_target_seqs $BLASTPLUS_MAX_TARGET_SEQS \\\n')
                    script_file_id.write( '                -max_hsps $BLASTPLUS_MAX_HSPS \\\n')
                    script_file_id.write( '                -qcov_hsp_perc $BLASTPLUS_QCOV_HSP_PERC \\\n')
                    script_file_id.write( '                -outfmt 5 \\\n')
                    if blastplus_other_parameters_blastp.upper() != 'NONE':
                        parameter_list = [x.strip() for x in blastplus_other_parameters_blastp.split(';')]
                        for parameter in parameter_list:
                            if parameter.find('=') > 0:
                                pattern = r'^--(.+)=(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                parameter_value = mo.group(2).strip()
                                script_file_id.write(f'                -{parameter_name} {parameter_value} \\\n')
                            else:
                                pattern = r'^--(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                script_file_id.write(f'                -{parameter_name} \\\n')
                    script_file_id.write(f'                -out ${current_code.upper()}_BLAST_XML\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error blastp $RC; fi\n')
                    script_file_id.write( '        echo "Alignment is done."\n')
                    script_file_id.write( '        conda deactivate\n')
                elif alignment_tool == xlib.get_diamond_name():
                    if i == 0:
                        script_file_id.write( '        if [[ -s $REIDENTIFIED_PEPTIDE_FILE ]]; then\n')
                    else:
                        script_file_id.write(f'        if [[ -s ${previus_code.upper()}_NON_ANNOTATED_PEPTIDE_FILE ]]; then\n')
                    script_file_id.write( '            source activate diamond\n')
                    script_file_id.write( '            /usr/bin/time \\\n')
                    script_file_id.write(f'                --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '                diamond blastp \\\n')
                    script_file_id.write( '                    --threads $THREADS \\\n')
                    script_file_id.write(f'                    --db ${current_code.upper()}_DIAMOND_DB_FILE \\\n')
                    if i == 0:
                        script_file_id.write( '                    --query $REIDENTIFIED_PEPTIDE_FILE \\\n')
                    else:
                        script_file_id.write(f'                    --query ${previus_code.upper()}_NON_ANNOTATED_PEPTIDE_FILE \\\n')
                    script_file_id.write( '                    --evalue $DIAMOND_EVALUE \\\n')
                    script_file_id.write( '                    --max-target-seqs $DIAMOND_MAX_TARGET_SEQS \\\n')
                    script_file_id.write( '                    --max-hsps $DIAMOND_MAX_HSPS \\\n')
                    if diamond_other_parameters_blastp.upper() != 'NONE':
                        parameter_list = [x.strip() for x in diamond_other_parameters_blastp.split(';')]
                        for parameter in parameter_list:
                            if parameter.find('=') > 0:
                                pattern = r'^--(.+)=(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                parameter_value = mo.group(2).strip()
                                script_file_id.write(f'                --{parameter_name} {parameter_value} \\\n')
                            else:
                                pattern = r'^--(.+)$'
                                mo = re.search(pattern, parameter)
                                parameter_name = mo.group(1).strip()
                                script_file_id.write(f'                --{parameter_name} \\\n')
                    script_file_id.write( '                    --outfmt 5 \\\n')
                    script_file_id.write(f'                    --out ${current_code.upper()}_BLAST_XML\n')
                    script_file_id.write( '            RC=$?\n')
                    script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error diamond-blastp $RC; fi\n')
                    script_file_id.write( '        else\n')
                    script_file_id.write(f'            touch ${current_code.upper()}_BLAST_XML\n')
                    script_file_id.write( '        fi\n')
                    script_file_id.write( '        echo "Alignment is done."\n')
                    script_file_id.write( '        conda deactivate\n')
                if len(database_list) == 1:
                    script_file_id.write( '        echo "Restoring sequence identifications in alignment file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            restore-ids.py \\\n')
                    script_file_id.write(f'                --in=${current_code.upper()}_BLAST_XML \\\n')
                    script_file_id.write( '                --format=XML \\\n')
                    script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                    script_file_id.write( '                --relationships2=$TOA_TRANSDECODER_RELATIONSHIP_FILE \\\n')
                    script_file_id.write( '                --out=$MERGED_BLAST_XML \\\n')
                    script_file_id.write( '                --verbose=N \\\n')
                    script_file_id.write( '                --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi\n')
                    script_file_id.write( '        echo "Identifications are restored."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function load_alignment_{current_code}_proteome\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/load_alignment_{current_code}_proteome.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "LOAD OF PEPTIDE ALIGNMENT TO {current_code.upper()} PROTEOME INTO TOA DATABASE"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Loading alignmnet data ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            load-blast-data.py \\\n')
                script_file_id.write( '                --db=$TOA_DB \\\n')
                script_file_id.write(f'                --dataset={current_code} \\\n')
                script_file_id.write( '                --format=5 \\\n')
                script_file_id.write(f'                --blast=${current_code.upper()}_BLAST_XML \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error load-blast-data.py $RC; fi\n')
                script_file_id.write( '        echo "Data are loaded."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'function annotate_peptides_{current_code}\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write(f'    STEP_STATUS=$STATUS_DIR/annotate_peptides_{current_code}.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "ANNOTATION OF PEPTIDES WITH {current_code.upper()}"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Annotating peptides ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            annotate-sequences.py \\\n')
                script_file_id.write( '                --db=$TOA_DB \\\n')
                script_file_id.write(f'                --dataset={current_code} \\\n')
                script_file_id.write(f'                --aligner={alignment_tool} \\\n')
                if i == 0:
                    script_file_id.write( '                --seqs=$REIDENTIFIED_PEPTIDE_FILE \\\n')
                else:
                    script_file_id.write(f'                --seqs=${previus_code.upper()}_NON_ANNOTATED_PEPTIDE_FILE \\\n')
                script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                script_file_id.write( '                --relationships2=$TOA_TRANSDECODER_RELATIONSHIP_FILE \\\n')
                if current_code == 'nr':
                    script_file_id.write(f'                --annotation=$NR_VIRIDIPLANTAE_ANNOTATION_FILE \\\n')
                    script_file_id.write(f'                --annotation2=$NR_CONTAMINATION_ANNOTATION_FILE \\\n')
                else:
                    script_file_id.write(f'                --annotation=${current_code.upper()}_ANNOTATION_FILE \\\n')
                    script_file_id.write(f'                --annotation2=NONE \\\n')
                script_file_id.write(f'                --nonann=${current_code.upper()}_NON_ANNOTATED_PEPTIDE_FILE \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error annotate-sequences.py $RC; fi\n')
                script_file_id.write( '        echo "Annotation is done."\n')
                if len(database_list) == 1:
                    if current_code == 'nr':
                        script_file_id.write( '        ANNOTATION_FILE_TMP=$NR_VIRIDIPLANTAE_ANNOTATION_FILE.tmp\n')
                        script_file_id.write( '        echo "Deleting the header record of $NR_VIRIDIPLANTAE_ANNOTATION_FILE ..."\n')
                        script_file_id.write( '        /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write( '            tail -n +2 $NR_VIRIDIPLANTAE_ANNOTATION_FILE > $ANNOTATION_FILE_TMP\n')
                    else:
                        script_file_id.write(f'        ANNOTATION_FILE_TMP=${current_code.upper()}_ANNOTATION_FILE.tmp\n')
                        script_file_id.write(f'        echo "Deleting the header record of ${current_code.upper()}_ANNOTATION_FILE ..."\n')
                        script_file_id.write( '        /usr/bin/time \\\n')
                        script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                        script_file_id.write(f'            tail -n +2 ${current_code.upper()}_ANNOTATION_FILE > $ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                    script_file_id.write( '        echo "Record is deleted."\n')
                    script_file_id.write( '        echo "Creating plant annotation file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            merge-annotation-files.py \\\n')
                    script_file_id.write( '            --file1=$ANNOTATION_FILE_TMP \\\n')
                    script_file_id.write(f'            --type1={database_type_dict[database_list[0]]} \\\n')
                    script_file_id.write( '            --file2=NONE \\\n')
                    script_file_id.write( '            --type2=NONE \\\n')
                    script_file_id.write( '            --operation=SAVE1 \\\n')
                    script_file_id.write( '            --mfile=$PLANT_ANNOTATION_FILE \\\n')
                    script_file_id.write( '            --header=Y \\\n')
                    script_file_id.write( '            --verbose=N \\\n')
                    script_file_id.write( '            --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                    script_file_id.write( '        echo "File is created."\n')
                    script_file_id.write( '        echo "Deleting temporal file  ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            rm $ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                    script_file_id.write( '        echo "File is deleted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            if len(database_list) > 1:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_alignment_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_alignment_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "MERGER OF ALIGNMENT FILES"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Merging alignment files ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            merge-xml-files.py \\\n')
                blast_xml_list = []
                for database_code in database_list:
                    blast_xml_list.append(f'${database_code.upper()}_BLAST_XML')
                script_file_id.write(f'                --list={",".join(blast_xml_list)} \\\n')
                script_file_id.write( '                --relationships=$TOA_TRANSCRIPTOME_RELATIONSHIP_FILE \\\n')
                script_file_id.write( '                --relationships2=$TOA_TRANSDECODER_RELATIONSHIP_FILE \\\n')
                script_file_id.write( '                --mfile=$MERGED_BLAST_XML \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-xml-files.py $RC; fi\n')
                script_file_id.write( '        echo "Files are merged."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_annotation_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_annotation_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "MERGER OF PLANT ANNOTATION FILES"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                for database_code in database_list2:
                    script_file_id.write(f'        {database_code.upper()}_ANNOTATION_FILE_TMP=${database_code.upper()}_ANNOTATION_FILE".tmp"\n')
                    script_file_id.write(f'        {database_code.upper()}_ANNOTATION_FILE_SORTED=${database_code.upper()}_ANNOTATION_FILE".sorted"\n')
                tmp_file_list = []
                for i in range(len(database_list2) - 2):
                    script_file_id.write(f'        MERGED_ANNOTATION_FILE_TMP{i + 1}=$MERGED_ANNOTATION_FILE".tmp{i + 1}"\n')
                    tmp_file_list.append(f'$MERGED_ANNOTATION_FILE_TMP{i + 1}')
                script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            tail -n +2 ${database_list2[0].upper()}_ANNOTATION_FILE > ${database_list2[0].upper()}_ANNOTATION_FILE_TMP\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write( '        echo "Record is deleted."\n')
                script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[0].upper()}_ANNOTATION_FILE_TMP > ${database_list2[0].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '        echo "Records are sorted."\n')
                script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            tail -n +2 ${database_list2[1].upper()}_ANNOTATION_FILE > ${database_list2[1].upper()}_ANNOTATION_FILE_TMP\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        echo "Record is deleted."\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[1].upper()}_ANNOTATION_FILE_TMP > ${database_list2[1].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '        echo "Records are sorted."\n')
                script_file_id.write(f'        echo "Merging `basename ${database_list2[0].upper()}_ANNOTATION_FILE` and `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            merge-annotation-files.py \\\n')
                script_file_id.write(f'                --file1=${database_list2[0].upper()}_ANNOTATION_FILE_SORTED \\\n')
                script_file_id.write(f'                --type1={database_type_dict[database_list[0]]} \\\n')
                script_file_id.write(f'                --file2=${database_list2[1].upper()}_ANNOTATION_FILE_SORTED \\\n')
                script_file_id.write(f'                --type2={database_type_dict[database_list[1]]} \\\n')
                script_file_id.write( '                --operation=1AND2 \\\n')
                if len(database_list2) > 2:
                    script_file_id.write( '                --mfile=$MERGED_ANNOTATION_FILE_TMP1 \\\n')
                    script_file_id.write( '                --header=N \\\n')
                else:
                    script_file_id.write( '                --mfile=$PLANT_ANNOTATION_FILE \\\n')
                    script_file_id.write( '                --header=Y \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                script_file_id.write( '        echo "Files are merged."\n')
                script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[0].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm ${database_list2[0].upper()}_ANNOTATION_FILE_TMP ${database_list2[0].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[1].upper()}_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm ${database_list2[1].upper()}_ANNOTATION_FILE_TMP ${database_list2[1].upper()}_ANNOTATION_FILE_SORTED\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                for i in range(2, len(database_list2)):
                    script_file_id.write(f'        echo "Deleting the header record of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            tail -n +2 ${database_list2[i].upper()}_ANNOTATION_FILE > ${database_list2[i].upper()}_ANNOTATION_FILE_TMP\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        echo "Record is deleted."\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                    script_file_id.write(f'        echo "Sorting data records of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            sort --field-separator=";" --key=2,5 < ${database_list2[i].upper()}_ANNOTATION_FILE_TMP > ${database_list2[i].upper()}_ANNOTATION_FILE_SORTED\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                    script_file_id.write( '        echo "Records are sorted."\n')
                    script_file_id.write(f'        echo "Adding annotation of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` to the merged annotation file ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write( '            merge-annotation-files.py \\\n')
                    script_file_id.write(f'                --file1=$MERGED_ANNOTATION_FILE_TMP{i - 1} \\\n')
                    script_file_id.write( '                --type1=MERGER \\\n')
                    script_file_id.write(f'                --file2=${database_list2[i].upper()}_ANNOTATION_FILE_SORTED \\\n')
                    script_file_id.write(f'                --type2={database_type_dict[database_list[i]]} \\\n')
                    script_file_id.write( '                --operation=1AND2 \\\n')
                    if i < len(database_list2) - 1:
                        script_file_id.write(f'                --mfile=$MERGED_ANNOTATION_FILE_TMP{i} \\\n')
                        script_file_id.write( '                --header=N \\\n')
                    else:
                        script_file_id.write( '                --mfile=$PLANT_ANNOTATION_FILE \\\n')
                        script_file_id.write( '                --header=Y \\\n')
                    script_file_id.write( '                --verbose=N \\\n')
                    script_file_id.write( '                --trace=N\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                    script_file_id.write( '        echo "Files are merged."\n')
                    script_file_id.write(f'        echo "Deleting temporal files of `basename ${database_list2[i].upper()}_ANNOTATION_FILE` ..."\n')
                    script_file_id.write( '        /usr/bin/time \\\n')
                    script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                    script_file_id.write(f'            rm ${database_list2[i].upper()}_ANNOTATION_FILE_TMP ${database_list2[i].upper()}_ANNOTATION_FILE_SORTED\n')
                    script_file_id.write( '        RC=$?\n')
                    script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                    script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write( '        echo "Deleting temporal annotation files ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write(f'            rm {" ".join(tmp_file_list)}\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '        echo "Files are deleted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function split_merged_plant_annotation_file\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/split_merged_plant_annotation_file.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "SPLIT OF MERGED PLANT ANNOTATION FILE"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Splitting merged file `basename $PLANT_ANNOTATION_FILE` ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            split-annotation-file.py \\\n')
                script_file_id.write( '                --annotation=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '                --type=MERGER \\\n')
                script_file_id.write( '                --header=Y \\\n')
                script_file_id.write( '                --rnum=250000 \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error split-annotation-file.py $RC; fi\n')
                script_file_id.write( '        echo "File is splitted."\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function calculate_annotation_stats\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $OUTPUT_DIR\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/calculate_annotation_stats.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "CALCULATE ANNOTATION STATISTICS"\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "Calculating stats ..."\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            calculate-annotation-stats.py \\\n')
            script_file_id.write( '                --db=$TOA_DB \\\n')
            script_file_id.write( '                --transcriptome=$TRANSCRIPTOME_FILE \\\n')
            script_file_id.write( '                --peptides=$PEPTIDE_FILE \\\n')
            script_file_id.write(f'                --dslist={",".join(database_list)} \\\n')
            script_file_id.write(f'                --nonannlist={",".join(non_annotation_file_list)} \\\n')
            if len(database_list) > 1:
                script_file_id.write( '                --annotation=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '                --type=MERGER \\\n')
            else:
                script_file_id.write(f'                --annotation=${database_list[0].upper()}_ANNOTATION_FILE \\\n')
                script_file_id.write(f'                --type={database_type_dict[database_list[0]]} \\\n')
            script_file_id.write( '                --stats=$ANNOTATION_STATS_FILE \\\n')
            script_file_id.write( '                --verbose=N \\\n')
            script_file_id.write( '                --trace=N\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error calculate-annotation-stats.py $RC; fi\n')
            script_file_id.write( '        echo "Stats are calculated."\n')
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
            process_name = f'TOA - {xlib.get_toa_process_pipeline_aminoacid_name()}'
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
            script_file_id.write( '\n')
            script_file_id.write( '# re-identify sequences of the transcriptome file\n')
            script_file_id.write( 'reidentify_transcript_sequences\n')
            script_file_id.write( '\n')
            script_file_id.write( '# extract the long open reading frames and predict coding regions\n')
            script_file_id.write( 'extract_orfs\n')
            script_file_id.write( 'predict_coding_regions\n')
            script_file_id.write( '\n')
            script_file_id.write( '# re-identify sequences of the peptide file\n')
            script_file_id.write( 'reidentify_peptide_sequences\n')
            script_file_id.write( '\n')
            for i in range(len(database_list)):
                current_code = database_list[i]
                previus_code = database_list[i - 1] if i > 0 else ''
                if i == 0:
                    script_file_id.write(f'# complete peptide sequences -> {current_code}\n')
                else:
                    script_file_id.write(f'# peptide sequences not annotated with {previus_code} -> {current_code}\n')
                script_file_id.write(f'align_peptides_{current_code}_proteome\n')
                script_file_id.write(f'load_alignment_{current_code}_proteome\n')
                script_file_id.write(f'annotate_peptides_{current_code}\n')
                script_file_id.write( '\n')
            if len(database_list) > 1:
                script_file_id.write( '# merged files\n')
                script_file_id.write( 'merge_alignment_files\n')
                script_file_id.write( 'merge_annotation_files\n')
                script_file_id.write( '# -- split_merged_plant_annotation_file\n')
                script_file_id.write( '\n')
            script_file_id.write( '# annotation statistics\n')
            script_file_id.write( 'calculate_annotation_stats\n')
            script_file_id.write( '\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_aminoacid_pipeline_script()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_aminoacid_pipeline_starter(current_run_dir):
    '''
    Build the starter of the script to process a amino acid pipeline.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_aminoacid_pipeline_starter())):
            os.makedirs(os.path.dirname(get_aminoacid_pipeline_starter()))
        with open(get_aminoacid_pipeline_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_aminoacid_pipeline_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')

    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_aminoacid_pipeline_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_script():
    '''
    Get the script path to process a amino acid pipeline.
    '''

    # assign the script path
    aminoacid_pipeline_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_pipeline_aminoacid_code()}-process.sh'

    # return the script path
    return aminoacid_pipeline_script

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_starter():
    '''
    Get the starter path to process a amino acid pipeline.
    '''

    # assign the starter path
    aminoacid_pipeline_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_pipeline_aminoacid_code()}-process-starter.sh'

    # return the starter path
    return aminoacid_pipeline_starter

#-------------------------------------------------------------------------------

def restart_pipeline_process(cluster_name, experiment_id, pipeline_type, pipeline_dataset_id, log, function=None):
    '''
    Restart a pipeline process from the last step ended OK.
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

    # get the starter
    if OK:

        # nucleotide pipelines
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            starter = get_nucleotide_pipeline_starter()

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            starter = get_aminoacid_pipeline_starter()

    # get the current run directory
    if OK:
        current_run_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, pipeline_dataset_id)

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(starter)} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(starter), log)

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

def create_annotation_merger_config_file(experiment_id=xlib.get_toa_result_pipeline_dir(),pipeline_dataset_id_1='toapipelineaa-170101-000000', pipeline_dataset_id_2='toapipelinent-170101-000000', merger_operation='1AND2'):
    '''
    Create FastQC config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the FastQC config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_annotation_merger_config_file())):
            os.makedirs(os.path.dirname(get_annotation_merger_config_file()))
        with open(get_annotation_merger_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the annotation merger parameters\n')
            file_id.write( '[annotation merger parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'pipeline_dataset_id_1 = {pipeline_dataset_id_1}', '# identification of the first pipeline dataset'))
            file_id.write( '{0:<50} {1}\n'.format(f'pipeline_dataset_id_2 = {pipeline_dataset_id_2}', '# identification of the second pipeline dataset'))
            file_id.write( '{0:<50} {1}\n'.format(f'merger_operation = {merger_operation}', f'# merger operation: {xlib.get_annotation_merger_operation_code_list_text()}'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_annotation_merger_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_annotation_merger_process(cluster_name, log, function=None):
    '''
    Run a annotation merger process.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the TOA config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_toa_name()} config file ...\n')
    OK = os.path.isfile(get_toa_config_file())
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write(f'*** ERROR: The {get_toa_config_file()} config file does not exist.\n')
        log.write('Please recreate this file.\n')
        OK = False

    # check the annotation merger config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Checking the {xlib.get_toa_process_merge_annotations_name()} config file ...\n')
        (OK, error_list) = check_annotation_merger_config_file(strict=True)
        if OK:
            log.write('The file is OK.\n')
        else:
            log.write('*** ERROR: The config file is not valid.\n')
            log.write('Please correct this file or recreate it.\n')

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

    # check the TOA is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_toa_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, _) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_toa_name()} is not installed.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_pipeline_dir(), xlib.get_toa_process_merge_annotations_code())

        # create current run directory
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        script = get_annotation_merger_script()
        log.write(f'Building the process script {script} ...\n')
        (OK, error_list) = build_annotation_merger_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')
            log.write('*** ERROR: The file could not be built.\n')

    # upload the script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {script} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(script)}'
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(script)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(script)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        starter = get_annotation_merger_starter()
        log.write(f'Building the process starter {starter} ...\n')
        (OK, error_list) = build_annotation_merger_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {starter} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(starter)}'
        (OK, error_list) = xssh.put_file(sftp_client, starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(starter)} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(starter)}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(starter)} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(starter), log)

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

def check_annotation_merger_config_file(strict):
    '''
    Check the FastQC config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        annotation_merger_option_dict = xlib.get_option_dict(get_annotation_merger_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in annotation_merger_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = annotation_merger_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

        # check section "annotation merger parameters"
        if 'annotation merger parameters' not in sections_list:
            error_list.append('*** ERROR: the section "annotation merger parameters" is not found.')
            OK = False
        else:

            # check section "annotation merger parameters" - key "pipeline_dataset_id_1"
            pipeline_dataset_id_1 = annotation_merger_option_dict.get('annotation merger parameters', {}).get('pipeline_dataset_id_1', not_found)
            if pipeline_dataset_id_1 == not_found:
                error_list.append('*** ERROR: the key "pipeline_dataset_id_1" is not found in the section "annotation merger parameters".')
                OK = False

            # check section "annotation merger parameters" - key "pipeline_dataset_id_2"
            pipeline_dataset_id_2 = annotation_merger_option_dict.get('annotation merger parameters', {}).get('pipeline_dataset_id_2', not_found)
            if pipeline_dataset_id_2 == not_found:
                error_list.append('*** ERROR: the key "pipeline_dataset_id_2" is not found in the section "annotation merger parameters".')
                OK = False

            # check if pipeline_dataset_id_1 and pipeline_dataset_id_2 are different 
            if pipeline_dataset_id_1 == pipeline_dataset_id_2:
                error_list.append('*** ERROR: the "pipeline_dataset_id_1" and "pipeline_dataset_id_2" values hasve to be different.')
                OK = False

            # check section "annotation merger parameters" - key "merger_operation"
            merger_operation = annotation_merger_option_dict.get('annotation merger parameters', {}).get('merger_operation', not_found)
            if merger_operation == not_found:
                error_list.append('*** ERROR: the key "merger_operation" is not found in the section "annotation merger parameters".')
                OK = False
            elif not xlib.check_code(merger_operation, xlib.get_annotation_merger_operation_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "merger_operation" has to be {xlib.get_annotation_merger_operation_code_list_text()}.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_toa_process_merge_annotations_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_annotation_merger_script(cluster_name, current_run_dir):
    '''
    Build the script to process a annotation merger.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # get the pipeline option dictionary
    pipeline_option_dict = xlib.get_option_dict(get_annotation_merger_config_file())

    # get the options
    pipeline_dataset_id_1 = pipeline_option_dict['annotation merger parameters']['pipeline_dataset_id_1']
    pipeline_dataset_id_2 = pipeline_option_dict['annotation merger parameters']['pipeline_dataset_id_2']
    merger_operation = pipeline_option_dict['annotation merger parameters']['merger_operation']

    # write the script
    if OK:
        try:
            if not os.path.exists(os.path.dirname(get_annotation_merger_script())):
                os.makedirs(os.path.dirname(get_annotation_merger_script()))
            with open(get_annotation_merger_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
                script_file_id.write( '#!/bin/bash\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '# output directory\n')
                script_file_id.write(f'OUTPUT_DIR={current_run_dir}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                with open(get_toa_config_file(), mode='r', encoding='iso-8859-1', newline='\n') as toa_config_file_id:
                    records = toa_config_file_id.readlines()
                    for record in records:
                        script_file_id.write(record)
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'MINICONDA_BIN_DIR={toa_config_dict["MINICONDA3_BIN_DIR"]}\n')
                script_file_id.write(f'TOA_DIR={toa_config_dict["TOA_DIR"]}\n')
                script_file_id.write( 'export PATH=$MINICONDA_BIN_DIR:$TOA_DIR:$PATH\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir -p $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'mkdir -p $STATS_DIR\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function init\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    INIT_DATETIME=`date +%s`\n')
                script_file_id.write( '    FORMATTED_INIT_DATETIME=`date "+%Y-%m-%d %H:%M:%S"`\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Script started at $FORMATTED_INIT_DATETIME."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function merge_annotation_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "MERGER OF PLANT ANNOTATION FILES"\n')
                script_file_id.write(f'    ANNOTATION_FILE_1={toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id_1}/`basename $PLANT_ANNOTATION_FILE`\n')
                script_file_id.write( '    ANNOTATION_FILE_1_TMP=$OUTPUT_DIR/pipeline1-`basename $PLANT_ANNOTATION_FILE`.tmp\n')
                script_file_id.write( '    ANNOTATION_FILE_1_SORTED=$OUTPUT_DIR/pipeline1-`basename $PLANT_ANNOTATION_FILE`.sorted\n')
                script_file_id.write(f'    ANNOTATION_FILE_2={toa_config_dict["RESULT_DIR"]}/{xlib.get_toa_result_pipeline_dir()}/{pipeline_dataset_id_2}/`basename $PLANT_ANNOTATION_FILE`\n')
                script_file_id.write( '    ANNOTATION_FILE_2_TMP=$OUTPUT_DIR/pipeline2-`basename $PLANT_ANNOTATION_FILE`.tmp\n')
                script_file_id.write( '    ANNOTATION_FILE_2_SORTED=$OUTPUT_DIR/pipeline2-`basename $PLANT_ANNOTATION_FILE`.sorted\n')
                script_file_id.write( '    echo "Deleting the header record of $ANNOTATION_FILE_1 ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        tail -n +2 $ANNOTATION_FILE_1 > $ANNOTATION_FILE_1_TMP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write( '    echo "Record is deleted."\n')
                script_file_id.write( '    echo "Sorting data records of $ANNOTATION_FILE_1 ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        sort --field-separator=";" --key=2,5 < $ANNOTATION_FILE_1_TMP > $ANNOTATION_FILE_1_SORTED\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '    echo "Records are sorted."\n')
                script_file_id.write( '    echo "Deleting the header record of $ANNOTATION_FILE_2 ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        tail -n +2 $ANNOTATION_FILE_2 > $ANNOTATION_FILE_2_TMP\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tail $RC; fi\n')
                script_file_id.write( '    echo "Record is deleted."\n')
                script_file_id.write( '    echo "Sorting data records of $ANNOTATION_FILE_2 ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        sort --field-separator=";" --key=2,5 < $ANNOTATION_FILE_2_TMP > $ANNOTATION_FILE_2_SORTED\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error sort $RC; fi\n')
                script_file_id.write( '    echo "Records are sorted."\n')
                script_file_id.write( '    echo "Merging annotation files ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        $MINICONDA_BIN_DIR/python3 $TOA_DIR/merge-annotation-files.py \\\n')
                script_file_id.write( '        --file1=$ANNOTATION_FILE_1_SORTED \\\n')
                script_file_id.write( '        --type1=MERGER \\\n')
                script_file_id.write( '        --file2=$ANNOTATION_FILE_2_SORTED \\\n')
                script_file_id.write( '        --type2=MERGER \\\n')
                script_file_id.write(f'        --operation={merger_operation} \\\n')
                script_file_id.write( '        --mfile=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '        --header=Y \\\n')
                script_file_id.write( '        --verbose=N \\\n')
                script_file_id.write( '        --trace=N\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi\n')
                script_file_id.write( '    echo "Files are merged."\n')
                script_file_id.write( '    echo "Deleting temporal files  ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '        rm $ANNOTATION_FILE_1_TMP $ANNOTATION_FILE_2_TMP $ANNOTATION_FILE_1_SORTED $ANNOTATION_FILE_2_SORTED\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
                script_file_id.write( '    echo "Files are deleted."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function calculate_annotation_stats\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $OUTPUT_DIR\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/calculate_annotation_stats.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "CALCULATE ANNOTATION STATISTICS"\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write( '        echo "Calculating stats ..."\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write( '            $MINICONDA_BIN_DIR/python3 $TOA_DIR/calculate-annotation-stats.py \\\n')
                script_file_id.write( '                --db=$TOA_DB \\\n')
                script_file_id.write( '                --transcriptome=NONE \\\n')
                script_file_id.write( '                --peptides=NONE \\\n')
                script_file_id.write( '                --dslist=NONE \\\n')
                script_file_id.write( '                --nonannlist=NONE \\\n')
                script_file_id.write( '                --annotation=$PLANT_ANNOTATION_FILE \\\n')
                script_file_id.write( '                --type=MERGER \\\n')
                script_file_id.write( '                --stats=$ANNOTATION_STATS_FILE \\\n')
                script_file_id.write( '                --verbose=N \\\n')
                script_file_id.write( '                --trace=N\n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error calculate-annotation-stats.py $RC; fi\n')
                script_file_id.write( '        echo "Stats are calculated."\n')
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
                process_name = f'TOA - {xlib.get_toa_process_merge_annotations_name()}'
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
                script_file_id.write( 'merge_annotation_files\n')
                script_file_id.write( 'calculate_annotation_stats\n')
                script_file_id.write( '\n')
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append(f'*** ERROR: The file {get_annotation_merger_script()} can not be created.')
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_annotation_merger_starter(current_run_dir):
    '''
    Build the starter of the script to process a annotation merger.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_annotation_merger_starter())):
            os.makedirs(os.path.dirname(get_annotation_merger_starter()))
        with open(get_annotation_merger_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_annotation_merger_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_annotation_merger_starter()} can not be created.')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_annotation_merger_config_file():
    '''
    Get the annotation merger config file path.
    '''

    # assign the annotation merger config file path
    annotation_merger_config_file = f'{xlib.get_config_dir()}/{xlib.get_toa_process_merge_annotations_code()}-config.txt'

    # return the annotation merger config file path
    return annotation_merger_config_file

#-------------------------------------------------------------------------------

def get_annotation_merger_script():
    '''
    Get the script path to process a annotation merger.
    '''

    # assign the script path
    annotation_merger_script = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_merge_annotations_code()}-process.sh'

    # return the script path
    return annotation_merger_script

#-------------------------------------------------------------------------------

def get_annotation_merger_starter():
    '''
    Get the starter path to process a annotation merger.
    '''

    # assign the starter path
    annotation_merger_starter = f'{xlib.get_temp_dir()}/{xlib.get_toa_process_merge_annotations_code()}-process-starter.sh'

    # return the starter path
    return annotation_merger_starter

#-------------------------------------------------------------------------------

def get_nucleotide_annotation_database_code_list():
    '''
    Get the code list of "nucleotide_annotation_database".
    '''

    return ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant', 'nt']

#-------------------------------------------------------------------------------

def get_nucleotide_annotation_database_code_list_text():
    '''
    Get the code list of "nucleotide_annotation_database" as text.
    '''

    return str(get_nucleotide_annotation_database_code_list()).strip('[]').replace('\'', '').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_aminoacid_annotation_database_code_list():
    '''
    Get the code list of "aminoacid_annotation_database".
    '''

    return ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant', 'nr']

#-------------------------------------------------------------------------------

def get_aminoacid_annotation_database_code_list_text():
    '''
    Get the code list of "aminoacid_annotation_database" as text.
    '''

    return str(get_aminoacid_annotation_database_code_list()).strip('[]').replace('\'', '').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_assembly_origin_code_list():
    '''
    Get the code list of "assembly_origin".
    '''

    return ['NGSCLOUD', 'EXTERNAL']

#-------------------------------------------------------------------------------

def get_assembly_origin_code_list_text():
    '''
    Get the code list of "assembly_origin" as text.
    '''

    return str(get_assembly_origin_code_list()).strip('[]').replace('\'', '').replace(',', ' or')

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
     print('This file contains functions related to the TOA (Tree-oriented Annotation) process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
