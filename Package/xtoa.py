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
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(cluster_app_dir, xlib.get_toa_name())
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

def install_toa(cluster_name, log, function=None):
    '''
    Install the TOA software in the cluster.
    '''

    # initialize the control variable
    OK = True

    # get the version and download URL of TOA
    (toa_version, toa_url, toa_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_toa_name())

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
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(cluster_app_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: There is not any volume mounted in the directory.\n')
            log.write('You have to link a volume in the mounting point {0} for the cluster {1}.\n'.format(cluster_app_dir, cluster_name))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_toa_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TOA installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_toa_installation_script()))
        (OK, error_list) = build_toa_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the TOA installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} ...\n'.format(get_toa_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_toa_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_toa_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TOA installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_toa_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_toa_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TOA installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_toa_installation_starter()))
        (OK, error_list) = build_toa_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the TOA installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} ...\n'.format(get_toa_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_toa_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_toa_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TOA installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_toa_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_toa_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the TOA installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_toa_installation_starter())))
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
    (toa_version, toa_url, toa_channel) = xconfiguration.get_bioinfo_app_data(xlib.get_toa_name())

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
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_toa_installation_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_toa_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_toa_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_toa_installation_script():
    '''
    Get the TOA installation script path in the local computer.
    '''

    # assign the script path
    toa_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_toa_name())

    # return the script path
    return toa_installation_script

#-------------------------------------------------------------------------------

def get_toa_installation_starter():
    '''
    Get the TOA installation starter path in the local computer.
    '''

    # assign the starter path
    toa_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_name())

    # return the starter path
    return toa_installation_starter

#-------------------------------------------------------------------------------

def create_toa_config_file(toa_dir=f'{xlib.get_cluster_app_dir()}/TOA/Package', miniconda3_bin_dir=f'{xlib.get_cluster_app_dir()}/Miniconda3/bin', db_dir=xlib.get_cluster_database_dir(), result_dir=xlib.get_cluster_result_dir()):
    '''
    Create the TOA config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the NGScloud config file
    toa_config_file = get_toa_config_file()

    # create the TOA config file
    if OK:
        try:
            if not os.path.exists(os.path.dirname(toa_config_file)):
                os.makedirs(os.path.dirname(toa_config_file))
            with open(toa_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
                file_id.write( '{0}\n'.format('# environment'))
                file_id.write( '{0}\n'.format('TOA_DIR={0}'.format(toa_dir)))
                file_id.write( '{0}\n'.format('MINICONDA3_BIN_DIR={0}'.format(miniconda3_bin_dir)))
                file_id.write( '{0}\n'.format('RESULT_DIR={0}'.format(result_dir)))
                file_id.write( '{0}\n'.format('DB_DIR={0}'.format(db_dir)))
                file_id.write( '{0}\n'.format('TOA_DB_DIR={0}/TOA'.format(db_dir)))
                file_id.write( '{0}\n'.format('DATA_DIR={0}/data'.format(db_dir)))
                file_id.write( '{0}\n'.format('PLAZA_DIR={0}/PLAZA'.format(db_dir)))
                file_id.write( '{0}\n'.format('NCBI_DIR={0}/NCBI'.format(db_dir)))
                file_id.write( '{0}\n'.format('INTERPRO_DIR={0}/InterPro'.format(db_dir)))
                file_id.write( '{0}\n'.format('GO_DIR={0}/GO'.format(db_dir)))
                file_id.write( '{0}\n'.format('EC_DIR={0}/EC'.format(db_dir)))
                file_id.write( '{0}\n'.format('KEGG_DIR={0}/KEGG'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# TOA database'))
                file_id.write( '{0}\n'.format('TOA_DB={0}/TOA/toa.db'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# basic data'))
                file_id.write( '{0}\n'.format('DATASET_FILE={0}/data/datasets.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('SPECIES_FILE={0}/data/species.txt'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# other basic data'))
                file_id.write( '{0}\n'.format('EC_IDS_FTP=ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'))
                file_id.write( '{0}\n'.format('EC_IDS_FILE={0}/EC/enzyme.dat'.format(db_dir)))
                file_id.write( '{0}\n'.format('KEGG_IDS_FTP=ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'))
                file_id.write( '{0}\n'.format('KEGG_IDS_FILE={0}/KEGG/ko_list.gz'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# Gymno PLAZA 1.0'))
                file_id.write( '{0}\n'.format('GYMNO_01_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Fasta/cds.csv.gz'))
                file_id.write( '{0}\n'.format('GYMNO_01_CDS_FILE={0}/PLAZA/gymno_01-cds.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Fasta/proteome.csv.gz'))
                file_id.write( '{0}\n'.format('GYMNO_01_PROTEOME_FILE={0}/PLAZA/gymno_01-proteome.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/Descriptions'))
                file_id.write( '{0}\n'.format('GYMNO_01_GENEDESC_DIR={0}/PLAZA/gymno_01-descriptions'.format(db_dir)))
                file_id.write( '{0}\n'.format("GYMNO_01_GENEDESC_FILE_PATTERN='gene_description.*.csv.gz'"))
                file_id.write( '{0}\n'.format('GYMNO_01_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/InterPro/interpro.csv.gz'))
                file_id.write( '{0}\n'.format('GYMNO_01_INTERPRO_FILE={0}/PLAZA/gymno_01-interpro.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/GO/go.csv.gz'))
                file_id.write( '{0}\n'.format('GYMNO_01_GO_FILE={0}/PLAZA/gymno_01-go.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/MapMan/mapman.csv.gz'))
                file_id.write( '{0}\n'.format('GYMNO_01_MAPMAN_FILE={0}/PLAZA/gymno_01-mapman.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_PROTEOME_DB_NAME=gymno_01_proteome'))
                file_id.write( '{0}\n'.format('GYMNO_01_PROTEOME_DB_DIR={0}/PLAZA/gymno_01_proteome_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_PROTEOME_DB_FILE={0}/PLAZA/gymno_01_proteome_db/gymno_01_proteome'.format(db_dir)))
                file_id.write( '{0}\n'.format('GYMNO_01_BLAST_XML=$OUTPUT_DIR/gymno_01-alignment.xml'))
                file_id.write( '{0}\n'.format('GYMNO_01_ANNOTATION_FILE=$OUTPUT_DIR/gymno_01-annotation.csv'))
                file_id.write( '{0}\n'.format('GYMNO_01_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/gymno_01-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('GYMNO_01_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/gymno_01-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# Dicots PLAZA 4.0'))
                file_id.write( '{0}\n'.format('DICOTS_04_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/cds.all_transcripts.fasta.gz'))
                file_id.write( '{0}\n'.format('DICOTS_04_CDS_FILE={0}/PLAZA/dicots_04-cds.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/proteome.all_transcripts.fasta.gz'))
                file_id.write( '{0}\n'.format('DICOTS_04_PROTEOME_FILE={0}/PLAZA/dicots_04-proteome.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Descriptions'))
                file_id.write( '{0}\n'.format('DICOTS_04_GENEDESC_DIR={0}/PLAZA/dicots_04-descriptions'.format(db_dir)))
                file_id.write( '{0}\n'.format("DICOTS_04_GENEDESC_FILE_PATTERN='gene_description.*.csv.gz'"))
                file_id.write( '{0}\n'.format('DICOTS_04_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/InterPro/interpro.csv.gz'))
                file_id.write( '{0}\n'.format('DICOTS_04_INTERPRO_FILE={0}/PLAZA/dicots_04-interpro.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/GO/go.csv.gz'))
                file_id.write( '{0}\n'.format('DICOTS_04_GO_FILE={0}/PLAZA/dicots_04-go.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/MapMan/mapman.csv.gz'))
                file_id.write( '{0}\n'.format('DICOTS_04_MAPMAN_FILE={0}/PLAZA/dicots_04-mapman.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_PROTEOME_DB_NAME=dicots_04_proteome'))
                file_id.write( '{0}\n'.format('DICOTS_04_PROTEOME_DB_DIR={0}/PLAZA/dicots_04_proteome_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_PROTEOME_DB_FILE={0}/PLAZA/dicots_04_proteome_db/dicots_04_proteome'.format(db_dir)))
                file_id.write( '{0}\n'.format('DICOTS_04_BLAST_XML=$OUTPUT_DIR/dicots_04-alignment.xml'))
                file_id.write( '{0}\n'.format('DICOTS_04_ANNOTATION_FILE=$OUTPUT_DIR/dicots_04-annotation.csv'))
                file_id.write( '{0}\n'.format('DICOTS_04_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/dicots_04-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('DICOTS_04_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/dicots_04-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# Monocots PLAZA 4.0'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_CDS_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Fasta/cds.all_transcripts.fasta.gz'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_CDS_FILE={0}/PLAZA/monocots_04-cds.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_PROTEOME_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Fasta/proteome.all_transcripts.fasta.gz'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_PROTEOME_FILE={0}/PLAZA/monocots_04-proteome.fasta.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_GENEDESC_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Descriptions'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_GENEDESC_DIR={0}/PLAZA/monocots_04-descriptions'.format(db_dir)))
                file_id.write( '{0}\n'.format("MONOCOTS_04_GENEDESC_FILE_PATTERN='gene_description.*.csv.gz'"))
                file_id.write( '{0}\n'.format('MONOCOTS_04_INTERPRO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/InterPro/interpro.csv.gz'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_INTERPRO_FILE={0}/PLAZA/monocots_04-interpro.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_GO_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/GO/go.csv.gz'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_GO_FILE={0}/PLAZA/monocots_04-go.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_MAPMAN_FTP=ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/MapMan/mapman.csv.gz'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_MAPMAN_FILE={0}/PLAZA/monocots_04-mapman.csv.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_PROTEOME_DB_NAME=monocots_04_proteome'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_PROTEOME_DB_DIR={0}/PLAZA/monocots_04_proteome_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_PROTEOME_DB_FILE={0}/PLAZA/monocots_04_proteome_db/monocots_04_proteome'.format(db_dir)))
                file_id.write( '{0}\n'.format('MONOCOTS_04_BLAST_XML=$OUTPUT_DIR/monocots_04-alignment.xml'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_ANNOTATION_FILE=$OUTPUT_DIR/monocots_04-annotation.csv'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/monocots_04-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('MONOCOTS_04_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/monocots_04-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI RefSeq Plant'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_FTP=ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_LOCAL={0}/NCBI/ftp.ncbi.nlm.nih.gov/refseq/release/plant/'.format(db_dir)))
                file_id.write( '{0}\n'.format("REFSEQ_PROTEIN_FILE_PATTERN='.protein.faa.gz'"))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_PROTEOME_FILE={0}/NCBI/refseq_plant-proteome.fasta'.format(db_dir)))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_PROTEOME_DB_NAME=refseq_plant_proteome'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_PROTEOME_DB_DIR={0}/NCBI/refseq_plant_proteome_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_PROTEOME_DB_FILE={0}/NCBI/refseq_plant_proteome_db/refseq_plant_proteome'.format(db_dir)))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_FILE_LIST={0}/NCBI/refseq_plant_proteome_db/files_list.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_BLAST_XML=$OUTPUT_DIR/refseq_plant-alignment.xml'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_ANNOTATION_FILE=$OUTPUT_DIR/refseq_plant-annotation.csv'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/refseq_plant-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('REFSEQ_PLANT_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/refseq_plant-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI BLAST database'))
                # -- file_id.write( '{0}\n'.format('BLAST_DATABASES_FTP=ftp://ftp.ncbi.nlm.nih.gov/blast/db/'))
                file_id.write( '{0}\n'.format('BLAST_DATABASES_FTP=ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/'))
                # -- file_id.write( '{0}\n'.format('BLAST_DATABASES_LOCAL={0}/NCBI/ftp.ncbi.nlm.nih.gov/blast/db/'.format(db_dir)))
                file_id.write( '{0}\n'.format('BLAST_DATABASES_LOCAL={0}/NCBI/ftp.ncbi.nlm.nih.gov/blast/db/v4/'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI BLAST database NT'))
                file_id.write( '{0}\n'.format('NT_DB_NAME=nt'))
                # -- file_id.write( '{0}\n'.format("NT_FILE_PATTERN='nt.*.tar.gz'"))
                file_id.write( '{0}\n'.format("NT_FILE_PATTERN='nt_v4.*.tar.gz'"))
                file_id.write( '{0}\n'.format('NT_DB_DIR={0}/NCBI/nt_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('NT_FILE_LIST={0}/NCBI/nt_db/files_list.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('NT_VIRIDIPLANTAE_BLAST_XML=$OUTPUT_DIR/nt-viridiplantae-alignment.xml'))
                file_id.write( '{0}\n'.format('NT_VIRIDIPLANTAE_ANNOTATION_FILE=$OUTPUT_DIR/nt-viridiplantae-annotation.csv'))
                file_id.write( '{0}\n'.format('NT_VIRIDIPLANTAE_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/nt-viridiplantae-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('REIDENTIFIED_NT_REMAINDER_BLAST_XML=$OUTPUT_DIR/reidentified-nt-remainder-alignment.xml'))
                file_id.write( '{0}\n'.format('NT_REMAINDER_BLAST_XML=$OUTPUT_DIR/nt-remainder-alignment.xml'))
                file_id.write( '{0}\n'.format('NT_REMAINDER_ANNOTATION_FILE=$OUTPUT_DIR/nt-remainder-annotation.csv'))
                file_id.write( '{0}\n'.format('NT_REMAINDER_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/nt-remainder-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('NR_REMAINDER_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/nr-remainder-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI Nucleotide GenInfo identifier lists'))
                file_id.write( '{0}\n'.format('NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST={0}/NCBI/nucleotide_viridiplantae.gi'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI BLAST database NR'))
                file_id.write( '{0}\n'.format('NR_DB_NAME=nr'))
                # -- file_id.write( '{0}\n'.format("NR_FILE_PATTERN='nr.*.tar.gz'"))
                file_id.write( '{0}\n'.format("NR_FILE_PATTERN='nr_v4.*.tar.gz'"))
                file_id.write( '{0}\n'.format('NR_DB_DIR={0}/NCBI/nr_db'.format(db_dir)))
                file_id.write( '{0}\n'.format('NR_FILE_LIST={0}/NCBI/nr_db/files_list.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('NR_VIRIDIPLANTAE_BLAST_XML=$OUTPUT_DIR/nr-viridiplantae-alignment.xml'))
                file_id.write( '{0}\n'.format('NR_VIRIDIPLANTAE_ANNOTATION_FILE=$OUTPUT_DIR/nr-viridiplantae-annotation.csv'))
                file_id.write( '{0}\n'.format('NR_VIRIDIPLANTAE_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/nr-viridiplantae-nonann-peptides.fasta'))
                file_id.write( '{0}\n'.format('REIDENTIFIED_NR_REMAINDER_BLAST_XML=$OUTPUT_DIR/reidentified-nr-remainder-alignment.xml'))
                file_id.write( '{0}\n'.format('NR_REMAINDER_BLAST_XML=$OUTPUT_DIR/nr-remainder-alignment.xml'))
                file_id.write( '{0}\n'.format('NR_REMAINDER_ANNOTATION_FILE=$OUTPUT_DIR/nr-remainder-annotation.csv'))
                file_id.write( '{0}\n'.format('NT_REMAINDER_NON_ANNOTATED_TRANSCRIPT_FILE=$OUTPUT_DIR/nt-remainder-nonann-transcripts.fasta'))
                file_id.write( '{0}\n'.format('NR_REMAINDER_NON_ANNOTATED_PEPTIDE_FILE=$OUTPUT_DIR/nr-remainder-nonann-peptides.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI Protein GenInfo identifier lists'))
                file_id.write( '{0}\n'.format('PROTEIN_VIRIDIPLANTAE_GI_LIST={0}/NCBI/protein_viridiplantae.gi'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# NCBI Gene'))
                file_id.write( '{0}\n'.format('GENE_GENE2GO_FTP=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'))
                file_id.write( '{0}\n'.format('GENE_GENE2GO_FILE={0}/NCBI/gene-gene2go.gz'.format(db_dir)))
                file_id.write( '{0}\n'.format('GENE_GENE2REFSEQ_FTP=ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'))
                file_id.write( '{0}\n'.format('GENE_GENE2REFSEQ_FILE={0}/NCBI/gene-gene2refseq.gz'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# InterPro'))
                file_id.write( '{0}\n'.format('INTERPRO_INTERPRO2GO_FTP=ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go'))
                file_id.write( '{0}\n'.format('INTERPRO_INTERPRO2GO_FILE={0}/InterPro/interpro2go'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# Gene Onlogolgy'))
                file_id.write( '{0}\n'.format('GO_ONTOLOGY_FTP=http://purl.obolibrary.org/obo/go.obo'))
                file_id.write( '{0}\n'.format('GO_ONTOLOGY_FILE={0}/GO/go.obo'.format(db_dir)))
                file_id.write( '{0}\n'.format('#GO_EC2GO_FTP=http://geneontology.org/external2go/ec2go'))
                file_id.write( '{0}\n'.format('GO_EC2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/ec2go'))
                file_id.write( '{0}\n'.format('GO_EC2GO_FILE={0}/GO/ec2go.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('#GO_KEGG2GO_FTP=http://geneontology.org/external2go/kegg2go'))
                file_id.write( '{0}\n'.format('GO_KEGG2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/kegg2go'))
                file_id.write( '{0}\n'.format('GO_KEGG2GO_FILE={0}/GO/kegg2go.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('#GO_METACYC2GO_FTP=http://geneontology.org/external2go/metacyc2go'))
                file_id.write( '{0}\n'.format('GO_METACYC2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/metacyc2go'))
                file_id.write( '{0}\n'.format('GO_METACYC2GO_FILE={0}/GO/metacyc2go.txt'.format(db_dir)))
                file_id.write( '{0}\n'.format('#GO_INTERPRO2GO_FTP=http://geneontology.org/external2go/interpro2go'))
                file_id.write( '{0}\n'.format('GO_INTERPRO2GO_FTP=https://build.berkeleybop.org/view/GO/job/Update%20external2go/lastSuccessfulBuild/artifact/external2go/interpro2go'))
                file_id.write( '{0}\n'.format('GO_INTERPRO2GO_FILE={0}/GO/interpro2go.txt'.format(db_dir)))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# other transcriptome files'))
                file_id.write( '{0}\n'.format('REIDENTIFIED_TRANSCRIPTOME_FILE=$OUTPUT_DIR/reidentified-transcriptome.fasta'))
                file_id.write( '{0}\n'.format('RELATIONSHIP_FILE=$OUTPUT_DIR/relationships.csv'))
                file_id.write( '{0}\n'.format('PURGED_TRANSCRIPTOME_FILE=$OUTPUT_DIR/purged-transcriptome.fasta'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# peptide sequence files'))
                file_id.write( '{0}\n'.format('TRANSDECODER_OUTPUT_DIR=$OUTPUT_DIR/transdecoder'))
                file_id.write( '{0}\n'.format('ORF_FILE=$TRANSDECODER_OUTPUT_DIR/longest_orfs.pep'))
                file_id.write( '{0}\n'.format('PEPTIDE_FILE=$OUTPUT_DIR/`basename "$TRANSCRIPTOME_FILE"`.transdecoder.pep'))
                file_id.write( '{0}\n'.format('REIDENTIFIED_PEPTIDE_FILE=$OUTPUT_DIR/reidentified-peptides.fasta'))
                file_id.write( '{0}\n'.format('RELATIONSHIP_FILE=$OUTPUT_DIR/relationships.csv'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# merger files'))
                file_id.write( '{0}\n'.format('MERGED_ANNOTATION_FILE=$OUTPUT_DIR/merged-annotation.csv'))
                file_id.write( '{0}\n'.format('PLANT_ANNOTATION_FILE=$OUTPUT_DIR/plant-annotation.csv'))
                file_id.write( '{0}\n'.format('REIDENTIFIED_PLANT_BLAST_XML=$OUTPUT_DIR/reidentified-plant-alignment.xml'))
                file_id.write( '{0}\n'.format('PLANT_BLAST_XML=$OUTPUT_DIR/plant-alignment.xml'))
                file_id.write( '\n')
                file_id.write( '{0}\n'.format('# statistics'))
                file_id.write( '{0}\n'.format('STATS_SUBDIR_NAME=stats'))
                file_id.write( '{0}\n'.format('STATS_DIR=$OUTPUT_DIR/stats'))
                file_id.write( '{0}\n'.format('STATS_BASE_NAME=stats'))
                file_id.write( '{0}\n'.format('ANNOTATION_STATS_FILE=$OUTPUT_DIR/stats/stats.csv'))
        except Exception as e:
            error_list.append(f'*** EXCEPTION: "{e}".')
            error_list.append('*** ERROR: The file {0} can not be created'.format(toa_config_file))
            OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_toa_config_file():
    '''
    Get the path of the TOA config file corresponding to the environment.
    '''

    # assign the config file
    toa_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_toa_code())

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
    except Exception as e:
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
            file_id.write( '{0}\n'.format('# This file contains the data of the genomic datasets used by Tree-oriented Annotation (TOA) software package.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# RECORD FORMAT: "dataset_id";"dataset_name";"repository_id";"ftp_adress"'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('"dicots_04";"Dicots PLAZA 4.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/"'))
            file_id.write( '{0}\n'.format('"gene";"Gene";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"'))
            file_id.write( '{0}\n'.format('"gymno_01";"Gymno PLAZA 1.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_gymno_01/"'))
            file_id.write( '{0}\n'.format('"monocots_04";"Monocots PLAZA 4.0";"plaza";"ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/"'))
            file_id.write( '{0}\n'.format('"nr_remainder";"nr remainder";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"'))
            file_id.write( '{0}\n'.format('"nr_viridiplantae";"nr Viridiplantae";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"'))
            file_id.write( '{0}\n'.format('"nt_remainder";"nt remainder";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"'))
            file_id.write( '{0}\n'.format('"nt_viridiplantae";"nt Viridiplantae";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/blast/db/"'))
            file_id.write( '{0}\n'.format('"refseq_plant";"RefSeq Plant";"ncbi";"ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/"'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_dataset_file()))
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
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_dataset_file()))
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
                    dataset_id = mo.group(1).strip()
                    dataset_name = mo.group(2).strip()
                    repository_id = mo.group(3).strip()
                    ftp_adress = mo.group(4).strip()
                except Exception as e:
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

            # read the next record
            record = dataset_file_id.readline()

    # close the file of datasets
    dataset_file_id.close()

    # warn that the file of datasets is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_dataset_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dataset_file():
    '''
    Get the dataset file path.
    '''

    # assign the dataset file path
    dataset_file = '{0}/datasets.txt'.format(xlib.get_config_dir())

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

    # create the file of restriction sites and write the default data
    try:
        if not os.path.exists(os.path.dirname(get_species_file())):
            os.makedirs(os.path.dirname(get_species_file()))
        with open(get_species_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# This file contains the data of specioes used by Tree-oriented Annotation (TOA) software package.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# RECORD FORMAT: "species_name";"plaza_id"'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('"Actinidia chinensis";"ach"'))
            file_id.write( '{0}\n'.format('"Amborella trichopoda";"atr"'))
            file_id.write( '{0}\n'.format('"Ananas comosus";"aco"'))
            file_id.write( '{0}\n'.format('"Arabidopsis lyrata";"aly"'))
            file_id.write( '{0}\n'.format('"Arabidopsis thaliana";"ath"'))
            file_id.write( '{0}\n'.format('"Arachis ipaensis";"aip"'))
            file_id.write( '{0}\n'.format('"Beta vulgaris";"bvu"'))
            file_id.write( '{0}\n'.format('"Brachypodium distachyon";"bdi"'))
            file_id.write( '{0}\n'.format('"Brassica oleracea";"bol"'))
            file_id.write( '{0}\n'.format('"Brassica rapa";"bra"'))
            file_id.write( '{0}\n'.format('"Cajanus cajan";"ccaj"'))
            file_id.write( '{0}\n'.format('"Capsella rubella";"cru"'))
            file_id.write( '{0}\n'.format('"Capsicum annuum";"can"'))
            file_id.write( '{0}\n'.format('"Carica papaya";"cpa"'))
            file_id.write( '{0}\n'.format('"Chenopodium quinoa";"cqu"'))
            file_id.write( '{0}\n'.format('"Chlamydomonas reinhardtii";"cre"'))
            file_id.write( '{0}\n'.format('"Cicer arietinum";"car"'))
            file_id.write( '{0}\n'.format('"Citrullus lanatus";"cla"'))
            file_id.write( '{0}\n'.format('"Citrus clementina";"ccl"'))
            file_id.write( '{0}\n'.format('"Coffea canephora";"ccan"'))
            file_id.write( '{0}\n'.format('"Corchorus olitorius";"col"'))
            file_id.write( '{0}\n'.format('"Cucumis melo";"cme"'))
            file_id.write( '{0}\n'.format('"Cucumis sativus L.";"csa"'))
            file_id.write( '{0}\n'.format('"Cycas micholitzii";"cmi"'))
            file_id.write( '{0}\n'.format('"Daucus carota";"dca"'))
            file_id.write( '{0}\n'.format('"Elaeis guineensis";"egu"'))
            file_id.write( '{0}\n'.format('"Erythranthe guttata";"egut"'))
            file_id.write( '{0}\n'.format('"Eucalyptus grandis";"egr"'))
            file_id.write( '{0}\n'.format('"Fragaria vesca";"fve"'))
            file_id.write( '{0}\n'.format('"Ginkgo biloba";"gbi"'))
            file_id.write( '{0}\n'.format('"Glycine max";"gma"'))
            file_id.write( '{0}\n'.format('"Gnetum montanum";"gmo"'))
            file_id.write( '{0}\n'.format('"Gossypium raimondii";"gra"'))
            file_id.write( '{0}\n'.format('"Hevea brasiliensis";"hbr"'))
            file_id.write( '{0}\n'.format('"Hordeum vulgare";"hvu"'))
            file_id.write( '{0}\n'.format('"Malus domestica";"mdo"'))
            file_id.write( '{0}\n'.format('"Manihot esculenta";"mes"'))
            file_id.write( '{0}\n'.format('"Marchantia polymorpha";"mpo"'))
            file_id.write( '{0}\n'.format('"Medicago truncatula";"mtr"'))
            file_id.write( '{0}\n'.format('"Micromonas commoda";"mco"'))
            file_id.write( '{0}\n'.format('"Musa acuminata";"mac"'))
            file_id.write( '{0}\n'.format('"Nelumbo nucifera";"nnu"'))
            file_id.write( '{0}\n'.format('"Oropetium thomaeum";"oth"'))
            file_id.write( '{0}\n'.format('"Oryza brachyantha";"obr"'))
            file_id.write( '{0}\n'.format('"Oryza sativa ssp. indica";"osaindica"'))
            file_id.write( '{0}\n'.format('"Oryza sativa ssp. japonica";"osa"'))
            file_id.write( '{0}\n'.format('"Petunia axillaris";"pax"'))
            file_id.write( '{0}\n'.format('"Phalaenopsis equestris";"peq"'))
            file_id.write( '{0}\n'.format('"Phyllostachys edulis";"ped"'))
            file_id.write( '{0}\n'.format('"Physcomitrella patens";"ppa"'))
            file_id.write( '{0}\n'.format('"Picea abies";"pab"'))
            file_id.write( '{0}\n'.format('"Picea glauca";"pgl"'))
            file_id.write( '{0}\n'.format('"Picea sitchensis";"psi"'))
            file_id.write( '{0}\n'.format('"Pinus pinaster";"ppi"'))
            file_id.write( '{0}\n'.format('"Pinus sylvestris";"psy"'))
            file_id.write( '{0}\n'.format('"Pinus taeda";"pta"'))
            file_id.write( '{0}\n'.format('"Populus trichocarpa";"ptr"'))
            file_id.write( '{0}\n'.format('"Prunus persica";"ppe"'))
            file_id.write( '{0}\n'.format('"Pseudotsuga menziesii";"pme"'))
            file_id.write( '{0}\n'.format('"Pyrus bretschneideri";"pbr"'))
            file_id.write( '{0}\n'.format('"Ricinus communis";"rco"'))
            file_id.write( '{0}\n'.format('"Schrenkiella parvula";"spa"'))
            file_id.write( '{0}\n'.format('"Selaginella moellendorffii";"smo"'))
            file_id.write( '{0}\n'.format('"Setaria italica";"sit"'))
            file_id.write( '{0}\n'.format('"Solanum lycopersicum";"sly"'))
            file_id.write( '{0}\n'.format('"Solanum tuberosum";"stu"'))
            file_id.write( '{0}\n'.format('"Sorghum bicolor";"sbi"'))
            file_id.write( '{0}\n'.format('"Spirodela polyrhiza";"spo"'))
            file_id.write( '{0}\n'.format('"Tarenaya hassleriana";"tha"'))
            file_id.write( '{0}\n'.format('"Taxus baccata";"tba"'))
            file_id.write( '{0}\n'.format('"Theobroma cacao";"tca"'))
            file_id.write( '{0}\n'.format('"Trifolium pratense";"tpr"'))
            file_id.write( '{0}\n'.format('"Triticum aestivum";"tae"'))
            file_id.write( '{0}\n'.format('"Utricularia gibba";"ugi"'))
            file_id.write( '{0}\n'.format('"Vigna radiata var. radiata";"vra"'))
            file_id.write( '{0}\n'.format('"Vitis vinifera";"vvi"'))
            file_id.write( '{0}\n'.format('"Zea mays";"zma"'))
            file_id.write( '{0}\n'.format('"Ziziphus jujuba";"zju"'))
            file_id.write( '{0}\n'.format('"Zostera marina";"zosmarina"'))
            file_id.write( '{0}\n'.format('"Zoysia japonica ssp. nagirizaki";"zjn"'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_species_file()))
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
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_species_file()))
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
                    species_name = mo.group(1).strip()
                    plaza_species_id = mo.group(2).strip()
                except Exception as e:
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

            # read the next record
            record = species_file_id.readline()

    # close the file of species
    species_file_id.close()

    # warn that the file of species is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_species_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_species_file():
    '''
    Get the species file path.
    '''

    # assign the species file path
    species_file = '{0}/species.txt'.format(xlib.get_config_dir())

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
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_toa_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_toa_name()))
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
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if process_type == xlib.get_toa_type_recreate():
            script = get_recreate_toa_database_script()
            log.write('Building the process script {0} ...\n'.format(script))
            (OK, error_list) = build_recreate_toa_database_script(cluster_name, current_run_dir)
        elif process_type == xlib.get_toa_type_rebuild(): 
            script = get_rebuild_toa_database_script()
            log.write('Building the process script {0} ...\n'.format(script))
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
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(script, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(script)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the script starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if process_type == xlib.get_toa_type_recreate():
            starter = get_recreate_toa_database_starter()
            log.write('Building the process starter {0} ...\n'.format(starter))
            (OK, error_list) = build_recreate_toa_database_starter(current_run_dir)
        elif process_type == xlib.get_toa_type_rebuild():
            starter = get_rebuild_toa_database_starter()
            log.write('Building the process starter {0} ...\n'.format(starter))
            (OK, error_list) = build_rebuild_toa_database_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(get_recreate_toa_database_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_recreate_toa_database_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_recreate_toa_database_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_recreate_toa_database_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_recreate_toa_database_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_recreate_toa_database_starter())))
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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function create_toa_database_dir'))
                script_file_id.write( '{\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Creating the database directory ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        mkdir --parents {0}'.format(toa_config_dict['TOA_DB_DIR'])))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-ncbi-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function recreate_toa_database'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        recreate-database.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-ncbi-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('create_toa_database_dir'))
                script_file_id.write( '{0}\n'.format('recreate_toa_database'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_recreate_toa_database_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_recreate_toa_database_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_recreate_toa_database_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_recreate_toa_database_script():
    '''
    Get the script path in the local computer to recreate a TOA database.
    '''

    # assign the script path
    recreate_toa_database_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_recreate_toa_database_code())

    # return the script path
    return recreate_toa_database_script

#-------------------------------------------------------------------------------

def get_recreate_toa_database_starter():
    '''
    Get the starter path in the local computer to recreate a TOA database.
    '''

    # assign the starter path
    recreate_toa_database_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_recreate_toa_database_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function rebuild_toa_database'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        rebuild-database.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-ncbi-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('rebuild_toa_database'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_rebuild_toa_database_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_rebuild_toa_database_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_rebuild_toa_database_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_rebuild_toa_database_script():
    '''
    Get the script path in the local computer to rebuild a TOA database.
    '''

    # assign the script path
    rebuild_toa_database_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_rebuild_toa_database_code())

    # return the script path
    return rebuild_toa_database_script

#-------------------------------------------------------------------------------

def get_rebuild_toa_database_starter():
    '''
    Get the starter path in the local computer to rebuild a TOA database.
    '''

    # assign the starter path
    rebuild_toa_database_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_rebuild_toa_database_code())

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
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_toa_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_toa_name()))
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

        # processes to build BLAST databases
        elif process_type == xlib.get_toa_type_build_blastdb():
            if genomic_database == xlib.get_toa_data_nt_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_blastdb_nt_code())
            elif genomic_database == xlib.get_toa_data_nr_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_blastdb_nr_code())

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code())
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                current_run_dir = xlib.get_cluster_current_run_dir(xlib.get_toa_result_database_dir(), xlib.get_toa_process_gilist_viridiplantae_protein_gi_code())

        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # create the data subdirectory from the database directory
    if OK:
        if process_type == xlib.get_toa_type_load_data() and genomic_database == xlib.get_toa_data_basic_data_code():
            log.write(f'{xlib.get_separator()}\n')
            log.write('Creating the TOA data directory ...\n')
            command = 'mkdir --parents {0}'.format(data_dir)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The directory path {0} is created.\n'.format(data_dir))
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
                log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_dataset_file(), data_dir))
                cluster_path = '{0}/{1}'.format(data_dir, os.path.basename(get_dataset_file()))
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
            log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_species_file(), data_dir))
            cluster_path = '{0}/{1}'.format(data_dir, os.path.basename(get_species_file()))
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
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_gymno_01_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_proteome_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_dicots_04_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_proteome_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_monocots_04_proteome_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_refseq_plant_code():
                script = get_refseq_plant_proteome_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_refseq_plant_proteome_script(cluster_name, current_run_dir)

        # processes to download functional annotations from a genomic database server
        elif process_type == xlib.get_toa_type_download_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                script = get_basic_data_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_basic_data_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                script = get_gymno_01_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_gymno_01_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_dicots_04_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_monocots_04_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                script = get_gene_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_gene_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                script = get_interpro_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_interpro_download_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                script = get_go_download_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_go_download_script(cluster_name, current_run_dir)

        # processes to load data of a genomic database into TOA database
        elif process_type == xlib.get_toa_type_load_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                script = get_basic_data_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_basic_data_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                script = get_gymno_01_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_gymno_01_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                script = get_dicots_04_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_dicots_04_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                script = get_monocots_04_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_monocots_04_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                script = get_gene_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_gene_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                script = get_interpro_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_interpro_load_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                script = get_go_load_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_go_load_script(cluster_name, current_run_dir)

        # processes to build BLAST databases
        elif process_type == xlib.get_toa_type_build_blastdb():
            if genomic_database == xlib.get_toa_data_nt_code():
                script = get_nt_blastdb_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_nt_blastdb_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_nr_code():
                script = get_nr_blastdb_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_nr_blastdb_script(cluster_name, current_run_dir)

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                script = get_viridiplantae_nucleotide_gi_gilist_script()
                log.write('Building the process script {0} ...\n'.format(script))
                (OK, error_list) = build_viridiplantae_nucleotide_gi_gilist_script(cluster_name, current_run_dir)
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                script = get_viridiplantae_protein_gi_gilist_script()
                log.write('Building the process script {0} ...\n'.format(script))
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
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(script, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(script)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
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
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_gymno_01_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_proteome_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_dicots_04_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_proteome_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_monocots_04_proteome_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_refseq_plant_code():
                starter = get_refseq_plant_proteome_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_refseq_plant_proteome_starter(current_run_dir)

        # processes to download functional annotations from a genomic database server
        elif process_type == xlib.get_toa_type_download_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                starter = get_basic_data_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_basic_data_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                starter = get_gymno_01_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_gymno_01_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_dicots_04_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_monocots_04_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                starter = get_gene_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_gene_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                starter = get_interpro_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_interpro_download_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                starter = get_go_download_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_go_download_starter(current_run_dir)

        # processes to load data of a genomic database into TOA database
        elif process_type == xlib.get_toa_type_load_data():
            if genomic_database == xlib.get_toa_data_basic_data_code():
                starter = get_basic_data_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_basic_data_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gymno_01_code():
                starter = get_gymno_01_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_gymno_01_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_dicots_04_code():
                starter = get_dicots_04_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_dicots_04_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_monocots_04_code():
                starter = get_monocots_04_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_monocots_04_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_gene_code():
                starter = get_gene_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_gene_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_interpro_code():
                starter = get_interpro_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_interpro_load_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_go_code():
                starter = get_go_load_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_go_load_starter(current_run_dir)

        # processes to build BLAST databases
        elif process_type == xlib.get_toa_type_build_blastdb():
            if genomic_database == xlib.get_toa_data_nt_code():
                starter = get_nt_blastdb_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_nt_blastdb_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_nr_code():
                starter = get_nr_blastdb_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_nr_blastdb_starter(current_run_dir)

        # processes to build GeneId identifier list
        elif process_type == xlib.get_toa_type_build_gilist():
            if genomic_database == xlib.get_toa_data_viridiplantae_nucleotide_gi_code():
                starter = get_viridiplantae_nucleotide_gi_gilist_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_viridiplantae_nucleotide_gi_gilist_starter(current_run_dir)
            elif genomic_database == xlib.get_toa_data_viridiplantae_protein_gi_code():
                starter = get_viridiplantae_protein_gi_gilist_starter()
                log.write('Building the process starter {0} ...\n'.format(starter))
                (OK, error_list) = build_viridiplantae_protein_gi_gilist_starter(current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(starter, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(starter))
        (OK, error_list) = xssh.put_file(sftp_client, starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(starter)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(starter))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(starter)))
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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$EC_DIR" ]; then mkdir --parents $EC_DIR; fi'))
                script_file_id.write( '{0}\n'.format('if [ ! -d "$KEGG_DIR" ]; then mkdir --parents $KEGG_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_basic_data'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Enzyme Commission (EC) ids ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document  $EC_IDS_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $EC_IDS_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading KEGG ids ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $KEGG_IDS_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $KEGG_IDS_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_basic_data'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_basic_data_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_basic_data_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_basic_data_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_basic_data_download_script():
    '''
    Get the script path in the local computer to download other basic data.
    '''

    # assign the script path
    basic_data_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_basic_data_code())

    # return the script path
    return basic_data_download_script

#-------------------------------------------------------------------------------

def get_basic_data_download_starter():
    '''
    Get the script path in the local computer to download other basic data.
    '''

    # assign the starter path
    basic_data_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_basic_data_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$TOA_DB_DIR" ]; then mkdir --parents $TOA_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function load_basic_data'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading basic data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-basic-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --datasets=$DATASET_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --species=$SPECIES_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --ecids=$EC_IDS_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --keggids=$KEGG_IDS_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_basic_data'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_basic_data_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_basic_data_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_basic_data_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_basic_data_load_script():
    '''
    Get the script path in the local computer to load basic data into TOA database.
    '''

    # assign the script path
    basic_data_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_basic_data_code())

    # return the script path
    return basic_data_load_script

#-------------------------------------------------------------------------------

def get_basic_data_load_starter():
    '''
    Get the starter path in the local computer to load basic data into TOA database.
    '''

    # assign the starter path
    basic_data_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_basic_data_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$GYMNO_01_PROTEOME_DB_DIR" ]; then mkdir --parents $GYMNO_01_PROTEOME_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_gymno01_proteome'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GYMNO_01_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GYMNO_01_PROTEOME_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Decompressing proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        gzip --decompress --force $GYMNO_01_PROTEOME_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gzip $RC; fi'))
                script_file_id.write( '{0}\n'.format("    GYMNO_01_PROTEOME_FILE=`echo $GYMNO_01_PROTEOME_FILE | sed 's/.gz//g'`"))
                script_file_id.write( '{0}\n'.format('    echo "File is decompressed."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating BLAST database ..."'))
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        makeblastdb \\'))
                script_file_id.write( '{0}\n'.format('            -title $GYMNO_01_PROTEOME_DB_NAME \\'))
                script_file_id.write( '{0}\n'.format('            -dbtype prot \\'))
                script_file_id.write( '{0}\n'.format('            -input_type fasta \\'))
                script_file_id.write( '{0}\n'.format('            -in $GYMNO_01_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            -out $GYMNO_01_PROTEOME_DB_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "BLAST database is generated."'))
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
                script_file_id.write( '{0}\n'.format('build_gymno01_proteome'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_gymno_01_proteome_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gymno_01_proteome_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gymno_01_proteome_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_proteome_script():
    '''
    Get the script path in the local computer to build the Gymno PLAZA 1.0 proteome.
    '''

    # assign the script path
    gymno_01_proteome_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_gymno_01_code())

    # return the script path
    return gymno_01_proteome_script

#-------------------------------------------------------------------------------

def get_gymno_01_proteome_starter():
    '''
    Get the starter path in the local computer to build the Gymno PLAZA 1.0 proteome.
    '''

    # assign the starter path
    gymno_01_proteome_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_gymno_01_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$GYMNO_01_GENEDESC_DIR" ]; then mkdir --parents $GYMNO_01_GENEDESC_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_gymno01_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading gene description files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=0 \\'))
                script_file_id.write( '{0}\n'.format('            --no-host-directories \\'))
                script_file_id.write( '{0}\n'.format('            --cut-dirs=4 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$GYMNO_01_GENEDESC_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$GYMNO_01_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $GYMNO_01_GENEDESC_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading InterPro file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GYMNO_01_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GYMNO_01_INTERPRO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GYMNO_01_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GYMNO_01_GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GYMNO_01_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GYMNO_01_MAPMAN_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_gymno01_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_gymno_01_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gymno_01_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gymno_01_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_download_script():
    '''
    Get the script path in the local computer to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # assign the script path
    gymno_01_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_gymno_01_code())

    # return the script path
    return gymno_01_download_script

#-------------------------------------------------------------------------------

def get_gymno_01_download_starter():
    '''
    Get the script path in the local computer to download the Gymno PLAZA 1.0 functional annotation.
    '''

    # assign the starter path
    gymno_01_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_gymno_01_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_gymno01_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-plaza-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --dataset=gymno_01 \\'))
                script_file_id.write( '{0}\n'.format('            --species=all \\'))
                script_file_id.write( '{0}\n'.format('            --genedesc=$GYMNO_01_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            --interpro=$GYMNO_01_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --go=$GYMNO_01_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --mapman=$GYMNO_01_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_gymno01_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_gymno_01_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gymno_01_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gymno_01_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gymno_01_load_script():
    '''
    Get the script path in the local computer to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # assign the script path
    gymno_01_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_gymno_01_code())

    # return the script path
    return gymno_01_load_script

#-------------------------------------------------------------------------------

def get_gymno_01_load_starter():
    '''
    Get the starter path in the local computer to load Gymno PLAZA 1.0 functional annotation into TOA database.
    '''

    # assign the starter path
    gymno_01_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_gymno_01_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$DICOTS_04_PROTEOME_DB_DIR" ]; then mkdir --parents $DICOTS_04_PROTEOME_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_dicots04_proteome'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $DICOTS_04_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $DICOTS_04_PROTEOME_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Decompressing proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        gzip --decompress --force $DICOTS_04_PROTEOME_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gzip $RC; fi'))
                script_file_id.write( '{0}\n'.format("    DICOTS_04_PROTEOME_FILE=`echo $DICOTS_04_PROTEOME_FILE | sed 's/.gz//g'`"))
                script_file_id.write( '{0}\n'.format('    echo "File is decompressed."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating BLAST database ..."'))
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        makeblastdb \\'))
                script_file_id.write( '{0}\n'.format('            -title $DICOTS_04_PROTEOME_DB_NAME \\'))
                script_file_id.write( '{0}\n'.format('            -dbtype prot \\'))
                script_file_id.write( '{0}\n'.format('            -input_type fasta \\'))
                script_file_id.write( '{0}\n'.format('            -in $DICOTS_04_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            -out $DICOTS_04_PROTEOME_DB_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "BLAST database is generated."'))
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
                script_file_id.write( '{0}\n'.format('build_dicots04_proteome'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_dicots_04_proteome_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_dicots_04_proteome_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_dicots_04_proteome_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_proteome_script():
    '''
    Get the script path in the local computer to build the Dicots PLAZA 4.0 proteome.
    '''

    # assign the script path
    dicots_04_proteome_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_dicots_04_code())

    # return the script path
    return dicots_04_proteome_script

#-------------------------------------------------------------------------------

def get_dicots_04_proteome_starter():
    '''
    Get the starter path in the local computer to build the Dicots PLAZA 4.0 proteome.
    '''

    # assign the starter path
    dicots_04_proteome_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_dicots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$DICOTS_04_GENEDESC_DIR" ]; then mkdir --parents $DICOTS_04_GENEDESC_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_dicots04_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading gene description files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=0 \\'))
                script_file_id.write( '{0}\n'.format('            --no-host-directories \\'))
                script_file_id.write( '{0}\n'.format('            --cut-dirs=4 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$DICOTS_04_GENEDESC_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$DICOTS_04_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $DICOTS_04_GENEDESC_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading InterPro file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $DICOTS_04_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $DICOTS_04_INTERPRO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $DICOTS_04_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $DICOTS_04_GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $DICOTS_04_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $DICOTS_04_MAPMAN_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_dicots04_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_dicots_04_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_dicots_04_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_dicots_04_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_download_script():
    '''
    Get the script path in the local computer to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # assign the script path
    dicots_04_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_dicots_04_code())

    # return the script path
    return dicots_04_download_script

#-------------------------------------------------------------------------------

def get_dicots_04_download_starter():
    '''
    Get the script path in the local computer to download the Dicots PLAZA 4.0 functional annotation.
    '''

    # assign the starter path
    dicots_04_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_dicots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_dicots04_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-plaza-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --dataset=dicots_04 \\'))
                script_file_id.write( '{0}\n'.format('            --species=all \\'))
                script_file_id.write( '{0}\n'.format('            --genedesc=$DICOTS_04_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            --interpro=$DICOTS_04_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --go=$DICOTS_04_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --mapman=$DICOTS_04_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_dicots04_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_dicots_04_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_dicots_04_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_dicots_04_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_dicots_04_load_script():
    '''
    Get the script path in the local computer to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the script path
    dicots_04_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_dicots_04_code())

    # return the script path
    return dicots_04_load_script

#-------------------------------------------------------------------------------

def get_dicots_04_load_starter():
    '''
    Get the starter path in the local computer to load Dicots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the starter path
    dicots_04_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_dicots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$MONOCOTS_04_PROTEOME_DB_DIR" ]; then mkdir --parents $MONOCOTS_04_PROTEOME_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_monocots04_proteome'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $MONOCOTS_04_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $MONOCOTS_04_PROTEOME_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Decompressing proteome file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        gzip --decompress --force $MONOCOTS_04_PROTEOME_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gzip $RC; fi'))
                script_file_id.write( '{0}\n'.format("    MONOCOTS_04_PROTEOME_FILE=`echo $MONOCOTS_04_PROTEOME_FILE | sed 's/.gz//g'`"))
                script_file_id.write( '{0}\n'.format('    echo "File is decompressed."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating BLAST database ..."'))
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        makeblastdb \\'))
                script_file_id.write( '{0}\n'.format('            -title $MONOCOTS_04_PROTEOME_DB_NAME \\'))
                script_file_id.write( '{0}\n'.format('            -dbtype prot \\'))
                script_file_id.write( '{0}\n'.format('            -input_type fasta \\'))
                script_file_id.write( '{0}\n'.format('            -in $MONOCOTS_04_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            -out $MONOCOTS_04_PROTEOME_DB_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "BLAST database is generated."'))
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
                script_file_id.write( '{0}\n'.format('build_monocots04_proteome'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_monocots_04_proteome_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_monocots_04_proteome_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_monocots_04_proteome_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_proteome_script():
    '''
    Get the script path in the local computer to build the Monocots PLAZA 4.0 proteome.
    '''

    # assign the script path
    monocots_04_proteome_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_monocots_04_code())

    # return the script path
    return monocots_04_proteome_script

#-------------------------------------------------------------------------------

def get_monocots_04_proteome_starter():
    '''
    Get the starter path in the local computer to build the Monocots PLAZA 4.0 proteome.
    '''

    # assign the starter path
    monocots_04_proteome_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_monocots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$MONOCOTS_04_GENEDESC_DIR" ]; then mkdir --parents $MONOCOTS_04_GENEDESC_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_monocots04_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading gene description files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=0 \\'))
                script_file_id.write( '{0}\n'.format('            --no-host-directories \\'))
                script_file_id.write( '{0}\n'.format('            --cut-dirs=4 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$MONOCOTS_04_GENEDESC_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$MONOCOTS_04_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $MONOCOTS_04_GENEDESC_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading InterPro file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $MONOCOTS_04_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $MONOCOTS_04_INTERPRO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $MONOCOTS_04_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $MONOCOTS_04_GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $MONOCOTS_04_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $MONOCOTS_04_MAPMAN_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_monocots04_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_monocots_04_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_monocots_04_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_monocots_04_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_download_script():
    '''
    Get the script path in the local computer to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # assign the script path
    monocots_04_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_monocots_04_code())

    # return the script path
    return monocots_04_download_script

#-------------------------------------------------------------------------------

def get_monocots_04_download_starter():
    '''
    Get the script path in the local computer to download the Monocots PLAZA 4.0 functional annotation.
    '''

    # assign the starter path
    monocots_04_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_monocots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_monocots04_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-plaza-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --dataset=monocots_04 \\'))
                script_file_id.write( '{0}\n'.format('            --species=all \\'))
                script_file_id.write( '{0}\n'.format('            --genedesc=$MONOCOTS_04_GENEDESC_DIR \\'))
                script_file_id.write( '{0}\n'.format('            --interpro=$MONOCOTS_04_INTERPRO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --go=$MONOCOTS_04_GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --mapman=$MONOCOTS_04_MAPMAN_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-plaza-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_monocots04_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_monocots_04_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_monocots_04_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_monocots_04_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_monocots_04_load_script():
    '''
    Get the script path in the local computer to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the script path
    monocots_04_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_monocots_04_code())

    # return the script path
    return monocots_04_load_script

#-------------------------------------------------------------------------------

def get_monocots_04_load_starter():
    '''
    Get the starter path in the local computer to load Monocots PLAZA 4.0 functional annotation into TOA database.
    '''

    # assign the starter path
    monocots_04_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_monocots_04_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$REFSEQ_PLANT_LOCAL" ]; then mkdir --parents $REFSEQ_PLANT_LOCAL; fi'))
                script_file_id.write( '{0}\n'.format('if [ ! -d "$REFSEQ_PLANT_PROTEOME_DB_DIR" ]; then mkdir --parents $REFSEQ_PLANT_PROTEOME_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_refseqplant_proteome'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading protein FASTA files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=1 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$REFSEQ_PROTEIN_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$NCBI_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $REFSEQ_PLANT_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Building proteome file ..."'))
                script_file_id.write( '{0}\n'.format('    > $REFSEQ_PLANT_PROTEOME_FILE'))
                script_file_id.write( '{0}\n'.format("    ls `echo $REFSEQ_PLANT_LOCAL/'*'$REFSEQ_PROTEIN_FILE_PATTERN` > $REFSEQ_PLANT_FILE_LIST"))
                script_file_id.write( '{0}\n'.format('    while read FILE_GZ; do'))
                script_file_id.write( '{0}\n'.format("        FILE_FASTA=`echo $FILE_GZ | sed 's/.gz//g'`"))
                script_file_id.write( '{0}\n'.format('        gzip --decompress --force  $FILE_GZ'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error gzip $RC; fi'))
                script_file_id.write( '{0}\n'.format('        cat $FILE_FASTA >> $REFSEQ_PLANT_PROTEOME_FILE'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error cat $RC; fi'))
                script_file_id.write( '{0}\n'.format('        rm -f $FILE_FASTA'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('    done < $REFSEQ_PLANT_FILE_LIST'))
                script_file_id.write( '{0}\n'.format('    echo "Proteome is buit."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating BLAST database ..."'))
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        makeblastdb \\'))
                script_file_id.write( '{0}\n'.format('            -title $REFSEQ_PLANT_PROTEOME_DB_NAME \\'))
                script_file_id.write( '{0}\n'.format('            -dbtype prot \\'))
                script_file_id.write( '{0}\n'.format('            -input_type fasta \\'))
                script_file_id.write( '{0}\n'.format('            -in $REFSEQ_PLANT_PROTEOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('            -out $REFSEQ_PLANT_PROTEOME_DB_FILE'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-basic-data.py $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "BLAST database is generated."'))
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
                script_file_id.write( '{0}\n'.format('build_refseqplant_proteome'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_refseq_plant_proteome_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_refseq_plant_proteome_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_refseq_plant_proteome_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_refseq_plant_proteome_script():
    '''
    Get the script path in the local computer to build the NCBI RefSeq Plant proteome.
    '''

    # assign the script path
    refseq_plant_proteome_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_refseq_plant_code())

    # return the script path
    return refseq_plant_proteome_script

#-------------------------------------------------------------------------------

def get_refseq_plant_proteome_starter():
    '''
    Get the starter path in the local computer to build the NCBI RefSeq Plant proteome.
    '''

    # assign the starter path
    refseq_plant_proteome_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_proteome_refseq_plant_code())

    # return the starter path
    return refseq_plant_proteome_starter

#-------------------------------------------------------------------------------

def build_nt_blastdb_script(cluster_name, current_run_dir):
    '''
    Build the script to build BLAST database NT.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nt_blastdb_script())):
            os.makedirs(os.path.dirname(get_nt_blastdb_script()))
        with open(get_nt_blastdb_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$NT_DB_DIR" ]; then mkdir --parents $NT_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_database_nt'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading nt database files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=1 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$NT_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$NCBI_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $BLAST_DATABASES_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Decompressing nt database files ..."'))
                script_file_id.write( '{0}\n'.format('    ls `echo $BLAST_DATABASES_LOCAL/$NT_FILE_PATTERN` > $NT_FILE_LIST'))
                script_file_id.write( '{0}\n'.format('    while read NT_FILE; do'))
                script_file_id.write( '{0}\n'.format('        tar --extract --gzip --file=$NT_FILE --directory=$NT_DB_DIR'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tar $RC; fi'))
                script_file_id.write( '{0}\n'.format('        rm -f $NT_FILE'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('    done < $NT_FILE_LIST'))
                script_file_id.write( '{0}\n'.format('    echo "Files are decompressed."'))
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
                process_name = f'TOA - {xlib.get_toa_process_blastdb_nt_name()}'
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
                script_file_id.write( '{0}\n'.format('build_database_nt'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_nt_blastdb_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nt_blastdb_starter(current_run_dir):
    '''
    Build the starter of the script to build BLAST database NT.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nt_blastdb_starter())):
            os.makedirs(os.path.dirname(get_nt_blastdb_starter()))
        with open(get_nt_blastdb_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_nt_blastdb_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_nt_blastdb_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nt_blastdb_script():
    '''
    Get the script path in the local computer to build BLAST database NT.
    '''

    # assign the script path
    nt_blastdb_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_blastdb_nt_code())

    # return the script path
    return nt_blastdb_script

#-------------------------------------------------------------------------------

def get_nt_blastdb_starter():
    '''
    Get the starter path in the local computer to build BLAST database NT.
    '''

    # assign the starter path
    nt_blastdb_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_blastdb_nt_code())

    # return the starter path
    return nt_blastdb_starter

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_viridiplantae_nucleotide_gilist'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Building Viridiplantae nucleotide GI list ..."'))
                script_file_id.write(f'        source activate {xlib.get_entrez_direct_code()}\n')
                script_file_id.write( '{0}\n'.format('    esearch -db nucleotide -query "Viridiplantae[Organism]" | efetch -format uid > $NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error entrez-direct $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "List is built."'))
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
                script_file_id.write( '{0}\n'.format('build_viridiplantae_nucleotide_gilist'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_viridiplantae_nucleotide_gi_gilist_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_viridiplantae_nucleotide_gi_gilist_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_viridiplantae_nucleotide_gi_gilist_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_viridiplantae_nucleotide_gi_gilist_script():
    '''
    Get the script path in the local computer to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # assign the script path
    viridiplantae_nucleotide_gi_gilist_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code())

    # return the script path
    return viridiplantae_nucleotide_gi_gilist_script

#-------------------------------------------------------------------------------

def get_viridiplantae_nucleotide_gi_gilist_starter():
    '''
    Get the starter path in the local computer to build the NCBI Nucleotide GenInfo viridiplantae identifier list.
    '''

    # assign the starter path
    viridiplantae_nucleotide_gi_gilist_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_gilist_viridiplantae_nucleotide_gi_code())

    # return the starter path
    return viridiplantae_nucleotide_gi_gilist_starter

#-------------------------------------------------------------------------------

def build_nr_blastdb_script(cluster_name, current_run_dir):
    '''
    Build the script to build BLAST database NR.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of TOA configuration.
    toa_config_dict = get_toa_config_dict()

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nr_blastdb_script())):
            os.makedirs(os.path.dirname(get_nr_blastdb_script()))
        with open(get_nr_blastdb_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$NR_DB_DIR" ]; then mkdir --parents $NR_DB_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_database_nr'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading nr database files ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --recursive \\'))
                script_file_id.write( '{0}\n'.format('            --level=1 \\'))
                script_file_id.write( '{0}\n'.format('            --accept=$NR_FILE_PATTERN \\'))
                script_file_id.write( '{0}\n'.format('            --directory-prefix=$NCBI_DIR \\'))
                script_file_id.write( '{0}\n'.format('            $BLAST_DATABASES_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Files are downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Decompressing nr database files ..."'))
                script_file_id.write( '{0}\n'.format('    ls `echo $BLAST_DATABASES_LOCAL/$NR_FILE_PATTERN` > $NR_FILE_LIST'))
                script_file_id.write( '{0}\n'.format('    while read NR_FILE; do'))
                script_file_id.write( '{0}\n'.format('        tar --extract --gzip --file=$NR_FILE --directory=$NR_DB_DIR'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tar $RC; fi'))
                script_file_id.write( '{0}\n'.format('        rm -f $NR_FILE'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('    done < $NR_FILE_LIST'))
                script_file_id.write( '{0}\n'.format('    echo "Files are decompressed."'))
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
                process_name = f'TOA - {xlib.get_toa_process_blastdb_nr_name()}'
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
                script_file_id.write( '{0}\n'.format('build_database_nr'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_nr_blastdb_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_nr_blastdb_starter(current_run_dir):
    '''
    Build the starter of the script to build BLAST database NR.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starter
    try:
        if not os.path.exists(os.path.dirname(get_nr_blastdb_starter())):
            os.makedirs(os.path.dirname(get_nr_blastdb_starter()))
        with open(get_nr_blastdb_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_nr_blastdb_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_nr_blastdb_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nr_blastdb_script():
    '''
    Get the script path in the local computer to build BLAST database NR.
    '''

    # assign the script path
    nr_blastdb_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_blastdb_nr_code())

    # return the script path
    return nr_blastdb_script

#-------------------------------------------------------------------------------

def get_nr_blastdb_starter():
    '''
    Get the starter path in the local computer to build BLAST database NR.
    '''

    # assign the starter path
    nr_blastdb_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_blastdb_nr_code())

    # return the starter path
    return nr_blastdb_starter

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function build_viridiplantae_protein_gilist'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Building Viridiplantae protein GI list ..."'))
                script_file_id.write(f'        source activate {xlib.get_entrez_direct_anaconda_code()}\n')
                script_file_id.write( '{0}\n'.format('    esearch -db protein -query "Viridiplantae[Organism]" | efetch -format uid > $PROTEIN_VIRIDIPLANTAE_GI_LIST'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error entrez-direct $RC; fi'))
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '{0}\n'.format('    echo "List is built."'))
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
                process_name = f'TOA - {get_toa_process_gilist_viridiplantae_nucleotide_gi_name()}'
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
                script_file_id.write( '{0}\n'.format('build_viridiplantae_protein_gilist'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_viridiplantae_protein_gi_gilist_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_viridiplantae_protein_gi_gilist_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_viridiplantae_protein_gi_gilist_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_viridiplantae_protein_gi_gilist_script():
    '''
    Get the script path in the local computer to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # assign the script path
    viridiplantae_protein_gi_gilist_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_gilist_viridiplantae_protein_gi_code())

    # return the script path
    return viridiplantae_protein_gi_gilist_script

#-------------------------------------------------------------------------------

def get_viridiplantae_protein_gi_gilist_starter():
    '''
    Get the starter path in the local computer to build the NCBI Protein GenInfo viridiplantae identifier list.
    '''

    # assign the starter path
    viridiplantae_protein_gi_gilist_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_gilist_viridiplantae_protein_gi_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$NCBI_DIR" ]; then mkdir --parents $NCBI_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_gene_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene to RefSeq file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GENE_GENE2REFSEQ_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GENE_GENE2REFSEQ_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading Gene Ontology file ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GENE_GENE2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GENE_GENE2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_gene_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_gene_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gene_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gene_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gene_download_script():
    '''
    Get the script path in the local computer to download the NCBI Gene functional annotation.
    '''

    # assign the script path
    gene_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_gene_code())

    # return the script path
    return gene_download_script

#-------------------------------------------------------------------------------

def get_gene_download_starter():
    '''
    Get the script path in the local computer to download the NCBI Gene functional annotation.
    '''

    # assign the starter path
    gene_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_gene_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_gene_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-ncbi-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --dataset=gene \\'))
                script_file_id.write( '{0}\n'.format('            --gene2refseq=$GENE_GENE2REFSEQ_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --gene2go=$GENE_GENE2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-ncbi-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_gene_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_gene_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gene_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gene_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gene_load_script():
    '''
    Get the script path in the local computer to load NCBI Gene functional annotation into TOA database.
    '''

    # assign the script path
    gene_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_gene_code())

    # return the script path
    return gene_load_script

#-------------------------------------------------------------------------------

def get_gene_load_starter():
    '''
    Get the starter path in the local computer to load NCBI Gene functional annotation into TOA database.
    '''

    # assign the starter path
    gene_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_gene_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$INTERPRO_DIR" ]; then mkdir --parents $INTERPRO_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_interpro_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading file of mappings of InterPro entries to Gene Ontology terms ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $INTERPRO_INTERPRO2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $INTERPRO_INTERPRO2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_interpro_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_interpro_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_interpro_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_interpro_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_interpro_download_script():
    '''
    Get the script path in the local computer to download the InterPro functional annotation.
    '''

    # assign the script path
    interpro_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_interpro_code())

    # return the script path
    return interpro_download_script

#-------------------------------------------------------------------------------

def get_interpro_download_starter():
    '''
    Get the script path in the local computer to download the InterPro functional annotation.
    '''

    # assign the starter path
    interpro_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_interpro_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_interpro_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-interpro-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --interpro2go=$INTERPRO_INTERPRO2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-interpro-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_interpro_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_interpro_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_interpro_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_interpro_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_interpro_load_script():
    '''
    Get the script path in the local computer to load InterPro functional annotation into TOA database.
    '''

    # assign the script path
    interpro_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_interpro_code())

    # return the script path
    return interpro_load_script

#-------------------------------------------------------------------------------

def get_interpro_load_starter():
    '''
    Get the starter path in the local computer to load InterPro functional annotation into TOA database.
    '''

    # assign the starter path
    interpro_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_interpro_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
                script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
                script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
                script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('if [ ! -d "$GO_DIR" ]; then mkdir --parents $GO_DIR; fi'))
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
                script_file_id.write( '{0}\n'.format('function download_go_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading ontology ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GO_ONTOLOGY_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GO_ONTOLOGY_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading GO to Enzyme Commission (EC) ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document  $GO_EC2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GO_EC2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading GO to KEGG ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GO_KEGG2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GO_KEGG2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading GO to MetaCyc ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GO_METACYC2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GO_METACYC2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Downloading GO to InterPro ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        wget \\'))
                script_file_id.write( '{0}\n'.format('            --quiet \\'))
                script_file_id.write( '{0}\n'.format('            --output-document $GO_INTERPRO2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            $GO_INTERPRO2GO_FTP'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error wget $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "File is downloaded."'))
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
                script_file_id.write( '{0}\n'.format('download_go_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_go_download_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_go_download_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_go_download_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_go_download_script():
    '''
    Get the script path in the local computer to download the Gene Ontology functional annotation.
    '''

    # assign the script path
    go_download_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_go_code())

    # return the script path
    return go_download_script

#-------------------------------------------------------------------------------

def get_go_download_starter():
    '''
    Get the script path in the local computer to download the Gene Ontology functional annotation.
    '''

    # assign the starter path
    go_download_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_download_go_code())

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
                script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
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
                script_file_id.write( '{0}\n'.format('function load_go_functional_annotation'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Loading functional annotation data into TOA database ..."'))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('        load-go-data.py \\'))
                script_file_id.write( '{0}\n'.format('            --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('            --ontology=$GO_ONTOLOGY_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --ec2go=$GO_EC2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --kegg2go=$GO_KEGG2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --metacyc2go=$GO_METACYC2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --interpro2go=$GO_INTERPRO2GO_FILE \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('            --trace=N'))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error load-go-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    echo "Data are loaded."'))
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
                script_file_id.write( '{0}\n'.format('load_go_functional_annotation'))
                script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_go_load_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_go_load_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_go_load_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_go_load_script():
    '''
    Get the script path in the local computer to load Gene Ontology functional annotation into TOA database.
    '''

    # assign the script path
    go_load_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_go_code())

    # return the script path
    return go_load_script

#-------------------------------------------------------------------------------

def get_go_load_starter():
    '''
    Get the starter path in the local computer to load Gene Ontology functional annotation into TOA database.
    '''

    # assign the starter path
    go_load_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_load_go_code())

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
    except:
        dicots_04_order = 0
    try:
        gymno_01_order = database_list.index('gymno_01') + 1
    except:
        gymno_01_order = 0
    try:
        monocots_04_order = database_list.index('monocots_04') + 1
    except:
        monocots_04_order = 0
    try:
        refseq_plant_order = database_list.index('refseq_plant') + 1
    except:
        refseq_plant_order = 0
    if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
        try:
            nt_viridiplantae_order = database_list.index('nt_viridiplantae') + 1
        except:
            nt_viridiplantae_order = 0
        try:
            nt_complete_order = database_list.index('nt_complete') + 1
        except:
            nt_complete_order = 0
    elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
        try:
            nr_viridiplantae_order = database_list.index('nr_viridiplantae') + 1
        except:
            nr_viridiplantae_order = 0
        try:
            nr_complete_order = database_list.index('nr_complete') + 1
        except:
            nr_complete_order = 0

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
            file_id.write( '{0}\n'.format('# The assembly files have to be located in:'))
            file_id.write( '{0}\n'.format('#    {0}/experiment_id/assembly_dataset_id  (when the assembly origin is NGSCLOUD)'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('#    {0}/reference_dataset_id  (when the assembly origin is EXTERNAL)'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, database_dataset_id, assembly_dataset_id, reference_dataset_id and trasncriptome_file names are fixed in the identification section.'))
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# You can consult the parameters of nucleotide pipeline (TOA package) and their meaning in "https://github.com/GGFHF/TOA".'))
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format('assembly_origin = {0}'.format(assembly_origin), '# assembly origin: {0}'.format(get_assembly_origin_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification; NONE when the assembly origin is EXTERNAL'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}; NONE when the assembly origin is EXTERNAL'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification; NONE when the assembly origin is EXTERNAL'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification; NONE when the assembly origin is NGSCLOUD'))
            file_id.write( '{0:<50} {1}\n'.format('transcriptome_file = {0}'.format(transcriptome_file), '# transcriptome file name; NONE when the assembly origin is NGSCLOUD'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the database parameters'))
            file_id.write( '{0}\n'.format('[database parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('gymno_01 = {0}'.format(gymno_01_order), '# order of Gymno PLAZA 1.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format('dicots_04 = {0}'.format(dicots_04_order), '# order of Dicots PLAZA 4.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format('monocots_04 = {0}'.format(monocots_04_order), '# order of Monocots PLAZA 4.0 in the annotation; 0 if it is not used'))
            file_id.write( '{0:<50} {1}\n'.format('refseq_plant = {0}'.format(refseq_plant_order), '# order of NCBI RefSeq Plant in the annotation; 0 if it is not used'))
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                file_id.write( '{0:<50} {1}\n'.format('nt_viridiplantae = {0}'.format(nt_viridiplantae_order), '# order of NCBI BLAST database NT (Viridiplantae) in the annotation; 0 if it is not used'))
                file_id.write( '{0:<50} {1}\n'.format('nt_complete = {0}'.format(nt_complete_order), '# order of NCBI BLAST database NT (complete) in the annotation, it has to be the last; 0 if it is not used'))
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                file_id.write( '{0:<50} {1}\n'.format('nr_viridiplantae = {0}'.format(nr_viridiplantae_order), '# order of NCBI BLAST database NR (Viridiplantae) in the annotation; 0 if it is not used'))
                file_id.write( '{0:<50} {1}\n'.format('nr_complete = {0}'.format(nr_complete_order), '# order of NCBI BLAST database NR (complete) in the annotation, it has to be the last; 0 if it is not used'))
                pass
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the BLAST parameters'))
            file_id.write( '{0}\n'.format('[BLAST parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('thread_number = 16', '# threads number'))
            file_id.write( '{0:<50} {1}\n'.format('e_value = 1E-6', '# expectation value (E-value) threshold for saving hits'))
            file_id.write( '{0:<50} {1}\n'.format('max_target_seqs = 20', '# maximum number of aligned sequences to keep'))
            file_id.write( '{0:<50} {1}\n'.format('max_hsps = 999999', '# maximum number of HSPs per subject sequence to save for each query'))
            file_id.write( '{0:<50} {1}\n'.format('qcov_hsp_perc = 0.0', '# alignments below the specified query coverage per HSPs are removed'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(config_file))
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
                error_list.append('*** ERROR: the key "assembly_origin" has to be {0}.'.format(get_assembly_origin_code_list_text()))
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
                    error_list.append('*** ERROR: the key "assembly_software" has to be {0} when the assembly origin is NGSCLOUD.'.format(get_assembly_software_code_list_text()))
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
                    error_list.append('*** ERROR: the key "assembly_software" has to start with {0} when the assembly origin is NGSCLOUD.'.format(get_assembly_software_code_list_text()))
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

            # check section "database parameters" - key "nt_viridiplantae"/"nr_viridiplantae"
            is_ok_nx_viridiplantae = False
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                nt_viridiplantae = pipeline_option_dict.get('database parameters', {}).get('nt_viridiplantae', not_found)
                if nt_viridiplantae == not_found:
                    error_list.append('*** ERROR: the key "nt_viridiplantae" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nt_viridiplantae, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nt_viridiplantae" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx_viridiplantae = True
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                nr_viridiplantae = pipeline_option_dict.get('database parameters', {}).get('nr_viridiplantae', not_found)
                if nr_viridiplantae == not_found:
                    error_list.append('*** ERROR: the key "nr_viridiplantae" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nr_viridiplantae, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nr_viridiplantae" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx_viridiplantae = True

            # check section "database parameters" - key "nt_complete"/"nr_complete"
            is_ok_nx_complete = False
            if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                nt_complete = pipeline_option_dict.get('database parameters', {}).get('nt_complete', not_found)
                if nt_complete == not_found:
                    error_list.append('*** ERROR: the key "nt_complete" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nt_complete, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nt_complete" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx_complete = True
            elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                nr_complete = pipeline_option_dict.get('database parameters', {}).get('nr_complete', not_found)
                if nr_complete == not_found:
                    error_list.append('*** ERROR: the key "nr_complete" is not found in the section "database parameters".')
                    OK = False
                elif not xlib.check_int(nr_complete, minimum=0, maximum=6):
                    error_list.append('*** ERROR: the key "nr_complete" has to be an integer number between 0 and 6.')
                    OK = False
                else:
                    is_ok_nx_complete = True

            # check the order of databases
            if is_ok_gymno_01 and is_ok_dicots_04 and is_ok_monocots_04 and is_ok_refseq_plant and is_ok_nx_viridiplantae and is_ok_nx_complete:
                if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
                    database_code_list = get_nucleotide_annotation_database_code_list()
                    database_order_list = [int(dicots_04), int(gymno_01), int(monocots_04), int(refseq_plant), int(nt_viridiplantae), int(nt_complete)]
                    last_database_code = 'nt_complete'
                elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
                    database_code_list = get_aminoacid_annotation_database_code_list()
                    database_order_list = [int(dicots_04), int(gymno_01), int(monocots_04), int(refseq_plant), int(nr_viridiplantae), int(nr_complete)]
                    last_database_code = 'nr_complete'
                (OK2, error_list2) = check_database_order(database_code_list, database_order_list, last_database_code)
                if not OK2:
                    OK = False
                    error_list = error_list + error_list2

        # check section "BLAST parameters"
        if 'BLAST parameters' not in sections_list:
            error_list.append('*** ERROR: the section "BLAST parameters" is not found.')
            OK = False
        else:

            # check section "BLAST parameters" - key "thread_number"
            thread_number = pipeline_option_dict.get('BLAST parameters', {}).get('thread_number', not_found)
            if thread_number == not_found:
                error_list.append('*** ERROR: the key "thread_number" is not found in the section "BLAST parameters".')
                OK = False
            elif not xlib.check_int(thread_number, minimum=1):
                error_list.append('*** ERROR: the key "thread_number" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BLAST parameters" - key "e_value"
            e_value = pipeline_option_dict.get('BLAST parameters', {}).get('e_value', not_found)
            if e_value == not_found:
                error_list.append('*** ERROR: the key "e_value" is not found in the section "BLAST parameters".')
                OK = False
            elif not xlib.check_float(e_value, minimum=0., mne=1E-12):
                error_list.append('*** ERROR: the key "e_value" has to be a float number greater than to 0.0.')
                OK = False

            # check section "BLAST parameters" - key "max_target_seqs"
            max_target_seqs = pipeline_option_dict.get('BLAST parameters', {}).get('max_target_seqs', not_found)
            if max_target_seqs == not_found:
                error_list.append('*** ERROR: the key "max_target_seqs" is not found in the section "BLAST parameters".')
                OK = False
            elif not xlib.check_int(max_target_seqs, minimum=1):
                error_list.append('*** ERROR: the key "max_target_seqsr" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BLAST parameters" - key "max_hsps"
            max_hsps = pipeline_option_dict.get('BLAST parameters', {}).get('max_hsps', not_found)
            if max_hsps == not_found:
                error_list.append('*** ERROR: the key "max_hsps" is not found in the section "BLAST parameters".')
                OK = False
            elif not xlib.check_int(max_hsps, minimum=1):
                error_list.append('*** ERROR: the key "max_hsps" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BLAST parameters" - key "qcov_hsp_perc"
            qcov_hsp_perc = pipeline_option_dict.get('BLAST parameters', {}).get('qcov_hsp_perc', not_found)
            if qcov_hsp_perc == not_found:
                error_list.append('*** ERROR: the key "qcov_hsp_perc" is not found in the section "BLAST parameters".')
                OK = False
            elif not xlib.check_float(qcov_hsp_perc, minimum=0., maximum=100., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "qcov_hsp_perc" has to be a float number greater than or equal to 0.0 and less than 100.0.')
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_toa_process_pipeline_nucleotide_name()))

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
            error_list.append('*** ERROR: the {0} order has to be the last one.'.format(last_database_code))
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
    nucleotide_pipeline_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_toa_process_pipeline_nucleotide_code())

    # return the nucleotide pipeline config file path
    return nucleotide_pipeline_config_file

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_config_file():
    '''
    Get the amino acid pipeline config file path.
    '''

    # assign the amino acid pipeline config file path
    aminoacid_pipeline_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_toa_process_pipeline_aminoacid_code())

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

    # build the sentence to set the cluster PATH to check requirements
    path_sentence = 'export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])

    # get the selected database list
    selected_database_list = get_selected_database_list(pipeline_type)

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the TOA config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_toa_name()))
    OK = os.path.isfile(get_toa_config_file())
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The {0} config file does not exist.\n')
        log.write('Please recreate this file.\n')
        OK = False

    # check the pipeline config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        if pipeline_type == xlib.get_toa_process_pipeline_nucleotide_code():
            log.write('Checking the {0} config file ...\n'.format(xlib.get_toa_process_pipeline_nucleotide_name()))
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            log.write('Checking the {0} config file ...\n'.format(xlib.get_toa_process_pipeline_aminoacid_name()))
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
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_toa_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_toa_name()))
            OK = False

    # check the basic data load in TOA databasep
    if OK:
        check_sentence = 'check-data-load.py --db={0} --group=basic'.format(toa_config_dict['TOA_DB'])
        command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: The basic data load in TOA database is wrong.\n')
            OK = False
        else:
            log.write('... basic data load in TOA database is OK ...\n')

    # check the Gymno PLAZA 1.0 proteome
    if OK:
        if 'gymno_01' in selected_database_list:
            command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['GYMNO_01_PROTEOME_DB_DIR'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Gymno PLAZA 1.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Gymno PLAZA 1.0 proteome is OK ...\n')

    # check the Gymno PLAZA 1.0 data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=gymno_01'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Gymno PLAZA 1.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Gymno PLAZA 1.0 load data in TOA database is OK ...\n')

    # check the Dicots PLAZA 4.0 proteome
    if OK:
        if 'dicots_04' in selected_database_list:
            command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['DICOTS_04_PROTEOME_DB_DIR'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Dicots PLAZA 4.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Dicots PLAZA 4.0 proteome is OK ...\n')

    # check the Dicots PLAZA 4.0 data load in TOA database
    if OK:
        if 'dicots_04' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=dicots_04'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Dicots PLAZA 4.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Dicots PLAZA 4.0 data load in TOA database is OK ...\n')

    # check the Monocots PLAZA 4.0 proteome
    if OK:
        if 'monocots_04' in selected_database_list:
            command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['MONOCOTS_04_PROTEOME_DB_DIR'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Monocots PLAZA 4.0 proteome is not built.\n')
                OK = False
            else:
                log.write('... Monocots PLAZA 4.0 proteome is OK ...\n')

    # check the Monocots PLAZA 4.0 data load in TOA database
    if OK:
        if 'monocots_04' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=monocots_04'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: Monocots PLAZA 4.0 load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... Monocots PLAZA 4.0 data load in TOA database is OK ...\n')

    # check the NCBI RefSeq Plant proteome
    if OK:
        if 'refseq_plant' in selected_database_list:
            command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['REFSEQ_PLANT_PROTEOME_DB_DIR'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI RefSeq Plant proteome is not built.\n')
                OK = False
            else:
                log.write('... NCBI RefSeq Plant proteome is OK ...\n')

    # check the NCBI BLAST database NT
    if OK:
        if 'nt_viridiplantae' in selected_database_list or 'nt_complete' in selected_database_list:
            command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['NT_DB_DIR'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI BLAST database NT is not built.\n')
                OK = False
            else:
                log.write('... NCBI BLAST database NT is OK ...\n')

    # check the NCBI Nucleotide GenInfo viridiplantae identifier list
    if OK:
        if 'nt_viridiplantae' in selected_database_list:
            command = '[ -f {0} ] && echo RC=0 || echo RC=1'.format(toa_config_dict['NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST'])
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI Nucleotide GenInfo viridiplantae identifier list is not built.\n')
                OK = False
            else:
                log.write('... NCBI Nucleotide GenInfo viridiplantae identifier list is OK ...\n')

    # check the NCBI Gene data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=gene'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: NCBI Gene load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... NCBI Gene data load in TOA database is OK ...\n')

    # check the InterPro data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=interpro'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if stdout[len(stdout) - 1] != 'RC=0':
                log.write('*** ERROR: InterPro load in TOA database is wrong.\n')
                OK = False
            else:
                log.write('... InterPro data load in TOA database is OK ...\n')

    # check the Gene Ontology data load in TOA database
    if OK:
        if 'gymno_01' in selected_database_list or 'dicots_04' in selected_database_list or 'monocots_04' in selected_database_list:
            check_sentence = 'check-data-load.py --db={0} --group=go'.format(toa_config_dict['TOA_DB'])
            command = '{0}; {1} && echo RC=0 || echo RC=1'.format(path_sentence, check_sentence)
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
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
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
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
            log.write('Building the process script {0} ...\n'.format(script))
            (OK, error_list) = build_nucleotide_pipeline_script(cluster_name, current_run_dir)

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            script = get_aminoacid_pipeline_script()
            log.write('Building the process script {0} ...\n'.format(script))
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
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(script, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, error_list) = xssh.put_file(sftp_client, script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(script)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(script))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
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
            log.write('Building the process starter {0} ...\n'.format(starter))
            (OK, error_list) = build_nucleotide_pipeline_starter(current_run_dir)

        # amino acid pipelines
        elif pipeline_type == xlib.get_toa_process_pipeline_aminoacid_code():
            starter = get_aminoacid_pipeline_starter()
            log.write('Building the process starter {0} ...\n'.format(starter))
            (OK, error_list) = build_aminoacid_pipeline_starter(current_run_dir)

        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the script starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(starter, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(starter))
        (OK, error_list) = xssh.put_file(sftp_client, starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the script starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(starter)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(starter))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(starter)))
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
    thread_number = pipeline_option_dict['BLAST parameters']['thread_number']
    e_value = pipeline_option_dict['BLAST parameters']['e_value']
    max_target_seqs = pipeline_option_dict['BLAST parameters']['max_target_seqs']
    max_hsps = pipeline_option_dict['BLAST parameters']['max_hsps']
    qcov_hsp_perc = pipeline_option_dict['BLAST parameters']['qcov_hsp_perc']

    # get the all selected database list
    all_database_list = get_selected_database_list(xlib.get_toa_process_pipeline_nucleotide_code())

    # change code "nt_complete" by "nt_remainder"
    if all_database_list[len(all_database_list) - 1] == 'nt_complete':
        all_database_list[len(all_database_list) - 1] = 'nt_remainder'

    # get the plant database list
    plant_database_list = all_database_list.copy()
    if 'nt_remainder' in plant_database_list:
        plant_database_list.remove('nt_remainder')

    # get the database type dictionary
    database_type_dict = {}
    for i in range(len(all_database_list)):
        if all_database_list[i] in ['gymno_01', 'dicots_04', 'monocots_04']:
            database_type_dict[all_database_list[i]] = 'PLAZA'
        elif all_database_list[i] == 'refseq_plant':
            database_type_dict[all_database_list[i]] = 'REFSEQ'
        elif all_database_list[i] in ['nt_viridiplantae', 'nt_remainder']:
            database_type_dict[all_database_list[i]] = 'NT'

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
    for database in all_database_list:
        non_annotation_file_list.append('${0}_NON_ANNOTATED_TRANSCRIPT_FILE'.format(database.upper()))

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_nucleotide_pipeline_script())):
            os.makedirs(os.path.dirname(get_nucleotide_pipeline_script()))
        with open(get_nucleotide_pipeline_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('# transcriptome file'))
            script_file_id.write( '{0}\n'.format('TRANSCRIPTOME_FILE={0}'.format(transcriptome_file)))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# BLAST parameters'))
            script_file_id.write( '{0}\n'.format('NUM_THREADS={0}'.format(thread_number)))
            script_file_id.write( '{0}\n'.format('E_VALUE={0}'.format(e_value)))
            script_file_id.write( '{0}\n'.format('MAX_TARGET_SEQS={0}'.format(max_target_seqs)))
            script_file_id.write( '{0}\n'.format('MAX_HSPS={0}'.format(max_hsps)))
            script_file_id.write( '{0}\n'.format('QCOV_HSP_PERC={0}'.format(qcov_hsp_perc)))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# output directory'))
            script_file_id.write( '{0}\n'.format('OUTPUT_DIR={0}'.format(current_run_dir)))
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
            script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATS_DIR'))
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
            script_file_id.write( '{0}\n'.format('    echo "TRANSCRIPTOME FILE: $TRANSCRIPTOME_FILE"'))
            script_file_id.write( '{0}\n'.format('    echo "ALIGNMENT DATASETS: {0}"'.format(','.join(all_database_list))))
            script_file_id.write( '{0}\n'.format('    echo "NUM_THREADS: $NUM_THREADS"'))
            script_file_id.write( '{0}\n'.format('    echo "E_VALUE: $E_VALUE"'))
            script_file_id.write( '{0}\n'.format('    echo "MAX_TARGET_SEQS: $MAX_TARGET_SEQS"'))
            script_file_id.write( '{0}\n'.format('    echo "MAX_HSPS: $MAX_HSPS"'))
            script_file_id.write( '{0}\n'.format('    echo "QCOV_HSP_PERC: $QCOV_HSP_PERC"'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function reidentify_sequences'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/reidentify_sequences.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "RE-INDENTIFY SEQUENCES OF THE TRANSCRIPTOME FILE"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write( '{0}\n'.format('        echo "Re-identifing sequences ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            reid-fasta-file.py \\'))
            script_file_id.write( '{0}\n'.format('                --fasta=$TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --out=$REIDENTIFIED_TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
            script_file_id.write( '{0}\n'.format('                --trace=N'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error reid-fasta-file.py $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "Sequences are re-identified."'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            for i in range(len(all_database_list)):
                current_code = all_database_list[i]
                previus_code = all_database_list[i - 1] if i > 0 else ''
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function align_transcripts_{0}_proteome'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/align_transcripts_{0}_proteome.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "ALIGNMENT OF TRANSCRIPTS TO {0} PROTEOME"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '{0}\n'.format('        echo "Aligning transcripts ..."'))
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    script_file_id.write( '{0}\n'.format('        export BLASTDB=${0}_PROTEOME_DB_DIR'.format(current_code.upper())))
                elif current_code in ['nt_viridiplantae', 'nt_remainder']:
                    script_file_id.write( '{0}\n'.format('        export BLASTDB=$NT_DB_DIR'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    script_file_id.write( '{0}\n'.format('            blastx \\'))
                elif current_code in ['nt_viridiplantae', 'nt_remainder']:
                    script_file_id.write( '{0}\n'.format('            blastn \\'))
                script_file_id.write( '{0}\n'.format('                -num_threads $NUM_THREADS \\'))
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    script_file_id.write( '{0}\n'.format('                -db ${0}_PROTEOME_DB_NAME \\'.format(current_code.upper())))
                elif current_code in ['nt_viridiplantae', 'nt_remainder']:
                    script_file_id.write( '{0}\n'.format('                -db $NT_DB_NAME \\'))
                if i == 0:
                    script_file_id.write( '{0}\n'.format('                -query $REIDENTIFIED_TRANSCRIPTOME_FILE \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                -query ${0}_NON_ANNOTATED_TRANSCRIPT_FILE \\'.format(previus_code.upper())))
                if current_code == 'nt_viridiplantae':
                    script_file_id.write( '{0}\n'.format('                -gilist $NUCLEOTIDE_VIRIDIPLANTAE_GI_LIST \\'))
                script_file_id.write( '{0}\n'.format('                -evalue $E_VALUE \\'))
                script_file_id.write( '{0}\n'.format('                -max_target_seqs $MAX_TARGET_SEQS \\'))
                script_file_id.write( '{0}\n'.format('                -max_hsps $MAX_HSPS \\'))
                script_file_id.write( '{0}\n'.format('                -qcov_hsp_perc $QCOV_HSP_PERC \\'))
                script_file_id.write( '{0}\n'.format('                -outfmt 5 \\'))
                if current_code == 'nt_remainder':
                    script_file_id.write( '{0}\n'.format('                -out $REIDENTIFIED_NT_REMAINDER_BLAST_XML'))
                else:
                    script_file_id.write( '{0}\n'.format('                -out ${0}_BLAST_XML'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error blastx $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Alignment is done."'))
                script_file_id.write( '{0}\n'.format('        conda deactivate'))
                if current_code == 'nt_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                    script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_NT_REMAINDER_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                    script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --out=$NT_REMAINDER_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                if len(plant_database_list) == 1 and current_code != 'nt_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Creating plant alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            cp ${0}_BLAST_XML $REIDENTIFIED_PLANT_BLAST_XML'.format(current_code.upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error cp $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "File is created."'))
                    script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in merged alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                    script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                    script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --out=$PLANT_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function load_alignment_{0}_proteome'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/load_alignment_{0}_proteome.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "LOAD OF TRANSCRIPT ALIGNMENT TO {0} PROTEOME INTO TOA DATABASE"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Loading alignmnet data ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            load-blast-data.py \\'))
                script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('                --dataset={0} \\'.format(current_code)))
                script_file_id.write( '{0}\n'.format('                --format=5 \\'))
                if current_code == 'nt_remainder':
                    script_file_id.write( '{0}\n'.format('                --blast=$REIDENTIFIED_NT_REMAINDER_BLAST_XML \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --blast=${0}_BLAST_XML \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error load-blast-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Data are loaded."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function annotate_transcripts_{0}'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/annotate_transcripts_{0}.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "ANNOTATION OF TRANSCRIPTS WITH {0}"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Annotating transcripts ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            annotate-sequences.py \\'))
                script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('                --dataset={0} \\'.format(current_code)))
                if i == 0:
                    script_file_id.write( '{0}\n'.format('                --seqs=$REIDENTIFIED_TRANSCRIPTOME_FILE \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --seqs=${0}_NON_ANNOTATED_TRANSCRIPT_FILE \\'.format(previus_code.upper())))
                script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --annotation=${0}_ANNOTATION_FILE \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --nonann=${0}_NON_ANNOTATED_TRANSCRIPT_FILE \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error annotate-sequences.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Annotation is done."'))
                if len(plant_database_list) == 1 and current_code != 'nt_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Creating plant annotation file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            cp ${0}_ANNOTATION_FILE $PLANT_ANNOTATION_FILE'.format(current_code.upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error cp $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "File is created."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            if len(plant_database_list) > 1:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function merge_plant_alignment_files'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/merge_plant_alignment_files.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "MERGER OF PLANT ALIGNMENT FILES"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Merging alignment files ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            merge-xml-files.py \\'))
                plant_blast_xml_list = []
                for database_code in plant_database_list:
                    plant_blast_xml_list.append('${0}_BLAST_XML'.format(database_code.upper()))
                script_file_id.write( '{0}\n'.format('                --list={0} \\'.format(','.join(plant_blast_xml_list))))
                script_file_id.write( '{0}\n'.format('                --relationships=NONE \\'))
                script_file_id.write( '{0}\n'.format('                --mfile=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-xml-files.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in merged alignment file ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --out=$PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function merge_plant_annotation_files'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/merge_plant_annotation_files.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "MERGER OF PLANT ANNOTATION FILES"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                for database_code in plant_database_list:
                    script_file_id.write( '{0}\n'.format('        {0}_ANNOTATION_FILE_TMP=${0}_ANNOTATION_FILE".tmp"'.format(database_code.upper())))
                    script_file_id.write( '{0}\n'.format('        {0}_ANNOTATION_FILE_SORTED=${0}_ANNOTATION_FILE".sorted"'.format(database_code.upper())))
                tmp_file_list = []
                for i in range(len(plant_database_list) - 2):
                    script_file_id.write( '{0}\n'.format('        MERGED_ANNOTATION_FILE_TMP{0}=$MERGED_ANNOTATION_FILE".tmp{0}"'.format(i + 1)))
                    tmp_file_list.append('$MERGED_ANNOTATION_FILE_TMP{0}'.format(i + 1))
                script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                script_file_id.write( '{0}\n'.format('        echo "Merging `basename ${0}_ANNOTATION_FILE` and `basename ${1}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper(), plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            merge-annotation-files.py \\'))
                script_file_id.write( '{0}\n'.format('                --file1=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('                --type1={0} \\'.format(database_type_dict[plant_database_list[0]])))
                script_file_id.write( '{0}\n'.format('                --file2=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('                --type2={0} \\'.format(database_type_dict[plant_database_list[1]])))
                if len(plant_database_list) > 2:
                    script_file_id.write( '{0}\n'.format('                --mfile=$MERGED_ANNOTATION_FILE_TMP1 \\'))
                    script_file_id.write( '{0}\n'.format('                --header=N \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --mfile=$PLANT_ANNOTATION_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                for i in range(2, len(plant_database_list)):
                    script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                    script_file_id.write( '{0}\n'.format('        echo "Adding annotation of `basename ${0}_ANNOTATION_FILE` to the merged annotation file ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            merge-annotation-files.py \\'))
                    script_file_id.write( '{0}\n'.format('                --file1=$MERGED_ANNOTATION_FILE_TMP{} \\'.format(i - 1)))
                    script_file_id.write( '{0}\n'.format('                --type1=MERGER \\'))
                    script_file_id.write( '{0}\n'.format('                --file2=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('                --type2={0} \\'.format(database_type_dict[plant_database_list[i]])))
                    if i < len(plant_database_list) - 1:
                        script_file_id.write( '{0}\n'.format('                --mfile=$MERGED_ANNOTATION_FILE_TMP{0} \\'.format(i)))
                        script_file_id.write( '{0}\n'.format('                --header=N \\'))
                    else:
                        script_file_id.write( '{0}\n'.format('                --mfile=$PLANT_ANNOTATION_FILE \\'))
                        script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                    script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal annotation files ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm {0}'.format(' '.join(tmp_file_list))))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function split_merged_plant_annotation_file'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/split_merged_plant_annotation_file.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "SPLIT OF MERGED PLANT ANNOTATION FILE"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Splitting merged file `basename $PLANT_ANNOTATION_FILE` ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            split-annotation-file.py \\'))
                script_file_id.write( '{0}\n'.format('                --annotation=$PLANT_ANNOTATION_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --type=MERGER \\'))
                script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                script_file_id.write( '{0}\n'.format('                --rnum=250000 \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error split-annotation-file.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "File is splitted."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            if all_database_list[0] != 'nt_remainder':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function purge_transcriptome'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/purge_transcriptome.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "PURGE TRANSCRIPTOME REMOVING NON-ANNOTATED TRANSCRIPTS"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Purging transcriptome files ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            merge-fasta-files.py \\'))
                script_file_id.write( '{0}\n'.format('                --file1=$REIDENTIFIED_TRANSCRIPTOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --file2=${0}_NON_ANNOTATED_TRANSCRIPT_FILE \\'.format(plant_database_list[len(plant_database_list) - 1].upper())))
                script_file_id.write( '{0}\n'.format('                --mfile=$PURGED_TRANSCRIPTOME_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --operation=1LESS2 \\'))
                script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-fasta-files.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Transcriptome is purged."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function calculate_annotation_stats'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/calculate_annotation_stats.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "CALCULATE ANNOTATION STATISTICS"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write( '{0}\n'.format('        echo "Calculating stats ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            calculate-annotation-stats.py \\'))
            script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
            script_file_id.write( '{0}\n'.format('                --transcriptome=$TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --peptides=NONE \\'))
            script_file_id.write( '{0}\n'.format('                --dslist={0} \\'.format(','.join(all_database_list))))
            script_file_id.write( '{0}\n'.format('                --nonannlist={0} \\'.format(','.join(non_annotation_file_list))))
            if len(plant_database_list) > 1:
                script_file_id.write( '{0}\n'.format('                --annotation=$PLANT_ANNOTATION_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --type=MERGER \\'))
            else:
                script_file_id.write( '{0}\n'.format('                --annotation=${0}_ANNOTATION_FILE \\'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('                --type={0} \\'.format(database_type_dict[plant_database_list[0]])))
            script_file_id.write( '{0}\n'.format('                --stats=$ANNOTATION_STATS_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
            script_file_id.write( '{0}\n'.format('                --trace=N'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error calculate-annotation-stats.py $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "Stats are calculated."'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
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
            script_file_id.write( '{0}\n'.format('# re-identify sequences of the transcriptome file'))
            script_file_id.write( '{0}\n'.format('reidentify_sequences'))
            script_file_id.write( '\n')
            for i in range(len(all_database_list)):
                current_code = all_database_list[i]
                previus_code = all_database_list[i - 1] if i > 0 else ''
                if i == 0:
                    script_file_id.write( '{0}\n'.format('# complete transcriptome -> {0}'.format(current_code)))
                else:
                    script_file_id.write( '{0}\n'.format('# transcripts not annotated with {0} -> {1}'.format(previus_code, current_code)))
                script_file_id.write( '{0}\n'.format('align_transcripts_{0}_proteome'.format(current_code)))
                script_file_id.write( '{0}\n'.format('load_alignment_{0}_proteome'.format(current_code)))
                script_file_id.write( '{0}\n'.format('annotate_transcripts_{0}'.format(current_code)))
                script_file_id.write( '\n')
            if len(plant_database_list) > 1:
                script_file_id.write( '{0}\n'.format('# merged plant files'))
                script_file_id.write( '{0}\n'.format('merge_plant_alignment_files'))
                script_file_id.write( '{0}\n'.format('merge_plant_annotation_files'))
                script_file_id.write( '{0}\n'.format('# -- split_merged_plant_annotation_file'))
                script_file_id.write( '\n')
            if all_database_list[0] != 'nt_remainder':
                script_file_id.write( '{0}\n'.format('# transcriptome with plant transcripts'))
                script_file_id.write( '{0}\n'.format('purge_transcriptome'))
                script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# annotation statistics'))
            script_file_id.write( '{0}\n'.format('calculate_annotation_stats'))
            script_file_id.write( '\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_nucleotide_pipeline_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_nucleotide_pipeline_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_nucleotide_pipeline_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_nucleotide_pipeline_script():
    '''
    Get the script path in the local computer to process a nucleotide pipeline.
    '''

    # assign the script path
    nucleotide_pipeline_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_pipeline_nucleotide_code())

    # return the script path
    return nucleotide_pipeline_script

#-------------------------------------------------------------------------------

def get_nucleotide_pipeline_starter():
    '''
    Get the starter path in the local computer to process a nucleotide pipeline.
    '''

    # assign the starter path
    nucleotide_pipeline_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_pipeline_nucleotide_code())

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
    thread_number = pipeline_option_dict['BLAST parameters']['thread_number']
    e_value = pipeline_option_dict['BLAST parameters']['e_value']
    max_target_seqs = pipeline_option_dict['BLAST parameters']['max_target_seqs']
    max_hsps = pipeline_option_dict['BLAST parameters']['max_hsps']
    qcov_hsp_perc = pipeline_option_dict['BLAST parameters']['qcov_hsp_perc']

    # get the all selected database list
    all_database_list = get_selected_database_list(xlib.get_toa_process_pipeline_aminoacid_code())

    # change code "nt_complete" by "nr_remainder"
    if all_database_list[len(all_database_list) - 1] == 'nr_complete':
        all_database_list[len(all_database_list) - 1] = 'nr_remainder'

    # get the plant database list
    plant_database_list = all_database_list.copy()
    if 'nr_remainder' in plant_database_list:
        plant_database_list.remove('nr_remainder')

    # get the database type dictionary
    database_type_dict = {}
    for i in range(len(all_database_list)):
        if all_database_list[i] in ['gymno_01', 'dicots_04', 'monocots_04']:
            database_type_dict[all_database_list[i]] = 'PLAZA'
        elif all_database_list[i] == 'refseq_plant':
            database_type_dict[all_database_list[i]] = 'REFSEQ'
        elif all_database_list[i] in ['nr_viridiplantae', 'nr_remainder']:
            database_type_dict[all_database_list[i]] = 'NR'

    # set the transcriptome file path
    if OK:
        if assembly_origin == 'NGSCLOUD':
            if assembly_software == xlib.get_soapdenovotrans_code():
                if assembly_type == 'CONTIGS':
                    transcriptome_file = '{0}/{1}-{2}.contig'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id), experiment_id, assembly_dataset_id)
                elif  assembly_type == 'SCAFFOLDS':
                    transcriptome_file = '{0}/{1}-{2}.scafSeq'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id), experiment_id, assembly_dataset_id)
            elif assembly_software == xlib.get_transabyss_code():
                transcriptome_file = '{0}/transabyss-final.fa'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
            elif assembly_software == xlib.get_trinity_code():
                transcriptome_file = '{0}/Trinity.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
            elif assembly_software == xlib.get_cd_hit_est_code():
                transcriptome_file = '{0}/clustered-transcriptome.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
            elif assembly_software == xlib.get_transcript_filter_code():
                transcriptome_file = '{0}/filtered-transcriptome.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
        elif assembly_origin == 'EXTERNAL':
            transcriptome_file = xlib.get_cluster_reference_file(reference_dataset_id, transcriptome_file)

    # get the non annotation file list
    non_annotation_file_list = []
    for database in all_database_list:
        non_annotation_file_list.append('${0}_NON_ANNOTATED_PEPTIDE_FILE'.format(database.upper()))

    # write the script
    try:
        if not os.path.exists(os.path.dirname(get_aminoacid_pipeline_script())):
            os.makedirs(os.path.dirname(get_aminoacid_pipeline_script()))
        with open(get_aminoacid_pipeline_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('# transcriptome file'))
            script_file_id.write( '{0}\n'.format('TRANSCRIPTOME_FILE={0}'.format(transcriptome_file)))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# BLAST parameters'))
            script_file_id.write( '{0}\n'.format('NUM_THREADS={0}'.format(thread_number)))
            script_file_id.write( '{0}\n'.format('E_VALUE={0}'.format(e_value)))
            script_file_id.write( '{0}\n'.format('MAX_TARGET_SEQS={0}'.format(max_target_seqs)))
            script_file_id.write( '{0}\n'.format('MAX_HSPS={0}'.format(max_hsps)))
            script_file_id.write( '{0}\n'.format('QCOV_HSP_PERC={0}'.format(qcov_hsp_perc)))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# output directory'))
            script_file_id.write( '{0}\n'.format('OUTPUT_DIR={0}'.format(current_run_dir)))
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
            script_file_id.write( '{0}\n'.format('export PATH={0}:{1}:$PATH'.format(toa_config_dict['MINICONDA3_BIN_DIR'], toa_config_dict['TOA_DIR'])))
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('mkdir --parents $STATS_DIR'))
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
            script_file_id.write( '{0}\n'.format('    echo "TRANSCRIPTOME FILE: $TRANSCRIPTOME_FILE"'))
            script_file_id.write( '{0}\n'.format('    echo "ALIGNMENT DATASETS: {0}"'.format(','.join(all_database_list))))
            script_file_id.write( '{0}\n'.format('    echo "NUM_THREADS: $NUM_THREADS"'))
            script_file_id.write( '{0}\n'.format('    echo "E_VALUE: $E_VALUE"'))
            script_file_id.write( '{0}\n'.format('    echo "MAX_TARGET_SEQS: $MAX_TARGET_SEQS"'))
            script_file_id.write( '{0}\n'.format('    echo "MAX_HSPS: $MAX_HSPS"'))
            script_file_id.write( '{0}\n'.format('    echo "QCOV_HSP_PERC: $QCOV_HSP_PERC"'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function extract_orfs'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/extract_orfs.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "EXTRACT THE LONG OPEN READING FRAMES"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_transdecoder_anaconda_code()}\n')
            script_file_id.write( '{0}\n'.format('        echo "Extracting ORFs ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            TransDecoder.LongOrfs \\'))
            script_file_id.write( '{0}\n'.format('                -t $TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                -m 100 \\'))
            script_file_id.write( '{0}\n'.format('                --output_dir $TRANSDECODER_OUTPUT_DIR'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error TransDecoder.LongOrfs $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "ORFs are extracted."'))
            script_file_id.write( '{0}\n'.format('        conda deactivate'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function predict_coding_regions'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/predict_coding_regions.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "PREDICT THE LIKELY CODING REGIONS"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_transdecoder_anaconda_code()}\n')
            script_file_id.write( '{0}\n'.format('        echo "Predicting codign regions ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            TransDecoder.Predict \\'))
            script_file_id.write( '{0}\n'.format('                -t $TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --output_dir $TRANSDECODER_OUTPUT_DIR'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error TransDecoder.Predict $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "Coding regions are predicted."'))
            script_file_id.write( '{0}\n'.format('        conda deactivate'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function reidentify_sequences'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/reidentify_sequences.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "RE-INDENTIFY SEQUENCES OF THE PEPTIDE FILE"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write( '{0}\n'.format('        echo "Re-identifing sequences ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            reid-fasta-file.py \\'))
            script_file_id.write( '{0}\n'.format('                --fasta=$PEPTIDE_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --out=$REIDENTIFIED_PEPTIDE_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
            script_file_id.write( '{0}\n'.format('                --trace=N'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error reid-fasta-file.py $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "Sequences are re-identified."'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            for i in range(len(all_database_list)):
                current_code = all_database_list[i]
                previus_code = all_database_list[i - 1] if i > 0 else ''
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function align_peptides_{0}_proteome'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/align_peptides_{0}_proteome.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "ALIGNMENT OF PEPTIDES TO {0} PROTEOME"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write(f'        source activate {xlib.get_blastplus_anaconda_code()}\n')
                script_file_id.write( '{0}\n'.format('        echo "Aligning peptides ..."'))
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    script_file_id.write( '{0}\n'.format('        export BLASTDB=${0}_PROTEOME_DB_DIR'.format(current_code.upper())))
                elif current_code in ['nr_viridiplantae', 'nr_remainder']:
                    script_file_id.write( '{0}\n'.format('        export BLASTDB=$NR_DB_DIR'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            blastp \\'))
                script_file_id.write( '{0}\n'.format('                -num_threads $NUM_THREADS \\'))
                if current_code in ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant']:
                    script_file_id.write( '{0}\n'.format('                -db ${0}_PROTEOME_DB_NAME \\'.format(current_code.upper())))
                elif current_code in ['nr_viridiplantae', 'nr_remainder']:
                    script_file_id.write( '{0}\n'.format('                -db $NR_DB_NAME \\'))
                if i == 0:
                    script_file_id.write( '{0}\n'.format('                -query $REIDENTIFIED_PEPTIDE_FILE \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                -query ${0}_NON_ANNOTATED_PEPTIDE_FILE \\'.format(previus_code.upper())))
                if current_code == 'nr_viridiplantae':
                    script_file_id.write( '{0}\n'.format('                -gilist $PROTEIN_VIRIDIPLANTAE_GI_LIST \\'))
                script_file_id.write( '{0}\n'.format('                -evalue $E_VALUE \\'))
                script_file_id.write( '{0}\n'.format('                -max_target_seqs $MAX_TARGET_SEQS \\'))
                script_file_id.write( '{0}\n'.format('                -max_hsps $MAX_HSPS \\'))
                script_file_id.write( '{0}\n'.format('                -qcov_hsp_perc $QCOV_HSP_PERC \\'))
                script_file_id.write( '{0}\n'.format('                -outfmt 5 \\'))
                if current_code == 'nr_remainder':
                    script_file_id.write( '{0}\n'.format('                -out $REIDENTIFIED_NR_REMAINDER_BLAST_XML'))
                else:
                    script_file_id.write( '{0}\n'.format('                -out ${0}_BLAST_XML'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error blastx $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Alignment is done."'))
                script_file_id.write( '{0}\n'.format('        conda deactivate'))
                if current_code == 'nr_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                    script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_NR_REMAINDER_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                    script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --out=$NR_REMAINDER_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                if len(plant_database_list) == 1 and current_code != 'nr_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Creating plant alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('           cp ${0}_BLAST_XML $REIDENTIFIED_PLANT_BLAST_XML'.format(current_code.upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error cp $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "File is created."'))
                    script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in merged alignment file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                    script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                    script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --out=$PLANT_BLAST_XML \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function load_alignment_{0}_proteome'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/load_alignment_{0}_proteome.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "LOAD OF PEPTIDE ALIGNMENT TO {0} PROTEOME INTO TOA DATABASE"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Loading alignmnet data ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            load-blast-data.py \\'))
                script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('                --dataset={0} \\'.format(current_code)))
                script_file_id.write( '{0}\n'.format('                --format=5 \\'))
                if current_code == 'nr_remainder':
                    script_file_id.write( '{0}\n'.format('                --blast=$REIDENTIFIED_NR_REMAINDER_BLAST_XML \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --blast=${0}_BLAST_XML \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error load-blast-data.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Data are loaded."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function annotate_peptides_{0}'.format(current_code)))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/annotate_peptides_{0}.ok'.format(current_code)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "ANNOTATION OF PEPTIDES WITH {0}"'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Annotating peptides ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            annotate-sequences.py \\'))
                script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
                script_file_id.write( '{0}\n'.format('                --dataset={0} \\'.format(current_code)))
                if i == 0:
                    script_file_id.write( '{0}\n'.format('                --seqs=$REIDENTIFIED_PEPTIDE_FILE \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --seqs=${0}_NON_ANNOTATED_PEPTIDE_FILE \\'.format(previus_code.upper())))
                script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --annotation=${0}_ANNOTATION_FILE \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --nonann=${0}_NON_ANNOTATED_PEPTIDE_FILE \\'.format(current_code.upper())))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error annotate-sequences.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Annotation is done."'))
                if len(plant_database_list) == 1 and current_code != 'nr_remainder':
                    script_file_id.write( '{0}\n'.format('        echo "Creating plant annotation file ..."'))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('            cp ${0}_ANNOTATION_FILE $PLANT_ANNOTATION_FILE'.format(current_code.upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error cp $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "File is created."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            if len(plant_database_list) > 1:
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function merge_plant_alignment_files'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/merge_plant_alignment_files.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "MERGER OF PLANT ALIGNMENT FILES"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Merging alignment files ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            merge-xml-files.py \\'))
                plant_blast_xml_list = []
                for database_code in plant_database_list:
                    plant_blast_xml_list.append('${0}_BLAST_XML'.format(database_code.upper()))
                script_file_id.write( '{0}\n'.format('                --list={0} \\'.format(','.join(plant_blast_xml_list))))
                script_file_id.write( '{0}\n'.format('                --relationships=NONE \\'))
                script_file_id.write( '{0}\n'.format('                --mfile=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-xml-files.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                script_file_id.write( '{0}\n'.format('        echo "Restoring sequence identifications in merged alignment file ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('            restore-ids.py \\'))
                script_file_id.write( '{0}\n'.format('                --in=$REIDENTIFIED_PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --format=XML \\'))
                script_file_id.write( '{0}\n'.format('                --relationships=$RELATIONSHIP_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --out=$PLANT_BLAST_XML \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error restore-ids.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Identifications are restored."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function merge_plant_annotation_files'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/merge_plant_annotation_files.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "MERGER OF PLANT ANNOTATION FILES"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                for database_code in plant_database_list:
                    script_file_id.write( '{0}\n'.format('        {0}_ANNOTATION_FILE_TMP=${0}_ANNOTATION_FILE".tmp"'.format(database_code.upper())))
                    script_file_id.write( '{0}\n'.format('        {0}_ANNOTATION_FILE_SORTED=${0}_ANNOTATION_FILE".sorted"'.format(database_code.upper())))
                tmp_file_list = []
                for i in range(len(plant_database_list) - 2):
                    script_file_id.write( '{0}\n'.format('        MERGED_ANNOTATION_FILE_TMP{0}=$MERGED_ANNOTATION_FILE".tmp{0}"'.format(i + 1)))
                    tmp_file_list.append('$MERGED_ANNOTATION_FILE_TMP{0}'.format(i + 1))
                script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                script_file_id.write( '{0}\n'.format('        echo "Merging `basename ${0}_ANNOTATION_FILE` and `basename ${1}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper(), plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            merge-annotation-files.py \\'))
                script_file_id.write( '{0}\n'.format('                --file1=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('                --type1={0} \\'.format(database_type_dict[plant_database_list[0]])))
                script_file_id.write( '{0}\n'.format('                --file2=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('                --type2={0} \\'.format(database_type_dict[plant_database_list[1]])))
                if len(plant_database_list) > 2:
                    script_file_id.write( '{0}\n'.format('                --mfile=$MERGED_ANNOTATION_FILE_TMP1 \\'))
                    script_file_id.write( '{0}\n'.format('                --header=N \\'))
                else:
                    script_file_id.write( '{0}\n'.format('                --mfile=$PLANT_ANNOTATION_FILE \\'))
                    script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[1].upper())))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                for i in range(2, len(plant_database_list)):
                    script_file_id.write( '{0}\n'.format('        echo "Deleting the header record of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            tail --lines=+2 ${0}_ANNOTATION_FILE > ${0}_ANNOTATION_FILE_TMP'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        echo "Records are deleted."'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error tail $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Sorting data records of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            sort --field-separator=";" --key=1,2 --key=4,4 --key=6,6 < ${0}_ANNOTATION_FILE_TMP > ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error sort $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Records are sorted."'))
                    script_file_id.write( '{0}\n'.format('        echo "Adding annotation of `basename ${0}_ANNOTATION_FILE` to the merged annotation file ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            merge-annotation-files.py \\'))
                    script_file_id.write( '{0}\n'.format('                --file1=$MERGED_ANNOTATION_FILE_TMP{} \\'.format(i - 1)))
                    script_file_id.write( '{0}\n'.format('                --type1=MERGER \\'))
                    script_file_id.write( '{0}\n'.format('                --file2=${0}_ANNOTATION_FILE_SORTED \\'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('                --type2={0} \\'.format(database_type_dict[plant_database_list[i]])))
                    if i < len(plant_database_list) - 1:
                        script_file_id.write( '{0}\n'.format('                --mfile=$MERGED_ANNOTATION_FILE_TMP{0} \\'.format(i)))
                        script_file_id.write( '{0}\n'.format('                --header=N \\'))
                    else:
                        script_file_id.write( '{0}\n'.format('                --mfile=$PLANT_ANNOTATION_FILE \\'))
                        script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                    script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                    script_file_id.write( '{0}\n'.format('                --trace=N'))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error merge-annotation-files.py $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Files are merged."'))
                    script_file_id.write( '{0}\n'.format('        echo "Deleting temporal files of `basename ${0}_ANNOTATION_FILE` ..."'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '{0}\n'.format('            rm ${0}_ANNOTATION_FILE_TMP ${0}_ANNOTATION_FILE_SORTED'.format(plant_database_list[i].upper())))
                    script_file_id.write( '{0}\n'.format('        RC=$?'))
                    script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                    script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        echo "Deleting temporal annotation files ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            rm {0}'.format(' '.join(tmp_file_list))))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "Files are deleted."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function split_merged_plant_annotation_file'))
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/split_merged_plant_annotation_file.ok'))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "SPLIT OF MERGED PLANT ANNOTATION FILE"'))
                script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
                script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
                script_file_id.write( '    else\n')
                script_file_id.write( '{0}\n'.format('        echo "Splitting merged file `basename $PLANT_ANNOTATION_FILE` ..."'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '{0}\n'.format('            split-annotation-file.py \\'))
                script_file_id.write( '{0}\n'.format('                --annotation=$PLANT_ANNOTATION_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --type=MERGER \\'))
                script_file_id.write( '{0}\n'.format('                --header=Y \\'))
                script_file_id.write( '{0}\n'.format('                --rnum=250000 \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
                script_file_id.write( '{0}\n'.format('                --trace=N'))
                script_file_id.write( '{0}\n'.format('        RC=$?'))
                script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error split-annotation-file.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('        echo "File is splitted."'))
                script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function calculate_annotation_stats'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '{0}\n'.format('    STEP_STATUS=$STATUS_DIR/calculate_annotation_stats.ok'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "CALCULATE ANNOTATION STATISTICS"'))
            script_file_id.write( '{0}\n'.format('    if [ -f $STEP_STATUS ]; then'))
            script_file_id.write( '{0}\n'.format('        echo "This step was previously run."'))
            script_file_id.write( '    else\n')
            script_file_id.write( '{0}\n'.format('        echo "Calculating stats ..."'))
            script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('            calculate-annotation-stats.py \\'))
            script_file_id.write( '{0}\n'.format('                --db=$TOA_DB \\'))
            script_file_id.write( '{0}\n'.format('                --transcriptome=$TRANSCRIPTOME_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --peptides=$PEPTIDE_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --dslist={0} \\'.format(','.join(all_database_list))))
            script_file_id.write( '{0}\n'.format('                --nonannlist={0} \\'.format(','.join(non_annotation_file_list))))
            if len(plant_database_list) > 1:
                script_file_id.write( '{0}\n'.format('                --annotation=$PLANT_ANNOTATION_FILE \\'))
                script_file_id.write( '{0}\n'.format('                --type=MERGER \\'))
            else:
                script_file_id.write( '{0}\n'.format('                --annotation=${0}_ANNOTATION_FILE \\'.format(plant_database_list[0].upper())))
                script_file_id.write( '{0}\n'.format('                --type={0} \\'.format(database_type_dict[plant_database_list[0]])))
            script_file_id.write( '{0}\n'.format('                --stats=$ANNOTATION_STATS_FILE \\'))
            script_file_id.write( '{0}\n'.format('                --verbose=N \\'))
            script_file_id.write( '{0}\n'.format('                --trace=N'))
            script_file_id.write( '{0}\n'.format('        RC=$?'))
            script_file_id.write( '{0}\n'.format('        if [ $RC -ne 0 ]; then manage_error calculate-annotation-stats.py $RC; fi'))
            script_file_id.write( '{0}\n'.format('        echo "Stats are calculated."'))
            script_file_id.write( '{0}\n'.format('        touch $STEP_STATUS'))
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
            script_file_id.write( '{0}\n'.format('# extract the long open reading frames'))
            script_file_id.write( '{0}\n'.format('extract_orfs'))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# extract the long open reading frames'))
            script_file_id.write( '{0}\n'.format('predict_coding_regions'))
            script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# re-identify sequences of the peptide file'))
            script_file_id.write( '{0}\n'.format('reidentify_sequences'))
            script_file_id.write( '\n')
            for i in range(len(all_database_list)):
                current_code = all_database_list[i]
                previus_code = all_database_list[i - 1] if i > 0 else ''
                if i == 0:
                    script_file_id.write( '{0}\n'.format('# complete peptide sequences -> {0}'.format(current_code)))
                else:
                    script_file_id.write( '{0}\n'.format('# peptide sequences not annotated with {0} -> {1}'.format(previus_code, current_code)))
                script_file_id.write( '{0}\n'.format('align_peptides_{0}_proteome'.format(current_code)))
                script_file_id.write( '{0}\n'.format('load_alignment_{0}_proteome'.format(current_code)))
                script_file_id.write( '{0}\n'.format('annotate_peptides_{0}'.format(current_code)))
                script_file_id.write( '\n')
            if len(plant_database_list) > 1:
                script_file_id.write( '{0}\n'.format('# merged plant files'))
                script_file_id.write( '{0}\n'.format('merge_plant_alignment_files'))
                script_file_id.write( '{0}\n'.format('merge_plant_annotation_files'))
                script_file_id.write( '{0}\n'.format('# -- split_merged_plant_annotation_file'))
                script_file_id.write( '\n')
            script_file_id.write( '{0}\n'.format('# annotation statistics'))
            script_file_id.write( '{0}\n'.format('calculate_annotation_stats'))
            script_file_id.write( '\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created.'.format(get_aminoacid_pipeline_script()))
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
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_aminoacid_pipeline_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_aminoacid_pipeline_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_script():
    '''
    Get the script path in the local computer to process a amino acid pipeline.
    '''

    # assign the script path
    aminoacid_pipeline_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_pipeline_aminoacid_code())

    # return the script path
    return aminoacid_pipeline_script

#-------------------------------------------------------------------------------

def get_aminoacid_pipeline_starter():
    '''
    Get the starter path in the local computer to process a amino acid pipeline.
    '''

    # assign the starter path
    aminoacid_pipeline_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_toa_process_pipeline_aminoacid_code())

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
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(starter)))
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

def get_nucleotide_annotation_database_code_list():
    '''
    Get the code list of "nucleotide_annotation_database".
    '''

    return ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant', 'nt_viridiplantae', 'nt_complete']

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

    return ['gymno_01', 'dicots_04', 'monocots_04', 'refseq_plant', 'nr_viridiplantae', 'nr_complete']

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
