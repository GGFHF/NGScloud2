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
                error_list.append('{0}\n'.format(error))
                OK = False

    # check the ddRADseqTools directory is created
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ddradseqtools_name())
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

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_ddradseqtools_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseqTools installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the installation script {0} ...\n'.format(get_ddradseqtools_installation_script()))
        (OK, error_list) = build_ddradseqtools_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the ddRADseqTools installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the installation script {0} in the directory {1} of the master ...\n'.format(get_ddradseqtools_installation_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseqtools_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseqTools installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseqTools installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ddradseqtools_installation_starter()))
        (OK, error_list) = build_ddradseqtools_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the ddRADseqTools installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_ddradseqtools_installation_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseqtools_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseqTools installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the ddRADseqTools installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_starter())))
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
    (ddradseqtools_version, ddradseqtools_url) = xconfiguration.get_bioinfo_app_data(xlib.get_ddradseqtools_name())

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
            script_file_id.write( '{0}\n'.format('function remove_ddradseqtools_directory'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing {0} directory ..."'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    if [ -d "{0}" ]; then'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('        rm -rf {0}'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('        echo "The directory is removed."'))
            script_file_id.write( '{0}\n'.format('    else'))
            script_file_id.write( '{0}\n'.format('        echo "The directory is not found."'))
            script_file_id.write( '{0}\n'.format('    fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function download_ddradseqtools_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Downloading the {0} installation file ..."'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            download_script = f'import requests; r = requests.get(\'{ddradseqtools_url}\') ; open(\'{xlib.get_ddradseqtools_name()}.zip\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $PYTHON3_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is downloaded."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function decompress_ddradseqtools_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Decompressing the {0} installation file ..."'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write(f'    unzip -u {xlib.get_ddradseqtools_name()}.zip\n')
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error tar $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is decompressed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function rename_ddradseqtools_directory'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Renaming the {0} directory ..."'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write(f'    mv {xlib.get_ddradseqtools_name()}-master {xlib.get_ddradseqtools_name()}\n')
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
            script_file_id.write( '{0}\n'.format('    chmod u+x {0}/Package/*.py'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo'))
            script_file_id.write( '{0}\n'.format('    echo "The directory is renamed."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function remove_ddradseqtools_installation_file'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Removing the {0} installation file ..."'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(xlib.get_cluster_app_dir())))
            script_file_id.write( '{0}\n'.format('    rm -f {0}.zip'.format(xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The file is removed."'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok('{0} installation'.format(xlib.get_ddradseqtools_name()), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} installation"'.format(xlib.get_project_name(), xlib.get_ddradseqtools_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong('{0} installation'.format(xlib.get_ddradseqtools_name()), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('remove_ddradseqtools_directory'))
            script_file_id.write( '{0}\n'.format('download_ddradseqtools_installation_file'))
            script_file_id.write( '{0}\n'.format('decompress_ddradseqtools_installation_file'))
            script_file_id.write( '{0}\n'.format('rename_ddradseqtools_directory'))
            script_file_id.write( '{0}\n'.format('set_execution_permissions'))
            script_file_id.write( '{0}\n'.format('remove_ddradseqtools_installation_file'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ddradseqtools_installation_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_ddradseqtools_installation_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ddradseqtools_installation_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ddradseqtools_installation_script():
    '''
    Get the ddRADseqTools installation path in the local computer.
    '''

    # assign the ddRADseqTools installation path
    ddradseqtools_installation_script = '{0}/{1}-installation.sh'.format(xlib.get_temp_dir(), xlib.get_ddradseqtools_name())

    # return the ddRADseqTools installation path
    return ddradseqtools_installation_script

#-------------------------------------------------------------------------------

def get_ddradseqtools_installation_starter():
    '''
    Get the ddRADseqTools installation starter path in the local computer.
    '''

    # assign the ddRADseqTools installation starter path
    ddradseqtools_installation_starter = '{0}/{1}-installation-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ddradseqtools_name())

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
            file_id.write( '{0}\n'.format('# This file contains the restriction sites recognition motifs and their cut sites'))
            file_id.write( '{0}\n'.format('# (a cut site is represented by * in the sequence).'))
            file_id.write( '{0}\n'.format('# New enzymes can be included at the end of the file'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# RECORD FORMAT: enzyme_id;restriction_site_seq(5\'->3\')'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('AatII;GACGT*C'))
            file_id.write( '{0}\n'.format('Acc65I;G*GTACC'))
            file_id.write( '{0}\n'.format('AclI;AA*CGTT'))
            file_id.write( '{0}\n'.format('AatII;GACGT*C'))
            file_id.write( '{0}\n'.format('Acc65I;G*GTACC'))
            file_id.write( '{0}\n'.format('AclI;AA*CGTT'))
            file_id.write( '{0}\n'.format('AfeI;AGC*GCT'))
            file_id.write( '{0}\n'.format('AflII;C*TTAAG'))
            file_id.write( '{0}\n'.format('AgeI;A*CCGGT'))
            file_id.write( '{0}\n'.format('ApaI;GGGCC*C'))
            file_id.write( '{0}\n'.format('ApaLI;G*TGCAC'))
            file_id.write( '{0}\n'.format('AscI;GG*CGCGCC'))
            file_id.write( '{0}\n'.format('AseI;AT*TAAT'))
            file_id.write( '{0}\n'.format('AsiSI;GCGAT*CGC'))
            file_id.write( '{0}\n'.format('AvrII;C*CTAGG'))
            file_id.write( '{0}\n'.format('BamHI;G*GATCC'))
            file_id.write( '{0}\n'.format('BclI;T*GATCA'))
            file_id.write( '{0}\n'.format('BglII;A*GATCT'))
            file_id.write( '{0}\n'.format('BmtI;GCTAG*C'))
            file_id.write( '{0}\n'.format('BsiWI;C*GTACG'))
            file_id.write( '{0}\n'.format('BspEI;T*CCGGA'))
            file_id.write( '{0}\n'.format('BspHI;T*CATGA'))
            file_id.write( '{0}\n'.format('BsrGI;T*GTACA'))
            file_id.write( '{0}\n'.format('BssHII;G*CGCGC'))
            file_id.write( '{0}\n'.format('BstBI;TT*CGAA'))
            file_id.write( '{0}\n'.format('BstZ17I;GTA*TAC'))
            file_id.write( '{0}\n'.format('ClaI;AT*CGAT'))
            file_id.write( '{0}\n'.format('Csp6I;G*TAC'))
            file_id.write( '{0}\n'.format('DraI;TTT*AAA'))
            file_id.write( '{0}\n'.format('EagI;C*GGCCG'))
            file_id.write( '{0}\n'.format('EcoRI;G*AATTC'))
            file_id.write( '{0}\n'.format('EcoRV;GAT*ATC'))
            file_id.write( '{0}\n'.format('FseI;GGCCGG*CC'))
            file_id.write( '{0}\n'.format('FspI;TGC*GCA'))
            file_id.write( '{0}\n'.format('HindIII;A*AGCTT'))
            file_id.write( '{0}\n'.format('HpaI;GTT*AAC'))
            file_id.write( '{0}\n'.format('KpnI;GGTAC*C'))
            file_id.write( '{0}\n'.format('MfeI;C*AATTG'))
            file_id.write( '{0}\n'.format('MluI;A*CGCGT'))
            file_id.write( '{0}\n'.format('MscI;TGG*CCA'))
            file_id.write( '{0}\n'.format('MseI;T*TAA'))
            file_id.write( '{0}\n'.format('NaeI;GCC*GGC'))
            file_id.write( '{0}\n'.format('NarI;GG*CGCC'))
            file_id.write( '{0}\n'.format('NcoI;C*CATGG'))
            file_id.write( '{0}\n'.format('NdeI;CA*TATG'))
            file_id.write( '{0}\n'.format('NgoMIV;G*CCGGC'))
            file_id.write( '{0}\n'.format('NheI;G*CTAGC'))
            file_id.write( '{0}\n'.format('NlaIII;CATG*'))
            file_id.write( '{0}\n'.format('NotI;GC*GGCCGC'))
            file_id.write( '{0}\n'.format('NruI;TCG*CGA'))
            file_id.write( '{0}\n'.format('NsiI;ATGCA*T'))
            file_id.write( '{0}\n'.format('PacI;TTAAT*TAA'))
            file_id.write( '{0}\n'.format('PciI;A*CATGT'))
            file_id.write( '{0}\n'.format('PmeI;GTTT*AAAC'))
            file_id.write( '{0}\n'.format('PmlI;CAC*GTG'))
            file_id.write( '{0}\n'.format('PsiI;TTA*TAA'))
            file_id.write( '{0}\n'.format('PspOMI;G*GGCCC'))
            file_id.write( '{0}\n'.format('PstI;CTGCA*G'))
            file_id.write( '{0}\n'.format('PvuI;CGAT*CG'))
            file_id.write( '{0}\n'.format('PvuII;CAG*CTG'))
            file_id.write( '{0}\n'.format('SacI;GAGCT*C'))
            file_id.write( '{0}\n'.format('SacII;CCGC*GG'))
            file_id.write( '{0}\n'.format('SalI;G*TCGAC'))
            file_id.write( '{0}\n'.format('SbfI;CCTGCA*GG'))
            file_id.write( '{0}\n'.format('ScaI;AGT*ACT'))
            file_id.write( '{0}\n'.format('SfoI;GGC*GCC'))
            file_id.write( '{0}\n'.format('SmaI;CCC*GGG'))
            file_id.write( '{0}\n'.format('SnaBI;TAC*GTA'))
            file_id.write( '{0}\n'.format('SpeI;A*CTAGT'))
            file_id.write( '{0}\n'.format('SphI;GCATG*C'))
            file_id.write( '{0}\n'.format('SspI;AAT*ATT'))
            file_id.write( '{0}\n'.format('StuI;AGG*CCT'))
            file_id.write( '{0}\n'.format('SwaI;ATTT*AAAT'))
            file_id.write( '{0}\n'.format('XbaI;T*CTAGA'))
            file_id.write( '{0}\n'.format('XhoI;C*TCGAG'))
            file_id.write( '{0}\n'.format('XmaI;C*CCGGG'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('AccI;GT*MKAC'))
            file_id.write( '{0}\n'.format('AccB1I;G*GYRCC'))
            file_id.write( '{0}\n'.format('AccB2I;RGCGC*Y'))
            file_id.write( '{0}\n'.format('AceI;G*CWGC'))
            file_id.write( '{0}\n'.format('AdeI;CACNNN*GTG'))
            file_id.write( '{0}\n'.format('MslI;CAYNN*NNRTG'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_restriction_site_file()))
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

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # set the pattern of the record of restriction site file (record format: enzyme_id;restriction_site_seq)
    pattern = r'^(.*);(.*)$'

    # open the file of restriction sites
    try:
        restriction_site_file_id = open(get_restriction_site_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_restriction_site_file()))
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
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the restriction site sequence is correct
                if not xlib.is_valid_sequence(restriction_site_seq, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: The sequence is invalid in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the enzyme identification is formed by alphanumeric characters
                if not xlib.is_name_valid(enzyme_id):
                    error_list.append('*** ERROR: The enzyme id has characters non-alphanumeric in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

            # read the next record
            record = restriction_site_file_id.readline()

    # close the file of restriction sites
    restriction_site_file_id.close()

    # warn that the file of restriction sites is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_restriction_site_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_restriction_site_file():
    '''
    Get the restriction site file path.
    '''

    # assign the restriction site file path
    restriction_site_file = '{0}/restrictionsites.txt'.format(xlib.get_config_dir())

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
                    enzyme_id_seq = '{0} ({1})'.format(enzyme_id, restriction_site_seq)
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
    except Exception as e:
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
            file_id.write( '{0}\n'.format('# This file contains end sequences integrated by indexes, degenerate nucleotides to indentify'))
            file_id.write( '{0}\n'.format('# the PCR duplicates (DBR), adapter and primer. A read has two ends: one where the adapter 1'))
            file_id.write( '{0}\n'.format('#  is and another where adapter 2 is.'))
            file_id.write( '{0}\n'.format('# '))
            file_id.write( '{0}\n'.format('# IND1_IND2_DBR technique: The sequence of the end corresponding to the adapter 1 includes'))
            file_id.write( '{0}\n'.format('# a index1 sequence (111..., each digit 1 represents a nucleotide of the index1), the sequence'))
            file_id.write( '{0}\n'.format('# of the end corresponing to the adapter 2 include a index2 sequence (222..., each digit 2'))
            file_id.write( '{0}\n'.format('# represents a nucleotide of the index2). A DBR sequence (333..., each digit 3 represents a'))
            file_id.write( '{0}\n'.format('# nucleotide of the DBR) has to be included in the end sequence of adapter 1 or adapter2.'))
            file_id.write( '{0}\n'.format('# '))
            file_id.write( '{0}\n'.format('# IND1_IND2 technique: The sequence of the end corresponding to the adapter 1 includes'))
            file_id.write( '{0}\n'.format('# a index1 sequence (111...), the sequence of the end corresponing to the adapter 2 include a'))
            file_id.write( '{0}\n'.format('# index2 sequence (222...). The DBR sequence (333...) is not considered.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# IND1_DBR technique: The sequence of the end corresponding to the adapter 1 includes a index1'))
            file_id.write( '{0}\n'.format('# sequence (111...) and a DBR sequence (333...).  The index2 sequence (222...) is not'))
            file_id.write( '{0}\n'.format('# considered.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# IND1 technique: The sequence of the end corresponding to the adapter 1 includes a index1'))
            file_id.write( '{0}\n'.format('# sequence (111...). The index2 sequence (222...) and the DBR sequence (333...) are not'))
            file_id.write( '{0}\n'.format('# considered.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The length of index1 (111...), index2 (222...) and DBR (333...) of ends used by a program'))
            file_id.write( '{0}\n'.format('# have to be equal to the value of the options index1len, index2len and dbrlen received by'))
            file_id.write( '{0}\n'.format('# that program.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# RECORD FORMAT: end_id;end_seq(5\'->3\')'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end01;AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT111111C'))
            file_id.write( '{0}\n'.format('end02;CAAGCAGAAGACGGCATACGAGAT3333222222GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end03;CAAGCAGAAGACGGCATACGAGAT111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end04;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end05;CAAGCAGAAGACGGCATACGAGAT3333111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end06;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end07;CAAGCAGAAGACGGCATACGAGAT3333111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end08;AATGATACGGCGACCACCGAGATCTACAC222222ACACTCTTTCCCTACACGACGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end09;CAAGCAGAAGACGGCATACGAGAT111111GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'))
            file_id.write( '{0}\n'.format('end10;AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATC'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end51;111111'))
            file_id.write( '{0}\n'.format('end52;'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end61;3333111111'))
            file_id.write( '{0}\n'.format('end62;'))
            file_id.write( '{0}\n'.format('end63;33331111111'))
            file_id.write( '{0}\n'.format('end64;'))
            file_id.write( '{0}\n'.format('end65;333331111111'))
            file_id.write( '{0}\n'.format('end66;'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end71;111111'))
            file_id.write( '{0}\n'.format('end72;222222'))
            file_id.write( '{0}\n'.format('end73;11111'))
            file_id.write( '{0}\n'.format('end74;22222'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end81;111111'))
            file_id.write( '{0}\n'.format('end82;3333222222'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('end91;3333111111'))
            file_id.write( '{0}\n'.format('end92;222222'))
            file_id.write( '{0}\n'.format('end93;33331111111'))
            file_id.write( '{0}\n'.format('end94;2222222'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_end_file()))
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

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # set the pattern of the record of end file (record format: end_id;end_seq)
    pattern = r'^(.*);(.*)$'

    # open the file of ends
    try:
        end_file_id = open(get_end_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_end_file()))
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
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the end sequence is correct
                if not xlib.is_valid_sequence(end_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=['1', '2', '3'], cut_tag_check=False):
                    error_list.append('*** ERROR: The sequence is invalid in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the end identification is formed by alphanumeric characters
                if not xlib.is_name_valid(end_id):
                    error_list.append('*** ERROR: The end id has characters non-alphanumeric in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

            # read the next record
            record = end_file_id.readline()

    # close the file of ends
    end_file_id.close()

    # warn that the file of ends is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_end_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_end_file():
    '''
    Get the end file path.
    '''

    # assign the end file path
    end_file = '{0}/ends.txt'.format(xlib.get_config_dir())

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
            file_id.write( '{0}\n'.format('# This file contains the sequences of each invidual: index1 sequences that are attached to the'))
            file_id.write( '{0}\n'.format('# 5\' end of the Watson strand and, optionally, the index2 sequences that are attached to the'))
            file_id.write( '{0}\n'.format('# 5\' end of the Crick strand.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# RECORD FORMAT: individual_id;replicated_individual_id or NONE;population_id;index1_seq(5\'->3\' in Watson strand);[index2_seq(5\'->3\' in Crick strand)]'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# These data corresponding to an example with 8 index1 (6bp) and 6 index2 (6bp) by index1 -> up 48 individuals.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0101;NONE;pop01;GGTCTT;ATCACG'))
            file_id.write( '{0}\n'.format('ind0102;NONE;pop01;GGTCTT;CGATGT'))
            file_id.write( '{0}\n'.format('ind0103;NONE;pop01;GGTCTT;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0104;NONE;pop01;GGTCTT;TGACCA'))
            file_id.write( '{0}\n'.format('ind0105;NONE;pop01;GGTCTT;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0103r;ind0103;pop01;GGTCTT;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0201;NONE;pop01;CTGGTT;ATCACG'))
            file_id.write( '{0}\n'.format('ind0202;NONE;pop01;CTGGTT;CGATGT'))
            file_id.write( '{0}\n'.format('ind0203;NONE;pop01;CTGGTT;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0204;NONE;pop01;CTGGTT;TGACCA'))
            file_id.write( '{0}\n'.format('ind0205;NONE;pop01;CTGGTT;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0204r;ind0204;pop01;CTGGTT;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0301;NONE;pop01;AAGATA;ATCACG'))
            file_id.write( '{0}\n'.format('ind0302;NONE;pop01;AAGATA;CGATGT'))
            file_id.write( '{0}\n'.format('ind0303;NONE;pop01;AAGATA;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0304;NONE;pop01;AAGATA;TGACCA'))
            file_id.write( '{0}\n'.format('ind0305;NONE;pop01;AAGATA;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0306;NONE;pop01;AAGATA;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0401;NONE;pop01;ACTTCC;ATCACG'))
            file_id.write( '{0}\n'.format('ind0402;NONE;pop01;ACTTCC;CGATGT'))
            file_id.write( '{0}\n'.format('ind0403;NONE;pop01;ACTTCC;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0404;NONE;pop01;ACTTCC;TGACCA'))
            file_id.write( '{0}\n'.format('ind0405;NONE;pop01;ACTTCC;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0406;NONE;pop01;ACTTCC;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0501;NONE;pop01;TTACGG;ATCACG'))
            file_id.write( '{0}\n'.format('ind0502;NONE;pop01;TTACGG;CGATGT'))
            file_id.write( '{0}\n'.format('ind0503;NONE;pop01;TTACGG;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0504;NONE;pop01;TTACGG;TGACCA'))
            file_id.write( '{0}\n'.format('ind0505;NONE;pop01;TTACGG;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0506;NONE;pop01;TTACGG;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0601;NONE;pop01;AACGAA;ATCACG'))
            file_id.write( '{0}\n'.format('ind0602;NONE;pop01;AACGAA;CGATGT'))
            file_id.write( '{0}\n'.format('ind0603;NONE;pop01;AACGAA;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0604;NONE;pop01;AACGAA;TGACCA'))
            file_id.write( '{0}\n'.format('ind0605;NONE;pop01;AACGAA;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0606;NONE;pop01;AACGAA;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0701;NONE;pop01;ATTCAT;ATCACG'))
            file_id.write( '{0}\n'.format('ind0702;NONE;pop01;ATTCAT;CGATGT'))
            file_id.write( '{0}\n'.format('ind0703;NONE;pop01;ATTCAT;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0704;NONE;pop01;ATTCAT;TGACCA'))
            file_id.write( '{0}\n'.format('ind0705;NONE;pop01;ATTCAT;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0706;NONE;pop01;ATTCAT;GCCAAT'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('ind0801;NONE;pop01;CCGACC;ATCACG'))
            file_id.write( '{0}\n'.format('ind0802;NONE;pop01;CCGACC;CGATGT'))
            file_id.write( '{0}\n'.format('ind0803;NONE;pop01;CCGACC;TTAGGC'))
            file_id.write( '{0}\n'.format('ind0804;NONE;pop01;CCGACC;TGACCA'))
            file_id.write( '{0}\n'.format('ind0805;NONE;pop01;CCGACC;ACAGTG'))
            file_id.write( '{0}\n'.format('ind0806;NONE;pop01;CCGACC;GCCAAT'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_individual_file()))
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

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # set the pattern of the record of individual file (record format: individual_id;replicated_individual_id;population_id;index1_seq;[index2_seq])
    pattern = r'^(.+);(.+);(.+);(.+);(.*)$'

    # open the file of individuals
    try:
        individual_file_id = open(get_individual_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_individual_file()))
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
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the index1 sequence is correct
                if not xlib.is_valid_sequence(index1_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                    error_list.append('*** ERROR: The index 1 sequence is invalid in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the index2 sequence is correct
                if index2_seq != '':
                    if not xlib.is_valid_sequence(index2_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                        error_list.append('*** ERROR: The index 2 sequence is invalid in the record "{0}".'.format(record.replace("\n", "")))
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
                    error_list.append('*** ERROR: The individual id has characters non-alphanumeric in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the replicated individual identification is formed by alphanumeric characters
                if not xlib.is_name_valid(replicated_individual_id):
                    error_list.append('*** ERROR: The replicated individual id has characters non-alphanumeric in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # check that the population identification is formed by alphanumeric characters
                if not xlib.is_name_valid(population_id):
                    error_list.append('*** ERROR: The population id has characters non-alphanumeric in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

            # read the next record
            record = individual_file_id.readline()

    # close the file of individuals
    individual_file_id.close()

    # warn that the file of individuals is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_individual_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_individual_file():
    '''
    Get the individual file path.
    '''

    # assign the individual file path
    individual_file = '{0}/individuals.txt'.format(xlib.get_config_dir())

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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference file has to be located in the cluster directory {0}/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The reference_dataset_id and reference_file_name are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of rsitesearch (ddRADseqTools package) and their meaning in https://github.com/GGFHF/.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the rsitesearch parameters'))
            file_id.write( '{0}\n'.format('[rsitesearch parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('enzyme1 = {0}'.format(enzyme1), '# id of 1st restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format('enzyme2 = {0}'.format(enzyme2), '# id of 2nd restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format('minfragsize = 101', '# lower loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format('maxfragsize = 300', '# upper loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format('fragstinterval = 25', '# interval length of fragment size'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_rsitesearch_config_file()))
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

    # get the rsitesearch option dictionary
    rsitesearch_option_dict = xlib.get_option_dict(get_rsitesearch_config_file())

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the file of restriction sites
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the file {0} ...\n'.format(get_restriction_site_file()))
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
        log.write('Checking the {0} config file ...\n'.format(xlib.get_rsitesearch_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the ddRADseqTools is installed
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ddradseqtools_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ddradseqtools_name()))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('simulation', xlib.get_rsitesearch_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the file of restriction sites to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_restriction_site_file(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_restriction_site_file()))
        (OK, error_list) = xssh.put_file(sftp_client, get_restriction_site_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the rsitesearch process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_rsitesearch_process_script()))
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
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_rsitesearch_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_rsitesearch_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_rsitesearch_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the rsitesearch process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsitesearch_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_rsitesearch_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the rsitesearch process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_rsitesearch_process_starter()))
        (OK, error_list) = build_rsitesearch_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the rsitesearch process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_rsitesearch_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_rsitesearch_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_rsitesearch_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the rsitesearch process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsitesearch_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_rsitesearch_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the rsitesearch process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsitesearch_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
                error_list.append('*** ERROR: The value maxfragsize value ({0}) is less than the minfragsize value ({1}).'.format(maxfragsize, minfragsize))
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
        error_list.append('The {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_rsitesearch_name()))

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
    rsfile = '{0}/{1}'.format(current_run_dir, os.path.basename(get_restriction_site_file()))

    fragsfile = '{0}/genome-fragments.txt'.format(current_run_dir)
    fragstfile = '{0}/genome-fragment-stats.txt'.format(current_run_dir)

    # write the rsitesearch process script
    if OK:
        try:
            if not os.path.exists(os.path.dirname(get_rsitesearch_process_script())):
                os.makedirs(os.path.dirname(get_rsitesearch_process_script()))
            with open(get_rsitesearch_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
                script_file_id.write( '#!/bin/bash\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( '{0}\n'.format('PYTHON3_PATH={0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('DDRADSEQTOOLS_PATH={0}/{1}/Package'.format(xlib.get_cluster_app_dir(), xlib.get_ddradseqtools_name())))
                script_file_id.write( '{0}\n'.format('export PATH=$PYTHON3_PATH:$DDRADSEQTOOLS_PATH:$PATH'))
                script_file_id.write( '{0}\n'.format('export MPLBACKEND="agg"'))
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
                script_file_id.write( '{0}\n'.format('function run_rsitesearch_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Analyzing enzymes of a ddRADseq experiment ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        rsitesearch.py \\'))
                script_file_id.write( '{0}\n'.format('            --genfile={0} \\'.format(genfile)))
                script_file_id.write( '{0}\n'.format('            --fragsfile={0} \\'.format(fragsfile)))
                script_file_id.write( '{0}\n'.format('            --rsfile={0} \\'.format(rsfile)))
                script_file_id.write( '{0}\n'.format('            --enzyme1={0} \\'.format(enzyme1)))
                script_file_id.write( '{0}\n'.format('            --enzyme2={0} \\'.format(enzyme2)))
                script_file_id.write( '{0}\n'.format('            --minfragsize={0} \\'.format(minfragsize)))
                script_file_id.write( '{0}\n'.format('            --maxfragsize={0} \\'.format(maxfragsize)))
                script_file_id.write( '{0}\n'.format('            --fragstfile={0} \\'.format(fragstfile)))
                script_file_id.write( '{0}\n'.format('            --fragstinterval={0} \\'.format(fragstinterval)))
                script_file_id.write( '{0}\n'.format('            --plot=YES \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=NO \\'))
                script_file_id.write( '{0}\n'.format('            --trace=NO'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rsitesearch.py $RC; fi'))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_rsitesearch_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_rsitesearch_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_rsitesearch_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_rsitesearch_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('run_rsitesearch_process'))
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be created.'.format(get_rsitesearch_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_rsitesearch_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_rsitesearch_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_rsitesearch_config_file():
    '''
    Get the rsitesearch config file path.
    '''

    # assign the rsitesearch config file path
    rsitesearch_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_rsitesearch_code())

    # return the rsitesearch config file path
    return rsitesearch_config_file

#-------------------------------------------------------------------------------

def get_rsitesearch_process_script():
    '''
    Get the rsitesearch process script path in the local computer.
    '''

    # assign the rsitesearch script path
    rsitesearch_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_rsitesearch_code())

    # return the rsitesearch script path
    return rsitesearch_process_script

#-------------------------------------------------------------------------------

def get_rsitesearch_process_starter():
    '''
    Get the rsitesearch process starter path in the local computer.
    '''

    # assign the rsitesearch process starter path
    rsitesearch_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_rsitesearch_code())

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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference file has to be located in the cluster directory {0}/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The reference_dataset_id and reference_file_name are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of ddRADseqTools programs and their meaning in https://github.com/GGFHF/.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format('rsfile = ./config/restrictionsites.txt', '# local path of the restriction sites file'))
            file_id.write( '{0:<50} {1}\n'.format('endsfile = ./config/ends.txt', '# local path oh the end sequences file'))
            file_id.write( '{0:<50} {1}\n'.format('individualsfile = ./config/individuals.txt', '# local path oh the end sequences file'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the ddRADseq simulation parameters'))
            file_id.write( '{0}\n'.format('[ddRADseq simulation parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('enzyme1 = {0}'.format(enzyme1), '# id of 1st restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format('enzyme2 = {0}'.format(enzyme2), '# id of 2nd restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format('minfragsize = 101', '# lower loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format('maxfragsize = 300', '# upper loci fragment size'))
            file_id.write( '{0:<50} {1}\n'.format('fragstinterval = 25', '# interval length of fragment size'))
            file_id.write( '{0:<50} {1}\n'.format('technique = IND1_IND2_DBR', '# technique: {0}'.format(get_technique_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('readtype = PE', '# read type: {0}'.format(get_read_type_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('index1len = 6', '# index sequence length in the adapter 1'))
            file_id.write( '{0:<50} {1}\n'.format('index2len = 6', '# index sequence length in the adapter 2 (it has to be 0 when technique is IND1)'))
            file_id.write( '{0:<50} {1}\n'.format('dbrlen = 4', '# DBR sequence length (it has to be 0 when technique is IND1 or IND1_IND2)'))
            file_id.write( '{0:<50} {1}\n'.format('wend = end91', '# code used in endsfile corresponding to the end where the adapter 1 is'))
            file_id.write( '{0:<50} {1}\n'.format('cend = end92', '# code used in endsfile corresponding to the end where the adapter 2 is'))
            file_id.write( '{0:<50} {1}\n'.format('locinum = 3000', '# loci number to sample'))
            file_id.write( '{0:<50} {1}\n'.format('readsnum = 300000', '# reads number'))
            file_id.write( '{0:<50} {1}\n'.format('minreadvar = 0.8', '# lower variation on reads number per locus (0.5 <= minreadvar <= 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('maxreadvar = 1.2', '# upper variation on reads number per locus (1.0 <= maxreadvar <= 1.5)'))
            file_id.write( '{0:<50} {1}\n'.format('insertlen = 100', '# read length, i. e. genome sequence length inserted in reads'))
            file_id.write( '{0:<50} {1}\n'.format('mutprob = 0.2', '# mutation probability (0.0 <= mutprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('locusmaxmut = 1', '# maximum mutations number by locus (1 <= locusmaxmut <= 5)'))
            file_id.write( '{0:<50} {1}\n'.format('indelprob = 0.1', '# insertion/deletion probability (0.0 <= indelprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('maxindelsize = 10', '# upper insertion/deletion size (1 <= maxindelsize < 30)'))
            file_id.write( '{0:<50} {1}\n'.format('dropout = 0.0', '# mutation probability in the enzyme recognition sites (0.0 <= dropout < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('pcrdupprob = 0.2', '# probability of loci bearing PCR duplicates (0.0 <= pcrdupprob < 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('pcrdistribution = MULTINOMIAL', '# distribution type to calculate the PCR duplicates number: {0}'.format(get_pcrdistribution_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('multiparam = 0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015', '# probability values to multinomial distribution with format prob1,prob2,...,probn (they have to sum 1.0)'))
            file_id.write( '{0:<50} {1}\n'.format('poissonparam = 1.0', '# lambda value of the Poisson distribution'))
            file_id.write( '{0:<50} {1}\n'.format('gcfactor = 0.2', '# weight factor of GC ratio in a locus with PCR duplicates (0.0 <= gcfactor < 1.0)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_ddradseq_simulation_config_file()))
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

    # get the ddRADseq simulation option dictionary
    ddradseq_simulation_option_dict = xlib.get_option_dict(get_ddradseq_simulation_config_file())

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the ddRADseq simulation config file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Checking the {0} config file ...\n'.format(xlib.get_ddradseq_simulation_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the ddRADseq simulation is installed
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_ddradseqtools_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ddradseqtools_name()))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('simulation', xlib.get_ddradseq_simulation_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the file of restriction sites to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_restriction_site_file(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_restriction_site_file()))
        (OK, error_list) = xssh.put_file(sftp_client, get_restriction_site_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the file of ends to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_end_file(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_end_file()))
        (OK, error_list) = xssh.put_file(sftp_client, get_end_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the file of individuals to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_individual_file(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_individual_file()))
        (OK, error_list) = xssh.put_file(sftp_client, get_individual_file(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # build the ddRADseq simulation process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_ddradseq_simulation_process_script()))
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
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_ddradseq_simulation_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseq_simulation_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseq simulation process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ddRADseq simulation process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ddradseq_simulation_process_starter()))
        (OK, error_list) = build_ddradseq_simulation_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the ddRADseq simulation process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_ddradseq_simulation_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ddradseq_simulation_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ddRADseq simulation process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the ddRADseq simulation process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
                error_list.append('*** ERROR: The value maxfragsize value ({0}) is less than the minfragsize value ({1}).'.format(maxfragsize, minfragsize))
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
                error_list.append('*** ERROR: the key "technique" has to be {0}.'.format(get_technique_code_list_text()))
                OK = False
            else:
                is_ok_technique = True

            # check section "ddRADseq simulation parameters" - key "format"
            format = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "ddRADseq simulation parameters" - key "readtype"
            readtype = ddradseq_simulation_option_dict.get('ddRADseq simulation parameters', {}).get('readtype', not_found)
            if readtype == not_found:
                error_list.append('*** ERROR: the key "readtype" is not found in the section "ddRADseq simulation parameters".')
                OK = False
            elif not xlib.check_code(readtype, get_read_type_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "readtype" has to be {0}.'.format(get_read_type_code_list_text()))
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
                error_list.append('*** ERROR: the key "pcrdistribution" has to be {0}.'.format(get_pcrdistribution_code_list_text()))
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
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_ddradseq_simulation_name()))

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
    rsfile = '{0}/{1}'.format(current_run_dir, os.path.basename(get_restriction_site_file()))
    endsfile = '{0}/{1}'.format(current_run_dir, os.path.basename(get_end_file()))
    individualsfile = '{0}/{1}'.format(current_run_dir, os.path.basename(get_individual_file()))
    fragsfile = '{0}/genome-fragments.txt'.format(current_run_dir)
    fragstfile = '{0}/genome-fragment-stats.txt'.format(current_run_dir)
    readsfile = '{0}/reads'.format(current_run_dir)
    readsfile1 = '{0}/reads-1'.format(current_run_dir)
    readsfile2 = '{0}/reads-2'.format(current_run_dir)
    clearfile = '{0}/reads-cleared'.format(current_run_dir)
    clearfile1 = '{0}/reads-cleared-1'.format(current_run_dir)
    clearfile2 = '{0}/reads-cleared-2'.format(current_run_dir)
    dupstfile = '{0}/pcrduplicates-stats'.format(current_run_dir)
    read_file_list = '{0}/reads-files.txt'.format(current_run_dir)
    demultiplexed_file = '{0}/demultiplexed-ind'.format(current_run_dir)

    # write the ddRADseq simulation process script
    if OK:
        try:
            if not os.path.exists(os.path.dirname(get_ddradseq_simulation_process_script())):
                os.makedirs(os.path.dirname(get_ddradseq_simulation_process_script()))
            with open(get_ddradseq_simulation_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
                script_file_id.write( '#!/bin/bash\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'SEP="#########################################"\n')
                script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
                script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
                script_file_id.write( '{0}\n'.format('PYTHON3_PATH={0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
                script_file_id.write( '{0}\n'.format('DDRADSEQTOOLS_PATH={0}/{1}/Package'.format(xlib.get_cluster_app_dir(), xlib.get_ddradseqtools_name())))
                script_file_id.write( '{0}\n'.format('export PATH=$PYTHON3_PATH:$DDRADSEQTOOLS_PATH:$PATH'))
                script_file_id.write( '{0}\n'.format('export MPLBACKEND="agg"'))
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
                script_file_id.write( '{0}\n'.format('function run_rsitesearch_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating fragments ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        rsitesearch.py \\'))
                script_file_id.write( '{0}\n'.format('            --genfile={0} \\'.format(genfile)))
                script_file_id.write( '{0}\n'.format('            --fragsfile={0} \\'.format(fragsfile)))
                script_file_id.write( '{0}\n'.format('            --rsfile={0} \\'.format(rsfile)))
                script_file_id.write( '{0}\n'.format('            --enzyme1={0} \\'.format(enzyme1)))
                script_file_id.write( '{0}\n'.format('            --enzyme2={0} \\'.format(enzyme2)))
                script_file_id.write( '{0}\n'.format('            --minfragsize={0} \\'.format(minfragsize)))
                script_file_id.write( '{0}\n'.format('            --maxfragsize={0} \\'.format(maxfragsize)))
                script_file_id.write( '{0}\n'.format('            --fragstfile={0} \\'.format(fragstfile)))
                script_file_id.write( '{0}\n'.format('            --fragstinterval={0} \\'.format(fragstinterval)))
                script_file_id.write( '{0}\n'.format('            --plot=YES \\'))
                script_file_id.write( '{0}\n'.format('            --verbose=NO \\'))
                script_file_id.write( '{0}\n'.format('            --trace=NO'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rsitesearch.py $RC; fi'))
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function run_simddradseq_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Generating reads ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        simddradseq.py \\'))
                script_file_id.write( '{0}\n'.format('            --fragsfile={0} \\'.format(fragsfile)))
                script_file_id.write( '{0}\n'.format('            --technique={0} \\'.format(technique)))
                script_file_id.write( '{0}\n'.format('            --format={0} \\'.format(format)))
                script_file_id.write( '{0}\n'.format('            --readsfile={0} \\'.format(readsfile)))
                script_file_id.write( '{0}\n'.format('            --readtype={0} \\'.format(readtype)))
                script_file_id.write( '{0}\n'.format('            --rsfile={0} \\'.format(rsfile)))
                script_file_id.write( '{0}\n'.format('            --enzyme1={0} \\'.format(enzyme1)))
                script_file_id.write( '{0}\n'.format('            --enzyme2={0} \\'.format(enzyme2)))
                script_file_id.write( '{0}\n'.format('            --endsfile={0} \\'.format(endsfile)))
                script_file_id.write( '{0}\n'.format('            --index1len={0} \\'.format(index1len)))
                script_file_id.write( '{0}\n'.format('            --index2len={0} \\'.format(index2len)))
                script_file_id.write( '{0}\n'.format('            --dbrlen={0} \\'.format(dbrlen)))
                script_file_id.write( '{0}\n'.format('            --wend={0} \\'.format(wend)))
                script_file_id.write( '{0}\n'.format('            --cend={0} \\'.format(cend)))
                script_file_id.write( '{0}\n'.format('            --individualsfile={0} \\'.format(individualsfile)))
                script_file_id.write( '{0}\n'.format('            --locinum={0} \\'.format(locinum)))
                script_file_id.write( '{0}\n'.format('            --readsnum={0} \\'.format(readsnum)))
                script_file_id.write( '{0}\n'.format('            --minreadvar={0} \\'.format(minreadvar)))
                script_file_id.write( '{0}\n'.format('            --maxreadvar={0} \\'.format(maxreadvar)))
                script_file_id.write( '{0}\n'.format('            --insertlen={0} \\'.format(insertlen)))
                script_file_id.write( '{0}\n'.format('            --mutprob={0} \\'.format(mutprob)))
                script_file_id.write( '{0}\n'.format('            --locusmaxmut={0} \\'.format(locusmaxmut)))
                script_file_id.write( '{0}\n'.format('            --indelprob={0} \\'.format(indelprob)))
                script_file_id.write( '{0}\n'.format('            --maxindelsize={0} \\'.format(maxindelsize)))
                script_file_id.write( '{0}\n'.format('            --dropout={0} \\'.format(dropout)))
                script_file_id.write( '{0}\n'.format('            --pcrdupprob={0} \\'.format(pcrdupprob)))
                script_file_id.write( '{0}\n'.format('            --pcrdistribution={0} \\'.format(pcrdistribution)))
                script_file_id.write( '{0}\n'.format('            --multiparam={0} \\'.format(multiparam)))
                script_file_id.write( '{0}\n'.format('            --poissonparam={0} \\'.format(poissonparam)))
                script_file_id.write( '{0}\n'.format('            --gcfactor={0} \\'.format(gcfactor)))
                script_file_id.write( '{0}\n'.format('            --verbose=NO \\'))
                script_file_id.write( '{0}\n'.format('            --trace=NO'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error simddradseq.py $RC; fi'))
                script_file_id.write( '}\n')
                if int(dbrlen) > 0:
                    script_file_id.write( '#-------------------------------------------------------------------------------\n')
                    script_file_id.write( '{0}\n'.format('function run_pcrdupremoval_process'))
                    script_file_id.write( '{\n')
                    script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                    script_file_id.write( '    echo "$SEP"\n')
                    script_file_id.write( '{0}\n'.format('    echo "Removing PCR duplicates ..."'))
                    script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                    script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                    script_file_id.write( '{0}\n'.format('        pcrdupremoval.py \\'))
                    script_file_id.write( '{0}\n'.format('            --format={0} \\'.format(format)))
                    script_file_id.write( '{0}\n'.format('            --readtype={0} \\'.format(readtype)))
                    if format == 'FASTA':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(readsfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(readsfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fasta \\'.format(readsfile2)))
                    elif format == 'FASTQ':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(readsfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(readsfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fastq \\'.format(readsfile2)))
                    script_file_id.write( '{0}\n'.format('            --clearfile={0} \\'.format(clearfile)))
                    script_file_id.write( '{0}\n'.format('            --dupstfile={0} \\'.format(dupstfile)))
                    script_file_id.write( '{0}\n'.format('            --plot=YES \\'))
                    script_file_id.write( '{0}\n'.format('            --verbose=NO \\'))
                    script_file_id.write( '{0}\n'.format('            --trace=NO'))
                    script_file_id.write( '{0}\n'.format('    RC=$?'))
                    script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error pcrdupremoval.py $RC; fi'))
                    script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function run_indsdemultiplexing_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Demultiplexing individuals ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        indsdemultiplexing.py \\'))
                script_file_id.write( '{0}\n'.format('            --technique={0} \\'.format(technique)))
                script_file_id.write( '{0}\n'.format('            --format={0} \\'.format(format)))
                script_file_id.write( '{0}\n'.format('            --readtype={0} \\'.format(readtype)))
                script_file_id.write( '{0}\n'.format('            --endsfile={0} \\'.format(endsfile)))
                script_file_id.write( '{0}\n'.format('            --index1len={0} \\'.format(index1len)))
                script_file_id.write( '{0}\n'.format('            --index2len={0} \\'.format(index2len)))
                script_file_id.write( '{0}\n'.format('            --dbrlen={0} \\'.format(dbrlen)))
                script_file_id.write( '{0}\n'.format('            --wend={0} \\'.format(wend)))
                script_file_id.write( '{0}\n'.format('            --cend={0} \\'.format(cend)))
                script_file_id.write( '{0}\n'.format('            --individualsfile={0} \\'.format(individualsfile)))
                if int(dbrlen) > 0:
                    if format == 'FASTA':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(clearfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(clearfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fasta \\'.format(clearfile2)))
                    elif format == 'FASTQ':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(clearfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(clearfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fastq \\'.format(clearfile2)))
                else:
                    if format == 'FASTA':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(readsfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fasta \\'.format(readsfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fasta \\'.format(readsfile2)))
                    elif format == 'FASTQ':
                        if readtype == 'SE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(readsfile)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2=NONE \\'))
                        elif readtype == 'PE':
                            script_file_id.write( '{0}\n'.format('            --readsfile1={0}.fastq \\'.format(readsfile1)))
                            script_file_id.write( '{0}\n'.format('            --readsfile2={0}.fastq \\'.format(readsfile2)))
                script_file_id.write( '{0}\n'.format('            --verbose=NO \\'))
                script_file_id.write( '{0}\n'.format('            --trace=NO'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error indsdemultiplexing.py $RC; fi'))
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function run_readstrim_process'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Trimming reads ..."'))
                if format == 'FASTA':
                    if readtype == 'SE':
                        script_file_id.write( '{0}\n'.format('    ls {0}*.fasta > {1}'.format(demultiplexed_file, read_file_list)))
                    elif readtype == 'PE':
                        script_file_id.write( '{0}\n'.format('    ls {0}*-1.fasta > {1}'.format(demultiplexed_file, read_file_list)))
                elif format == 'FASTQ':
                    if readtype == 'SE':
                        script_file_id.write( '{0}\n'.format('    ls {0}*.fastq > {1}'.format(demultiplexed_file, read_file_list)))
                    elif readtype == 'PE':
                        script_file_id.write( '{0}\n'.format('    ls {0}*-1.fastq > {1}'.format(demultiplexed_file, read_file_list)))
                script_file_id.write( '{0}\n'.format('    while read FILE_1; do'))
                script_file_id.write( '{0}\n'.format('        if [[ $FILE_1 =~ .*errors.* ]]; then continue; fi'))
                script_file_id.write( '{0}\n'.format('        echo "$SEP"'))
                script_file_id.write( '{0}\n'.format('        echo "... file $FILE_1 ..."'))
                if readtype == 'SE':
                    script_file_id.write( '{0}\n'.format('        FILE_2=NONE'))
                elif readtype == 'PE':
                    if format == 'FASTA':
                        script_file_id.write( '{0}\n'.format('        FILE_2=`echo $FILE_1 | sed "s/-1.fasta/-2.fasta/g"`'))
                    elif format == 'FASTQ':
                        script_file_id.write( '{0}\n'.format('        FILE_2=`echo $FILE_1 | sed "s/-1.fastq/-2.fastq/g"`'))
                if format == 'FASTA':
                    script_file_id.write( '{0}\n'.format('        FILE_TRIMMED=`echo $FILE_1 | sed "s/-1.fasta/-trimmed/g"`'))
                elif format == 'FASTQ':
                    script_file_id.write( '{0}\n'.format('        FILE_TRIMMED=`echo $FILE_1 | sed "s/-1.fastq/-trimmed/g"`'))
                script_file_id.write( '{0}\n'.format('        /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('            --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('            readstrim.py \\'))
                script_file_id.write( '{0}\n'.format('                --technique={0} \\'.format(technique)))
                script_file_id.write( '{0}\n'.format('                --format={0} \\'.format(format)))
                script_file_id.write( '{0}\n'.format('                --readtype={0} \\'.format(readtype)))
                script_file_id.write( '{0}\n'.format('                --endsfile={0} \\'.format(endsfile)))
                script_file_id.write( '{0}\n'.format('                --index1len={0} \\'.format(index1len)))
                script_file_id.write( '{0}\n'.format('                --index2len={0} \\'.format(index2len)))
                script_file_id.write( '{0}\n'.format('                --dbrlen={0} \\'.format(dbrlen)))
                script_file_id.write( '{0}\n'.format('                --wend={0} \\'.format(wend)))
                script_file_id.write( '{0}\n'.format('                --cend={0} \\'.format(cend)))
                script_file_id.write( '{0}\n'.format('                --readsfile1=$FILE_1 \\'))
                if readtype == 'SE':
                    script_file_id.write( '{0}\n'.format('                --readsfile2=NONE \\'))
                elif readtype == 'PE':
                    script_file_id.write( '{0}\n'.format('                --readsfile2=$FILE_2 \\'))
                script_file_id.write( '{0}\n'.format('                --trimfile=$FILE_TRIMMED \\'))
                script_file_id.write( '{0}\n'.format('                --verbose=NO \\'))
                script_file_id.write( '{0}\n'.format('                --trace=NO'))
                script_file_id.write( '{0}\n'.format('            RC=$?'))
                script_file_id.write( '{0}\n'.format('            if [ $RC -ne 0 ]; then manage_error readstrim.py $RC; fi'))
                script_file_id.write( '{0}\n'.format('    done < {0}'.format(read_file_list)))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_ddradseq_simulation_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_ddradseq_simulation_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_ddradseq_simulation_name())))
                script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_ddradseq_simulation_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('run_rsitesearch_process'))
                script_file_id.write( '{0}\n'.format('run_simddradseq_process'))
                if int(dbrlen) > 0:
                    script_file_id.write( '{0}\n'.format('run_pcrdupremoval_process'))
                script_file_id.write( '{0}\n'.format('run_indsdemultiplexing_process'))
                script_file_id.write( '{0}\n'.format('run_readstrim_process'))
                script_file_id.write( 'end\n')
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be created.'.format(get_ddradseq_simulation_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_ddradseq_simulation_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ddradseq_simulation_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_config_file():
    '''
    Get the ddRADseq simulation config file path.
    '''

    # assign the ddRADseq simulation config file path
    ddradseq_simulation_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_ddradseq_simulation_code())

    # return the ddRADseq simulation config file path
    return ddradseq_simulation_config_file

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_process_script():
    '''
    Get the ddRADseq simulation process script path in the local computer.
    '''

    # assign the ddRADseq simulation script path
    ddradseq_simulation_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_ddradseq_simulation_code())

    # return the ddRADseq simulation script path
    return ddradseq_simulation_process_script

#-------------------------------------------------------------------------------

def get_ddradseq_simulation_process_starter():
    '''
    Get the ddRADseq simulation process starter path in the local computer.
    '''

    # assign the ddRADseq simulation process starter path
    ddradseq_simulation_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ddradseq_simulation_code())

    # return the ddRADseq simulation starter path
    return ddradseq_simulation_process_starter

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

if __name__ == '__main__':
     print('This file contains functions related to the ddRADseqTools process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
