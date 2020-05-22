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
This file contains functions related to the Transrate process used in both console
mode and gui mode.
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

def is_installed_transrate(cluster_name, passed_connection, ssh_client):
    '''
    Check if Transrate is installed.
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

    # check the Transrate directory is created
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_transrate_name()} ] && echo RC=0 || echo RC=1'
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

def install_transrate(cluster_name, log, function=None):
    '''
    Install the Transrate software in the cluster.
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
            log.write(f'You have to link a volume in the mounting point {xlib.get_cluster_app_dir()} for the template {cluster_name}.\n')
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Installation requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir('installation', xlib.get_transrate_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Transrate installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the installation script {get_transrate_installation_script()} ...\n')
        (OK, error_list) = build_transrate_installation_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Transrate installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the installation script {get_transrate_installation_script()} in the directory {current_run_dir} of the master ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transrate_installation_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transrate_installation_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Transrate installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transrate_installation_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transrate_installation_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Transrate installation starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_transrate_installation_starter()} ...\n')
        (OK, error_list) = build_transrate_installation_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Transrate installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_transrate_installation_starter()} in the directory {current_run_dir} of the master ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transrate_installation_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transrate_installation_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Transrate installation starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transrate_installation_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transrate_installation_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Transrate installation
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_transrate_installation_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_transrate_installation_starter()), log)

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

def build_transrate_installation_script(cluster_name, current_run_dir):
    '''
    Build the Transrate installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the version and download URL of Transrate
    (transrate_version, transrate_url) = xconfiguration.get_bioinfo_app_data(xlib.get_transrate_name())

    # write the Transrate installation script
    try:
        if not os.path.exists(os.path.dirname(get_transrate_installation_script())):
            os.makedirs(os.path.dirname(get_transrate_installation_script()))
        with open(get_transrate_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write(f'PYTHON3_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
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
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( f'function remove_transrate_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing {xlib.get_transrate_name()} directory ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    if [ -d "{xlib.get_transrate_name()}" ]; then\n')
            script_file_id.write(f'        rm -rf {xlib.get_transrate_name()}\n')
            script_file_id.write( '        echo "The directory is removed."\n')
            script_file_id.write( '    else\n')
            script_file_id.write( '        echo "The directory is not found."\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function download_transrate_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Downloading the {xlib.get_transrate_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            download_script = f'import requests; r = requests.get(\'{transrate_url}\') ; open(\'{xlib.get_transrate_name()}.tar.gz\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $PYTHON3_PATH/python3 -c "{download_script}"\n')
            script_file_id.write(f'    RC=$?\n')
            script_file_id.write(f'    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write(f'    echo "The file is downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function decompress_transrate_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Decompressing the {xlib.get_transrate_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    tar -xzvf {xlib.get_transrate_name()}.tar.gz\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write( '    echo "The file is decompressed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_transrate_directory\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Renaming the {xlib.get_transrate_name()} directory ..."\n'.format(''.format()))
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    mv {xlib.get_transrate_code()}-{transrate_version}-linux-x86_64 {xlib.get_transrate_name()}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
            script_file_id.write( '    echo "The directory is renamed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function remove_transrate_installation_file\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Removing the {xlib.get_transrate_name()} installation file ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    rm -f {xlib.get_transrate_name()}.tar.gz\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "The file is removed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_dependencies\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    echo "Installing {xlib.get_transrate_name()} dependencies ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_transrate_name()}\n')
            script_file_id.write( '    ./transrate --install-deps ref\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error transrate $RC; fi\n')
            script_file_id.write( '    echo "Dependencies are installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function rename_librt\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Renaming library librt.so.1 ..."\n')
            script_file_id.write(f'    cd {xlib.get_cluster_app_dir()}/{xlib.get_transrate_name()}/bin\n')
            script_file_id.write( '    mv librt.so.1 librt.so.1.wrong\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then mv $RC; fi\n')
            script_file_id.write( '    echo "librt.so.1 is renamed."\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {xlib.get_rnaquast_name()} installation"\n')
            message = xlib.get_mail_message_ok(f'{xlib.get_rnaquast_name()} installation', cluster_name)
            script_file_id.write(f'    MESSAGE="{message}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {xlib.get_rnaquast_name()} installation"\n')
            message = xlib.get_mail_message_wrong(f'{xlib.get_rnaquast_name()} installation', cluster_name)
            script_file_id.write(f'    MESSAGE="{message}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
            script_file_id.write( '    exit 3\n')
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
            script_file_id.write( 'remove_transrate_directory\n')
            script_file_id.write( 'download_transrate_installation_file\n')
            script_file_id.write( 'decompress_transrate_installation_file\n')
            script_file_id.write( 'rename_transrate_directory\n')
            script_file_id.write( 'remove_transrate_installation_file\n')
            script_file_id.write( 'install_dependencies\n')
            script_file_id.write( 'rename_librt\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_transrate_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transrate_installation_starter(current_run_dir):
    '''
    Build the starter of the Transrate installation.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Transrate installation starter
    try:
        if not os.path.exists(os.path.dirname(get_transrate_installation_starter())):
            os.makedirs(os.path.dirname(get_transrate_installation_starter()))
        with open(get_transrate_installation_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_transrate_installation_script())} &>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_transrate_installation_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transrate_installation_script():
    '''
    Get the Transrate installation path in the local computer.
    '''

    # assign the Transrate installation path
    transrate_installation_script = f'{xlib.get_temp_dir()}/{xlib.get_transrate_name()}-installation.sh'

    # return the Transrate installation path
    return transrate_installation_script

#-------------------------------------------------------------------------------

def get_transrate_installation_starter():
    '''
    Get the Transrate installation starter path in the local computer.
    '''

    # assign the Transrate installation starter path
    transrate_installation_starter = f'{xlib.get_temp_dir()}/{xlib.get_transrate_name()}-installation-starter.sh'

    # return the Transrate installation starter path
    return transrate_installation_starter

#-------------------------------------------------------------------------------

def create_transrate_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='GCF_000001735.3_TAIR10_genomic.fna', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq'], assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create Transrate config file with the default options. It is necessary
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

    # create the Transrate config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_transrate_config_file())):
            os.makedirs(os.path.dirname(get_transrate_config_file()))
        with open(get_transrate_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id'))
            file_id.write( '{0}\n'.format(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id'))
            file_id.write( '{0}\n'.format(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id'))
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file_name, read_dataset_id and assembly_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of Transrate and their meaning in http://hibberdlab.com/transrate/.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# WARNING: The files have to be decompressed.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the Transrate parameters'))
            file_id.write( '{0}\n'.format('[Transrate parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('loglevel = INFO', f'# log level: {get_loglevel_code_list_text()}'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first library.'))
                file_id.write( '{0}\n'.format(f'[library-{i + 1}]'))
                file_id.write( '{0:<50} {1}\n'.format(f'read_file_1 = {os.path.basename(file_1_list[i])}', '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format(f'read_file_2 = {os.path.basename(file_2_list[i])}', '# name of the - strand reads file in PE read type or NONE in SE case'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more libraries, you have to repeat the section library-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_transrate_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_transrate_process(cluster_name, log, function=None):
    '''
    Run a Transrate process.
    '''

    # initialize the control variable
    OK = True

    # get the Transrate option dictionary
    transrate_option_dict = xlib.get_option_dict(get_transrate_config_file())

    # get the experiment identification
    experiment_id = transrate_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the Transrate config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_transrate_name()} config file ...\n')
    (OK, error_list) = check_transrate_config_file(strict=True)
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

    # check the Transrate is installed
    if OK:
        command = f'[ -d {xlib.get_cluster_app_dir()}/{xlib.get_transrate_name()} ] && echo RC=0 || echo RC=1'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write(f'*** ERROR: {xlib.get_transrate_name()} is not installed.\n')
            OK = False
        # -- (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_transrate_bioconda_code(), cluster_name, True, ssh_client)
        # -- if OK:
        # --     if not is_installed:
        # --         log.write(f'*** ERROR: {xlib.get_busco_name()} is not installed.\n')
        # --         OK = False
        # -- else:
        # --     log.write(f'*** ERROR: The verification of {xlib.get_busco_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_transrate_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Transrate process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_transrate_process_script()} ...\n')
        (OK, error_list) = build_transrate_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Transrate process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_transrate_process_script()} in the directory {current_run_dir} of the master ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transrate_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transrate_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Transrate process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transrate_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transrate_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Transrate process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_transrate_process_starter()} ...\n')
        (OK, error_list) = build_transrate_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Transrate process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_transrate_process_starter()} in the directory {current_run_dir} of the master ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_transrate_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_transrate_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Transrate process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_transrate_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_transrate_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Transrate process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_transrate_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_transrate_process_starter()), log)

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

def check_transrate_config_file(strict):
    '''
    Check the Transrate config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        transrate_option_dict = xlib.get_option_dict(get_transrate_config_file())
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in transrate_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = transrate_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = transrate_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = transrate_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False
            elif reference_file.find('.gz') != -1:
                error_list.append('*** ERROR: the key "reference_file" has to be a decompressed file (.gz).')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = transrate_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = transrate_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = transrate_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" has to start with {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = transrate_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
                    OK = False

        # check section "Transrate parameters"
        if 'Transrate parameters' not in sections_list:
            error_list.append('*** ERROR: the section "Transrate parameters" is not found.')
            OK = False
        else:

            # check section "Transrate parameters" - key "threads"
            threads = transrate_option_dict.get('Transrate parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "rnaQUAST parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "Transrate parameters" - key "loglevel"
            loglevel = transrate_option_dict.get('Transrate parameters', {}).get('loglevel', not_found)
            if loglevel == not_found:
                error_list.append('*** ERROR: the key "loglevel" is not found in the section "Transrate parameters".')
                OK = False
            elif not xlib.check_code(loglevel, get_loglevel_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "loglevel" has to be {get_loglevel_code_list_text()}.')
                OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "format"
            format = transrate_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = transrate_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "read_type" has to be {get_read_type_code_list_text()}.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'Transrate parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = transrate_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False
                    elif read_file_1.find('.gz') != -1:
                        error_list.append(f'*** ERROR: the key "read_file_1" in the section "{section}" has to be a decompressed file (.gz).')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = transrate_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False
                    elif read_file_2.find('.gz') != -1:
                        error_list.append(f'*** ERROR: the key "read_file_2" in the section "{section}" has to be a decompressed file (.gz).')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_transrate_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transrate_process_script(cluster_name, current_run_dir):
    '''
    Build the current Transrate process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the Transrate option dictionary
    transrate_option_dict = xlib.get_option_dict(get_transrate_config_file())

    # get the options
    experiment_id = transrate_option_dict['identification']['experiment_id']
    reference_dataset_id = transrate_option_dict['identification']['reference_dataset_id']
    reference_file = transrate_option_dict['identification']['reference_file']
    assembly_software = transrate_option_dict['identification']['assembly_software']
    read_dataset_id = transrate_option_dict['identification']['read_dataset_id']
    assembly_dataset_id = transrate_option_dict['identification']['assembly_dataset_id']
    assembly_type = transrate_option_dict['identification']['assembly_type']
    threads = transrate_option_dict['Transrate parameters']['threads']
    loglevel = transrate_option_dict['Transrate parameters']['loglevel']
    read_type = transrate_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in transrate_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build library files
    files1 = ''
    files2 = ''
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = transrate_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            files1 += read_file_1 + ','
            if read_type == 'PE':
                read_file_2 = transrate_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                files2 += read_file_2 + ','
    files1 = files1[:len(files1) - 1]
    if read_type == 'PE':
        files2 = files2[:len(files2) - 1]

    # set the reference file path
    if reference_dataset_id.upper() != 'NONE':
        reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

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

    # write the Transrate process script
    try:
        if not os.path.exists(os.path.dirname(get_transrate_process_script())):
            os.makedirs(os.path.dirname(get_transrate_process_script()))
        with open(get_transrate_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write(f'TRANSRATE_PATH={xlib.get_cluster_app_dir()}/{xlib.get_transrate_name()}\n')
            # -- script_file_id.write(f'TRANSRATE_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_transrate_bioconda_code()}/bin\n')
            script_file_id.write( 'export PATH=$TRANSRATE_PATH:$PATH\n')
            # -- script_file_id.write(f'source {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin/activate {xlib.get_transrate_bioconda_code()}\n')
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
            script_file_id.write(f'    echo "HOST_IP: $HOST_IP - HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_transrate_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write( '        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\\n')
            script_file_id.write( '        transrate \\\n')
            script_file_id.write(f'            --threads={threads} \\\n')
            script_file_id.write(f'            --assembly={transcriptome_file} \\\n')
            if reference_dataset_id.upper() != 'NONE':
                script_file_id.write(f'            --reference={reference_file} \\\n')
            script_file_id.write(f'            --left={read_file_1} \\\n')
            if read_type.upper() == 'PE':
                script_file_id.write(f'            --right={read_file_2} \\\n')
            script_file_id.write(f'            --loglevel={loglevel.lower()} \\\n')
            script_file_id.write(f'            --output={current_run_dir}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error transrate $RC; fi\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {xlib.get_transrate_name()} process"\n')
            script_file_id.write(f'    MESSAGE="{xlib.get_mail_message_ok(xlib.get_transrate_name(), cluster_name)}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: {xlib.get_transrate_name()} process"\n')
            script_file_id.write(f'    MESSAGE="{xlib.get_mail_message_wrong(xlib.get_transrate_name(), cluster_name)}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '    touch $SCRIPT_STATUS_WRONG\n')
            script_file_id.write( '    exit 3\n')
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
            script_file_id.write( 'run_transrate_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_transrate_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_transrate_process_starter(current_run_dir):
    '''
    Build the starter of the current Transrate process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Transrate process starter
    try:
        if not os.path.exists(os.path.dirname(get_transrate_process_starter())):
            os.makedirs(os.path.dirname(get_transrate_process_starter()))
        with open(get_transrate_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_transrate_process_script())} &>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_transrate_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_transrate_config_file():
    '''
    Get the Transrate config file path.
    '''

    # assign the Transrate config file path
    transrate_config_file = f'{xlib.get_config_dir()}/{xlib.get_transrate_code()}-config.txt'

    # return the Transrate config file path
    return transrate_config_file

#-------------------------------------------------------------------------------

def get_transrate_process_script():
    '''
    Get the Transrate process script path in the local computer.
    '''

    # assign the Transrate script path
    transrate_process_script = f'{xlib.get_temp_dir()}/{xlib.get_transrate_code()}-process.sh'

    # return the Transrate script path
    return transrate_process_script

#-------------------------------------------------------------------------------

def get_transrate_process_starter():
    '''
    Get the Transrate process starter path in the local computer.
    '''

    # assign the Transrate process starter path
    transrate_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_transrate_code()}-process-starter.sh'

    # return the Transrate starter path
    return transrate_process_starter

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
    
def get_loglevel_code_list():
    '''
    Get the code list of "loglevel".
    '''

    return ['ERROR', 'INFO', 'WARN', 'DEBUG']

#-------------------------------------------------------------------------------
    
def get_loglevel_code_list_text():
    '''
    Get the code list of "loglevel" as text.
    '''

    return str(get_loglevel_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_format_code_list():
    '''
    Get the code list of "format".
    '''

    return ['FASTQ']

#-------------------------------------------------------------------------------
    
def get_format_code_list_text():
    '''
    Get the code list of "format" as text.
    '''

    return 'FASTQ (only this format is allowed)'

#-------------------------------------------------------------------------------
    
def get_read_type_code_list():
    '''
    Get the code list of "read_type".
    '''

    return ['SE', 'PE']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "read_type" as text.
    '''

    return 'SE (single-end) or PE (pair-end)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the Transrate process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
