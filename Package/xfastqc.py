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
This file contains functions related to the FastQC process used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys

import xbioinfoapp
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_fastqc_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), file_list=['rnaseq-a_1.fastq', 'rnaseq-a_2.fastq']):
    '''
    Create FastQC config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the FastQC config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_fastqc_config_file())):
            os.makedirs(os.path.dirname(get_fastqc_config_file()))
        with open(get_fastqc_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of FastQC and their meaning in "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the FastQC parameters\n')
            file_id.write( '[FastQC parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use'))
            for i in range(len(file_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '# This section has the information of the first file.\n')
                file_id.write(f'[file-{i + 1}]\n')
                file_id.write( '{0:<50} {1}\n'.format(f'file_name = {os.path.basename(file_list[i])}', '# read file name'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '# If there are more files, you have to repeat the section file-1 with the data of each file.\n')
                    file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_fastqc_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_fastqc_process(cluster_name, log, function=None):
    '''
    Run a FastQC process.
    '''

    # initialize the control variable
    OK = True

    # get the FastQC option dictionary
    fastqc_option_dict = xlib.get_option_dict(get_fastqc_config_file())

    # get the experiment identification
    experiment_id = fastqc_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the FastQC config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_fastqc_name()} config file ...\n')
    (OK, error_list) = check_fastqc_config_file(strict=True)
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

    # check the FastQC is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_fastqc_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_fastqc_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_fastqc_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_fastqc_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the FastQC process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_fastqc_process_script()} ...\n')
        (OK, error_list) = build_fastqc_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the FastQC process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_fastqc_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_fastqc_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_fastqc_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the FastQC process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_fastqc_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_fastqc_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the FastQC process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_fastqc_process_starter()} ...\n')
        (OK, error_list) = build_fastqc_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the FastQC process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_fastqc_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_fastqc_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_fastqc_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the FastQC process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_fastqc_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_fastqc_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the FastQC process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_fastqc_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_fastqc_process_starter()), log)

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

def check_fastqc_config_file(strict):
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
        fastqc_option_dict = xlib.get_option_dict(get_fastqc_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in fastqc_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = fastqc_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = fastqc_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "FastQC parameters"
        if 'FastQC parameters' not in sections_list:
            error_list.append('*** ERROR: the section "FastQC parameters" is not found.')
            OK = False
        else:

            # check section "FastQC parameters" - key "threads"
            threads = fastqc_option_dict.get('FastQC parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

        # check section "file-1"
        if 'file-1' not in sections_list:
            error_list.append('*** ERROR: the section "file-1" is not found.')
            OK = False

        # check all sections "file-n"
        for section in sections_list:

            if section not in ['identification', 'FastQC parameters']:

                # check than the section identification is like file-n 
                if not re.match('^file-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "file-n" - key "file_name"
                    file_name = fastqc_option_dict.get(section, {}).get('file_name', not_found)
                    if file_name == not_found:
                        error_list.append(f'*** ERROR: the key "file_name" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.is_valid_path(file_name, 'linux'):
                        error_list.append(f'*** ERROR: the file {file_name} in the key "file_name" of the section "{section}" has a non valid file name.')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_fastqc_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_fastqc_process_script(cluster_name, current_run_dir):
    '''
    Build the current FastQC process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the FastQC option dictionary
    fastqc_option_dict = xlib.get_option_dict(get_fastqc_config_file())

    # get the options
    experiment_id = fastqc_option_dict['identification']['experiment_id']
    read_dataset_id = fastqc_option_dict['identification']['read_dataset_id']
    threads = fastqc_option_dict['FastQC parameters']['threads']

    # get the sections list
    sections_list = []
    for section in fastqc_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build the file name list
    file_name_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^file-[0-9]+$', section):
            file_name = fastqc_option_dict[section]['file_name']
            file_name_list.append(file_name)

    # write the FastQC process script
    try:
        if not os.path.exists(os.path.dirname(get_fastqc_process_script())):
            os.makedirs(os.path.dirname(get_fastqc_process_script()))
        with open(get_fastqc_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function run_fastqc_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_fastqc_anaconda_code()}\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    fastqc --version\n')
            for file_name in file_name_list:
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '        fastqc \\\n')
                script_file_id.write(f'            {xlib.get_cluster_read_file(experiment_id, read_dataset_id, file_name)} \\\n')
                script_file_id.write(f'            --threads={threads} \\\n')
                script_file_id.write(f'            --outdir={current_run_dir}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error fastqc $RC; fi\n')
            script_file_id.write( '    conda deactivate\n')
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
            process_name = f'{xlib.get_fastqc_name()} process'
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
            script_file_id.write( 'run_fastqc_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_fastqc_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_fastqc_process_starter(current_run_dir):
    '''
    Build the starter of the current FastQC process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the FastQC process starter
    try:
        if not os.path.exists(os.path.dirname(get_fastqc_process_starter())):
            os.makedirs(os.path.dirname(get_fastqc_process_starter()))
        with open(get_fastqc_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_fastqc_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_fastqc_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_fastqc_config_file():
    '''
    Get the FastQC config file path.
    '''

    # assign the FastQC config file path
    fastqc_config_file = f'{xlib.get_config_dir()}/{xlib.get_fastqc_code()}-config.txt'

    # return the FastQC config file path
    return fastqc_config_file

#-------------------------------------------------------------------------------

def get_fastqc_process_script():
    '''
    Get the FastQC process script path in the local computer.
    '''

    # assign the FastQC script path
    fastqc_process_script = f'{xlib.get_temp_dir()}/{xlib.get_fastqc_code()}-process.sh'

    # return the FastQC script path
    return fastqc_process_script

#-------------------------------------------------------------------------------

def get_fastqc_process_starter():
    '''
    Get the FastQC process starter path in the local computer.
    '''

    # assign the FastQC process starter path
    fastqc_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_fastqc_code()}-process-starter.sh'

    # return the FastQC starter path
    return fastqc_process_starter

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the FastQC process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
