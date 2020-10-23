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
This file contains functions related to the starcode process used in both console
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

def create_starcode_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create starcode config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the starcode config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_starcode_config_file())):
            os.makedirs(os.path.dirname(get_starcode_config_file()))
        with open(get_starcode_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters and trimming sets of starcode and their meaning in "https://github.com/gui11aume/starcode".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the starcode parameters.\n')
            file_id.write( '[starcode parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('distance = AUTO', '# maximum Levenshtein distance for clustering (AUTO is computed as min(8, 2 + [median seq length]/30))'))
            file_id.write( '{0:<50} {1}\n'.format('spheres = NO', f'# perform sphere clustering algorithm instead of the default message passing algorithm: {get_spheres_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('cluster_ratio = 5', '# minimum sequence count ratio to cluster two matching sequences (it only applies to message passing algorithm)'))
            file_id.write( '{0:<50} {1}\n'.format('non_redundant = NO', f'# remove redundant sequences from the output: {get_non_redundant_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('print_clusters = NO', f'# add a third column to the output containing the sequences associated with each cluster: {get_print_clusters_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('seq_id = NO', f'# show the clustered sequence numbers (1-based) following the original input order: {get_seq_id_code_list_text()}'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'concatenate_files = NO', f'# concatenate files building one file (SE) or two files (PE): {get_concatenate_files_code_list_text()}'))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '# This section has the information of the first library.\n')
                file_id.write(f'[library-{i + 1}]\n')
                file_id.write( '{0:<50} {1}\n'.format(f'read_file_1 = {os.path.basename(file_1_list[i])}', '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format( 'read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format(f'read_file_2 = {os.path.basename(file_2_list[i])}', '# name of the - strand reads file in PE read type or NONE in SE case'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '# If there are more libraries, you have to repeat the section library-1 with the data of each file.\n')
                    file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_starcode_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_starcode_process(cluster_name, log, function=None):
    '''
    Run a starcode process.
    '''

    # initialize the control variable
    OK = True

    # get the starcode option dictionary
    starcode_option_dict = xlib.get_option_dict(get_starcode_config_file())

    # get the experiment identification
    experiment_id = starcode_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the starcode config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_starcode_name()} config file ...\n')
    (OK, error_list) = check_starcode_config_file(strict=True)
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

    # check the starcode is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_starcode_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_starcode_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_starcode_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_starcode_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the starcode process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_starcode_process_script()} ...\n')
        (OK, error_list) = build_starcode_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the starcode process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_starcode_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_starcode_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_starcode_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the starcode process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_starcode_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_starcode_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the starcode process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_starcode_process_starter()} ...\n')
        (OK, error_list) = build_starcode_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the starcode process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_starcode_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_starcode_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_starcode_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the starcode process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_starcode_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_starcode_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the starcode process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_starcode_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_starcode_process_starter()), log)

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

def check_starcode_config_file(strict):
    '''
    Check the starcode config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        starcode_option_dict = xlib.get_option_dict(get_starcode_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in starcode_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = starcode_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = starcode_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "starcode parameters"
        if 'starcode parameters' not in sections_list:
            error_list.append('*** ERROR: the section "starcode parameters" is not found.')
            OK = False
        else:

            # check section "starcode parameters" - key "threads"
            threads = starcode_option_dict.get('starcode parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "starcode parameters" - key "distance"
            distance = starcode_option_dict.get('starcode parameters', {}).get('distance', not_found)
            if distance == not_found:
                error_list.append('*** ERROR: the key "distance" is not found in the section "starcode parameters".')
                OK = False
            elif distance.upper() != 'AUTO' and not xlib.check_int(distance, minimum=1):
                error_list.append('*** ERROR: the key "distance" has to be an integer number greater than or equal to 1 or AUTO.')
                OK = False

            # check section "starcode parameters" - key "spheres"
            spheres = starcode_option_dict.get('starcode parameters', {}).get('spheres', not_found)
            if spheres == not_found:
                error_list.append('*** ERROR: the key "spheres" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_code(spheres, get_spheres_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "spheres" has to be {get_spheres_code_list_text()}.')
                OK = False

            # check section "starcode parameters" - key "cluster_ratio"
            cluster_ratio = starcode_option_dict.get('starcode parameters', {}).get('cluster_ratio', not_found)
            if cluster_ratio == not_found:
                error_list.append('*** ERROR: the key "cluster_ratio" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_int(cluster_ratio, minimum=1):
                error_list.append('*** ERROR: the key "cluster_ratio" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "starcode parameters" - key "non_redundant"
            non_redundant = starcode_option_dict.get('starcode parameters', {}).get('non_redundant', not_found)
            if non_redundant == not_found:
                error_list.append('*** ERROR: the key "non_redundant" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_code(non_redundant, get_non_redundant_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "non_redundant" has to be {get_non_redundant_code_list_text()}.')
                OK = False

            # check section "starcode parameters" - key "print_clusters"
            print_clusters = starcode_option_dict.get('starcode parameters', {}).get('print_clusters', not_found)
            if print_clusters == not_found:
                error_list.append('*** ERROR: the key "print_clusters" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_code(print_clusters, get_print_clusters_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "print_clusters" has to be {get_print_clusters_code_list_text()}.')
                OK = False

            # check section "starcode parameters" - key "seq_id"
            seq_id = starcode_option_dict.get('starcode parameters', {}).get('seq_id', not_found)
            if seq_id == not_found:
                error_list.append('*** ERROR: the key "seq_id" is not found in the section "starcode parameters".')
                OK = False
            elif not xlib.check_code(seq_id, get_seq_id_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "seq_id" has to be {get_seq_id_code_list_text()}.')
                OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "format"
            format = starcode_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = starcode_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "read_type" has to be {get_read_type_code_list_text()}.')
                OK = False

            # check section "library" - key "concatenate_files"
            concatenate_files = starcode_option_dict.get('library', {}).get('concatenate_files', not_found)
            if concatenate_files == not_found:
                error_list.append('*** ERROR: the key "concatenate_files" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(concatenate_files, get_concatenate_files_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "concatenate_files" has to be {get_concatenate_files_code_list_text()}.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'starcode parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "readsfile1"
                    read_file_1 = starcode_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = starcode_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_starcode_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_starcode_process_script(cluster_name, current_run_dir):
    '''
    Build the current starcode process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the starcode option dictionary
    starcode_option_dict = xlib.get_option_dict(get_starcode_config_file())

    # get the options
    experiment_id = starcode_option_dict['identification']['experiment_id']
    read_dataset_id = starcode_option_dict['identification']['read_dataset_id']
    threads = starcode_option_dict['starcode parameters']['threads']
    distance = starcode_option_dict['starcode parameters']['distance'].upper()
    spheres = starcode_option_dict['starcode parameters']['spheres'].upper()
    cluster_ratio = starcode_option_dict['starcode parameters']['cluster_ratio']
    non_redundant = starcode_option_dict['starcode parameters']['non_redundant'].upper()
    print_clusters = starcode_option_dict['starcode parameters']['print_clusters'].upper()
    seq_id = starcode_option_dict['starcode parameters']['seq_id'].upper()
    format = starcode_option_dict['library']['format'].upper()
    read_type = starcode_option_dict['library']['read_type'].upper()
    concatenate_files = starcode_option_dict['library']['concatenate_files'].upper()

    # get the sections list
    sections_list = []
    for section in starcode_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # get the input read directory
    input_read_dir = xlib.get_cluster_experiment_read_dataset_dir(experiment_id, read_dataset_id)

    # build the file name list
    file_name_1_list = []
    file_name_2_list = []
    for section in sections_list:
        if re.match('^library-[0-9]+$', section):
            read_file_1 = starcode_option_dict[section]['read_file_1']
            file_name_1_list.append(f'{input_read_dir}/{read_file_1}')
            if read_type == 'PE':
                read_file_2 = starcode_option_dict[section]['read_file_2']
                file_name_2_list.append(f'{input_read_dir}/{read_file_2}')

    # write the starcode process script
    try:
        if not os.path.exists(os.path.dirname(get_starcode_process_script())):
            os.makedirs(os.path.dirname(get_starcode_process_script()))
        with open(get_starcode_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STARCODE_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_starcode_anaconda_code()}/bin\n')
            script_file_id.write( 'PATH=$STARCODE_PATH:$PATH\n')
            script_file_id.write(f'cd {xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write(f'source activate {xlib.get_starcode_anaconda_code()}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'STATUS_DIR={xlib.get_status_dir(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_OK={xlib.get_status_ok(current_run_dir)}\n')
            script_file_id.write(f'SCRIPT_STATUS_WRONG={xlib.get_status_wrong(current_run_dir)}\n')
            script_file_id.write( 'mkdir --parents $STATUS_DIR\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_OK ]; then rm $SCRIPT_STATUS_OK; fi\n')
            script_file_id.write( 'if [ -f $SCRIPT_STATUS_WRONG ]; then rm $SCRIPT_STATUS_WRONG; fi\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'CURRENT_DIR={current_run_dir}\n')
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
            if concatenate_files == 'NO':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function run_initial_starcode_process\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Running initial Starcode process"\n')
                script_file_id.write( '    starcode --version\n')
                for i in range(len(file_name_1_list)):
                        script_file_id.write( '    echo "$SEP"\n')
                        script_file_id.write( '    /usr/bin/time \\\n')
                        script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                        script_file_id.write( '        starcode \\\n')
                        script_file_id.write(f'            --threads={threads} \\\n')
                        if distance != 'AUTO':
                            script_file_id.write(f'            --distance {distance} \\\n')
                        if spheres == 'YES':
                            script_file_id.write( '            --spheres \\\n')
                        if spheres == 'NO':
                            script_file_id.write(f'            --cluster-ratio {cluster_ratio} \\\n')
                        if non_redundant == 'YES':
                            script_file_id.write( '            --non-redundant \\\n')
                        if print_clusters == 'YES':
                            script_file_id.write( '            --print-clusters \\\n')
                        if seq_id == 'YES':
                            script_file_id.write( '            --seq-id \\\n')
                        if read_type == 'SE':
                            script_file_id.write(f'            --input {file_name_1_list[i]} \\\n')
                        elif read_type == 'PE':
                            script_file_id.write(f'            --input1 {file_name_1_list[i]} \\\n')
                            script_file_id.write(f'            --input2 {file_name_2_list[i]} \\\n')
                        if read_type == 'PE' and non_redundant == 'YES':
                            script_file_id.write(f'            --output1 $CURRENT_DIR/starcode-initial-output-1-{i}.txt \\\n')
                            script_file_id.write(f'            --output2 $CURRENT_DIR/starcode-initial-output-2-{i}.txt\n')
                        else:
                            script_file_id.write(f'            --output $CURRENT_DIR/starcode-initial-output-{i}.txt\n')
                        script_file_id.write( '    RC=$?\n')
                        script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error starcode $RC; fi\n')
                script_file_id.write( '    echo "Starcode is processed"\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function concatenate_intermediate_starcode_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                if read_type == 'SE':
                    script_file_id.write( '    echo "Concatenate the intermediate Starcode files ..."\n')
                    output_file_list = []
                    for i in range(len(file_name_1_list)):
                        output_file_list.append(f'$CURRENT_DIR/starcode-initial-output-{i}.txt')
                    script_file_id.write(f'    cat {" ".join(output_file_list)} > starcode-initial-output-concatenated.txt\n')
                    script_file_id.write( '    echo "The concatenation is done."\n')
                elif read_type == 'PE':
                    pass
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_intermediate_fasta_file\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                if read_type == 'SE':
                    script_file_id.write( '    echo "Building the intermediate FASTA file ..."\n')
                    script_file_id.write( '    i=0\n')
                    script_file_id.write(f'    FASTA_FILE=$CURRENT_DIR/starcode-intermediate.fasta\n')
                    script_file_id.write( '    touch $FASTA_FILE\n')
                    script_file_id.write(f'    for STARCODE_RECORD in $(cut -f 1 $CURRENT_DIR/starcode-initial-output-concatenated.txt); do\n')
                    script_file_id.write( '        ((i++))\n')
                    script_file_id.write( '        echo ">seq$i" >>$FASTA_FILE \n')
                    script_file_id.write( '        echo "$STARCODE_RECORD" >>$FASTA_FILE \n')
                    script_file_id.write( '    done\n')
                    script_file_id.write( '    echo "The files are built."\n')
                elif read_type == 'PE':
                    pass
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function run_final_starcode_process\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                if read_type == 'SE':
                    script_file_id.write( '    echo "Running final Starcode process"\n')
                    script_file_id.write( '    echo "$SEP"\n')
                    script_file_id.write( '    /usr/bin/time \\\n')
                    script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                    script_file_id.write( '        starcode \\\n')
                    script_file_id.write(f'            --threads={threads} \\\n')
                    if distance != 'AUTO':
                        script_file_id.write(f'            --distance {distance} \\\n')
                    if spheres == 'YES':
                        script_file_id.write( '            --spheres \\\n')
                    if spheres == 'NO':
                        script_file_id.write(f'            --cluster-ratio {cluster_ratio} \\\n')
                    if non_redundant == 'YES':
                        script_file_id.write( '            --non-redundant \\\n')
                    if print_clusters == 'YES':
                        script_file_id.write( '            --print-clusters \\\n')
                    if seq_id == 'YES':
                        script_file_id.write( '            --seq-id \\\n')
                    script_file_id.write(f'            --input starcode-intermediate.fasta \\\n')
                    script_file_id.write(f'            --output $CURRENT_DIR/starcode-final-output.txt\n')
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error starcode $RC; fi\n')
                    script_file_id.write( '    echo "Starcode is processed"\n')
                elif read_type == 'PE':
                    pass
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_final_fasta_file\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                if read_type == 'SE':
                    script_file_id.write( '    echo "Building the final FASTA file ..."\n')
                    script_file_id.write( '    i=0\n')
                    script_file_id.write( '    FASTA_FILE=$CURRENT_DIR/starcode.fasta\n')
                    script_file_id.write( '    touch $FASTA_FILE\n')
                    script_file_id.write( '    for STARCODE_RECORD in $(cut -f 1 $CURRENT_DIR/starcode-final-output.txt); do\n')
                    script_file_id.write( '        ((i++))\n')
                    script_file_id.write( '        echo ">seq$i" >>$FASTA_FILE \n')
                    script_file_id.write( '        echo "$STARCODE_RECORD" >>$FASTA_FILE \n')
                    script_file_id.write( '    done\n')
                    script_file_id.write( '    echo "The file is built."\n')
                elif read_type == 'PE':
                    pass
                script_file_id.write( '}\n')
            elif concatenate_files == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function concatenate_input_files\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Concatenate the files of the library ..."\n')
                if format == 'FASTQ':
                    concatenated_library_1 = f'$CURRENT_DIR/concatenated_library_1.fastq'
                elif format == 'FASTA':
                    concatenated_library_1 = f'$CURRENT_DIR/concatenated_library_1.fasta'
                script_file_id.write(f'    cat {" ".join(file_name_1_list)} > {concatenated_library_1}\n')
                file_name_1_list = [concatenated_library_1]
                if read_type == 'PE':
                    if format == 'FASTQ':
                        concatenated_library_2 = f'$CURRENT_DIR/concatenated_library_2.fastq'
                    elif format == 'FASTA':
                        concatenated_library_2 = f'$CURRENT_DIR/concatenated_library_2.fasta'
                    script_file_id.write(f'    cat {" ".join(file_name_2_list)} > {concatenated_library_2}\n')
                    file_name_2_list = [concatenated_library_2]
                script_file_id.write( '    echo "The concatenation is done."\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function run_starcode_process\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Running Starcode process"\n')
                script_file_id.write( '    starcode --version\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '        starcode \\\n')
                script_file_id.write(f'            --threads={threads} \\\n')
                if distance != 'AUTO':
                    script_file_id.write(f'            --distance {distance} \\\n')
                if spheres == 'YES':
                    script_file_id.write( '            --spheres \\\n')
                if spheres == 'NO':
                    script_file_id.write(f'            --cluster-ratio {cluster_ratio} \\\n')
                if non_redundant == 'YES':
                    script_file_id.write( '            --non-redundant \\\n')
                if print_clusters == 'YES':
                    script_file_id.write( '            --print-clusters \\\n')
                if seq_id == 'YES':
                    script_file_id.write( '            --seq-id \\\n')
                if read_type == 'SE':
                    script_file_id.write(f'            --input {concatenated_library_1} \\\n')
                elif read_type == 'PE':
                    script_file_id.write(f'            --input1 {concatenated_library_1} \\\n')
                    script_file_id.write(f'            --input2 {concatenated_library_2} \\\n')
                if read_type == 'PE' and non_redundant == 'YES':
                    script_file_id.write(f'            --output1 $CURRENT_DIR/starcode-output-1.txt \\\n')
                    script_file_id.write(f'            --output2 $CURRENT_DIR/starcode-output-2.txt\n')
                else:
                    script_file_id.write(f'            --output $CURRENT_DIR/starcode-output.txt\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error starcode $RC; fi\n')
                script_file_id.write( '    echo "Starcode is processed"\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_fasta_file\n')
                script_file_id.write( '{\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                if read_type == 'SE':
                    script_file_id.write( '    echo "Building the final FASTA file ..."\n')
                    script_file_id.write( '    i=0\n')
                    script_file_id.write( '    FASTA_FILE=$CURRENT_DIR/starcode.fasta\n')
                    script_file_id.write( '    touch $FASTA_FILE\n')
                    script_file_id.write( '    for STARCODE_RECORD in $(cut -f 1 $CURRENT_DIR/starcode-output.txt); do\n')
                    script_file_id.write( '        ((i++))\n')
                    script_file_id.write( '        echo ">seq$i" >>$FASTA_FILE \n')
                    script_file_id.write( '        echo "$STARCODE_RECORD" >>$FASTA_FILE \n')
                    script_file_id.write( '    done\n')
                    script_file_id.write( '    echo "The file is built."\n')
                elif read_type == 'PE':
                    pass
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
            process_name = f'{xlib.get_starcode_name()} process'
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
            if concatenate_files == 'NO':
                script_file_id.write( 'run_initial_starcode_process\n')
                script_file_id.write( 'concatenate_intermediate_starcode_files\n')
                script_file_id.write( 'build_intermediate_fasta_file\n')
                script_file_id.write( 'run_final_starcode_process\n')
                script_file_id.write( 'build_final_fasta_file\n')
            elif concatenate_files == 'YES':
                script_file_id.write( 'concatenate_input_files\n')
                script_file_id.write( 'run_starcode_process\n')
                script_file_id.write( 'build_fasta_file\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_starcode_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_starcode_process_starter(current_run_dir):
    '''
    Build the starter of the current starcode process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the starcode process starter
    try:
        if not os.path.exists(os.path.dirname(get_starcode_process_starter())):
            os.makedirs(os.path.dirname(get_starcode_process_starter()))
        with open(get_starcode_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_starcode_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_starcode_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_starcode_config_file():
    '''
    Get the starcode config file path.
    '''

    # assign the starcode config file path
    starcode_config_file = f'{xlib.get_config_dir()}/{xlib.get_starcode_code()}-config.txt'

    # return the starcode config file path
    return starcode_config_file

#-------------------------------------------------------------------------------

def get_starcode_process_script():
    '''
    Get the starcode process script path in the local computer.
    '''

    # assign the starcode script path
    starcode_process_script = f'{xlib.get_temp_dir()}/{xlib.get_starcode_code()}-process.sh'

    # return the starcode script path
    return starcode_process_script

#-------------------------------------------------------------------------------

def get_starcode_process_starter():
    '''
    Get the starcode process starter path in the local computer.
    '''

    # assign the starcode process starter path
    starcode_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_starcode_code()}-process-starter.sh'

    # return the starcode starter path
    return starcode_process_starter

#-------------------------------------------------------------------------------
    
def get_spheres_code_list():
    '''
    Get the code list of "spheres".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_spheres_code_list_text():
    '''
    Get the code list of "spheres" as text.
    '''

    return str(get_spheres_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_non_redundant_code_list():
    '''
    Get the code list of "non_redundant".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_non_redundant_code_list_text():
    '''
    Get the code list of "non_redundant" as text.
    '''

    return str(get_non_redundant_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_print_clusters_code_list():
    '''
    Get the code list of "print_clusters".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_print_clusters_code_list_text():
    '''
    Get the code list of "print_clusters" as text.
    '''

    return str(get_print_clusters_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_seq_id_code_list():
    '''
    Get the code list of "seq_id".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_seq_id_code_list_text():
    '''
    Get the code list of "seq_id" as text.
    '''

    return str(get_seq_id_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_concatenate_files_code_list():
    '''
    Get the code list of "concatenate_files".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_concatenate_files_code_list_text():
    '''
    Get the code list of "concatenate_files" as text.
    '''

    return str(get_concatenate_files_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
    Get the code list of "read_type".
    '''

    # -- return ['SE', 'PE']
    return ['SE']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "read_type" as text.
    '''

    # -- return 'SE (single-end) or PE (pair-end)'
    return 'SE (single-end)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the starcode process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
