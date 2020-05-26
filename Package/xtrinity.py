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
This file contains functions related to the Trinity process used in both
console mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import subprocess
import sys

import xbioinfoapp
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_trinity_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create Trinity config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the Trinity config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_trinity_config_file())):
            os.makedirs(os.path.dirname(get_trinity_config_file()))
        with open(get_trinity_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of Trinity and their meaning in "https://github.com/trinityrnaseq/trinityrnaseq/wiki".\n')
            file_id.write( '#\n')
            file_id.write( '# There are two formats to set an option:\n')
            file_id.write( '#\n')
            file_id.write( '#    option = value                             <- if the option supports a single value\n')
            file_id.write( '#\n')
            file_id.write( '#    option = value-1, value-2, ..., value-n    <- if the option supports a values list\n')
            file_id.write( '#\n')
            file_id.write( '# In section "Trinity parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of Trinity and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --max_reads_per_graph=200000; --min_glue=2\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information that identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the Trinity parameters\n')
            file_id.write( '[Trinity parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_memory = 60', '# suggested maximum memory in GiB to use by Trinity where limiting can be enabled'))
            file_id.write( '{0:<50} {1}\n'.format( 'kmer = 25', '# value or values list of kmer size: maximum, 32.'))
            file_id.write( '{0:<50} {1}\n'.format( 'min_kmer_cov = 1', '# minimum count for Kmers to be assembled by Inchworm'))
            file_id.write( '{0:<50} {1}\n'.format( 'bfly_heap_space_max = 4', '# java maximum heap space setting in GiB'))
            file_id.write( '{0:<50} {1}\n'.format( 'bfly_calculate_cpu = YES', f'# calculate CPUs based on 0.8 of max_memory divided by heap space setting for Butterfly: {get_bfly_calculate_cpu_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'normalized_reads = NO', f'# use normalized reads: {get_normalized_reads_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
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
        error_list.append(f'*** ERROR: The file {get_trinity_config_file()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_trinity_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in Trinity config file.
    '''

    # initialize the control variable
    OK = True

    # get the Trinity option dictionary
    trinity_option_dict = xlib.get_option_dict(get_trinity_config_file())

    # get the experiment identification
    experiment_id = trinity_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the Trinity config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_trinity_name()} config file ...\n')
    (OK, error_list) = check_trinity_config_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file is not valid.\n')
        log.write('Please, correct this file or recreate it.\n')

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

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_trinity_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_trinity_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_trinity_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # for each kmer value, build the process, copy it the cluster and run it
    if OK:

        # get the kmer list
        kmer = trinity_option_dict['Trinity parameters']['kmer']
        kmer_value_list = xlib.split_literal_to_integer_list(kmer)
        
        # for each kmer value, do the tasks
        i = 1
        for kmer_value in kmer_value_list:

            # determine the run directory in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Determining the run directory for kmer {kmer_value} in the cluster ...\n')
            if i > 1:
                current_run_dir = f'{xlib.get_cluster_current_run_dir(experiment_id, xlib.get_trinity_code())}-{i}'
            else:
                current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_trinity_code())
            command = f'mkdir --parents {current_run_dir}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write(f'The directory path is {current_run_dir}.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            i += 1

            # build the Trinity process script
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Building the process script {get_trinity_process_script()} ...\n')
            (OK, error_list) = build_trinity_process_script(cluster_name, current_run_dir, kmer_value)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')
                break

            # upload the process script to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the process script {get_trinity_process_script()} to the directory {current_run_dir} ...\n')
            cluster_path = f'{current_run_dir}/{os.path.basename(get_trinity_process_script())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_trinity_process_script(), cluster_path)
            if OK:
                log.write('The file id uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process script in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_trinity_process_script())} ...\n')
            command = f'chmod u+x {current_run_dir}/{os.path.basename(get_trinity_process_script())}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # build the process starter
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Building the process starter {get_trinity_process_starter()} ...\n')
            (OK, error_list) = build_trinity_process_starter(current_run_dir)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('***ERROR: The file could not be built.\n')
                break

            # upload the process starter to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the process starter {get_trinity_process_starter()} to the directory {current_run_dir} ...\n')
            cluster_path = f'{current_run_dir}/{os.path.basename(get_trinity_process_starter())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_trinity_process_starter(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process starter in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_trinity_process_starter())} ...\n')
            command = f'chmod u+x {current_run_dir}/{os.path.basename(get_trinity_process_starter())}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # submit the process
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_trinity_process_starter())} ...\n')
            OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_trinity_process_starter()), log)

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

def check_trinity_config_file(strict):
    '''
    Check the Trinity config file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        trinity_option_dict = xlib.get_option_dict(get_trinity_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in trinity_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = trinity_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = trinity_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "Trinity parameters"
        if 'Trinity parameters' not in sections_list:
            error_list.append('*** ERROR: the section "Trinity parameters" is not found.')
            OK = False
        else:

            # check section "Trinity parameters" - key "ncpu"
            ncpu = trinity_option_dict.get('Trinity parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "Trinity parameters" - key "max_memory"
            max_memory = trinity_option_dict.get('Trinity parameters', {}).get('max_memory', not_found)
            if max_memory == not_found:
                error_list.append('*** ERROR: the key "max_memory" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_int(max_memory, minimum=1):
                error_list.append('*** ERROR: the key "max_memory" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "Trinity parameters" - key "kmer"
            kmer = trinity_option_dict.get('Trinity parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section Trinity parameters".')
                OK = False
            else:
                kmer_list = xlib.split_literal_to_integer_list(kmer)
                if kmer_list == []:
                    error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 1 and 32.')
                    OK = False
                else:
                    for kmer_item in kmer_list:
                        if not xlib.check_int(kmer_item, minimum=1, maximum=32):
                            error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 1 and 32.')
                            OK = False
                            break

            # check section "Trinity parameters" - key "min_kmer_cov"
            min_kmer_cov = trinity_option_dict.get('Trinity parameters', {}).get('min_kmer_cov', not_found)
            if min_kmer_cov == not_found:
                error_list.append('*** ERROR: the key "min_kmer_cov" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_int(min_kmer_cov, minimum=1):
                error_list.append('*** ERROR: the key "min_kmer_cov" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "Trinity parameters" - key "bfly_heap_space_max"
            bfly_heap_space_max = trinity_option_dict.get('Trinity parameters', {}).get('bfly_heap_space_max', not_found)
            if bfly_heap_space_max == not_found:
                error_list.append('*** ERROR: the key "bfly_heap_space_max" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_int(bfly_heap_space_max, minimum=1):
                error_list.append('*** ERROR: the key "min_kmer_cov" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "Trinity parameters" - key "bfly_calculate_cpu"
            bfly_calculate_cpu = trinity_option_dict.get('Trinity parameters', {}).get('bfly_calculate_cpu', not_found)
            if bfly_calculate_cpu == not_found:
                error_list.append('*** ERROR: the key "bfly_calculate_cpu" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_code(bfly_calculate_cpu, get_bfly_calculate_cpu_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "bfly_calculate_cpu" has to be {get_bfly_calculate_cpu_code_list_text()}.')
                OK = False

            # check section "Trinity parameters" - key "normalized_reads"
            normalized_reads = trinity_option_dict.get('Trinity parameters', {}).get('normalized_reads', not_found)
            if normalized_reads == not_found:
                error_list.append('*** ERROR: the key "normalized_reads" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_code(normalized_reads, get_normalized_reads_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "normalized_reads" has to be {get_normalized_reads_code_list_text()}.')
                OK = False

            # check section "Trinity parameters" - key "other_parameters"
            not_allowed_parameters_list = ['no_version_check', 'seqType', 'left', 'right', 'single', 'CPU', 'max_memory', 'KMER_SIZE', 'bflyHeapSpaceMax', 'bflyCalculateCPU', 'no_normalize_reads', 'output', 'full_cleanup', 'genome_guided_bam']
            other_parameters = trinity_option_dict.get('Trinity parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "Trinity parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "format"
            format = trinity_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = trinity_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'Trinity parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = trinity_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = trinity_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the Trinity config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_trinity_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_trinity_process_script(cluster_name, current_run_dir, kmer_value):
    '''
    Build the current Trinity process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    trinity_option_dict = xlib.get_option_dict(get_trinity_config_file())

    # get the options
    experiment_id = trinity_option_dict['identification']['experiment_id']
    read_dataset_id = trinity_option_dict['identification']['read_dataset_id']
    ncpu = trinity_option_dict['Trinity parameters']['ncpu']
    max_memory = trinity_option_dict['Trinity parameters']['max_memory']
    bfly_heap_space_max = trinity_option_dict['Trinity parameters']['bfly_heap_space_max']
    bfly_calculate_cpu = trinity_option_dict['Trinity parameters']['bfly_calculate_cpu']
    normalized_reads = trinity_option_dict['Trinity parameters']['normalized_reads']
    other_parameters = trinity_option_dict['Trinity parameters']['other_parameters']
    format = 'fq' if trinity_option_dict['library']['format'].upper() == 'FASTQ' else 'fa'
    read_type = trinity_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in trinity_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build library files
    files1 = ''
    files2 = ''
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = trinity_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            files1 += read_file_1 + ','
            if read_type.upper() == 'PE':
                read_file_2 = trinity_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                files2 += read_file_2 + ','
    files1 = files1[:len(files1) - 1]
    if read_type.upper() == 'PE':
        files2 = files2[:len(files2) - 1]

    # write the Trinity process script
    try:
        if not os.path.exists(os.path.dirname(get_trinity_process_script())):
            os.makedirs(os.path.dirname(get_trinity_process_script()))
        with open(get_trinity_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'ulimit -s unlimited\n')
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
            script_file_id.write( 'function run_trinity_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_trinity_anaconda_code()}\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    Trinity --no_version_check --version\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '        Trinity \\\n')
            script_file_id.write( '            --no_version_check \\\n')
            script_file_id.write(f'            --CPU {ncpu} \\\n')
            # -- script_file_id.write(f'            --KMER_SIZE {kmer_value} \\\n')
            script_file_id.write(f'            --seqType {format} \\\n')
            if read_type.upper() == 'PE':
                script_file_id.write(f'            --left {files1} \\\n')
                script_file_id.write(f'            --right {files2} \\\n')
            else:
                script_file_id.write(f'            --single {files1} \\\n')
            script_file_id.write(f'            --max_memory {max_memory}G \\\n')
            script_file_id.write(f'            --bflyHeapSpaceMax {bfly_heap_space_max}G \\\n')
            if bfly_calculate_cpu.upper() == 'YES':
                script_file_id.write( '            --bflyCalculateCPU \\\n')
            if normalized_reads.upper() == 'NO':
                script_file_id.write( '            --no_normalize_reads \\\n')
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for parameter in parameter_list:
                    if parameter.find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write(f'            --{parameter_name} {parameter_value} \\\n')
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        script_file_id.write(f'            --{parameter_name} \\\n')
            script_file_id.write(f'            --output {current_run_dir}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error Trinity $RC; fi\n')
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
            process_name = f'{xlib.get_trinity_name()} process'
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
            script_file_id.write( 'run_trinity_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_trinity_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_trinity_process_starter(current_run_dir):
    '''
    Build the starter of the current Trinity process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Trinity process starter
    try:
        if not os.path.exists(os.path.dirname(get_trinity_process_starter())):
            os.makedirs(os.path.dirname(get_trinity_process_starter()))
        with open(get_trinity_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_trinity_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_trinity_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_trinity_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a Trinity process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_trinity_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_trinity_process_starter()), log)

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

def get_trinity_config_file():
    '''
    Get the Trinity config file path.
    '''

    # assign the Trinity config file path
    trinity_config_file = f'{xlib.get_config_dir()}/{xlib.get_trinity_code()}-config.txt'

    # return the Trinity config file path
    return trinity_config_file

#-------------------------------------------------------------------------------

def get_trinity_process_script():
    '''
    Get the Trinity process script path in the local computer.
    '''

    # assign the Trinity script path
    trinity_process_script = f'{xlib.get_temp_dir()}/{xlib.get_trinity_code()}-process.sh'

    # return the Trinity script path
    return trinity_process_script

#-------------------------------------------------------------------------------

def get_trinity_process_starter():
    '''
    Get the Trinity process starter path in the local computer.
    '''

    # assign the Trinity process starter path
    trinity_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_trinity_code()}-process-starter.sh'

    # return the Trinity starter path
    return trinity_process_starter

#-------------------------------------------------------------------------------

def create_ggtrinity_config_file(experiment_id='exp001', alignment_dataset_id='star-170101-235959'):
    '''
    Create Genome-guided Trinity config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the alignment software
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

    # create the Genome-guided Trinity config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_ggtrinity_config_file())):
            os.makedirs(os.path.dirname(get_ggtrinity_config_file()))
        with open(get_ggtrinity_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and alignment_dataset_id names are fixed in the identification section.\n')
            file_id.write(f'# The alignment file has to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/alignment_dataset_id\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of Trinity in https://github.com/trinityrnaseq/trinityrnaseq/wiki\n')
            file_id.write( '#\n')
            file_id.write( '# In section "Genome-guided Trinity parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of Trinity and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --max_reads_per_graph=200000; --min_glue=2\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information that identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_software = {alignment_software}', f'# alignment software: {get_alignment_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_dataset_id = {alignment_dataset_id}', '# alignment dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the Genome-guided Trinity parameters\n')
            file_id.write( '[Genome-guided Trinity parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_memory = 10', '# suggested maximum memory in GiB to use by Trinity where limiting can be enabled'))
            file_id.write( '{0:<50} {1}\n'.format( 'genome_guided_max_intron = 10000', '# maximum allowed intron length'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ggtrinity_config_file()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_ggtrinity_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in Genome-guided Trinity config file.
    '''

    # initialize the control variable
    OK = True

    # get the Genome-guided Trinity option dictionary
    ggtrinity_option_dict = xlib.get_option_dict(get_ggtrinity_config_file())

    # get the experiment identification
    experiment_id = ggtrinity_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the Genome-guided Trinity config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_ggtrinity_name()} config file ...\n')
    (OK, error_list) = check_ggtrinity_config_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file is not valid.\n')
        log.write('Please, correct this file or recreate it.\n')

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

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_trinity_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_ggtrinity_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_ggtrinity_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_ggtrinity_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Genome-guided Trinity process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_ggtrinity_process_script()} ...\n')
        (OK, error_list) = build_ggtrinity_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Genome-guided Trinity process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_ggtrinity_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ggtrinity_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ggtrinity_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Genome-guided Trinity process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ggtrinity_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ggtrinity_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Genome-guided Trinity process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_ggtrinity_process_starter()} ...\n')
        (OK, error_list) = build_ggtrinity_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Genome-guided Trinity process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_ggtrinity_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_ggtrinity_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_ggtrinity_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Genome-guided Trinity process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_ggtrinity_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_ggtrinity_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Genome-guided Trinity process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ggtrinity_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ggtrinity_process_starter()), log)

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

def check_ggtrinity_config_file(strict):
    '''
    Check the Genome-guided Trinity configu file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        ggtrinity_option_dict = xlib.get_option_dict(get_ggtrinity_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in ggtrinity_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = ggtrinity_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "alignment_software"
            alignment_software = ggtrinity_option_dict.get('identification', {}).get('alignment_software', not_found)
            if alignment_software == not_found:
                error_list.append(f'*** ERROR: the key "alignment_software" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_code(alignment_software, get_alignment_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "alignment_software" has to be {get_alignment_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "alignment_dataset_id"
            alignment_dataset_id = ggtrinity_option_dict.get('identification', {}).get('alignment_dataset_id', not_found)
            if alignment_dataset_id == not_found:
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_startswith(alignment_dataset_id, get_alignment_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" has to start with {get_alignment_software_code_list_text()}.')
                OK = False

        # check section "Genome-guided Trinity parameters"
        if 'Genome-guided Trinity parameters' not in sections_list:
            error_list.append('*** ERROR: the section "Genome-guided Trinity parameters" is not found.')
            OK = False
        else:

            # check section "Genome-guided Trinity parameters" - key "ncpu"
            ncpu = ggtrinity_option_dict.get('Genome-guided Trinity parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "Genome-guided Trinity parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')

            # check section "Genome-guided Trinity parameters" - key "max_memory"
            max_memory = ggtrinity_option_dict.get('Genome-guided Trinity parameters', {}).get('max_memory', not_found)
            if max_memory == not_found:
                error_list.append('*** ERROR: the key "max_memory" is not found in the section "Genome-guided Trinity parameters".')
                OK = False
            elif not xlib.check_int(max_memory, minimum=1):
                error_list.append('*** ERROR: the key "max_memory" has to be an integer number greater than or equal to 1.')

            # check section "Genome-guided Trinity parameters" - key "genome_guided_max_intron"
            genome_guided_max_intron = ggtrinity_option_dict.get('Genome-guided Trinity parameters', {}).get('genome_guided_max_intron', not_found)
            if genome_guided_max_intron == not_found:
                error_list.append('*** ERROR: the key "genome_guided_max_intron" is not found in the section "Genome-guided Trinity parameters".')
                OK = False
            elif not xlib.check_int(genome_guided_max_intron, minimum=1):
                error_list.append('*** ERROR: the key "genome_guided_max_intron" has to be an integer number greater than or equal to 1.')

            # check section "Genome-guided Trinity parameters" - key "other_parameters"
            not_allowed_parameters_list = ['no_version_check', 'CPU', 'max_memory', 'genome_guided_bam', 'genome_guided_max_intron', 'output', 'full_cleanup']
            other_parameters = ggtrinity_option_dict.get('Genome-guided Trinity parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "Genome-guided Trinity parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the Genome-guided Trinity config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_ggtrinity_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ggtrinity_process_script(cluster_name, current_run_dir):
    '''
    Build the current Genome-guided Trinity process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    ggtrinity_option_dict = xlib.get_option_dict(get_ggtrinity_config_file())

    # get the options
    experiment_id = ggtrinity_option_dict['identification']['experiment_id']
    alignment_software = ggtrinity_option_dict['identification']['alignment_software']
    alignment_dataset_id = ggtrinity_option_dict['identification']['alignment_dataset_id']
    ncpu = ggtrinity_option_dict['Genome-guided Trinity parameters']['ncpu']
    max_memory = ggtrinity_option_dict['Genome-guided Trinity parameters']['max_memory']
    genome_guided_max_intron = ggtrinity_option_dict['Genome-guided Trinity parameters']['genome_guided_max_intron']
    other_parameters = ggtrinity_option_dict['Genome-guided Trinity parameters']['other_parameters']

    # set the alignment file paths
    if alignment_software == xlib.get_bowtie2_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/alignment.sam'
        bam_files = '$BAM_DIR/alignment.bam'
    elif alignment_software == xlib.get_gsnap_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/*-split.concordant_uniq'
        bam_files = '$BAM_DIR/*-split.concordant_uniq.bam'
    elif alignment_software == xlib.get_hisat2_code():
        sam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/alignment.sam'
        bam_files = '$BAM_DIR/alignment.bam'
    elif alignment_software == xlib.get_star_code():
        bam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/*-Aligned.sortedByCoord.out.bam'
    elif alignment_software == xlib.get_tophat_code():
        bam_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/accepted_hits.bam'

    # set the merged bam file path
    if alignment_software == xlib.get_bowtie2_code():
        merged_bam_file = bam_files
    elif alignment_software == xlib.get_gsnap_code():
        merged_bam_file = '$BAM_DIR/merged_bam_file.bam'
    elif alignment_software == xlib.get_hisat2_code():
        merged_bam_file = bam_files
    elif alignment_software == xlib.get_star_code():
        merged_bam_file = '$BAM_DIR/merged_bam_file.bam'
    elif alignment_software == xlib.get_tophat_code():
        merged_bam_file = bam_files
 
    # write the Genome-guided Trinity process script
    try:
        if not os.path.exists(os.path.dirname(get_ggtrinity_process_script())):
            os.makedirs(os.path.dirname(get_ggtrinity_process_script()))
        with open(get_ggtrinity_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'ulimit -s unlimited\n')
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
            script_file_id.write(f'BAM_DIR={current_run_dir}/BAM\n')
            script_file_id.write( 'if [ ! -d "$BAM_DIR" ]; then mkdir --parents $BAM_DIR; fi\n')
            script_file_id.write(f'MERGED_BAM_FILE={merged_bam_file}\n')
            script_file_id.write(f'SORTED_BAM_FILE=$BAM_DIR/sorted_bam_file.bam\n')
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
            if alignment_software in [xlib.get_gsnap_code(), xlib.get_hisat2_code()]:
                script_file_id.write( 'function convert_sam2bam\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/convert_sam2bam.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Converting SAM files to BAM format ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write(f'        source activate {xlib.get_samtools_anaconda_code()}\n')
                script_file_id.write(f'        ls {sam_files} > sam-files.txt\n')
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
            if alignment_software in [xlib.get_gsnap_code(), xlib.get_star_code()]:
                script_file_id.write( 'function merge_bam_files\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    STEP_STATUS=$STATUS_DIR/merge_bam_files.ok\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Merging BAM files ..."\n')
                script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
                script_file_id.write( '        echo "This step was previously run."\n')
                script_file_id.write( '    else\n')
                script_file_id.write(f'        FILENUM=`ls -l {bam_files} | grep -v ^d | wc -l`\n')
                script_file_id.write( '        echo "File number: $FILENUM"\n')
                script_file_id.write( '        if [ "$FILENUM" -gt 1 ]; then\n')
                script_file_id.write(f'            ls {bam_files} > bam-files.txt\n')
                script_file_id.write( '            BAM_FILE_LIST=""\n')
                script_file_id.write( '            while read BAM_FILE; do\n')
                script_file_id.write( '                BAM_FILE_LIST=`echo "$BAM_FILE_LIST $BAM_FILE"`\n')
                script_file_id.write( '            done < bam-files.txt\n')
                script_file_id.write(f'            source activate {xlib.get_samtools_anaconda_code()}\n')
                script_file_id.write( '            samtools merge $MERGED_BAM_FILE $BAM_FILE_LIST\n')
                script_file_id.write( '            RC=$?\n')
                script_file_id.write( '            if [ $RC -ne 0 ]; then manage_error samtools-merge $RC; fi\n')
                script_file_id.write( '            conda deactivate\n')
                script_file_id.write( '            echo "BAM files are merged."\n')
                script_file_id.write( '        else\n')
                script_file_id.write(f'           MERGED_BAM_FILE=`ls {bam_files}`\n')
                script_file_id.write( '           echo "There only one BAM file. The merger is not done."\n')
                script_file_id.write( '        fi\n')
                script_file_id.write( '        touch $STEP_STATUS\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function sort_merged_bam_file\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/sort_merged_bam_file.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Sorting the [merged] BAM file ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_samtools_anaconda_code()}\n')
            script_file_id.write( '        samtools sort -o $SORTED_BAM_FILE -O bam $MERGED_BAM_FILE\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error samtools-merger $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "BAM file is sorted."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_ggtrinity_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_ggtrinity_process.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Assembling reads from BAM file ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_trinity_anaconda_code()}\n')
            script_file_id.write( '        Trinity --no_version_check --version\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '            Trinity \\\n')
            script_file_id.write( '                --no_version_check \\\n')
            script_file_id.write(f'                --CPU {ncpu} \\\n')
            script_file_id.write(f'                --max_memory {max_memory}G \\\n')
            script_file_id.write( '                --genome_guided_bam $SORTED_BAM_FILE \\\n')
            script_file_id.write(f'                --genome_guided_max_intron {genome_guided_max_intron} \\\n')
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for parameter in parameter_list:
                    if parameter.find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write(f'            --{parameter_name} {parameter_value} \\\n')
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        script_file_id.write(f'            --{parameter_name} \\\n')
            script_file_id.write(f'                --output {current_run_dir}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error Trinity $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Reads are assembled."\n')
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
            process_name = f'{xlib.get_ggtrinity_name()} process'
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
            if alignment_software in [xlib.get_gsnap_code(), xlib.get_hisat2_code()]:
                script_file_id.write( 'convert_sam2bam\n')
            if alignment_software in [xlib.get_gsnap_code(), xlib.get_star_code()]:
                script_file_id.write( 'merge_bam_files\n')
            script_file_id.write( 'sort_merged_bam_file\n')
            script_file_id.write( 'run_ggtrinity_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ggtrinity_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ggtrinity_process_starter(current_run_dir):
    '''
    Build the starter of the current Genome-guided Trinity process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the Genome-guided Trinity process starter
    try:
        if not os.path.exists(os.path.dirname(get_ggtrinity_process_starter())):
            os.makedirs(os.path.dirname(get_ggtrinity_process_starter()))
        with open(get_ggtrinity_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_ggtrinity_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_ggtrinity_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_ggtrinity_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a Genome-guided Trinity process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_ggtrinity_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ggtrinity_process_starter()), log)

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

def get_ggtrinity_config_file():
    '''
    Get the Genome-guided Trinity config file path.
    '''

    # assign the Genome-guided Trinity config file path
    ggtrinity_config_file = f'{xlib.get_config_dir()}/{xlib.get_ggtrinity_code()}-config.txt'

    # return the Genome-guided Trinity config file path
    return ggtrinity_config_file

#-------------------------------------------------------------------------------

def get_ggtrinity_process_script():
    '''
    Get the Genome-guided Trinity process script path in the local computer.
    '''

    # assign the Genome-guided Trinity script path
    ggtrinity_process_script = f'{xlib.get_temp_dir()}/{xlib.get_ggtrinity_code()}-process.sh'

    # return the Genome-guided Trinity script path
    return ggtrinity_process_script

#-------------------------------------------------------------------------------

def get_ggtrinity_process_starter():
    '''
    Get the Genome-guided Trinity process starter path in the local computer.
    '''

    # assign the Genome-guided Trinity process starter path
    ggtrinity_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_ggtrinity_code()}-process-starter.sh'

    # return the Genome-guided Trinity starter path
    return ggtrinity_process_starter

#-------------------------------------------------------------------------------

def create_insilico_read_normalization_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create insilico_read_normalization config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the Trinity config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_insilico_read_normalization_config_file())):
            os.makedirs(os.path.dirname(get_insilico_read_normalization_config_file()))
        with open(get_insilico_read_normalization_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of insilico_read_normalization (Trinity package) and their meaning in "https://github.com/trinityrnaseq/trinityrnaseq/wiki".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "insilico_read_normalization parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of insilico_read_normalization and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --pairs_together; --PARALLEL_STATS\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information that identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the insilico_read_normalization parameters\n')
            file_id.write( '[insilico_read_normalization parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'kmer = 25', '# K-MER size'))
            file_id.write( '{0:<50} {1}\n'.format( 'ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'jm = 10', '# maximum memory in GiB to use for k-mer counting by jellyfish'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_cov = 30', '# targeted maximum coverage for reads'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
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
        error_list.append(f'*** ERROR: The file {get_insilico_read_normalization_config_file()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_insilico_read_normalization_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in insilico_read_normalization config file.
    '''

    # initialize the control variable
    OK = True

    # get the insilico_read_normalization option dictionary
    insilico_read_normalization_option_dict = xlib.get_option_dict(get_insilico_read_normalization_config_file())

    # get the experiment identification
    experiment_id = insilico_read_normalization_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the insilico_read_normalization config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_insilico_read_normalization_name()} config file ...\n')
    (OK, error_list) = check_insilico_read_normalization_config_file(strict=True)
    if OK:
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The config file is not valid.\n')
        log.write('Please, correct this file or recreate it.\n')

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

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_trinity_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_trinity_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_trinity_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_insilico_read_normalization_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the insilico_read_normalization process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_insilico_read_normalization_process_script()} ...\n')
        (OK, error_list) = build_insilico_read_normalization_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the insilico_read_normalization process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_insilico_read_normalization_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_insilico_read_normalization_process_script(), cluster_path)
        if OK:
            log.write('The file id uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the insilico_read_normalization process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set on.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the insilico_read_normalization process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_insilico_read_normalization_process_starter()} ...\n')
        (OK, error_list) = build_insilico_read_normalization_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the insilico_read_normalization process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_insilico_read_normalization_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_insilico_read_normalization_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the insilico_read_normalization process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set on.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the insilico_read_normalization process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter()), log)

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

def check_insilico_read_normalization_config_file(strict):
    '''
    Check the insilico_read_normalization configu file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        insilico_read_normalization_option_dict = xlib.get_option_dict(get_insilico_read_normalization_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in insilico_read_normalization_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = insilico_read_normalization_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = insilico_read_normalization_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "insilico_read_normalization parameters"
        if 'insilico_read_normalization parameters' not in sections_list:
            error_list.append('*** ERROR: the section "insilico_read_normalization parameters" is not found.')
            OK = False
        else:

            # check section "insilico_read_normalization parameters" - key "kmer"
            kmer = insilico_read_normalization_option_dict.get('insilico_read_normalization parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section "insilico_read_normalization parameters".')
                OK = False
            elif not xlib.check_int(kmer, minimum=1):
                error_list.append('*** ERROR: the key "kmer" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "insilico_read_normalization parameters" - key "ncpu"
            ncpu = insilico_read_normalization_option_dict.get('insilico_read_normalization parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "insilico_read_normalization parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "insilico_read_normalization parameters" - key "jm"
            jm = insilico_read_normalization_option_dict.get('insilico_read_normalization parameters', {}).get('jm', not_found)
            if jm == not_found:
                error_list.append('*** ERROR: the key "jm" is not found in the section "insilico_read_normalization parameters".')
                OK = False
            elif not xlib.check_int(jm, minimum=1):
                error_list.append('*** ERROR: the key "jm" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "insilico_read_normalization parameters" - key "max_cov"
            max_cov = insilico_read_normalization_option_dict.get('insilico_read_normalization parameters', {}).get('max_cov', not_found)
            if max_cov == not_found:
                error_list.append('*** ERROR: the key "max_cov" is not found in the section "insilico_read_normalization parameters".')
                OK = False
            elif not xlib.check_int(max_cov, minimum=1):
                error_list.append('*** ERROR: the key "max_cov" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "insilico_read_normalization parameters" - key "other_parameters"
            not_allowed_parameters_list = ['seqType', 'left', 'right', 'left_list', 'right_list', 'single', 'single_list', 'KMER_SIZE', 'CPU', 'JM', 'max_cov', 'no_cleanup', 'output']
            other_parameters = insilico_read_normalization_option_dict.get('insilico_read_normalization parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "insilico_read_normalization parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "format"
            format = insilico_read_normalization_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = insilico_read_normalization_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'insilico_read_normalization parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = insilico_read_normalization_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = insilico_read_normalization_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the insilico_read_normalization config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_insilico_read_normalization_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_insilico_read_normalization_process_script(cluster_name, current_run_dir):
    '''
    Build the current insilico_read_normalization process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    insilico_read_normalization_option_dict = xlib.get_option_dict(get_insilico_read_normalization_config_file())

    # get the options
    experiment_id = insilico_read_normalization_option_dict['identification']['experiment_id']
    read_dataset_id = insilico_read_normalization_option_dict['identification']['read_dataset_id']
    kmer = insilico_read_normalization_option_dict['insilico_read_normalization parameters']['kmer']
    ncpu = insilico_read_normalization_option_dict['insilico_read_normalization parameters']['ncpu']
    jm = insilico_read_normalization_option_dict['insilico_read_normalization parameters']['jm']
    max_cov = insilico_read_normalization_option_dict['insilico_read_normalization parameters']['max_cov']
    other_parameters = insilico_read_normalization_option_dict['insilico_read_normalization parameters']['other_parameters']
    format = 'fq' if insilico_read_normalization_option_dict['library']['format'].upper() == 'FASTQ' else 'fa'
    read_type = insilico_read_normalization_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in insilico_read_normalization_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build library files
    files1 = ''
    files2 = ''
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = insilico_read_normalization_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            files1 += read_file_1 + ','
            if read_type.upper() == 'PE':
                read_file_2 = insilico_read_normalization_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                files2 += read_file_2 + ','
    files1 = files1[:len(files1) - 1]
    if read_type.upper() == 'PE':
        files2 = files2[:len(files2) - 1]

    # set the output normalised read directory
    normalised_read_dir = xlib.get_cluster_experiment_read_dataset_dir(experiment_id, os.path.basename(current_run_dir))

    # write the insilico_read_normalization process script
    try:
        if not os.path.exists(os.path.dirname(get_insilico_read_normalization_process_script())):
            os.makedirs(os.path.dirname(get_insilico_read_normalization_process_script()))
        with open(get_insilico_read_normalization_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function run_insilico_read_normalization_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_insilico_read_normalization_process.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Normalizing read files dataset ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_trinity_anaconda_code()}\n')
            script_file_id.write( '        Trinity --no_version_check --version\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '            insilico_read_normalization.pl \\\n')
            script_file_id.write(f'                --CPU {ncpu} \\\n')
            script_file_id.write(f'                --KMER_SIZE {kmer} \\\n')
            script_file_id.write(f'                --seqType {format} \\\n')
            if read_type.upper() == 'PE':
                script_file_id.write(f'            --left {files1} \\\n')
                script_file_id.write(f'            --right {files2} \\\n')
            else:
                script_file_id.write(f'            --single {files1} \\\n')
            script_file_id.write(f'                --JM {jm}G \\\n')
            script_file_id.write(f'                --max_cov {max_cov} \\\n')
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
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
            script_file_id.write(f'                --output {current_run_dir}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error insilico_read_normalization.pl $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "Files are normalized."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function move_normalized_reads\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/move_normalized_reads.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Moving normalized files to read dataset ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        mkdir --parents {normalised_read_dir}\n')
            script_file_id.write(f'        rm left.norm.fq right.norm.fq\n')
            script_file_id.write(f'        mv *.fq {normalised_read_dir}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then mv $RC; fi\n')
            script_file_id.write( '        echo "Files are moved."\n')
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
            process_name = f'{xlib.get_insilico_read_normalization_name()} process'
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
            script_file_id.write( 'run_insilico_read_normalization_process\n')
            script_file_id.write( 'move_normalized_reads\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_insilico_read_normalization_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_insilico_read_normalization_process_starter(current_run_dir):
    '''
    Build the starter of the current insilico_read_normalization process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the insilico_read_normalization process starter
    try:
        if not os.path.exists(os.path.dirname(get_insilico_read_normalization_process_starter())):
            os.makedirs(os.path.dirname(get_insilico_read_normalization_process_starter()))
        with open(get_insilico_read_normalization_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_insilico_read_normalization_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_insilico_read_normalization_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a insilico_read_normalization process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_insilico_read_normalization_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter()), log)

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

def get_insilico_read_normalization_config_file():
    '''
    Get the insilico_read_normalization config file path.
    '''

    # assign the insilico_read_normalization config file path
    insilico_read_normalization_config_file = f'{xlib.get_config_dir()}/{xlib.get_insilico_read_normalization_code()}-config.txt'

    # return the insilico_read_normalization config file path
    return insilico_read_normalization_config_file

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_process_script():
    '''
    Get the insilico_read_normalization process script path in the local computer.
    '''

    # assign the insilico_read_normalization script path
    insilico_read_normalization_process_script = f'{xlib.get_temp_dir()}/{xlib.get_insilico_read_normalization_code()}-process.sh'

    # return the insilico_read_normalization script path
    return insilico_read_normalization_process_script

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_process_starter():
    '''
    Get the insilico_read_normalization process starter path in the local computer.
    '''

    # assign the insilico_read_normalization process starter path
    insilico_read_normalization_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_insilico_read_normalization_code()}-process-starter.sh'

    # return the insilico_read_normalization starter path
    return insilico_read_normalization_process_starter

#-------------------------------------------------------------------------------
    
def get_bfly_calculate_cpu_code_list():
    '''
    Get the code list of "bfly_calculate_cpu".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_bfly_calculate_cpu_code_list_text():
    '''
    Get the code list of "bfly_calculate_cpu" as text.
    '''

    return str(get_bfly_calculate_cpu_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_normalized_reads_code_list():
    '''
    Get the code list of "normalized_reads".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_normalized_reads_code_list_text():
    '''
    Get the code list of "normalized_reads" as text.
    '''

    return str(get_normalized_reads_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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

    return ['SE', 'PE']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "read_type" as text.
    '''

    return 'SE (single-end) or PE (pair-end)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the Trinity process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
