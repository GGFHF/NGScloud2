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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of Trinity and their meaning in https://github.com/trinityrnaseq/trinityrnaseq/wiki.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# There are two formats to set an option:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value                             <- if the option supports a single value'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value-1, value-2, ..., value-n    <- if the option supports a values list'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "Trinity parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of Trinity and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --max_reads_per_graph=200000; --min_glue=2'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the Trinity parameters'))
            file_id.write( '{0}\n'.format('[Trinity parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format('max_memory = 60', '# suggested maximum memory in GiB to use by Trinity where limiting can be enabled'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = 25', '# value or values list of kmer size: maximum, 32.'))
            file_id.write( '{0:<50} {1}\n'.format('min_kmer_cov = 1', '# minimum count for Kmers to be assembled by Inchworm'))
            file_id.write( '{0:<50} {1}\n'.format('bfly_heap_space_max = 4', '# java maximum heap space setting in GiB'))
            file_id.write( '{0:<50} {1}\n'.format('bfly_calculate_cpu = YES', '# calculate CPUs based on 0.8 of max_memory divided by heap space setting for Butterfly: {0}'.format(get_bfly_calculate_cpu_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('normalized_reads = NO', '# use normalized reads: {0}'.format(get_normalized_reads_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# read type: {0}'.format(get_read_type_code_list_text())))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first library.'))
                file_id.write( '{0}\n'.format('[library-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('read_file_1 = {0}'.format(os.path.basename(file_1_list[i])), '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = {0}'.format(os.path.basename(file_2_list[i])), '# name of the - strand reads file in PE read type or NONE in SE case'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more libraries, you have to repeat the section library-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_trinity_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_trinity_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_trinity_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_trinity_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_trinity_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # for each kmer value, build the process, copy it the cluster and run it
    if OK:

        # get the kmer list
        kmer = trinity_option_dict['Trinity parameters']['kmer']
        kmer_list = xlib.split_literal_to_integer_list(kmer)
        
        # for each kmer value, do the tasks
        i = 1
        for kmer_value in kmer_list:

            # determine the run directory in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Determining the run directory for kmer {0} in the cluster ...\n'.format(kmer_value))
            if i > 1:
                current_run_dir = '{0}-{1}'.format(xlib.get_cluster_current_run_dir(experiment_id, xlib.get_trinity_code()), i)
            else:
                current_run_dir = '{0}'.format(xlib.get_cluster_current_run_dir(experiment_id, xlib.get_trinity_code()))
            command = f'mkdir --parents {current_run_dir}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The directory path is {0}.\n'.format(current_run_dir))
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            i += 1

            # build the Trinity process script
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process script {0} ...\n'.format(get_trinity_process_script()))
            (OK, error_list) = build_trinity_process_script(cluster_name, current_run_dir, kmer_value)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')
                break

            # upload the process script to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_trinity_process_script(), current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_trinity_process_script()))
            (OK, error_list) = xssh.put_file(sftp_client, get_trinity_process_script(), cluster_path)
            if OK:
                log.write('The file id uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process script in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_trinity_process_script())))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_trinity_process_script()))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # build the process starter
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process starter {0} ...\n'.format(get_trinity_process_starter()))
            (OK, error_list) = build_trinity_process_starter(current_run_dir)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('***ERROR: The file could not be built.\n')
                break

            # upload the process starter to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_trinity_process_starter(), current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_trinity_process_starter()))
            (OK, error_list) = xssh.put_file(sftp_client, get_trinity_process_starter(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process starter in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_trinity_process_starter())))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_trinity_process_starter()))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # submit the process
            log.write(f'{xlib.get_separator()}\n')
            log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_trinity_process_starter())))
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
                error_list.append('*** ERROR: the key "bfly_calculate_cpu" has to be {0}.'.format(get_bfly_calculate_cpu_code_list_text()))
                OK = False

            # check section "Trinity parameters" - key "normalized_reads"
            normalized_reads = trinity_option_dict.get('Trinity parameters', {}).get('normalized_reads', not_found)
            if normalized_reads == not_found:
                error_list.append('*** ERROR: the key "normalized_reads" is not found in the section "Trinity parameters".')
                OK = False
            elif not xlib.check_code(normalized_reads, get_normalized_reads_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "normalized_reads" has to be {0}.'.format(get_normalized_reads_code_list_text()))
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = trinity_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "read_type" has to be {0}.'.format(get_read_type_code_list_text()))
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
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = trinity_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = trinity_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the Trinity config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_trinity_name()))

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
            script_file_id.write( '{0}\n'.format('ulimit -s unlimited'))
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('TRINITY_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_trinity_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$TRINITY_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_trinity_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_trinity_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    Trinity --no_version_check --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        Trinity \\'))
            script_file_id.write( '{0}\n'.format('            --no_version_check \\'))
            script_file_id.write( '{0}\n'.format('            --CPU {0} \\'.format(ncpu)))
            script_file_id.write( '{0}\n'.format('            --KMER_SIZE {0} \\'.format(kmer_value)))
            script_file_id.write( '{0}\n'.format('            --seqType {0} \\'.format(format)))
            if read_type.upper() == 'PE':
                script_file_id.write( '{0}\n'.format('            --left {0} \\'.format(files1)))
                script_file_id.write( '{0}\n'.format('            --right {0} \\'.format(files2)))
            else:
                script_file_id.write( '{0}\n'.format('            --single {0} \\'.format(files1)))
            script_file_id.write( '{0}\n'.format('            --max_memory {0}G \\'.format(max_memory)))
            script_file_id.write( '{0}\n'.format('            --bflyHeapSpaceMax {0}G \\'.format(bfly_heap_space_max)))
            if bfly_calculate_cpu.upper() == 'YES':
                script_file_id.write( '{0}\n'.format('            --bflyCalculateCPU \\'))
            if normalized_reads.upper() == 'NO':
                script_file_id.write( '{0}\n'.format('            --no_normalize_reads \\'))
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for parameter in parameter_list:
                    if parameter.find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} {1} \\'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
            script_file_id.write( '{0}\n'.format('            --output {0}'.format(current_run_dir)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error Trinity $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_trinity_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_trinity_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_trinity_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_trinity_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_trinity_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_trinity_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_trinity_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_trinity_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def determine_trinity_cluster():
    '''
    Determine the cluster to the current Trinity experiment.
    '''

    # initialize the template and cluster names
    template_name = ''
    cluster_name = ''

    # ...

    # return the template and cluster names
    return (template_name, cluster_name)

#-------------------------------------------------------------------------------

def get_trinity_config_file():
    '''
    Get the Trinity config file path.
    '''

    # assign the Trinity config file path
    trinity_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_trinity_code())

    # return the Trinity config file path
    return trinity_config_file

#-------------------------------------------------------------------------------

def get_trinity_process_script():
    '''
    Get the Trinity process script path in the local computer.
    '''

    # assign the Trinity script path
    trinity_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_trinity_code())

    # return the Trinity script path
    return trinity_process_script

#-------------------------------------------------------------------------------

def get_trinity_process_starter():
    '''
    Get the Trinity process starter path in the local computer.
    '''

    # assign the Trinity process starter path
    trinity_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_trinity_code())

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
    if alignment_dataset_id.startswith(xlib.get_star_code()):
        alignment_software = xlib.get_star_code()
    elif alignment_dataset_id.startswith(xlib.get_tophat_code()):
        alignment_software = xlib.get_tophat_code()

    # create the Genome-guided Trinity config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_ggtrinity_config_file())):
            os.makedirs(os.path.dirname(get_ggtrinity_config_file()))
        with open(get_ggtrinity_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and alignment_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('# The alignment file has to be located in the cluster directory {0}/experiment_id/alignment_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of Trinity in https://github.com/trinityrnaseq/trinityrnaseq/wiki.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "Genome-guided Trinity parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of Trinity and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --max_reads_per_graph=200000; --min_glue=2'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('alignment_software = {0}'.format(alignment_software), '# alignment software: {0}'.format(get_alignment_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('alignment_dataset_id = {0}'.format(alignment_dataset_id), '# alignment dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the Genome-guided Trinity parameters'))
            file_id.write( '{0}\n'.format('[Genome-guided Trinity parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format('max_memory = 10', '# suggested maximum memory in GiB to use by Trinity where limiting can be enabled'))
            file_id.write( '{0:<50} {1}\n'.format('genome_guided_max_intron = 10000', '# maximum allowed intron length'))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ggtrinity_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_ggtrinity_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_trinity_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ggtrinity_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_ggtrinity_name()))

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
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Genome-guided Trinity process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_ggtrinity_process_script()))
        (OK, error_list) = build_ggtrinity_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the Genome-guided Trinity process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_ggtrinity_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ggtrinity_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ggtrinity_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Genome-guided Trinity process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ggtrinity_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ggtrinity_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the Genome-guided Trinity process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ggtrinity_process_starter()))
        (OK, error_list) = build_ggtrinity_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the Genome-guided Trinity process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_ggtrinity_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ggtrinity_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ggtrinity_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the Genome-guided Trinity process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ggtrinity_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ggtrinity_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the Genome-guided Trinity process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ggtrinity_process_starter())))
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
                error_list.append('*** ERROR: the key "alignment_software" is not found in the section "{0}".'.format(section))
                OK = False
            elif not xlib.check_code(alignment_software, get_alignment_software_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "alignment_software" has to be {0}.'.format(get_alignment_software_code_list_text()))
                OK = False

            # check section "identification" - key "alignment_dataset_id"
            alignment_dataset_id = ggtrinity_option_dict.get('identification', {}).get('alignment_dataset_id', not_found)
            if alignment_dataset_id == not_found:
                error_list.append('*** ERROR: the key "alignment_dataset_id" is not found in the section "{0}".'.format(section))
                OK = False
            elif not xlib.check_startswith(alignment_dataset_id, get_alignment_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "alignment_dataset_id" has to start with {0}.'.format(get_alignment_software_code_list_text()))
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
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_ggtrinity_name()))

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

    # set the alignment file path
    if alignment_software == xlib.get_star_code():
        alignment_file = '{0}/starAligned.sortedByCoord.out.bam'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id))
    elif alignment_software == xlib.get_tophat_code():
        alignment_file = '{0}/accepted_hits.bam'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id))

    # write the Genome-guided Trinity process script
    try:
        if not os.path.exists(os.path.dirname(get_ggtrinity_process_script())):
            os.makedirs(os.path.dirname(get_ggtrinity_process_script()))
        with open(get_ggtrinity_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('ulimit -s unlimited'))
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('TRINITY_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_trinity_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$TRINITY_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_trinity_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_ggtrinity_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    Trinity --no_version_check --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Assembling from genome-aligned reads ..."'))
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        Trinity \\'))
            script_file_id.write( '{0}\n'.format('            --no_version_check \\'))
            script_file_id.write( '{0}\n'.format('            --CPU {0} \\'.format(ncpu)))
            script_file_id.write( '{0}\n'.format('            --max_memory {0}G \\'.format(max_memory)))
            script_file_id.write( '{0}\n'.format('            --genome_guided_bam {0} \\'.format(alignment_file)))
            script_file_id.write( '{0}\n'.format('            --genome_guided_max_intron {0} \\'.format(genome_guided_max_intron)))
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for parameter in parameter_list:
                    if parameter.find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} {1} \\'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
            script_file_id.write( '{0}\n'.format('            --output {0}'.format(current_run_dir)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error Trinity $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_ggtrinity_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_ggtrinity_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_ggtrinity_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_ggtrinity_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_ggtrinity_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ggtrinity_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_ggtrinity_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ggtrinity_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ggtrinity_config_file():
    '''
    Get the Genome-guided Trinity config file path.
    '''

    # assign the Genome-guided Trinity config file path
    ggtrinity_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_ggtrinity_code())

    # return the Genome-guided Trinity config file path
    return ggtrinity_config_file

#-------------------------------------------------------------------------------

def get_ggtrinity_process_script():
    '''
    Get the Genome-guided Trinity process script path in the local computer.
    '''

    # assign the Genome-guided Trinity script path
    ggtrinity_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_ggtrinity_code())

    # return the Genome-guided Trinity script path
    return ggtrinity_process_script

#-------------------------------------------------------------------------------

def get_ggtrinity_process_starter():
    '''
    Get the Genome-guided Trinity process starter path in the local computer.
    '''

    # assign the Genome-guided Trinity process starter path
    ggtrinity_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ggtrinity_code())

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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of insilico_read_normalization (Trinity package) and their meaning in https://github.com/trinityrnaseq/trinityrnaseq/wiki.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "insilico_read_normalization parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of insilico_read_normalization and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --pairs_together; --PARALLEL_STATS'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the insilico_read_normalization parameters'))
            file_id.write( '{0}\n'.format('[insilico_read_normalization parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = 25', '# K-MER size'))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of CPUs for use'))
            file_id.write( '{0:<50} {1}\n'.format('jm = 10', '# maximum memory in GiB to use for k-mer counting by jellyfish'))
            file_id.write( '{0:<50} {1}\n'.format('max_cov = 30', '# targeted maximum coverage for reads'))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# read type: {0}'.format(get_read_type_code_list_text())))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first library.'))
                file_id.write( '{0}\n'.format('[library-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('read_file_1 = {0}'.format(os.path.basename(file_1_list[i])), '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = {0}'.format(os.path.basename(file_2_list[i])), '# name of the - strand reads file in PE read type or NONE in SE case'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more libraries, you have to repeat the section library-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_insilico_read_normalization_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_insilico_read_normalization_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check Trinity is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_trinity_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_trinity_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_trinity_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = '{0}'.format(xlib.get_cluster_current_run_dir(experiment_id, xlib.get_insilico_read_normalization_code()))
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the insilico_read_normalization process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_insilico_read_normalization_process_script()))
        (OK, error_list) = build_insilico_read_normalization_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the insilico_read_normalization process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_insilico_read_normalization_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_insilico_read_normalization_process_script(), cluster_path)
        if OK:
            log.write('The file id uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the insilico_read_normalization process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set on.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the insilico_read_normalization process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_insilico_read_normalization_process_starter()))
        (OK, error_list) = build_insilico_read_normalization_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the insilico_read_normalization process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_insilico_read_normalization_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_insilico_read_normalization_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the insilico_read_normalization process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set on.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the insilico_read_normalization process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_starter())))
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = insilico_read_normalization_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "read_type" has to be {0}.'.format(get_read_type_code_list_text()))
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
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = insilico_read_normalization_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = insilico_read_normalization_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the insilico_read_normalization config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_insilico_read_normalization_name()))

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
    read_file_1_list = []
    read_file_2_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = insilico_read_normalization_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            read_file_1_list.append(read_file_1)
            if read_type.upper() == 'PE':
                read_file_2 = insilico_read_normalization_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                read_file_2_list.append(read_file_2)

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
            script_file_id.write( '{0}\n'.format('TRINITY_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_trinity_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$TRINITY_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_trinity_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_insilico_read_normalization_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    Trinity --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        insilico_read_normalization.pl \\'))
            script_file_id.write( '{0}\n'.format('            --seqType {0} \\'.format(format)))
            if read_type.upper() == 'PE':
                if len(read_file_1_list) == 1:
                    script_file_id.write( '{0}\n'.format('            --left {0} \\'.format(read_file_1_list[0])))
                    script_file_id.write( '{0}\n'.format('            --right {0} \\'.format(read_file_2_list[0])))
                else:
                    for i in range(len(read_file_1_list)):
                        script_file_id.write( '{0}\n'.format('            --left_list {0} \\'.format(read_file_1_list[i])))
                        script_file_id.write( '{0}\n'.format('            --right_list {0} \\'.format(read_file_2_list[i])))
            else:
                if len(read_file_1_list) == 1:
                    script_file_id.write( '{0}\n'.format('            --single {0} \\'.format(read_file_1_list[0])))
                else:
                    for i in range(len(read_file_1_list)):
                        script_file_id.write( '{0}\n'.format('            --single_list {0} \\'.format(read_file_1_list[0])))
            script_file_id.write( '{0}\n'.format('            --KMER_SIZE {0} \\'.format(kmer)))
            script_file_id.write( '{0}\n'.format('            --CPU {0} \\'.format(ncpu)))
            script_file_id.write( '{0}\n'.format('            --JM {0}G \\'.format(jm)))
            script_file_id.write( '{0}\n'.format('            --max_cov {0} \\'.format(max_cov)))
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for parameter in parameter_list:
                    if parameter.find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} {1} \\'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter)
                        parameter_name = mo.group(1).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
            script_file_id.write( '{0}\n'.format('            --output {0}'.format(current_run_dir)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error insilico_read_normalization.pl $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_insilico_read_normalization_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_insilico_read_normalization_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_insilico_read_normalization_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_insilico_read_normalization_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_insilico_read_normalization_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_insilico_read_normalization_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_insilico_read_normalization_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_insilico_read_normalization_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_config_file():
    '''
    Get the insilico_read_normalization config file path.
    '''

    # assign the insilico_read_normalization config file path
    insilico_read_normalization_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_insilico_read_normalization_code())

    # return the insilico_read_normalization config file path
    return insilico_read_normalization_config_file

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_process_script():
    '''
    Get the insilico_read_normalization process script path in the local computer.
    '''

    # assign the insilico_read_normalization script path
    insilico_read_normalization_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_insilico_read_normalization_code())

    # return the insilico_read_normalization script path
    return insilico_read_normalization_process_script

#-------------------------------------------------------------------------------

def get_insilico_read_normalization_process_starter():
    '''
    Get the insilico_read_normalization process starter path in the local computer.
    '''

    # assign the insilico_read_normalization process starter path
    insilico_read_normalization_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_insilico_read_normalization_code())

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

    return [xlib.get_star_code(), xlib.get_tophat_code()]

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list_text():
    '''
    Get the code list of "alignment_software" as text.
    '''

    return '{0} ({1}) or {2} ({3})'.format(xlib.get_star_code(), xlib.get_star_name(), xlib.get_tophat_code(), xlib.get_tophat_name())

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
