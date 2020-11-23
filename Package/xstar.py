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
This file contains functions related to the STAR process used in both console
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

def create_star_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', gtf_file='Arabidopsis_thaliana.TAIR10.36.gtf', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create STAR config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the STAR config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_star_config_file())):
            os.makedirs(os.path.dirname(get_star_config_file()))
        with open(get_star_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference and GTF files have to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id, reference_dataset_id, reference_file, gtf_file and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of STAR and their meaning in "https://github.com/alexdobin/STAR".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "STAR parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of STAR and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --outSAMattributes=All; --limitGenomeGenerateRAM=48000000000\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format(f'gtf_file = {gtf_file}', '# GTF file name'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the STAR parameters\n')
            file_id.write( '[STAR parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'index_building = YES', f'# index building when a reference is used: {get_index_building_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'limit_genome_generate_ram = 31000000000', '# maximum available RAM (in bytes) for genome index generation'))
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'two_pass_mode = NONE', f'# 2-pass mapping mode: {get_two_pass_mode_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'two_pass_1_readsn = -1', '# number of reads to process for the 1st step; use -1 to map all reads in the first step'))
            file_id.write( '{0:<50} {1}\n'.format( 'out_filter_multimap_nmax = 20', '# maximun number of multiple alignments allowed for a read'))
            file_id.write( '{0:<50} {1}\n'.format( 'out_reads_unmapped = FASTX', f'# output of unmapped and partially mapped reads in separate files: {get_out_reads_unmapped_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'cleanup = NO', f'# leave intermediate files: {get_cleanup_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the information of the library (only one library is allowed)\n')
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
        error_list.append(f'*** ERROR: The file {get_star_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_star_process(cluster_name, log, function=None):
    '''
    Run a STAR process.
    '''

    # initialize the control variable
    OK = True

    # get the STAR option dictionary
    star_option_dict = xlib.get_option_dict(get_star_config_file())

    # get the experiment identification
    experiment_id = star_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the STAR config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_star_name()} config file ...\n')
    (OK, error_list) = check_star_config_file(strict=True)
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

    # check STAR is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_star_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_star_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_star_name()} installation could not be performed.\n')

    # check SAMtools is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_samtools_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_samtools_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_samtools_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_star_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the STAR process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_star_process_script()} ...\n')
        (OK, error_list) = build_star_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the STAR process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_star_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_star_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_star_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the STAR process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_star_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_star_process_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the STAR process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_star_process_starter()} ...\n')
        (OK, error_list) = build_star_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the STAR process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_star_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_star_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_star_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the STAR process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_star_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_star_process_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the STAR process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_star_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_star_process_starter()), log)

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

def check_star_config_file(strict):
    '''
    Check the STAR config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        star_option_dict = xlib.get_option_dict(get_star_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in star_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = star_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = star_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = star_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "gtf_file"
            gtf_file = star_option_dict.get('identification', {}).get('gtf_file', not_found)
            if gtf_file == not_found:
                error_list.append('*** ERROR: the key "gtf_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = star_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "STAR parameters"
        if 'STAR parameters' not in sections_list:
            error_list.append('*** ERROR: the section "STAR parameters" is not found.')
            OK = False
        else:

            # check section "STAR parameters" - key "index_building"
            index_building = star_option_dict.get('STAR parameters', {}).get('index_building', not_found)
            if index_building == not_found:
                error_list.append('*** ERROR: the key "index_building" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(index_building, get_index_building_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "index_building" has to be {get_index_building_code_list_text()}.')
                OK = False

            # check section "STAR parameters" - key "limit_genome_generate_ram"
            limit_genome_generate_ram = star_option_dict.get('STAR parameters', {}).get('limit_genome_generate_ram', not_found)
            if limit_genome_generate_ram == not_found:
                error_list.append('*** ERROR: the key "limit_genome_generate_ram" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_int(limit_genome_generate_ram, minimum=1):
                error_list.append('*** ERROR: the key "limit_genome_generate_ram" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "STAR parameters" - key "threads"
            threads = star_option_dict.get('STAR parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "STAR parameters" - key "out_filter_multimap_nmax"
            out_filter_multimap_nmax = star_option_dict.get('STAR parameters', {}).get('out_filter_multimap_nmax', not_found)
            if out_filter_multimap_nmax == not_found:
                error_list.append('*** ERROR: the key "out_filter_multimap_nmax" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_int(out_filter_multimap_nmax, minimum=1):
                error_list.append('*** ERROR: the key "out_filter_multimap_nmax" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "STAR parameters" - key "two_pass_mode"
            two_pass_mode = star_option_dict.get('STAR parameters', {}).get('two_pass_mode', not_found)
            if two_pass_mode == not_found:
                error_list.append('*** ERROR: the key "two_pass_mode" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(two_pass_mode, get_two_pass_mode_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "two_pass_mode" has to be {get_two_pass_mode_code_list_text()}.')
                OK = False

            # check section "STAR parameters" - key "two_pass_1_readsn"
            two_pass_1_readsn = star_option_dict.get('STAR parameters', {}).get('two_pass_1_readsn', not_found)
            if two_pass_1_readsn == not_found:
                error_list.append('*** ERROR: the key "two_pass_1_readsn" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_int(two_pass_1_readsn, minimum=1) and not xlib.check_int(two_pass_1_readsn, minimum=-1, maximum=-1):
                error_list.append('*** ERROR: the key "two_pass_1_readsn" has to be an integer number greater than or equal to 1, or -1 (all redas).')
                OK = False

            # check section "STAR parameters" - key "out_reads_unmapped"
            out_reads_unmapped = star_option_dict.get('STAR parameters', {}).get('out_reads_unmapped', not_found)
            if out_reads_unmapped == not_found:
                error_list.append('*** ERROR: the key "out_reads_unmapped" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(out_reads_unmapped, get_out_reads_unmapped_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "out_reads_unmapped" has to be {get_out_reads_unmapped_code_list_text()}.')
                OK = False

            # check section "STAR parameters" - key "cleanup"
            cleanup = star_option_dict.get('STAR parameters', {}).get('cleanup', not_found)
            if cleanup == not_found:
                error_list.append('*** ERROR: the key "cleanup" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(cleanup, get_cleanup_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "cleanup" has to be {get_cleanup_code_list_text()}.')
                OK = False

            # check section "STAR parameters" - key "other_parameters"
            not_allowed_parameters_list = ['runMode', 'runThreadN', 'limitGenomeGenerateRAM', 'genomeDir', 'readFilesCommand', 'readFilesIn', 'sjdbGTFfile', 'twopassMode', 'twopass1readsN', 'quantMode', 'outFilterMultimapNmax', 'outSAMunmapped', 'outFileNamePrefix', 'outTmpKeep', 'outSAMtype', 'outSAMattributes']
            other_parameters = star_option_dict.get('STAR parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "STAR parameters".')
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
            format = star_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = star_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'STAR parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = star_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = star_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_star_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_star_process_script(cluster_name, current_run_dir):
    '''
    Build the current STAR process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the STAR option dictionary
    star_option_dict = xlib.get_option_dict(get_star_config_file())

    # get the options
    experiment_id = star_option_dict['identification']['experiment_id']
    reference_dataset_id = star_option_dict['identification']['reference_dataset_id']
    reference_file = star_option_dict['identification']['reference_file']
    gtf_file = star_option_dict['identification']['gtf_file']
    read_dataset_id = star_option_dict['identification']['read_dataset_id']
    index_building = star_option_dict['STAR parameters']['index_building']
    limit_genome_generate_ram = star_option_dict['STAR parameters']['limit_genome_generate_ram']
    threads = star_option_dict['STAR parameters']['threads']
    two_pass_mode = star_option_dict['STAR parameters']['two_pass_mode']
    two_pass_1_readsn = star_option_dict['STAR parameters']['two_pass_1_readsn']
    out_filter_multimap_nmax = star_option_dict['STAR parameters']['out_filter_multimap_nmax']
    out_reads_unmapped = star_option_dict['STAR parameters']['out_reads_unmapped']
    cleanup = star_option_dict['STAR parameters']['cleanup']
    other_parameters = star_option_dict['STAR parameters']['other_parameters']
    read_type = star_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in star_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build read file lists
    read_file_1_list = []
    read_file_2_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = star_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            read_file_1_list.append(read_file_1)
            if read_type.upper() == 'PE':
                read_file_2 = star_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                read_file_2_list.append(read_file_2)

    # set the reference file path
    reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

    # set the gtf file path
    gtf_file = xlib.get_cluster_reference_file(reference_dataset_id, gtf_file)

    # set the STAR indexes directory
    (reference_file_name, _) = os.path.splitext(reference_file)
    star_indexes_dir = f'{reference_file_name}_star_indexes'

    # write the STAR process script
    try:
        if not os.path.exists(os.path.dirname(get_star_process_script())):
            os.makedirs(os.path.dirname(get_star_process_script()))
        with open(get_star_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write(f'CURRENT_DIR={current_run_dir}\n')
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
            script_file_id.write( 'function print_star_version\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_star_anaconda_code()}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    STAR --version\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '}\n')
            if index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function create_star_indexes\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    source activate {xlib.get_star_anaconda_code()}\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Creating indexes ..."\n')
                script_file_id.write(f'    rm -rf {star_indexes_dir}\n')
                script_file_id.write(f'    mkdir --parents {star_indexes_dir}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        STAR \\\n')
                script_file_id.write( '            --runMode genomeGenerate \\\n')
                script_file_id.write(f'            --limitGenomeGenerateRAM {limit_genome_generate_ram} \\\n')
                script_file_id.write(f'            --runThreadN {threads} \\\n')
                script_file_id.write(f'            --genomeDir {star_indexes_dir} \\\n')
                script_file_id.write(f'            --genomeFastaFiles {reference_file}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error STAR $RC; fi\n')
                script_file_id.write( '    echo "Indexes are built."\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_star_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_star_anaconda_code()}\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            for i in range(len(read_file_1_list)):
                # set the library name
                if read_file_1.endswith('.gz'):
                    (library_name, _) = os.path.splitext(os.path.basename(read_file_1_list[i][:-3]))
                    gz_compression = True
                else:
                    (library_name, _) = os.path.splitext(os.path.basename(read_file_1_list[i]))
                    gz_compression = False
                # write the instructions for the library
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "Mapping reads of {library_name} ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        STAR \\\n')
                script_file_id.write( '            --runMode alignReads \\\n')
                script_file_id.write(f'            --runThreadN {threads} \\\n')
                script_file_id.write(f'            --genomeDir {star_indexes_dir} \\\n')
                if gz_compression:
                    script_file_id.write( '            --readFilesCommand gzip -c \\\n')
                if read_type.upper() == 'SE':
                    script_file_id.write(f'            --readFilesIn {read_file_1_list[i]} \\\n')
                elif read_type.upper() == 'PE':
                    script_file_id.write(f'            --readFilesIn {read_file_1_list[i]} {read_file_2_list[i]} \\\n')
                script_file_id.write(f'            --sjdbGTFfile {gtf_file} \\\n')
                if two_pass_mode.upper() == 'NONE':
                    script_file_id.write( '            --twopassMode None \\\n')
                elif two_pass_mode.upper() == 'BASIC':
                    script_file_id.write( '            --twopassMode Basic \\\n')
                    script_file_id.write(f'            --twopass1readsN {two_pass_1_readsn} \\\n')
                script_file_id.write( '            --quantMode TranscriptomeSAM \\\n')
                if out_reads_unmapped.upper() == 'NONE':
                    script_file_id.write( '            --outReadsUnmapped None \\\n')
                elif out_reads_unmapped.upper() == 'FASTX':
                    script_file_id.write( '            --outReadsUnmapped Fastx \\\n')
                script_file_id.write(f'            --outFilterMultimapNmax {out_filter_multimap_nmax} \\\n')
                script_file_id.write( '            --outSAMunmapped None \\\n')
                script_file_id.write(f'            --outFileNamePrefix "{current_run_dir}/{library_name}-" \\\n')
                if cleanup.upper() == 'NO':
                    script_file_id.write( '            --outTmpKeep Yes \\\n')
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
                script_file_id.write( '            --outSAMtype SAM \\\n')
                script_file_id.write( '            --outSAMattributes NH HI AS nM XS\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error STAR $RC; fi\n')
                script_file_id.write( '    echo "Reads are mapped."\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function convert_sam2bam\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Converting SAM files to BAM format ..."\n')
            script_file_id.write(f'    source activate {xlib.get_samtools_anaconda_code()}\n')
            script_file_id.write( '    ls *.sam > sam-files.txt\n')
            script_file_id.write( '    while read FILE_SAM; do\n')
            script_file_id.write( '        FILE_BAM=`basename $FILE_SAM | sed "s|.sam|.bam|g"`\n')
            script_file_id.write( '        echo "Converting file $FILE_SAM to BAM format ..."\n')
            script_file_id.write(f'        samtools view --threads {threads} -b -S -o $FILE_BAM $FILE_SAM\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error samtools-view $RC; fi\n')
            script_file_id.write( '        echo "$FILE_BAM is created."\n')
            script_file_id.write( '        echo "Compressing $FILE_SAM ..."\n')
            script_file_id.write( '        gzip $FILE_SAM\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error gzip $RC; fi\n')
            script_file_id.write( '        echo "$FILE_SAM is compressed."\n')
            script_file_id.write( '    done < sam-files.txt\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function sort_and_index_bam_files\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    source activate {xlib.get_samtools_anaconda_code()}\n')
            script_file_id.write( '    ls *-Aligned.out.bam > bam-files.txt\n')
            script_file_id.write( '    while read FILE_BAM; do\n')
            script_file_id.write( '        FILE_SORTED_BAM=`basename $FILE_BAM | sed "s|.bam|.sorted.bam|g"`\n')
            script_file_id.write( '        echo "Sorting and indexing $FILE_BAM ..."\n')
            script_file_id.write(f'        samtools sort --threads {threads} $FILE_BAM -o $FILE_SORTED_BAM\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error samtools-sort $RC; fi\n')
            script_file_id.write(f'        samtools index -@ {threads} $FILE_SORTED_BAM\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error samtools-index $RC; fi\n')
            script_file_id.write( '        echo "$FILE_SORTED_BAM is created."\n')
            script_file_id.write( '        echo "Deleting file $FILE_BAM ..."\n')
            script_file_id.write( '        rm -f $FILE_BAM\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '        echo "$FILE_BAM is deleted."\n')
            script_file_id.write( '    done < bam-files.txt\n')
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
            process_name = f'{xlib.get_star_name()} process'
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
            script_file_id.write( 'print_star_version\n')
            if index_building.upper() == 'YES':
                script_file_id.write( 'create_star_indexes\n')
            script_file_id.write( 'run_star_process\n')
            script_file_id.write( 'convert_sam2bam\n')
            script_file_id.write( 'sort_and_index_bam_files\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_star_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_star_process_starter(current_run_dir):
    '''
    Build the starter of the current STAR process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the STAR process starter
    try:
        if not os.path.exists(os.path.dirname(get_star_process_starter())):
            os.makedirs(os.path.dirname(get_star_process_starter()))
        with open(get_star_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_star_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_star_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_star_config_file():
    '''
    Get the STAR config file path.
    '''

    # assign the STAR config file path
    star_config_file = f'{xlib.get_config_dir()}/{xlib.get_star_code()}-config.txt'

    # return the STAR config file path
    return star_config_file

#-------------------------------------------------------------------------------

def get_star_process_script():
    '''
    Get the STAR process script path in the local computer.
    '''

    # assign the STAR script path
    star_process_script = f'{xlib.get_temp_dir()}/{xlib.get_star_code()}-process.sh'

    # return the STAR script path
    return star_process_script

#-------------------------------------------------------------------------------

def get_star_process_starter():
    '''
    Get the STAR process starter path in the local computer.
    '''

    # assign the STAR process starter path
    star_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_star_code()}-process-starter.sh'

    # return the STAR starter path
    return star_process_starter

#-------------------------------------------------------------------------------
    
def get_index_building_code_list():
    '''
    Get the code list of "index_building".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_index_building_code_list_text():
    '''
    Get the code list of "index_building" as text.
    '''

    return 'YES (built indexes) or NO (old indexes will be used)'

#-------------------------------------------------------------------------------
    
def get_two_pass_mode_code_list():
    '''
    Get the code list of "two_pass_mode".
    '''

    return ['NONE', 'BASIC']

#-------------------------------------------------------------------------------
    
def get_two_pass_mode_code_list_text():
    '''
    Get the code list of "two_pass_mode" as text.
    '''

    return 'NONE (1-pass mapping) or BASIC (basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly)'

#-------------------------------------------------------------------------------
    
def get_out_reads_unmapped_code_list():
    '''
    Get the code list of "out_reads_unmapped".
    '''

    return ['NONE', 'FASTX']

#-------------------------------------------------------------------------------
    
def get_out_reads_unmapped_code_list_text():
    '''
    Get the code list of "out_reads_unmapped" as text.
    '''

    return 'NONE (no output) or FASTX (output in separate fasta/fastq files)'

#-------------------------------------------------------------------------------
    
def get_cleanup_code_list():
    '''
    Get the code list of "cleanup".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_cleanup_code_list_text():
    '''
    Get the code list of "cleanup" as text.
    '''

    return str(get_cleanup_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
     print('This file contains functions related to the STAR process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
