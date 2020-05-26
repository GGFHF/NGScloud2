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
This file contains functions related to the SOAPdenovo2 process used in
both console mode and gui mode.
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

def create_soapdenovo2_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create SOAPdenovo2 config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the SOAPdenovo2 config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_soapdenovo2_config_file())):
            os.makedirs(os.path.dirname(get_soapdenovo2_config_file()))
        with open(get_soapdenovo2_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id and read_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of SOAPdenovo2 and their meaning in "http://soap.genomics.org.cn/SOAPdenovo2.html".\n')
            file_id.write( '#\n')
            file_id.write( '# There are two formats to set an option:\n')
            file_id.write( '#\n')
            file_id.write( '#    option = value                             <- if the option supports a single value\n')
            file_id.write( '#\n')
            file_id.write( '#    option = value-1, value-2, ..., value-n    <- if the option supports a values list\n')
            file_id.write( '#\n')
            file_id.write( '# WARNING: The files have to be decompressed.\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information that identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the SOAPdenovo2 parameters\n')
            file_id.write( '[SOAPdenovo2 parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'version = 63', f'# SOAPdenovo2 version: {get_version_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'ncpu = 4', '# number of cpu for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'init_memory_assumption = 0', '# memory assumption initialized to avoid further reallocation in GiB'))
            file_id.write( '{0:<50} {1}\n'.format( 'kmer = 25', '# value or values list of kmer size: minimum, 13; maximum: version value.'))
            file_id.write( '{0:<50} {1}\n'.format( 'kmer_freq_cutoff = 0', '# kmers with frequency no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format( 'edge_cov_cutoff = 2', '# edges with coverage no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format( 'resolve_repeats = NO', f'# resolve repeats by reads: {get_resolve_repeats_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'merge_level = 1', '# strength of merging similar sequences during contiging: minimum, 0; maximum, 3.'))
            file_id.write( '{0:<50} {1}\n'.format( 'filter = 0', '# weight to filter arc when linearize two edges'))
            file_id.write( '{0:<50} {1}\n'.format( 'merge_bubble = NO', f'# merge clean bubble before iterate: {get_merge_bubble_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'srkgf = NO', f'# output gap related redas for SRkgf to fill gap: {get_srkgf_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'fill = NO', f'# fill gaps in scaffolds: {get_fill_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'unmask = NO', f'# un-mask contigs with high/low coverage before scaffolding: {get_unmask_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'keep_contigs = NO', f'# keep contigs weakly connected to other contigs in scaffold: {get_keep_contigs_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'gap_len_diff = 50', '# allowed length difference between estimated and filled gap'))
            file_id.write( '{0:<50} {1}\n'.format( 'min_contig_len = 27', '# shortest contig for scaffolding'))
            file_id.write( '{0:<50} {1}\n'.format( 'min_contig_cvg = 0.1', '# minimum contig coverage'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_contig_cvg = 2.0', '# maximum contig coverage'))
            file_id.write( '{0:<50} {1}\n'.format( 'insert_size_upper_bound = 1.5', '# upper bound of insert size for large insert size'))
            file_id.write( '{0:<50} {1}\n'.format( 'bubble_coverage = 0.6', '# lower coverage in bubble structure'))
            file_id.write( '{0:<50} {1}\n'.format( 'genome_size = 0', '# genome size for statistics'))
            file_id.write( '{0:<50} {1}\n'.format( 'visualization = NO', f'# output visualization information of assembly: {get_visualization_code_list_text()}'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'max_rd_len = 100', '# any read longer than the value will be cut to this length'))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '# This section has the information of the first library.\n')
                file_id.write(f'[library-{i + 1}]\n')
                file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
                file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
                file_id.write( '{0:<50} {1}\n'.format(f'read_file_1 = {os.path.basename(file_1_list[i])}', '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format( 'read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format(f'read_file_2 = {os.path.basename(file_2_list[i])}', '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                file_id.write( '{0:<50} {1}\n'.format( 'avg_ins = 200', '# average insert size'))
                file_id.write( '{0:<50} {1}\n'.format( 'reverse_seq = 0', f'# sequence direction: {get_reverse_seq_code_list_text()}'))
                file_id.write( '{0:<50} {1}\n'.format( 'asm_flags = 3', f'# in which part(s) the reads are used: {get_asm_flags_code_list_text()}'))
                file_id.write( '{0:<50} {1}\n'.format( 'rd_len_cutof = 100', '# reads will be cut to this length'))
                file_id.write( '{0:<50} {1}\n'.format( 'pair_num_cutoff = 3', '# cutoff value of pair number for a reliable connection between two contigs or pre-scaffolds (minimum: 3)'))
                file_id.write( '{0:<50} {1}\n'.format( 'map_len = 32', '# minimum aligned length between a read and a contigs for a reliable read location (minimum: 32)'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '# If there are more libraries, you have to repeat the section library-1 with the data of each file.\n')
                    file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_soapdenovo2_config_file()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_soapdenovo2_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in SOAPdenovo2 config file.
    '''

    # initialize the control variable
    OK = True

    # get the SOAPdenovo2 option dictionary
    soapdenovo2_option_dict = xlib.get_option_dict(get_soapdenovo2_config_file())

    # get the experiment identification
    experiment_id = soapdenovo2_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the SOAPdenovo2 configuration file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_soapdenovo2_name()} config file ...\n')
    (OK, error_list) = check_soapdenovo2_config_file(strict=True)
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

    # check the SOAPdenovo2 is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_soapdenovo2_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_soapdenovo2_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_soapdenovo2_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # build the process configuration file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process config file {get_soapdenovo2_process_config_file()} ...\n')
        (OK, error_list) = build_soapdenovo2_process_config_file()
        if OK:
            log.write('The file is built.\n')
        else:
            log.write('*** ERROR: The file could not be built.\n')

    # for each kmer value, build the process, copy it the cluster and run it
    if OK:

        # get the kmer list
        kmer = soapdenovo2_option_dict['SOAPdenovo2 parameters']['kmer']
        kmer_value_list = xlib.split_literal_to_integer_list(kmer)
        
        # for each kmer value, do the tasks
        i = 1
        for kmer_value in kmer_value_list:

            # determine the run directory in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Determining the run directory for kmer {kmer_value} in the cluster is being determined ...\n')
            if i > 1:
                current_run_dir = f'{xlib.get_cluster_current_run_dir(experiment_id, xlib.get_soapdenovo2_code())}-{i}'
            else:
                current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_soapdenovo2_code())
            command = f'mkdir --parents {current_run_dir}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write(f'The directory path is {current_run_dir}.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            i += 1

            # build the SOAPdenovo2 process script
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Building the process script {get_soapdenovo2_process_script()} ...\n')
            (OK, error_list) = build_soapdenovo2_process_script(cluster_name, current_run_dir, kmer_value)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')
                break

            # upload the process configuration file to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the process config file {get_soapdenovo2_process_config_file()} to the directory {current_run_dir} ...\n')
            cluster_path = f'{current_run_dir}/{os.path.basename(get_soapdenovo2_process_config_file())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_soapdenovo2_process_config_file(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # upload the process script to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the process script {get_soapdenovo2_process_script()} to the directory {current_run_dir} ...\n')
            cluster_path = f'{current_run_dir}/{os.path.basename(get_soapdenovo2_process_script())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_soapdenovo2_process_script(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process script in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_soapdenovo2_process_script())} ...\n')
            command = f'chmod u+x {current_run_dir}/{os.path.basename(get_soapdenovo2_process_script())}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # build the process starter
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Building the process starter {get_soapdenovo2_process_starter()} ...\n')
            (OK, error_list) = build_soapdenovo2_process_starter(current_run_dir)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('***ERROR: The file could not be built.\n')
                break

            # upload the process starter to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Uploading the process starter {get_soapdenovo2_process_starter()} to the directory {current_run_dir} ...\n')
            cluster_path = f'{current_run_dir}/{os.path.basename(get_soapdenovo2_process_starter())}'
            (OK, error_list) = xssh.put_file(sftp_client, get_soapdenovo2_process_starter(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process starter in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_soapdenovo2_process_starter())} ...\n')
            command = f'chmod u+x {current_run_dir}/{os.path.basename(get_soapdenovo2_process_starter())}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # submit the process
            log.write(f'{xlib.get_separator()}\n')
            log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_soapdenovo2_process_starter())} ...\n')
            OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_soapdenovo2_process_starter()), log)

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

def check_soapdenovo2_config_file(strict):
    '''
    Check the SOAPdenovo2 config file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'

    # get the option dictionary
    try:
        soapdenovo2_option_dict = xlib.get_option_dict(get_soapdenovo2_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in soapdenovo2_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = soapdenovo2_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the sectfillion "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = soapdenovo2_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "SOAPdenovo2 parameters"
        if 'SOAPdenovo2 parameters' not in sections_list:
            error_list.append('*** ERROR: the section "SOAPdenovo2 parameters" is not found.')
            OK = False
        else:

            # check section "SOAPdenovo2 parameters" - key "version"
            version = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('version', not_found)
            if version == not_found:
                error_list.append('*** ERROR: the key "version" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(version, get_version_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "version" has to be {get_version_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "ncpu"
            ncpu = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "init_memory_assumption"
            init_memory_assumption = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('init_memory_assumption', not_found)
            if init_memory_assumption == not_found:
                error_list.append('*** ERROR: the key "init_memory_assumption" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(init_memory_assumption, minimum=0):
                error_list.append('*** ERROR: the key "init_memory_assumption" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "kmer"
            kmer = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            else:
                kmer_list = xlib.split_literal_to_integer_list(kmer)
                if kmer_list == []:
                    error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 13 and the version value.')
                    OK = False
                else:
                    for kmer_item in kmer_list:
                        if not xlib.check_int(kmer_item, minimum=13, maximum=int(version)):
                            error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 13 and the version value.')
                            OK = False
                            break

            # check section "SOAPdenovo2 parameters" - key "kmer_freq_cutoff"
            kmer_freq_cutoff = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('kmer_freq_cutoff', not_found)
            if kmer_freq_cutoff == not_found:
                error_list.append('*** ERROR: the key "kmer_freq_cutoff" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(kmer_freq_cutoff, minimum=0):
                error_list.append('*** ERROR: the key "kmer_freq_cutoff" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "edge_cov_cutoff"
            edge_cov_cutoff = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('edge_cov_cutoff', not_found)
            if edge_cov_cutoff == not_found:
                error_list.append('*** ERROR: the key "edge_cov_cutoff" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(edge_cov_cutoff, minimum=0):
                error_list.append('*** ERROR: the key "edge_cov_cutoff" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "resolve_repeats"
            resolve_repeats = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('resolve_repeats', not_found)
            if resolve_repeats == not_found:
                error_list.append('*** ERROR: the key "resolve_repeats" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(resolve_repeats, get_resolve_repeats_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "resolve_repeats" has to be {get_resolve_repeats_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "merge_level"
            merge_level = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('merge_level', not_found)
            if merge_level == not_found:
                error_list.append('*** ERROR: the key "merge_level" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(merge_level, minimum=0, maximum=3):
                error_list.append('*** ERROR: the key "merge_level" has to be an integer number between 0 and 3.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "filter"
            filter = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('filter', not_found)
            if filter == not_found:
                error_list.append('*** ERROR: the key "filter" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(filter, minimum=0):
                error_list.append('*** ERROR: the key "filter" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "merge_bubble"
            merge_bubble = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('merge_bubble', not_found)
            if merge_bubble == not_found:
                error_list.append('*** ERROR: the key "merge_bubble" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(merge_bubble, get_merge_bubble_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "merge_bubble" has to be {get_merge_bubble_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "srkgf"
            srkgf = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('srkgf', not_found)
            if srkgf == not_found:
                error_list.append('*** ERROR: the key "srkgf" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(srkgf, get_srkgf_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "srkgf" has to be {get_srkgf_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "fill"
            fill = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('fill', not_found)
            if fill == not_found:
                error_list.append('*** ERROR: the key "fill" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(fill, get_fill_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "fill" has to be {get_fill_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "unmask"
            unmask = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('unmask', not_found)
            if unmask == not_found:
                error_list.append('*** ERROR: the key "unmask" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(unmask, get_unmask_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "unmask" has to be {get_unmask_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "keep_contigs"
            keep_contigs = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('keep_contigs', not_found)
            if keep_contigs == not_found:
                error_list.append('*** ERROR: the key "keep_contigs" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(keep_contigs, get_keep_contigs_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "keep_contigs" has to be {get_keep_contigs_code_list_text()}.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "gap_len_diff"
            gap_len_diff = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('gap_len_diff', not_found)
            if gap_len_diff == not_found:
                error_list.append('*** ERROR: the key "gap_len_diff" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(gap_len_diff, minimum=0):
                    error_list.append('*** ERROR: the key "gap_len_diff" has to be an integer number greater than or equal to 0.')
                    OK = False

            # check section "SOAPdenovo2 parameters" - key "min_contig_len"
            min_contig_len = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('min_contig_len', not_found)
            if min_contig_len == not_found:
                error_list.append('*** ERROR: the key "min_contig_len" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(min_contig_len, minimum=0):
                error_list.append('*** ERROR: the key "min_contig_len" is not an integer value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "min_contig_cvg"
            min_contig_cvg = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('min_contig_cvg', not_found)
            if min_contig_cvg == not_found:
                error_list.append('*** ERROR: the key "min_contig_cvg" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_float(min_contig_cvg, minimum=0., mne=0.):
                error_list.append('*** ERROR: the key "min_contig_cvg" is not a float value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "max_contig_cvg"
            max_contig_cvg = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('max_contig_cvg', not_found)
            if max_contig_cvg == not_found:
                error_list.append('*** ERROR: the key "max_contig_cvg" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_float(max_contig_cvg, minimum=0., mne=0.):
                error_list.append('*** ERROR: the key "max_contig_cvg" is not a float value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "insert_size_upper_bound"
            insert_size_upper_bound = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('insert_size_upper_bound', not_found)
            if insert_size_upper_bound == not_found:
                error_list.append('*** ERROR: the key "insert_size_upper_bound" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_float(insert_size_upper_bound, minimum=0., mne=0.):
                error_list.append('*** ERROR: the key "insert_size_upper_bound" is not a float value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "bubble_coverage"
            bubble_coverage = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('bubble_coverage', not_found)
            if bubble_coverage == not_found:
                error_list.append('*** ERROR: the key "bubble_coverage" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_float(bubble_coverage, minimum=0., mne=0.):
                error_list.append('*** ERROR: the key "bubble_coverage" is not a float value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "genome_size"
            genome_size = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('genome_size', not_found)
            if genome_size == not_found:
                error_list.append('*** ERROR: the key "genome_size" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(genome_size, minimum=0):
                error_list.append('*** ERROR: the key "genome_size" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "visualization"
            visualization = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('visualization', not_found)
            if visualization == not_found:
                error_list.append('*** ERROR: the key "visualization" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(visualization, get_visualization_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "visualization" has to be {get_visualization_code_list_text()}.')
                OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "max_rd_len"
            max_rd_len = soapdenovo2_option_dict.get('library', {}).get('max_rd_len', not_found)
            if max_rd_len == not_found:
                error_list.append('*** ERROR: the key "max_rd_len" is not found in the section "library".')
                OK = False
            elif not xlib.check_int(max_rd_len, minimum=0):
                error_list.append('*** ERROR: the key "max_rd_len" has to be an integer number greater than or equal to 0.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'SOAPdenovo2 parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "format"
                    format = soapdenovo2_option_dict.get(section, {}).get('format', not_found)
                    if format == not_found:
                        error_list.append(f'*** ERROR: the key "format" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                        error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                        OK = False

                    # check section "library-n" - key "read_type"
                    read_type = soapdenovo2_option_dict.get(section, {}).get('read_type', not_found)
                    if read_type == not_found:
                        error_list.append(f'*** ERROR: the key "read_type" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                        error_list.append(f'*** ERROR: the key "read_type" has to be {get_read_type_code_list_text()}.')
                        OK = False
                    elif read_type == 'SP' and format != 'FASTA':
                        error_list.append('*** ERROR: if read type is SP, the format has to be FASTA.')
                        OK = False

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = soapdenovo2_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False
                    elif read_file_1.find('.gz') != -1:
                        error_list.append(f'*** ERROR: the key "read_file_1" in the section "{section}" has to be a decompressed file.')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = soapdenovo2_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False
                    elif read_file_2.find('.gz') != -1:
                        error_list.append(f'*** ERROR: the key "read_file_2" in the section "{section}" has to be a decompressed file.')
                        OK = False

                    # check section "library" - key "avg_ins"
                    avg_ins = soapdenovo2_option_dict.get(section, {}).get('avg_ins', not_found)
                    if avg_ins == not_found:
                        error_list.append(f'*** ERROR: the key "avg_ins" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_int(avg_ins, minimum=0):
                        error_list.append(f'*** ERROR: the key "avg_ins" in the section "{section}" has to be an integer number greater than or equal to 0.')
                        OK = False

                    # check section "library-n" - key "reverse_seq"
                    reverse_seq = soapdenovo2_option_dict.get(section, {}).get('reverse_seq', not_found)
                    if reverse_seq == not_found:
                        error_list.append(f'*** ERROR: the key "reverse_seq" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_code(reverse_seq, get_reverse_seq_code_list(), case_sensitive=False):
                        error_list.append(f'*** ERROR: the key "reverse_seq" has to be {get_reverse_seq_code_list_text()}.')
                        OK = False

                    # check section "library-n" - key "asm_flags"
                    asm_flags = soapdenovo2_option_dict.get(section, {}).get('asm_flags', not_found)
                    if asm_flags == not_found:
                        error_list.append(f'*** ERROR: the key "asm_flags" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_code(asm_flags, get_asm_flags_code_list(), case_sensitive=False):
                        error_list.append(f'*** ERROR: the key "asm_flags" has to be {get_asm_flags_code_list_text()}.')
                        OK = False

                    # check section "library-n" - key "rd_len_cutof"
                    rd_len_cutof = soapdenovo2_option_dict.get(section, {}).get('rd_len_cutof', not_found)
                    if rd_len_cutof == not_found:
                        error_list.append(f'*** ERROR: the key "rd_len_cutof" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_int(rd_len_cutof, minimum=0):
                        error_list.append(f'*** ERROR: the key "rd_len_cutof" in the section "{section}" has to be an integer number greater than or equal to 1.')
                        OK = False

                    # check section "library-n" - key "pair_num_cutoff"
                    pair_num_cutoff = soapdenovo2_option_dict.get(section, {}).get('pair_num_cutoff', not_found)
                    if pair_num_cutoff == not_found:
                        error_list.append(f'*** ERROR: the key "pair_num_cutoff" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_int(pair_num_cutoff, minimum=0):
                        error_list.append(f'*** ERROR: the key "pair_num_cutoff" in the section "{section}" has to be an integer number greater than or equal to 1.')
                        OK = False

                    # check section "library-n" - key "map_len"
                    map_len = soapdenovo2_option_dict.get(section, {}).get('map_len', not_found)
                    if map_len == not_found:
                        error_list.append(f'*** ERROR: the key "map_len" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_int(map_len, minimum=32):
                        error_list.append(f'*** ERROR: the key "map_len" in the section "{section}" has to be an integer number greater than or equal to 32.')
                        OK = False

    # warn that the SOAPdenovo2 config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_soapdenovo2_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovo2_process_script(cluster_name, current_run_dir, kmer_value):
    '''
    Build the current SOAPdenovo2 process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    soapdenovo2_options_dict = xlib.get_option_dict(get_soapdenovo2_config_file())

    # get the options
    experiment_id = soapdenovo2_options_dict['identification']['experiment_id']
    version = soapdenovo2_options_dict['SOAPdenovo2 parameters']['version']
    ncpu = soapdenovo2_options_dict['SOAPdenovo2 parameters']['ncpu']
    init_memory_assumption = soapdenovo2_options_dict['SOAPdenovo2 parameters']['init_memory_assumption']
    kmer_freq_cutoff = soapdenovo2_options_dict['SOAPdenovo2 parameters']['kmer_freq_cutoff']
    edge_cov_cutoff = soapdenovo2_options_dict['SOAPdenovo2 parameters']['edge_cov_cutoff']
    resolve_repeats = soapdenovo2_options_dict['SOAPdenovo2 parameters']['resolve_repeats'].upper()
    merge_level = soapdenovo2_options_dict['SOAPdenovo2 parameters']['merge_level']
    filter = soapdenovo2_options_dict['SOAPdenovo2 parameters']['filter']
    merge_bubble = soapdenovo2_options_dict['SOAPdenovo2 parameters']['merge_bubble'].upper()
    srkgf = soapdenovo2_options_dict['SOAPdenovo2 parameters']['srkgf'].upper()
    fill = soapdenovo2_options_dict['SOAPdenovo2 parameters']['fill'].upper()
    unmask = soapdenovo2_options_dict['SOAPdenovo2 parameters']['unmask'].upper()
    keep_contigs = soapdenovo2_options_dict['SOAPdenovo2 parameters']['keep_contigs'].upper()
    gap_len_diff = soapdenovo2_options_dict['SOAPdenovo2 parameters']['gap_len_diff']
    min_contig_len = soapdenovo2_options_dict['SOAPdenovo2 parameters']['min_contig_len']
    min_contig_cvg = soapdenovo2_options_dict['SOAPdenovo2 parameters']['min_contig_cvg']
    max_contig_cvg = soapdenovo2_options_dict['SOAPdenovo2 parameters']['max_contig_cvg']
    insert_size_upper_bound = soapdenovo2_options_dict['SOAPdenovo2 parameters']['insert_size_upper_bound']
    bubble_coverage = soapdenovo2_options_dict['SOAPdenovo2 parameters']['bubble_coverage']
    genome_size = soapdenovo2_options_dict['SOAPdenovo2 parameters']['genome_size']
    visualization = soapdenovo2_options_dict['SOAPdenovo2 parameters']['visualization'].upper()

    # write the SOAPdenovo2 process script
    try:
        if not os.path.exists(os.path.dirname(get_soapdenovo2_process_script())):
            os.makedirs(os.path.dirname(get_soapdenovo2_process_script()))
        with open(get_soapdenovo2_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'ulimit -s unlimited\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
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
            script_file_id.write( 'function run_soapdenovo2_pregraph\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_soapdenovo2_pregraph.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running pregraph process ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_soapdenovo2_anaconda_code()}\n')
            script_file_id.write(f'        cd {current_run_dir}\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write(f'            SOAPdenovo-{version}mer pregraph \\\n')
            script_file_id.write(f'                -p {ncpu} \\\n')
            script_file_id.write(f'                -s {current_run_dir}/{os.path.basename(get_soapdenovo2_process_config_file())} \\\n')
            script_file_id.write(f'                -a {init_memory_assumption} \\\n')
            script_file_id.write(f'                -K {kmer_value} \\\n')
            script_file_id.write(f'                -d {kmer_freq_cutoff} \\\n')
            if resolve_repeats == 'YES':
                script_file_id.write( '                -R \\\n')
            script_file_id.write(f'                -o {experiment_id}-{os.path.basename(current_run_dir)}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write(f'        if [ $RC -ne 0 ]; then manage_error SOAPdenovo-{version}mer $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "pregraph process is ended."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_soapdenovo2_contig\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_soapdenovo2_contig.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running contig process ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_soapdenovo2_anaconda_code()}\n')
            script_file_id.write(f'        cd {current_run_dir}\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write(f'            SOAPdenovo-{version}mer contig \\\n')
            script_file_id.write(f'                -p {ncpu} \\\n')
            script_file_id.write(f'                -s {current_run_dir}/{os.path.basename(get_soapdenovo2_process_config_file())} \\\n')
            script_file_id.write(f'                -D {edge_cov_cutoff} \\\n')
            if resolve_repeats == 'YES':
                script_file_id.write( '                -R \\\n')
            script_file_id.write(f'                -M {merge_level} \\\n')
            script_file_id.write(f'                -e {filter} \\\n')
            if merge_bubble == 'YES':
                script_file_id.write( '                -E \\\n')
            script_file_id.write(f'                -g {experiment_id}-{os.path.basename(current_run_dir)}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write(f'        if [ $RC -ne 0 ]; then manage_error SOAPdenovo-{version}mer $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "contig process is ended."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_soapdenovo2_map\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_soapdenovo2_map.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running map process ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_soapdenovo2_anaconda_code()}\n')
            script_file_id.write(f'        cd {current_run_dir}\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write(f'            SOAPdenovo-{version}mer map \\\n')
            script_file_id.write(f'                -p {ncpu} \\\n')
            script_file_id.write(f'                -s {current_run_dir}/{os.path.basename(get_soapdenovo2_process_config_file())} \\\n')
            script_file_id.write(f'                -K {kmer_value} \\\n')
            if srkgf == 'YES':
                script_file_id.write( '                -f \\\n')
            script_file_id.write(f'                -g {experiment_id}-{os.path.basename(current_run_dir)}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write(f'        if [ $RC -ne 0 ]; then manage_error SOAPdenovo-{version}mer $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "map process is ended."\n')
            script_file_id.write( '        touch $STEP_STATUS\n')
            script_file_id.write( '    fi\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_soapdenovo2_scaff\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    STEP_STATUS=$STATUS_DIR/run_soapdenovo2_scaff.ok\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Running scaff process ..."\n')
            script_file_id.write( '    if [ -f $STEP_STATUS ]; then\n')
            script_file_id.write( '        echo "This step was previously run."\n')
            script_file_id.write( '    else\n')
            script_file_id.write(f'        source activate {xlib.get_soapdenovo2_anaconda_code()}\n')
            script_file_id.write(f'        cd {current_run_dir}\n')
            script_file_id.write( '        /usr/bin/time \\\n')
            script_file_id.write(f'            --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write(f'            SOAPdenovo-{version}mer scaff \\\n')
            script_file_id.write(f'                -p {ncpu} \\\n')
            if fill == 'YES':
                script_file_id.write( '                -F \\\n')
            if unmask == 'YES':
                script_file_id.write( '                -u \\\n')
            if keep_contigs == 'YES':
                script_file_id.write( '                -w \\\n')
            if visualization == 'YES':
                script_file_id.write(f'                -V {visualization} \\\n')
            script_file_id.write(f'                -G {gap_len_diff} \\\n')
            script_file_id.write(f'                -L {min_contig_len} \\\n')
            script_file_id.write(f'                -c {min_contig_cvg} \\\n')
            script_file_id.write(f'                -C {max_contig_cvg} \\\n')
            script_file_id.write(f'                -b {insert_size_upper_bound} \\\n')
            script_file_id.write(f'                -B {bubble_coverage} \\\n')
            script_file_id.write(f'                -N {genome_size} \\\n')
            script_file_id.write(f'                -g {experiment_id}-{os.path.basename(current_run_dir)}\n')
            script_file_id.write( '        RC=$?\n')
            script_file_id.write(f'        if [ $RC -ne 0 ]; then manage_error SOAPdenovo-{version}mer $RC; fi\n')
            script_file_id.write( '        conda deactivate\n')
            script_file_id.write( '        echo "scaff process is ended."\n')
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
            process_name = f'{xlib.get_soapdenovo2_name()} process'
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
            script_file_id.write( 'run_soapdenovo2_pregraph\n')
            script_file_id.write( 'run_soapdenovo2_contig\n')
            script_file_id.write( 'run_soapdenovo2_map\n')
            script_file_id.write( 'run_soapdenovo2_scaff\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_soapdenovo2_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovo2_process_starter(current_run_dir):
    '''
    Build the starter of the current SOAPdenovo2 process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the SOAPdenovo2 process starter
    try:
        if not os.path.exists(os.path.dirname(get_soapdenovo2_process_starter())):
            os.makedirs(os.path.dirname(get_soapdenovo2_process_starter()))
        with open(get_soapdenovo2_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_soapdenovo2_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_soapdenovo2_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def restart_soapdenovo2_process(cluster_name, experiment_id, result_dataset_id, log, function=None):
    '''
    Restart a SOAPdenovo2 process from the last step ended OK.
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
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_soapdenovo2_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_soapdenovo2_process_starter()), log)

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

def build_soapdenovo2_process_config_file():
    '''
    Build the SOAPdenovo2 process config file to the current SOPAdenovo2 experiment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    soapdenovo2_option_dict = xlib.get_option_dict(get_soapdenovo2_config_file())

    # get the sections list
    sections_list = []
    for section in soapdenovo2_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # get the experiment identification and read dataset identification
    experiment_id = soapdenovo2_option_dict['identification']['experiment_id']
    read_dataset_id = soapdenovo2_option_dict['identification']['read_dataset_id']

    # write the SOAPdenovo2 process configuration file
    try:
        if not os.path.exists(os.path.dirname(get_soapdenovo2_process_config_file())):
            os.makedirs(os.path.dirname(get_soapdenovo2_process_config_file()))
        with open(get_soapdenovo2_process_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#maximal read length\n')
            file_id.write(f'max_rd_len={soapdenovo2_option_dict["library"]["max_rd_len"]}\n')
            for section in sections_list:
                if re.match('^library-[0-9]+$', section):
                    format = soapdenovo2_option_dict[section]['format'].upper()
                    read_type = soapdenovo2_option_dict[section]['read_type'].upper()
                    read_file_1 = soapdenovo2_option_dict[section]['read_file_1']
                    read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
                    if read_type == 'PE':
                        read_file_2 = soapdenovo2_option_dict[section]['read_file_2']
                        read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                    file_id.write( '[LIB]\n')
                    file_id.write( '# average insert size\n')
                    file_id.write(f'avg_ins={soapdenovo2_option_dict[section]["avg_ins"]}\n')
                    file_id.write( '# if sequence needs to be reversed\n')
                    file_id.write(f'reverse_seq={soapdenovo2_option_dict[section]["reverse_seq"]}\n')
                    file_id.write( '# in which part(s) the reads are used\n')
                    file_id.write(f'asm_flags={soapdenovo2_option_dict[section]["asm_flags"]}\n')
                    file_id.write( '# maximal read length in this lib\n')
                    file_id.write(f'rd_len_cutof={soapdenovo2_option_dict[section]["rd_len_cutof"]}\n')
                    file_id.write( '# cutoff of pair number for a reliable connection (at least 3 for short insert size)\n')
                    file_id.write(f'pair_num_cutoff={soapdenovo2_option_dict[section]["pair_num_cutoff"]}\n')
                    file_id.write( '# minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\n')
                    file_id.write(f'map_len={soapdenovo2_option_dict[section]["map_len"]}\n')
                    if format == 'FASTA' and read_type == 'SE':
                        file_id.write( '# fasta file for single reads\n')
                        file_id.write(f'f={read_file_1}\n')
                    elif format == 'FASTA' and read_type == 'PE':
                        file_id.write( '# fasta file for read 1\n')
                        file_id.write(f'f1={read_file_1}\n')
                        file_id.write( '# fastq file for read 2 always follows fastq file for read 1\n')
                        file_id.write(f'f2={read_file_2}\n')
                    elif format == 'FASTA' and read_type == 'SP':
                        file_id.write( '# a single fasta file for paired reads\n')
                        file_id.write(f'p={read_file_1}\n')
                    elif format == 'FASTQ' and read_type == 'SE':
                        file_id.write( '# fastq file for single reads\n')
                        file_id.write(f'q={read_file_1}\n')
                    elif format == 'FASTQ' and read_type == 'PE':
                        file_id.write( '# fastq file for read 1\n')
                        file_id.write(f'q1={read_file_1}\n')
                        file_id.write( '# fastq file for read 2 always follows fastq file for read 1\n')
                        file_id.write(f'q2={read_file_2}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_soapdenovo2_process_config_file()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_soapdenovo2_config_file():
    '''
    Get the SOAPdenovo2 config file path.
    '''

    # assign the SOAPdenovo2 config file path
    soapdenovo2_config_file = f'{xlib.get_config_dir()}/{xlib.get_soapdenovo2_code()}-config.txt'

    # return the SOAPdenovo2 config file path
    return soapdenovo2_config_file

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_script():
    '''
    Get the SOAPdenovo2 process script path in the local computer.
    '''

    # assign the SOAPdenovo2 script path
    soapdenovo2_process_script = f'{xlib.get_temp_dir()}/{xlib.get_soapdenovo2_code()}-process.sh'

    # return the SOAPdenovo2 script path
    return soapdenovo2_process_script

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_starter():
    '''
    Get the SOAPdenovo2 process starter path in the local computer.
    '''

    # assign the SOAPdenovo2 process starter path
    soapdenovo2_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_soapdenovo2_code()}-process-starter.sh'

    # return the SOAPdenovo2 starter path
    return soapdenovo2_process_starter

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_config_file():
    '''
    Get the SOAPdenovo2 process config file path in the local computer.
    '''

    # assign the SOAPdenovo2 process config file path
    soapdenovo2_process_config_file = f'{xlib.get_temp_dir()}/{xlib.get_soapdenovo2_code()}-process-config.txt'

    # return the SOAPdenovo2 process config file path
    return soapdenovo2_process_config_file

#-------------------------------------------------------------------------------
    
def get_version_code_list():
    '''
    Get the code list of "version".
    '''

    return ['63', '127']

#-------------------------------------------------------------------------------
    
def get_version_code_list_text():
    '''
    Get the code list of "version" as text.
    '''

    return '63 (SOAPdenovo-63mer) or 127 (SOAPdenovo-127mer)'

#-------------------------------------------------------------------------------
    
def get_resolve_repeats_code_list():
    '''
    Get the code list of "resolve_repeats".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_resolve_repeats_code_list_text():
    '''
    Get the code list of "resolve_repeats" as text.
    '''

    return str(get_resolve_repeats_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_merge_bubble_code_list():
    '''
    Get the code list of "merge_bubble".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_merge_bubble_code_list_text():
    '''
    Get the code list of "merge_bubble" as text.
    '''

    return str(get_merge_bubble_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_srkgf_code_list():
    '''
    Get the code list of "srkgf".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_srkgf_code_list_text():
    '''
    Get the code list of "srkgf" as text.
    '''

    return str(get_srkgf_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_fill_code_list():
    '''
    Get the code list of "fill".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_fill_code_list_text():
    '''
    Get the code list of "fill" as text.
    '''

    return str(get_fill_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_unmask_code_list():
    '''
    Get the code list of "unmask".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_unmask_code_list_text():
    '''
    Get the code list of "unmask" as text.
    '''

    return str(get_unmask_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_keep_contigs_code_list():
    '''
    Get the code list of "keep_contigs".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_keep_contigs_code_list_text():
    '''
    Get the code list of "keep_contigs" as text.
    '''

    return str(get_keep_contigs_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_visualization_code_list():
    '''
    Get the code list of "visualization".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_visualization_code_list_text():
    '''
    Get the code list of "visualization" as text.
    '''

    return str(get_visualization_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_reverse_seq_code_list():
    '''
    Get the code list of "reverse_seq".
    '''

    return ['0', '1']

#-------------------------------------------------------------------------------
    
def get_reverse_seq_code_list_text():
    '''
    Get the code list of "reverse_seq" as text.
    '''

    return '0 (forward-reverse) or 1 (reverse-forward)'

#-------------------------------------------------------------------------------
    
def get_asm_flags_code_list():
    '''
    Get the code list of "asm_flags".
    '''

    return ['1', '2', '3', '4']

#-------------------------------------------------------------------------------
    
def get_asm_flags_code_list_text():
    '''
    Get the code list of "asm_flags" as text.
    '''

    return '1 (only contig assembly) or 2 (only scaffolod assembly) or 3 (both contig and scffold assembly) or 4 (only gap closure)'

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

    return ['SE', 'PE', 'SP']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "read_type" as text.
    '''

    return 'SE (single-end) or PE (pair-end) or SP (single fasta file for paired reads)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the SOAPdenovo2 process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
