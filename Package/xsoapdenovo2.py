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

    # get the SOAPdenovo2 config file path
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # create the SOAPdenovo2 config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(soapdenovo2_config_file)):
            os.makedirs(os.path.dirname(soapdenovo2_config_file))
        with open(soapdenovo2_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of SOAPdenovo2 and their meaning in http://soap.genomics.org.cn/SOAPdenovo2.html.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# There are two formats to set an option:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value                             <- if the option supports a single value'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value-1, value-2, ..., value-n    <- if the option supports a values list'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# WARNING: The files have to be decompressed.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the SOAPdenovo2 parameters'))
            file_id.write( '{0}\n'.format('[SOAPdenovo2 parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('version = 63', '# SOAPdenovo2 version: {0}'.format(get_version_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of cpu for use'))
            file_id.write( '{0:<50} {1}\n'.format('init_memory_assumption = 0', '# memory assumption initialized to avoid further reallocation in GiB'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = 25', '# value or values list of kmer size: minimum, 13; maximum: version value.'))
            file_id.write( '{0:<50} {1}\n'.format('kmer_freq_cutoff = 0', '# kmers with frequency no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format('edge_cov_cutoff = 2', '# edges with coverage no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format('resolve_repeats = NO', '# resolve repeats by reads: {0}'.format(get_resolve_repeats_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('merge_level = 1', '# strength of merging similar sequences during contiging: minimum, 0; maximum, 3.'))
            file_id.write( '{0:<50} {1}\n'.format('filter = 0', '# weight to filter arc when linearize two edges'))
            file_id.write( '{0:<50} {1}\n'.format('mapping_kmer = 25', '# kmer size used for mapping read to contig: minimum, 13; maximum: 63.'))
            file_id.write( '{0:<50} {1}\n'.format('keep_read = NO', '# keep available read(*.read): {0}'.format(get_keep_read_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('merge_bubble = NO', '# merge clean bubble before iterate: {0}'.format(get_merge_bubble_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('srkgf = NO', '# output gap related redas for SRkgf to fill gap: {0}'.format(get_srkgf_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('fill = NO', '# fill gaps in scaffolds: {0}'.format(get_fill_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('unmask = NO', '# un-mask contigs with high/low coverage before scaffolding: {0}'.format(get_unmask_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('keep_contigs = NO', '# keep contigs weakly connected to other contigs in scaffold: {0}'.format(get_keep_contigs_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('gap_len_diff = 50', '# allowed length difference between estimated and filled gap'))
            file_id.write( '{0:<50} {1}\n'.format('min_contig_len = 27', '# shortest contig for scaffolding'))
            file_id.write( '{0:<50} {1}\n'.format('min_contig_cvg = 0.1', '# minimum contig coverage'))
            file_id.write( '{0:<50} {1}\n'.format('max_contig_cvg = 2.0', '# maximum contig coverage'))
            file_id.write( '{0:<50} {1}\n'.format('insert_size_upper_bound = 1.5', '# upper bound of insert size for large insert size'))
            file_id.write( '{0:<50} {1}\n'.format('bubble_coverage = 0.6', '# lower coverage in bubble structure'))
            file_id.write( '{0:<50} {1}\n'.format('genome_size = 0', '# genome size for statistics'))
            file_id.write( '{0:<50} {1}\n'.format('visualization = NO', '# output visualization information of assembly: {0}'.format(get_visualization_code_list_text())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('max_rd_len = 100', '# any read longer than the value will be cut to this length'))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first library.'))
                file_id.write( '{0}\n'.format('[library-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# read type: {0}'.format(get_read_type_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('read_file_1 = {0}'.format(os.path.basename(file_1_list[i])), '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = {0}'.format(os.path.basename(file_2_list[i])), '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                file_id.write( '{0:<50} {1}\n'.format('avg_ins = 200', '# average insert size'))
                file_id.write( '{0:<50} {1}\n'.format('reverse_seq = 0', '# sequence direction: {0}'.format(get_reverse_seq_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('asm_flags = 3', '# in which part(s) the reads are used: {0}'.format(get_asm_flags_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('rd_len_cutof = 100', '# reads will be cut to this length'))
                file_id.write( '{0:<50} {1}\n'.format('pair_num_cutoff = 3', '# cutoff value of pair number for a reliable connection between two contigs or pre-scaffolds (minimum: 3)'))
                file_id.write( '{0:<50} {1}\n'.format('map_len = 32', '# minimum aligned length between a read and a contigs for a reliable read location (minimum: 32)'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more libraries, you have to repeat the section library-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovo2_config_file))
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

    # get the SOAPdenovo2 code and name
    soapdenovo2_code = xlib.get_soapdenovo2_code()
    soapdenovo2_name = xlib.get_soapdenovo2_name()

    # get the SOAPdenovo2 config file
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # get the SOAPdenovo2 option dictionary
    soapdenovo2_option_dict = xlib.get_option_dict(soapdenovo2_config_file)

    # get the experiment identification
    experiment_id = soapdenovo2_option_dict['identification']['experiment_id']

    # get the SOAPdenovo2 process script path in the local computer
    soapdenovo2_process_script = get_soapdenovo2_process_script()

    # get the SOAPdenovo2 process starter path in the local computer
    soapdenovo2_process_starter = get_soapdenovo2_process_starter()

    # get the SOAPdenovo2 process config file path in the local computer
    soapdenovo2_process_config_file = get_soapdenovo2_process_config_file()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the SOAPdenovo2 configuration file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(soapdenovo2_name))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the SOAPdenovo2 is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_soapdenovo2_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(soapdenovo2_name))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(soapdenovo2_name))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # build the process configuration file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process config file {0} ...\n'.format(soapdenovo2_process_config_file))
        (OK, error_list) = build_soapdenovo2_process_config_file()
        if OK:
            log.write('The file is built.\n')
        else:
            log.write('*** ERROR: The file could not be built.\n')

    # for each kmer value, build the process, copy it the cluster and run it
    if OK:

        # get the kmer list
        kmer = soapdenovo2_option_dict['SOAPdenovo2 parameters']['kmer']
        kmer_list = xlib.split_literal_to_integer_list(kmer)
        
        # for each kmer value, do the tasks
        i = 1
        for kmer_value in kmer_list:

            # determine the run directory in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Determining the run directory for kmer {0} in the cluster is being determined ...\n'.format(kmer_value))
            if i > 1:
                current_run_dir = '{0}-{1}'.format(xlib.get_cluster_current_run_dir(experiment_id, soapdenovo2_code), i)
            else:
                current_run_dir = '{0}'.format(xlib.get_cluster_current_run_dir(experiment_id, soapdenovo2_code))
            command = f'mkdir --parents {current_run_dir}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The directory path is {0}.\n'.format(current_run_dir))
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            i += 1

            # build the SOAPdenovo2 process script
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process script {0} ...\n'.format(soapdenovo2_process_script))
            (OK, error_list) = build_soapdenovo2_process_script(cluster_name, current_run_dir, kmer_value)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')
                break

            # upload the process configuration file to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process config file {0} to the directory {1} of the master ...\n'.format(soapdenovo2_process_config_file, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovo2_process_config_file))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovo2_process_config_file, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # upload the process script to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(soapdenovo2_process_script, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovo2_process_script))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovo2_process_script, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process script in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovo2_process_script)))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(soapdenovo2_process_script))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # build the process starter
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process starter {0} ...\n'.format(soapdenovo2_process_starter))
            (OK, error_list) = build_soapdenovo2_process_starter(current_run_dir)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('***ERROR: The file could not be built.\n')
                break

            # upload the process starter to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(soapdenovo2_process_starter, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovo2_process_starter))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovo2_process_starter, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process starter in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovo2_process_starter)))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(soapdenovo2_process_starter))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # submit the process
            log.write(f'{xlib.get_separator()}\n')
            log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovo2_process_starter)))
            OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(soapdenovo2_process_starter), log)

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

    # get the SOAPdenovo2 name
    soapdenovo2_name = xlib.get_soapdenovo2_name()

    # get the SOAPdenovo2 configuration file path
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # get the option dictionary
    try:
        soapdenovo2_option_dict = xlib.get_option_dict(soapdenovo2_config_file)
    except Exception as e:
        error_list.append('*** ERROR: The dictionary with options could  not be retrieved from the config file.')
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
                error_list.append('*** ERROR: the key "version" has to be {0}.'.format(get_version_code_list_text()))
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
                error_list.append('*** ERROR: the key "resolve_repeats" has to be {0}.'.format(get_resolve_repeats_code_list_text()))
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

            # check section "SOAPdenovo2 parameters" - key "mapping_kmer"
            mapping_kmer = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('mapping_kmer', not_found)
            if mapping_kmer == not_found:
                error_list.append('*** ERROR: the key "mapping_kmer" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_int(mapping_kmer, minimum=13, maximum=63):
                error_list.append('*** ERROR: the key "mapping_kmer" has to be an integer number between 13 and 63.')
                OK = False

            # check section "SOAPdenovo2 parameters" - key "keep_read"
            keep_read = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('keep_read', not_found)
            if keep_read == not_found:
                error_list.append('*** ERROR: the key "keep_read" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(keep_read, get_keep_read_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "keep_read" has to be {0}.'.format(get_keep_read_code_list_text()))
                OK = False

            # check section "SOAPdenovo2 parameters" - key "merge_bubble"
            merge_bubble = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('merge_bubble', not_found)
            if merge_bubble == not_found:
                error_list.append('*** ERROR: the key "merge_bubble" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(merge_bubble, get_merge_bubble_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "merge_bubble" has to be {0}.'.format(get_merge_bubble_code_list_text()))
                OK = False

            # check section "SOAPdenovo2 parameters" - key "srkgf"
            srkgf = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('srkgf', not_found)
            if srkgf == not_found:
                error_list.append('*** ERROR: the key "srkgf" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(srkgf, get_srkgf_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "srkgf" has to be {0}.'.format(get_srkgf_code_list_text()))
                OK = False

            # check section "SOAPdenovo2 parameters" - key "fill"
            fill = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('fill', not_found)
            if fill == not_found:
                error_list.append('*** ERROR: the key "fill" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(fill, get_fill_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "fill" has to be {0}.'.format(get_fill_code_list_text()))
                OK = False

            # check section "SOAPdenovo2 parameters" - key "unmask"
            unmask = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('unmask', not_found)
            if unmask == not_found:
                error_list.append('*** ERROR: the key "unmask" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(unmask, get_unmask_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "unmask" has to be {0}.'.format(get_unmask_code_list_text()))
                OK = False

            # check section "SOAPdenovo2 parameters" - key "keep_contigs"
            keep_contigs = soapdenovo2_option_dict.get('SOAPdenovo2 parameters', {}).get('keep_contigs', not_found)
            if keep_contigs == not_found:
                error_list.append('*** ERROR: the key "keep_contigs" is not found in the section "SOAPdenovo2 parameters".')
                OK = False
            elif not xlib.check_code(keep_contigs, get_keep_contigs_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "keep_contigs" has to be {0}.'.format(get_keep_contigs_code_list_text()))
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
                error_list.append('*** ERROR: the key "visualization" has to be {0}.'.format(get_visualization_code_list_text()))
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
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "format"
                    format = soapdenovo2_option_dict.get(section, {}).get('format', not_found)
                    if format == not_found:
                        error_list.append('*** ERROR: the key "format" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                        OK = False

                    # check section "library-n" - key "read_type"
                    read_type = soapdenovo2_option_dict.get(section, {}).get('read_type', not_found)
                    if read_type == not_found:
                        error_list.append('*** ERROR: the key "read_type" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "read_type" has to be {0}.'.format(get_read_type_code_list_text()))
                        OK = False
                    elif read_type == 'SP' and format != 'FASTA':
                        error_list.append('*** ERROR: if read type is SP, the format has to be FASTA.')
                        OK = False

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = soapdenovo2_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False
                    elif read_file_1.find('.gz') != -1:
                        error_list.append('*** ERROR: the key "read_file_1" in the section "{0}" has to be a decompressed file.'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = soapdenovo2_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False
                    elif read_file_2.find('.gz') != -1:
                        error_list.append('*** ERROR: the key "read_file_2" in the section "{0}" has to be a decompressed file.'.format(section))
                        OK = False

                    # check section "library" - key "avg_ins"
                    avg_ins = soapdenovo2_option_dict.get(section, {}).get('avg_ins', not_found)
                    if avg_ins == not_found:
                        error_list.append('*** ERROR: the key "avg_ins" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(avg_ins, minimum=0):
                        error_list.append('*** ERROR: the key "avg_ins" in the section "{0}" has to be an integer number greater than or equal to 0.'.format(section))
                        OK = False

                    # check section "library-n" - key "reverse_seq"
                    reverse_seq = soapdenovo2_option_dict.get(section, {}).get('reverse_seq', not_found)
                    if reverse_seq == not_found:
                        error_list.append('*** ERROR: the key "reverse_seq" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(reverse_seq, get_reverse_seq_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "reverse_seq" has to be {0}.'.format(get_reverse_seq_code_list_text()))
                        OK = False

                    # check section "library-n" - key "asm_flags"
                    asm_flags = soapdenovo2_option_dict.get(section, {}).get('asm_flags', not_found)
                    if asm_flags == not_found:
                        error_list.append('*** ERROR: the key "asm_flags" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(asm_flags, get_asm_flags_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "asm_flags" has to be {0}.'.format(get_asm_flags_code_list_text()))
                        OK = False

                    # check section "library-n" - key "rd_len_cutof"
                    rd_len_cutof = soapdenovo2_option_dict.get(section, {}).get('rd_len_cutof', not_found)
                    if rd_len_cutof == not_found:
                        error_list.append('*** ERROR: the key "rd_len_cutof" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(rd_len_cutof, minimum=0):
                        error_list.append('*** ERROR: the key "rd_len_cutof" in the section "{0}" has to be an integer number greater than or equal to 1.'.format(section))
                        OK = False

                    # check section "library-n" - key "pair_num_cutoff"
                    pair_num_cutoff = soapdenovo2_option_dict.get(section, {}).get('pair_num_cutoff', not_found)
                    if pair_num_cutoff == not_found:
                        error_list.append('*** ERROR: the key "pair_num_cutoff" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(pair_num_cutoff, minimum=0):
                        error_list.append('*** ERROR: the key "pair_num_cutoff" in the section "{0}" has to be an integer number greater than or equal to 1.'.format(section))
                        OK = False

                    # check section "library-n" - key "map_len"
                    map_len = soapdenovo2_option_dict.get(section, {}).get('map_len', not_found)
                    if map_len == not_found:
                        error_list.append('*** ERROR: the key "map_len" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(map_len, minimum=32):
                        error_list.append('*** ERROR: the key "map_len" in the section "{0}" has to be an integer number greater than or equal to 32.'.format(section))
                        OK = False

    # warn that the SOAPdenovo2 config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(soapdenovo2_name))

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

    # get the SOAPdenovo2 config file path
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # get the option dictionary
    soapdenovo2_options_dict = xlib.get_option_dict(soapdenovo2_config_file)

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
    mapping_kmer = soapdenovo2_options_dict['SOAPdenovo2 parameters']['mapping_kmer']
    keep_read = soapdenovo2_options_dict['SOAPdenovo2 parameters']['keep_read'].upper()
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

    # get the SOAPdenovo2 process config file name
    soapdenovo2_process_config_file = get_soapdenovo2_process_config_file()

    # get the SOAPdenovo2 process script name
    soapdenovo2_process_script = get_soapdenovo2_process_script()

    # write the SOAPdenovo2 process script
    try:
        if not os.path.exists(os.path.dirname(soapdenovo2_process_script)):
            os.makedirs(os.path.dirname(soapdenovo2_process_script))
        with open(soapdenovo2_process_script, mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('ulimit -s unlimited'))
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('SOAPDENOVO2_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_soapdenovo2_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$SOAPDENOVO2_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_soapdenovo2_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_soapdenovo2_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        SOAPdenovo-{0}mer all \\'.format(version)))
            script_file_id.write( '{0}\n'.format('            -s {0}/{1} \\'.format(current_run_dir, os.path.basename(soapdenovo2_process_config_file))))
            script_file_id.write( '{0}\n'.format('            -o {0}-{1} \\'.format(experiment_id, os.path.basename(current_run_dir))))
            script_file_id.write( '{0}\n'.format('            -p {0} \\'.format(ncpu)))
            script_file_id.write( '{0}\n'.format('            -a {0} \\'.format(init_memory_assumption)))
            script_file_id.write( '{0}\n'.format('            -K {0} \\'.format(kmer_value)))
            script_file_id.write( '{0}\n'.format('            -d {0} \\'.format(kmer_freq_cutoff)))
            script_file_id.write( '{0}\n'.format('            -D {0} \\'.format(edge_cov_cutoff)))
            if resolve_repeats == 'YES':
                script_file_id.write( '{0}\n'.format('            -R \\'))

            script_file_id.write( '{0}\n'.format('            -M {0} \\'.format(merge_level)))
            script_file_id.write( '{0}\n'.format('            -e {0} \\'.format(filter)))
            script_file_id.write( '{0}\n'.format('            -k {0} \\'.format(mapping_kmer)))
            if keep_read == 'YES':
                script_file_id.write( '{0}\n'.format('            -r \\'))
            if merge_bubble == 'YES':
                script_file_id.write( '{0}\n'.format('            -E \\'))
            if srkgf == 'YES':
                script_file_id.write( '{0}\n'.format('            -f \\'))
            if fill == 'YES':
                script_file_id.write( '{0}\n'.format('            -F \\'))
            if unmask == 'YES':
                script_file_id.write( '{0}\n'.format('            -u \\'))
            if keep_contigs == 'YES':
                script_file_id.write( '{0}\n'.format('            -w \\'))
            if visualization == 'YES':
                script_file_id.write( '{0}\n'.format('            -V {0} \\'.format(visualization)))
            script_file_id.write( '{0}\n'.format('            -G {0} \\'.format(gap_len_diff)))
            script_file_id.write( '{0}\n'.format('            -L {0} \\'.format(min_contig_len)))
            script_file_id.write( '{0}\n'.format('            -c {0} \\'.format(min_contig_cvg)))
            script_file_id.write( '{0}\n'.format('            -C {0} \\'.format(max_contig_cvg)))
            script_file_id.write( '{0}\n'.format('            -b {0} \\'.format(insert_size_upper_bound)))
            script_file_id.write( '{0}\n'.format('            -B {0} \\'.format(bubble_coverage)))
            script_file_id.write( '{0}\n'.format('            -N {0}'.format(genome_size)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error SOAPdenovo-{0}mer $RC; fi'.format(version)))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_soapdenovo2_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_soapdenovo2_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_soapdenovo2_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_soapdenovo2_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_soapdenovo2_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovo2_process_script))
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

    # get the SOAPdenovo2 process starter path in local
    soapdenovo2_process_starter = get_soapdenovo2_process_starter()

    # get the SOAPdenovo2 process script path in local
    soapdenovo2_process_script = get_soapdenovo2_process_script()

    # get the log file name
    log_file = xlib.get_cluster_log_file()

    # write the SOAPdenovo2 process starter
    try:
        if not os.path.exists(os.path.dirname(soapdenovo2_process_starter)):
            os.makedirs(os.path.dirname(soapdenovo2_process_starter))
        with open(soapdenovo2_process_starter, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(soapdenovo2_process_script), log_file)))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovo2_process_starter))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovo2_process_config_file():
    '''
    Build the SOAPdenovo2 process config file to the current SOPAdenovo2 experiment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the SOAPdenovo2 config file path
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # get the SOAPdenovo2 process config file path
    soapdenovo2_process_config_file = get_soapdenovo2_process_config_file()

    # get the option dictionary
    soapdenovo2_option_dict = xlib.get_option_dict(soapdenovo2_config_file)

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
        if not os.path.exists(os.path.dirname(soapdenovo2_process_config_file)):
            os.makedirs(os.path.dirname(soapdenovo2_process_config_file))
        with open(soapdenovo2_process_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#maximal read length'))
            file_id.write( '{0}\n'.format('max_rd_len={0}'.format(soapdenovo2_option_dict['library']['max_rd_len'])))
            for section in sections_list:
                if re.match('^library-[0-9]+$', section):
                    format = soapdenovo2_option_dict[section]['format'].upper()
                    read_type = soapdenovo2_option_dict[section]['read_type'].upper()
                    read_file_1 = soapdenovo2_option_dict[section]['read_file_1']
                    read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
                    if read_type == 'PE':
                        read_file_2 = soapdenovo2_option_dict[section]['read_file_2']
                        read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                    file_id.write( '{0}\n'.format('[LIB]'))
                    file_id.write( '{0}\n'.format('# average insert size'))
                    file_id.write( '{0}\n'.format('avg_ins={0}'.format(soapdenovo2_option_dict[section]['avg_ins'])))
                    file_id.write( '{0}\n'.format('# if sequence needs to be reversed'))
                    file_id.write( '{0}\n'.format('reverse_seq={0}'.format(soapdenovo2_option_dict[section]['reverse_seq'])))
                    file_id.write( '{0}\n'.format('# in which part(s) the reads are used'))
                    file_id.write( '{0}\n'.format('asm_flags={0}'.format(soapdenovo2_option_dict[section]['asm_flags'])))
                    file_id.write( '{0}\n'.format('# maximal read length in this lib'))
                    file_id.write( '{0}\n'.format('rd_len_cutof={0}'.format(soapdenovo2_option_dict[section]['rd_len_cutof'])))
                    file_id.write( '{0}\n'.format('# cutoff of pair number for a reliable connection (at least 3 for short insert size)'))
                    file_id.write( '{0}\n'.format('pair_num_cutoff={0}'.format(soapdenovo2_option_dict[section]['pair_num_cutoff'])))
                    file_id.write( '{0}\n'.format('# minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)'))
                    file_id.write( '{0}\n'.format('map_len={0}'.format(soapdenovo2_option_dict[section]['map_len'])))
                    if format == 'FASTA' and read_type == 'SE':
                        file_id.write( '{0}\n'.format('# fasta file for single reads'))
                        file_id.write( '{0}\n'.format('f={0}'.format(read_file_1)))
                    elif format == 'FASTA' and read_type == 'PE':
                        file_id.write( '{0}\n'.format('# fasta file for read 1'))
                        file_id.write( '{0}\n'.format('f1={0}'.format(read_file_1)))
                        file_id.write( '{0}\n'.format('# fastq file for read 2 always follows fastq file for read 1'))
                        file_id.write( '{0}\n'.format('f2={0}'.format(read_file_2)))
                    elif format == 'FASTA' and read_type == 'SP':
                        file_id.write( '{0}\n'.format('# a single fasta file for paired reads'))
                        file_id.write( '{0}\n'.format('p={0}'.format(read_file_1)))
                    elif format == 'FASTQ' and read_type == 'SE':
                        file_id.write( '{0}\n'.format('# fastq file for single reads'))
                        file_id.write( '{0}\n'.format('q={0}'.format(read_file_1)))
                    elif format == 'FASTQ' and read_type == 'PE':
                        file_id.write( '{0}\n'.format('# fastq file for read 1'))
                        file_id.write( '{0}\n'.format('q1={0}'.format(read_file_1)))
                        file_id.write( '{0}\n'.format('# fastq file for read 2 always follows fastq file for read 1'))
                        file_id.write( '{0}\n'.format('q2={0}'.format(read_file_2)))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovo2_process_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def determine_soapdenovo2_cluster():
    '''
    Determine the cluster to the current SOAPdenovo2 experiment.
    '''

    # initialize the template and cluster names
    template_name = ''
    cluster_name = ''

    # get the SOAPdenovo2 config file path
    soapdenovo2_config_file = get_soapdenovo2_config_file()

    # get the option dictionary
    soapdenovo2_options_dict = xlib.get_option_dict(soapdenovo2_config_file)

    # determine the template and cluster names
    template_name = soapdenovo2_options_dict['performance']['template_name']
    if template_name == 'AUTO': 
        template_name = 'cl-t1.micro'
        cluster_name = template_name
    else: 
        cluster_name = template_name

    # return the template and cluster names
    return (template_name, cluster_name)

#-------------------------------------------------------------------------------

def get_soapdenovo2_config_file():
    '''
    Get the SOAPdenovo2 config file path.
    '''

    # assign the SOAPdenovo2 config file path
    soapdenovo2_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_soapdenovo2_code())

    # return the SOAPdenovo2 config file path
    return soapdenovo2_config_file

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_script():
    '''
    Get the SOAPdenovo2 process script path in the local computer.
    '''

    # assign the SOAPdenovo2 script path
    soapdenovo2_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_soapdenovo2_code())

    # return the SOAPdenovo2 script path
    return soapdenovo2_process_script

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_starter():
    '''
    Get the SOAPdenovo2 process starter path in the local computer.
    '''

    # assign the SOAPdenovo2 process starter path
    soapdenovo2_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_soapdenovo2_code())

    # return the SOAPdenovo2 starter path
    return soapdenovo2_process_starter

#-------------------------------------------------------------------------------

def get_soapdenovo2_process_config_file():
    '''
    Get the SOAPdenovo2 process config file path in the local computer.
    '''

    # assign the SOAPdenovo2 process config file path
    soapdenovo2_process_config_file = '{0}/{1}-process-config.txt'.format(xlib.get_temp_dir(), xlib.get_soapdenovo2_code())

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
    
def get_keep_read_code_list():
    '''
    Get the code list of "keep_read".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_keep_read_code_list_text():
    '''
    Get the code list of "keep_read" as text.
    '''

    return str(get_keep_read_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
