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
This file contains functions related to the HISAT2 process used in both console
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

def create_hisat2_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', gtf_file='Arabidopsis_thaliana.TAIR10.36.gtf', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create HISAT2 config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the HISAT2 config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_hisat2_config_file())):
            os.makedirs(os.path.dirname(get_hisat2_config_file()))
        with open(get_hisat2_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference and GTF files have to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n'.format(''.format()))
            file_id.write( '# The experiment_id, reference_dataset_id, reference_file, assembly_dataset_id and read_dataset_id are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of HISAT2 and their meaning in "https://ccb.jhu.edu/software/hisat2/".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "HISAT2 parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of HISAT2 and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --new-summary; --score-min=L,0,-0.2\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format(f'gtf_file = {gtf_file}', '# GTF file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the HISAT2 parameters\n')
            file_id.write( '[HISAT2 parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'index_building = YES', f'# index building: {get_index_building_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'large_index = YES', f'# a large index is force, even if the reference is less than ~ 4 billion nucleotides long: {get_large_index_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'dta_cufflinks = YES', f'# alignments tailored specifically for Cufflinks: {get_dta_cufflinks_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'min_mp = 2', '# minimum mismatch penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_mp = 6', '# maximum mismatch penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'no_softclip = NO', f'# disallow soft-clipping: {get_no_softclip_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'min_sp = 1', '# minimum penalty for soft-clipping per base'))
            file_id.write( '{0:<50} {1}\n'.format( 'max_sp = 2', '# maximum penalty for soft-clipping per base'))
            file_id.write( '{0:<50} {1}\n'.format( 'np = 1', '# penalty for positions where the read, reference, or both, contain an ambiguous character such as N'))
            file_id.write( '{0:<50} {1}\n'.format( 'open_rdg = 5', '# read gap open penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'extend_rdg = 3', '# read gap extend penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'open_rfg = 5', '# reference gap open penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'extend_rfg = 3', '# reference gap extend penalty'))
            file_id.write( '{0:<50} {1}\n'.format( 'orientation = FR', f'# orientation of paired-end reads: {get_orientation_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'quality-score = 33', f'# FASTQ quality score: {get_quality_score_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'format = FASTQ', f'# format: {get_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_type = {read_type}', f'# read type: {get_read_type_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'library_concatenation = NO', f'# {get_library_concatenation_code_list_text()}'))
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
        error_list.append(f'*** ERROR: The file {get_hisat2_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_hisat2_process(cluster_name, log, function=None):
    '''
    Run a HISAT2 process.
    '''

    # initialize the control variable
    OK = True

    # get the HISAT2 option dictionary
    hisat2_option_dict = xlib.get_option_dict(get_hisat2_config_file())

    # get the experiment identification
    experiment_id = hisat2_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the HISAT2 config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_hisat2_name()} config file ...\n')
    (OK, error_list) = check_hisat2_config_file(strict=True)
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

    # check the HISAT2 is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_hisat2_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_hisat2_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_hisat2_name()} installation could not be performed.\n')

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
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_hisat2_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the HISAT2 process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_hisat2_process_script()} ...\n')
        (OK, error_list) = build_hisat2_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the HISAT2 process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_hisat2_process_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_hisat2_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_hisat2_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the HISAT2 process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_hisat2_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_hisat2_process_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the HISAT2 process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_hisat2_process_starter()} ...\n')
        (OK, error_list) = build_hisat2_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the HISAT2 process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_hisat2_process_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_hisat2_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_hisat2_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the HISAT2 process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_hisat2_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_hisat2_process_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the HISAT2 process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_hisat2_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_hisat2_process_starter()), log)

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

def check_hisat2_config_file(strict):
    '''
    Check the HISAT2 config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        hisat2_option_dict = xlib.get_option_dict(get_hisat2_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in hisat2_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = hisat2_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = hisat2_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = hisat2_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "gtf_file"
            gtf_file = hisat2_option_dict.get('identification', {}).get('gtf_file', not_found)
            if gtf_file == not_found:
                error_list.append('*** ERROR: the key "gtf_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = hisat2_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "HISAT2 parameters"
        if 'HISAT2 parameters' not in sections_list:
            error_list.append('*** ERROR: the section "HISAT2 parameters" is not found.')
            OK = False
        else:

            # check section "HISAT2 parameters" - key "index_building"
            index_building = hisat2_option_dict.get('HISAT2 parameters', {}).get('index_building', not_found)
            if index_building == not_found:
                error_list.append('*** ERROR: the key "index_building" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(index_building, get_index_building_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "index_building" has to be {get_index_building_code_list_text()}.')
                OK = False

            # check section "HISAT2 parameters" - key "large_index"
            large_index = hisat2_option_dict.get('HISAT2 parameters', {}).get('large_index', not_found)
            if large_index == not_found:
                error_list.append('*** ERROR: the key "large_index" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(large_index, get_large_index_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "large_index" has to be {get_large_index_code_list_text()}.')
                OK = False

            # check section "HISAT2 parameters" - key "threads"
            threads = hisat2_option_dict.get('HISAT2 parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "HISAT2 parameters" - key "dta_cufflinks"
            dta_cufflinks = hisat2_option_dict.get('HISAT2 parameters', {}).get('dta_cufflinks', not_found)
            if dta_cufflinks == not_found:
                error_list.append('*** ERROR: the key "dta_cufflinks" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(dta_cufflinks, get_dta_cufflinks_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "dta_cufflinks" has to be {get_dta_cufflinks_code_list_text()}.')
                OK = False

            # check section "HISAT2 parameters" - key "min_mp"
            min_mp = hisat2_option_dict.get('HISAT2 parameters', {}).get('min_mp', not_found)
            is_ok_min_mp = False
            if min_mp == not_found:
                error_list.append('*** ERROR: the key "min_mp" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(min_mp, minimum=0):
                error_list.append('*** ERROR: the key "min_mp" has to be an integer number greater than or equal to 0.')
                OK = False
            else:
                is_ok_min_mp = True

            # check section "HISAT2 parameters" - key "max_mp"
            max_mp = hisat2_option_dict.get('HISAT2 parameters', {}).get('max_mp', not_found)
            is_ok_max_mp = False
            if max_mp == not_found:
                error_list.append('*** ERROR: the key "max_mp" is not found in the section "HISAT2 parameters".')
                OK = False
                max_mp = 99999999999999
            elif not xlib.check_int(max_mp, minimum=0):
                error_list.append('*** ERROR: the key "max_mp" has to be an integer number greater than or equal to 0.')
                OK = False
            else:
                is_ok_max_mp = True

            # check if max_mp value is greater than or equal than min_mp value
            if is_ok_min_mp and is_ok_max_mp and int(max_mp) < int(min_mp):
                error_list.append(f'*** ERROR: The value max_mp value ({max_mp}) is less than the min_mp value ({min_mp}).')
                OK = False

            # check section "HISAT2 parameters" - key "no_softclip"
            no_softclip = hisat2_option_dict.get('HISAT2 parameters', {}).get('no_softclip', not_found)
            if no_softclip == not_found:
                error_list.append('*** ERROR: the key "no_softclip" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(no_softclip, get_no_softclip_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "no_softclip" has to be {get_no_softclip_code_list_text()}.')
                OK = False

            # check section "HISAT2 parameters" - key "min_sp"
            min_sp = hisat2_option_dict.get('HISAT2 parameters', {}).get('min_sp', not_found)
            is_ok_min_sp = False
            if min_sp == not_found:
                error_list.append('*** ERROR: the key "min_sp" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(min_sp, minimum=0):
                error_list.append('*** ERROR: the key "min_sp" has to be an integer number greater than or equal to 0.')
                OK = False
            else:
                is_ok_min_sp = True

            # check section "HISAT2 parameters" - key "max_sp"
            max_sp = hisat2_option_dict.get('HISAT2 parameters', {}).get('max_sp', not_found)
            is_ok_max_sp = False
            if max_sp == not_found:
                error_list.append('*** ERROR: the key "max_sp" is not found in the section "HISAT2 parameters".')
                OK = False
                max_sp = 99999999999999
            elif not xlib.check_int(max_sp, minimum=0):
                error_list.append('*** ERROR: the key "max_sp" has to be an integer number greater than or equal to 0.')
                OK = False
            else:
                is_ok_max_sp = True

            # check if max_sp value is greater than or equal than min_sp value
            if is_ok_min_sp and is_ok_max_sp and int(max_sp) < int(min_sp):
                error_list.append(f'*** ERROR: The value max_sp value ({max_sp}) is less than the min_sp value ({min_sp}).')
                OK = False

            # check section "HISAT2 parameters" - key "np"
            np = hisat2_option_dict.get('HISAT2 parameters', {}).get('np', not_found)
            if np == not_found:
                error_list.append('*** ERROR: the key "np" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(np, minimum=0):
                error_list.append('*** ERROR: the key "np" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "HISAT2 parameters" - key "open_rdg"
            open_rdg = hisat2_option_dict.get('HISAT2 parameters', {}).get('open_rdg', not_found)
            if open_rdg == not_found:
                error_list.append('*** ERROR: the key "open_rdg" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(open_rdg, minimum=0):
                error_list.append('*** ERROR: the key "open_rdg" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "HISAT2 parameters" - key "extend_rdg"
            extend_rdg = hisat2_option_dict.get('HISAT2 parameters', {}).get('extend_rdg', not_found)
            if extend_rdg == not_found:
                error_list.append('*** ERROR: the key "extend_rdg" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(extend_rdg, minimum=0):
                error_list.append('*** ERROR: the key "extend_rdg" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "HISAT2 parameters" - key "open_rfg"
            open_rfg = hisat2_option_dict.get('HISAT2 parameters', {}).get('open_rfg', not_found)
            if open_rfg == not_found:
                error_list.append('*** ERROR: the key "open_rfg" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(open_rfg, minimum=0):
                error_list.append('*** ERROR: the key "open_rfg" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "HISAT2 parameters" - key "extend_rfg"
            extend_rfg = hisat2_option_dict.get('HISAT2 parameters', {}).get('extend_rfg', not_found)
            if extend_rfg == not_found:
                error_list.append('*** ERROR: the key "extend_rfg" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(extend_rfg, minimum=0):
                error_list.append('*** ERROR: the key "extend_rfg" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "HISAT2 parameters" - key "orientation"
            orientation = hisat2_option_dict.get('HISAT2 parameters', {}).get('orientation', not_found)
            if orientation == not_found:
                error_list.append('*** ERROR: the key "orientation" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(orientation, get_orientation_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "orientation" has to be {get_orientation_code_list_text()}.')
                OK = False

            # check section "HISAT2 parameters" - key "quality-zero-score"
            quality_score = hisat2_option_dict.get('HISAT2 parameters', {}).get('quality-score', not_found)
            if quality_score == not_found:
                error_list.append('*** ERROR: the key "quality-score" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(quality_score):
                error_list.append('*** ERROR: the key "quality-score" has to be an integer number.')
                OK = False

            # check section "HISAT2 parameters" - key "other_parameters"
            not_allowed_parameters_list = ['threads', 'dta-cufflinks', 'qseq', 'phred33', 'phred64', 'mp', 'sp', 'no-softclip', 'np', 'rdg', 'rfg', 'known-splicesite-infile', 'time', 'un', 'un-gz', 'un-bz2', 'al', 'al-gz', 'al-bz2', 'un-conc', 'un-conc-gz', 'un-conc-bz2', 'al-conc', 'al-conc-gz', 'al-conc-bz2', 'quiet', 'summary-file']
            other_parameters = hisat2_option_dict.get('HISAT2 parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "HISAT2 parameters".')
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
            format = hisat2_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = hisat2_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "read_type" has to be {get_read_type_code_list_text()}.')
                OK = False

            # check section "library" - key "library_concatenation"
            library_concatenation = hisat2_option_dict.get('library', {}).get('library_concatenation', not_found)
            if library_concatenation == not_found:
                error_list.append('*** ERROR: the key "library_concatenation" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(library_concatenation, get_library_concatenation_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "library_concatenation" has to be {get_library_concatenation_code_list_text()}.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'HISAT2 parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = hisat2_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = hisat2_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_hisat2_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_hisat2_process_script(cluster_name, current_run_dir):
    '''
    Build the current HISAT2 process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the HISAT2 option dictionary
    hisat2_option_dict = xlib.get_option_dict(get_hisat2_config_file())

    # get the options
    experiment_id = hisat2_option_dict['identification']['experiment_id']
    reference_dataset_id = hisat2_option_dict['identification']['reference_dataset_id']
    reference_file = hisat2_option_dict['identification']['reference_file']
    gtf_file = hisat2_option_dict['identification']['gtf_file']
    read_dataset_id = hisat2_option_dict['identification']['read_dataset_id']
    index_building = hisat2_option_dict['HISAT2 parameters']['index_building']
    large_index = hisat2_option_dict['HISAT2 parameters']['large_index']
    threads = hisat2_option_dict['HISAT2 parameters']['threads']
    dta_cufflinks = hisat2_option_dict['HISAT2 parameters']['dta_cufflinks']
    min_mp = hisat2_option_dict['HISAT2 parameters']['min_mp']
    max_mp = hisat2_option_dict['HISAT2 parameters']['max_mp']
    no_softclip = hisat2_option_dict['HISAT2 parameters']['no_softclip']
    min_sp = hisat2_option_dict['HISAT2 parameters']['min_sp']
    max_sp = hisat2_option_dict['HISAT2 parameters']['max_sp']
    np = hisat2_option_dict['HISAT2 parameters']['np']
    open_rdg = hisat2_option_dict['HISAT2 parameters']['open_rdg']
    extend_rdg = hisat2_option_dict['HISAT2 parameters']['extend_rdg']
    open_rfg = hisat2_option_dict['HISAT2 parameters']['open_rfg']
    extend_rfg = hisat2_option_dict['HISAT2 parameters']['extend_rfg']
    orientation = hisat2_option_dict['HISAT2 parameters']['orientation']
    other_parameters = hisat2_option_dict['HISAT2 parameters']['other_parameters']
    format = hisat2_option_dict['library']['format']
    read_type = hisat2_option_dict['library']['read_type']
    library_concatenation = hisat2_option_dict['library']['library_concatenation']

    # get the sections list
    sections_list = []
    for section in hisat2_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build read file lists
    read_file_1_list = []
    read_file_2_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = hisat2_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            read_file_1_list.append(read_file_1)
            if read_type.upper() == 'PE':
                read_file_2 = hisat2_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                read_file_2_list.append(read_file_2)
    if library_concatenation.upper() == 'YES':
        read_file_1_list = [",".join(read_file_1_list)]
        if read_type.upper() == 'PE':
            read_file_2_list = [",".join(read_file_2_list)]

    # set the cluster reference dataset directory
    cluster_reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)

    # set the cluster reference file
    cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

    # set the cluster GTF file
    if gtf_file.upper() != 'NONE':
        cluster_gtf_file = xlib.get_cluster_reference_file(reference_dataset_id, gtf_file)
    else:
        cluster_gtf_file = 'NONE'

    # set the cluster splice site file
    cluster_splice_site_file = xlib.get_cluster_reference_file(reference_dataset_id, 'hisat2_splice_site_file.txt')

    # set the cluster exon file
    cluster_exon_file = xlib.get_cluster_reference_file(reference_dataset_id, 'hisat2_exon_file.txt')

    # set the directory and basename of the index HISAT2
    (reference_file_name, _) = os.path.splitext(reference_file)
    hisat2_index_dir = f'{cluster_reference_dataset_dir}/{reference_file_name}-hisat2_indexes'
    hisat2_index_basename = 'hisat2_indexes'

    # write the GSMAP process script
    try:
        if not os.path.exists(os.path.dirname(get_hisat2_process_script())):
            os.makedirs(os.path.dirname(get_hisat2_process_script()))
        with open(get_hisat2_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function print_hisat2_version\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_hisat2_anaconda_code()}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    hisat2 --version\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '}\n')
            if reference_dataset_id.upper() != 'NONE' and index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_splice_site_file\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    source activate {xlib.get_hisat2_anaconda_code()}\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building splice site file ..."\n')
                script_file_id.write(f'    echo -n > {cluster_splice_site_file}\n')
                script_file_id.write(f'    if [ "{cluster_gtf_file}" != "NONE" ]; then\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            hisat2_extract_splice_sites.py \\\n')
                script_file_id.write(f'                {cluster_gtf_file} \\\n')
                script_file_id.write(f'                >> {cluster_splice_site_file} \n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error hisat2_extract_splice_sites.py $RC; fi\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '    echo "The file is built."\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_exon_file\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    source activate {xlib.get_hisat2_anaconda_code()}\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building exon file ..."\n')
                script_file_id.write(f'    echo -n > {cluster_exon_file} \n')
                script_file_id.write(f'    if [ "{cluster_gtf_file}" != "NONE" ]; then\n')
                script_file_id.write( '        /usr/bin/time \\\n')
                script_file_id.write(f'            --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '            hisat2_extract_exons.py \\\n')
                script_file_id.write(f'                {cluster_gtf_file} \\\n')
                script_file_id.write(f'                >> {cluster_exon_file} \n')
                script_file_id.write( '        RC=$?\n')
                script_file_id.write( '        if [ $RC -ne 0 ]; then manage_error hisat2_extract_exons.py $RC; fi\n')
                script_file_id.write( '    fi\n')
                script_file_id.write( '    echo "The file is built."\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '}\n')
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_hisat2_indexes\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    source activate {xlib.get_hisat2_anaconda_code()}\n')
                script_file_id.write( '    cd $CURRENT_DIR\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building indexes ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        hisat2-build \\\n')
                script_file_id.write(f'            -p {threads} \\\n')
                script_file_id.write( '            -f \\\n')
                if large_index.upper() == 'YES':
                    script_file_id.write( '            --large-index \\\n')
                if cluster_gtf_file != 'NONE':
                    script_file_id.write(f'            --ss {cluster_splice_site_file} \\\n')
                    script_file_id.write(f'            --exon {cluster_exon_file} \\\n')
                script_file_id.write(f'            {cluster_reference_file} \\\n')
                script_file_id.write(f'            {hisat2_index_basename}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error hisat2-build $RC; fi\n')
                script_file_id.write(f'    mkdir --parents {hisat2_index_dir}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mkdir $RC; fi\n')
                script_file_id.write(f'    mv -f {hisat2_index_basename}.* {hisat2_index_dir}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mv $RC; fi\n')
                script_file_id.write( '    echo "Indexes are built."\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_hisat2_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_hisat2_anaconda_code()}\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            for i in range(len(read_file_1_list)):
                # set file names for the library
                if library_concatenation.upper() == 'YES':
                    library_name = 'concatenated_libraries'
                    alignment_file = 'alignment.sam'
                    un_gz = 'unpairednotaligned.fastq.gz'
                    al_gz = 'unpairedaligned.fastq.gz'
                    un_conc_gz = 'pairednotaligned.fastq.gz'
                    al_conc_gz = 'pairednotaligned.fastq.gz'
                    summary_file = 'summary.txt'
                else:
                    if read_file_1.endswith('.gz'):
                        (library_name, _) = os.path.splitext(os.path.basename(read_file_1_list[i][:-3]))
                    else:
                        (library_name, _) = os.path.splitext(os.path.basename(read_file_1_list[i]))
                    alignment_file = f'{library_name}-alignment.sam'
                    un_gz = f'{library_name}-unpairednotaligned.fastq.gz'
                    al_gz = f'{library_name}-unpairedaligned.fastq.gz'
                    un_conc_gz = f'{library_name}-pairednotaligned.fastq.gz'
                    al_conc_gz = f'{library_name}-pairednotaligned.fastq.gz'
                    summary_file = f'{library_name}-summary.txt'
                # write the instructions for the library
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    echo "Mapping reads of {library_name} ..."\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
                script_file_id.write( '        hisat2 \\\n')
                script_file_id.write(f'            --threads {threads} \\\n')
                if dta_cufflinks.upper() == 'YES':
                    script_file_id.write( '            --dta-cufflinks \\\n')
                script_file_id.write(f'            --mp {max_mp},{min_mp} \\\n')
                if no_softclip.upper() != 'YES':
                    script_file_id.write( '            --no-softclip \\\n')
                else:
                    script_file_id.write(f'            --sp {max_sp},{min_sp} \\\n')
                script_file_id.write(f'            --np {np} \\\n')
                script_file_id.write(f'            --rdg {open_rdg},{extend_rdg} \\\n')
                script_file_id.write(f'            --rfg {open_rfg},{extend_rfg} \\\n')
                script_file_id.write(f'            --{orientation.lower()} \\\n')
                if other_parameters.upper() != 'NONE':
                    parameter_list = [x.strip() for x in other_parameters.split(';')]
                    for j in range(len(parameter_list)):
                        if parameter_list[j].find('=') > 0:
                            pattern = r'^--(.+)=(.+)$'
                            mo = re.search(pattern, parameter_list[j])
                            parameter_name = mo.group(1).strip()
                            parameter_value = mo.group(2).strip()
                            script_file_id.write(f'            --{parameter_name} {parameter_value} \\\n')
                        else:
                            pattern = r'^--(.+)$'
                            mo = re.search(pattern, parameter_list[j])
                            parameter_name = mo.group(1).strip()
                            script_file_id.write(f'            --{parameter_name} \\\n')
                script_file_id.write(f'            --known-splicesite-infile {cluster_splice_site_file} \\\n')
                script_file_id.write(f'            -x {hisat2_index_dir}/{hisat2_index_basename} \\\n')
                if format.upper() == 'FASTQ':
                    script_file_id.write( '            -q \\\n')
                elif format.upper() == 'FASTA':
                    script_file_id.write( '            -f \\\n')
                if read_type.upper() == 'SE':
                    script_file_id.write(f'            -U {read_file_1_list[i]}\n')
                elif read_type.upper() == 'PE':
                    script_file_id.write(f'            -1 {read_file_1_list[i]} \\\n')
                    script_file_id.write(f'            -2 {read_file_2_list[i]} \\\n')
                script_file_id.write(f'            -S {alignment_file} \\\n')
                script_file_id.write(f'            --un-gz {un_gz} \\\n')
                script_file_id.write(f'            --al-gz {al_gz} \\\n')
                script_file_id.write(f'            --un-conc-gz {un_conc_gz} \\\n')
                script_file_id.write(f'            --al-conc-gz {al_conc_gz} \\\n')
                script_file_id.write(f'            --summary-file {summary_file} \\\n')
                script_file_id.write( '            --time\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error hisat2 $RC; fi\n')
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
            script_file_id.write( '    ls *.bam > bam-files.txt\n')
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
            process_name = f'{xlib.get_hisat2_name()} process'
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
            script_file_id.write( 'print_hisat2_version\n')
            if reference_dataset_id.upper() != 'NONE' and index_building.upper() == 'YES':
                script_file_id.write( 'build_splice_site_file\n')
                script_file_id.write( 'build_exon_file\n')
                script_file_id.write( 'build_hisat2_indexes\n')
            script_file_id.write( 'run_hisat2_process\n')
            script_file_id.write( 'convert_sam2bam\n')
            script_file_id.write( 'sort_and_index_bam_files\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_hisat2_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_hisat2_process_starter(current_run_dir):
    '''
    Build the starter of the current HISAT2 process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the HISAT2 process starter
    try:
        if not os.path.exists(os.path.dirname(get_hisat2_process_starter())):
            os.makedirs(os.path.dirname(get_hisat2_process_starter()))
        with open(get_hisat2_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_hisat2_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_hisat2_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_hisat2_config_file():
    '''
    Get the HISAT2 config file path.
    '''

    # assign the HISAT2 config file path
    hisat2_config_file = f'{xlib.get_config_dir()}/{xlib.get_hisat2_code()}-config.txt'

    # return the HISAT2 config file path
    return hisat2_config_file

#-------------------------------------------------------------------------------

def get_hisat2_process_script():
    '''
    Get the HISAT2 process script path in the local computer.
    '''

    # assign the HISAT2 script path
    hisat2_process_script = f'{xlib.get_temp_dir()}/{xlib.get_hisat2_code()}-process.sh'

    # return the HISAT2 script path
    return hisat2_process_script

#-------------------------------------------------------------------------------

def get_hisat2_process_starter():
    '''
    Get the HISAT2 process starter path in the local computer.
    '''

    # assign the HISAT2 process starter path
    hisat2_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_hisat2_code()}-process-starter.sh'

    # return the HISAT2 starter path
    return hisat2_process_starter

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

def get_dta_cufflinks_code_list():
    '''
    Get the code list of "dta_cufflinks".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------

def get_dta_cufflinks_code_list_text():
    '''
    Get the code list of "dta_cufflinks" as text.
    '''

    return str(get_no_softclip_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_no_softclip_code_list():
    '''
    Get the code list of "no_softclip".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------

def get_no_softclip_code_list_text():
    '''
    Get the code list of "no_softclip" as text.
    '''

    return str(get_no_softclip_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_large_index_code_list():
    '''
    Get the code list of "large_index".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------

def get_large_index_code_list_text():
    '''
    Get the code list of "large_index" as text.
    '''

    return str(get_large_index_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

def get_orientation_code_list():
    '''
    Get the code list of "orientation".
    '''

    return ['FR', 'RF', 'FF']

#-------------------------------------------------------------------------------

def get_orientation_code_list_text():
    '''
    Get the code list of "orientation" as text.
    '''

    return 'FR (fwd-rev, or typical Illumina) or RF (rev-fwd, for circularized inserts) or FF (fwd-fwd, same strand)'

#-------------------------------------------------------------------------------

def get_quality_score_code_list():
    '''
    Get the code list of "quality_score".
    '''

    return ['33', '64']

#-------------------------------------------------------------------------------

def get_quality_score_code_list_text():
    '''
    Get the code list of "quality_score" as text.
    '''

    return '33 (Phred+33) or 64 (Phred+64)'

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

def get_library_concatenation_code_list():
    '''
    Get the code list of "library_concatenation".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------

def get_library_concatenation_code_list_text():
    '''
    Get the code list of "library_concatenation" as text.
    '''

    return 'YES (map concatanated libraries) or NO (map each library separately)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the HISAT2 process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
