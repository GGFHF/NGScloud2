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
This file contains functions related to the GAMP-GSNAP process used in both console
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

def create_gmap_config_file(experiment_id='exp001', reference_dataset_id='NONE', reference_file='NONE', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create GMAP config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # set the app
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

    # create the GMAP config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_gmap_config_file())):
            os.makedirs(os.path.dirname(get_gmap_config_file()))
        with open(get_gmap_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write( '# The experiment_id, reference_dataset_id, reference_file and assembly_dataset_id are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of GMAP and their meaning in "http://research-pub.gene.com/gmap/".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "GMAP parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of GMAP and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --no-chimeras; --canonical-mode=2\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the GMAP parameters\n')
            file_id.write( '[GMAP parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'index_building = YES', f'# index building: {get_index_building_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'gmap_version = gmapl', f'# GMAP version: {get_gmap_version_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'kmer = NONE', '# kmer size to use in genome database or NONE (the program will find the highest available kmer size in the genome database)'))
            file_id.write( '{0:<50} {1}\n'.format( 'sampling = NONE', '# Sampling to use in genome database or NONE (the program will find the smallest available sampling value in the genome database within selected k-mer size)'))
            file_id.write( '{0:<50} {1}\n'.format( 'input-buffer-size = 1000', '# size of input buffer'))
            file_id.write( '{0:<50} {1}\n'.format( 'output-buffer-size = 1000', '# size of buffer size in queries for output thread'))
            file_id.write( '{0:<50} {1}\n'.format( 'prunelevel = 0', f'# pruning level: {get_prunelevel_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'format = COMPRESS', f'# format for output: {get_gmap_output_format_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gmap_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_gmap_process(cluster_name, log, function=None):
    '''
    Run a GMAP process.
    '''

    # initialize the control variable
    OK = True

    # get the GMAP option dictionary
    gmap_option_dict = xlib.get_option_dict(get_gmap_config_file())

    # get the experiment identification
    experiment_id = gmap_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the GMAP config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_gmap_name()} config file ...\n')
    (OK, error_list) = check_gmap_config_file(strict=True)
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

    # check the GMAP-GSNAP is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_gmap_gsnap_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_gmap_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_gmap_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_gmap_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GMAP process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_gmap_process_script()} ...\n')
        (OK, error_list) = build_gmap_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the GMAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_gmap_process_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_gmap_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_gmap_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GMAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_gmap_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_gmap_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GMAP process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_gmap_process_starter()} ...\n')
        (OK, error_list) = build_gmap_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the GMAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_gmap_process_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_gmap_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_gmap_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GMAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_gmap_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_gmap_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the GMAP process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_gmap_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_gmap_process_starter()), log)

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

def check_gmap_config_file(strict):
    '''
    Check the GMAP config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        gmap_option_dict = xlib.get_option_dict(get_gmap_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in gmap_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = gmap_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = gmap_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = gmap_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = gmap_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = gmap_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" hast to start with {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = gmap_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
                    OK = False

        # check section "GMAP parameters"
        if 'GMAP parameters' not in sections_list:
            error_list.append('*** ERROR: the section "GMAP parameters" is not found.')
            OK = False
        else:

            # check section "GMAP parameters" - key "index_building"
            index_building = gmap_option_dict.get('GMAP parameters', {}).get('index_building', not_found)
            if index_building == not_found:
                error_list.append('*** ERROR: the key "index_building" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(index_building, get_index_building_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "index_building" has to be {get_index_building_code_list_text()}.')
                OK = False

            # check section "GMAP parameters" - key "version"
            gmap_version = gmap_option_dict.get('GMAP parameters', {}).get('gmap_version', not_found)
            if gmap_version == not_found:
                error_list.append('*** ERROR: the key "gmap_version" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(gmap_version, get_gmap_version_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "gmap_version" has to be {get_gmap_version_code_list_text()}.')
                OK = False

            # check section "GMAP parameters" - key "threads"
            threads = gmap_option_dict.get('GMAP parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GMAP parameters" - key "kmer"
            kmer = gmap_option_dict.get('GMAP parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section "GMAP parameters".')
                OK = False
            elif kmer.upper() != 'NONE' and not xlib.check_int(kmer, minimum=1, maximum=30):
                error_list.append('*** ERROR: the key "kmer" has to be an integer number between 1 and 16 or NONE.')
                OK = False

            # check section "GMAP parameters" - key "sampling"
            sampling = gmap_option_dict.get('GMAP parameters', {}).get('sampling', not_found)
            if sampling == not_found:
                error_list.append('*** ERROR: the key "sampling" is not found in the section "GMAP parameters".')
                OK = False
            elif sampling.upper() != 'NONE' and not xlib.check_int(sampling, minimum=1):
                error_list.append('*** ERROR: the key "sampling" has to be an integer number greater than or equal to 1 or NONE.')
                OK = False

            # check section "GMAP parameters" - key "input-buffer-size"
            input_buffer_size = gmap_option_dict.get('GMAP parameters', {}).get('input-buffer-size', not_found)
            if input_buffer_size == not_found:
                error_list.append('*** ERROR: the key "input-buffer-size" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_int(input_buffer_size, minimum=1):
                error_list.append('*** ERROR: the key "input-buffer-size" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GMAP parameters" - key "output-buffer-size"
            output_buffer_size = gmap_option_dict.get('GMAP parameters', {}).get('output-buffer-size', not_found)
            if output_buffer_size == not_found:
                error_list.append('*** ERROR: the key "output-buffer-size" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_int(output_buffer_size, minimum=1):
                error_list.append('*** ERROR: the key "output-buffer-size" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GMAP parameters" - key "prunelevel"
            prunelevel = gmap_option_dict.get('GMAP parameters', {}).get('prunelevel', not_found)
            if prunelevel == not_found:
                error_list.append('*** ERROR: the key "prunelevel" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(prunelevel, get_prunelevel_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "prunelevel" has to be {get_prunelevel_code_list_text()}.')
                OK = False

            # check section "GMAP parameters" - key "format"
            format = gmap_option_dict.get('GMAP parameters', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(format, get_gmap_output_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_gmap_output_format_code_list_text()}.')
                OK = False

            # check section "GMAP parameters" - key "other_parameters"
            not_allowed_parameters_list = ['nthreads', 'kmer', 'sampling', 'input-buffer-size', 'output-buffer-size', 'prunelevel', 'compress', 'summary', 'align', 'format', 'ordered']
            other_parameters = gmap_option_dict.get('GMAP parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "GMAP parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_gmap_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gmap_process_script(cluster_name, current_run_dir):
    '''
    Build the current GMAP process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the GMAP option dictionary
    gmap_option_dict = xlib.get_option_dict(get_gmap_config_file())

    # get the options
    experiment_id = gmap_option_dict['identification']['experiment_id']
    reference_dataset_id = gmap_option_dict['identification']['reference_dataset_id']
    reference_file = gmap_option_dict['identification']['reference_file']
    assembly_software = gmap_option_dict['identification']['assembly_software']
    assembly_dataset_id = gmap_option_dict['identification']['assembly_dataset_id']
    assembly_type = gmap_option_dict['identification']['assembly_type']
    index_building = gmap_option_dict['GMAP parameters']['index_building']
    gmap_version = gmap_option_dict['GMAP parameters']['gmap_version']
    threads = gmap_option_dict['GMAP parameters']['threads']
    kmer = gmap_option_dict['GMAP parameters']['kmer']
    sampling = gmap_option_dict['GMAP parameters']['sampling']
    input_buffer_size = gmap_option_dict['GMAP parameters']['input-buffer-size']
    output_buffer_size = gmap_option_dict['GMAP parameters']['output-buffer-size']
    prunelevel = gmap_option_dict['GMAP parameters']['prunelevel']
    format = gmap_option_dict['GMAP parameters']['format']
    other_parameters = gmap_option_dict['GMAP parameters']['other_parameters']

    # set the cluster reference dataset directory
    cluster_reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)

    # set the cluster reference file
    cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

    # set the GMAP database name
    reference_file_name, reference_file_extension = os.path.splitext(reference_file)
    gmap_database = f'{reference_file_name}-gmap_database'

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

    # set the output file path
    output_file = 'gmap_output_fasta'

    # write the GMAP process script
    try:
        if not os.path.exists(os.path.dirname(get_gmap_process_script())):
            os.makedirs(os.path.dirname(get_gmap_process_script()))
        with open(get_gmap_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'GMAP_GSNAP_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_gmap_gsnap_anaconda_code()}/bin\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( 'PATH=$GMAP_GSNAP_PATH:$MINICONDA3_BIN_PATH:$PATH\n')
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
            if index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_gmap_database\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    source activate {xlib.get_gmap_gsnap_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '        gmap_build \\\n')
                script_file_id.write(f'            --dir={cluster_reference_dataset_dir} \\\n')
                script_file_id.write(f'            --db={gmap_database} \\\n')
                if kmer.upper() != 'NONE':
                    script_file_id.write(f'            --kmer={kmer} \\\n')
                script_file_id.write(f'            {cluster_reference_file}\n'.format(''.format()))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error gmap_build $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_gmap_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    gmap --version\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write(f'    source activate {xlib.get_gmap_gsnap_anaconda_code()}\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write(f'        {gmap_version.lower()} \\\n')
            script_file_id.write(f'            --nthreads={threads} \\\n')
            script_file_id.write(f'            --dir={cluster_reference_dataset_dir} \\\n')
            script_file_id.write(f'            --db={gmap_database} \\\n')
            if kmer.upper() != 'NONE':
                script_file_id.write(f'            --kmer={kmer} \\\n')
            if sampling.upper() != 'NONE':
                script_file_id.write(f'            --sampling={sampling} \\\n')
            script_file_id.write(f'            --input-buffer-size={input_buffer_size} \\\n')
            script_file_id.write(f'            --output-buffer-size={output_buffer_size} \\\n')
            script_file_id.write(f'            --prunelevel={prunelevel} \\\n')
            if format.upper() == 'COMPRESS':
                script_file_id.write( '            --compress \\\n')
            elif format.upper() == 'SUMMARY':
                script_file_id.write( '            --summary \\\n')
            elif format.upper() == 'ALIGN':
                script_file_id.write( '            --align \\\n')
            else:
                script_file_id.write(f'            --format={format.lower()} \\\n')
            script_file_id.write( '            --ordered \\\n')
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for i in range(len(parameter_list)):
                    if parameter_list[i].find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write(f'            --{parameter_name}={parameter_value} \\\n')
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        script_file_id.write(f'            --{parameter_name} \\\n')
            script_file_id.write(f'            {transcriptome_file} \\\n')
            script_file_id.write(f'            > {output_file}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error gmap $RC; fi\n')
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
            process_name = f'{xlib.get_gmap_name()} process'
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
            if index_building.upper() == 'YES':
                script_file_id.write(f'build_gmap_database\n')
            script_file_id.write( 'run_gmap_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gmap_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gmap_process_starter(current_run_dir):
    '''
    Build the starter of the current GMAP process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the GMAP process starter
    try:
        if not os.path.exists(os.path.dirname(get_gmap_process_starter())):
            os.makedirs(os.path.dirname(get_gmap_process_starter()))
        with open(get_gmap_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gmap_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gmap_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gmap_config_file():
    '''
    Get the GMAP config file path.
    '''

    # assign the GMAP config file path
    gmap_config_file = f'{xlib.get_config_dir()}/{xlib.get_gmap_code()}-config.txt'

    # return the GMAP config file path
    return gmap_config_file

#-------------------------------------------------------------------------------

def get_gmap_process_script():
    '''
    Get the GMAP process script path in the local computer.
    '''

    # assign the GMAP script path
    gmap_process_script = f'{xlib.get_temp_dir()}/{xlib.get_gmap_code()}-process.sh'

    # return the GMAP script path
    return gmap_process_script

#-------------------------------------------------------------------------------

def get_gmap_process_starter():
    '''
    Get the GMAP process starter path in the local computer.
    '''

    # assign the GMAP process starter path
    gmap_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_gmap_code()}-process-starter.sh'

    # return the GMAP starter path
    return gmap_process_starter

#-------------------------------------------------------------------------------

def create_gsnap_config_file(experiment_id='exp001', reference_dataset_id='NONE', reference_file='NONE', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create GSNAP config file with the default options. It is necessary
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
    elif assembly_dataset_id.startswith(xlib.get_soapdenovo2_code()):
        assembly_software = xlib.get_soapdenovo2_code()
    elif assembly_dataset_id.startswith(xlib.get_starcode_code()):
        assembly_software = xlib.get_starcode_code()
    elif assembly_dataset_id.upper() == 'NONE':
        assembly_software = 'NONE'

    # create the GSNAP config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_gsnap_config_file())):
            os.makedirs(os.path.dirname(get_gsnap_config_file()))
        with open(get_gsnap_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '# The experiment_id, reference_dataset_id, reference_file, assembly_dataset_id and read_dataset_id are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of GSNAP and their meaning in "http://research-pub.gene.com/gmap/".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "GSNAP parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of GSNAP and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --allow-pe-name-mismatch; --max-gmap-improvement=5\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_extended_assembly_software_code_list_text()}; or NONE if a reference is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification or NONE if a reference is used'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}'.format(), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the GSNAP parameters\n')
            file_id.write( '[GSNAP parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'index_building = YES', f'# index building when a reference is used: {get_index_building_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'gsnap_version = gsnapl', f'# GSNAP version: {get_gsnap_version_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'kmer = NONE', '# kmer size to use in genome database or NONE (the program will find the highest available kmer size in the genome database)'))
            file_id.write( '{0:<50} {1}\n'.format( 'sampling = NONE', '# Sampling to use in genome database or NONE (the program will find the smallest available sampling value in the genome database within selected k-mer size)'))
            file_id.write( '{0:<50} {1}\n'.format( 'input-buffer-size = 1000', '# size of input buffer'))
            file_id.write( '{0:<50} {1}\n'.format( 'output-buffer-size = 1000', '# size of buffer size in queries for output thread'))
            file_id.write( '{0:<50} {1}\n'.format( 'max-mismatches = NONE', '# maximum number of mismatches allowed or NONE (the program will calculate the ultrafast level)'))
            file_id.write( '{0:<50} {1}\n'.format( 'indel-endlength = 4', '# minimum length at end required for indel alignments'))
            file_id.write( '{0:<50} {1}\n'.format( 'orientation = FR', f'# orientation of paired-end reads: {get_orientation_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'quality-zero-score = 33', '# FASTQ quality scores are zero at this ASCII value'))
            file_id.write( '{0:<50} {1}\n'.format( 'quality-print-shift = 0', '# shift FASTQ quality scores by this amount in output'))
            file_id.write( '{0:<50} {1}\n'.format( 'format = SAM', f'# format for output: {get_gsnap_output_format_code_list_text()}'))
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
        error_list.append(f'*** ERROR: The file {get_gsnap_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_gsnap_process(cluster_name, log, function=None):
    '''
    Run a GSNAP process.
    '''

    # initialize the control variable
    OK = True

    # get the GSNAP option dictionary
    gsnap_option_dict = xlib.get_option_dict(get_gsnap_config_file())

    # get the experiment identification
    experiment_id = gsnap_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the GSNAP config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_gsnap_name()} config file ...\n')
    (OK, error_list) = check_gsnap_config_file(strict=True)
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

    # check the GMAP-GSNAP is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_gmap_gsnap_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_gmap_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_gsnap_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_gsnap_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GSNAP process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_gsnap_process_script()} ...\n')
        (OK, error_list) = build_gsnap_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the GSNAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_gsnap_process_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_gsnap_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_gsnap_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GSNAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_gsnap_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_gsnap_process_script())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GSNAP process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_gsnap_process_starter()} ...\n')
        (OK, error_list) = build_gsnap_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the GSNAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_gsnap_process_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_gsnap_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_gsnap_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GSNAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_gsnap_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_gsnap_process_starter())}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the GSNAP process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_gsnap_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_gsnap_process_starter()), log)

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

def check_gsnap_config_file(strict):
    '''
    Check the GSNAP config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        gsnap_option_dict = xlib.get_option_dict(get_gsnap_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in gsnap_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = gsnap_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = gsnap_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            is_ok_reference_dataset_id = False
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False
            else:
                is_ok_reference_dataset_id = True

            # check section "identification" - key "reference_file"
            reference_file = gsnap_option_dict.get('identification', {}).get('reference_file', not_found)
            is_ok_reference_file = False
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False
            else:
                is_ok_reference_file = True

            # check that "reference_file" has to be NONE if "reference_dataset_id" is NONE
            if is_ok_reference_dataset_id and is_ok_reference_file and reference_dataset_id.upper() == 'NONE' and reference_file.upper() != 'NONE':
                error_list.append('*** ERROR: "reference_file" has to be NONE if "reference_dataset_id" is NONE.')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = gsnap_option_dict.get('identification', {}).get('assembly_software', not_found)
            is_ok_assembly_software = False
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif assembly_software.upper() != 'NONE' and not xlib.check_code(assembly_software, get_extended_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_extended_assembly_software_code_list_text()}; or NONE if a reference is used.')
                OK = False
            else:
                is_ok_assembly_software = True

            # check that "assembly_software" has to be NONE if "reference_dataset_id" is not NONE, and vice versa
            if is_ok_reference_dataset_id and is_ok_assembly_software and (reference_dataset_id.upper() == 'NONE' and assembly_software.upper() == 'NONE' or reference_dataset_id.upper() != 'NONE' and assembly_software.upper() != 'NONE'):
                error_list.append('*** ERROR: "assembly_software" has to be NONE if "reference_dataset_id" is not NONE, and vice versa.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = gsnap_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            is_ok_assembly_dataset_id = False
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.upper() != 'NONE' and not xlib.check_startswith(assembly_dataset_id, get_extended_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" does not have to start with {get_extended_assembly_software_code_list_text()}.')
                OK = False
            else:
                is_ok_assembly_dataset_id = True

            # check that "assembly_dataset_id" has to be NONE if "assembly_software" is NONE
            if is_ok_assembly_software and is_ok_assembly_dataset_id and assembly_software.upper() == 'NONE' and assembly_dataset_id.upper() != 'NONE':
                error_list.append('*** ERROR: "assembly_dataset_id" has to be NONE if "assembly_software" is NONE.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = gsnap_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif (assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or assembly_dataset_id.startswith(xlib.get_soapdenovo2_code())) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not (assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) or assembly_dataset_id.startswith(xlib.get_soapdenovo2_code())) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} and {xlib.get_soapdenovo2_name()}; or NONE in any other case.')
                    OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = gsnap_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "GSNAP parameters"
        if 'GSNAP parameters' not in sections_list:
            error_list.append('*** ERROR: the section "GSNAP parameters" is not found.')
            OK = False
        else:

            # check section "GSNAP parameters" - key "index_building"
            index_building = gsnap_option_dict.get('GSNAP parameters', {}).get('index_building', not_found)
            if index_building == not_found:
                error_list.append('*** ERROR: the key "index_building" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_code(index_building, get_index_building_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "index_building" has to be {get_index_building_code_list_text()}.')
                OK = False

            # check section "GSNAP parameters" - key "version"
            gsnap_version = gsnap_option_dict.get('GSNAP parameters', {}).get('gsnap_version', not_found)
            if gsnap_version == not_found:
                error_list.append('*** ERROR: the key "gsnap_version" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_code(gsnap_version, get_gsnap_version_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "gsnap_version" has to be {get_gsnap_version_code_list_text()}.')
                OK = False

            # check section "GSNAP parameters" - key "threads"
            threads = gsnap_option_dict.get('GSNAP parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GSNAP parameters" - key "kmer"
            kmer = gsnap_option_dict.get('GSNAP parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section "GSNAP parameters".')
                OK = False
            elif kmer.upper() != 'NONE' and not xlib.check_int(kmer, minimum=1, maximum=30):
                error_list.append('*** ERROR: the key "kmer" has to be an integer number between 1 and 16 or NONE.')
                OK = False

            # check section "GSNAP parameters" - key "sampling"
            sampling = gsnap_option_dict.get('GSNAP parameters', {}).get('sampling', not_found)
            if sampling == not_found:
                error_list.append('*** ERROR: the key "sampling" is not found in the section "GSNAP parameters".')
                OK = False
            elif sampling.upper() != 'NONE' and not xlib.check_int(sampling, minimum=1):
                error_list.append('*** ERROR: the key "sampling" has to be an integer number greater than or equal to 1 or NONE.')
                OK = False

            # check section "GSNAP parameters" - key "input-buffer-size"
            input_buffer_size = gsnap_option_dict.get('GSNAP parameters', {}).get('input-buffer-size', not_found)
            if input_buffer_size == not_found:
                error_list.append('*** ERROR: the key "input-buffer-size" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(input_buffer_size, minimum=1):
                error_list.append('*** ERROR: the key "input-buffer-size" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GSNAP parameters" - key "output-buffer-size"
            output_buffer_size = gsnap_option_dict.get('GSNAP parameters', {}).get('output-buffer-size', not_found)
            if output_buffer_size == not_found:
                error_list.append('*** ERROR: the key "output-buffer-size" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(output_buffer_size, minimum=1):
                error_list.append('*** ERROR: the key "output-buffer-size" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GSNAP parameters" - key "max-mismatches"
            max_mismatches = gsnap_option_dict.get('GSNAP parameters', {}).get('max-mismatches', not_found)
            if max_mismatches == not_found:
                error_list.append('*** ERROR: the key "max-mismatches" is not found in the section "GSNAP parameters".')
                OK = False
            elif  max_mismatches.upper() != 'NONE' and not xlib.check_float(max_mismatches, minimum=0., mne=0.):
                error_list.append('*** ERROR: the key "max-mismatches" has to be a float number greater than 0.0 or NONE.')
                OK = False

            # check section "GSNAP parameters" - key "indel-endlength"
            indel_endlength = gsnap_option_dict.get('GSNAP parameters', {}).get('indel-endlength', not_found)
            if indel_endlength == not_found:
                error_list.append('*** ERROR: the key "indel-endlength" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(indel_endlength, minimum=1):
                error_list.append('*** ERROR: the key "indel-endlength" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "GSNAP parameters" - key "orientation"
            orientation = gsnap_option_dict.get('GSNAP parameters', {}).get('orientation', not_found)
            if orientation == not_found:
                error_list.append('*** ERROR: the key "orientation" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_code(orientation, get_orientation_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "orientation" has to be {get_orientation_code_list_text()}.')
                OK = False

            # check section "GSNAP parameters" - key "quality-zero-score"
            quality_zero_score = gsnap_option_dict.get('GSNAP parameters', {}).get('quality-zero-score', not_found)
            if quality_zero_score == not_found:
                error_list.append('*** ERROR: the key "quality-zero-score" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(quality_zero_score):
                error_list.append('*** ERROR: the key "quality-zero-score" has to be an integer number.')
                OK = False

            # check section "GSNAP parameters" - key "quality-print-shift"
            quality_print_shift = gsnap_option_dict.get('GSNAP parameters', {}).get('quality-print-shift', not_found)
            if quality_print_shift == not_found:
                error_list.append('*** ERROR: the key "quality-print-shift" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_int(quality_print_shift):
                error_list.append('*** ERROR: the key "quality-print-shift" has to be an integer number.')
                OK = False

            # check section "GSNAP parameters" - key "format"
            format = gsnap_option_dict.get('GSNAP parameters', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_code(format, get_gsnap_output_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_gsnap_output_format_code_list_text()}.')
                OK = False

            # check section "GSNAP parameters" - key "other_parameters"
            not_allowed_parameters_list = ['nthreads', 'kmer', 'sampling', 'input-buffer-size', 'output-buffer-size', 'max-mismatches', 'indel-endlength', 'orientation', 'force-single-end', 'quality-protocol', 'quality-zero-score', 'quality-print-shift', 'format', 'ordered', 'split-output', 'failed-input']
            other_parameters = gsnap_option_dict.get('GSNAP parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "GSNAP parameters".')
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
            format = gsnap_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = gsnap_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'GSNAP parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = gsnap_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = gsnap_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_gsnap_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gsnap_process_script(cluster_name, current_run_dir):
    '''
    Build the current GSNAP process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the GSNAP option dictionary
    gsnap_option_dict = xlib.get_option_dict(get_gsnap_config_file())

    # get the options
    experiment_id = gsnap_option_dict['identification']['experiment_id']
    reference_dataset_id = gsnap_option_dict['identification']['reference_dataset_id']
    reference_file = gsnap_option_dict['identification']['reference_file']
    assembly_software = gsnap_option_dict['identification']['assembly_software']
    assembly_dataset_id = gsnap_option_dict['identification']['assembly_dataset_id']
    assembly_type = gsnap_option_dict['identification']['assembly_type']
    read_dataset_id = gsnap_option_dict['identification']['read_dataset_id']
    index_building = gsnap_option_dict['GSNAP parameters']['index_building']
    gsnap_version = gsnap_option_dict['GSNAP parameters']['gsnap_version']
    threads = gsnap_option_dict['GSNAP parameters']['threads']
    kmer = gsnap_option_dict['GSNAP parameters']['kmer']
    sampling = gsnap_option_dict['GSNAP parameters']['sampling']
    input_buffer_size = gsnap_option_dict['GSNAP parameters']['input-buffer-size']
    output_buffer_size = gsnap_option_dict['GSNAP parameters']['output-buffer-size']
    max_mismatches = gsnap_option_dict['GSNAP parameters']['max-mismatches']
    indel_endlength = gsnap_option_dict['GSNAP parameters']['indel-endlength']
    orientation = gsnap_option_dict['GSNAP parameters']['orientation']
    quality_zero_score = gsnap_option_dict['GSNAP parameters']['quality-zero-score']
    quality_print_shift = gsnap_option_dict['GSNAP parameters']['quality-print-shift']
    format = gsnap_option_dict['GSNAP parameters']['format']
    other_parameters = gsnap_option_dict['GSNAP parameters']['other_parameters']
    read_type = gsnap_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in gsnap_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build read file lists
    read_file_1_list = []
    read_file_2_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = gsnap_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            read_file_1_list.append(read_file_1)
            if read_type.upper() == 'PE':
                read_file_2 = gsnap_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                read_file_2_list.append(read_file_2)

    # set the cluster reference dataset directory
    if reference_dataset_id.upper() != 'NONE':
        cluster_reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)
    else:
        cluster_reference_dataset_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)

    # set the cluster reference file
    if reference_dataset_id.upper() != 'NONE':
        cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    else:
        if assembly_software in [xlib.get_soapdenovotrans_code(), xlib.get_soapdenovo2_code()]:
            if assembly_type == 'CONTIGS':
                cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.contig'
            elif  assembly_type == 'SCAFFOLDS':
                cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/{experiment_id}-{assembly_dataset_id}.scafSeq'
        elif assembly_software == xlib.get_transabyss_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/transabyss-final.fa'
        elif assembly_software == xlib.get_trinity_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity.fasta'
        elif assembly_software == xlib.get_ggtrinity_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/Trinity-GG.fasta'
        elif assembly_software == xlib.get_cd_hit_est_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/clustered-transcriptome.fasta'
        elif assembly_software == xlib.get_transcript_filter_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/filtered-transcriptome.fasta'
        elif assembly_software == xlib.get_starcode_code():
            cluster_reference_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id)}/starcode.fasta'

    # set the GMAP database name
    if reference_file.upper() != 'NONE':
        reference_file_name, reference_file_extension = os.path.splitext(reference_file)
        gmap_database = f'{reference_file_name}-gmap_database'
    else:
        gmap_database = 'pseudogenome-gmap_database'

    # write the GSMAP process script
    try:
        if not os.path.exists(os.path.dirname(get_gsnap_process_script())):
            os.makedirs(os.path.dirname(get_gsnap_process_script()))
        with open(get_gsnap_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write(f'GMAP_GSNAP_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/envs/{xlib.get_gmap_gsnap_anaconda_code()}/bin\n')
            script_file_id.write(f'MINICONDA3_BIN_PATH={xlib.get_cluster_app_dir()}/{xlib.get_miniconda3_name()}/bin\n')
            script_file_id.write( 'PATH=$GMAP_GSNAP_PATH:$MINICONDA3_BIN_PATH:$PATH\n')
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
            if index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( 'function build_gmap_database\n')
                script_file_id.write( '{\n')
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '    echo "Building database ..."\n')
                script_file_id.write(f'    source activate {xlib.get_gmap_gsnap_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write( '        gmap_build \\\n')
                script_file_id.write(f'            --dir={cluster_reference_dataset_dir} \\\n')
                script_file_id.write(f'            --genomedb={gmap_database} \\\n')
                if kmer.upper() != 'NONE':
                    script_file_id.write(f'            --kmer={kmer} \\\n')
                script_file_id.write(f'            {cluster_reference_file}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error gmap_build $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
                script_file_id.write( '    echo "Database is built."\n')
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_gsnap_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Mapping reads ..."\n')
            for i in range(len(read_file_1_list)):
                # set gunzip, bunzip2, split_output and failed_input values
                gunzip = False
                bunzip2 = False
                if read_file_1_list[i].endswith('.gz'):
                    basename = read_file_1_list[i][:-3]
                    gunzip = True
                elif read_file_1_list[i].endswith('.bz2'):
                    bunzip2 = True
                    basename = read_file_1_list[i][:-4]
                else:
                    basename = read_file_1_list[i]
                position = basename[::-1].find('.')
                if position > -1:
                    basename = f'{current_run_dir}/{os.path.basename(basename[:len(basename)-position-1])}'
                else:
                    basename = f'{current_run_dir}/{os.path.basename(basename)}'
                split_output = f'{basename}-split'
                failed_input = f'{basename}-failed'
                # write the gsnap run instructions
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write(f'    source activate {xlib.get_gmap_gsnap_anaconda_code()}\n')
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
                script_file_id.write(f'        {gsnap_version.lower()} \\\n')
                script_file_id.write(f'            --nthreads={threads} \\\n')
                script_file_id.write(f'            --dir={cluster_reference_dataset_dir} \\\n')
                script_file_id.write(f'            --db={gmap_database} \\\n')
                if kmer.upper() != 'NONE':
                    script_file_id.write(f'            --kmer={kmer} \\\n')
                if sampling.upper() != 'NONE':
                    script_file_id.write(f'            --sampling={sampling} \\\n')
                script_file_id.write(f'            --input-buffer-size={input_buffer_size} \\\n')
                script_file_id.write(f'            --output-buffer-size={output_buffer_size} \\\n')
                if max_mismatches.upper() != 'NONE':
                    script_file_id.write(f'            --max-mismatches={max_mismatches} \\\n')
                script_file_id.write(f'            --indel-endlength={indel_endlength} \\\n')
                if read_type.upper() == 'PE':
                    script_file_id.write(f'            --orientation={orientation.upper()} \\\n')
                script_file_id.write(f'            --quality-zero-score={quality_zero_score} \\\n')
                script_file_id.write(f'            --quality-print-shift={quality_print_shift} \\\n')
                script_file_id.write(f'            --format={format.lower()} \\\n')
                script_file_id.write( '            --ordered \\\n')
                script_file_id.write(f'            --split-output={split_output} \\\n')
                script_file_id.write(f'            --failed-input={failed_input} \\\n')
                if other_parameters.upper() != 'NONE':
                    parameter_list = [x.strip() for x in other_parameters.split(';')]
                    for i in range(len(parameter_list)):
                        if parameter_list[i].find('=') > 0:
                            pattern = r'^--(.+)=(.+)$'
                            mo = re.search(pattern, parameter_list[i])
                            parameter_name = mo.group(1).strip()
                            parameter_value = mo.group(2).strip()
                            script_file_id.write(f'            --{parameter_name}={parameter_value} \\\n')
                        else:
                            pattern = r'^--(.+)$'
                            mo = re.search(pattern, parameter_list[i])
                            parameter_name = mo.group(1).strip()
                            script_file_id.write(f'            --{parameter_name} \\\n')
                if gunzip:
                    script_file_id.write( '            --gunzip \\\n')
                if bunzip2:
                    script_file_id.write( '            --bunzip2 \\\n')
                if read_type.upper() == 'SE':
                    script_file_id.write(f'            {read_file_1_list[i]}\n')
                if read_type.upper() == 'PE':
                    script_file_id.write(f'            {read_file_1_list[i]} \\\n')
                    script_file_id.write(f'            {read_file_2_list[i]}\n')
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error gsnap $RC; fi\n')
                script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '    echo "Reads are mapped."\n')
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
            process_name = f'{xlib.get_gsnap_name()} process'
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
            if index_building.upper() == 'YES':
                script_file_id.write( 'build_gmap_database\n')
            script_file_id.write( 'run_gsnap_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gsnap_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gsnap_process_starter(current_run_dir):
    '''
    Build the starter of the current GSNAP process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the GSNAP process starter
    try:
        if not os.path.exists(os.path.dirname(get_gsnap_process_starter())):
            os.makedirs(os.path.dirname(get_gsnap_process_starter()))
        with open(get_gsnap_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_gsnap_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_gsnap_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gsnap_config_file():
    '''
    Get the GSNAP config file path.
    '''

    # assign the GSNAP config file path
    gsnap_config_file = f'{xlib.get_config_dir()}/{xlib.get_gsnap_code()}-config.txt'

    # return the GSNAP config file path
    return gsnap_config_file

#-------------------------------------------------------------------------------

def get_gsnap_process_script():
    '''
    Get the GSNAP process script path in the local computer.
    '''

    # assign the GSNAP script path
    gsnap_process_script = f'{xlib.get_temp_dir()}/{xlib.get_gsnap_code()}-process.sh'

    # return the GSNAP script path
    return gsnap_process_script

#-------------------------------------------------------------------------------

def get_gsnap_process_starter():
    '''
    Get the GSNAP process starter path in the local computer.
    '''

    # assign the GSNAP process starter path
    gsnap_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_gsnap_code()}-process-starter.sh'

    # return the GSNAP starter path
    return gsnap_process_starter

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

def get_extended_assembly_software_code_list():
    '''
    Get the code list of "assembly_software".
    '''

    return [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(),  xlib.get_transcript_filter_code(), xlib.get_soapdenovo2_code(), xlib.get_starcode_name()]

#-------------------------------------------------------------------------------

def get_extended_assembly_software_code_list_text():
    '''
    Get the code list of "assembly_software" as text.
    '''

    return f'{xlib.get_soapdenovotrans_code()} ({xlib.get_soapdenovotrans_name()}) or {xlib.get_transabyss_code()} ({xlib.get_transabyss_name()}) or {xlib.get_trinity_code()} ({xlib.get_trinity_name()}) or {xlib.get_ggtrinity_code()} ({xlib.get_ggtrinity_name()}) or {xlib.get_cd_hit_est_code()} ({xlib.get_cd_hit_est_name()}) or {xlib.get_transcript_filter_code()} ({xlib.get_transcript_filter_name()}) or {xlib.get_soapdenovo2_code()} ({xlib.get_soapdenovo2_name()}) or {xlib.get_starcode_code()} ({xlib.get_starcode_name()})'

#-------------------------------------------------------------------------------
    
def get_gmap_version_code_list():
    '''
    Get the code list of "gmap_version".
    '''

    return ['gmap', 'gmapl']

#-------------------------------------------------------------------------------
    
def get_gmap_version_code_list_text():
    '''
    Get the code list of "gmap_version" as text.
    '''

    return 'gmap (small genones) or gmapl (large genomes of more than 2^32 bp in total length)'

#-------------------------------------------------------------------------------
    
def get_gsnap_version_code_list():
    '''
    Get the code list of "gsnap_version".
    '''

    return ['gsnap', 'gsnapl']

#-------------------------------------------------------------------------------
    
def get_gsnap_version_code_list_text():
    '''
    Get the code list of "gsnap_version" as text.
    '''

    return 'gsnap (small genones) or gsnapl (large genomes of more than 2^32 bp in total length)'

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
    
def get_prunelevel_code_list():
    '''
    Get the code list of "prunelevel".
    '''

    return ['0', '1', '2', '3']

#-------------------------------------------------------------------------------
    
def get_prunelevel_code_list_text():
    '''
    Get the code list of "prunelevel" as text.
    '''

    return '0 (no pruning) or 1 (poor seqs) or 2 (repetitive seqs) or 3 (poor and repetitive)'

#-------------------------------------------------------------------------------
    
def get_gmap_output_format_code_list():
    '''
    Get the code list of "format".
    '''

    return ['COMPRESS', 'SUMMARY', 'ALIGN', 'PLS', 'GFF3_GENE', 'SPLICESITES', 'INTRONS', 'MAP_EXONS', 'MAP_RANGES', 'COORDS']

#-------------------------------------------------------------------------------
    
def get_gmap_output_format_code_list_text():
    '''
    Get the code list of "format" as text.
    '''

    return str(get_gmap_output_format_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
    
def get_quality_protocol_code_list():
    '''
    Get the code list of "quality_protocol".
    '''

    return ['ILLUMINA', 'SANGER']

#-------------------------------------------------------------------------------
    
def get_quality_protocol_code_list_text():
    '''
    Get the code list of "quality_protocol" as text.
    '''

    return 'ILLUMINA (ASCII 64-126) or SANGER (ASCII 33-126)'

#-------------------------------------------------------------------------------
    
def get_gsnap_output_format_code_list():
    '''
    Get the code list of "format".
    '''

    return ['SAM', 'M8']

#-------------------------------------------------------------------------------
    
def get_gsnap_output_format_code_list_text():
    '''
    Get the code list of "format" as text.
    '''

    return str(get_gsnap_output_format_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
     print('This file contains functions related to the GMAP-GSNAP process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
