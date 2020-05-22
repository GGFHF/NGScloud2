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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference file has to be located in the cluster directory {0}/experiment_id/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file and assembly_dataset_id are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of GMAP and their meaning in http://research-pub.gene.com/gmap/.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "GMAP parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of GMAP and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --no-chimeras; --canonical-mode=2'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification or NONE'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_type = {0}'.format(assembly_type), '# assembly type: CONTIGS or SCAFFOLDS in {0}; NONE in any other case'.format(xlib.get_soapdenovotrans_name())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the GMAP parameters'))
            file_id.write( '{0}\n'.format('[GMAP parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('index_building = YES', '# index building: {0}'.format(get_index_building_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('gmap_version = gmapl', '# GMAP version: {0}'.format(get_gmap_version_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = NONE', '# kmer size to use in genome database or NONE (the program will find the highest available kmer size in the genome database)'))
            file_id.write( '{0:<50} {1}\n'.format('sampling = NONE', '# Sampling to use in genome database or NONE (the program will find the smallest available sampling value in the genome database within selected k-mer size)'))
            file_id.write( '{0:<50} {1}\n'.format('input-buffer-size = 1000', '# size of input buffer'))
            file_id.write( '{0:<50} {1}\n'.format('output-buffer-size = 1000', '# size of buffer size in queries for output thread'))
            file_id.write( '{0:<50} {1}\n'.format('prunelevel = 0', '# pruning level: {0}'.format(get_prunelevel_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('format = COMPRESS', '# format for output: {0}'.format(get_gmap_output_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_gmap_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_gmap_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the GMAP-GSNAP is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_gmap_gsnap_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_gmap_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_gmap_name()))

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
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GMAP process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_gmap_process_script()))
        (OK, error_list) = build_gmap_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the GMAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} in the directory {1} of the master ...\n'.format(get_gmap_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_gmap_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_gmap_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GMAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gmap_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_gmap_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GMAP process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_gmap_process_starter()))
        (OK, error_list) = build_gmap_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the GMAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_gmap_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_gmap_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_gmap_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GMAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gmap_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_gmap_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the GMAP process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gmap_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
                error_list.append('*** ERROR: the key "assembly_software" has to be {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = gmap_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "assembly_dataset_id" hast to start with {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = gmap_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {0} or NONE in any other case.'.format(xlib.get_soapdenovotrans_name()))
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
                error_list.append('*** ERROR: the key "index_building" has to be {0}.'.format(get_index_building_code_list_text()))
                OK = False

            # check section "GMAP parameters" - key "version"
            gmap_version = gmap_option_dict.get('GMAP parameters', {}).get('gmap_version', not_found)
            if gmap_version == not_found:
                error_list.append('*** ERROR: the key "gmap_version" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(gmap_version, get_gmap_version_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "gmap_version" has to be {0}.'.format(get_gmap_version_code_list_text()))
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
                error_list.append('*** ERROR: the key "prunelevel" has to be {0}.'.format(get_prunelevel_code_list_text()))
                OK = False

            # check section "GMAP parameters" - key "format"
            format = gmap_option_dict.get('GMAP parameters', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "GMAP parameters".')
                OK = False
            elif not xlib.check_code(format, get_gmap_output_format_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_gmap_output_format_code_list_text()))
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
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_gmap_name()))

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
    gmap_database = '{0}-gmap_database'.format(reference_file_name)

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
    output_file = 'gmap_output_{0}.txt'.format(format.lower())

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
            script_file_id.write( '{0}\n'.format('GMAP_GSNAP_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_gmap_gsnap_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$GMAP_GSNAP_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_gmap_gsnap_bioconda_code())))
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
            if index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function build_gmap_database'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        gmap_build \\'))
                script_file_id.write( '{0}\n'.format('            --dir={0} \\'.format(cluster_reference_dataset_dir)))
                script_file_id.write( '{0}\n'.format('            --db={0} \\'.format(gmap_database)))
                if kmer.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --kmer={0} \\'.format(kmer)))
                script_file_id.write( '{0}\n'.format('            {0}'.format(cluster_reference_file)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gmap_build $RC; fi'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_gmap_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    gmap --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        {0} \\'.format(gmap_version.lower())))
            script_file_id.write( '{0}\n'.format('            --nthreads={0} \\'.format(threads)))
            script_file_id.write( '{0}\n'.format('            --dir={0} \\'.format(cluster_reference_dataset_dir)))
            script_file_id.write( '{0}\n'.format('            --db={0} \\'.format(gmap_database)))
            if kmer.upper() != 'NONE':
                script_file_id.write( '{0}\n'.format('            --kmer={0} \\'.format(kmer)))
            if sampling.upper() != 'NONE':
                script_file_id.write( '{0}\n'.format('            --sampling={0} \\'.format(sampling)))
            script_file_id.write( '{0}\n'.format('            --input-buffer-size={0} \\'.format(input_buffer_size)))
            script_file_id.write( '{0}\n'.format('            --output-buffer-size={0} \\'.format(output_buffer_size)))
            script_file_id.write( '{0}\n'.format('            --prunelevel={0} \\'.format(prunelevel)))
            if format.upper() == 'COMPRESS':
                script_file_id.write( '{0}\n'.format('            --compress \\'))
            elif format.upper() == 'SUMMARY':
                script_file_id.write( '{0}\n'.format('            --summary \\'))
            elif format.upper() == 'ALIGN':
                script_file_id.write( '{0}\n'.format('            --align \\'))
            else:
                script_file_id.write( '{0}\n'.format('            --format={0} \\'.format(format.lower())))
            script_file_id.write( '{0}\n'.format('            --ordered \\'))
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for i in range(len(parameter_list)):
                    if parameter_list[i].find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write( '{0}\n'.format('            --{0}={1} \\'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(transcriptome_file)))
            script_file_id.write( '{0}\n'.format('            > {0}'.format(output_file)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gmap $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_gmap_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_gmap_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_gmap_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_gmap_name(), cluster_name))))
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
            if index_building.upper() == 'YES':
                script_file_id.write( '{0}\n'.format('build_gmap_database'))
            script_file_id.write( '{0}\n'.format('run_gmap_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gmap_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_gmap_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gmap_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gmap_config_file():
    '''
    Get the GMAP config file path.
    '''

    # assign the GMAP config file path
    gmap_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_gmap_code())

    # return the GMAP config file path
    return gmap_config_file

#-------------------------------------------------------------------------------

def get_gmap_process_script():
    '''
    Get the GMAP process script path in the local computer.
    '''

    # assign the GMAP script path
    gmap_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_gmap_code())

    # return the GMAP script path
    return gmap_process_script

#-------------------------------------------------------------------------------

def get_gmap_process_starter():
    '''
    Get the GMAP process starter path in the local computer.
    '''

    # assign the GMAP process starter path
    gmap_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_gmap_code())

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
    elif assembly_dataset_id.upper() == 'NONE':
        assembly_software = 'NONE'

    # create the GSNAP config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_gsnap_config_file())):
            os.makedirs(os.path.dirname(get_gsnap_config_file()))
        with open(get_gsnap_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference file has to be located in the cluster directory {0}/experiment_id/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file, assembly_dataset_id and read_dataset_id are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of GSNAP and their meaning in http://research-pub.gene.com/gmap/.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "GSNAP parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of GSNAP and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --allow-pe-name-mismatch; --max-gmap-improvement=5'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name or NONE if an assembly is used'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}; or NONE if a reference is used'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification or NONE if a reference is used'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_type = {0}'.format(assembly_type), '# assembly type: CONTIGS or SCAFFOLDS in {0}; NONE in any other case'.format(xlib.get_soapdenovotrans_name())))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the GSNAP parameters'))
            file_id.write( '{0}\n'.format('[GSNAP parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('index_building = YES', '# index building when a reference is used: {0}'.format(get_index_building_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('gsnap_version = gsnapl', '# GSNAP version: {0}'.format(get_gsnap_version_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = NONE', '# kmer size to use in genome database or NONE (the program will find the highest available kmer size in the genome database)'))
            file_id.write( '{0:<50} {1}\n'.format('sampling = NONE', '# Sampling to use in genome database or NONE (the program will find the smallest available sampling value in the genome database within selected k-mer size)'))
            file_id.write( '{0:<50} {1}\n'.format('input-buffer-size = 1000', '# size of input buffer'))
            file_id.write( '{0:<50} {1}\n'.format('output-buffer-size = 1000', '# size of buffer size in queries for output thread'))
            file_id.write( '{0:<50} {1}\n'.format('max-mismatches = NONE', '# maximum number of mismatches allowed or NONE (the program will calculate the ultrafast level)'))
            file_id.write( '{0:<50} {1}\n'.format('indel-endlength = 4', '# minimum length at end required for indel alignments'))
            file_id.write( '{0:<50} {1}\n'.format('orientation = FR', '# orientation of paired-end reads: {0}'.format(get_orientation_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('quality-zero-score = 33', '# FASTQ quality scores are zero at this ASCII value'))
            file_id.write( '{0:<50} {1}\n'.format('quality-print-shift = 0', '# shift FASTQ quality scores by this amount in output'))
            file_id.write( '{0:<50} {1}\n'.format('format = SAM', '# format for output: {0}'.format(get_gsnap_output_format_code_list_text())))
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
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_gsnap_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_gsnap_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the GMAP-GSNAP is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_gmap_gsnap_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_gmap_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_gsnap_name()))

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
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GSNAP process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_gsnap_process_script()))
        (OK, error_list) = build_gsnap_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the GSNAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} in the directory {1} of the master ...\n'.format(get_gsnap_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_gsnap_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_gsnap_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GSNAP process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gsnap_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_gsnap_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the GSNAP process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_gsnap_process_starter()))
        (OK, error_list) = build_gsnap_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the GSNAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_gsnap_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_gsnap_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_gsnap_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the GSNAP process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gsnap_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_gsnap_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the GSNAP process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_gsnap_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
            elif assembly_software.upper() != 'NONE' and not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "assembly_software" has to be {0}; or NONE if a reference is used.'.format(get_assembly_software_code_list_text()))
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
            elif assembly_dataset_id.upper() != 'NONE' and not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "assembly_dataset_id" does not have to start with {0}.'.format(get_assembly_software_code_list_text()))
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
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {0} or NONE in any other case.'.format(xlib.get_soapdenovotrans_name()))
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
                error_list.append('*** ERROR: the key "index_building" has to be {0}.'.format(get_index_building_code_list_text()))
                OK = False

            # check section "GSNAP parameters" - key "version"
            gsnap_version = gsnap_option_dict.get('GSNAP parameters', {}).get('gsnap_version', not_found)
            if gsnap_version == not_found:
                error_list.append('*** ERROR: the key "gsnap_version" is not found in the section "GSNAP parameters".')
                OK = False
            elif not xlib.check_code(gsnap_version, get_gsnap_version_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "gsnap_version" has to be {0}.'.format(get_gsnap_version_code_list_text()))
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
            elif  max_mismatches.upper() != 'NONE' and not xlib.check_float(max_mismatches, minimum=0., maximum=1., mne=0., mxe=0.):
                error_list.append('*** ERROR: the key "max-mismatches" has to be a float numberbetween 0.0 and 1.0 or NONE.')
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
                error_list.append('*** ERROR: the key "orientation" has to be {0}.'.format(get_orientation_code_list_text()))
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_gsnap_output_format_code_list_text()))
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = gsnap_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'GSNAP parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = gsnap_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = gsnap_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_gsnap_name()))

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
        cluster_reference_dataset_dir = current_run_dir

    # set the cluster reference file
    if reference_dataset_id.upper() != 'NONE':
        cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    else:
        if assembly_software == xlib.get_soapdenovotrans_code():
            if assembly_type == 'CONTIGS':
                cluster_reference_file = '{0}/{1}-{2}.contig'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id), experiment_id, assembly_dataset_id)
            elif  assembly_type == 'SCAFFOLDS':
                cluster_reference_file = '{0}/{1}-{2}.scafSeq'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id), experiment_id, assembly_dataset_id)
        elif assembly_software == xlib.get_transabyss_code():
            cluster_reference_file = '{0}/transabyss-final.fa'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
        elif assembly_software == xlib.get_trinity_code():
            cluster_reference_file = '{0}/Trinity.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
        elif assembly_software == xlib.get_ggtrinity_code():
            cluster_reference_file = '{0}/Trinity-GG.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
        elif assembly_software == xlib.get_cd_hit_est_code():
            cluster_reference_file = '{0}/clustered-transcriptome.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))
        elif assembly_software == xlib.get_transcript_filter_code():
            cluster_reference_file = '{0}/filtered-transcriptome.fasta'.format(xlib.get_cluster_experiment_result_dataset_dir(experiment_id, assembly_dataset_id))

    # set the GMAP database name
    # -- gmap_database = 'gmap_database'
    reference_file_name, reference_file_extension = os.path.splitext(reference_file)
    gmap_database = '{0}-gmap_database'.format(reference_file_name)

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
            script_file_id.write( '{0}\n'.format('GMAP_GSNAP_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_gmap_gsnap_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$GMAP_GSNAP_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_gmap_gsnap_bioconda_code())))
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
            if reference_dataset_id.upper() != 'NONE' and index_building.upper() == 'YES':
                script_file_id.write( '#-------------------------------------------------------------------------------\n')
                script_file_id.write( '{0}\n'.format('function build_gmap_database'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        gmap_build \\'))
                script_file_id.write( '{0}\n'.format('            --dir={0} \\'.format(cluster_reference_dataset_dir)))
                script_file_id.write( '{0}\n'.format('            --db={0} \\'.format(gmap_database)))
                if kmer.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --kmer={0} \\'.format(kmer)))
                script_file_id.write( '{0}\n'.format('            {0}'.format(cluster_reference_file)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gmap_build $RC; fi'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_gsnap_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    gsnap --version'))
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
                    basename = '{0}/{1}'.format(current_run_dir, os.path.basename(basename[:len(basename)-position-1]))
                else:
                    basename = '{0}/{1}'.format(current_run_dir, os.path.basename(basename))
                split_output = '{0}-split.txt'.format(basename)
                failed_input = '{0}-failed.txt'.format(basename)
                # write the gsnap run instructions
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        {0} \\'.format(gsnap_version.lower())))
                script_file_id.write( '{0}\n'.format('            --nthreads={0} \\'.format(threads)))
                if reference_dataset_id.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --dir={0} \\'.format(cluster_reference_dataset_dir)))
                    script_file_id.write( '{0}\n'.format('            --db={0} \\'.format(gmap_database)))
                if kmer.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --kmer={0} \\'.format(kmer)))
                if sampling.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --sampling={0} \\'.format(sampling)))
                script_file_id.write( '{0}\n'.format('            --input-buffer-size={0} \\'.format(input_buffer_size)))
                script_file_id.write( '{0}\n'.format('            --output-buffer-size={0} \\'.format(output_buffer_size)))
                if max_mismatches.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --max-mismatches={0} \\'.format(max_mismatches)))
                script_file_id.write( '{0}\n'.format('            --indel-endlength={0} \\'.format(indel_endlength)))
                if read_type.upper() == 'PE':
                    script_file_id.write( '{0}\n'.format('            --orientation={0} \\'.format(orientation.upper())))
                script_file_id.write( '{0}\n'.format('            --quality-zero-score={0} \\'.format(quality_zero_score)))
                script_file_id.write( '{0}\n'.format('            --quality-print-shift={0} \\'.format(quality_print_shift)))
                script_file_id.write( '{0}\n'.format('            --format={0} \\'.format(format.lower())))
                script_file_id.write( '{0}\n'.format('            --ordered \\'))
                script_file_id.write( '{0}\n'.format('            --split-output={0} \\'.format(split_output)))
                script_file_id.write( '{0}\n'.format('            --failed-input={0} \\'.format(failed_input)))
                if other_parameters.upper() != 'NONE':
                    parameter_list = [x.strip() for x in other_parameters.split(';')]
                    for i in range(len(parameter_list)):
                        if parameter_list[i].find('=') > 0:
                            pattern = r'^--(.+)=(.+)$'
                            mo = re.search(pattern, parameter_list[i])
                            parameter_name = mo.group(1).strip()
                            parameter_value = mo.group(2).strip()
                            script_file_id.write( '{0}\n'.format('            --{0}={1} \\'.format(parameter_name, parameter_value)))
                        else:
                            pattern = r'^--(.+)$'
                            mo = re.search(pattern, parameter_list[i])
                            parameter_name = mo.group(1).strip()
                            script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
                if gunzip:
                    script_file_id.write( '{0}\n'.format('            --gunzip \\'))
                if bunzip2:
                    script_file_id.write( '{0}\n'.format('            --bunzip2 \\'))
                if read_type.upper() == 'SE':
                    script_file_id.write( '{0}\n'.format('            {0}'.format(read_file_1_list[i])))
                if read_type.upper() == 'PE':
                    script_file_id.write( '{0}\n'.format('            {0} \\'.format(read_file_1_list[i])))
                    script_file_id.write( '{0}\n'.format('            {0}'.format(read_file_2_list[i])))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gsnap $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_gsnap_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_gsnap_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_gsnap_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_gsnap_name(), cluster_name))))
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
            if reference_dataset_id.upper() != 'NONE' and index_building.upper() == 'YES':
                script_file_id.write( '{0}\n'.format('build_gmap_database'))
            script_file_id.write( '{0}\n'.format('run_gsnap_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gsnap_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_gsnap_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gsnap_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gsnap_config_file():
    '''
    Get the GSNAP config file path.
    '''

    # assign the GSNAP config file path
    gsnap_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_gsnap_code())

    # return the GSNAP config file path
    return gsnap_config_file

#-------------------------------------------------------------------------------

def get_gsnap_process_script():
    '''
    Get the GSNAP process script path in the local computer.
    '''

    # assign the GSNAP script path
    gsnap_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_gsnap_code())

    # return the GSNAP script path
    return gsnap_process_script

#-------------------------------------------------------------------------------

def get_gsnap_process_starter():
    '''
    Get the GSNAP process starter path in the local computer.
    '''

    # assign the GSNAP process starter path
    gsnap_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_gsnap_code())

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

    return '{0} ({1}) or {2} ({3}) or {4} ({5}) or {6} ({7}) or {8} ({9}) or {10} ({11})'.format(xlib.get_soapdenovotrans_code(), xlib.get_soapdenovotrans_name(), xlib.get_transabyss_code(), xlib.get_transabyss_name(), xlib.get_trinity_code(), xlib.get_trinity_name(), xlib.get_ggtrinity_code(), xlib.get_ggtrinity_name(), xlib.get_cd_hit_est_code(), xlib.get_cd_hit_est_name(), xlib.get_transcript_filter_code(), xlib.get_transcript_filter_name())

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
