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
This file contains functions related to the TopHat process used in both
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

def create_tophat_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', annotation_file='Arabidopsis_thaliana.TAIR10.36.gtf', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create TopHat config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the TopHat config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_tophat_config_file())):
            os.makedirs(os.path.dirname(get_tophat_config_file()))
        with open(get_tophat_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file, gtf_file and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# You can consult the parameters of TopHat and their meaning in "http://ccb.jhu.edu/software/tophat/"'))
            file_id.write( '{0}\n'.format('# and the ones of Trinity in "https://github.com/trinityrnaseq/trinityrnaseq/wiki".'))
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# In section "TopHat parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of TopHat and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('#    other_parameters = --fusion-search; --fusion-read-mismatches=3'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_file = {reference_file}', '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format('annotation_file = {0}'.format(annotation_file), '# reference annotation file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the TopHat parameters'))
            file_id.write( '{0}\n'.format('[TopHat parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('read_mismatches  = 2', '# final read alignments having more than these many mismatches are discarded'))
            file_id.write( '{0:<50} {1}\n'.format('read_gap_length = 2', '# final read alignments having more than these many total length of gaps are discarded')) 
            file_id.write( '{0:<50} {1}\n'.format('read_edit_dist = 2', '# final read alignments having more than these many edit distance are discarded'))
            file_id.write( '{0:<50} {1}\n'.format('library_type = FR-UNSTRANDED', '# library type: {0}'.format(get_library_type_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '# This section has the global information of all libraries.\n')
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# {0}'.format(get_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# SE (single-end) or PE (paired-end)'))
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
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_tophat_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_tophat_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in TopHat config file.
    '''

    # initialize the control variable
    OK = True

    # get the TopHat option dictionary
    tophat_option_dict = xlib.get_option_dict(get_tophat_config_file())

    # get the experiment identification
    experiment_id = tophat_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the TopHat config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_tophat_name()))
    (OK, error_list) = check_tophat_config_file(strict=True)
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

    # check TopHat is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_tophat_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_tophat_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_tophat_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_tophat_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TopHat process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_tophat_process_script()))
        (OK, error_list) = build_tophat_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the TopHat process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(get_tophat_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_tophat_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_tophat_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TopHat process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_tophat_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_tophat_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the TopHat process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_tophat_process_starter()))
        (OK, error_list) = build_tophat_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the TopHat process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(get_tophat_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_tophat_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_tophat_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the TopHat process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_tophat_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_tophat_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the TopHat process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_tophat_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_tophat_process_starter()), log)

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

def check_tophat_config_file(strict):
    '''
    Check the TopHat configu file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        tophat_option_dict = xlib.get_option_dict(get_tophat_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in tophat_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = tophat_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = tophat_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = tophat_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "annotation_file"
            annotation_file = tophat_option_dict.get('identification', {}).get('annotation_file', not_found)
            if annotation_file == not_found:
                error_list.append('*** ERROR: the key "annotation_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = tophat_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "TopHat parameters"
        if 'TopHat parameters' not in sections_list:
            error_list.append('*** ERROR: the section "TopHat parameters" is not found.')
            OK = False
        else:

            # check section "TopHat parameters" - key "threads"
            threads = tophat_option_dict.get('TopHat parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "TopHat parameters" - key "read_mismatches"
            read_mismatches = tophat_option_dict.get('TopHat parameters', {}).get('read_mismatches', not_found)
            if read_mismatches == not_found:
                error_list.append('*** ERROR: the key "read_mismatches" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_int(read_mismatches, minimum=0):
                error_list.append('*** ERROR: the key "read_mismatches" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "TopHat parameters" - key "read_gap_length"
            read_gap_length = tophat_option_dict.get('TopHat parameters', {}).get('read_gap_length', not_found)
            if read_gap_length == not_found:
                error_list.append('*** ERROR: the key "read_gap_length" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_int(read_gap_length, minimum=0):
                error_list.append('*** ERROR: the key "read_gap_length" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "TopHat parameters" - key "read_edit_dist"
            read_edit_dist = tophat_option_dict.get('TopHat parameters', {}).get('read_edit_dist', not_found)
            if read_edit_dist == not_found:
                error_list.append('*** ERROR: the key "read_edit_dist" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_int(read_edit_dist, minimum=0):
                error_list.append('*** ERROR: the key "read_edit_dist" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "TopHat parameters" - key "library_type"
            library_type = tophat_option_dict.get('TopHat parameters', {}).get('library_type', not_found).lower()
            if library_type == not_found:
                error_list.append('*** ERROR: the key "library_type" is not found in the section "TopHat parameters".')
                OK = False
            elif not xlib.check_code(library_type, get_library_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "library_type" has to be {get_library_type_code_list_text()}.')
                OK = False

            # check section "TopHat parameters" - key "other_parameters"
            not_allowed_parameters_list = ['num-threads', 'read-mismatches', 'read-gap-length', 'read-edit-dist', 'library-type', 'GTF', 'output-dir', 'bowtie1']
            other_parameters = tophat_option_dict.get('TopHat parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "TopHat parameters".')
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
            format = tophat_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "format" has to be {get_format_code_list_text()}.')
                OK = False

            # check section "library" - key "read_type"
            read_type = tophat_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'TopHat parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = tophat_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_1" is not found in the section "{section}"')
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = tophat_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append(f'*** ERROR: the key "read_file_2" is not found in the section "{section}"')
                        OK = False

    # warn that the TopHat config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_tophat_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_tophat_process_script(cluster_name, current_run_dir):
    '''
    Build the current TopHat process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    tophat_option_dict = xlib.get_option_dict(get_tophat_config_file())

    # get the options
    experiment_id = tophat_option_dict['identification']['experiment_id']
    reference_dataset_id = tophat_option_dict['identification']['reference_dataset_id']
    reference_file = tophat_option_dict['identification']['reference_file']
    annotation_file = tophat_option_dict['identification']['annotation_file']
    read_dataset_id = tophat_option_dict['identification']['read_dataset_id']
    threads = tophat_option_dict['TopHat parameters']['threads']
    read_mismatches = tophat_option_dict['TopHat parameters']['read_mismatches']
    read_gap_length = tophat_option_dict['TopHat parameters']['read_gap_length']
    read_edit_dist = tophat_option_dict['TopHat parameters']['read_edit_dist']
    library_type = tophat_option_dict['TopHat parameters']['library_type']
    other_parameters = tophat_option_dict['TopHat parameters']['other_parameters']
    format = tophat_option_dict['library']['format'].upper()
    read_type = tophat_option_dict['library']['read_type']

    # get the sections list
    sections_list = []
    for section in tophat_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build library files
    files1 = ''
    files2 = ''
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = tophat_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            files1 += read_file_1 + ','
            if read_type.upper() == 'PE':
                read_file_2 = tophat_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                files2 += read_file_2 + ','
    files1 = files1[:len(files1) - 1]
    if read_type.upper() == 'PE':
        files2 = files2[:len(files2) - 1]

    # set the reference file path
    reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

    # set the annotation file path
    annotation_file = xlib.get_cluster_reference_file(reference_dataset_id, annotation_file)

    # set the Bowtie2 index dir and index base
    (reference_file_name, reference_file_extension) = os.path.splitext(reference_file)
    bowtie2_index_dir = '{0}_bowtie2_indexes'.format(reference_file_name)
    bowtie2_index_base = '{0}/{1}'.format(bowtie2_index_dir, os.path.basename(reference_file_name))

    # write the TopHat process script
    try:
        if not os.path.exists(os.path.dirname(get_tophat_process_script())):
            os.makedirs(os.path.dirname(get_tophat_process_script()))
        with open(get_tophat_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('TOPHAT_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_tophat_anaconda_code())))
            script_file_id.write( '{0}\n'.format('TRINITY_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_trinity_anaconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$TOPHAT_PATH:$TRINITY_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_tophat_anaconda_code())))
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
            script_file_id.write( '{0}\n'.format('function create_bowtie2_indexes'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Creating Bowtie2 indexes ..."'))
            script_file_id.write( '{0}\n'.format('    rm -rf {0}'.format(bowtie2_index_dir)))
            script_file_id.write( '{0}\n'.format('    mkdir --parents {0}'.format(bowtie2_index_dir)))
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('        bowtie2-build \\'))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(reference_file)))
            script_file_id.write( '{0}\n'.format('            {0}'.format(bowtie2_index_base)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error bowtie2 $RC; fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_tophat_process'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    tophat --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('        tophat \\'))
            script_file_id.write( '{0}\n'.format('            --num-threads {0} \\'.format(threads)))
            script_file_id.write( '{0}\n'.format('            --read-mismatches {0} \\'.format(read_mismatches)))
            script_file_id.write( '{0}\n'.format('            --read-gap-length {0} \\'.format(read_gap_length)))
            script_file_id.write( '{0}\n'.format('            --read-edit-dist {0} \\'.format(read_edit_dist)))
            script_file_id.write( '{0}\n'.format('            --library-type {0} \\'.format(library_type)))
            script_file_id.write( '{0}\n'.format('            --GTF {0} \\'.format(annotation_file)))
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
            script_file_id.write( '{0}\n'.format('            --output-dir {0} \\'.format(current_run_dir)))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(bowtie2_index_base)))
            if read_type.upper() == 'PE':
                script_file_id.write( '{0}\n'.format('            {0} \\'.format(files1)))
                script_file_id.write( '{0}\n'.format('            {0}'.format(files2)))
            else:
                script_file_id.write( '{0}\n'.format('            {0}'.format(files1)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error tophat $RC; fi'))
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
            process_name = f'{xlib.get_tophat_name()} process'
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
            script_file_id.write( '{0}\n'.format('create_bowtie2_indexes'))
            script_file_id.write( '{0}\n'.format('run_tophat_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_tophat_process_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_tophat_process_starter(current_run_dir):
    '''
    Build the starter of the current TopHat process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the TopHat process starter
    try:
        if not os.path.exists(os.path.dirname(get_tophat_process_starter())):
            os.makedirs(os.path.dirname(get_tophat_process_starter()))
        with open(get_tophat_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_tophat_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_tophat_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_tophat_config_file():
    '''
    Get the TopHat config file path.
    '''

    # assign the TopHat config file path
    tophat_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_tophat_code())

    # return the TopHat config file path
    return tophat_config_file

#-------------------------------------------------------------------------------

def get_tophat_process_script():
    '''
    Get the TopHat process script path in the local computer.
    '''

    # assign the TopHat script path
    tophat_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_tophat_code())

    # return the TopHat script path
    return tophat_process_script

#-------------------------------------------------------------------------------

def get_tophat_process_starter():
    '''
    Get the TopHat process starter path in the local computer.
    '''

    # assign the TopHat process starter path
    tophat_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_tophat_code())

    # return the TopHat starter path
    return tophat_process_starter

#-------------------------------------------------------------------------------
    
def get_library_type_code_list():
    '''
    Get the code list of "library_type".
    '''

    return ['FR-FIRSTSTRAND', 'FR-SECONDSTRAND', 'FR-UNSTRANDED']

#-------------------------------------------------------------------------------
    
def get_library_type_code_list_text():
    '''
    Get the code list of "library_type" as text.
    '''

    return str(get_library_type_code_list()).strip('[]').replace('\'','').replace(',', ' or')

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
     print('This file contains functions related to the TopHat process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
