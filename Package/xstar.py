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
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference and GTF files have to be located in the cluster directory {0}/experiment_id/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file, gtf_file and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of STAR and their meaning in https://github.com/alexdobin/STAR'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "STAR parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of STAR and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --outSAMattributes=All; --limitGenomeGenerateRAM=48000000000'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format('gtf_file = {0}'.format(gtf_file), '# GTF file name'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the STAR parameters'))
            file_id.write( '{0}\n'.format('[STAR parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('index_building = YES', '# index building when a reference is used: {0}'.format(get_index_building_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('two_pass_mode = NONE', '# 2-pass mapping mode: {0}'.format(get_two_pass_mode_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('two_pass_1_readsn = -1', '# number of reads to process for the 1st step; use -1 to map all reads in the first step'))
            file_id.write( '{0:<50} {1}\n'.format('out_filter_multimap_nmax = 20', '# maximun number of multiple alignments allowed for a read'))
            file_id.write( '{0:<50} {1}\n'.format('out_reads_unmapped = FASTX', '# output of unmapped and partially mapped reads in separate files: {0}'.format(get_out_reads_unmapped_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('cleanup = NO', '# leave intermediate files: {0}'.format(get_cleanup_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information of the library (only one library is allowed)'))
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
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_star_config_file()))
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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_star_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check STAR is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_star_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_star_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_star_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_star_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the STAR process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_star_process_script()))
        (OK, error_list) = build_star_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the STAR process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_star_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_star_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_star_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the STAR process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_star_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_star_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the STAR process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_star_process_starter()))
        (OK, error_list) = build_star_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the STAR process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_star_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_star_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_star_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the STAR process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_star_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_star_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the STAR process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_star_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
                error_list.append('*** ERROR: the key "index_building" has to be {0}.'.format(get_index_building_code_list_text()))
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
                error_list.append('*** ERROR: the key "two_pass_mode" has to be {0}.'.format(get_two_pass_mode_code_list_text()))
                OK = False

            # check section "STAR parameters" - key "two_pass_1_readsn"
            two_pass_1_readsn = star_option_dict.get('STAR parameters', {}).get('two_pass_1_readsn', not_found)
            if two_pass_1_readsn == not_found:
                error_list.append('*** ERROR: the key "two_pass_1_readsn" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1) and not xlib.check_int(threads, minimum=-1, maximum=-1):
                error_list.append('*** ERROR: the key "two_pass_1_readsn" has to be an integer number greater than or equal to 1, or -1 (all redas).')
                OK = False

            # check section "STAR parameters" - key "out_reads_unmapped"
            out_reads_unmapped = star_option_dict.get('STAR parameters', {}).get('out_reads_unmapped', not_found)
            if out_reads_unmapped == not_found:
                error_list.append('*** ERROR: the key "out_reads_unmapped" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(out_reads_unmapped, get_out_reads_unmapped_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "out_reads_unmapped" has to be {0}.'.format(get_out_reads_unmapped_code_list_text()))
                OK = False

            # check section "STAR parameters" - key "cleanup"
            cleanup = star_option_dict.get('STAR parameters', {}).get('cleanup', not_found)
            if cleanup == not_found:
                error_list.append('*** ERROR: the key "cleanup" is not found in the section "STAR parameters".')
                OK = False
            elif not xlib.check_code(cleanup, get_cleanup_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "cleanup" has to be {0}.'.format(get_cleanup_code_list_text()))
                OK = False

            # check section "STAR parameters" - key "other_parameters"
            not_allowed_parameters_list = ['runMode', 'runThreadN', 'genomeDir', 'readFilesCommand', 'readFilesIn', 'sjdbGTFfile', 'twopassMode', 'twopass1readsN', 'quantMode', 'outFilterMultimapNmax', 'outSAMunmapped', 'outFileNamePrefix', 'outTmpKeep', 'outSAMtype']
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = star_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'STAR parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = star_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = star_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_star_name()))

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
    (reference_file_name, reference_file_extension) = os.path.splitext(reference_file)
    star_indexes_dir = '{0}_star_indexes'.format(reference_file_name)

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
            script_file_id.write( '{0}\n'.format('STAR_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_star_bioconda_code())))
            script_file_id.write( '{0}\n'.format('TRINITY_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_trinity_bioconda_code())))
            script_file_id.write( '{0}\n'.format('BOWTIE2_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_bowtie2_bioconda_code())))
            script_file_id.write( '{0}\n'.format('export PATH=$STAR_PATH:$TRINITY_PATH:$BOWTIE2_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_star_bioconda_code())))
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
                script_file_id.write( '{0}\n'.format('function create_star_indexes'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Creating indexes ..."'))
                script_file_id.write( '{0}\n'.format('    rm -rf {0}'.format(star_indexes_dir)))
                script_file_id.write( '{0}\n'.format('    mkdir --parents {0}'.format(star_indexes_dir)))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        STAR \\'))
                script_file_id.write( '{0}\n'.format('            --runMode genomeGenerate \\'))
                script_file_id.write( '{0}\n'.format('            --runThreadN {0} \\'.format(threads)))
                script_file_id.write( '{0}\n'.format('            --genomeDir {0} \\'.format(star_indexes_dir)))
                script_file_id.write( '{0}\n'.format('            --genomeFastaFiles {0}'.format(reference_file)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error STAR $RC; fi'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_star_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    STAR --version'))
            for i in range(len(read_file_1_list)):
                # get the read file paths
                # set the name of read file 1
                if read_file_1.endswith('.gz'):
                    (read_file_1_name, read_file_1_extension) = os.path.splitext(os.path.basename(read_file_1_list[i][:-3]))
                    gz_compression = True
                else:
                    (read_file_1_name, read_file_1_extension) = os.path.splitext(os.path.basename(read_file_1_list[i]))
                    gz_compression = False
                # write the STAR run instructions
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Assembling reads ..."'))
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        STAR \\'))
                script_file_id.write( '{0}\n'.format('            --runMode alignReads \\'))
                script_file_id.write( '{0}\n'.format('            --runThreadN {0} \\'.format(threads)))
                script_file_id.write( '{0}\n'.format('            --genomeDir {0} \\'.format(star_indexes_dir)))
                if gz_compression:
                    script_file_id.write( '{0}\n'.format('            --readFilesCommand gzip -c \\'))
                if read_type.upper() == 'SE':
                    script_file_id.write( '{0}\n'.format('            --readFilesIn {0} \\'.format(read_file_1_list[i])))
                elif read_type.upper() == 'PE':
                    script_file_id.write( '{0}\n'.format('            --readFilesIn {0} {1} \\'.format(read_file_1_list[i], read_file_2_list[i])))
                script_file_id.write( '{0}\n'.format('            --sjdbGTFfile {0} \\'.format(gtf_file)))
                if two_pass_mode.upper() == 'NONE':
                    script_file_id.write( '{0}\n'.format('            --twopassMode None \\'))
                elif two_pass_mode.upper() == 'BASIC':
                    script_file_id.write( '{0}\n'.format('            --twopassMode Basic \\'))
                    script_file_id.write( '{0}\n'.format('            --twopass1readsN {0} \\'.format(two_pass_1_readsn)))
                script_file_id.write( '{0}\n'.format('            --quantMode TranscriptomeSAM \\'))
                if out_reads_unmapped.upper() == 'NONE':
                    script_file_id.write( '{0}\n'.format('            --outReadsUnmapped None \\'))
                elif out_reads_unmapped.upper() == 'FASTX':
                    script_file_id.write( '{0}\n'.format('            --outReadsUnmapped Fastx \\'))
                script_file_id.write( '{0}\n'.format('            --outFilterMultimapNmax {0} \\'.format(out_filter_multimap_nmax)))
                script_file_id.write( '{0}\n'.format('            --outSAMunmapped None \\'))
                script_file_id.write( '{0}\n'.format('            --outFileNamePrefix "{0}/{1}-" \\'.format(current_run_dir, read_file_1_name)))
                if cleanup.upper() == 'NO':
                    script_file_id.write( '{0}\n'.format('            --outTmpKeep Yes \\'))
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
                script_file_id.write( '{0}\n'.format('            --outSAMtype BAM SortedByCoordinate'))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error STAR $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_star_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_star_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_star_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_star_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('create_star_indexes'))
            script_file_id.write( '{0}\n'.format('run_star_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_star_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_star_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_star_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_star_config_file():
    '''
    Get the STAR config file path.
    '''

    # assign the STAR config file path
    star_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_star_code())

    # return the STAR config file path
    return star_config_file

#-------------------------------------------------------------------------------

def get_star_process_script():
    '''
    Get the STAR process script path in the local computer.
    '''

    # assign the STAR script path
    star_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_star_code())

    # return the STAR script path
    return star_process_script

#-------------------------------------------------------------------------------

def get_star_process_starter():
    '''
    Get the STAR process starter path in the local computer.
    '''

    # assign the STAR process starter path
    star_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_star_code())

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
