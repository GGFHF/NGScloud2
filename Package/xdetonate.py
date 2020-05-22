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
This file contains functions related to the DETONATE process used in
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

def create_rsem_eval_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq'], assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create RSEM-EVAL config file with the default options. It is necessary
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

    # create the RSEM-EVAL config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_rsem_eval_config_file())):
            os.makedirs(os.path.dirname(get_rsem_eval_config_file()))
        with open(get_rsem_eval_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, read_dataset_id and assembly_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of RSEM-EVAL (DETONATE package) and their meaning in http://deweylab.biostat.wisc.edu/detonate/.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the assembly result dataset.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_type = {0}'.format(assembly_type), '# assembly type: CONTIGS or SCAFFOLDS in {0}; NONE in any other case'.format(xlib.get_soapdenovotrans_name())))
            file_id.write( '\n')
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the RSEM-EVAL parameters'))
            file_id.write( '{0}\n'.format('[RSEM-EVAL parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('num_threads = 2', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('bowtie2_mismatch_rate = 0.1', '# maximum mismatch rate allowed (Bowtie 2 parameter)'))
            file_id.write( '{0:<50} {1}\n'.format('keep_intermediate_files = NO', '# keep temporary files generated: {0}'.format(get_keep_intermediate_file_code_list_text())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# read type: {0}'.format(get_read_type_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('length = 200', '# average read length in SE read type or average fragment length in PE read type'))
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
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_rsem_eval_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_rsem_eval_process(cluster_name, log, function=None):
    '''
    Run a RSEM-EVAL process.
    '''

    # initialize the control variable
    OK = True

    # get the RSEM-EVAL option dictionary
    rsem_eval_option_dict = xlib.get_option_dict(get_rsem_eval_config_file())

    # get the experiment identification
    experiment_id = rsem_eval_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the RSEM-EVAL config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_rsem_eval_name()))
    (OK, error_list) = check_rsem_eval_config_file(strict=True)
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

    # check the DETONATE is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_detonate_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_detonate_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_detonate_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_rsem_eval_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the RSEM-EVAL process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_rsem_eval_process_script()))
        (OK, error_list) = build_rsem_eval_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the RSEM-EVAL process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_rsem_eval_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_rsem_eval_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_rsem_eval_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RSEM-EVAL process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsem_eval_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_rsem_eval_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the RSEM-EVAL process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_rsem_eval_process_starter()))
        (OK, error_list) = build_rsem_eval_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the RSEM-EVAL process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_rsem_eval_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_rsem_eval_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_rsem_eval_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the RSEM-EVAL process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsem_eval_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_rsem_eval_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the RSEM-EVAL process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_rsem_eval_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_rsem_eval_process_starter()), log)

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

def check_rsem_eval_config_file(strict):
    '''
    Check the RSEM-EVAL config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        rsem_eval_option_dict = xlib.get_option_dict(get_rsem_eval_config_file())
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in rsem_eval_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = rsem_eval_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = rsem_eval_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = rsem_eval_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "assembly_software" has to be {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = rsem_eval_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "assembly_dataset_id" has to start with {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = rsem_eval_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {0} or NONE in any other case.'.format(xlib.get_soapdenovotrans_name()))
                    OK = False

        # check section "RSEM-EVAL parameters"
        if 'RSEM-EVAL parameters' not in sections_list:
            error_list.append('*** ERROR: the section "RSEM-EVAL parameters" is not found.')
            OK = False
        else:

            # check section "RSEM-EVAL parameters" - key "num_threads"
            num_threads = rsem_eval_option_dict.get('RSEM-EVAL parameters', {}).get('num_threads', not_found)
            if num_threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "RSEM-EVAL parameters".')
                OK = False
            elif not xlib.check_int(num_threads, minimum=1):
                error_list.append('*** ERROR: the key "num_threads" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "RSEM-EVAL parameters" - key "bowtie2_mismatch_rate"
            bowtie2_mismatch_rate = rsem_eval_option_dict.get('RSEM-EVAL parameters', {}).get('bowtie2_mismatch_rate', not_found)
            if bowtie2_mismatch_rate == not_found:
                error_list.append('*** ERROR: the key "bowtie2_mismatch_rate" is not found in the section "RSEM-EVAL parameters".')
                OK = False
            elif not xlib.check_float(bowtie2_mismatch_rate, minimum=0., maximum=1.):
                error_list.append('*** ERROR: the key "bowtie2_mismatch_rate" has to be a float number between 0.0 and 1.0.')
                OK = False

            # check section "RSEM-EVAL parameters" - key "keep_intermediate_files"
            keep_intermediate_files = rsem_eval_option_dict.get('RSEM-EVAL parameters', {}).get('keep_intermediate_files', not_found)
            if keep_intermediate_files == not_found:
                error_list.append('*** ERROR: the key "keep_intermediate_files" is not found in the section "RSEM-EVAL parameters".')
                OK = False
            elif not xlib.check_code(keep_intermediate_files, get_keep_intermediate_file_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "keep_intermediate_files" has to be {0}.'.format(get_keep_intermediate_file_code_list_text()))
                OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "format"
            format = rsem_eval_option_dict.get('library', {}).get('format', not_found)
            if format == not_found:
                error_list.append('*** ERROR: the key "format" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = rsem_eval_option_dict.get('library', {}).get('read_type', not_found)
            if read_type == not_found:
                error_list.append('*** ERROR: the key "read_type" is not found in the section "library".')
                OK = False
            elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "read_type" has to be {0}.'.format(get_read_type_code_list_text()))
                OK = False

            # check section "library" - key "length"
            length = rsem_eval_option_dict.get('library', {}).get('length', not_found)
            if length == not_found:
                error_list.append('*** ERROR: the key "length" is not found in the section "library".')
                OK = False
            elif not xlib.check_int(length, minimum=1):
                error_list.append('*** ERROR: the key "length" has to be an integer number greater than or equal to 1.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'RSEM-EVAL parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = rsem_eval_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = rsem_eval_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_rsem_eval_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_rsem_eval_process_script(cluster_name, current_run_dir):
    '''
    Build the current RSEM-EVAL process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    rsem_eval_option_dict = xlib.get_option_dict(get_rsem_eval_config_file())

    # get the options
    experiment_id = rsem_eval_option_dict['identification']['experiment_id']
    read_dataset_id = rsem_eval_option_dict['identification']['read_dataset_id']
    assembly_software = rsem_eval_option_dict['identification']['assembly_software']
    assembly_dataset_id = rsem_eval_option_dict['identification']['assembly_dataset_id']
    assembly_type = rsem_eval_option_dict['identification']['assembly_type']
    num_threads = rsem_eval_option_dict['RSEM-EVAL parameters']['num_threads']
    bowtie2_mismatch_rate = rsem_eval_option_dict['RSEM-EVAL parameters']['bowtie2_mismatch_rate']
    keep_intermediate_files = rsem_eval_option_dict['RSEM-EVAL parameters']['keep_intermediate_files']
    format = rsem_eval_option_dict['library']['format'].upper()
    read_type = rsem_eval_option_dict['library']['read_type'].upper()
    length = rsem_eval_option_dict['library']['length']

    # get the sections list
    sections_list = []
    for section in rsem_eval_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build library files
    files1 = ''
    files2 = ''
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^library-[0-9]+$', section):
            read_file_1 = rsem_eval_option_dict[section]['read_file_1']
            read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
            files1 += read_file_1 + ','
            if read_type == 'PE':
                read_file_2 = rsem_eval_option_dict[section]['read_file_2']
                read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                files2 += read_file_2 + ','
    files1 = files1[:len(files1) - 1]
    if read_type == 'PE':
        files2 = files2[:len(files2) - 1]

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

    # set the distribution file path
    distribution_file = '{0}/distribution.txt'.format(current_run_dir)

    # set the temporaly directory path
    temp_dir = '{0}/temp'.format(current_run_dir)

    # write the RSEM-EVAL process script
    try:
        if not os.path.exists(os.path.dirname(get_rsem_eval_process_script())):
            os.makedirs(os.path.dirname(get_rsem_eval_process_script()))
        with open(get_rsem_eval_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('DETONATE_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_detonate_bioconda_code())))
            script_file_id.write( '{0}\n'.format('BOWTIE2_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_bowtie2_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$DETONATE_PATH:$BOWTIE2_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_detonate_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_rsem_eval_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    rsem-eval-calculate-score --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Running rsem-eval-estimate-transcript-length-distribution ... "'))
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        rsem-eval-estimate-transcript-length-distribution \\'))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(transcriptome_file)))
            script_file_id.write( '{0}\n'.format('            {0}'.format(distribution_file)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rsem-eval-estimate-transcript-length-distribution $RC; fi'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Running rsem-eval-calculate-score ... "'))
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        rsem-eval-calculate-score \\'))
            script_file_id.write( '{0}\n'.format('            --num-threads {0} \\'.format(num_threads)))
            script_file_id.write( '{0}\n'.format('            --bowtie2 \\'))
            script_file_id.write( '{0}\n'.format('            --bowtie2-mismatch-rate {0} \\'.format(bowtie2_mismatch_rate)))
            script_file_id.write( '{0}\n'.format('            --transcript-length-parameters {0} \\'.format(distribution_file)))
            script_file_id.write( '{0}\n'.format('            --temporary-folder {0} \\'.format(temp_dir)))
            if keep_intermediate_files.upper() == 'YES':
                script_file_id.write( '{0}\n'.format('            --keep-intermediate-files \\'))
            if format == 'FASTA':
                script_file_id.write( '{0}\n'.format('            --no-qualities \\'))
            if read_type == 'PE':
                script_file_id.write( '{0}\n'.format('            --paired-end {0} {1} \\'.format(files1, files2)))
            else:
                script_file_id.write( '{0}\n'.format('            {0} \\'.format(files1)))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(transcriptome_file)))
            script_file_id.write( '{0}\n'.format('            {0} \\'.format(current_run_dir)))
            script_file_id.write( '{0}\n'.format('            {0}'.format(length)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then rsem-eval-calculate-score $RC; fi'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function move_result_files'))
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Moving result files ... "'))
            script_file_id.write( '{0}\n'.format('    mv ../{0}.stat/* .'.format(os.path.basename(current_run_dir))))
            script_file_id.write( '{0}\n'.format('    rm -fr ../{0}.stat'.format(os.path.basename(current_run_dir))))
            script_file_id.write( '{0}\n'.format('    mv ../{0}.* .'.format(os.path.basename(current_run_dir))))
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
            script_file_id.write( '{0}\n'.format('    echo "FILTERING_DATA - ASSEMBLY_SOFTWARE: {0}"'.format(assembly_software)))
            script_file_id.write( '{0}\n'.format('    echo "FILTERING_DATA - ASSEMBLY_DATASET_ID: {0}"'.format(assembly_dataset_id)))
            script_file_id.write( '{0}\n'.format('    echo "FILTERING_DATA - ASSEMBLY_TYPE: {0}"'.format(assembly_type)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    RECIPIENT={0}'.format(xconfiguration.get_contact_data())))
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_rsem_eval_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_rsem_eval_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_rsem_eval_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_rsem_eval_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_rsem_eval_process'))
            script_file_id.write( '{0}\n'.format('move_result_files'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_rsem_eval_process_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_rsem_eval_process_starter(current_run_dir):
    '''
    Build the starter of the current RSEM-EVAL process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the RSEM-EVAL process starter
    try:
        if not os.path.exists(os.path.dirname(get_rsem_eval_process_starter())):
            os.makedirs(os.path.dirname(get_rsem_eval_process_starter()))
        with open(get_rsem_eval_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_rsem_eval_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_rsem_eval_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_rsem_eval_config_file():
    '''
    Get the RSEM-EVAL config file path.
    '''

    # assign the RSEM-EVAL config file path
    rsem_eval_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_rsem_eval_code())

    # return the RSEM-EVAL config file path
    return rsem_eval_config_file

#-------------------------------------------------------------------------------

def get_rsem_eval_process_script():
    '''
    Get the RSEM-EVAL process script path in the local computer.
    '''

    # assign the RSEM-EVAL script path
    rsem_eval_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_rsem_eval_code())

    # return the RSEM-EVAL script path
    return rsem_eval_process_script

#-------------------------------------------------------------------------------

def get_rsem_eval_process_starter():
    '''
    Get the RSEM-EVAL process starter path in the local computer.
    '''

    # assign the RSEM-EVAL process starter path
    rsem_eval_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_rsem_eval_code())

    # return the RSEM-EVAL starter path
    return rsem_eval_process_starter

#-------------------------------------------------------------------------------

def create_ref_eval_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type = 'PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq'], assembly_dataset_id='sndt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create REF-EVAL config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    #...

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_ref_eval_process(cluster_name, log, function=None):
    '''
    Run a REF-EVAL process.
    '''

    # initialize the control variable
    OK = True

    # get the REF-EVAL option dictionary
    ref_eval_option_dict = xlib.get_option_dict(get_ref_eval_config_file())

    # get the experiment identification
    experiment_id = ref_eval_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the REF-EVAL config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_ref_eval_name()))
    (OK, error_list) = check_ref_eval_config_file(strict=True)
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

    # check the DETONATE is installed
    if OK:
        command = '[ -d {0}/{1} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_app_dir(), xlib.get_detonate_name())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_detonate_name()))
            OK = False

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_ref_eval_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the REF-EVAL process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_ref_eval_process_script()))
        (OK, error_list) = build_ref_eval_process_script(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the REF-EVAL process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_ref_eval_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ref_eval_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ref_eval_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the REF-EVAL process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ref_eval_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ref_eval_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the REF-EVAL process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ref_eval_process_starter()))
        (OK, error_list) = build_ref_eval_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the REF-EVAL process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_ref_eval_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ref_eval_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ref_eval_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the REF-EVAL process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ref_eval_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ref_eval_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the REF-EVAL process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ref_eval_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ref_eval_process_starter()), log)

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

def check_ref_eval_config_file(strict):
    '''
    Check the RSEM-EVAL config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # ...

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ref_eval_process_script(current_run_dir):
    '''
    Build the current RSEM-EVAL process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # ...

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ref_eval_process_starter(current_run_dir):
    '''
    Build the starter of the current RSEM-EVAL process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the RSEM-EVAL process starter
    try:
        if not os.path.exists(os.path.dirname(get_ref_eval_process_starter())):
            os.makedirs(os.path.dirname(get_ref_eval_process_starter()))
        with open(get_ref_eval_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_ref_eval_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ref_eval_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ref_eval_config_file():
    '''
    Get the REF-EVAL config file path.
    '''

    # assign the REF-EVAL config file path
    ref_eval_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_ref_eval_code())

    # return the REF-EVAL config file path
    return ref_eval_config_file

#-------------------------------------------------------------------------------

def get_ref_eval_process_script():
    '''
    Get the REF-EVAL process script path in the local computer.
    '''

    # assign the REF-EVAL script path
    ref_eval_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_ref_eval_code())

    # return the REF-EVAL script path
    return ref_eval_process_script

#-------------------------------------------------------------------------------

def get_ref_eval_process_starter():
    '''
    Get the REF-EVAL process starter path in the local computer.
    '''

    # assign the REF-EVAL process starter path
    ref_eval_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ref_eval_code())

    # return the REF-EVAL starter path
    return ref_eval_process_starter

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
    
def get_keep_intermediate_file_code_list():
    '''
    Get the code list of "keep_intermediate_files".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_keep_intermediate_file_code_list_text():
    '''
    Get the code list of "keep_intermediate_files" as text.
    '''

    return str(get_keep_intermediate_file_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the DETONATE process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
