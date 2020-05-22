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
This file contains functions related to the CD-HIT process used in
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

def create_cd_hit_est_config_file(experiment_id='exp001', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create CD-HIT-EST config file with the default options. It is necessary
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

    # create the CD-HIT-EST config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_cd_hit_est_config_file())):
            os.makedirs(os.path.dirname(get_cd_hit_est_config_file()))
        with open(get_cd_hit_est_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and assembly_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of CD-HIT-EST (CD-HIT package) and their meaning in http://weizhong-lab.ucsd.edu/cd-hit/.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "CD-HIT-EST parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of CD-HIT-EST and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --aS=0.9; --U=10'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the assembly result dataset.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_software = {0}'.format(assembly_software), '# assembly software: {0}'.format(get_assembly_software_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('assembly_dataset_id = {0}'.format(assembly_dataset_id), '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_type = {0}'.format(assembly_type), '# assembly type: CONTIGS or SCAFFOLDS in {0}; NONE in any other case'.format(xlib.get_soapdenovotrans_name())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the CD-HIT-EST parameters'))
            file_id.write( '{0}\n'.format('[CD-HIT-EST parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use; with 0, all CPUs will be used'))
            file_id.write( '{0:<50} {1}\n'.format('memory_limit = 0', '# memory limit (in MB) for the program; 0 for unlimitted'))
            file_id.write( '{0:<50} {1}\n'.format('seq_identity_threshold = 0.9', '# sequence identity threshold'))
            file_id.write( '{0:<50} {1}\n'.format('word_length = 8', '# word length'))
            file_id.write( '{0:<50} {1}\n'.format('mask = NX', '# masking letters (e.g. -mask NX, to mask out both "N" and "X")'))
            file_id.write( '{0:<50} {1}\n'.format('match = 2', '# matching score (1 for T-U and N-N)'))
            file_id.write( '{0:<50} {1}\n'.format('mismatch = -2', '# mismatching score'))
            file_id.write( '{0:<50} {1}\n'.format('other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_cd_hit_est_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_cd_hit_est_process(cluster_name, log, function=None):
    '''
    Run a CD-HIT-EST process.
    '''

    # initialize the control variable
    OK = True

    # get the CD-HIT-EST option dictionary
    cd_hit_est_option_dict = xlib.get_option_dict(get_cd_hit_est_config_file())

    # get the experiment identification
    experiment_id = cd_hit_est_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the CD-HIT-EST config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_cd_hit_est_name()))
    (OK, error_list) = check_cd_hit_est_config_file(strict=True)
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

    # check the CD-HIT is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_cd_hit_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_cd_hit_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_cd_hit_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_cd_hit_est_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the CD-HIT-EST process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_cd_hit_est_process_script()))
        (OK, error_list) = build_cd_hit_est_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the CD-HIT-EST process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(get_cd_hit_est_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_cd_hit_est_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the CD-HIT-EST process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the CD-HIT-EST process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_cd_hit_est_process_starter()))
        (OK, error_list) = build_cd_hit_est_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the CD-HIT-EST process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(get_cd_hit_est_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_cd_hit_est_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the CD-HIT-EST process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the CD-HIT-EST process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_cd_hit_est_process_starter()), log)

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

def check_cd_hit_est_config_file(strict):
    '''
    check the CD-HIT-EST config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        cd_hit_est_option_dict = xlib.get_option_dict(get_cd_hit_est_config_file())
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in cd_hit_est_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = cd_hit_est_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = cd_hit_est_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "assembly_software" has to be {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = cd_hit_est_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append('*** ERROR: the key "assembly_dataset_id" has to start with {0}.'.format(get_assembly_software_code_list_text()))
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = cd_hit_est_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append('*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {0} or NONE in any other case.'.format(xlib.get_soapdenovotrans_name()))
                    OK = False

        # check section "CD-HIT-EST parameters"
        if 'CD-HIT-EST parameters' not in sections_list:
            error_list.append('*** ERROR: the section "CD-HIT-EST parameters" is not found.')
            OK = False
        else:

            # check section "CD-HIT-EST parameters" - key "threads"
            threads = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=0):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "CD-HIT-EST parameters" - key "memory_limit"
            memory_limit = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('memory_limit', not_found)
            if memory_limit == not_found:
                error_list.append('*** ERROR: the key "memory_limit" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_int(memory_limit, minimum=0):
                error_list.append('*** ERROR: the key "memory_limit" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "CD-HIT-EST parameters" - key "seq_identity_threshold"
            seq_identity_threshold = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('seq_identity_threshold', not_found)
            if seq_identity_threshold == not_found:
                error_list.append('*** ERROR: the key "seq_identity_threshold" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_float(seq_identity_threshold, minimum=0., maximum=1.):
                error_list.append('*** ERROR: the key "seq_identity_threshold" has to be a float number between 0.0 and 1.0.')
                OK = False

            # check section "CD-HIT-EST parameters" - key "word_length"
            word_length = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('word_length', not_found)
            if word_length == not_found:
                error_list.append('*** ERROR: the key "word_length" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_int(word_length, minimum=1):
                    error_list.append('*** ERROR: the key "word_length" has to be an integer number greater than or equal to 1.')
                    OK = False

            # check section "CD-HIT-EST parameters" - key "mask"
            mask = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('mask', not_found)
            if mask == not_found:
                error_list.append('*** ERROR: the key "mask" is not found in the section "CD-HIT-EST parameters".')
                OK = False

            # check section "CD-HIT-EST parameters" - key "match"
            match = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('match', not_found)
            if match == not_found:
                error_list.append('*** ERROR: the key "match" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_int(match):
                error_list.append('*** ERROR: the key "match" has to be an integer number.')
                OK = False

            # check section "CD-HIT-EST parameters" - key "mismatch"
            mismatch = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('mismatch', not_found)
            if mismatch == not_found:
                error_list.append('*** ERROR: the key "mismatch" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif not xlib.check_int(mismatch):
                error_list.append('*** ERROR: the key "mismatch" has to be an integer number.')
                OK = False

            # check section "CD-HIT-EST parameters" - key "other_parameters"
            not_allowed_parameters_list = ['T', 'M', 'c', 'n', 'mask', 'match', 'mismatch']
            other_parameters = cd_hit_est_option_dict.get('CD-HIT-EST parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "CD-HIT-EST parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_cd_hit_est_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_cd_hit_est_process_script(cluster_name, current_run_dir):
    '''
    Build the current CD-HIT-EST process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the option dictionary
    cd_hit_est_option_dict = xlib.get_option_dict(get_cd_hit_est_config_file())

    # get the options
    experiment_id = cd_hit_est_option_dict['identification']['experiment_id']
    assembly_software = cd_hit_est_option_dict['identification']['assembly_software']
    assembly_dataset_id = cd_hit_est_option_dict['identification']['assembly_dataset_id']
    assembly_type = cd_hit_est_option_dict['identification']['assembly_type']
    threads = cd_hit_est_option_dict['CD-HIT-EST parameters']['threads']
    memory_limit = cd_hit_est_option_dict['CD-HIT-EST parameters']['memory_limit']
    seq_identity_threshold = cd_hit_est_option_dict['CD-HIT-EST parameters']['seq_identity_threshold']
    word_length = cd_hit_est_option_dict['CD-HIT-EST parameters']['word_length']
    mask = cd_hit_est_option_dict['CD-HIT-EST parameters']['mask']
    match = cd_hit_est_option_dict['CD-HIT-EST parameters']['match']
    mismatch = cd_hit_est_option_dict['CD-HIT-EST parameters']['mismatch']
    other_parameters = cd_hit_est_option_dict['CD-HIT-EST parameters']['other_parameters']

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
    if OK:
        output_file = '{0}/clustered-transcriptome.fasta'.format(current_run_dir)

    # write the CD-HIT-EST process script
    try:
        if not os.path.exists(os.path.dirname(get_cd_hit_est_process_script())):
            os.makedirs(os.path.dirname(get_cd_hit_est_process_script()))
        with open(get_cd_hit_est_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('CDHIT_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_cd_hit_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$CDHIT_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_cd_hit_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_cd_hit_est_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Running {0} process ..."'.format(xlib.get_cd_hit_est_name())))
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        cd-hit-est \\'))
            script_file_id.write( '{0}\n'.format('            -T {0} \\'.format(threads)))
            script_file_id.write( '{0}\n'.format('            -M {0} \\'.format(memory_limit)))
            script_file_id.write( '{0}\n'.format('            -i {0} \\'.format(transcriptome_file)))
            script_file_id.write( '{0}\n'.format('            -c {0} \\'.format(seq_identity_threshold)))
            script_file_id.write( '{0}\n'.format('            -n {0} \\'.format(word_length)))
            script_file_id.write( '{0}\n'.format('            -mask {0} \\'.format(mask)))
            script_file_id.write( '{0}\n'.format('            -match {0} \\'.format(match)))
            script_file_id.write( '{0}\n'.format('            -mismatch {0} \\'.format(mismatch)))
            if other_parameters.upper() == 'NONE':
                script_file_id.write( '{0}\n'.format('            -o {0}'.format(output_file)))
            else:
                script_file_id.write( '{0}\n'.format('            -o {0} \\'.format(output_file)))
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for i in range(len(parameter_list)):
                    if parameter_list[i].find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        if i < len(parameter_list) - 1:
                            script_file_id.write( '{0}\n'.format('            -{0} {1} \\'.format(parameter_name, parameter_value)))
                        else:
                            script_file_id.write( '{0}\n'.format('            -{0} {1}'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        if i < len(parameter_list):
                            script_file_id.write( '{0}\n'.format('            -{0} \\'.format(parameter_name)))
                        else:
                            script_file_id.write( '{0}\n'.format('            -{0}'.format(parameter_name)))
                    i += 1
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error cd-hit-est $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_cd_hit_est_name())))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_cd_hit_est_name())))
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
            script_file_id.write( '{0}\n'.format('run_cd_hit_est_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_cd_hit_est_process_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_cd_hit_est_process_starter(current_run_dir):
    '''
    Build the starter of the current CD-HIT-EST process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the CD-HIT-EST process starter
    try:
        if not os.path.exists(os.path.dirname(get_cd_hit_est_process_starter())):
            os.makedirs(os.path.dirname(get_cd_hit_est_process_starter()))
        with open(get_cd_hit_est_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_cd_hit_est_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_cd_hit_est_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_cd_hit_est_config_file():
    '''
    Get the CD-HIT-EST config file path.
    '''

    # assign the CD-HIT-EST config file path
    cd_hit_est_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_cd_hit_est_code())

    # return the CD-HIT-EST config file path
    return cd_hit_est_config_file

#-------------------------------------------------------------------------------

def get_cd_hit_est_process_script():
    '''
    Get the CD-HIT-EST process script path in the local computer.
    '''

    # assign the CD-HIT-EST script path
    cd_hit_est_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_cd_hit_est_code())

    # return the CD-HIT-EST script path
    return cd_hit_est_process_script

#-------------------------------------------------------------------------------

def get_cd_hit_est_process_starter():
    '''
    Get the CD-HIT-EST process starter path in the local computer.
    '''

    # assign the CD-HIT-EST process starter path
    cd_hit_est_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_cd_hit_est_code())

    # return the CD-HIT-EST starter path
    return cd_hit_est_process_starter

#-------------------------------------------------------------------------------
    
def get_assembly_software_code_list():
    '''
    Get the code list of "assembly_software".
    '''

    return [xlib.get_soapdenovotrans_code(), xlib.get_transabyss_code(), xlib.get_trinity_code(), xlib.get_ggtrinity_code(), xlib.get_cd_hit_est_code(),  xlib.get_transcript_filter_code()]

#-------------------------------------------------------------------------------
    
def get_assembly_software_code_list_text():
    '''
    Get the code list of "assembly_software".
    '''

    return '{0} ({1}) or {2} ({3}) or {4} ({5}) or {6} ({7}) or {8} ({9}) or {10} ({11})'.format(xlib.get_soapdenovotrans_code(), xlib.get_soapdenovotrans_name(), xlib.get_transabyss_code(), xlib.get_transabyss_name(), xlib.get_trinity_code(), xlib.get_trinity_name(), xlib.get_ggtrinity_code(), xlib.get_ggtrinity_name(), xlib.get_cd_hit_est_code(), xlib.get_cd_hit_est_name(), xlib.get_transcript_filter_code(), xlib.get_transcript_filter_name())

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the CD-HIT process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
