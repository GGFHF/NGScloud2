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
This file contains functions related to the eXpress process used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys
import urllib

import xbioinfoapp
import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_express_config_file(experiment_id='exp001', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS', alignment_dataset_id='bowtie2-170101-235959'):
    '''
    Create eXpress config file with the default options. It is necessary
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

    # set the alignment software
    if alignment_dataset_id.startswith(xlib.get_bowtie2_code()):
        alignment_software = xlib.get_bowtie2_code()

    # create the eXpress config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_express_config_file())):
            os.makedirs(os.path.dirname(get_express_config_file()))
        with open(get_express_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write(f'# The alignment file has to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/alignment_dataset_id\n')
            file_id.write( '# The experiment_id, assembly_dataset_id and alignment_dataset_id are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of eXpress and their meaning in "https://pachterlab.github.io/eXpress/".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "eXpress parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of Cufflinks and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --calc-covar; --forget-param=0.6\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_software = {alignment_software}', f'# alignment software: {get_alignment_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'alignment_dataset_id = {alignment_dataset_id}', '# alignment dataset identification'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the eXpress parameters\n')
            file_id.write( '[eXpress parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'frag-len-mean = 200', '# mean fragment length'))
            file_id.write( '{0:<50} {1}\n'.format( 'frag-len-stddev = 60', '# fragment length standard deviation'))
            file_id.write( '{0:<50} {1}\n'.format( 'library_type = NONE', f'# library type: {get_library_type_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'max-indel-size = 0', '# maximum allowed size of a single indel'))
            file_id.write( '{0:<50} {1}\n'.format( 'no-bias-correct = NO', f'# if YES, eXpress will not measure and account for sequence-specific biases: {get_no_bias_correct_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'no-error-model = NO', f'# if YES, eXpress will not measure and account for errors in alignments: {get_no_error_model_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_express_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_express_process(cluster_name, log, function=None):
    '''
    Run a eXpress process.
    '''

    # initialize the control variable
    OK = True

    # get the eXpress option dictionary
    express_option_dict = xlib.get_option_dict(get_express_config_file())

    # get the experiment identification
    experiment_id = express_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the eXpress config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_express_name()} config file ...\n')
    (OK, error_list) = check_express_config_file(strict=True)
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

    # check eXpress is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_express_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_express_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_express_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_express_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the eXpress process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_express_process_script()} ...\n')
        (OK, error_list) = build_express_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the eXpress process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_express_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_express_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_express_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the eXpress process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_express_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_express_process_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the eXpress process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_express_process_starter()} ...\n')
        (OK, error_list) = build_express_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the eXpress process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_express_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_express_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_express_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the eXpress process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_express_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_express_process_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the eXpress process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_express_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_express_process_starter()), log)

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

def check_express_config_file(strict):
    '''
    Check the eXpress config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        express_option_dict = xlib.get_option_dict(get_express_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in express_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = express_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = express_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = express_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" has to start with {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = express_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
                    OK = False

            # check section "identification" - key "alignment_software"
            alignment_software = express_option_dict.get('identification', {}).get('alignment_software', not_found)
            if alignment_software == not_found:
                error_list.append(f'*** ERROR: the key "alignment_software" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_code(alignment_software, get_alignment_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "alignment_software" has to be {get_alignment_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "alignment_dataset_id"
            alignment_dataset_id = express_option_dict.get('identification', {}).get('alignment_dataset_id', not_found)
            if alignment_dataset_id == not_found:
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" is not found in the section "{section}".')
                OK = False
            elif not xlib.check_startswith(alignment_dataset_id, get_alignment_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "alignment_dataset_id" has to start with {get_alignment_software_code_list_text()}.')
                OK = False

        # check section "eXpress parameters"
        if 'eXpress parameters' not in sections_list:
            error_list.append('*** ERROR: the section "eXpress parameters" is not found.')
            OK = False
        else:

            # check section "express parameters" - key "frag-len-mean"
            frag_len_mean = express_option_dict.get('eXpress parameters', {}).get('frag-len-mean', not_found)
            if frag_len_mean == not_found:
                error_list.append('*** ERROR: the key "frag-len-mean" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_int(frag_len_mean, minimum=1):
                error_list.append('*** ERROR: the key "frag-len-mean" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "express parameters" - key "frag-len-stddev"
            frag_len_stddev = express_option_dict.get('eXpress parameters', {}).get('frag-len-stddev', not_found)
            if frag_len_stddev == not_found:
                error_list.append('*** ERROR: the key "frag-len-stddev" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_int(frag_len_stddev, minimum=1):
                error_list.append('*** ERROR: the key "frag-len-stddev" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "eXpress parameters" - key "library_type"
            library_type = express_option_dict.get('eXpress parameters', {}).get('library_type', not_found)
            if library_type == not_found:
                error_list.append('*** ERROR: the key "library_type" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_code(library_type, get_library_type_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "library_type" has to be {get_library_type_code_list_text()}.')
                OK = False

            # check section "eXpress parameters" - key "max-indel-size"
            max_indel_size = express_option_dict.get('eXpress parameters', {}).get('max-indel-size', not_found)
            if max_indel_size == not_found:
                error_list.append('*** ERROR: the key "max-indel-size" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_int(max_indel_size, minimum=0):
                error_list.append('*** ERROR: the key "max-indel-size" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "eXpress parameters" - key "no-bias-correct"
            no_bias_correct = express_option_dict.get('eXpress parameters', {}).get('no-bias-correct', not_found)
            if no_bias_correct == not_found:
                error_list.append('*** ERROR: the key "no-bias-correct" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_code(no_bias_correct, get_no_bias_correct_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "no-bias-correct" has to be {get_no_bias_correct_code_list_text()}.')
                OK = False

            # check section "eXpress parameters" - key "no-error-model"
            no_error_model = express_option_dict.get('eXpress parameters', {}).get('no-error-model', not_found)
            if no_error_model == not_found:
                error_list.append('*** ERROR: the key "no-error-model" is not found in the section "eXpress parameters".')
                OK = False
            elif not xlib.check_code(no_error_model, get_no_error_model_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "no-error-model" has to be {get_no_error_model_code_list_text()}.')
                OK = False

            # check section "eXpress parameters" - key "other_parameters"
            not_allowed_parameters_list = ['no-update-check', 'frag-len-mean', 'frag-len-stddev', 'max-indel-size', 'fr-stranded', 'rf-stranded', 'f-stranded', 'r-stranded', 'no-bias-correct', 'no-error-model', 'output-dir']
            other_parameters = express_option_dict.get('eXpress parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "eXpress parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_express_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_express_process_script(cluster_name, current_run_dir):
    '''
    Build the current eXpress process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the eXpress option dictionary
    express_option_dict = xlib.get_option_dict(get_express_config_file())

    # get the options
    experiment_id = express_option_dict['identification']['experiment_id']
    assembly_software = express_option_dict['identification']['assembly_software']
    assembly_dataset_id = express_option_dict['identification']['assembly_dataset_id']
    assembly_type = express_option_dict['identification']['assembly_type']
    alignment_software = express_option_dict['identification']['alignment_software']
    alignment_dataset_id = express_option_dict['identification']['alignment_dataset_id']
    frag_len_mean = express_option_dict['eXpress parameters']['frag-len-mean']
    frag_len_stddev = express_option_dict['eXpress parameters']['frag-len-stddev']
    library_type = express_option_dict['eXpress parameters']['library_type']
    max_indel_size = express_option_dict['eXpress parameters']['max-indel-size']
    no_bias_correct = express_option_dict['eXpress parameters']['no-bias-correct']
    no_error_model = express_option_dict['eXpress parameters']['no-error-model']
    other_parameters = express_option_dict['eXpress parameters']['other_parameters']

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

    # set the alignment file path
    if alignment_software == xlib.get_bowtie2_code():
        alignment_file = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id)}/alignment.sam'

    # write the eXpress process script
    try:
        if not os.path.exists(os.path.dirname(get_express_process_script())):
            os.makedirs(os.path.dirname(get_express_process_script()))
        with open(get_express_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function run_express_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_express_anaconda_code()}\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '        express \\\n')
            script_file_id.write( '            --no-update-check \\\n')
            script_file_id.write(f'            --frag-len-mean {frag_len_mean} \\\n')
            script_file_id.write(f'            --frag-len-stddev {frag_len_stddev} \\\n')
            if library_type.lower() == 'fr-stranded':
                script_file_id.write( '            --fr-stranded \\\n')
            elif library_type.lower() == 'rf-stranded':
                script_file_id.write( '            --rf-stranded \\\n')
            elif library_type.lower() == 'f-stranded':
                script_file_id.write( '            --f-stranded \\\n')
            elif library_type.lower() == 'r-stranded':
                script_file_id.write( '            --r-stranded \\\n')
            script_file_id.write(f'            --max-indel-size {max_indel_size} \\\n')
            if no_bias_correct.upper() == 'YES':
                script_file_id.write( '            --no-bias-correct \\\n')
            if no_error_model.upper() == 'YES':
                script_file_id.write( '            --no-error-model \\\n')
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
            script_file_id.write(f'            --output-dir {current_run_dir} \\\n')
            script_file_id.write(f'            {transcriptome_file} \\\n')
            script_file_id.write(f'            {alignment_file}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error express $RC; fi\n')
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
            process_name = f'{xlib.get_express_name()} process'
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
            script_file_id.write( 'run_express_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_express_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_express_process_starter(current_run_dir):
    '''
    Build the starter of the current eXpress process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the eXpress process starter
    try:
        if not os.path.exists(os.path.dirname(get_express_process_starter())):
            os.makedirs(os.path.dirname(get_express_process_starter()))
        with open(get_express_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_express_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_express_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_express_config_file():
    '''
    Get the eXpress config file path.
    '''

    # assign the eXpress config file path
    express_config_file = f'{xlib.get_config_dir()}/{xlib.get_express_code()}-config.txt'

    # return the eXpress config file path
    return express_config_file

#-------------------------------------------------------------------------------

def get_express_process_script():
    '''
    Get the eXpress process script path in the local computer.
    '''

    # assign the eXpress script path
    express_process_script = f'{xlib.get_temp_dir()}/{xlib.get_express_code()}-process.sh'

    # return the eXpress script path
    return express_process_script

#-------------------------------------------------------------------------------

def get_express_process_starter():
    '''
    Get the eXpress process starter path in the local computer.
    '''

    # assign the eXpress process starter path
    express_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_express_code()}-process-starter.sh'

    # return the eXpress starter path
    return express_process_starter

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
    
def get_alignment_software_code_list():
    '''
    Get the code list of "alignment_software".
    '''

    return [xlib.get_bowtie2_code()]

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list_text():
    '''
    Get the code list of "alignment_software" as text.
    '''

    return f'{xlib.get_bowtie2_code()} ({xlib.get_bowtie2_name()}))'

#-------------------------------------------------------------------------------
    
def get_library_type_code_list():
    '''
    Get the code list of "library_type".
    '''

    return ['FR-STRANDED', 'RF-STRANDED', 'F-STRANDED', 'R-STRANDED', 'NONE']

#-------------------------------------------------------------------------------
    
def get_library_type_code_list_text():
    '''
    Get the code list of "library_type" as text.
    '''

    return str(get_library_type_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_no_bias_correct_code_list():
    '''
    Get the code list of "no-bias-correct".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_no_bias_correct_code_list_text():
    '''
    Get the code list of "no-bias-correct" as text.
    '''

    return str(get_no_bias_correct_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_no_error_model_code_list():
    '''
    Get the code list of "no-error-model".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_no_error_model_code_list_text():
    '''
    Get the code list of "no-error-model" as text.
    '''

    return str(get_no_error_model_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the eXpress process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
