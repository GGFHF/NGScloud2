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
This file contains functions related to the HTSeq process used in both console
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

def create_htseq_count_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', annotation_file='Arabidopsis_thaliana.TAIR10.36.gtf', alignment_dataset_id_list=['star-170101-235959','tophat-170101-235959']):
    '''
    Create htseq-count config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the htseq-count config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_htseq_count_config_file())):
            os.makedirs(os.path.dirname(get_htseq_count_config_file()))
        with open(get_htseq_count_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference and annotation files have to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write( '# The experiment_id, reference_dataset_id, and annotation_file names are fixed in the identification section.\n')
            file_id.write(f'# The alignment files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/alignment_dataset_id\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of htseq_count and their meaning in "https://github.com/simon-anders/htseq".\n')
            file_id.write( '#\n')
            file_id.write( '# In section "htseq-count parameters", the key "other_parameters" allows you to input additional parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of TopHat and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    other_parameters = --max-reads-in-buffer=30000000\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'reference_dataset_id = {reference_dataset_id}', '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'annotation_file = {annotation_file}', '# reference annotation file name'))
            for i in range(len(alignment_dataset_id_list)):
                # set the alignment software
                alignment_dataset_id = alignment_dataset_id_list[i]
                if alignment_dataset_id.startswith(xlib.get_bowtie2_code()):
                    alignment_software = xlib.get_bowtie2_code()
                elif alignment_dataset_id.startswith(xlib.get_gsnap_code()):
                    alignment_software = xlib.get_gsnap_code()
                elif alignment_dataset_id.startswith(xlib.get_hisat2_code()):
                    alignment_software = xlib.get_hisat2_code()
                elif alignment_dataset_id.startswith(xlib.get_star_code()):
                    alignment_software = xlib.get_star_code()
                elif alignment_dataset_id.startswith(xlib.get_tophat_code()):
                    alignment_software = xlib.get_tophat_code()
                # write the alignment dataset section
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '# This section has the information of the first alignment dataset.\n')
                file_id.write(f'[alignment-dataset-{i + 1}]\n')
                file_id.write( '{0:<50} {1}\n'.format(f'alignment_software = {alignment_software}', f'# alignment software: {get_alignment_software_code_list_text()}'))
                file_id.write( '{0:<50} {1}\n'.format(f'alignment_dataset_id = {alignment_dataset_id}', '# alignment dataset identification'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '# If there are more alignment datasets, you have to repeat the section alignment-dataset-1 with the data of each dataset.\n')
                    file_id.write( '# The section identification has to be alignment-dataset-n (n is an integer not repeated)\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the htseq-count parameters\n')
            file_id.write( '[htseq-count parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format( 'nprocesses = 4', '# number of CPU for use'))
            file_id.write( '{0:<50} {1}\n'.format( 'stranded = YES', f'# whether the data is from a strand-specific assay: {get_stranded_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'minaqual = 10', '# skip all reads with MAPQ (5th column in BAM file) alignment quality lower than the given minimum value'))
            file_id.write( '{0:<50} {1}\n'.format( 'type = exon', '# feature type (3rd column in GTF file) to be used, all features of other type are ignored'))
            file_id.write( '{0:<50} {1}\n'.format( 'idattr = gene_id', '# GTF file feature id used to identity the counts in the output table'))
            file_id.write( '{0:<50} {1}\n'.format( 'mode = UNION', f'# mode to handle reads overlapping more than one feature: {get_mode_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'nonunique = NONE', f'# Mode to handle reads that align to or are assigned to more than one feature in the overlap mode of choice: {get_nonunique_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format( 'other_parameters = NONE', '# additional parameters to the previous ones or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_htseq_count_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_htseq_count_process(cluster_name, log, function=None):
    '''
    Run a htseq-count process.
    '''

    # initialize the control variable
    OK = True

    # get the htseq-count option dictionary
    htseq_count_option_dict = xlib.get_option_dict(get_htseq_count_config_file())

    # get the experiment identification
    experiment_id = htseq_count_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the htseq-count config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_htseq_count_name()} config file ...\n')
    (OK, error_list) = check_htseq_count_config_file(strict=True)
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

    # check HTSeq is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_htseq_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_htseq_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_htseq_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_htseq_count_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the htseq-count process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_htseq_count_process_script()} ...\n')
        (OK, error_list) = build_htseq_count_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the htseq-count process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_htseq_count_process_script()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_htseq_count_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_htseq_count_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the htseq-count process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_htseq_count_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_htseq_count_process_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the htseq-count process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_htseq_count_process_starter()} ...\n')
        (OK, error_list) = build_htseq_count_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the htseq-count process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_htseq_count_process_starter()} to the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_htseq_count_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_htseq_count_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the htseq-count process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_htseq_count_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_htseq_count_process_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the htseq-count process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_htseq_count_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_htseq_count_process_starter()), log)

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

def check_htseq_count_config_file(strict):
    '''
    Check the htseq-count config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        htseq_count_option_dict = xlib.get_option_dict(get_htseq_count_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in htseq_count_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = htseq_count_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = htseq_count_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "annotation_file"
            annotation_file = htseq_count_option_dict.get('identification', {}).get('annotation_file', not_found)
            if annotation_file == not_found:
                error_list.append('*** ERROR: the key "annotation_file" is not found in the section "identification".')
                OK = False
            elif os.path.splitext(annotation_file)[1] not in ['.gtf', '.gff']:
                error_list.append('*** ERROR: the key "annotation_file" has to be a file name with .gtf/.gff extension.')
                OK = False

        # check section "alignment-dataset-1"
        if 'alignment-dataset-1' not in sections_list:
            error_list.append('*** ERROR: the section "alignment-dataset-1" is not found.')
            OK = False

        # check all sections "alignment-dataset-n"
        for section in sections_list:

            if section not in ['identification', 'htseq-count parameters']:

                # check than the section identification is like alignment-dataset-n 
                if not re.match('^alignment-dataset-[0-9]+$', section):
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "alignment-dataset-n" - key "alignment_software"
                    alignment_software = htseq_count_option_dict.get(section, {}).get('alignment_software', not_found)
                    if alignment_software == not_found:
                        error_list.append(f'*** ERROR: the key "alignment_software" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_code(alignment_software, get_alignment_software_code_list(), case_sensitive=False):
                        error_list.append(f'*** ERROR: the key "alignment_software" has to be {get_alignment_software_code_list_text()}.')
                        OK = False

                    # check section "alignment-dataset-n" - key "alignment_dataset_id"
                    alignment_dataset_id = htseq_count_option_dict.get(section, {}).get('alignment_dataset_id', not_found)
                    if alignment_dataset_id == not_found:
                        error_list.append(f'*** ERROR: the key "alignment_dataset_id" is not found in the section "{section}".')
                        OK = False
                    elif not xlib.check_startswith(alignment_dataset_id, get_alignment_software_code_list(), case_sensitive=True):
                        error_list.append(f'*** ERROR: the key "alignment_dataset_id" has to start with {get_alignment_software_code_list_text()}.')
                        OK = False

        # check section "htseq-count parameters"
        if 'htseq-count parameters' not in sections_list:
            error_list.append('*** ERROR: the section "htseq-count parameters" is not found.')
            OK = False
        else:

            # check section "htseq-count parameters" - key "nprocesses"
            nprocesses = htseq_count_option_dict.get('htseq-count parameters', {}).get('nprocesses', not_found)
            if nprocesses == not_found:
                error_list.append('*** ERROR: the key "nprocesses" is not found in the section "htseq-count parameters".')
                OK = False
            elif not xlib.check_int(nprocesses, minimum=1):
                error_list.append('*** ERROR: the key "nprocesses" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "htseq-count parameters" - key "stranded"
            stranded = htseq_count_option_dict.get('htseq-count parameters', {}).get('stranded', not_found)
            if stranded == not_found:
                error_list.append('*** ERROR: the key "stranded" is not found in the section "htseq-count parameters".')
                OK = False
            elif not xlib.check_code(stranded, get_stranded_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "stranded" has to be {get_stranded_code_list_text()}.')
                OK = False

            # check section "htseq-count parameters" - key "minaqual"
            minaqual = htseq_count_option_dict.get('htseq-count parameters', {}).get('minaqual', not_found)
            if minaqual == not_found:
                error_list.append('*** ERROR: the key "minaqual" is not found in the section "htseq-count parameters".')
                OK = False
            elif not xlib.check_int(minaqual):
                error_list.append('*** ERROR: the key "minaqual" has to be an integer number.')
                OK = False

            # check section "htseq-count parameters" - key "type"
            type = htseq_count_option_dict.get('htseq-count parameters', {}).get('type', not_found)
            if type == not_found:
                error_list.append('*** ERROR: the key "type" is not found in the section "htseq-count parameters".')
                OK = False

            # check section "htseq-count parameters" - key "idattr"
            idattr = htseq_count_option_dict.get('htseq-count parameters', {}).get('idattr', not_found)
            if idattr == not_found:
                error_list.append('*** ERROR: the key "idattr" is not found in the section "htseq-count parameters".')
                OK = False

            # check section "htseq-count parameters" - key "mode"
            mode = htseq_count_option_dict.get('htseq-count parameters', {}).get('mode', not_found)
            if mode == not_found:
                error_list.append('*** ERROR: the key "mode" is not found in the section "htseq-count parameters".')
                OK = False
            elif not xlib.check_code(mode, get_mode_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "mode" has to be {get_mode_code_list_text()}.')
                OK = False

            # check section "htseq-count parameters" - key "nonunique"
            nonunique = htseq_count_option_dict.get('htseq-count parameters', {}).get('nonunique', not_found)
            if nonunique == not_found:
                error_list.append('*** ERROR: the key "nonunique" is not found in the section "htseq-count parameters".')
                OK = False
            elif not xlib.check_code(nonunique, get_nonunique_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "nonunique" has to be {get_nonunique_code_list_text()}.')
                OK = False

            # check section "htseq-count parameters" - key "other_parameters"
            not_allowed_parameters_list = ['nprocesses', 'format', 'stranded', 'minaqual', 'type', 'idattr', 'mode', 'nonunique', 'quiet']
            other_parameters = htseq_count_option_dict.get('htseq-count parameters', {}).get('other_parameters', not_found)
            if other_parameters == not_found:
                error_list.append('*** ERROR: the key "other_parameters" is not found in the section "htseq-count parameters".')
                OK = False
            elif other_parameters.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(other_parameters, "other_parameters", not_allowed_parameters_list)
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_htseq_count_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_htseq_count_process_script(cluster_name, current_run_dir):
    '''
    Build the current htseq-count process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the htseq-count option dictionary
    htseq_count_option_dict = xlib.get_option_dict(get_htseq_count_config_file())

    # get the options
    experiment_id = htseq_count_option_dict['identification']['experiment_id']
    reference_dataset_id = htseq_count_option_dict['identification']['reference_dataset_id']
    annotation_file = htseq_count_option_dict['identification']['annotation_file']
    nprocesses = htseq_count_option_dict['htseq-count parameters']['nprocesses']
    stranded = htseq_count_option_dict['htseq-count parameters']['stranded']
    minaqual = htseq_count_option_dict['htseq-count parameters']['minaqual']
    type = htseq_count_option_dict['htseq-count parameters']['type']
    idattr = htseq_count_option_dict['htseq-count parameters']['idattr']
    mode = htseq_count_option_dict['htseq-count parameters']['mode']
    nonunique = htseq_count_option_dict['htseq-count parameters']['nonunique']
    other_parameters = htseq_count_option_dict['htseq-count parameters']['other_parameters']

    # get the sections list
    sections_list = []
    for section in htseq_count_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build alignment dataset identification list
    alignment_software_list = []
    alignment_dataset_id_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^alignment-dataset-[0-9]+$', section):
            alignment_software_list.append(htseq_count_option_dict[section]['alignment_software'])
            alignment_dataset_id_list.append(htseq_count_option_dict[section]['alignment_dataset_id'])

    # set the annotation file path
    annotation_file = xlib.get_cluster_reference_file(reference_dataset_id, annotation_file)

    # write the htseq-count process script
    try:
        if not os.path.exists(os.path.dirname(get_htseq_count_process_script())):
            os.makedirs(os.path.dirname(get_htseq_count_process_script()))
        with open(get_htseq_count_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function print_htseq_count_version\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_htseq_anaconda_code()}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    # -- htseq-count --version\n')
            script_file_id.write( '    conda deactivate\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_htseq_count_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_htseq_anaconda_code()}\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Counting reads ..."\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '        htseq-count \\\n')
            # -- script_file_id.write(f'            --nprocesses={nprocesses} \\\n')
            script_file_id.write( '            --format=bam \\\n')
            script_file_id.write(f'            --stranded={stranded.lower()} \\\n')
            script_file_id.write(f'            --minaqual={minaqual} \\\n')
            script_file_id.write(f'            --type={type} \\\n')
            script_file_id.write(f'            --idattr={idattr} \\\n')
            script_file_id.write(f'            --mode={mode.lower()} \\\n')
            script_file_id.write(f'            --nonunique={nonunique.lower()} \\\n')
            script_file_id.write( '            --quiet \\\n')
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
            for i in range(len(alignment_dataset_id_list)):
                alignment_files = f'{xlib.get_cluster_experiment_result_dataset_dir(experiment_id, alignment_dataset_id_list[i])}/*.sorted.bam'
                script_file_id.write(f'            {alignment_files} \\\n')
            script_file_id.write(f'            {annotation_file} \\\n')
            script_file_id.write(f'            > read-count.txt\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error htseq-count $RC; fi\n')
            script_file_id.write( '    echo "Reads are counted."\n')
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
            process_name = f'{xlib.get_htseq_count_name()} process'
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
            script_file_id.write( 'print_htseq_count_version\n')
            script_file_id.write( 'run_htseq_count_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_htseq_count_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_htseq_count_process_starter(current_run_dir):
    '''
    Build the starter of the current htseq-count process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the htseq-count process starter
    try:
        if not os.path.exists(os.path.dirname(get_htseq_count_process_starter())):
            os.makedirs(os.path.dirname(get_htseq_count_process_starter()))
        with open(get_htseq_count_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_htseq_count_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_htseq_count_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_htseq_count_config_file():
    '''
    Get the htseq-count config file path.
    '''

    # assign the htseq-count config file path
    htseq_count_config_file = f'{xlib.get_config_dir()}/{xlib.get_htseq_count_code()}-config.txt'

    # return the htseq-count config file path
    return htseq_count_config_file

#-------------------------------------------------------------------------------

def get_htseq_count_process_script():
    '''
    Get the htseq-count process script path in the local computer.
    '''

    # assign the htseq-count script path
    htseq_count_process_script = f'{xlib.get_temp_dir()}/{xlib.get_htseq_count_code()}-process.sh'

    # return the htseq-count script path
    return htseq_count_process_script

#-------------------------------------------------------------------------------

def get_htseq_count_process_starter():
    '''
    Get the htseq-count process starter path in the local computer.
    '''

    # assign the htseq-count process starter path
    htseq_count_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_htseq_count_code()}-process-starter.sh'

    # return the htseq-count starter path
    return htseq_count_process_starter

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list():
    '''
    Get the code list of "alignment_software".
    '''

    return [xlib.get_bowtie2_code(), xlib.get_gsnap_code(), xlib.get_hisat2_code(), xlib.get_star_code(), xlib.get_tophat_code()]

#-------------------------------------------------------------------------------
    
def get_alignment_software_code_list_text():
    '''
    Get the code list of "alignment_software" as text.
    '''

    return f'{xlib.get_bowtie2_code()} ({xlib.get_bowtie2_name()}) or {xlib.get_gsnap_code()} ({xlib.get_gsnap_name()}) or {xlib.get_hisat2_code()} ({xlib.get_hisat2_name()}) or {xlib.get_star_code()} ({xlib.get_star_name()}) or {xlib.get_tophat_code()} ({xlib.get_tophat_name()})'

#-------------------------------------------------------------------------------
    
def get_stranded_code_list():
    '''
    Get the code list of "stranded".
    '''

    return ['YES', 'NO', 'REVERSE']

#-------------------------------------------------------------------------------
    
def get_stranded_code_list_text():
    '''
    Get the code list of "stranded".
    '''

    return str(get_stranded_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_mode_code_list():
    '''
    Get the code list of "mode".
    '''

    return ['UNION', 'INTERSECTION-STRICT', 'INTERSECTION-NONEMPTY']

#-------------------------------------------------------------------------------
    
def get_mode_code_list_text():
    '''
    Get the code list of "mode".
    '''

    return str(get_mode_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_nonunique_code_list():
    '''
    Get the code list of "nonunique".
    '''

    return ['NONE', 'ALL', 'FRACTION', 'RANDOM']

#-------------------------------------------------------------------------------
    
def get_nonunique_code_list_text():
    '''
    Get the code list of "nonunique".
    '''

    return str(get_nonunique_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the HTSeq process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
