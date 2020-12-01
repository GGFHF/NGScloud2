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
This file contains functions related to the BUSCO process used in both console
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

def create_busco_config_file(experiment_id='exp001', assembly_dataset_id='sdnt-170101-235959', assembly_type='CONTIGS'):
    '''
    Create BUSCO config file with the default options. It is necessary
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

    # create the BUSCO config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_busco_config_file())):
            os.makedirs(os.path.dirname(get_busco_config_file()))
        with open(get_busco_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write(f'# The reference file has to be located in the cluster directory {xlib.get_cluster_reference_dir()}/experiment_id/reference_dataset_id\n')
            file_id.write(f'# The assembly files have to be located in the cluster directory {xlib.get_cluster_result_dir()}/experiment_id/assembly_dataset_id\n')
            file_id.write( '# The experiment_id and assembly_dataset_id names are fixed in the identification section.\n')
            file_id.write( '#\n')
            file_id.write( '# In section "BUSCO parameters", the key "augustus_options" allows you to input additional August parameters in the format:\n')
            file_id.write( '#\n')
            file_id.write( '#    augustus_options = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]\n')
            file_id.write( '#\n')
            file_id.write( '# parameter-i is a parameter name of Augustus and value-i a valid value of parameter-i, e.g.\n')
            file_id.write( '#\n')
            file_id.write( '#    augustus_options = --translation_table=6 --progress=true\n')
            file_id.write( '#\n')
            file_id.write( '# You can consult the parameters of BUSCO and their meaning in "http://busco.ezlab.org/"\n')
            file_id.write( '# and the ones of August in "http://bioinf.uni-greifswald.de/augustus/".\n')
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_software = {assembly_software}', f'# assembly software: {get_assembly_software_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_dataset_id = {assembly_dataset_id}', '# assembly dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format(f'assembly_type = {assembly_type}', f'# assembly type: CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()}; NONE in any other case'))
            file_id.write( '\n')
            file_id.write( '# This section has the information to set the BUSCO parameters\n')
            file_id.write( '[BUSCO parameters]\n')
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of threads/cores for use'))
            file_id.write( '{0:<50} {1}\n'.format('lineage_data_url = https://busco-data.ezlab.org/v4/data/lineages/viridiplantae_odb10.2020-09-10.tar.gz', '# the url of lineage data file that will be used'))
            file_id.write( '{0:<50} {1}\n'.format('mode = TRAN', f'# mode: {get_mode_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('evalue = 1E-03', '# E-value cutoff for BLAST searches'))
            file_id.write( '{0:<50} {1}\n'.format('limit = 3', '# number of candidate regions to consider'))
            file_id.write( '{0:<50} {1}\n'.format('species = NONE', '# identifier of existing Augustus species gene finding parameters or NONE'))
            file_id.write( '{0:<50} {1}\n'.format('long = NO', f'# Augustus optimization mode for self-training: {get_long_code_list_text()}'))
            file_id.write( '{0:<50} {1}\n'.format('augustus_options = NONE', '# additional parameters to August or NONE'))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_busco_config_file()} can not be recreated')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_busco_process(cluster_name, log, function=None):
    '''
    Run a BUSCO process.
    '''

    # initialize the control variable
    OK = True

    # get the BUSCO option dictionary
    busco_option_dict = xlib.get_option_dict(get_busco_config_file())

    # get the experiment identification
    experiment_id = busco_option_dict['identification']['experiment_id']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the BUSCO config file
    log.write(f'{xlib.get_separator()}\n')
    log.write(f'Checking the {xlib.get_busco_name()} config file ...\n')
    (OK, error_list) = check_busco_config_file(strict=True)
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

    # check BUSCO is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_anaconda_package(xlib.get_busco_anaconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write(f'*** ERROR: {xlib.get_busco_name()} is not installed.\n')
                OK = False
        else:
            log.write(f'*** ERROR: The verification of {xlib.get_busco_name()} installation could not be performed.\n')

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_busco_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the BUSCO process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process script {get_busco_process_script()} ...\n')
        (OK, error_list) = build_busco_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the BUSCO process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process script {get_busco_process_script()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_busco_process_script())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_busco_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the BUSCO process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_busco_process_script())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_busco_process_script())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the BUSCO process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the process starter {get_busco_process_starter()} ...\n')
        (OK, error_list) = build_busco_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the busco process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the process starter {get_busco_process_starter()} in the directory {current_run_dir} ...\n')
        cluster_path = f'{current_run_dir}/{os.path.basename(get_busco_process_starter())}'
        (OK, error_list) = xssh.put_file(sftp_client, get_busco_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the BUSCO process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Setting on the run permision of {current_run_dir}/{os.path.basename(get_busco_process_starter())} ...\n')
        command = f'chmod u+x {current_run_dir}/{os.path.basename(get_busco_process_starter())}'
        (OK, _, _) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the BUSCO process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the process script {current_run_dir}/{os.path.basename(get_busco_process_starter())} ...\n')
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_busco_process_starter()), log)

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

def check_busco_config_file(strict):
    '''
    Check the BUSCO config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        busco_option_dict = xlib.get_option_dict(get_busco_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in busco_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = busco_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "assembly_software"
            assembly_software = busco_option_dict.get('identification', {}).get('assembly_software', not_found)
            if assembly_software == not_found:
                error_list.append('*** ERROR: the key "assembly_software" is not found in the section "identification".')
                OK = False
            elif not xlib.check_code(assembly_software, get_assembly_software_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "assembly_software" has to be {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_dataset_id"
            assembly_dataset_id = busco_option_dict.get('identification', {}).get('assembly_dataset_id', not_found)
            if assembly_dataset_id == not_found:
                error_list.append('*** ERROR: the key "assembly_dataset_id" is not found in the section "identification".')
                OK = False
            elif not xlib.check_startswith(assembly_dataset_id, get_assembly_software_code_list(), case_sensitive=True):
                error_list.append(f'*** ERROR: the key "assembly_dataset_id" has to start with {get_assembly_software_code_list_text()}.')
                OK = False

            # check section "identification" - key "assembly_type"
            assembly_type = busco_option_dict.get('identification', {}).get('assembly_type', not_found)
            if assembly_type == not_found:
                error_list.append('*** ERROR: the key "assembly_type" is not found in the section "identification".')
                OK = False
            elif assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() not in ['CONTIGS', 'SCAFFOLDS'] or \
                not assembly_dataset_id.startswith(xlib.get_soapdenovotrans_code()) and assembly_type.upper() != 'NONE':
                    error_list.append(f'*** ERROR: the key "assembly_type" has to be CONTIGS or SCAFFOLDS in {xlib.get_soapdenovotrans_name()} or NONE in any other case.')
                    OK = False

        # check section "BUSCO parameters"
        if 'BUSCO parameters' not in sections_list:
            error_list.append('*** ERROR: the section "BUSCO parameters" is not found.')
            OK = False
        else:

            # check section "BUSCO parameters" - key "ncpu"
            ncpu = busco_option_dict.get('BUSCO parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "BUSCO parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BUSCO parameters" - key "lineage_data_url"
            lineage_data_url = busco_option_dict.get('BUSCO parameters', {}).get('lineage_data_url', not_found)
            if lineage_data_url == not_found:
                error_list.append('*** ERROR: the key "lineage_data_url" is not found in the section "BUSCO parameters"')
                OK = False
            else:
                try:
                    urllib.request.urlopen(lineage_data_url)
                except Exception as e:
                    error_list.append(f'*** EXCEPTION: "{e}".')
                    error_list.append('*** ERROR: the key "lineage_data_url" has to be a reachable address.')
                    OK = False

            # check section "BUSCO parameters" - key "mode"
            mode = busco_option_dict.get('BUSCO parameters', {}).get('mode', not_found)
            if mode == not_found:
                error_list.append('*** ERROR: the key "mode" is not found in the section "BUSCO parameters".')
                OK = False
            elif not xlib.check_code(mode, get_mode_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "mode" has to be {get_mode_code_list_text()}.')
                OK = False

            # check section "BUSCO parameters" - key "evalue"
            evalue = busco_option_dict.get('BUSCO parameters', {}).get('evalue', not_found)
            if evalue == not_found:
                error_list.append('*** ERROR: the key "evalue" is not found in the section "BUSCO parameters".')
                OK = False
            elif not xlib.check_float(evalue, minimum=0., mne=1E-12):
                error_list.append('*** ERROR: the key "evalue" has to be a float number greater than 0.')
                OK = False

            # check section "BUSCO parameters" - key "limit"
            limit = busco_option_dict.get('BUSCO parameters', {}).get('limit', not_found)
            if limit == not_found:
                error_list.append('*** ERROR: the key "limit" is not found in the section "BUSCO parameters".')
                OK = False
            elif not xlib.check_int(limit, minimum=1):
                error_list.append('*** ERROR: the key "limit" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "BUSCO parameters" - key "species"
            species = busco_option_dict.get('BUSCO parameters', {}).get('species', not_found)
            if species == not_found:
                error_list.append('*** ERROR: the key "species" is not found in the section "BUSCO parameters"')
                OK = False

            # check section "BUSCO parameters" - key "long"
            long = busco_option_dict.get('BUSCO parameters', {}).get('long', not_found)
            if long == not_found:
                error_list.append('*** ERROR: the key "long" is not found in the section "BUSCO parameters".')
                OK = False
            elif not xlib.check_code(long, get_long_code_list(), case_sensitive=False):
                error_list.append(f'*** ERROR: the key "long" has to be {get_long_code_list_text()}.')
                OK = False

            # check section "BUSCO parameters" - key "augustus_options"
            augustus_options = busco_option_dict.get('BUSCO parameters', {}).get('augustus_options', not_found)
            if augustus_options == not_found:
                error_list.append('*** ERROR: the key "augustus_options" is not found in the section "BUSCO parameters".')
                OK = False
            elif augustus_options.upper() != 'NONE':
                (OK, error_list2) = xlib.check_parameter_list(augustus_options, "augustus_options", [])
                error_list = error_list + error_list2

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append(f'\nThe {xlib.get_busco_name()} config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_busco_process_script(cluster_name, current_run_dir):
    '''
    Build the current BUSCO process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the BUSCO option dictionary
    busco_option_dict = xlib.get_option_dict(get_busco_config_file())

    # get the options
    experiment_id = busco_option_dict['identification']['experiment_id']
    assembly_software = busco_option_dict['identification']['assembly_software']
    assembly_dataset_id = busco_option_dict['identification']['assembly_dataset_id']
    assembly_type = busco_option_dict['identification']['assembly_type']
    ncpu = busco_option_dict['BUSCO parameters']['ncpu']
    lineage_data_url = busco_option_dict['BUSCO parameters']['lineage_data_url']
    mode = busco_option_dict['BUSCO parameters']['mode'].lower()
    evalue = busco_option_dict['BUSCO parameters']['evalue']
    limit = busco_option_dict['BUSCO parameters']['limit']
    species = busco_option_dict['BUSCO parameters']['species']
    long = busco_option_dict['BUSCO parameters']['long'].upper()
    augustus_options = busco_option_dict['BUSCO parameters']['augustus_options'].upper()

    # get the file and name from the lineage data url
    lineage_data_file = lineage_data_url.split("/")[-1]
    # -- lineage_data = lineage_data_file[:lineage_data_file.find('.tar.gz')]
    point_pos = lineage_data_file.find('.')
    lineage_data = lineage_data_file[:point_pos]

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

    # write the BUSCO process script
    try:
        if not os.path.exists(os.path.dirname(get_busco_process_script())):
            os.makedirs(os.path.dirname(get_busco_process_script()))
        with open(get_busco_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
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
            script_file_id.write( 'function download_lineage_data\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Downloading lineage data ..."\n')
            download_script = f'import requests; r = requests.get(\'{lineage_data_url}\') ; open(\'{lineage_data_file}\' , \'wb\').write(r.content)'
            script_file_id.write(f'    $MINICONDA3_BIN_PATH/python3 -c "{download_script}"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error download_script $RC; fi\n')
            script_file_id.write(f'    tar -xzvf ./{lineage_data_file}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
            script_file_id.write(f'    rm ./{lineage_data_file}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error rm $RC; fi\n')
            script_file_id.write( '    echo "Lineage data are downloaded."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function run_busco_process\n')
            script_file_id.write( '{\n')
            script_file_id.write(f'    source activate {xlib.get_busco_anaconda_code()}\n')
            script_file_id.write( '    cd $CURRENT_DIR\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Assessing the transcriptome quality ..."\n')
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format(separator=False)}" \\\n')
            script_file_id.write( '        busco \\\n')
            script_file_id.write(f'            --cpu={ncpu} \\\n')
            script_file_id.write(f'            --lineage_dataset=./{lineage_data} \\\n')
            script_file_id.write(f'            --mode={mode} \\\n')
            script_file_id.write(f'            --evalue={evalue} \\\n')
            script_file_id.write(f'            --limit={limit} \\\n')
            if species.upper() != 'NONE':
                script_file_id.write(f'            --species={species} \\\n')
            if long == 'YES':
                script_file_id.write( '            --long \\\n')
            if augustus_options.upper() != 'NONE':
                script_file_id.write(f'            --august_options="{augustus_options}" \\\n')
            script_file_id.write(f'            --in={transcriptome_file} \\\n')
            script_file_id.write(f'            --out={os.path.basename(current_run_dir)}\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error run_BUSCO.py $RC; fi\n')
            script_file_id.write( '    echo "The assessment is done."\n')
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
            process_name = f'{xlib.get_busco_name()} process'
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
            script_file_id.write( 'download_lineage_data\n')
            script_file_id.write( 'run_busco_process\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_busco_process_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_busco_process_starter(current_run_dir):
    '''
    Build the starter of the current BUSCO process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the BUSCO process starter
    try:
        if not os.path.exists(os.path.dirname(get_busco_process_starter())):
            os.makedirs(os.path.dirname(get_busco_process_starter()))
        with open(get_busco_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write(f'{current_run_dir}/{os.path.basename(get_busco_process_script())} &>>{current_run_dir}/{xlib.get_cluster_log_file()}\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append(f'*** ERROR: The file {get_busco_process_starter()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_busco_config_file():
    '''
    Get the BUSCO config file path.
    '''

    # assign the BUSCO config file path
    busco_config_file = f'{xlib.get_config_dir()}/{xlib.get_busco_code()}-config.txt'

    # return the BUSCO config file path
    return busco_config_file

#-------------------------------------------------------------------------------

def get_busco_process_script():
    '''
    Get the BUSCO process script path in the local computer.
    '''

    # assign the BUSCO script path
    busco_process_script = f'{xlib.get_temp_dir()}/{xlib.get_busco_code()}-process.sh'

    # return the BUSCO script path
    return busco_process_script

#-------------------------------------------------------------------------------

def get_busco_process_starter():
    '''
    Get the BUSCO process starter path in the local computer.
    '''

    # assign the BUSCO process starter path
    busco_process_starter = f'{xlib.get_temp_dir()}/{xlib.get_busco_code()}-process-starter.sh'

    # return the BUSCO starter path
    return busco_process_starter

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
    
def get_mode_code_list():
    '''
    Get the code list of "mode".
    '''

    return ['GENO', 'TRAN', 'PROT']

#-------------------------------------------------------------------------------
    
def get_mode_code_list_text():
    '''
    Get the code list of "mode".
    '''

    return 'GENO (genome assemblies, DNA) or TRAN (transcriptome assemblies, DNA) or PROT (annotated gene sets, proteins)'

#-------------------------------------------------------------------------------
    
def get_long_code_list():
    '''
    Get the code list of "long".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_long_code_list_text():
    '''
    Get the code list of "long" as text.
    '''

    return str(get_long_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the BUSCO process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
