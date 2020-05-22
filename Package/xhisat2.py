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

def create_hisat2_config_file(experiment_id='exp001', reference_dataset_id='Athaliana', reference_file='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa', splice_site_file='NONE', exon_file='NONE', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create HISAT2 config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the HISAT2 config file and write the default options
    if True:
    #try:
        if not os.path.exists(os.path.dirname(get_hisat2_config_file())):
            os.makedirs(os.path.dirname(get_hisat2_config_file()))
        with open(get_hisat2_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The reference file has to be located in the cluster directory {0}/experiment_id/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The assembly files have to be located in the cluster directory {0}/experiment_id/assembly_dataset_id'.format(xlib.get_cluster_result_dir())))
            file_id.write( '{0}\n'.format('# The read files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file, assembly_dataset_id and read_dataset_id are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of HISAT2 and their meaning in https://ccb.jhu.edu/software/hisat2/.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# In section "HISAT2 parameters", the key "other_parameters" allows you to input additional parameters in the format:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --parameter-1[=value-1][; --parameter-2[=value-2][; ...; --parameter-n[=value-n]]]'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# parameter-i is a parameter name of HISAT2 and value-i a valid value of parameter-i, e.g.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    other_parameters = --new-summary; --score-min=L,0,-0.2'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name'))
            file_id.write( '{0:<50} {1}\n'.format('splice_site_file = {0}'.format(splice_site_file), '# splice site file name'))
            file_id.write( '{0:<50} {1}\n'.format('exon_file = {0}'.format(exon_file), '# exon file name'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the HISAT2 parameters'))
            file_id.write( '{0}\n'.format('[HISAT2 parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('index_building = YES', '# index building when a reference is used: {0}'.format(get_index_building_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('large_index = NO', '# a large index is force, even if the reference is less than ~ 4 billion nucleotides long: {0}'.format(get_large_index_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('threads = 4', '# number of threads for use'))
            file_id.write( '{0:<50} {1}\n'.format('min_mp = 2', '# minimum mismatch penalty'))
            file_id.write( '{0:<50} {1}\n'.format('max_mp = 6', '# maximum mismatch penalty'))
            file_id.write( '{0:<50} {1}\n'.format('no_softclip = NO', '# disallow soft-clipping: {0}'.format(get_no_softclip_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('min_sp = 1', '# minimum penalty for soft-clipping per base'))
            file_id.write( '{0:<50} {1}\n'.format('max_sp = 2', '# maximum penalty for soft-clipping per base'))
            file_id.write( '{0:<50} {1}\n'.format('np = 1', '# penalty for positions where the read, reference, or both, contain an ambiguous character such as N'))
            file_id.write( '{0:<50} {1}\n'.format('open_rdg = 5', '# read gap open penalty'))
            file_id.write( '{0:<50} {1}\n'.format('extend_rdg = 3', '# read gap extend penalty'))
            file_id.write( '{0:<50} {1}\n'.format('open_rfg = 5', '# reference gap open penalty'))
            file_id.write( '{0:<50} {1}\n'.format('extend_rfg = 3', '# reference gap extend penalty'))
            file_id.write( '{0:<50} {1}\n'.format('orientation = FR', '# orientation of paired-end reads: {0}'.format(get_orientation_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('quality-score = 33', '# FASTQ quality score: {}'.format(get_quality_score_code_list_text())))
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
    #except Exception as e:
    #    error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_hisat2_config_file()))
    #    OK = False

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
    log.write('Checking the {0} config file ...\n'.format(xlib.get_hisat2_name()))
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
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the HISAT2 is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_hisat2_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_hisat2_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_hisat2_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_hisat2_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory path is {0}.\n'.format(current_run_dir))
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the HISAT2 process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_hisat2_process_script()))
        (OK, error_list) = build_hisat2_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the HISAT2 process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} in the directory {1} of the master ...\n'.format(get_hisat2_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_hisat2_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_hisat2_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the HISAT2 process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_hisat2_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_hisat2_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the HISAT2 process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_hisat2_process_starter()))
        (OK, error_list) = build_hisat2_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the HISAT2 process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} in the directory {1} of the master ...\n'.format(get_hisat2_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_hisat2_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_hisat2_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the HISAT2 process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_hisat2_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_hisat2_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the HISAT2 process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_hisat2_process_starter())))
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
        error_list.append('*** ERROR: The syntax is WRONG.')
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
                error_list.append('*** ERROR: the key "index_building" has to be {0}.'.format(get_index_building_code_list_text()))
                OK = False

            # check section "HISAT2 parameters" - key "large_index"
            large_index = hisat2_option_dict.get('HISAT2 parameters', {}).get('large_index', not_found)
            if large_index == not_found:
                error_list.append('*** ERROR: the key "large_index" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(large_index, get_large_index_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "large_index" has to be {0}.'.format(get_large_index_code_list_text()))
                OK = False

            # check section "HISAT2 parameters" - key "threads"
            threads = hisat2_option_dict.get('HISAT2 parameters', {}).get('threads', not_found)
            if threads == not_found:
                error_list.append('*** ERROR: the key "threads" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_int(threads, minimum=1):
                error_list.append('*** ERROR: the key "threads" has to be an integer number greater than or equal to 1.')
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
                error_list.append('*** ERROR: The value max_mp value ({0}) is less than the min_mp value ({1}).'.format(max_mp, min_mp))
                OK = False

            # check section "HISAT2 parameters" - key "no_softclip"
            no_softclip = hisat2_option_dict.get('HISAT2 parameters', {}).get('no_softclip', not_found)
            if no_softclip == not_found:
                error_list.append('*** ERROR: the key "no_softclip" is not found in the section "HISAT2 parameters".')
                OK = False
            elif not xlib.check_code(no_softclip, get_no_softclip_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "no_softclip" has to be {0}.'.format(get_no_softclip_code_list_text()))
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
                error_list.append('*** ERROR: The value max_sp value ({0}) is less than the min_sp value ({1}).'.format(max_sp, min_sp))
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
                error_list.append('*** ERROR: the key "orientation" has to be {0}.'.format(get_orientation_code_list_text()))
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
            not_allowed_parameters_list = ['nthreads', 'qseq', 'phred33', 'phred64', 'mp', 'sp', 'no-softclip', 'np', 'rdg', 'rfg', 'known-splicesite-infile', 'time', 'un', 'un-gz', 'un-bz2', 'al', 'al-gz', 'al-bz2', 'un-conc', 'un-conc-gz', 'un-conc-bz2', 'al-conc', 'al-conc-gz', 'al-conc-bz2', 'quiet', 'summary-file']
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
                error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                OK = False

            # check section "library" - key "read_type"
            read_type = hisat2_option_dict.get('library', {}).get('read_type', not_found)
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

            if section not in ['identification', 'HISAT2 parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = hisat2_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = hisat2_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_hisat2_name()))

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
    splice_site_file = hisat2_option_dict['identification']['splice_site_file']
    exon_file = hisat2_option_dict['identification']['exon_file']
    read_dataset_id = hisat2_option_dict['identification']['read_dataset_id']
    index_building = hisat2_option_dict['HISAT2 parameters']['index_building']
    large_index = hisat2_option_dict['HISAT2 parameters']['large_index']
    threads = hisat2_option_dict['HISAT2 parameters']['threads']
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

    # set the cluster reference dataset directory
    cluster_reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)

    # set the cluster reference file
    cluster_reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)

    # set the cluster splice site file
    if splice_site_file.upper() != 'NONE':
        cluster_splice_site_file = xlib.get_cluster_reference_file(reference_dataset_id, splice_site_file)

    # set the cluster exon file
    if exon_file.upper() != 'NONE':
        cluster_exon_file = xlib.get_cluster_reference_file(reference_dataset_id, exon_file)

    # set the directory and basename of the index HISAT2
    reference_file_name, reference_file_extension = os.path.splitext(reference_file)
    hisat2_index_dir = '{0}/{1}-hisat2_indexes'.format(cluster_reference_dataset_dir, reference_file_name)
    hisat2_index_basename = 'hisat2_indexes'

    # set the alignment file
    alignment_file = 'alignment.sam'

    # set the file of unpaired reads that fail to align
    un_gz = 'unpairednotaligned.fastq.gz'

    # set the file of unpaired reads that align at least once
    al_gz = 'unpairedaligned.fastq.gz'

    # set the file of paired-end reads that fail to align concordantly
    un_conc_gz = 'pairednotaligned.fastq.gz'

    # set the file of paired-end reads that align concordantly at least once
    al_conc_gz = 'pairednotaligned.fastq.gz'

    # set the summary file
    summary_file = 'summary.txt'

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
            script_file_id.write( '{0}\n'.format('HISAT2_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_hisat2_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$HISAT2_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_hisat2_bioconda_code())))
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
                script_file_id.write( '{0}\n'.format('function build_hisat2_index'))
                script_file_id.write( '{\n')
                script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
                script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        hisat2-build \\'))
                script_file_id.write( '{0}\n'.format('            -p {0} \\'.format(threads)))
                script_file_id.write( '{0}\n'.format('            -f \\'))
                if large_index.upper() == 'YES':
                    script_file_id.write( '{0}\n'.format('            --large-index \\'))
                if splice_site_file.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --ss {0} \\'.format(cluster_splice_site_file)))
                if exon_file.upper() != 'NONE':
                    script_file_id.write( '{0}\n'.format('            --exon {0} \\'.format(cluster_exon_file)))
                script_file_id.write( '{0}\n'.format('            {0} \\'.format(cluster_reference_file)))
                script_file_id.write( '{0}\n'.format('            {0}'.format(hisat2_index_basename)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error hisat2-build $RC; fi'))
                script_file_id.write( '{0}\n'.format('    mkdir --parents {0}'.format(hisat2_index_dir)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error mkdir $RC; fi'))
                script_file_id.write( '{0}\n'.format('    mv -f {0}.* {1}'.format(hisat2_index_basename, hisat2_index_dir)))
                script_file_id.write( '{0}\n'.format('    RC=$?'))
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error mv $RC; fi'))
                script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_hisat2_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    hisat2 --version'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        hisat2 \\'))
            script_file_id.write( '{0}\n'.format('            --threads {0} \\'.format(threads)))
            script_file_id.write( '{0}\n'.format('            --mp {0},{1} \\'.format(max_mp, min_mp)))
            if no_softclip.upper() != 'YES':
                script_file_id.write( '{0}\n'.format('            --no-softclip \\'))
            else:
                script_file_id.write( '{0}\n'.format('            --sp {0},{1} \\'.format(max_sp, min_sp)))
            script_file_id.write( '{0}\n'.format('            --np {0} \\'.format(np)))
            script_file_id.write( '{0}\n'.format('            --rdg {0},{1} \\'.format(open_rdg, extend_rdg)))
            script_file_id.write( '{0}\n'.format('            --rfg {0},{1} \\'.format(open_rfg, extend_rfg)))
            script_file_id.write( '{0}\n'.format('            --{0} \\'.format(orientation.lower())))
            if other_parameters.upper() != 'NONE':
                parameter_list = [x.strip() for x in other_parameters.split(';')]
                for i in range(len(parameter_list)):
                    if parameter_list[i].find('=') > 0:
                        pattern = r'^--(.+)=(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        parameter_value = mo.group(2).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} {1} \\'.format(parameter_name, parameter_value)))
                    else:
                        pattern = r'^--(.+)$'
                        mo = re.search(pattern, parameter_list[i])
                        parameter_name = mo.group(1).strip()
                        script_file_id.write( '{0}\n'.format('            --{0} \\'.format(parameter_name)))
            if splice_site_file.upper() != 'NONE':
                script_file_id.write( '{0}\n'.format('            --known-splicesite-infile {0} \\'.format(cluster_splice_site_file)))
            script_file_id.write( '{0}\n'.format('            -x {0}/{1} \\'.format(hisat2_index_dir, hisat2_index_basename)))
            if format.upper() == 'FASTQ':
                script_file_id.write( '{0}\n'.format('            -q \\'))
            elif format.upper() == 'FASTA':
                script_file_id.write( '{0}\n'.format('            -f \\'))
            if read_type.upper() == 'SE':
                script_file_id.write( '{0}\n'.format('            -U {0}'.format(','.join(read_file_1_list))))
            elif read_type.upper() == 'PE':
                script_file_id.write( '{0}\n'.format('            -1 {0} \\'.format(','.join(read_file_1_list))))
                script_file_id.write( '{0}\n'.format('            -2 {0} \\'.format(','.join(read_file_2_list))))
            script_file_id.write( '{0}\n'.format('            -S {0} \\'.format(alignment_file)))
            script_file_id.write( '{0}\n'.format('            --un-gz {0} \\'.format(un_gz)))
            script_file_id.write( '{0}\n'.format('            --al-gz {0} \\'.format(al_gz)))
            script_file_id.write( '{0}\n'.format('            --un-conc-gz {0} \\'.format(un_conc_gz)))
            script_file_id.write( '{0}\n'.format('            --al-conc-gz {0} \\'.format(al_conc_gz)))
            script_file_id.write( '{0}\n'.format('            --summary-file {0} \\'.format(summary_file)))
            script_file_id.write( '{0}\n'.format('            --time'))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error hisat2 $RC; fi'))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_hisat2_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_hisat2_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_hisat2_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_hisat2_name(), cluster_name))))
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
                script_file_id.write( '{0}\n'.format('build_hisat2_index'))
            script_file_id.write( '{0}\n'.format('run_hisat2_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_hisat2_process_script()))
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
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(get_hisat2_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_hisat2_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_hisat2_config_file():
    '''
    Get the HISAT2 config file path.
    '''

    # assign the HISAT2 config file path
    hisat2_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_hisat2_code())

    # return the HISAT2 config file path
    return hisat2_config_file

#-------------------------------------------------------------------------------

def get_hisat2_process_script():
    '''
    Get the HISAT2 process script path in the local computer.
    '''

    # assign the HISAT2 script path
    hisat2_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_hisat2_code())

    # return the HISAT2 script path
    return hisat2_process_script

#-------------------------------------------------------------------------------

def get_hisat2_process_starter():
    '''
    Get the HISAT2 process starter path in the local computer.
    '''

    # assign the HISAT2 process starter path
    hisat2_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_hisat2_code())

    # return the HISAT2 starter path
    return hisat2_process_starter

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
     print('This file contains functions related to the HISAT2 process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
