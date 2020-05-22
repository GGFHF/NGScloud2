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
This file contains functions related to the SOAPdenovo-Trans process used in
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

def create_soapdenovotrans_config_file(experiment_id='exp001', read_dataset_id=xlib.get_uploaded_read_dataset_name(), read_type='PE', file_1_list=['rnaseq-a_1.fastq'], file_2_list=['rnaseq-a_2.fastq']):
    '''
    Create SOAPdenovo-Trans config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the SOAPdenovo-Trans config file path
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # create the SOAPdenovo-Trans config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(soapdenovotrans_config_file)):
            os.makedirs(os.path.dirname(soapdenovotrans_config_file))
        with open(soapdenovotrans_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current run.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# The files have to be located in the cluster directory {0}/experiment_id/read_dataset_id'.format(xlib.get_cluster_read_dir())))
            file_id.write( '{0}\n'.format('# The experiment_id and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# You can consult the parameters of SOAPdenovo-Trans and their meaning in http://soap.genomics.org.cn/SOAPdenovo-Trans.html.'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# There are two formats to set an option:'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value                             <- if the option supports a single value'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('#    option = value-1, value-2, ..., value-n    <- if the option supports a values list'))
            file_id.write( '{0}\n'.format('#'))
            file_id.write( '{0}\n'.format('# WARNING: The files have to be decompressed.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information that identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('read_dataset_id = {0}'.format(read_dataset_id), '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the SOAPdenovo-Trans parameters'))
            file_id.write( '{0}\n'.format('[SOAPdenovo-Trans parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('version = 31', '# SOAPdenovo-Trans version: {0}'.format(get_version_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of cpu for use'))
            file_id.write( '{0:<50} {1}\n'.format('kmer = 25', '# value or values list of kmer size: minimum, 13; maximum: version value.'))
            file_id.write( '{0:<50} {1}\n'.format('kmer_freq_cutoff = 0', '# kmers with frequency no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format('edge_cov_cutoff = 2', '# edges with coverage no larger than the value will be deleted'))
            file_id.write( '{0:<50} {1}\n'.format('rpkm = NO', '# output assembly RPKM statistics: {0}'.format(get_rpkm_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('srkgf = NO', '# output gap related redas for SRkgf to fill gap: {0}'.format(get_srkgf_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('scaffold = NO', '# scaffold structure exists: {0}'.format(get_scaffold_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('fill = NO', '# fill gaps in scaffolds: {0}'.format(get_fill_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('merge_level = 1', '# strength of merging similar sequences during contiging: minimum, 0; maximum, 3.'))
            file_id.write( '{0:<50} {1}\n'.format('min_contig_len = 100', '# shortest contig for scaffolding'))
            file_id.write( '{0:<50} {1}\n'.format('locus_max_output = 5', '# output the number of transcripts no more than the value in one locus'))
            file_id.write( '{0:<50} {1}\n'.format('gap_len_diff = 50', '# allowed length difference between estimated and filled gap'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the global information of all libraries.'))
            file_id.write( '{0}\n'.format('[library]'))
            file_id.write( '{0:<50} {1}\n'.format('max_rd_len = 100', '# any read longer than the value will be cut to this length'))
            for i in range(len(file_1_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first library.'))
                file_id.write( '{0}\n'.format('[library-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('format = FASTQ', '# format: {0}'.format(get_format_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('read_type = {0}'.format(read_type), '# read type: {0}'.format(get_read_type_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('read_file_1 = {0}'.format(os.path.basename(file_1_list[i])), '# name of the read file in SE read type or the + strand read file in PE case'))
                if read_type == 'SE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = NONE', '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                elif read_type == 'PE':
                    file_id.write( '{0:<50} {1}\n'.format('read_file_2 = {0}'.format(os.path.basename(file_2_list[i])), '# name of the - strand reads file in PE read type or NONE in SE and SP case'))
                file_id.write( '{0:<50} {1}\n'.format('avg_ins = 200', '# average insert size'))
                file_id.write( '{0:<50} {1}\n'.format('reverse_seq = 0', '# sequence direction: {0}'.format(get_reverse_seq_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('asm_flags = 3', '# in which part(s) the reads are used: {0}'.format(get_asm_flags_code_list_text())))
                file_id.write( '{0:<50} {1}\n'.format('rd_len_cutof = 100', '# maximum read length'))
                file_id.write( '{0:<50} {1}\n'.format('map_len = 32', '# minimum aligned length to contigs for a reliable read location (minimum: 32)'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more libraries, you have to repeat the section library-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovotrans_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_soapdenovotrans_process(cluster_name, log, function=None):
    '''
    Run an experiment corresponding to the options in SOAPdenovo-Trans config file.
    '''

    # initialize the control variable
    OK = True

    # get the SOAPdenovo-Trans code and name
    soapdenovotrans_code = xlib.get_soapdenovotrans_code()
    soapdenovotrans_name = xlib.get_soapdenovotrans_name()

    # get the SOAPdenovo-Trans config file
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # get the SOAPdenovo-Trans option dictionary
    soapdenovotrans_option_dict = xlib.get_option_dict(soapdenovotrans_config_file)

    # get the experiment identification
    experiment_id = soapdenovotrans_option_dict['identification']['experiment_id']

    # get the SOAPdenovo-Trans process script path in the local computer
    soapdenovotrans_process_script = get_soapdenovotrans_process_script()

    # get the SOAPdenovo-Trans process starter path in the local computer
    soapdenovotrans_process_starter = get_soapdenovotrans_process_starter()

    # get the SOAPdenovo-Trans process config file path in the local computer
    soapdenovotrans_process_config_file = get_soapdenovotrans_process_config_file()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the SOAPdenovo-Trans configuration file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(soapdenovotrans_name))
    (OK, error_list) = check_soapdenovotrans_config_file(strict=True)
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

    # check the SOAPdenovo-Trans is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_bioconda_package(xlib.get_soapdenovotrans_bioconda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(soapdenovotrans_name))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(soapdenovotrans_name))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # build the process configuration file
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process config file {0} ...\n'.format(soapdenovotrans_process_config_file))
        (OK, error_list) = build_soapdenovotrans_process_config_file()
        if OK:
            log.write('The file is built.\n')
        else:
            log.write('*** ERROR: The file could not be built.\n')

    # for each kmer value, build the process, copy it the cluster and run it
    if OK:

        # get the kmer list
        kmer = soapdenovotrans_option_dict['SOAPdenovo-Trans parameters']['kmer']
        kmer_list = xlib.split_literal_to_integer_list(kmer)
        
        # for each kmer value, do the tasks
        i = 1
        for kmer_value in kmer_list:

            # determine the run directory in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Determining the run directory for kmer {0} in the cluster is being determined ...\n'.format(kmer_value))
            if i > 1:
                current_run_dir = '{0}-{1}'.format(xlib.get_cluster_current_run_dir(experiment_id, soapdenovotrans_code), i)
            else:
                current_run_dir = '{0}'.format(xlib.get_cluster_current_run_dir(experiment_id, soapdenovotrans_code))
            command = f'mkdir --parents {current_run_dir}'
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The directory path is {0}.\n'.format(current_run_dir))
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')
            i += 1

            # build the SOAPdenovo-Trans process script
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process script {0} ...\n'.format(soapdenovotrans_process_script))
            (OK, error_list) = build_soapdenovotrans_process_script(cluster_name, current_run_dir, kmer_value)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')
                break

            # upload the process configuration file to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process config file {0} to the directory {1} of the master ...\n'.format(soapdenovotrans_process_config_file, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_config_file))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovotrans_process_config_file, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # upload the process script to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process script {0} to the directory {1} of the master ...\n'.format(soapdenovotrans_process_script, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_script))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovotrans_process_script, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process script in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovotrans_process_script)))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_script))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # build the process starter
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the process starter {0} ...\n'.format(soapdenovotrans_process_starter))
            (OK, error_list) = build_soapdenovotrans_process_starter(current_run_dir)
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('***ERROR: The file could not be built.\n')
                break

            # upload the process starter to the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the process starter {0} to the directory {1} of the master ...\n'.format(soapdenovotrans_process_starter, current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_starter))
            (OK, error_list) = xssh.put_file(sftp_client, soapdenovotrans_process_starter, cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')
                break

            # set run permision to the process starter in the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovotrans_process_starter)))
            command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_starter))
            (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
            if OK:
                log.write('The run permision is set on.\n')
            else:
                log.write(f'*** ERROR: Wrong command ---> {command}\n')

            # submit the process
            log.write(f'{xlib.get_separator()}\n')
            log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(soapdenovotrans_process_starter)))
            OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(soapdenovotrans_process_starter), log)

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

def check_soapdenovotrans_config_file(strict):
    '''
    Check the SOAPdenovo-Trans config file checking the all the options have right values.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'

    # get the SOAPdenovo-Trans name
    soapdenovotrans_name = xlib.get_soapdenovotrans_name()

    # get the SOAPdenovo-Trans configuration file path
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # get the option dictionary
    try:
        soapdenovotrans_option_dict = xlib.get_option_dict(soapdenovotrans_config_file)
    except Exception as e:
        error_list.append('*** ERROR: The dictionary with options could  not be retrieved from the config file.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in soapdenovotrans_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = soapdenovotrans_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the sectfillion "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            run_id = soapdenovotrans_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if run_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "SOAPdenovo-Trans parameters"
        if 'SOAPdenovo-Trans parameters' not in sections_list:
            error_list.append('*** ERROR: the section "SOAPdenovo-Trans parameters" is not found.')
            OK = False
        else:

            # check section "SOAPdenovo-Trans parameters" - key "version"
            version = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('version', not_found)
            if version == not_found:
                error_list.append('*** ERROR: the key "version" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_code(version, get_version_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "version" has to be {0}.'.format(get_version_code_list_text()))
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "ncpu"
            ncpu = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=1):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 1.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "kmer"
            kmer = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('kmer', not_found)
            if kmer == not_found:
                error_list.append('*** ERROR: the key "kmer" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            else:
                kmer_list = xlib.split_literal_to_integer_list(kmer)
                if kmer_list == []:
                    error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 13 and the version value.')
                    OK = False
                else:
                    for kmer_item in kmer_list:
                        if not xlib.check_int(kmer_item, minimum=13, maximum=int(version)):
                            error_list.append('*** ERROR: the key "kmer" has to be an integer number, or an integer number list, between 13 and the version value.')
                            OK = False
                            break

            # check section "SOAPdenovo-Trans parameters" - key "kmer_freq_cutoff"
            kmer_freq_cutoff = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('kmer_freq_cutoff', not_found)
            if kmer_freq_cutoff == not_found:
                error_list.append('*** ERROR: the key "kmer_freq_cutoff" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(kmer_freq_cutoff, minimum=0):
                error_list.append('*** ERROR: the key "kmer_freq_cutoff" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "edge_cov_cutoff"
            edge_cov_cutoff = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('edge_cov_cutoff', not_found)
            if edge_cov_cutoff == not_found:
                error_list.append('*** ERROR: the key "edge_cov_cutoff" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(edge_cov_cutoff, minimum=0):
                error_list.append('*** ERROR: the key "edge_cov_cutoff" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "rpkm"
            rpkm = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('rpkm', not_found)
            if rpkm == not_found:
                error_list.append('*** ERROR: the key "rpkm" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_code(rpkm, get_rpkm_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "rpkm" has to be {0}.'.format(get_rpkm_code_list_text()))
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "srkgf"
            srkgf = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('srkgf', not_found)
            if srkgf == not_found:
                error_list.append('*** ERROR: the key "srkgf" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_code(srkgf, get_srkgf_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "srkgf" has to be {0}.'.format(get_srkgf_code_list_text()))
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "scaffold"
            scaffold = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('scaffold', not_found)
            if scaffold == not_found:
                error_list.append('*** ERROR: the key "scaffold" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_code(scaffold, get_scaffold_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "scaffold" has to be {0}.'.format(get_scaffold_code_list_text()))
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "fill"
            fill = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('fill', not_found)
            if fill == not_found:
                error_list.append('*** ERROR: the key "fill" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_code(fill, get_fill_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "fill" has to be {0}.'.format(get_fill_code_list_text()))
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "merge_level"
            merge_level = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('merge_level', not_found)
            if merge_level == not_found:
                error_list.append('*** ERROR: the key "merge_level" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(merge_level, minimum=0, maximum=3):
                error_list.append('*** ERROR: the key "merge_level" has to be an integer number between 0 and 3.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "min_contig_len"
            min_contig_len = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('min_contig_len', not_found)
            if min_contig_len == not_found:
                error_list.append('*** ERROR: the key "min_contig_len" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(min_contig_len, minimum=0):
                error_list.append('*** ERROR: the key "min_contig_len" is not an integer value greater or equal to 0.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "locus_max_output"
            locus_max_output = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('locus_max_output', not_found)
            if locus_max_output == not_found:
                error_list.append('*** ERROR: the key "locus_max_output" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(locus_max_output, minimum=0):
                error_list.append('*** ERROR: the key "locus_max_output" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "SOAPdenovo-Trans parameters" - key "gap_len_diff"
            gap_len_diff = soapdenovotrans_option_dict.get('SOAPdenovo-Trans parameters', {}).get('gap_len_diff', not_found)
            if gap_len_diff == not_found:
                error_list.append('*** ERROR: the key "gap_len_diff" is not found in the section "SOAPdenovo-Trans parameters".')
                OK = False
            elif not xlib.check_int(gap_len_diff, minimum=0):
                    error_list.append('*** ERROR: the key "gap_len_diff" has to be an integer number greater than or equal to 0.')
                    OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "max_rd_len"
            max_rd_len = soapdenovotrans_option_dict.get('library', {}).get('max_rd_len', not_found)
            if max_rd_len == not_found:
                error_list.append('*** ERROR: the key "max_rd_len" is not found in the section "library".')
                OK = False
            elif not xlib.check_int(max_rd_len, minimum=0):
                error_list.append('*** ERROR: the key "max_rd_len" has to be an integer number greater than or equal to 0.')
                OK = False

        # check section "library-1"
        if 'library-1' not in sections_list:
            error_list.append('*** ERROR: the section "library-1" is not found.')
            OK = False

        # check all sections "library-n"
        for section in sections_list:

            if section not in ['identification', 'SOAPdenovo-Trans parameters', 'library']:

                # check than the section identification is like library-n 
                if not re.match('^library-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "library-n" - key "format"
                    format = soapdenovotrans_option_dict.get(section, {}).get('format', not_found)
                    if format == not_found:
                        error_list.append('*** ERROR: the key "format" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(format, get_format_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "format" has to be {0}.'.format(get_format_code_list_text()))
                        OK = False

                    # check section "library-n" - key "read_type"
                    read_type = soapdenovotrans_option_dict.get(section, {}).get('read_type', not_found)
                    if read_type == not_found:
                        error_list.append('*** ERROR: the key "read_type" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(read_type, get_read_type_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "read_type" has to be {0}.'.format(get_read_type_code_list_text()))
                        OK = False
                    elif read_type == 'SP' and format != 'FASTA':
                        error_list.append('*** ERROR: if read type is SP, the format has to be FASTA.')
                        OK = False

                    # check section "library-n" - key "read_file_1"
                    read_file_1 = soapdenovotrans_option_dict.get(section, {}).get('read_file_1', not_found)
                    if read_file_1 == not_found:
                        error_list.append('*** ERROR: the key "read_file_1" is not found in the section "{0}"'.format(section))
                        OK = False
                    elif read_file_1.find('.gz') != -1:
                        error_list.append('*** ERROR: the key "read_file_1" in the section "{0}" has to be a decompressed file.'.format(section))
                        OK = False

                    # check section "library-n" - key "read_file_2"
                    read_file_2 = soapdenovotrans_option_dict.get(section, {}).get('read_file_2', not_found)
                    if read_file_2 == not_found:
                        error_list.append('*** ERROR: the key "read_file_2" is not found in the section "{0}"'.format(section))
                        OK = False
                    elif read_file_2.find('.gz') != -1:
                        error_list.append('*** ERROR: the key "read_file_2" in the section "{0}" has to be a decompressed file.'.format(section))
                        OK = False

                    # check section "library" - key "avg_ins"
                    avg_ins = soapdenovotrans_option_dict.get(section, {}).get('avg_ins', not_found)
                    if avg_ins == not_found:
                        error_list.append('*** ERROR: the key "avg_ins" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(avg_ins, minimum=0):
                        error_list.append('*** ERROR: the key "avg_ins" in the section "{0}" has to be an integer number greater than or equal to 0.'.format(section))
                        OK = False

                    # check section "library-n" - key "reverse_seq"
                    reverse_seq = soapdenovotrans_option_dict.get(section, {}).get('reverse_seq', not_found)
                    if reverse_seq == not_found:
                        error_list.append('*** ERROR: the key "reverse_seq" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(reverse_seq, get_reverse_seq_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "reverse_seq" has to be {0}.'.format(get_reverse_seq_code_list_text()))
                        OK = False

                    # check section "library-n" - key "asm_flags"
                    asm_flags = soapdenovotrans_option_dict.get(section, {}).get('asm_flags', not_found)
                    if asm_flags == not_found:
                        error_list.append('*** ERROR: the key "asm_flags" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_code(asm_flags, get_asm_flags_code_list(), case_sensitive=False):
                        error_list.append('*** ERROR: the key "asm_flags" has to be {0}.'.format(get_asm_flags_code_list_text()))
                        OK = False

                    # check section "library-n" - key "rd_len_cutof"
                    rd_len_cutof = soapdenovotrans_option_dict.get(section, {}).get('rd_len_cutof', not_found)
                    if rd_len_cutof == not_found:
                        error_list.append('*** ERROR: the key "rd_len_cutof" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(rd_len_cutof, minimum=0):
                        error_list.append('*** ERROR: the key "rd_len_cutof" in the section "{0}" has to be an integer number greater than or equal to 1.'.format(section))
                        OK = False

                    # check section "library-n" - key "map_len"
                    map_len = soapdenovotrans_option_dict.get(section, {}).get('map_len', not_found)
                    if map_len == not_found:
                        error_list.append('*** ERROR: the key "map_len" is not found in the section "{0}".'.format(section))
                        OK = False
                    elif not xlib.check_int(map_len, minimum=32):
                        error_list.append('*** ERROR: the key "map_len" in the section "{0}" has to be an integer number greater than or equal to 32.'.format(section))
                        OK = False

    # warn that the SOAPdenovo-Trans config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(soapdenovotrans_name))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovotrans_process_script(cluster_name, current_run_dir, kmer_value):
    '''
    Build the current SOAPdenovo-Trans process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the SOAPdenovo-Trans config file path
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # get the option dictionary
    soapdenovotrans_options_dict = xlib.get_option_dict(soapdenovotrans_config_file)

    # get the options
    experiment_id = soapdenovotrans_options_dict['identification']['experiment_id']
    version = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['version']
    rpkm = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['rpkm'].upper()
    srkgf = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['srkgf'].upper()
    scaffold = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['scaffold'].upper()
    fill = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['fill'].upper()
    ncpu = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['ncpu']
    kmer_freq_cutoff = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['kmer_freq_cutoff']
    edge_cov_cutoff = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['edge_cov_cutoff']
    merge_level = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['merge_level']
    min_contig_len = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['min_contig_len']
    locus_max_output = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['locus_max_output']
    gap_len_diff = soapdenovotrans_options_dict['SOAPdenovo-Trans parameters']['gap_len_diff']

    # get the SOAPdenovo-Trans process config file name
    soapdenovotrans_process_config_file = get_soapdenovotrans_process_config_file()

    # get the SOAPdenovo-Trans process script name
    soapdenovotrans_process_script = get_soapdenovotrans_process_script()

    # write the SOAPdenovo-Trans process script
    try:
        if not os.path.exists(os.path.dirname(soapdenovotrans_process_script)):
            os.makedirs(os.path.dirname(soapdenovotrans_process_script))
        with open(soapdenovotrans_process_script, mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('ulimit -s unlimited'))
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( '{0}\n'.format('SOAPDENOVOTRANS_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_soapdenovotrans_bioconda_code())))
            script_file_id.write( '{0}\n'.format('PATH=$SOAPDENOVOTRANS_PATH:$PATH'))
            script_file_id.write( '{0}\n'.format('cd {0}/{1}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('source activate {0}'.format(xlib.get_soapdenovotrans_bioconda_code())))
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
            script_file_id.write( '{0}\n'.format('function run_soapdenovotrans_process'))
            script_file_id.write( '{\n')
            script_file_id.write( '{0}\n'.format('    cd {0}'.format(current_run_dir)))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    /usr/bin/time \\'))
            script_file_id.write( '{0}\n'.format('        --format="$SEP\\nElapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
            script_file_id.write( '{0}\n'.format('        SOAPdenovo-Trans-{0}mer all \\'.format(version)))
            script_file_id.write( '{0}\n'.format('            -s {0}/{1} \\'.format(current_run_dir, os.path.basename(soapdenovotrans_process_config_file))))
            script_file_id.write( '{0}\n'.format('            -o {0}-{1} \\'.format(experiment_id, os.path.basename(current_run_dir))))
            script_file_id.write( '{0}\n'.format('            -p {0} \\'.format(ncpu)))
            script_file_id.write( '{0}\n'.format('            -K {0} \\'.format(kmer_value)))
            script_file_id.write( '{0}\n'.format('            -d {0} \\'.format(kmer_freq_cutoff)))
            script_file_id.write( '{0}\n'.format('            -e {0} \\'.format(edge_cov_cutoff)))
            if rpkm == 'YES':
                script_file_id.write( '{0}\n'.format('            -R \\'))
            if srkgf == 'YES':
                script_file_id.write( '{0}\n'.format('            -f \\'))
            if scaffold == 'YES':
                script_file_id.write( '{0}\n'.format('            -S \\'))
            if fill == 'YES':
                script_file_id.write( '{0}\n'.format('            -F \\'))
            script_file_id.write( '{0}\n'.format('            -M {0} \\'.format(merge_level)))
            script_file_id.write( '{0}\n'.format('            -L {0} \\'.format(min_contig_len)))
            script_file_id.write( '{0}\n'.format('            -t {0} \\'.format(locus_max_output)))
            script_file_id.write( '{0}\n'.format('            -G {0}'.format(gap_len_diff)))
            script_file_id.write( '{0}\n'.format('    RC=$?'))
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error SOAPdenovo-Trans-{0}mer $RC; fi'.format(version)))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_soapdenovotrans_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_ok(xlib.get_soapdenovotrans_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('    SUBJECT="{0}: {1} process"'.format(xlib.get_project_name(), xlib.get_soapdenovotrans_name())))
            script_file_id.write( '{0}\n'.format('    MESSAGE="{0}"'.format(xlib.get_mail_message_wrong(xlib.get_soapdenovotrans_name(), cluster_name))))
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
            script_file_id.write( '{0}\n'.format('run_soapdenovotrans_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovotrans_process_script))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovotrans_process_starter(current_run_dir):
    '''
    Build the starter of the current SOAPdenovo-Trans process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the SOAPdenovo-Trans process starter path in local
    soapdenovotrans_process_starter = get_soapdenovotrans_process_starter()

    # get the SOAPdenovo-Trans process script path in local
    soapdenovotrans_process_script = get_soapdenovotrans_process_script()

    # get the log file name
    log_file = xlib.get_cluster_log_file()

    # write the SOAPdenovo-Trans process starter
    try:
        if not os.path.exists(os.path.dirname(soapdenovotrans_process_starter)):
            os.makedirs(os.path.dirname(soapdenovotrans_process_starter))
        with open(soapdenovotrans_process_starter, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#!/bin/bash'))
            file_id.write( '{0}\n'.format('#-------------------------------------------------------------------------------'))
            file_id.write( '{0}\n'.format('{0}/{1} &>{0}/{2}'.format(current_run_dir, os.path.basename(soapdenovotrans_process_script), log_file)))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovotrans_process_starter))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_soapdenovotrans_process_config_file():
    '''
    Build the SOAPdenovo-Trans process config file to the current SOPAdenovo-Trans experiment.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the SOAPdenovo-Trans config file path
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # get the SOAPdenovo-Trans process config file path
    soapdenovotrans_process_config_file = get_soapdenovotrans_process_config_file()

    # get the option dictionary
    soapdenovotrans_option_dict = xlib.get_option_dict(soapdenovotrans_config_file)

    # get the sections list
    sections_list = []
    for section in soapdenovotrans_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # get the experiment identification and read dataset identification
    experiment_id = soapdenovotrans_option_dict['identification']['experiment_id']
    read_dataset_id = soapdenovotrans_option_dict['identification']['read_dataset_id']

    # write the SOAPdenovo-Trans process configuration file
    try:
        if not os.path.exists(os.path.dirname(soapdenovotrans_process_config_file)):
            os.makedirs(os.path.dirname(soapdenovotrans_process_config_file))
        with open(soapdenovotrans_process_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('#maximal read length'))
            file_id.write( '{0}\n'.format('max_rd_len={0}'.format(soapdenovotrans_option_dict['library']['max_rd_len'])))
            for section in sections_list:
                if re.match('^library-[0-9]+$', section):
                    format = soapdenovotrans_option_dict[section]['format'].upper()
                    read_type = soapdenovotrans_option_dict[section]['read_type'].upper()
                    read_file_1 = soapdenovotrans_option_dict[section]['read_file_1']
                    read_file_1 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_1)
                    if read_type == 'PE':
                        read_file_2 = soapdenovotrans_option_dict[section]['read_file_2']
                        read_file_2 = xlib.get_cluster_read_file(experiment_id, read_dataset_id, read_file_2)
                    file_id.write( '{0}\n'.format('[LIB]'))
                    file_id.write( '{0}\n'.format('# maximal read length in this lib'))
                    file_id.write( '{0}\n'.format('rd_len_cutof={0}'.format(soapdenovotrans_option_dict[section]['rd_len_cutof'])))
                    file_id.write( '{0}\n'.format('# average insert size'))
                    file_id.write( '{0}\n'.format('avg_ins={0}'.format(soapdenovotrans_option_dict[section]['avg_ins'])))
                    file_id.write( '{0}\n'.format('# if sequence needs to be reversed'))
                    file_id.write( '{0}\n'.format('reverse_seq={0}'.format(soapdenovotrans_option_dict[section]['reverse_seq'])))
                    file_id.write( '{0}\n'.format('# in which part(s) the reads are used'))
                    file_id.write( '{0}\n'.format('asm_flags={0}'.format(soapdenovotrans_option_dict[section]['asm_flags'])))
                    file_id.write( '{0}\n'.format('# minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)'))
                    file_id.write( '{0}\n'.format('map_len={0}'.format(soapdenovotrans_option_dict[section]['map_len'])))
                    if format == 'FASTA' and read_type == 'SE':
                        file_id.write( '{0}\n'.format('# fasta file for single reads'))
                        file_id.write( '{0}\n'.format('f={0}'.format(read_file_1)))
                    elif format == 'FASTA' and read_type == 'PE':
                        file_id.write( '{0}\n'.format('# fasta file for read 1'))
                        file_id.write( '{0}\n'.format('f1={0}'.format(read_file_1)))
                        file_id.write( '{0}\n'.format('# fastq file for read 2 always follows fastq file for read 1'))
                        file_id.write( '{0}\n'.format('f2={0}'.format(read_file_2)))
                    elif format == 'FASTA' and read_type == 'SP':
                        file_id.write( '{0}\n'.format('# a single fasta file for paired reads'))
                        file_id.write( '{0}\n'.format('p={0}'.format(read_file_1)))
                    elif format == 'FASTQ' and read_type == 'SE':
                        file_id.write( '{0}\n'.format('# fastq file for single reads'))
                        file_id.write( '{0}\n'.format('q={0}'.format(read_file_1)))
                    elif format == 'FASTQ' and read_type == 'PE':
                        file_id.write( '{0}\n'.format('# fastq file for read 1'))
                        file_id.write( '{0}\n'.format('q1={0}'.format(read_file_1)))
                        file_id.write( '{0}\n'.format('# fastq file for read 2 always follows fastq file for read 1'))
                        file_id.write( '{0}\n'.format('q2={0}'.format(read_file_2)))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(soapdenovotrans_process_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def determine_soapdenovotrans_cluster():
    '''
    Determine the cluster to the current SOPAdenovo-Trans experiment.
    '''

    # initialize the template and cluster names
    template_name = ''
    cluster_name = ''

    # get the SOAPdenovo-Trans config file path
    soapdenovotrans_config_file = get_soapdenovotrans_config_file()

    # get the option dictionary
    soapdenovotrans_options_dict = xlib.get_option_dict(soapdenovotrans_config_file)

    # determine the template and cluster names
    template_name = soapdenovotrans_options_dict['performance']['template_name']
    if template_name == 'AUTO': 
        template_name = 'cl-t1.micro'
        cluster_name = template_name
    else: 
        cluster_name = template_name

    # return the template and cluster names
    return (template_name, cluster_name)

#-------------------------------------------------------------------------------

def get_soapdenovotrans_config_file():
    '''
    Get the SOAPdenovo-Trans config file path.
    '''

    # assign the SOAPdenovo-Trans config file path
    soapdenovotrans_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_soapdenovotrans_code())

    # return the SOAPdenovo-Trans config file path
    return soapdenovotrans_config_file

#-------------------------------------------------------------------------------

def get_soapdenovotrans_process_script():
    '''
    Get the SOAPdenovo-Trans process script path in the local computer.
    '''

    # assign the SOAPdenovo-Trans script path
    soapdenovotrans_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_soapdenovotrans_code())

    # return the SOAPdenovo-Trans script path
    return soapdenovotrans_process_script

#-------------------------------------------------------------------------------

def get_soapdenovotrans_process_starter():
    '''
    Get the SOAPdenovo-Trans process starter path in the local computer.
    '''

    # assign the SOAPdenovo-Trans process starter path
    soapdenovotrans_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_soapdenovotrans_code())

    # return the SOAPdenovo-Trans starter path
    return soapdenovotrans_process_starter

#-------------------------------------------------------------------------------

def get_soapdenovotrans_process_config_file():
    '''
    Get the SOAPdenovo-Trans process config file path in the local computer.
    '''

    # assign the SOAPdenovo-Trans process config file path
    soapdenovotrans_process_config_file = '{0}/{1}-process-config.txt'.format(xlib.get_temp_dir(), xlib.get_soapdenovotrans_code())

    # return the SOAPdenovo-Trans process config file path
    return soapdenovotrans_process_config_file

#-------------------------------------------------------------------------------
    
def get_version_code_list():
    '''
    Get the code list of "version".
    '''

    return ['31', '127']

#-------------------------------------------------------------------------------
    
def get_version_code_list_text():
    '''
    Get the code list of "version" as text.
    '''

    return '31 (SOAPdenovo-Trans-31mer) or 127 (SOAPdenovo-Trans-127mer)'

#-------------------------------------------------------------------------------
    
def get_rpkm_code_list():
    '''
    Get the code list of "rpkm".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_rpkm_code_list_text():
    '''
    Get the code list of "rpkm" as text.
    '''

    return str(get_rpkm_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_srkgf_code_list():
    '''
    Get the code list of "srkgf".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_srkgf_code_list_text():
    '''
    Get the code list of "srkgf" as text.
    '''

    return str(get_srkgf_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_scaffold_code_list():
    '''
    Get the code list of "scaffold".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_scaffold_code_list_text():
    '''
    Get the code list of "scaffold" as text.
    '''

    return str(get_scaffold_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_fill_code_list():
    '''
    Get the code list of "fill".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_fill_code_list_text():
    '''
    Get the code list of "fill" as text.
    '''

    return str(get_fill_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_reverse_seq_code_list():
    '''
    Get the code list of "reverse_seq".
    '''

    return ['0', '1']

#-------------------------------------------------------------------------------
    
def get_reverse_seq_code_list_text():
    '''
    Get the code list of "reverse_seq" as text.
    '''

    return '0 (forward-reverse) or 1 (reverse-forward)'

#-------------------------------------------------------------------------------
    
def get_asm_flags_code_list():
    '''
    Get the code list of "asm_flags".
    '''

    return ['1', '2', '3']

#-------------------------------------------------------------------------------
    
def get_asm_flags_code_list_text():
    '''
    Get the code list of "asm_flags" as text.
    '''

    return '1 (only contig assembly) or 2 (only scaffolod assembly) or 3 (both contig and scffold assembly)'

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

    return ['SE', 'PE', 'SP']

#-------------------------------------------------------------------------------
    
def get_read_type_code_list_text():
    '''
    Get the code list of "read_type" as text.
    '''

    return 'SE (single-end) or PE (pair-end) or SP (single fasta file for paired reads)'

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the SOAPdenovo-Trans process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
