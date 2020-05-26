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
This file contains functions related to the ipyrad process used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys

import xbioinfoapp
import xconfiguration
import xddradseqtools
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_ipyrad_config_file(experiment_id='exp001', assembly_method='DENOVO', reference_dataset_id='NONE', reference_file='NONE', datatype='DDRAD', enzyme1='EcoRI', enzyme2='MseI',  read_dataset_id=xlib.get_uploaded_read_dataset_name(), file_pattern='rnaseq-a_1.fastq', are_data_demultiplexed='YES'):
    '''
    Create ipyrad config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the ipyrad config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_ipyrad_config_file())):
            os.makedirs(os.path.dirname(get_ipyrad_config_file()))
        with open(get_ipyrad_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '# You must review the information of this file and update the values with the corresponding ones to the current run.\n')
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# The reference and GTF files have to be located in the cluster directory {0}/experiment_id/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write(f'# The read files have to be located in the cluster directory {xlib.get_cluster_read_dir()}/experiment_id/read_dataset_id\n')
            file_id.write( '{0}\n'.format('# The experiment_id, reference_dataset_id, reference_file_name and read_dataset_id names are fixed in the identification section.'))
            file_id.write( '#\n')
            file_id.write( '{0}\n'.format('# You can consult the parameters of ipyrad and their meaning in "https://ipyrad.readthedocs.io/".'))
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification or NONE'))
            file_id.write( '{0:<50} {1}\n'.format('reference_file = {0}'.format(reference_file), '# reference file name or NONE'))
            file_id.write( '{0:<50} {1}\n'.format(f'read_dataset_id = {read_dataset_id}', '# read dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the ipyrad parameters'))
            file_id.write( '{0}\n'.format('[ipyrad parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('last_step = 7', '# last step to do'))
            file_id.write( '{0:<50} {1}\n'.format('ncpu = 4', '# number of CPUs to use; with 0, all CPUs will be used'))
            file_id.write( '{0:<50} {1}\n'.format('assembly_method = {0}'.format(assembly_method), '# assembly method: {0}'.format(get_assembly_method_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('datatype = {0}'.format(datatype), '# data type: {0}'.format(get_datatype_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('enzyme1 = {0}'.format(enzyme1), '# id of 1st restriction enzyme used in rsfile or its restriction site sequence'))
            file_id.write( '{0:<50} {1}\n'.format('enzyme2 = {0}'.format(enzyme2), '# id of 2nd restriction enzyme used in rsfile or its restriction site sequence (in DDRAD or PAIRDDRAD) or NONE (in any other case)'))
            file_id.write( '{0:<50} {1}\n'.format('max_low_qual_bases = 5', '# maximum low quality base calls (Q<20) in a read'))
            file_id.write( '{0:<50} {1}\n'.format('phred_qscore_offset = 33', '# phred Q score offset'))
            file_id.write( '{0:<50} {1}\n'.format('mindepth_statistical = 6', '# minimum depth for statistical base calling'))
            file_id.write( '{0:<50} {1}\n'.format('mindepth_majrule = 6', '# minimum depth for majority-rule base calling'))
            file_id.write( '{0:<50} {1}\n'.format('maxdepth = 10000', '# maximum cluster depth within samples'))
            file_id.write( '{0:<50} {1}\n'.format('clust_threshold = 0.85', '# clustering threshold for de novo assembly'))
            file_id.write( '{0:<50} {1}\n'.format('max_barcode_mismatch = 0', '# maximum number of allowable mismatches in barcodes'))
            file_id.write( '{0:<50} {1}\n'.format('filter_adapters = 0', '# filter for adapters/primers: {0}'.format(get_filter_adapters_code_list_text())))
            file_id.write( '{0:<50} {1}\n'.format('filter_min_trim_len = 35', '# minimum length of reads after adapter trim'))
            file_id.write( '{0:<50} {1}\n'.format('max_alleles_consens = 2', '# maximum alleles per site in consensus sequences'))
            file_id.write( '{0:<50} {1}\n'.format('max_ns_consens = 5,5', '# maximum Ns (uncalled bases) in consensus (R1;R2)'))
            file_id.write( '{0:<50} {1}\n'.format('max_hs_consens = 8,8', '# maximum Hs (heterozygotes) in consensus (R1, R2)'))
            file_id.write( '{0:<50} {1}\n'.format('min_samples_locus = 4', '# minimum sample number per locus for output'))
            file_id.write( '{0:<50} {1}\n'.format('max_snps_locus = 20,20', '# maximum SNP number per locus (R1, R2)'))
            file_id.write( '{0:<50} {1}\n'.format('max_indels_locus = 8,8', '# maximum indel number per locus (R1, R2)'))
            file_id.write( '{0:<50} {1}\n'.format('max_shared_hs_locus = 0.5', '# maximum heterozygous site number per locus (R1, R2)'))
            file_id.write( '{0:<50} {1}\n'.format('trim_reads = 0,0,0,0', '# trim raw read edges (R1>, <R1, R2>, <R2)'))
            file_id.write( '{0:<50} {1}\n'.format('trim_loci = 0,0,0,0', '# trim locus edges (R1>, <R1, R2>, <R2)'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information of the library (only one library is allowed)'))
            file_id.write( '[library]\n')
            file_id.write( '{0:<50} {1}\n'.format('file_pattern = {0}'.format(os.path.basename(file_pattern)), '# name of the read file in SE read type or the + strand read file in PE case'))
            file_id.write( '{0:<50} {1}\n'.format('data_demultiplexed = {0}'.format(are_data_demultiplexed), '# are data demultiplexed?: {0}'.format(get_data_demultiplexed_code_list_text())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_ipyrad_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_ipyrad_process(cluster_name, log, function=None):
    '''
    Run a ipyrad process.
    '''

    # initialize the control variable
    OK = True

    # get the ipyrad option dictionary
    ipyrad_option_dict = xlib.get_option_dict(get_ipyrad_config_file())

    # get the experiment identification and the indicator of data demultiplexed
    experiment_id = ipyrad_option_dict['identification']['experiment_id']
    data_demultiplexed = ipyrad_option_dict['library']['data_demultiplexed']

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the ipyrad config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(xlib.get_ipyrad_name()))
    (OK, error_list) = check_ipyrad_config_file(strict=True)
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

    # check ipyrad is installed
    if OK:
        (OK, error_list, is_installed) = xbioinfoapp.is_installed_conda_package(2, xlib.get_ipyrad_conda_code(), xlib.get_ipyrad_conda_code(), cluster_name, True, ssh_client)
        if OK:
            if not is_installed:
                log.write('*** ERROR: {0} is not installed.\n'.format(xlib.get_ipyrad_name()))
                OK = False
        else:
            log.write('*** ERROR: The verification of {0} installation could not be performed.\n'.format(xlib.get_ipyrad_name()))

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, xlib.get_ipyrad_code())
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the barcode file
    if OK:
        if data_demultiplexed == 'NO':
            log.write(f'{xlib.get_separator()}\n')
            log.write('Building the file {0} ...\n'.format(get_barcode_file()))
            (OK, error_list) = build_barcode_file()
            if OK:
                log.write('The file is built.\n')
            if not OK:
                log.write('*** ERROR: The file could not be built.\n')

    # upload the barcode file
    if OK:
        if data_demultiplexed == 'NO':
            log.write(f'{xlib.get_separator()}\n')
            log.write('Uploading the file {0} to the directory {1} ...\n'.format(get_barcode_file(), current_run_dir))
            cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_barcode_file()))
            (OK, error_list) = xssh.put_file(sftp_client, get_barcode_file(), cluster_path)
            if OK:
                log.write('The file is uploaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')

    # build the ipyrad process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(get_ipyrad_process_script()))
        (OK, error_list) = build_ipyrad_process_script(cluster_name, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the ipyrad process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(get_ipyrad_process_script(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ipyrad_process_script()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ipyrad_process_script(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ipyrad process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ipyrad_process_script())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ipyrad_process_script()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the ipyrad process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(get_ipyrad_process_starter()))
        (OK, error_list) = build_ipyrad_process_starter(current_run_dir)
        if OK:
            log.write('The file is built.\n')
        if not OK:
            log.write('***ERROR: The file could not be built.\n')

    # upload the ipyrad process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(get_ipyrad_process_starter(), current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(get_ipyrad_process_starter()))
        (OK, error_list) = xssh.put_file(sftp_client, get_ipyrad_process_starter(), cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the ipyrad process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ipyrad_process_starter())))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(get_ipyrad_process_starter()))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the ipyrad process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(get_ipyrad_process_starter())))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(get_ipyrad_process_starter()), log)

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

def check_ipyrad_config_file(strict):
    '''
    Check the ipyrad config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the dictionary of restriction enzymes
    (is_OK_restriction_enzyme_dict, enzyme_dict_error_list, restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()
    if enzyme_dict_error_list != []:
        error_list = error_list + enzyme_dict_error_list + ['\n']

    # get the enzime identification list
    if is_OK_restriction_enzyme_dict:
        enzyme_id_list = list(restriction_enzyme_dict.keys())

    # get the dictionary of individuals
    (is_OK_individual_dict, individual_dict_error_list, individual_dict) = xddradseqtools.get_individual_dict()
    if individual_dict_error_list != []:
        error_list = error_list + individual_dict_error_list + ['\n']

    # get the option dictionary
    ipyrad_option_dict = xlib.get_option_dict(get_ipyrad_config_file())
    try:
        ipyrad_option_dict = xlib.get_option_dict(get_ipyrad_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in ipyrad_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = ipyrad_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = ipyrad_option_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "reference_file"
            reference_file = ipyrad_option_dict.get('identification', {}).get('reference_file', not_found)
            if reference_file == not_found:
                error_list.append('*** ERROR: the key "reference_file" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "read_dataset_id"
            read_dataset_id = ipyrad_option_dict.get('identification', {}).get('read_dataset_id', not_found)
            if read_dataset_id == not_found:
                error_list.append('*** ERROR: the key "read_dataset_id" is not found in the section "identification".')
                OK = False

        # check section "ipyrad parameters"
        if 'ipyrad parameters' not in sections_list:
            error_list.append('*** ERROR: the section "ipyrad parameters" is not found.')
            OK = False
        else:

            # check section "ipyrad parameters" - key "last_step"
            last_step = ipyrad_option_dict.get('ipyrad parameters', {}).get('last_step', not_found)
            if last_step == not_found:
                error_list.append('*** ERROR: the key "last_step" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(last_step, minimum=1, maximum=7):
                error_list.append('*** ERROR: the key "last_step" has to be an integer number between 1 and 7.')
                OK = False

            # check section "ipyrad parameters" - key "ncpu"
            ncpu = ipyrad_option_dict.get('ipyrad parameters', {}).get('ncpu', not_found)
            if ncpu == not_found:
                error_list.append('*** ERROR: the key "ncpu" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(ncpu, minimum=0):
                error_list.append('*** ERROR: the key "ncpu" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "assembly_method"
            assembly_method = ipyrad_option_dict.get('ipyrad parameters', {}).get('assembly_method', not_found)
            if assembly_method == not_found:
                error_list.append('*** ERROR: the key "assembly_method" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_code(assembly_method, get_assembly_method_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "assembly_method" has to be {0}.'.format(get_assembly_method_code_list_text()))
                OK = False

            # check section "ipyrad parameters" - key "datatype"
            datatype = ipyrad_option_dict.get('ipyrad parameters', {}).get('datatype', not_found)
            is_ok_datatype = False
            if datatype == not_found:
                error_list.append('*** ERROR: the key "datatype" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_code(datatype, get_datatype_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "datatype" has to be {0}.'.format(get_datatype_code_list_text()))
                OK = False
            else:
                is_ok_datatype = True

            # check section "ipyrad parameters" - key "enzyme1"
            enzyme1 = ipyrad_option_dict.get('ipyrad parameters', {}).get('enzyme1', not_found)
            is_ok_enzyme1 = False
            if enzyme1 == not_found:
                error_list.append('*** ERROR: the key "enzyme1" is not found in the section "ipyrad parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if enzyme1 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme1, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    error_list.append('*** ERROR: the key "enzyme1" has to be a restriction enzyme id or a restriction site seq.')
                    OK = False
                else:
                    if enzyme1 in enzyme_id_list:
                        enzyme1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
                    else:
                        enzyme1_seq = enzyme1
                    is_ok_enzyme1 = True

            # check section "ipyrad parameters" - key "enzyme2"
            enzyme2 = ipyrad_option_dict.get('ipyrad parameters', {}).get('enzyme2', not_found)
            is_ok_enzyme2 = False
            if enzyme2 == not_found:
                error_list.append('*** ERROR: the key "enzyme2" is not found in the section "ipyrad parameters".')
                OK = False
            elif is_OK_restriction_enzyme_dict:
                if is_ok_datatype and datatype in ['DDRAD', 'PAIRDDRAD']:
                    if enzyme2 not in enzyme_id_list and not xlib.is_valid_sequence(enzyme2, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                        error_list.append('*** ERROR: the key "enzyme2" has to be a restriction enzyme id or a restriction site seq (in DDRAD or PAIRDDRAD) or NONE (in any other case).')
                        OK = False
                    else:
                        if enzyme2 in enzyme_id_list:
                            enzyme2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
                        else:
                            enzyme2_seq = enzyme2
                        is_ok_enzyme2 = True
                else:
                    if enzyme2.upper() != 'NONE':
                        error_list.append('*** ERROR: the key "enzyme2" has to be a restriction enzyme id or a restriction site seq (in DDRAD or PAIRDDRAD) or NONE (in any other case).')
                        OK = False
                    else:
                        is_ok_enzyme2 = True

            # check that enzyme1 has to be different from enzyme2
            if is_ok_datatype and is_ok_enzyme1 and is_ok_enzyme2 and datatype in ['DDRAD', 'PAIRDDRAD'] and enzyme1_seq == enzyme2_seq:
                error_list.append('*** ERROR: Both enzymes have the same sequence.')
                OK = False

            # check section "ipyrad parameters" - key "max_low_qual_bases"
            max_low_qual_bases = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_low_qual_bases', not_found)
            if max_low_qual_bases == not_found:
                error_list.append('*** ERROR: the key "max_low_qual_bases" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(max_low_qual_bases, minimum=0):
                error_list.append('*** ERROR: the key "max_low_qual_bases" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "phred_qscore_offset"
            phred_qscore_offset = ipyrad_option_dict.get('ipyrad parameters', {}).get('phred_qscore_offset', not_found)
            if phred_qscore_offset == not_found:
                error_list.append('*** ERROR: the key "phred_qscore_offset" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(phred_qscore_offset, minimum=0):
                error_list.append('*** ERROR: the key "phred_qscore_offset" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "mindepth_statistical"
            mindepth_statistical = ipyrad_option_dict.get('ipyrad parameters', {}).get('mindepth_statistical', not_found)
            if mindepth_statistical == not_found:
                error_list.append('*** ERROR: the key "mindepth_statistical" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(mindepth_statistical, minimum=0):
                error_list.append('*** ERROR: the key "mindepth_statistical" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "mindepth_majrule"
            mindepth_majrule = ipyrad_option_dict.get('ipyrad parameters', {}).get('mindepth_majrule', not_found)
            if mindepth_majrule == not_found:
                error_list.append('*** ERROR: the key "mindepth_majrule" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(mindepth_majrule, minimum=0):
                error_list.append('*** ERROR: the key "mindepth_majrule" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "maxdepth"
            maxdepth = ipyrad_option_dict.get('ipyrad parameters', {}).get('maxdepth', not_found)
            if maxdepth == not_found:
                error_list.append('*** ERROR: the key "maxdepth" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(maxdepth, minimum=0):
                error_list.append('*** ERROR: the key "maxdepth" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad simulation parameters" - key "clust_threshold"
            clust_threshold = ipyrad_option_dict.get('ipyrad parameters', {}).get('clust_threshold', not_found)
            if clust_threshold == not_found:
                error_list.append('*** ERROR: the key "clust_threshold" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_float(clust_threshold, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "clust_threshold" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ipyrad parameters" - key "max_barcode_mismatch"
            max_barcode_mismatch = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_barcode_mismatch', not_found)
            if max_barcode_mismatch == not_found:
                error_list.append('*** ERROR: the key "max_barcode_mismatch" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(max_barcode_mismatch, minimum=0):
                error_list.append('*** ERROR: the key "max_barcode_mismatch" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "filter_adapters"
            filter_adapters = ipyrad_option_dict.get('ipyrad parameters', {}).get('filter_adapters', not_found)
            if filter_adapters == not_found:
                error_list.append('*** ERROR: the key "filter_adapters" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_code(filter_adapters, get_filter_adapters_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "filter_adapters" has to be {0}.'.format(get_filter_adapters_code_list_text()))
                OK = False

            # check section "ipyrad parameters" - key "filter_min_trim_len"
            filter_min_trim_len = ipyrad_option_dict.get('ipyrad parameters', {}).get('filter_min_trim_len', not_found)
            if filter_min_trim_len == not_found:
                error_list.append('*** ERROR: the key "filter_min_trim_len" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(filter_min_trim_len, minimum=0):
                error_list.append('*** ERROR: the key "filter_min_trim_len" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "max_alleles_consens"
            max_alleles_consens = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_alleles_consens', not_found)
            if max_alleles_consens == not_found:
                error_list.append('*** ERROR: the key "max_alleles_consens" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(max_alleles_consens, minimum=0):
                error_list.append('*** ERROR: the key "max_alleles_consens" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "max_ns_consens"
            max_ns_consens = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_ns_consens', not_found)
            if max_ns_consens == not_found:
                error_list.append('*** ERROR: the key "max_ns_consens" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                max_ns_consens_list = xlib.split_literal_to_integer_list(max_ns_consens)
                if len(max_ns_consens_list) != 2 or min(max_ns_consens_list) < 0:
                    error_list.append('*** ERROR: the key "max_ns_consens" has to be two integer numbers greater than or equal to 0 separated by comma.')
                    OK = False

            # check section "ipyrad parameters" - key "max_hs_consens"
            max_hs_consens = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_hs_consens', not_found)
            if max_hs_consens == not_found:
                error_list.append('*** ERROR: the key "max_hs_consens" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                max_hs_consens_list = xlib.split_literal_to_integer_list(max_hs_consens)
                if len(max_hs_consens_list) != 2 or min(max_hs_consens_list) < 0:
                    error_list.append('*** ERROR: the key "max_hs_consens" has to be two integer numbers greater than or equal to 0 separated by comma.')
                    OK = False

            # check section "ipyrad parameters" - key "min_samples_locus"
            min_samples_locus = ipyrad_option_dict.get('ipyrad parameters', {}).get('min_samples_locus', not_found)
            if min_samples_locus == not_found:
                error_list.append('*** ERROR: the key "min_samples_locus" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_int(min_samples_locus, minimum=0):
                error_list.append('*** ERROR: the key "min_samples_locus" has to be an integer number greater than or equal to 0.')
                OK = False

            # check section "ipyrad parameters" - key "max_snps_locus"
            max_snps_locus = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_snps_locus', not_found)
            if max_snps_locus == not_found:
                error_list.append('*** ERROR: the key "max_snps_locus" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                max_snps_locus_list = xlib.split_literal_to_integer_list(max_snps_locus)
                if len(max_snps_locus_list) != 2 or min(max_snps_locus_list) < 1:
                    error_list.append('*** ERROR: the key "max_snps_locus" has to be two integer numbers greater than or equal to 1 separated by comma.')
                    OK = False

            # check section "ipyrad parameters" - key "max_indels_locus"
            max_indels_locus = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_indels_locus', not_found)
            if max_indels_locus == not_found:
                error_list.append('*** ERROR: the key "max_indels_locus" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                max_indels_locus_list = xlib.split_literal_to_integer_list(max_indels_locus)
                if len(max_indels_locus_list) != 2 or min(max_indels_locus_list) < 1:
                    error_list.append('*** ERROR: the key "max_indels_locus" has to be two integer numbers greater than or equal to 1 separated by comma.')
                    OK = False

            # check section "ipyrad simulation parameters" - key "max_shared_hs_locus"
            max_shared_hs_locus = ipyrad_option_dict.get('ipyrad parameters', {}).get('max_shared_hs_locus', not_found)
            if max_shared_hs_locus == not_found:
                error_list.append('*** ERROR: the key "max_shared_hs_locus" is not found in the section "ipyrad parameters".')
                OK = False
            elif not xlib.check_float(max_shared_hs_locus, minimum=0., maximum=1., mne=0., mxe=1E-12):
                error_list.append('*** ERROR: the key "max_shared_hs_locus" has to be a float number greater than or equal to 0.0 and less than 1.0.')
                OK = False

            # check section "ipyrad parameters" - key "trim_reads"
            trim_reads = ipyrad_option_dict.get('ipyrad parameters', {}).get('trim_reads', not_found)
            if trim_reads == not_found:
                error_list.append('*** ERROR: the key "trim_reads" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                trim_reads_list = xlib.split_literal_to_integer_list(trim_reads)
                if len(trim_reads_list) != 4 or min(trim_reads_list) < 0:
                    error_list.append('*** ERROR: the key "trim_reads" has to be four integer numbers greater than or equal to 0 separated by comma.')
                    OK = False

            # check section "ipyrad parameters" - key "trim_loci"
            trim_loci = ipyrad_option_dict.get('ipyrad parameters', {}).get('trim_loci', not_found)
            if trim_loci == not_found:
                error_list.append('*** ERROR: the key "trim_loci" is not found in the section "ipyrad parameters".')
                OK = False
            else:
                trim_loci_list = xlib.split_literal_to_integer_list(trim_loci)
                if len(trim_loci_list) != 4 or min(trim_loci_list) < 0:
                    error_list.append('*** ERROR: the key "trim_loci" has to be four integer numbers greater than or equal to 0 separated by comma.')
                    OK = False

        # check section "library"
        if 'library' not in sections_list:
            error_list.append('*** ERROR: the section "library" is not found.')
            OK = False
        else:

            # check section "library" - key "file_pattern"
            file_pattern = ipyrad_option_dict.get('library', {}).get('file_pattern', not_found)
            is_ok_file_pattern = False
            if file_pattern == not_found:
                error_list.append('*** ERROR: the key "file_pattern" is not found in the section "{0}"'.format(section))
                OK = False
            else:
                is_ok_file_pattern = True

            # check section "library" - key "data_demultiplexed"
            data_demultiplexed = ipyrad_option_dict.get('library', {}).get('data_demultiplexed', not_found)
            if data_demultiplexed == not_found:
                error_list.append('*** ERROR: the key "data_demultiplexed" is not found in the section "{0}"'.format(section))
                OK = False
            elif not xlib.check_code(data_demultiplexed, get_data_demultiplexed_code_list(), case_sensitive=False):
                error_list.append('*** ERROR: the key "data_demultiplexed" has to be {0}.'.format(get_data_demultiplexed_code_list_text()))
                OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_ipyrad_name()))

    # return the result of the control variables and the error list
    return (OK and is_OK_restriction_enzyme_dict and is_OK_individual_dict, error_list)

#-------------------------------------------------------------------------------

def build_barcode_file():
    '''
    Build the barcodes file from individual file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # set the pattern of the record of individual file (record format: individual_id;replicated_individual_id;population_id;index1_seq;[index2_seq])
    pattern = r'^(.+);(.+);(.+);(.+);(.*)$'

    # open the file of individuals
    try:
        individual_file_id = open(xddradseqtools.get_individual_file(), mode='r', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(xddradseqtools.get_individual_file()))
        OK = False

    # open the file of barcodes
    try:
        barcode_file_id = open(get_barcode_file(), mode='w', encoding='iso-8859-1', newline='\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be opened.'.format(get_barcode_file()))
        OK = False

    # for very individual, write his barcode
    if OK:

        # read the first record
        record = individual_file_id.readline()
        is_there_index2 = None

        # while there are records
        while record != '':

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try:
                    mo = re.search(pattern, record)
                    individual_id = mo.group(1).strip()
                    replicated_individual_id = mo.group(2).strip()
                    population_id = mo.group(3).strip()
                    index1_seq = mo.group(4).strip()
                    index2_seq = mo.group(5).strip()
                except Exception as e:
                    error_list.append('*** ERROR: There is a format error in the record "{0}".'.format(record.replace("\n", "")))
                    OK = False
                    break

                # write the barcode
                if index2_seq == '':
                    barcode_file_id.write( '{0}\t{1}\n'.format(individual_id, index1_seq))
                else:
                    barcode_file_id.write( '{0}\t{1}\t{2}\n'.format(individual_id, index1_seq, index2_seq))

            # read the next record
            record = individual_file_id.readline()

    # close the file of individuals
    individual_file_id.close()

    # close the file of barcodes
    barcode_file_id.close()

    # warn that the file of individuals is not valid if there are any errors
    if not OK:
        error_list.append('\nThe file {0} is not valid. Please, correct this file or recreate it.'.format(get_individual_file()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ipyrad_process_script(cluster_name, current_run_dir):
    '''
    Build the current ipyrad process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the dictionary of restriction enzymes
    (is_OK_restriction_enzyme_dict, enzyme_dict_error_list, restriction_enzyme_dict) = xddradseqtools.get_restriction_enzyme_dict()

    # get the enzime identification list
    enzyme_id_list = list(restriction_enzyme_dict.keys())

    # get the ipyrad option dictionary
    ipyrad_option_dict = xlib.get_option_dict(get_ipyrad_config_file())

    # get the options
    experiment_id = ipyrad_option_dict['identification']['experiment_id']
    reference_dataset_id = ipyrad_option_dict['identification']['reference_dataset_id']
    reference_file = ipyrad_option_dict['identification']['reference_file']
    read_dataset_id = ipyrad_option_dict['identification']['read_dataset_id']
    last_step = ipyrad_option_dict['ipyrad parameters']['last_step']
    ncpu = ipyrad_option_dict['ipyrad parameters']['ncpu']
    assembly_method = ipyrad_option_dict['ipyrad parameters']['assembly_method'].lower()
    datatype = ipyrad_option_dict['ipyrad parameters']['datatype'].lower()
    enzyme1 = ipyrad_option_dict['ipyrad parameters']['enzyme1']
    enzyme2 = ipyrad_option_dict['ipyrad parameters']['enzyme2']
    max_low_qual_bases = ipyrad_option_dict['ipyrad parameters']['max_low_qual_bases']
    phred_qscore_offset = ipyrad_option_dict['ipyrad parameters']['phred_qscore_offset']
    mindepth_statistical = ipyrad_option_dict['ipyrad parameters']['mindepth_statistical']
    mindepth_majrule = ipyrad_option_dict['ipyrad parameters']['mindepth_majrule']
    maxdepth = ipyrad_option_dict['ipyrad parameters']['maxdepth']
    clust_threshold = ipyrad_option_dict['ipyrad parameters']['clust_threshold']
    max_barcode_mismatch = ipyrad_option_dict['ipyrad parameters']['max_barcode_mismatch']
    filter_adapters = ipyrad_option_dict['ipyrad parameters']['filter_adapters']
    filter_min_trim_len = ipyrad_option_dict['ipyrad parameters']['filter_min_trim_len']
    max_alleles_consens = ipyrad_option_dict['ipyrad parameters']['max_alleles_consens']
    max_ns_consens = ipyrad_option_dict['ipyrad parameters']['max_ns_consens']
    max_hs_consens = ipyrad_option_dict['ipyrad parameters']['max_hs_consens']
    min_samples_locus = ipyrad_option_dict['ipyrad parameters']['min_samples_locus']
    max_snps_locus = ipyrad_option_dict['ipyrad parameters']['max_snps_locus']
    max_indels_locus = ipyrad_option_dict['ipyrad parameters']['max_indels_locus']
    max_shared_hs_locus = ipyrad_option_dict['ipyrad parameters']['max_shared_hs_locus']
    trim_reads = ipyrad_option_dict['ipyrad parameters']['trim_reads']
    trim_loci = ipyrad_option_dict['ipyrad parameters']['trim_loci']
    file_pattern = ipyrad_option_dict['library']['file_pattern']
    data_demultiplexed = ipyrad_option_dict['library']['data_demultiplexed'].upper()

    # set the parameter file
    parameter_file = 'params-{0}.txt'.format(experiment_id)

    # build steps
    firt_step = 1 if data_demultiplexed == 'NO' else 2
    steps = ''
    for i in range(1, int(last_step) + 1):
        steps += '{0}'.format(i)

    # get the file pattern path
    file_pattern = xlib.get_cluster_read_file(experiment_id, read_dataset_id, file_pattern.replace('.*','*'))

    # set the reference file path
    if reference_file.upper() != 'NONE':
        reference_file = xlib.get_cluster_reference_file(reference_dataset_id, reference_file)
    else:
        reference_file = ''

    # get the sequence of enzyme1
    if enzyme1 in enzyme_id_list:
        ressite1_seq = restriction_enzyme_dict[enzyme1]['restriction_site_seq']
    else:
        ressite1_seq = enzyme1

    # get the restriction_overhag of enzyme1
    cutsite1 = ressite1_seq.find('*')
    ressite1_lcut_seq = ressite1_seq[:cutsite1]
    ressite1_rcut_seq = ressite1_seq[cutsite1 + 1:]
    if len(ressite1_lcut_seq) >= len(ressite1_rcut_seq):
        restriction_overhang1 = xlib.get_reverse_complementary_sequence(ressite1_lcut_seq)
    else:
        restriction_overhang1 = ressite1_rcut_seq

    # get the sequence of enzyme2
    if enzyme2.upper() != 'NONE':
        if enzyme2 in enzyme_id_list:
            ressite2_seq = restriction_enzyme_dict[enzyme2]['restriction_site_seq']
        else:
            ressite2_seq = enzyme2

    # get the restriction_overhag of enzyme2
    if enzyme2.upper() != 'NONE':
        cutsite2 = ressite2_seq.find('*')
        ressite2_lcut_seq = ressite2_seq[:cutsite2]
        ressite2_rcut_seq = ressite2_seq[cutsite2 + 1:]
        if len(ressite2_lcut_seq) >= len(ressite2_rcut_seq):
            restriction_overhang2 = ressite1_lcut_seq
        else:
            restriction_overhang2 = xlib.get_reverse_complementary_sequence(ressite2_rcut_seq)
    else:
        restriction_overhang2 = ''

    # build the restriction overhang
    restriction_overhang = '{0},{1}'.format(restriction_overhang1, restriction_overhang2)

    # write the ipyrad process script
    try:
        if not os.path.exists(os.path.dirname(get_ipyrad_process_script())):
            os.makedirs(os.path.dirname(get_ipyrad_process_script()))
        with open(get_ipyrad_process_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('IPYRAD_PATH={0}/{1}/envs/{2}/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name(), xlib.get_ipyrad_code())))
            script_file_id.write( '{0}\n'.format('CUTADAPT_PATH={0}/{1}/envs/py27/bin'.format(xlib.get_cluster_app_dir(), xlib.get_miniconda3_name())))
            script_file_id.write( '{0}\n'.format('export PATH=$IPYRAD_PATH:$CUTADAPT_PATH:$PATH'))
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
            script_file_id.write( '{0}\n'.format('function create_parameter_file'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Creating parameter file ..."'))
            script_file_id.write( '{0}\n'.format('    ipyrad \\'))
            script_file_id.write( '{0}\n'.format('        -n {0}'.format(experiment_id)))
            script_file_id.write( '{0}\n'.format('    echo "The file is created."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function update_parameters'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Updating the parameters ..."'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_00="## \[0\] \[assembly_name\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_00/c\\{0}   $PARAMETRO_00" {1}'.format(experiment_id, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_01="## \[1\] \[project_dir\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_01/c\\{0}   $PARAMETRO_01" {1}'.format(current_run_dir, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            if data_demultiplexed == 'NO':
                script_file_id.write( '{0}\n'.format('    PARAMETRO_02="## \[2\] \[raw_fastq_path\]"'))
                script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_02/c\\{0}   $PARAMETRO_02" {1}'.format(file_pattern, parameter_file)))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
                script_file_id.write( '{0}\n'.format('    PARAMETRO_03="## \[3\] \[barcodes_path\]"'))
                script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_03/c\\{0}/{1}   $PARAMETRO_03" {2}'.format(current_run_dir, os.path.basename(get_barcode_file()), parameter_file)))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            elif data_demultiplexed == 'YES':
                script_file_id.write( '{0}\n'.format('    PARAMETRO_04="## \[4\] \[sorted_fastq_path\]"'))
                script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_04/c\\{0}   $PARAMETRO_04" {1}'.format(file_pattern, parameter_file)))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_05="## \[5\] \[assembly_method\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_05/c\\{0}   $PARAMETRO_05" {1}'.format(assembly_method, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_06="## \[6\] \[reference_sequence\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_06/c\\{0}   $PARAMETRO_06" {1}'.format(reference_file, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_07="## \[7\] \[datatype\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_07/c\\{0}   $PARAMETRO_07" {1}'.format(datatype, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_08="## \[8\] \[restriction_overhang\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_08/c\\{0}   $PARAMETRO_08" {1}'.format(restriction_overhang, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_09="## \[9\] \[max_low_qual_bases\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_09/c\\{0}   $PARAMETRO_09" {1}'.format(max_low_qual_bases, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_10="## \[10\] \[phred_Qscore_offset\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_10/c\\{0}   $PARAMETRO_10" {1}'.format(phred_qscore_offset, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_11="## \[11\] \[mindepth_statistical\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_11/c\\{0}   $PARAMETRO_11" {1}'.format(mindepth_statistical, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_12="## \[12\] \[mindepth_majrule\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_12/c\\{0}   $PARAMETRO_12" {1}'.format(mindepth_majrule, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_13="## \[13\] \[maxdepth\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_13/c\\{0}   $PARAMETRO_13" {1}'.format(maxdepth, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_14="## \[14\] \[clust_threshold\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_14/c\\{0}   $PARAMETRO_14" {1}'.format(clust_threshold, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_15="## \[15\] \[max_barcode_mismatch\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_15/c\\{0}   $PARAMETRO_15" {1}'.format(max_barcode_mismatch, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_16="## \[16\] \[filter_adapters\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_16/c\\{0}   $PARAMETRO_16" {1}'.format(filter_adapters, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_17="## \[17\] \[filter_min_trim_len\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_17/c\\{0}   $PARAMETRO_17" {1}'.format(filter_min_trim_len, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_18="## \[18\] \[max_alleles_consens\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_18/c\\{0}   $PARAMETRO_18" {1}'.format(max_alleles_consens, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_19="## \[19\] \[max_Ns_consens\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_19/c\\{0}   $PARAMETRO_19" {1}'.format(max_ns_consens, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_20="## \[20\] \[max_Hs_consens\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_20/c\\{0}   $PARAMETRO_20" {1}'.format(max_hs_consens, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_21="## \[21\] \[min_samples_locus\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_21/c\\{0}   $PARAMETRO_21" {1}'.format(min_samples_locus, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_22="## \[22\] \[max_SNPs_locus\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_22/c\\{0}   $PARAMETRO_22" {1}'.format(max_snps_locus, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_23="## \[23\] \[max_Indels_locus\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_23/c\\{0}   $PARAMETRO_23" {1}'.format(max_indels_locus, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_24="## \[24\] \[max_shared_Hs_locus\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_24/c\\{0}   $PARAMETRO_24" {1}'.format(max_shared_hs_locus, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_25="## \[25\] \[trim_reads\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_25/c\\{0}   $PARAMETRO_25" {1}'.format(trim_reads, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_26="## \[26\] \[trim_loci\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_26/c\\{0}   $PARAMETRO_26" {1}'.format(trim_loci, parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    PARAMETRO_27="## \[27\] \[output_formats\]"'))
            script_file_id.write( '{0}\n'.format('    sed -i "/$PARAMETRO_27/c\\{0}   $PARAMETRO_27" {1}'.format('*', parameter_file)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error sed $RC; fi'))
            script_file_id.write( '{0}\n'.format('    echo "The parameters are updated."'))
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( '{0}\n'.format('function run_ipyrad_process'))
            script_file_id.write( '{\n')
            script_file_id.write(f'    cd {current_run_dir}\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    ipyrad -v'))
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '{0}\n'.format('    echo "Assembling reads ..."'))
            script_file_id.write( '    /usr/bin/time \\\n')
            script_file_id.write(f'        --format="{xlib.get_time_output_format()}" \\\n')
            script_file_id.write( '{0}\n'.format('        ipyrad \\'))
            script_file_id.write( '{0}\n'.format('            -q \\'))
            script_file_id.write( '{0}\n'.format('            -c 0 \\'))
            script_file_id.write( '{0}\n'.format('            -p {0} \\'.format(parameter_file)))
            script_file_id.write( '{0}\n'.format('            -s {0}'.format(steps)))
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error ipyrad $RC; fi'))
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
            process_name = f'{xlib.get_ipyrad_name()} process'
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
            script_file_id.write( '{0}\n'.format('create_parameter_file'))
            script_file_id.write( '{0}\n'.format('update_parameters'))
            script_file_id.write( '{0}\n'.format('run_ipyrad_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ipyrad_process_script()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_ipyrad_process_starter(current_run_dir):
    '''
    Build the starter of the current ipyrad process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # write the ipyrad process starter
    try:
        if not os.path.exists(os.path.dirname(get_ipyrad_process_starter())):
            os.makedirs(os.path.dirname(get_ipyrad_process_starter()))
        with open(get_ipyrad_process_starter(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_ipyrad_process_script()), xlib.get_cluster_log_file())))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_ipyrad_process_starter()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_ipyrad_config_file():
    '''
    Get the ipyrad config file path.
    '''

    # assign the ipyrad config file path
    ipyrad_config_file = '{0}/{1}-config.txt'.format(xlib.get_config_dir(), xlib.get_ipyrad_code())

    # return the ipyrad config file path
    return ipyrad_config_file

#-------------------------------------------------------------------------------

def get_barcode_file():
    '''
    Get the ipyrad barcodes file path.
    '''

    # assign the ipyrad config file path
    ipyrad_barcodes_file = '{0}/{1}-barcodes.txt'.format(xlib.get_temp_dir(), xlib.get_ipyrad_code())

    # return the ipyrad config file path
    return ipyrad_barcodes_file

#-------------------------------------------------------------------------------

def get_ipyrad_process_script():
    '''
    Get the ipyrad process script path in the local computer.
    '''

    # assign the ipyrad script path
    ipyrad_process_script = '{0}/{1}-process.sh'.format(xlib.get_temp_dir(), xlib.get_ipyrad_code())

    # return the ipyrad script path
    return ipyrad_process_script

#-------------------------------------------------------------------------------

def get_ipyrad_process_starter():
    '''
    Get the ipyrad process starter path in the local computer.
    '''

    # assign the ipyrad process starter path
    ipyrad_process_starter = '{0}/{1}-process-starter.sh'.format(xlib.get_temp_dir(), xlib.get_ipyrad_code())

    # return the ipyrad starter path
    return ipyrad_process_starter

#-------------------------------------------------------------------------------
    
def get_assembly_method_code_list():
    '''
    Get the code list of "assembly_method".
    '''

    return ['DENOVO', 'REFERENCE', 'DENOVO+REFERENCE', 'DENOVO-REFERENCE']

#-------------------------------------------------------------------------------
    
def get_assembly_method_code_list_text():
    '''
    Get the code list of "assembly_method" as text.
    '''

    return str(get_assembly_method_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_datatype_code_list():
    '''
    Get the code list of "datatype".
    '''

    return ['RAD', 'DDRAD', 'PAIRDDRAD', 'GBS', 'PAIRGBS', '2BRAD', 'PAIR3RAD']

#-------------------------------------------------------------------------------
    
def get_datatype_code_list_text():
    '''
    Get the code list of "datatype" as text.
    '''

    return str(get_datatype_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------
    
def get_filter_adapters_code_list():
    '''
    Get the code list of "filter_adapters".
    '''

    return ['0', '1', '2']

#-------------------------------------------------------------------------------
    
def get_filter_adapters_code_list_text():
    '''
    Get the code list of "filter_adapters" as text.
    '''

    return '0 (no adapter filtering), 1 (filter based on quality scores) or 2 (strict filter for adapters)'

#-------------------------------------------------------------------------------
    
def get_data_demultiplexed_code_list():
    '''
    Get the code list of "data_demultiplexed".
    '''

    return ['YES', 'NO']

#-------------------------------------------------------------------------------
    
def get_data_demultiplexed_code_list_text():
    '''
    Get the code list of "data_demultiplexed" as text.
    '''

    return str(get_data_demultiplexed_code_list()).strip('[]').replace('\'','').replace(',', ' or')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the ipyrad process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
