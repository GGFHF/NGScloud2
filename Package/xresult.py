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
This file contains functions related to result datasets used in both console mode
and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import pathlib
import re
import subprocess
import sys

import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_result_transfer_config_file(experiment_id='exp001', result_dataset_id='trinity-160629-151313', status='uncompressed', selected_file_list=['./Trinity.fasta'], local_dir='./results'):
    '''
    Create o recreate the result transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the result transfer config file path
    result_transfer_config_file = get_result_transfer_config_file()

    # create the result transfer config file to download result files of a run
    try:
        if not os.path.exists(os.path.dirname(result_transfer_config_file)):
            os.makedirs(os.path.dirname(result_transfer_config_file))
        with open(result_transfer_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current transfer.'))
            file_id.write( '{0}\n'.format('# The files will be copied from the cluster directory {0}/experiment_id/result_dataset_id.'.format(xlib.get_cluster_result_dir())))
            file_id.write( '\n')
            file_id.write( '# This section has the information identifies the experiment.\n')
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format(f'experiment_id = {experiment_id}', '# experiment identification'))
            file_id.write( '{0:<50} {1}\n'.format('result_dataset_id = {0}'.format(result_dataset_id), '# run identification'))
            file_id.write( '{0:<50} {1}\n'.format('status = {0}'.format(status), '# result dataset status (it has to be always {0})'.format(status)))
            file_id.write( '{0:<50} {1}\n'.format('local_dir = {0}'.format(local_dir), '# local path where the file will be download'))
            if status == 'uncompressed':
                for i in range(len(selected_file_list)):
                    file_id.write( '\n')
                    if i == 0:
                        file_id.write( '{0}\n'.format('# This section has the information of the first result file.'))
                    file_id.write( '{0}\n'.format('[file-{0}]'.format(i + 1)))
                    file_id.write( '{0:<50} {1}\n'.format('dataset_subdirectory = {0}'.format(os.path.dirname(selected_file_list[i])), '# subdirectory of result dataset'))
                    file_id.write( '{0:<50} {1}\n'.format('file_name = {0}'.format(os.path.basename(selected_file_list[i])), '# result file name'))
                    if i == 0:
                        file_id.write( '\n')
                        file_id.write( '{0}\n'.format('# If there are more files, you have to repeat the section file-1 with the data of each file.'))
                        file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(result_transfer_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def download_result_dataset(cluster_name, log, function=None):
    '''
    Download the result dataset of a run from the cluster.
    '''
    
    # initialize the control variable
    OK = True

    # get the read transfer config file
    result_transfer_config_file = get_result_transfer_config_file()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # get and check the result transfer config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the result transfer config file ...\n')
    if check_result_transfer_config_file(strict=True):
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The result transfer config file is not valid.\n')
        log.write('Please correct this file or recreate the config files.\n')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name)
        for error in error_list:
            log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        sftp_client = xssh.create_sftp_client(ssh_transport)

    # get the option dictionary
    if OK:
        result_transfer_options_dict = xlib.get_option_dict(result_transfer_config_file)

    # download the result dataset
    if OK:

        # get the sections list
        sections_list = []
        for section in result_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # get the experiment identification, run identification and local directory from the section "identification"
        experiment_id = result_transfer_options_dict['identification']['experiment_id']
        result_dataset_id = result_transfer_options_dict['identification']['result_dataset_id']
        status = result_transfer_options_dict['identification']['status'].lower()
        local_dir = result_transfer_options_dict['identification']['local_dir']

        # download files when the status is uncompressed
        if status == 'uncompressed':

            # for each section "file-n"
            for section in sections_list:

                # check than the section identification is like file-n 
                if re.match('^file-[0-9]+$', section):

                    # get the dataset subdirectory and file name
                    dataset_subdirectory = result_transfer_options_dict[section]['dataset_subdirectory']
                    file_name = result_transfer_options_dict[section]['file_name']

                    # check if the dataset subdirectory is created
                    pathlib.Path(os.path.normpath('{0}/{1}'.format(local_dir, dataset_subdirectory))).mkdir(parents=True, exist_ok=True)

                    # assign the cluster path and local path
                    cluster_path = '{0}/{1}/{2}/{3}/{4}'.format(xlib.get_cluster_result_dir(), experiment_id, result_dataset_id, dataset_subdirectory, file_name)
                    local_path = os.path.normpath('{0}/{1}/{2}'.format(local_dir, dataset_subdirectory, file_name))

                    # download the result file from the cluster
                    log.write(f'{xlib.get_separator()}\n')
                    log.write('Downloading the file {0} to {1} ...\n'.format(cluster_path, local_dir))
                    (OK, error_list) = xssh.get_file(sftp_client, cluster_path, local_path)
                    if OK:
                        log.write('The file has been downloaded.\n')
                    else:
                        for error in error_list:
                            log.write(f'{error}\n')
                        break

        # download files when the status is compressed
        elif status == 'compressed':

            # assign the cluster path and local path
            cluster_path = '{0}/{1}/{2}'.format(xlib.get_cluster_result_dir(), experiment_id, result_dataset_id)
            local_path = '{0}/{1}'.format(local_dir, result_dataset_id)

            # download the result file from the cluster
            log.write(f'{xlib.get_separator()}\n')
            log.write('Downloading the file {0} to {1} ...\n'.format(cluster_path, local_dir))
            (OK, error_list) = xssh.get_file(sftp_client, cluster_path, local_path)
            if OK:
                log.write('The file has been downloaded.\n')
            else:
                for error in error_list:
                    log.write(f'{error}\n')

    # close the SSH transport connection
    if OK:
        xssh.close_ssh_transport_connection(ssh_transport)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

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

def check_result_transfer_config_file(strict):
    '''
    Check the result transfer config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the result transfer config file
    result_transfer_config_file = get_result_transfer_config_file()


    # get the option dictionary
    try:
        result_transfer_options_dict = xlib.get_option_dict(result_transfer_config_file)
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in result_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = result_transfer_options_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False
            elif not experiment_id.isidentifier():
                error_list.append('*** ERROR: the key "experiment_id" value in the section "identification" has some non-alphanumeric characters')
                OK = False

            # check section "identification" - key "result_dataset_id"

            result_dataset_id = result_transfer_options_dict.get('identification', {}).get('result_dataset_id', not_found)
            if result_dataset_id == not_found:
                error_list.append('*** ERROR: the key "result_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "status"
            status = result_transfer_options_dict.get('identification', {}).get('status', not_found).lower()
            if status == not_found:
                error_list.append('*** ERROR: the key "status" is not found in the section "identification".')
                OK = False
            else:
                if status not in ['compressed', 'uncompressed']:
                    error_list.append('*** ERROR: the key "status" value in the section "identification" has to be uncompressed or compressed.')
                    OK = False
                    status == 'WRONG'

            # check section "identification" - key "local_dir"
            local_dir = result_transfer_options_dict.get('identification', {}).get('local_dir', not_found)
            if local_dir == not_found:
                error_list.append('*** ERROR: the key "local_dir" is not found in the section "identification".')
                OK = False
            elif not os.path.isdir(local_dir):
                error_list.append('*** ERROR: the key "local_id" value in the section "identification" is a non existing directory path.')
                OK = False

        # check section "file-1"
        if status == 'uncompressed':
            if 'file-1' not in sections_list:
                error_list.append('*** ERROR: the section "file-1" is not found.')
                OK = False

        # check all sections "file-n"
        if status == 'uncompressed':
            for section in sections_list:

                if section not in ['identification']:

                    # check than the section identification is like file-n 
                    if not re.match('^file-[0-9]+$', section):
                        error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                        OK = False

                    else:

                        # check section "file-n" - key "dataset_subdirectory"
                        dataset_subdirectory = result_transfer_options_dict.get(section, {}).get('dataset_subdirectory', not_found)
                        if dataset_subdirectory == not_found:
                            error_list.append('*** ERROR: the key "dataset_subdirectory" is not found in the section "{0}".'.format(section))
                            OK = False
                        elif not xlib.is_valid_path(dataset_subdirectory, 'linux'):
                            error_list.append('*** ERROR: the file {0} in the key "dataset_subdirectory" of the section "{1}" has a non valid file name.'.format(dataset_subdirectory, section))
                            OK = False

                        # check section "file-n" - key "file_name"
                        file_name = result_transfer_options_dict.get(section, {}).get('file_name', not_found)
                        if file_name == not_found:
                            error_list.append('*** ERROR: the key "file_name" is not found in the section "{0}".'.format(section))
                            OK = False
                        elif not xlib.is_valid_path(file_name, 'linux'):
                            error_list.append('*** ERROR: the file {0} in the key "file_name" of the section "{1}" has a non valid file name.'.format(file_name, section))
                            OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe result transfer config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_result_dataset_dict(cluster_name, experiment_id, status, passed_connection, ssh_client):
    '''
    Get a dictionary with the result datasets of an experiment in the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the result directory in the cluster
    cluster_result_dir = xlib.get_cluster_result_dir()

    # initialize the dictionary of the result datasets
    result_dataset_dict = {}

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the result directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(cluster_result_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the result directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the cluster {1}.\n'.format(cluster_result_dir, cluster_name))
            OK = False

    # get the dictionary of the result datasets
    if OK:
        if status == 'uncompressed':
            command = 'cd  {0}/{1}; for list in `ls`; do ls -ld $list | grep -v ^- > /dev/null && echo $list; done;'.format(cluster_result_dir, experiment_id)
        elif status == 'compressed':
            command = 'cd {0}/{1}; for list in `ls`; do ls -ld $list | grep -v ^d > /dev/null && echo $list; done;'.format(cluster_result_dir, experiment_id)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            if status == 'uncompressed':
                input_pattern = '{0}-(.+)-(.+)'
                output_pattern = '{0} ({1} {2})'
            elif status == 'compressed':
                input_pattern = '{0}-(.+)-(.+).tar.gz'
                output_pattern = '{0} ({1} {2}) [compressed]'
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    result_dataset_id = line
                    if result_dataset_id.startswith(xlib.get_bowtie2_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_bowtie2_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_bowtie2_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_busco_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_busco_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_busco_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_cd_hit_est_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_cd_hit_est_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_cd_hit_est_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_cutadapt_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_cutadapt_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_cutadapt_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_cuffdiff_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_cuffdiff_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_cuffdiff_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_cufflinks_cuffmerge_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_cufflinks_cuffmerge_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_cufflinks_cuffmerge_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_cuffquant_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_cuffquant_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_cuffquant_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_ddradseq_simulation_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_ddradseq_simulation_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_ddradseq_simulation_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_fastqc_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_fastqc_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_fastqc_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_ggtrinity_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_ggtrinity_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_ggtrinity_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_gmap_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_gmap_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_gmap_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_gsnap_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_gsnap_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_gsnap_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_gzip_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_gzip_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_gzip_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_hisat2_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_hisat2_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_hisat2_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_htseq_count_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_htseq_count_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_htseq_count_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_insilico_read_normalization_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_insilico_read_normalization_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_insilico_read_normalization_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_ipyrad_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_ipyrad_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_ipyrad_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_kallisto_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_kallisto_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_kallisto_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_quast_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_quast_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_quast_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_ref_eval_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_ref_eval_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_ref_eval_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_rnaquast_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_rnaquast_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_rnaquast_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_rsem_eval_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_rsem_eval_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_rsem_eval_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_rsitesearch_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_rsitesearch_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_rsitesearch_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_soapdenovo2_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_soapdenovo2_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_soapdenovo2_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_soapdenovotrans_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_soapdenovotrans_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_soapdenovotrans_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_star_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_star_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_star_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_starcode_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_starcode_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_starcode_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_toa_process_pipeline_aminoacid_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_toa_process_pipeline_aminoacid_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_toa_process_pipeline_aminoacid_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_toa_process_pipeline_nucleotide_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_toa_process_pipeline_nucleotide_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_toa_process_pipeline_nucleotide_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_tophat_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_tophat_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_tophat_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_transabyss_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_transabyss_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_transabyss_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_transcript_filter_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_transcript_filter_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_transcript_filter_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_transcriptome_blastx_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_transcriptome_blastx_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_transcriptome_blastx_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_transrate_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_transrate_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_transrate_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_trimmomatic_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_trimmomatic_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_trimmomatic_name(), date, time)
                    elif result_dataset_id.startswith(xlib.get_trinity_code()+'-'):
                        mo = re.match(input_pattern.format(xlib.get_trinity_code()), result_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        result_dataset_name = output_pattern.format(xlib.get_trinity_name(), date, time)
                    else:
                        result_dataset_name = result_dataset_id
                    result_dataset_dict[result_dataset_id] = {'result_dataset_id': result_dataset_id, 'result_dataset_name': result_dataset_name}

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and dictionary of the result datasets
    return (OK, error_list, result_dataset_dict)

#-------------------------------------------------------------------------------

def get_result_dataset_name_list(cluster_name, experiment_id, status, app_list, passed_connection=False, ssh_client=None):
    '''
    Get a list of the result dataset names of an experiment in the cluster.
    '''

    # initialize the list of the result dataset names
    result_dataset_name_list = []

    # get the dictionary of the result datasets
    (OK, error_list, result_dataset_dict) = get_result_dataset_dict(cluster_name, experiment_id, status, passed_connection, ssh_client)

    # build the list of the result dataset names
    for result_dataset_id in result_dataset_dict.keys():
        for app in app_list:
            if app == xlib.get_all_applications_selected_code() or result_dataset_id.startswith(f'{app}-'):
                result_dataset_name_list.append(result_dataset_dict[result_dataset_id]['result_dataset_name'])
                break

    # sort the list of the result dataset names
    if result_dataset_name_list != []:
        result_dataset_name_list.sort()

    # return the control variable, error list and list of the result dataset names
    return (OK, error_list, result_dataset_name_list)

#-------------------------------------------------------------------------------

def get_result_dataset_id(cluster_name, experiment_id, result_dataset_name, status, passed_connection=False, ssh_client=None):
    '''
    Get the result dataset identification from the result dataset name.
    '''

    # initialize the control variable
    result_dataset_id_found = None

    # get the dictionary of the result datasets
    (OK, error_list, result_dataset_dict) = get_result_dataset_dict(cluster_name, experiment_id, status, passed_connection, ssh_client)

    # search the result dataset identification
    if OK:
        for result_dataset_id in result_dataset_dict.keys():
            if result_dataset_dict[result_dataset_id]['result_dataset_name'] == result_dataset_name:
                result_dataset_id_found = result_dataset_id
                break

    # return the control variable, error list and result dataset identification
    return (OK, error_list, result_dataset_id_found)

#-------------------------------------------------------------------------------

def get_result_transfer_config_file():
    '''
    Get the transfer config file path of the results of a run.
    '''

    # assign the result transfer config file
    result_transfer_config_file = '{0}/{1}'.format(xlib.get_config_dir(), 'result-transfer-config.txt')

    # return the result transfer config file
    return result_transfer_config_file

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the result datasets used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
