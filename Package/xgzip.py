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
This file contains functions related to the gzip process used in both console
mode and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import sys

import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_gzip_config_file(action='compress', dataset_type='read', experiment_id='exp001', dataset_id=xlib.get_uploaded_read_dataset_name(), file_list=['default']):
    '''
    Create gzip config file with the default options. It is necessary
    update the options in each run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # adapt the values by default to the dataset type when the gzip config file is created
    if file_list[0] == 'default':
        if dataset_type == 'reference':
            experiment_id = None
            dataset_id = 'Athaliana'
            file_list = ['./GCF_000001735.3_TAIR10_genomic.fna']
        elif dataset_type == 'database':
            experiment_id = None
            dataset_id = 'RefSeq_Plan_Protein'
            file_list = ['./RefSeq_Plan_Protein.faa', './RefSeq_Plan_Protein.pal']
        elif dataset_type == 'read':
            file_list = ['./rnaseq-a_1.fastq', './rnaseq-a_2.fastq']
        elif dataset_type == 'result':
            dataset_id = 'fastqc-170224-134133'
            file_list = ['./fixed_reads_end1_fastqc.html', './fixed_reads_end2_fastqc.html']
        elif dataset_type in ['whole-result']:
            dataset_id = 'fastqc-170224-134133'
            file_list = [None]

    # create the gzip config file and write the default options
    try:
        if not os.path.exists(os.path.dirname(get_gzip_config_file(dataset_type))):
            os.makedirs(os.path.dirname(get_gzip_config_file(dataset_type)))
        with open(get_gzip_config_file(dataset_type), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# This section has the information identifies the dataset.'))
            file_id.write( '[identification]\n')
            value = 'NONE' if dataset_type == 'reference' else experiment_id
            comment = 'It has to be always NONE' if dataset_type == 'reference' else 'experiment identification'
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(value), '# {0}'.format(comment)))
            file_id.write( '{0:<50} {1}\n'.format('dataset_type = {0}'.format(dataset_type), '# dataset type (it has to be always {0})'.format(dataset_type)))
            file_id.write( '{0:<50} {1}\n'.format('dataset_id = {0}'.format(dataset_id), '# dataset identification'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information to set the gzip parameters'))
            file_id.write( '{0}\n'.format('[gzip parameters]'))
            file_id.write( '{0:<50} {1}\n'.format('action = {0}'.format(action), '# action: compress or decompress'))
            if dataset_type in ['reference', 'database', 'read', 'result']:
                for i in range(len(file_list)):
                    file_id.write( '\n')
                    if i == 0:
                        file_id.write( '{0}\n'.format('# This section has the information of the first file.'))
                    file_id.write( '{0}\n'.format('[file-{0}]'.format(i + 1)))
                    file_id.write( '{0:<50} {1}\n'.format('dataset_subdirectory = {0}'.format(os.path.dirname(file_list[i])), '# subdirectory of {0} dataset'.format(dataset_type)))
                    file_id.write( '{0:<50} {1}\n'.format('file_name = {0}'.format(os.path.basename(file_list[i])), '# {0} file name'.format(dataset_type)))
                    if i == 0:
                        file_id.write( '\n')
                        file_id.write( '{0}\n'.format('# If there are more files, you have to repeat the section file-1 with the data of each file.'))
                        file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be recreated'.format(get_gzip_config_file(dataset_type)))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def run_gzip_process(cluster_name, dataset_type, log, function=None):
    '''
    Run a gzip process.
    '''

    # initialize the control variable
    OK = True

    # get the gzip code and name
    gzip_code = xlib.get_gzip_code()
    gzip_name = xlib.get_gzip_name()

    # get the gzip option dictionary
    gzip_option_dict = xlib.get_option_dict(get_gzip_config_file(dataset_type))

    # get the experiment identification
    experiment_id = gzip_option_dict['identification']['experiment_id']

    # get the gzip process script path in the local computer
    gzip_process_script = get_gzip_process_script(dataset_type)

    # get the gzip process starter path in the local computer
    gzip_process_starter = get_gzip_process_starter(dataset_type)

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the gzip config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the {0} config file ...\n'.format(gzip_name))
    (OK, error_list) = check_gzip_config_file(dataset_type, strict=True)
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

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # determine the run directory in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Determining the run directory in the cluster ...\n')
        if dataset_type == 'reference':
            current_run_dir = xlib.get_cluster_current_run_dir('reference', gzip_code)
        elif dataset_type == 'database':
            current_run_dir = xlib.get_cluster_current_run_dir('database', gzip_code)
        else:
            current_run_dir = xlib.get_cluster_current_run_dir(experiment_id, gzip_code)
        command = f'mkdir --parents {current_run_dir}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write(f'The directory path is {current_run_dir}.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the gzip process script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process script {0} ...\n'.format(gzip_process_script))
        (OK, error_list) = build_gzip_process_script(cluster_name, dataset_type, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            log.write('*** ERROR: The file could not be built.\n')

    # upload the gzip process script to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process script {0} to the directory {1} ...\n'.format(gzip_process_script, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(gzip_process_script))
        (OK, error_list) = xssh.put_file(sftp_client, gzip_process_script, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the gzip process script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(gzip_process_script)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(gzip_process_script))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # build the gzip process starter
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the process starter {0} ...\n'.format(gzip_process_starter))
        (OK, error_list) = build_gzip_process_starter(dataset_type, current_run_dir)
        if OK:
            log.write('The file is built.\n')
        else:
            log.write('***ERROR: The file could not be built.\n')

    # upload the gzip process starter to the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the process starter {0} to the directory {1} ...\n'.format(gzip_process_starter, current_run_dir))
        cluster_path = '{0}/{1}'.format(current_run_dir, os.path.basename(gzip_process_starter))
        (OK, error_list) = xssh.put_file(sftp_client, gzip_process_starter, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the gzip process starter in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision of {0}/{1} ...\n'.format(current_run_dir, os.path.basename(gzip_process_starter)))
        command = 'chmod u+x {0}/{1}'.format(current_run_dir, os.path.basename(gzip_process_starter))
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the gzip process
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the process script {0}/{1} ...\n'.format(current_run_dir, os.path.basename(gzip_process_starter)))
        OK = xssh.submit_script(cluster_name, ssh_client, current_run_dir, os.path.basename(gzip_process_starter), log)

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

def check_gzip_config_file(dataset_type, strict):
    '''
    Check the gzip config file of a run.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        gzip_option_dict = xlib.get_option_dict(get_gzip_config_file(dataset_type))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in gzip_option_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = gzip_option_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False
            elif dataset_type == 'reference' and experiment_id.upper() != 'NONE':
                error_list.append('*** ERROR: the key "experiment_id" has to be always NONE')
                OK = False

            # check section "identification" - key "dataset_type"
            dataset_type_2 = gzip_option_dict.get('identification', {}).get('dataset_type', not_found)
            if dataset_type_2 == not_found:
                error_list.append('*** ERROR: the key "dataset_type" is not found in the section "identification".')
                OK = False
            else:
                if dataset_type in ['reference', 'read']:
                    if dataset_type_2.lower() != dataset_type:
                        error_list.append('*** ERROR: the key "dataset_type" has to be {0}.'.format(dataset_type))
                        OK = False
                elif dataset_type == 'result':
                    if dataset_type_2.lower() not in ['result', 'whole-result']:
                        error_list.append('*** ERROR: the key "dataset_type" has to be result or whole-result.')
                        OK = False

            # check section "identification" - key "dataset_id"
            dataset_id = gzip_option_dict.get('identification', {}).get('dataset_id', not_found)
            if dataset_id == not_found:
                error_list.append('*** ERROR: the key "dataset_id" is not found in the section "identification".')
                OK = False

        # check section "gzip parameters"
        if 'gzip parameters' not in sections_list:
            error_list.append('*** ERROR: the section "gzip parameters" is not found.')
            OK = False
        else:

            # check section "gzip parameters" - key "action"
            action = gzip_option_dict.get('gzip parameters', {}).get('action', not_found)
            if action == not_found:
                error_list.append('*** ERROR: the key "action" is not found in the section "gzip parameters".')
                OK = False
            else:
                if action.lower() not in ['compress', 'decompress']:
                    error_list.append('*** ERROR: the key "action" has to be compress or decompress.')
                    OK = False

        # check section "file-1"
        if dataset_type_2.lower() in ['reference', 'database', 'read', 'result']:
            if 'file-1' not in sections_list:
                error_list.append('*** ERROR: the section "file-1" is not found.')
                OK = False

        # check all sections "file-n"
        if dataset_type_2.lower() in ['reference', 'database', 'read', 'result']:
            for section in sections_list:

                if section not in ['identification', 'gzip parameters']:

                    # check than the section identification is like file-n 
                    if not re.match('^file-[0-9]+$', section):
                        error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                        OK = False

                    else:

                        # check section "file-n" - key "dataset_subdirectory"
                        dataset_subdirectory = gzip_option_dict.get(section, {}).get('dataset_subdirectory', not_found)
                        if dataset_subdirectory == not_found:
                            error_list.append('*** ERROR: the key "dataset_subdirectory" is not found in the section "{0}".'.format(section))
                            OK = False
                        elif not xlib.is_valid_path(dataset_subdirectory, 'linux'):
                            error_list.append('*** ERROR: the file {0} in the key "dataset_subdirectory" of the section "{1}" has a non valid file name.'.format(dataset_subdirectory, section))
                            OK = False

                        # check section "file-n" - key "file_name"
                        file_name = gzip_option_dict.get(section, {}).get('file_name', not_found)
                        if file_name == not_found:
                            error_list.append('*** ERROR: the key "file_name" is not found in the section "{0}".'.format(section))
                            OK = False
                        elif not xlib.is_valid_path(file_name, 'linux'):
                            error_list.append('*** ERROR: the file {0} in the key "file_name" of the section "{1}" has a non valid file name.'.format(file_name, section))
                            OK = False

    # warn that the results config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe {0} config file is not valid. Please, correct this file or recreate it.'.format(xlib.get_gzip_name()))

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gzip_process_script(cluster_name, dataset_type, current_run_dir):
    '''
    Build the current gzip process script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the gzip option dictionary
    gzip_option_dict = xlib.get_option_dict(get_gzip_config_file(dataset_type))

    # get the options
    experiment_id = gzip_option_dict['identification']['experiment_id']
    dataset_type_2 = gzip_option_dict['identification']['dataset_type']
    dataset_id = gzip_option_dict['identification']['dataset_id']
    action = gzip_option_dict['gzip parameters']['action']

    # get the sections list
    sections_list = []
    for section in gzip_option_dict.keys():
        sections_list.append(section)
    sections_list.sort()

    # build the dataset subdirectory and file name lists
    dataset_subdirectory_list = []
    file_name_list = []
    for section in sections_list:
        # if the section identification is like library-n
        if re.match('^file-[0-9]+$', section):
            dataset_subdirectory = gzip_option_dict[section]['dataset_subdirectory']
            dataset_subdirectory_list.append(dataset_subdirectory)
            file_name = gzip_option_dict[section]['file_name']
            file_name_list.append(file_name)

    # get the dataset directory
    if dataset_type_2 == 'reference':
        dataset_dir = xlib.get_cluster_reference_dataset_dir(dataset_id)
    elif dataset_type_2 == 'database':
        dataset_dir = xlib.get_cluster_database_dataset_dir(dataset_id)
    elif dataset_type_2 == 'read':
        dataset_dir = xlib.get_cluster_experiment_read_dataset_dir(experiment_id, dataset_id)
    elif dataset_type_2 == 'result':
        dataset_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, dataset_id)
    elif dataset_type_2 == 'whole-result':
        dataset_dir = xlib.get_cluster_experiment_result_dataset_dir(experiment_id, dataset_id)

    # write the gzip process script
    try:
        if not os.path.exists(os.path.dirname(get_gzip_process_script(dataset_type_2))):
            os.makedirs(os.path.dirname(get_gzip_process_script(dataset_type_2)))
        with open(get_gzip_process_script(dataset_type_2), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
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
            script_file_id.write( '{0}\n'.format('function run_gzip_process'))
            script_file_id.write( '{\n')
            if dataset_type_2 in ['reference', 'database', 'read', 'result']:
                script_file_id.write(f'    cd {current_run_dir}\n')
                for i in range(len(dataset_subdirectory_list)):
                    script_file_id.write( '    echo "$SEP"\n')
                    script_file_id.write( '{0}\n'.format('    echo "Compressing/decompressing {0}/{1}/{2} ..."'.format(dataset_dir, dataset_subdirectory_list[i], file_name_list[i])))
                    script_file_id.write( '    /usr/bin/time \\\n')
                    script_file_id.write( '{0}\n'.format('        --format="Elapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                    if action == 'compress':
                        script_file_id.write( '{0}\n'.format('        gzip {0}/{1}/{2}'.format(dataset_dir, dataset_subdirectory_list[i], file_name_list[i])))
                    elif action == 'decompress':
                        script_file_id.write( '{0}\n'.format('        gzip --decompress {0}/{1}/{2}'.format(dataset_dir, dataset_subdirectory_list[i], file_name_list[i])))
                    script_file_id.write( '    RC=$?\n')
                    script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error gzip $RC; fi'))
            elif dataset_type_2 == 'whole-result':
                script_file_id.write(f'    cd {current_run_dir}\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Compressing/decompressing {0} ..."'.format(dataset_dir)))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '{0}\n'.format('        --format="Elapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                if action == 'compress':
                    script_file_id.write( '{0}\n'.format('        tar --create --gzip --verbose --file={0}.tar.gz {0}'.format(dataset_dir)))
                elif action == 'decompress':
                    script_file_id.write( '{0}\n'.format('        tar --extract --gzip --verbose --file={0} --directory=/'.format(dataset_dir)))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error tar $RC; fi\n')
                script_file_id.write( '    echo "$SEP"\n')
                script_file_id.write( '{0}\n'.format('    echo "Removing {0} ..."'.format(dataset_dir)))
                script_file_id.write( '    /usr/bin/time \\\n')
                script_file_id.write( '{0}\n'.format('        --format="Elapsed real time (s): %e\\nCPU time in kernel mode (s): %S\\nCPU time in user mode (s): %U\\nPercentage of CPU: %P\\nMaximum resident set size(Kb): %M\\nAverage total memory use (Kb):%K" \\'))
                script_file_id.write( '{0}\n'.format('        rm -rf {0}'.format(dataset_dir)))
                script_file_id.write( '    RC=$?\n')
                script_file_id.write( '{0}\n'.format('    if [ $RC -ne 0 ]; then manage_error rm $RC; fi'))
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
            process_name = f'{xlib.get_gzip_name()} process'
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
            script_file_id.write( '{0}\n'.format('run_gzip_process'))
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gzip_process_script(dataset_type_2)))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def build_gzip_process_starter(dataset_type, current_run_dir):
    '''
    Build the starter of the current gzip process.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the log file name
    log_file = xlib.get_cluster_log_file()

    # write the gzip process starter
    try:
        if not os.path.exists(os.path.dirname(get_gzip_process_starter(dataset_type))):
            os.makedirs(os.path.dirname(get_gzip_process_starter(dataset_type)))
        with open(get_gzip_process_starter(dataset_type), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '#!/bin/bash\n')
            file_id.write( '#-------------------------------------------------------------------------------\n')
            file_id.write( '{0}\n'.format('{0}/{1} &>>{0}/{2}'.format(current_run_dir, os.path.basename(get_gzip_process_script(dataset_type)), log_file)))
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_gzip_process_starter(dataset_type)))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_gzip_config_file(dataset_type):
    '''
    Get the gzip config file path.
    '''

    # assign the gzip config file path
    name = 'result' if dataset_type == 'whole-result' else dataset_type
    gzip_config_file = '{0}/{1}-{2}-config.txt'.format(xlib.get_config_dir(), xlib.get_gzip_code(), name)

    # return the gzip config file path
    return gzip_config_file

#-------------------------------------------------------------------------------

def get_gzip_process_script(dataset_type):
    '''
    Get the gzip process script path in the local computer.
    '''

    # assign the gzip script path
    name = 'result' if dataset_type == 'whole-result' else dataset_type
    gzip_process_script = '{0}/{1}-{2}-process.sh'.format(xlib.get_temp_dir(), xlib.get_gzip_code(), name)

    # return the gzip script path
    return gzip_process_script

#-------------------------------------------------------------------------------

def get_gzip_process_starter(dataset_type):
    '''
    Get the gzip process starter path in the local computer.
    '''

    # assign the gzip process starter path
    name = 'result' if dataset_type == 'whole-result' else dataset_type
    gzip_process_starter = '{0}/{1}-{2}-starter.sh'.format(xlib.get_temp_dir(), xlib.get_gzip_code(), name)

    # return the gzip starter path
    return gzip_process_starter

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the gzip process used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
