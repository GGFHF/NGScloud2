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
This file contains functions related to read datasets used in both console mode
and gui mode.
'''

#-------------------------------------------------------------------------------

import os
import re
import subprocess
import sys

import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def create_read_transfer_config_file(experiment_id='exp001', local_dir='./data', selected_file_list=['rnaseq-1.fastq']):
    '''
    Create o recreate the read transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the read transfer config file path
    read_transfer_config_file = get_read_transfer_config_file()

    # create the read transfer config file to upload RNA-seq read files
    try:
        if not os.path.exists(os.path.dirname(read_transfer_config_file)):
            os.makedirs(os.path.dirname(read_transfer_config_file))
        with open(read_transfer_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current transfer.'))
            file_id.write( '{0}\n'.format('# The files will be copied in the cluster directory {0}/experiment_id/{1}.'.format(xlib.get_cluster_read_dir(), xlib.get_uploaded_read_dataset_name())))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the experiment.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('experiment_id = {0}'.format(experiment_id), '# experiment identification'))
            for i in range(len(selected_file_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first read file.'))
                file_id.write( '{0}\n'.format('[file-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('local_path = {0}/{1}'.format(local_dir, selected_file_list[i]), '# local path of a read file'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more files, you have to repeat the section file-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(read_transfer_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def upload_read_dataset(cluster_name, log, function=None):
    '''
    Upload the read dataset to the cluster.
    '''
    
    # initialize the control variable
    OK = True

    # get the read transfer config file
    read_transfer_config_file = get_read_transfer_config_file()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # get and check the read transfer config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the read transfer config file ...\n')
    if check_read_transfer_config_file(strict=True):
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The read transfer config file is not valid.\n')
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
        read_transfer_options_dict = xlib.get_option_dict(read_transfer_config_file)

    # get the experiment identification and create the experiment reads directory
    if OK:

        # get the experiment identification
        experiment_id = read_transfer_options_dict['identification']['experiment_id']

        # Get the directory of read and results datasets of the experiment
        cluster_experiment_reads_dir = xlib.get_cluster_experiment_read_dataset_dir(experiment_id, xlib.get_uploaded_read_dataset_name())
        cluster_experiment_result_dir = xlib.get_cluster_experiment_result_dir(experiment_id)

        # create the experiment reads directory
        log.write(f'{xlib.get_separator()}\n')
        log.write('The reads directory {0} in the cluster is being created ...\n'.format(cluster_experiment_reads_dir))
        command = 'mkdir --parents {0}'.format(cluster_experiment_reads_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory is created.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

        # create the experiment run result directory
        log.write(f'{xlib.get_separator()}\n')
        log.write('The run result directory {0} in the cluster is being created ...\n'.format(cluster_experiment_result_dir))
        command = 'mkdir --parents {0}'.format(cluster_experiment_result_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory is created.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # upload the read dataset
    if OK:

        # get the sections list
        sections_list = []
        for section in read_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # for each section "file-n"
        for section in sections_list:

            # check than the section identification is like file-n 
            if re.match('^file-[0-9]+$', section):

                # get local path and cluster directory
                local_path = read_transfer_options_dict[section]['local_path']

                # upload the reference file to the cluster
                log.write(f'{xlib.get_separator()}\n')
                log.write('The file {0} is being uploaded to {1} ...\n'.format(local_path, cluster_experiment_reads_dir))
                cluster_path = '{0}/{1}'.format(cluster_experiment_reads_dir, os.path.basename(local_path))
                (OK, error_list) = xssh.put_file(sftp_client, local_path, cluster_path)
                if OK:
                    log.write('The file has been uploaded.\n')
                else:
                    for error in error_list:
                        log.write(f'{error}\n')
                    break

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

def check_read_transfer_config_file(strict):
    '''
    Check the read transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the read transfer config file
    read_transfer_config_file = get_read_transfer_config_file()

    # get the option dictionary
    read_transfer_options_dict = xlib.get_option_dict(read_transfer_config_file)
    try:
        read_transfer_options_dict = xlib.get_option_dict(read_transfer_config_file)
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in read_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "experiment_id"
            experiment_id = read_transfer_options_dict.get('identification', {}).get('experiment_id', not_found)
            if experiment_id == not_found:
                error_list.append('*** ERROR: the key "experiment_id" is not found in the section "identification".')
                OK = False

        # check section "file-1"
        if 'file-1' not in sections_list:
            error_list.append('*** ERROR: the section "file-1" is not found.')
            OK = False

        # check all sections "file-n"
        for section in sections_list:

            if section not in ['identification']:

                # check than the section identification is like file-n 
                if not re.match('^file-[0-9]+$', section):
                    error_list.append('*** ERROR: the section "{0}" has a wrong identification.'.format(section))
                    OK = False

                else:

                    # check section "file-n" - key "local_path"
                    local_path = read_transfer_options_dict.get(section, {}).get('local_path', not_found)
                    if local_path == not_found:
                        error_list.append('*** ERROR: the key "local_path" is not found in the section "{0}".'.format(section))
                        OK = False
                    else:
                        try:
                            open(local_path, mode='r').close()
                        except FileNotFoundError:
                            if strict:
                                error_list.append('*** ERROR: the file {0} in the key "local_path" of the section "{1}" does not exist or it is not accessible.'.format(local_path, section))
                                OK = False
                            else:
                                error_list.append('*** WARNING: the file {0} in the key "local_path" of the section "{1}" does not exist or it is not accessible.'.format(local_path, section))
                        except OSError:
                            error_list.append('*** ERROR: the file name "{0}" in the key "local_path" of the section "{1}" is not correct.'.format(local_path, section))
                            OK = False


    # warn that the reads config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe read transfer config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_read_dataset_dict(cluster_name, experiment_id, passed_connection=False, ssh_client=None):
    '''
    Get a dictoinary with the read datasets of an experiment in the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the read directory in the cluster
    cluster_read_dir = xlib.get_cluster_read_dir()

    # initialize the dictionary of the read datasets
    read_dataset_dict = {}

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the read directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(cluster_read_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the read directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(cluster_read_dir, cluster_name))
            OK = False

    # get the dictionary of the read datasets
    if OK:
        command = 'ls {0}/{1}'.format(cluster_read_dir, experiment_id)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    read_dataset_id = line
                    if read_dataset_id == xlib.get_uploaded_read_dataset_name():
                        read_dataset_name = 'uploaded reads'
                    elif read_dataset_id.startswith(xlib.get_cutadapt_code()):
                        mo = re.match('{0}-(.+)-(.+)'.format(xlib.get_cutadapt_code()), read_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        read_dataset_name = '{0} ({1} {2})'.format(xlib.get_cutadapt_name(), date, time)
                    elif read_dataset_id.startswith(xlib.get_insilico_read_normalization_code()):
                        mo = re.match('{0}-(.+)-(.+)'.format(xlib.get_insilico_read_normalization_code()), read_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        read_dataset_name = '{0} ({1} {2})'.format(xlib.get_insilico_read_normalization_name(), date, time)
                    elif read_dataset_id.startswith(xlib.get_trimmomatic_code()):
                        mo = re.match('{0}-(.+)-(.+)'.format(xlib.get_trimmomatic_code()), read_dataset_id)
                        date = mo.group(1)
                        time = mo.group(2)
                        read_dataset_name = '{0} ({1} {2})'.format(xlib.get_trimmomatic_name(), date, time)
                    else:
                        read_dataset_name = 'xxx'
                    read_dataset_dict[read_dataset_id] = {'read_dataset_id': read_dataset_id, 'read_dataset_name': read_dataset_name}

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and dictionary of the read datasets
    return (OK, error_list, read_dataset_dict)

#-------------------------------------------------------------------------------

def get_read_dataset_name_list(cluster_name, experiment_id, passed_connection=False, ssh_client=None):
    '''
    Get a list of the read dataset names of an experiment in the cluster.
    '''

    # initialize the list of the read dataset names
    read_dataset_name_list = []

    # get the dictionary of the read datasets
    (OK, error_list, read_dataset_dict) = get_read_dataset_dict(cluster_name, experiment_id, passed_connection, ssh_client)

    # build the list of the read dataset names
    for read_dataset_id in read_dataset_dict.keys():
        read_dataset_name_list.append(read_dataset_dict[read_dataset_id]['read_dataset_name'])

    # sort the list of the read dataset names
    if read_dataset_name_list != []:
        if 'uploaded reads' in read_dataset_name_list:
            read_dataset_name_list.remove('uploaded reads')
            if read_dataset_name_list != []:
                read_dataset_name_list = ['uploaded reads'] + sorted(read_dataset_name_list)
            else:
                read_dataset_name_list = ['uploaded reads']
        else:
            read_dataset_name_list.sort()

    # return the control variable, error list and list of the read dataset names
    return (OK, error_list, read_dataset_name_list)

#-------------------------------------------------------------------------------

def get_read_dataset_id(cluster_name, experiment_id, read_dataset_name, passed_connection=False, ssh_client=None):
    '''
    Get the read dataset identification from the read dataset name.
    '''

    # initialize the read dataset identification
    read_dataset_id_found = None

    # get the dictionary of the read datasets
    (OK, error_list, read_dataset_dict) = get_read_dataset_dict(cluster_name, experiment_id, passed_connection, ssh_client)

    # search the read dataset identification
    if OK:
        for read_dataset_id in read_dataset_dict.keys():
            if read_dataset_dict[read_dataset_id]['read_dataset_name'] == read_dataset_name:
                read_dataset_id_found = read_dataset_id
                break

    # return the control variable, error list and read dataset identification
    return (OK, error_list, read_dataset_id_found)

#-------------------------------------------------------------------------------

def get_read_transfer_config_file():
    '''
    Get the read transfer config file path.
    '''

    # assign the read transfer config file
    read_transfer_config_file = '{0}/{1}'.format(xlib.get_config_dir(), 'read-transfer-config.txt')

    # return the read transfer config file
    return read_transfer_config_file

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the read datasets used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
