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
This file contains functions related to reference datasets used in both console mode
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

def create_reference_transfer_config_file(local_dir='./data', selected_file_list=['GCF_000146045.2_R64_genomic.fna.gz'], reference_dataset_id = 'Scerevisiae'):
    '''
    Create o recreate the reference transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the reference transfer config file path
    reference_transfer_config_file = get_reference_transfer_config_file()

    # create the reference config config file to upload reference files
    try:
        if not os.path.exists(os.path.dirname(reference_transfer_config_file)):
            os.makedirs(os.path.dirname(reference_transfer_config_file))
        with open(reference_transfer_config_file, mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current transfer.'))
            file_id.write( '{0}\n'.format('# The reference files have to be located in the cluster directory {0}/reference_dataset_id'.format(xlib.get_cluster_reference_dir())))
            file_id.write( '{0}\n'.format('# The reference_dataset_id is fixed in the identification section.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the reference dataset.'))
            file_id.write( '{0}\n'.format('[identification]'))
            file_id.write( '{0:<50} {1}\n'.format('reference_dataset_id = {0}'.format(reference_dataset_id), '# reference dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('local_dir = {0}'.format(local_dir), '# local directory of reference files'))
            for i in range(len(selected_file_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first reference file.'))
                file_id.write( '{0}\n'.format('[file-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('file_name = {0}'.format(selected_file_list[i]), '# reference file name'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more files, you have to repeat the section file-1 with the data of each file.'))
                    file_id.write( '{0}\n'.format('# The section identification has to be library-n (n is an integer not repeated)'))
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(reference_transfer_config_file))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def upload_reference_dataset(cluster_name, log, function=None):
    '''
    Upload the reference dataset to the cluster.
    '''

    # initialize the control variable
    OK = True

    # get the reference transfer config file
    reference_transfer_config_file = get_reference_transfer_config_file()

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the reference transfer config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the reference transfer config file ...\n')
    if check_reference_transfer_config_file(strict=True):
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The reference transfer config file is not valid.\n')
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

    # upload the reference dataset
    if OK:

        # get the option dictionary
        reference_transfer_options_dict = xlib.get_option_dict(get_reference_transfer_config_file())

        # get the reference dataset identification and the local directory of the reference files
        reference_dataset_id = reference_transfer_options_dict['identification']['reference_dataset_id']
        local_dir = reference_transfer_options_dict['identification']['local_dir']

        # set the cluster reference directory
        cluster_reference_dir = '{0}/{1}'.format(xlib.get_cluster_reference_dir(), reference_dataset_id)

        # create the data directory in the cluster
        log.write(f'{xlib.get_separator()}\n')
        log.write('The reference directory {0} in the cluster is being created ...\n'.format(cluster_reference_dir))
        command = 'mkdir --parents {0}'.format(cluster_reference_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory is created.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

        # get the sections list
        sections_list = []
        for section in reference_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # for each section "file-n"
        for section in sections_list:

            # check than the section identification is like file-n 
            if re.match('^file-[0-9]+$', section):

                # get the file name
                file_name = reference_transfer_options_dict[section]['file_name']

                # set the local path and cluster path
                local_path = '{0}/{1}'.format(local_dir, file_name)
                cluster_path = '{0}/{1}'.format(cluster_reference_dir, file_name)

                # upload the reference file in the cluster
                log.write(f'{xlib.get_separator()}\n')
                log.write('The file {0} is being uploaded to {1} ...\n'.format(file_name, cluster_reference_dir))
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

def check_reference_transfer_config_file(strict):
    '''
    Check the reference transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the reference transfer config file
    reference_transfer_config_file = get_reference_transfer_config_file()

    # get the option dictionary
    try:
        reference_transfer_options_dict = xlib.get_option_dict(reference_transfer_config_file)
    except Exception as e:
        error_list.append('*** ERROR: The syntax is WRONG.')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in reference_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "reference_dataset_id"
            reference_dataset_id = reference_transfer_options_dict.get('identification', {}).get('reference_dataset_id', not_found)
            if reference_dataset_id == not_found:
                error_list.append('*** ERROR: the key "reference_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "local_dir"
            local_dir = reference_transfer_options_dict.get('identification', {}).get('local_dir', not_found)
            if local_dir == not_found:
                error_list.append('*** ERROR: the key "local_dir" is not found in the section "identification".')
                OK = False
            else:
                if not os.path.isdir(local_dir):
                    error_list.append('*** ERROR: {0} is not a directory or does not exist.'.format(local_dir))
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

                    # check section "file-n" - key "file_name"
                    file_name = reference_transfer_options_dict.get(section, {}).get('file_name', not_found)
                    if file_name == not_found:
                        error_list.append('*** ERROR: the key "file_name" is not found in the section "{0}".'.format(section))
                        OK = False
                    else:
                        if not os.path.isfile(os.path.join(local_dir, file_name)):
                            error_list.append('*** ERROR: the file {0} in the key "file_name" does not exist or is not accessible in the local directory {1}.'.format(file_name, local_dir))
                            OK = False

    # warn that the reference config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe reference transfer config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_reference_dataset_dict(cluster_name, passed_connection=False, ssh_client=None):
    '''
    Get a dictionary with the reference datasets in the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the reference directory in the cluster
    cluster_reference_dir = xlib.get_cluster_reference_dir()

    # initialize the dictionary of the reference datasets
    reference_dataset_dict = {}

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the app directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(cluster_reference_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the reference directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(cluster_reference_dir, cluster_name))
            OK = False

    # build the dictionary of the reference datasets
    if OK:
        command = 'ls {0}'.format(cluster_reference_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    reference_dataset_id = line
                    reference_dataset_name = reference_dataset_id
                    reference_dataset_dict[reference_dataset_id] = {'reference_dataset_id': reference_dataset_id, 'reference_dataset_name': reference_dataset_name}

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and dictionary of the reference datasets
    return (OK, error_list, reference_dataset_dict)

#-------------------------------------------------------------------------------

def get_reference_dataset_name_list(cluster_name, passed_connection=False, ssh_client=None):
    '''
    Get a list of the reference dataset names in the cluster.
    '''

    # initialize the list of the reference dataset names
    reference_dataset_name_list = []

    # get the dictionary of the reference datasets
    (OK, error_list, reference_dataset_dict) = get_reference_dataset_dict(cluster_name, passed_connection, ssh_client)

    # build the list of the reference dataset names
    for reference_dataset_id in reference_dataset_dict.keys():
        reference_dataset_name_list.append(reference_dataset_dict[reference_dataset_id]['reference_dataset_name'])

    # sort the list of the reference dataset names
    if reference_dataset_name_list != []:
        reference_dataset_name_list.sort()

    # return the control variable, error list and list of the reference dataset names
    return (OK, error_list, reference_dataset_name_list)

#-------------------------------------------------------------------------------

def get_reference_dataset_id(cluster_name, reference_dataset_name, passed_connection=False, ssh_client=None):
    '''
    Get the reference dataset identification from the reference dataset name.
    '''

    # initialize the control variable
    reference_dataset_id_found = None

    # get the dictionary of the reference datasets
    (OK, error_list, reference_dataset_dict) = get_reference_dataset_dict(cluster_name, passed_connection, ssh_client)

    # search the reference dataset identification
    if OK:
        for reference_dataset_id in reference_dataset_dict.keys():
            if reference_dataset_dict[reference_dataset_id]['reference_dataset_name'] == reference_dataset_name:
                reference_dataset_id_found = reference_dataset_id
                break

    # return the control variable, error list and reference dataset identification
    return (OK, error_list, reference_dataset_id_found)

#-------------------------------------------------------------------------------

def get_reference_file_name_list(cluster_name, reference_dataset_id, passed_connection=False, ssh_client=None):
    '''
    Get a list of the reference file names in a reference dataset of the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the reference directory in the cluster
    cluster_reference_dir = xlib.get_cluster_reference_dir()

    # initialize the dictionary of the reference datasets
    reference_file_name_list = []

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the reference directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(cluster_reference_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the reference directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the template {1}.\n'.format(cluster_reference_dir, cluster_name))
            OK = False

    # get the reference dataset directory
    reference_dataset_dir = xlib.get_cluster_reference_dataset_dir(reference_dataset_id)

    # build the list of the reference file name of the reference dataset
    if OK:
        command = 'find {0} -maxdepth 1 -type f'.format(reference_dataset_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    reference_file_name_list.append(os.path.basename(line))

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and list of the reference file names
    return (OK, error_list, reference_file_name_list)

#-------------------------------------------------------------------------------

def get_reference_transfer_config_file():
    '''
    Get the reference transfer config file path.
    '''

    # assign the reference transfer config file
    reference_transfer_config_file = '{0}/{1}'.format(xlib.get_config_dir(), 'reference-transfer-config.txt')

    # return the reference transfer config file
    return reference_transfer_config_file

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the reference datasets used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
