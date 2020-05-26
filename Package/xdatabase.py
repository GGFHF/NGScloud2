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
This file contains functions related to databases used in both console mode
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

def create_database_transfer_config_file(local_dir='./data', selected_file_list=['RefSeq_Plan_Protein.faa', 'RefSeq_Plan_Protein.pal'], database_dataset_id = 'RefSeq_Plan_Protein'):
    '''
    Create o recreate the database transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # create the config config file to upload database files
    try:
        if not os.path.exists(os.path.dirname(get_database_transfer_config_file())):
            os.makedirs(os.path.dirname(get_database_transfer_config_file()))
        with open(get_database_transfer_config_file(), mode='w', encoding='iso-8859-1', newline='\n') as file_id:
            file_id.write( '{0}\n'.format('# You must review the information of this file and update the values with the corresponding ones to the current transfer.'))
            file_id.write( '{0}\n'.format('# The database files have to be located in the cluster directory {0}/database_dataset_id'.format(xlib.get_cluster_database_dir())))
            file_id.write( '{0}\n'.format('# The database_dataset_id is fixed in the identification section.'))
            file_id.write( '\n')
            file_id.write( '{0}\n'.format('# This section has the information identifies the database dataset.'))
            file_id.write( '[identification]\n')
            file_id.write( '{0:<50} {1}\n'.format('database_dataset_id = {0}'.format(database_dataset_id), '# database dataset identification'))
            file_id.write( '{0:<50} {1}\n'.format('local_dir = {0}'.format(local_dir), '# local directory of database files'))
            for i in range(len(selected_file_list)):
                file_id.write( '\n')
                if i == 0:
                    file_id.write( '{0}\n'.format('# This section has the information of the first database file.'))
                file_id.write( '{0}\n'.format('[file-{0}]'.format(i + 1)))
                file_id.write( '{0:<50} {1}\n'.format('file_name = {0}'.format(selected_file_list[i]), '# database file name'))
                if i == 0:
                    file_id.write( '\n')
                    file_id.write( '{0}\n'.format('# If there are more files, you have to repeat the section file-1 with the data of each file.'))
                    file_id.write( '# The section identification has to be library-n (n is an integer not repeated)\n')
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be created'.format(get_database_transfer_config_file()))
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def upload_database_dataset(cluster_name, log, function=None):
    '''
    Upload the database dataset to the cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the log window does not have to be closed
    if not isinstance(log, xlib.DevStdOut):
        log.write('This process might take several minutes. Do not close this window, please wait!\n')

    # check the database transfer config file
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking the database transfer config file ...\n')
    if check_database_transfer_config_file(strict=True):
        log.write('The file is OK.\n')
    else:
        log.write('*** ERROR: The database transfer config file is not valid.\n')
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

    # upload the database dataset
    if OK:

        # get the option dictionary
        database_transfer_options_dict = xlib.get_option_dict(get_database_transfer_config_file())

        # get the database dataset identification and the local directory of the database files
        database_dataset_id = database_transfer_options_dict['identification']['database_dataset_id']
        local_dir = database_transfer_options_dict['identification']['local_dir']

        # set the cluster database directory
        cluster_database_dir = '{0}/{1}'.format(xlib.get_cluster_database_dir(), database_dataset_id)

        # create the data directory in the cluster
        log.write(f'{xlib.get_separator()}\n')
        log.write('The database directory {0} in the cluster is being created ...\n'.format(cluster_database_dir))
        command = 'mkdir --parents {0}'.format(cluster_database_dir)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The directory is created.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

        # get the sections list
        sections_list = []
        for section in database_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # for each section "file-n"
        for section in sections_list:

            # check than the section identification is like file-n 
            if re.match('^file-[0-9]+$', section):

                # get the file name
                file_name = database_transfer_options_dict[section]['file_name']

                # set the local path and cluster path
                local_path = '{0}/{1}'.format(local_dir, file_name)
                cluster_path = '{0}/{1}'.format(cluster_database_dir, file_name)

                # upload the database file to the cluster
                log.write(f'{xlib.get_separator()}\n')
                log.write('The file {0} is being uploaded to {1} ...\n'.format(file_name, cluster_database_dir))
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

def check_database_transfer_config_file(strict):
    '''
    Check the database transfer config file.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # intitialize variable used when value is not found
    not_found = '***NOTFOUND***'.upper()

    # get the option dictionary
    try:
        database_transfer_options_dict = xlib.get_option_dict(get_database_transfer_config_file())
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
        error_list.append('*** ERROR: The option dictionary could not be built from the config file')
        OK = False
    else:

        # get the sections list
        sections_list = []
        for section in database_transfer_options_dict.keys():
            sections_list.append(section)
        sections_list.sort()

        # check section "identification"
        if 'identification' not in sections_list:
            error_list.append('*** ERROR: the section "identification" is not found.')
            OK = False
        else:

            # check section "identification" - key "database_dataset_id"
            database_dataset_id = database_transfer_options_dict.get('identification', {}).get('database_dataset_id', not_found)
            if database_dataset_id == not_found:
                error_list.append('*** ERROR: the key "database_dataset_id" is not found in the section "identification".')
                OK = False

            # check section "identification" - key "local_dir"
            local_dir = database_transfer_options_dict.get('identification', {}).get('local_dir', not_found)
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
                    error_list.append(f'*** ERROR: the section "{section}" has a wrong identification.')
                    OK = False

                else:

                    # check section "file-n" - key "file_name"
                    file_name = database_transfer_options_dict.get(section, {}).get('file_name', not_found)
                    if file_name == not_found:
                        error_list.append('*** ERROR: the key "file_name" is not found in the section "{0}".'.format(section))
                        OK = False
                    else:
                        if not os.path.isfile(os.path.join(local_dir, file_name)):
                            error_list.append('*** ERROR: the file {0} in the key "file_name" does not exist or is not accessible in the local directory {1}.'.format(file_name, local_dir))
                            OK = False

    # warn that the database config file is not valid if there are any errors
    if not OK:
        error_list.append('\nThe database transfer config file is not valid. Please, correct this file or recreate it.')

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_database_dataset_dict(cluster_name, passed_connection=False, ssh_client=None):
    '''
    Get a dictionary with the database datasets in the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of the database datasets
    database_dataset_dict = {}

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the database directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_database_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the database directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the cluster {1}.\n'.format(xlib.get_cluster_database_dir(), cluster_name))
            OK = False

    # build the dictionary of the database datasets
    if OK:
        command = 'ls {0}'.format(xlib.get_cluster_database_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    database_dataset_id = line
                    database_dataset_name = database_dataset_id
                    database_dataset_dict[database_dataset_id] = {'database_dataset_id': database_dataset_id, 'database_dataset_name': database_dataset_name}

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and dictionary of the database datasets
    return (OK, error_list, database_dataset_dict)

#-------------------------------------------------------------------------------

def get_database_dataset_name_list(cluster_name, passed_connection=False, ssh_client=None):
    '''
    Get a list of the database dataset names in the cluster.
    '''

    # initialize the list of the database dataset names
    database_dataset_name_list = []

    # get the dictionary of the database datasets
    (OK, error_list, database_dataset_dict) = get_database_dataset_dict(cluster_name, passed_connection, ssh_client)

    # build the list of the database dataset names
    for database_dataset_id in database_dataset_dict.keys():
        database_dataset_name_list.append(database_dataset_dict[database_dataset_id]['database_dataset_name'])

    # sort the list of the database dataset names
    if database_dataset_name_list != []:
        database_dataset_name_list.sort()

    # return the control variable, error list and list of the database dataset names
    return (OK, error_list, database_dataset_name_list)

#-------------------------------------------------------------------------------

def get_database_dataset_id(cluster_name, database_dataset_name, passed_connection=False, ssh_client=None):
    '''
    Get the database dataset identification from the database dataset name.
    '''

    # initialize the control variable
    database_dataset_id_found = None

    # get the dictionary of the database datasets
    (OK, error_list, database_dataset_dict) = get_database_dataset_dict(cluster_name, passed_connection, ssh_client)

    # search the database dataset identification
    if OK:
        for database_dataset_id in database_dataset_dict.keys():
            if database_dataset_dict[database_dataset_id]['database_dataset_name'] == database_dataset_name:
                database_dataset_id_found = database_dataset_id
                break

    # return the control variable, error list and database dataset identification
    return (OK, error_list, database_dataset_id_found)

#-------------------------------------------------------------------------------

def get_database_file_name_list(cluster_name, database_dataset_id, file_type, passed_connection=False, ssh_client=None):
    '''
    Get a list of the database file names in a database dataset of the cluster.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the dictionary of the database datasets
    database_file_name_list = []

    # create the SSH client connection
    if not passed_connection:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)

    # check the app directory is created
    if OK:
        command = '[ -d {0} ] && echo RC=0 || echo RC=1'.format(xlib.get_cluster_database_dir())
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if stdout[len(stdout) - 1] != 'RC=0':
            error_list.append('*** ERROR: There is not any volume mounted in the database directory.\n')
            error_list.append('You have to link a volume in the mounting point {0} for the cluster {1}.\n'.format(xlib.get_cluster_database_dir(), cluster_name))
            OK = False

    # build the list of the database file name of the database dataset
    if OK:
        if file_type == 'all':
            command = 'ls {0}'.format(xlib.get_cluster_database_dataset_dir(database_dataset_id))
        else:
            command = 'ls {0}/*{1}'.format(xlib.get_cluster_database_dataset_dir(database_dataset_id), file_type)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            for line in stdout:
                line = line.rstrip('\n')
                if line != 'lost+found':
                    database_file_name_list.append(os.path.basename(line))

    # close the SSH client connection
    if OK and not passed_connection:
        xssh.close_ssh_client_connection(ssh_client)

    # return the control variable, error list and list of the database file names
    return (OK, error_list, database_file_name_list)

#-------------------------------------------------------------------------------

def get_database_transfer_config_file():
    '''
    Get the database transfer config file path.
    '''

    # assign the database transfer config file
    database_transfer_config_file = '{0}/{1}'.format(xlib.get_config_dir(), 'database-transfer-config.txt')

    # return the database transfer config file
    return database_transfer_config_file

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains functions related to the databases used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
