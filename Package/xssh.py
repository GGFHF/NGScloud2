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
This file contains the functions related to the SSH used in both console mode and gui mode.
'''
#-------------------------------------------------------------------------------

import io
import os
import re

import paramiko

import xconfiguration
import xcluster
import xec2
import sys

#-------------------------------------------------------------------------------

def create_ssh_client_connection(cluster_name, node_name=None, user='root'):
    '''
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the SSH client object
    ssh_client = None

    # set the node name when it is None
    if node_name is None:
        if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native():
            node_name = 'instance'
        elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
            node_name = 'master'

    # initialize the port
    port = 22

    # get the keypair file
    keypair_file = xconfiguration.get_keypair_file()

    # get the public dns nameof the master node in the cluster
    public_dns_name = xec2.get_node_public_dns_name(cluster_name, node_name)

    # create the SSH client object
    ssh_client = paramiko.SSHClient()

    # get the RSA key text from the corresponding file
    try:
        file = open(keypair_file,'r')
        records = file.read()
        records_inmemory = io.StringIO(records)
        rsakey = paramiko.RSAKey.from_private_key(records_inmemory)
    except Exception as e:
        error_list.append('*** ERROR: The file {0} can not be read.'.format(keypair_file))
        OK = False

    # accept auto-accept unknown keys
    if OK:
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    
    # start the connection
    if OK:
        try:
            ssh_client.connect(public_dns_name, port, user, pkey=rsakey)
        except Exception as e:
            error_list.append('*** ERROR: {0} can not be connected.'.format(public_dns_name))
            OK = False

    # return the control variable, the error list and the SSH client object
    return (OK, error_list, ssh_client)

#-------------------------------------------------------------------------------

def execute_cluster_command(ssh_client, command):
    '''
    '''

    # initialize the control variable
    OK = True

    # execute the command in the ssh client
    (stdin, stdout, stderr) = ssh_client.exec_command(command)

    # flush the stdout internal buffer
    stdout.flush()

    # initialize the string lines list corresponding to stdout
    stdout_string_lines_list = []

    # build a string lines list corresponding to the stdout
    for stdout_bytes_line in stdout.read().splitlines():
        # non-ASCII caracters are replaces by one blank space
        stdout_bytes_line = re.sub(b'[^\x00-\x7F]+', b' ', stdout_bytes_line)
        # create a string from the bytes literal
        stdout_string_line = stdout_bytes_line.decode('utf-8')
        # add the string line to the string lines list
        stdout_string_lines_list.append(stdout_string_line)

    # flush the stderr internal buffer
    stderr.flush()

    # initialize the string lines corresponding to the stderr
    stderr_string_lines_list = []

    # build a string lines list corresponding to the stderr and set False to OK variable if there are any lines
    stderr_bytes_lines_list = stderr.read().splitlines()
    if stderr_bytes_lines_list != []:
        OK = False
        for stderr_bytes_line in stderr_bytes_lines_list:
            # non-ASCII caracters are replaces by one blank space
            stderr_bytes_line = re.sub(b'[^\x00-\x7F]+', b' ', stderr_bytes_line)
            # create a string from the bytes literal
            stderr_string_line = stderr_bytes_line.decode('utf-8')
            # add the string line to the string lines list
            stderr_string_lines_list.append(stderr_string_line)

    # return the control variableThis file contains the functions related to the SSH used in both console mode and gui mode.
    return (OK, stdout_string_lines_list, stderr_string_lines_list)

#-------------------------------------------------------------------------------

def close_ssh_client_connection(ssh_client):
    '''
    '''

    # close the connection of the SSH client object
    ssh_client.close()

#-------------------------------------------------------------------------------

def create_ssh_transport_connection(cluster_name, node_name=None, user='root'):
    '''
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # initialize the SSH transport objet
    ssh_transport = None

    # set the node name when it is None
    if node_name is None:
        if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native():
            node_name = 'instance'
        elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
            node_name = 'master'

    # initialize the port
    port = 22

    # get the keypair file
    keypair_file = xconfiguration.get_keypair_file()

    # get the public dns nameof the master node in the cluster
    public_dns_name = xec2.get_node_public_dns_name(cluster_name, node_name)

    # create the SSH transport object
    try:
        ssh_transport = paramiko.Transport((public_dns_name, port))
    except Exception as e:
        error_list.append('*** ERROR: {0} can not be connected.'.format(public_dns_name))
        OK = False

    # get the RSA key text from the corresponding file
    if OK:
        try:
            file = open(keypair_file,'r')
            records = file.read()
            records_inmemory = io.StringIO(records)
            rsakey = paramiko.RSAKey.from_private_key(records_inmemory)
        except Exception as e:
            error_list.append('*** ERROR: The file {0} can not be read.'.format(keypair_file))
            OK = False
    
    # start the connection
    if OK:
        try:
            ssh_transport.connect(username=user, pkey=rsakey)
        except Exception as e:
            error_list.append('*** ERROR: User/pkey is not valid in {0}.'.format(public_dns_name))
            OK = False

    # return the control variable, the error list and the SSH transport objet
    return (OK, error_list, ssh_transport)

#-------------------------------------------------------------------------------

def create_sftp_client(ssh_transport):
    '''
    '''

    # create the SFTP client object
    sftp_client = paramiko.SFTPClient.from_transport(ssh_transport)

    # return the SFTP client object
    return sftp_client

#-------------------------------------------------------------------------------

def put_file(sftp_client, local_path, cluster_path):
    '''
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # upload the local file to the cluster
    #try:
    sftp_client.put(local_path, cluster_path)
    #except Exception as e:
    #   error_list.append('*** ERROR: It is not possible to upload the local file {0} to cluster file {1}'.format(local_path, cluster_path))
    #   OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_file(sftp_client, cluster_path, local_path):
    '''
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # download the cluster file to the local machine
    try:
        sftp_client.get(cluster_path, local_path)
    except Exception as e:
       error_list.append('*** ERROR: It is not possible to download the cluster file {0} to local file {1}'.format(cluster_path, local_path))
       OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def close_ssh_transport_connection(ssh_transport):
    '''
    '''

    # close the connection of the SSH client object
    ssh_transport.close()

#-------------------------------------------------------------------------------

def submit_script(cluster_name, ssh_client, current_run_dir, script, log):
    '''
    '''

    # initialize the control variable
    OK = True

    # submit the script starter
    if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native():
        command = f'nohup {current_run_dir}/{script} &>/dev/null &'
    elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
        sge_env = xcluster.get_sge_env()
        command = f'{sge_env}; qsub -V -b n -cwd {current_run_dir}/{script}'
    (OK, stdout, stderr) = execute_cluster_command(ssh_client, command)
    if OK:
        log.write('The script is submitted.\n')
        for line in stdout:
            log.write('{0}\n'.format(line))
    else:
        log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the functions related to the SSH used in both console mode and gui mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
