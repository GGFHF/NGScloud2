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
This file contains the functions related to the node operation used in both console
mode and gui mode.
'''
#-------------------------------------------------------------------------------

import os
import subprocess
import sys

import xconfiguration
import xec2
import xlib
import xssh

#-------------------------------------------------------------------------------

def add_node(cluster_name, node_name, log, function=None):
    '''
    Add a node in a cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the requirements are being verified 
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking process requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name, node_name='master')
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).\n'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the node is not running
    if OK:
        pass

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # add node
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Adding node {0} in cluster {1} using StarCluster ...\n'.format(node_name, cluster_name))
        log.write('\n')
        command = '{0} addnode {1} --alias={2}'.format(xlib.get_starcluster(), cluster_name, node_name)
        rc = xlib.run_command(command, log)
        log.write('\n')
        if rc == 0:
            log.write('The node is added.\n')
        else:
            log.write('*** ERROR: Return code {0} in command -> {1}\n'.format(rc, command))
            OK = False

    # install infrastructure software in the node
    if OK:
        OK = install_node_infrastructure_software(cluster_name, node_name, log)

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

def remove_node(cluster_name, node_name, log, function=None):
    '''
    Remove a node in a cluster.
    '''

    # initialize the control variable
    OK = True

    # warn that the requirements are being verified 
    log.write(f'{xlib.get_separator()}\n')
    log.write('Checking process requirements ...\n')

    # check the master is running
    if OK:
        (master_state_code, master_state_name) = xec2.get_node_state(cluster_name, node_name='master')
        if master_state_code != 16:
            log.write('*** ERROR: The cluster {0} is not running. Its state is {1} ({2}).'.format(cluster_name, master_state_code, master_state_name))
            OK = False

    # check the node is running
    if OK:
        pass

    # warn that the requirements are OK 
    if OK:
        log.write('Process requirements are OK.\n')

    # remove node
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Removing node {0} in cluster {1} using StarCluster ...\n'.format(node_name, cluster_name))
        log.write('\n')
        command = '{0} removenode --confirm {1} --alias={2}'.format(xlib.get_starcluster(), cluster_name, node_name)
        rc = xlib.run_command(command, log)
        log.write('\n')
        if rc == 0:
            log.write('The node is removed.\n')
        else:
            log.write('*** ERROR: Return code {0} in command -> {1}\n'.format(rc, command))
            OK = False

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

def install_node_infrastructure_software(cluster_name, node_name, log):
    '''
    Install infrastructure software in a node.
    '''

    # initialize the control variable
    OK = True

    # get the infrastructure software installation script path in local compute
    local_script_path = get_infrastructure_software_installation_script()

    # set the infrastructure software installation script path in node
    node_script_path = './{0}'.format(os.path.basename(local_script_path))

    # set the infrastructure software installation log path in node
    node_log_path = node_script_path[:node_script_path.find('.sh')] + '.log'

    # build the infrastructure software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Building the infrastructure software installation script {0} ...\n'.format(local_script_path))
        (OK, error_list) = build_infrastructure_software_installation_script(cluster_name)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH client to node {0} ...\n'.format(node_name))
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name, node_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SSH transport to node {0} ...\n'.format(node_name))
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name, node_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Connecting the SFTP client to node {0} ...\n'.format(node_name))
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # upload the infraestructe software installation script to the node
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Uploading the infraestructe software installation script to the node {0} ...\n'.format(node_name))
        (OK, error_list) = xssh.put_file(sftp_client, local_script_path, node_script_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # set run permision to the infraestructe software installation script in the cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Setting on the run permision ...\n')
        command = 'chmod u+x {0}'.format(node_script_path)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the infraestructe software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Submitting the infraestructe software installation script in the node {0} ...\n'.format(node_name))
        command = '{0} &>{1} &'.format(node_script_path, node_log_path)
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The script is submitted.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

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

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def install_node_infrastructure_software_ant(cluster_name, node_name, log):
    '''
    Install infrastructure software in a node.
    '''

    # initialize the control variable
    OK = True

    # create the SSH client connection
    log.write(f'{xlib.get_separator()}\n')
    log.write('Connecting the SSH client in node {0} ...\n'.format(node_name))
    (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name, node_name)
    if OK:
        log.write('The SSH client is connected.\n')
    else:
        for error in error_list:
            log.write(f'{error}\n')

    # update file /etc/apt/sources.list
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Updating the file /etc/apt/sources.list in node {0} ...\n'.format(node_name))
        command = 'sed -i "s/us-east-1\.ec2\.archive\.ubuntu\.com/old-releases\.ubuntu\.com/g" /etc/apt/sources.list'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        log.write('command: {0}\n'.format(command))
        log.write('OK: {0}\n'.format(OK))
        log.write('stdout: {0}\n'.format(stdout))
        log.write('stderr: {0}\n'.format(stderr))
        OK = True
    if OK:
        command = 'sed -i "s/security\.ubuntu\.com/old-releases\.ubuntu\.com/g" /etc/apt/sources.list'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        log.write('command: {0}\n'.format(command))
        log.write('OK: {0}\n'.format(OK))
        log.write('stdout: {0}\n'.format(stdout))
        log.write('stderr: {0}\n'.format(stderr))
        OK = True
    if OK:
        command = 'export DEBIAN_FRONTEND=noninteractive; apt-get update'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        log.write('command: {0}\n'.format(command))
        log.write('OK: {0}\n'.format(OK))
        log.write('stdout: {0}\n'.format(stdout))
        log.write('stderr: {0}\n'.format(stderr))
        OK = True
    if OK:
        log.write('The file is updated.\n')
    else:
        log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # install libtbb2
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Installing the package libtbb2 in node {0} (it can be a very slow process, please be patient) ...\n'.format(node_name))
        command = 'export DEBIAN_FRONTEND=noninteractive; apt-get --assume-yes --force-yes install libtbb2'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        log.write('command: {0}\n'.format(command))
        log.write('OK: {0}\n'.format(OK))
        log.write('stdout: {0}\n'.format(stdout))
        log.write('stderr: {0}\n'.format(stderr))
        OK = True
        if OK:
            log.write('The package is installed.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    ## install OpenJDK8
    #if OK:
    #    log.write(f'{xlib.get_separator()}\n')
    #    log.write('Installing the package OpenJDK8 in node {0} (it can be a very slow process, please be patient) ...\n'.format(node_name))
    #    command = 'echo "deb http://ppa.launchpad.net/openjdk-r/ppa/ubuntu trusty main" | tee /etc/apt/sources.list.d/openjdk-r.list'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    command = 'echo "deb-src http://ppa.launchpad.net/openjdk-r/ppa/ubuntu trusty main" | tee -a /etc/apt/sources.list.d/openjdk-r.list'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    #command = 'export DEBIAN_FRONTEND=noninteractive; apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EEA14886'
    #    command = 'export DEBIAN_FRONTEND=noninteractive; apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EB9B1D8886F44E2A'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    command = 'export DEBIAN_FRONTEND=noninteractive; apt-get update'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    command = 'export DEBIAN_FRONTEND=noninteractive; apt-get --assume-yes --force-yes install openjdk-8-jdk'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    command = 'export DEBIAN_FRONTEND=noninteractive; update-java-alternatives --set java-1.8.0-openjdk-amd64'
    #    (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    #    log.write('command: {0}\n'.format(command))
    #    log.write('OK: {0}\n'.format(OK))
    #    log.write('stdout: {0}\n'.format(stdout))
    #    log.write('stderr: {0}\n'.format(stderr))
    #    OK = True
    #if OK:
    #    log.write('OpenJDK8 is installed.\n')
    #else:
    #    log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # install Oracle Java 8
    # -- if OK:
    # --     log.write(f'{xlib.get_separator()}\n')
    # --     log.write('Installing Oracle Java 8 in node {0} (it can be a very slow process, please be patient) ...\n'.format(node_name))
    # --     command = 'echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee /etc/apt/sources.list.d/webupd8team-java.list'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     command = 'echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     command = 'export DEBIAN_FRONTEND=noninteractive; apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EEA14886'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     command = 'export DEBIAN_FRONTEND=noninteractive; apt-get update'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     command = 'echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | debconf-set-selections'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     command = 'export DEBIAN_FRONTEND=noninteractive; apt-get --assume-yes --force-yes install openjdk-8-jdk'
    # --     (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
    # --     log.write('command: {0}\n'.format(command))
    # --     log.write('OK: {0}\n'.format(OK))
    # --     log.write('stdout: {0}\n'.format(stdout))
    # --     log.write('stderr: {0}\n'.format(stderr))
    # --     OK = True
    # -- if OK:
    # --     log.write('Java 8 is installed.\n')
    # -- else:
    # --     log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # close the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write('Closing the SSH client connection ...\n')
        xssh.close_ssh_client_connection(ssh_client)
        log.write('The connection is closed.\n')

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def build_infrastructure_software_installation_script(cluster_name):
    '''
    Build the infrastructure software installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the dataset structure and NGScloud_volume
    dataset_structure = ngscloud_options_dict['dataset info']['dataset_structure']
    ngscloud_volume = ngscloud_options_dict['dataset info']['ngscloud_volume']

    # write the infrastructure software installation script
    try:
        if not os.path.exists(os.path.dirname(get_infrastructure_software_installation_script())):
            os.makedirs(os.path.dirname(get_infrastructure_software_installation_script()))
        with open(get_infrastructure_software_installation_script(), mode='w', encoding='iso-8859-1', newline='\n') as script_file_id:
            script_file_id.write( '#!/bin/bash\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'SEP="#########################################"\n')
            script_file_id.write( 'export HOST_IP=`curl --silent checkip.amazonaws.com`\n')
            script_file_id.write( 'export HOST_ADDRESS="ec2-${HOST_IP//./-}-compute-1.amazonaws.com"\n')
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
            script_file_id.write( 'function create_dataset_structure\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Creating the dataset structure ..."\n')
            script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_app_dir()}\n')
            script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_database_dir()}\n')
            script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_read_dir()}\n')
            script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_reference_dir()}\n')
            script_file_id.write(f'    sudo mkdir --parents {xlib.get_cluster_result_dir()}\n')
            script_file_id.write( '    echo "The dataset structure is created."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function fix_source_list\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Fixing file /etc/apt/sources.list ..."\n')
            script_file_id.write( '    sed -i "s/us-east-1.ec2.archive.ubuntu.com/old-releases.ubuntu.com/g" /etc/apt/sources.list\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '    sed -i "s/security.ubuntu.com/old-releases.ubuntu\.com/g" /etc/apt/sources.list\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error sed $RC; fi\n')
            script_file_id.write( '    apt-get update\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "The file is fixed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_libtbb2\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package libtbb2 ..."\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    apt-get --assume-yes install libtbb2\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_mailutils\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package mailutils ..."\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    HOST_IP=`curl checkip.amazonaws.com`\n')
            script_file_id.write( '    HOST_IP2=`echo "${HOST_IP//./-}"`\n')
            script_file_id.write( '    HOST_ADDRESS="ec2-${HOST_IP2}-compute-1.amazonaws.com"\n')
            script_file_id.write( '    echo "HOST_IP: $HOST_IP   HOST_ADDRESS: $HOST_ADDRESS"\n')
            script_file_id.write( '    debconf-set-selections <<< "postfix postfix/mailname string $HOST_ADDRESS"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    debconf-set-selections <<< "postfix postfix/main_mailer_type string \'Internet Site\'"\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error debconf-set-selections $RC; fi\n')
            script_file_id.write( '    apt-get --assume-yes install mailutils\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function create_swapfile\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Creating a file which will be used for swap ..."\n')
            script_file_id.write( '    sudo dd if=/dev/zero of=/swapfile bs=1024 count=2097152\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error dd $RC; fi\n')
            script_file_id.write( '    sudo chmod 600 /swapfile\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error chmod $RC; fi\n')
            script_file_id.write( '    sudo mkswap /swapfile\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error mkswap $RC; fi\n')
            script_file_id.write( '    sudo swapon /swapfile\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error swapon $RC; fi\n')
            script_file_id.write( '    sudo echo "/swapfile swap swap defaults 0 0" >> /etc/fstab\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error echo $RC; fi\n')
            script_file_id.write( '    echo\n')
            script_file_id.write( '    echo "The file is created."\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: Infrastructure Software Installation"\n')
            message = xlib.get_mail_message_ok('Infrastructure Software Installation', cluster_name)
            script_file_id.write(f'    MESSAGE="{message}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject "$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
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
            script_file_id.write(f'    RECIPIENT={xconfiguration.get_contact_data()}\n')
            script_file_id.write(f'    SUBJECT="{xlib.get_project_name()}: Infrastructure software installation"\n')
            message = xlib.get_mail_message_wrong('Infrastructure Software Installation', cluster_name)
            script_file_id.write(f'    MESSAGE="{message}"\n')
            script_file_id.write( '    mail --append "Content-type: text/html;" --append "FROM:root@NGScloud2" --subject="$SUBJECT" "$RECIPIENT" <<< "$MESSAGE"\n')
            script_file_id.write( '    exit 3\n')
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
            if dataset_structure == xconfiguration.get_dataset_structure_singlevolume() and ngscloud_volume != '':
                script_file_id.write( 'create_dataset_structure\n')
            script_file_id.write( 'fix_source_list\n')
            script_file_id.write( 'install_libtbb2\n')
            script_file_id.write( 'install_mailutils\n')
            script_file_id.write( 'create_swapfile\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** ERROR: The file {get_infrastructure_software_installation_script()} can not be created')
        OK = False

    # return the control variable and the error list
    return (OK, error_list)

#-------------------------------------------------------------------------------

def get_infrastructure_software_installation_script():
    '''
    Get the infrastructure software installation script path in the local computer.
    '''

    # assign infrastructure software installation script path
    infrastructure_software_installation_script = '{0}/{1}'.format(xlib.get_temp_dir(), 'infrastructure_software_installation.sh')

    # return the infrastructure software installation script path
    return infrastructure_software_installation_script

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the functions related to the node operation used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
