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
import urllib

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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).\n')
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
        log.write(f'Adding node {node_name} in cluster {cluster_name} using StarCluster ...\n')
        log.write('\n')
        command = f'{xlib.get_starcluster()} addnode {cluster_name} --alias={node_name}'
        rc = xlib.run_command(command, log)
        log.write('\n')
        if rc == 0:
            log.write('The node is added.\n')
        else:
            log.write(f'*** ERROR: Return code {rc} in command -> {command}\n')
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
            log.write(f'*** ERROR: The cluster {cluster_name} is not running. Its state is {master_state_code} ({master_state_name}).')
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
        log.write(f'Removing node {node_name} in cluster {cluster_name} using StarCluster ...\n')
        log.write('\n')
        command = f'{xlib.get_starcluster()} removenode --confirm {cluster_name} --alias={node_name}'
        rc = xlib.run_command(command, log)
        log.write('\n')
        if rc == 0:
            log.write('The node is removed.\n')
        else:
            log.write(f'*** ERROR: Return code {rc} in command -> {command}\n')
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
    node_script_path = f'./{os.path.basename(local_script_path)}'

    # set the infrastructure software installation log path in node
    node_log_path = node_script_path[:node_script_path.find('.sh')] + '.log'

    # build the infrastructure software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Building the infrastructure software installation script {local_script_path} ...\n')
        (OK, error_list) = build_infrastructure_software_installation_script(cluster_name)
        if OK:
            log.write('The file is built.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH client connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Connecting the SSH client to node {node_name} ...\n')
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name, node_name)
        if OK:
            log.write('The SSH client is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SSH transport connection
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Connecting the SSH transport to node {node_name} ...\n')
        (OK, error_list, ssh_transport) = xssh.create_ssh_transport_connection(cluster_name, node_name)
        if OK:
            log.write('The SSH transport is connected.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # create the SFTP client 
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Connecting the SFTP client to node {node_name} ...\n')
        sftp_client = xssh.create_sftp_client(ssh_transport)
        log.write('The SFTP client is connected.\n')

    # download the AWSCLI2 compressed file to local computer
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Downloading the {xlib.get_awscli_name()} compressed file to local computer ...\n')
        local_path = f'{xlib.get_temp_dir()}/{xlib.get_awscli_name()}.zip'
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))
        try:
            urllib.request.urlretrieve(xlib.get_awscli_url(), local_path)
        except Exception as e:
            log.write(f'*** EXCEPTION: "{e}".')
            log.write(f'*** ERROR: The file {xlib.get_awscli_url()} can not be downloaded.\n')
            OK = False
        else:
            log.write('The file is downloaded.\n')

    # upload the AWSCLI2 compressed file to cluster
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the {xlib.get_awscli_name()} compressed file to the cluster ...\n')
        cluster_path = f'./{os.path.basename(local_path)}'
        (OK, error_list) = xssh.put_file(sftp_client, local_path, cluster_path)
        if OK:
            log.write('The file is uploaded.\n')
        else:
            for error in error_list:
                log.write(f'{error}\n')

    # upload the infraestructe software installation script to the node
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Uploading the infraestructe software installation script to the node {node_name} ...\n')
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
        command = f'chmod u+x {node_script_path}'
        (OK, stdout, stderr) = xssh.execute_cluster_command(ssh_client, command)
        if OK:
            log.write('The run permision is set.\n')
        else:
            log.write(f'*** ERROR: Wrong command ---> {command}\n')

    # submit the infraestructe software installation script
    if OK:
        log.write(f'{xlib.get_separator()}\n')
        log.write(f'Submitting the infraestructe software installation script in the node {node_name} ...\n')
        command = f'{node_script_path} &>{node_log_path} &'
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

def build_infrastructure_software_installation_script(cluster_name):
    '''
    Build the infrastructure software installation script.
    '''

    # initialize the control variable and the error list
    OK = True
    error_list = []

    # get the connetion data
    (user_id, access_key_id, secret_access_key)  = xconfiguration.get_basic_aws_data()

    # get the old region and user identification
    current_region_name = xconfiguration.get_current_region_name()

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the dataset structure and NGScloud_volume
    dataset_structure = ngscloud_options_dict['dataset info']['dataset_structure']

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
            script_file_id.write( 'export AWS_CONFIG_FILE=/home/ubuntu/.aws/config\n')
            script_file_id.write( 'export AWS_SHARED_CREDENTIALS_FILE=/home/ubuntu/.aws/credentials\n')
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
            if dataset_structure in [xconfiguration.get_dataset_structure_singlevolume(), xconfiguration.get_dataset_structure_none()]:
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
            script_file_id.write( 'function install_awscli\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the AWS CLI ..."\n')
            script_file_id.write(f'    unzip {xlib.get_awscli_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then unzip $RC; fi\n')
            script_file_id.write( '    sudo ./aws/install\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then install $RC; fi\n')
            script_file_id.write( '    rm -rf aws\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then rm $RC; fi\n')
            script_file_id.write(f'    rm {xlib.get_awscli_name()}.zip\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then rm $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function setup_aws\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Setting up AWS ..."\n')
            script_file_id.write( '    UBUNTU_AWS_DIR=/home/ubuntu/.aws\n')
            script_file_id.write( '    mkdir --parents $UBUNTU_AWS_DIR\n')
            script_file_id.write(f'    CONFIG_FILE=$UBUNTU_AWS_DIR/config\n')
            script_file_id.write( '    echo "[default]" > $CONFIG_FILE\n')
            script_file_id.write(f'    echo "region = {current_region_name}" >> $CONFIG_FILE\n')
            script_file_id.write( '    CREDENTIALS_FILE=$UBUNTU_AWS_DIR/credentials\n')
            script_file_id.write( '    echo "[default]" > $CREDENTIALS_FILE\n')
            script_file_id.write(f'    echo "aws_access_key_id = {access_key_id}" >> $CREDENTIALS_FILE\n')
            script_file_id.write(f'    echo "aws_secret_access_key = {secret_access_key}" >> $CREDENTIALS_FILE\n')
            script_file_id.write( '    sudo echo "AWS is set up."\n')
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
            script_file_id.write( 'function install_xorg\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package xorg ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install xorg\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
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
            script_file_id.write( 'function install_libxt6\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package libxt6 ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install libxt6\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_parallel\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package parallel ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install parallel\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function install_texlive\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Installing the package texlive ..."\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-latex-base\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-fonts-recommended\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-fonts-extra\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    sudo apt-get --assume-yes install texlive-latex-extra\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "The package is installed."\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            script_file_id.write( 'function uninstall_mysql\n')
            script_file_id.write( '{\n')
            script_file_id.write( '    echo "$SEP"\n')
            script_file_id.write( '    echo "Uninstalling MySQL ..."\n')
            script_file_id.write( '    sudo apt-get purge --auto-remove --assume-yes mysql-client mysql-client-5.5 mysql-client-core-5.5 mysql-common mysql-server mysql-server-5.5 mysql-server-core-5.5\n')
            script_file_id.write( '    RC=$?\n')
            script_file_id.write( '    if [ $RC -ne 0 ]; then manage_error apt-get $RC; fi\n')
            script_file_id.write( '    echo "MySQL is uninstalled."\n')
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
            script_file_id.write( '    send_mail ok\n')
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
            script_file_id.write( '    exit 3\n')
            script_file_id.write( '}\n')
            script_file_id.write( '#-------------------------------------------------------------------------------\n')
            process_name = 'Infrastructure software installation'
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
            if dataset_structure in [xconfiguration.get_dataset_structure_singlevolume(), xconfiguration.get_dataset_structure_none()]:
                script_file_id.write( 'create_dataset_structure\n')
            script_file_id.write( 'install_awscli\n')
            script_file_id.write( 'setup_aws\n')
            script_file_id.write( 'fix_source_list\n')
            script_file_id.write( 'install_xorg\n')
            script_file_id.write( 'install_libtbb2\n')
            script_file_id.write( 'install_libxt6\n')
            script_file_id.write( 'install_parallel\n')
            script_file_id.write( 'install_texlive\n')
            script_file_id.write( 'uninstall_mysql\n')
            # -- script_file_id.write( 'create_swapfile\n')
            script_file_id.write( 'end\n')
    except Exception as e:
        error_list.append(f'*** EXCEPTION: "{e}".')
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
    infrastructure_software_installation_script = f'{xlib.get_temp_dir()}/infrastructure_software_installation.sh'

    # return the infrastructure software installation script path
    return infrastructure_software_installation_script

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the functions related to the node operation used in both console mode and gui mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
