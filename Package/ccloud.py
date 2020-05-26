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
This file contains the functions related to forms corresponding to Cloud Control
menu items in console mode.
'''

#-------------------------------------------------------------------------------

import os
import sys

import cinputs
import clib
import xbowtie2
import xbusco
import xcdhit
import xcluster
import xconfiguration
import xcufflinks
import xcutadapt
import xdatabase
import xddradseqtools
import xdetonate
import xec2
import xexpress
import xfastqc
import xgmap
import xgzip
import xhisat2
import xhtseq
import xinstance
import xipyrad
import xkallisto
import xlib
import xngshelper
import xnode
import xquast
import xraddesigner
import xread
import xreference
import xresult
import xrnaquast
import xsoapdenovo2
import xsoapdenovotrans
import xssh
import xstar
import xstarcode
import xtoa
import xtophat
import xtransabyss
import xtransrate
import xtrimmomatic
import xtrinity
import xvolume

#-------------------------------------------------------------------------------

def form_set_environment():
    '''
    Set the environment.
    '''

    # print headers
    clib.clear_screen()
    clib.print_headers_without_environment('Set environment')
    # -- print(f'function name: {sys._getframe().f_code.co_name}')

    # initialize the environment and the input environment
    xconfiguration.environment = ''
    environment = ''

    # get the current environments list
    environments_list = xconfiguration.get_environments_list()

    # print the available region names
    if environments_list != []:
        environments_list_text = str(environments_list).strip('[]').replace('\'', '')
        print(f'Current environments list: {environments_list_text} ...')
        input_text = '... Enter the environment name: '
    else:
        print('Currently there is not any environment recorded.')
        input_text = 'Enter a new environment name: '

    # input and check the environment
    while xconfiguration.environment == '':
        xconfiguration.environment = input(input_text)
        if xconfiguration.environment not in environments_list:
            print(xlib.get_separator())
            anwser = input(f'{xconfiguration.environment} is not a recorded environment. Do you like to record it? (Y/N): ')
            if anwser not in ['Y', 'y']:
                xconfiguration.environment = ''
            else:
                (OK, error_list) = xconfiguration.add_environment(xconfiguration.environment)
                if not OK:
                    for error in error_list:
                        print(error)
                    raise xlib.ProgramException('C002')

    # check if it is necesary to create the NGScloud config file corresponding to the environment
    if not xconfiguration.is_ngscloud_config_file_created():

        print(xlib.get_separator())
        print('Creating the config files ...')

        # create the NGScloud config file
        form_create_ngscloud_config_file(is_menu_call=False)

        # create the key pairs directory
        if not os.path.exists(xlib.get_keypairs_dir()):
            os.makedirs(xlib.get_keypairs_dir())

        # create the Bowtie2 config file
        (OK, error_list) = xbowtie2.create_bowtie2_config_file()

        # create the BUSCO config file
        (OK, error_list) = xbusco.create_busco_config_file()

        # create the CD-HIT-EST config file
        (OK, error_list) = xcdhit.create_cd_hit_est_config_file()

        # create the Cuffdiff config file
        (OK, error_list) = xcufflinks.create_cuffdiff_config_file()

        # create the Cufflinks-Cuffmerge config file
        (OK, error_list) = xcufflinks.create_cufflinks_cuffmerge_config_file()

        # create the Cuffquant config file
        (OK, error_list) = xcufflinks.create_cuffquant_config_file()

        # create the cutadapt config file
        (OK, error_list) = xcutadapt.create_cutadapt_config_file()

        # create the ddRADseq simulation config file
        (OK, error_list) = xddradseqtools.create_ddradseq_simulation_config_file()

        # create the eXpress config file
        (OK, error_list) = xexpress.create_express_config_file()

        # create the FastQC config file
        (OK, error_list) = xfastqc.create_fastqc_config_file()

        # create the Genome-guided Trinity config file
        (OK, error_list) = xtrinity.create_ggtrinity_config_file()

        # create the GMAP config file
        (OK, error_list) = xgmap.create_gmap_config_file()

        # create the GSNAP config file
        (OK, error_list) = xgmap.create_gsnap_config_file()

        # create the HISAT2 config file
        (OK, error_list) = xhisat2.create_hisat2_config_file()

        # create the htseq-count config file
        (OK, error_list) = xhtseq.create_htseq_count_config_file()

        # create the insilico_read_normalization config file
        (OK, error_list) = xtrinity.create_insilico_read_normalization_config_file()

        # create the ipyrad config file
        (OK, error_list) = xipyrad.create_ipyrad_config_file()

        # create the kallisto config file
        (OK, error_list) = xkallisto.create_kallisto_config_file()

        # create the QUAST config file
        (OK, error_list) = xquast.create_quast_config_file()

        # create the RADdesigner config file
        (OK, error_list) = xraddesigner.create_raddesigner_config_file()

        # create the REF-EVAL config file
        (OK, error_list) = xdetonate.create_ref_eval_config_file()

        # create the rnaQUAST config file
        (OK, error_list) = xrnaquast.create_rnaquast_config_file()

        # create the RSEM-EVAL config file
        (OK, error_list) = xdetonate.create_rsem_eval_config_file()

        # create the rsitesearch config file
        (OK, error_list) = xddradseqtools.create_rsitesearch_config_file()

        # create the SOAPdenovo2 config file
        (OK, error_list) = xsoapdenovo2.create_soapdenovo2_config_file()

        # create the SOAPdenovo-Trans config file
        (OK, error_list) = xsoapdenovotrans.create_soapdenovotrans_config_file()

        # create the STAR config file
        (OK, error_list) = xstar.create_star_config_file()

        # create the starcode config file
        (OK, error_list) = xstarcode.create_starcode_config_file()

        # create the TopHat config file
        (OK, error_list) = xtophat.create_tophat_config_file()

        # create the Trans-ABySS config file
        (OK, error_list) = xtransabyss.create_transabyss_config_file()

        # create the transcript-filter config file
        (OK, error_list) = xngshelper.create_transcript_filter_config_file()

        # create the transcriptome-blastx config file
        (OK, error_list) = xngshelper.create_transcriptome_blastx_config_file()

        # create the Transrate config file
        (OK, error_list) = xtransrate.create_transrate_config_file()

        # create the Trimmomatic config file
        (OK, error_list) = xtrimmomatic.create_trimmomatic_config_file()

        # create the Trinity config file
        (OK, error_list) = xtrinity.create_trinity_config_file()

        # create the Variant calling config file
        (OK, error_list) = xddradseqtools.create_variant_calling_config_file()

        # create the RAD-seq data files
        (OK, error_list) = xddradseqtools.create_end_file()
        (OK, error_list) = xddradseqtools.create_individual_file()
        (OK, error_list) = xddradseqtools.create_restriction_site_file()
        (OK, error_list) = xngshelper.create_vcf_sample_file()
        (OK, error_list) = xraddesigner.create_condition_file()

        # create the TOA config and data files
        (OK, error_list) = xtoa.create_toa_config_file()
        (OK, error_list) = xtoa.create_dataset_file()
        (OK, error_list) = xtoa.create_species_file()
        (OK, error_list) = xtoa.create_pipeline_config_file(pipeline_type=xlib.get_toa_process_pipeline_nucleotide_code())
        (OK, error_list) = xtoa.create_pipeline_config_file(pipeline_type=xlib.get_toa_process_pipeline_aminoacid_code())

        # create the transfer config files
        (OK, error_list) = xreference.create_reference_transfer_config_file()
        (OK, error_list) = xdatabase.create_database_transfer_config_file()
        (OK, error_list) = xread.create_read_transfer_config_file()
        (OK, error_list) = xresult.create_result_transfer_config_file(status='uncompressed')

        # create the gzip config files
        (OK, error_list) = xgzip.create_gzip_config_file(dataset_type='reference')
        (OK, error_list) = xgzip.create_gzip_config_file(dataset_type='database')
        (OK, error_list) = xgzip.create_gzip_config_file(dataset_type='read')
        (OK, error_list) = xgzip.create_gzip_config_file(dataset_type='result')

    # set the environment variables corresponding to the NGScloud config file, the AWS access key identification,
    # AWS secret access key and the current region name
    print(xlib.get_separator())
    print('Setting the environment variables ...')
    xconfiguration.set_environment_variables()
    print('The environment variables are set.')

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_create_ngscloud_config_file(is_menu_call):
    '''
    Create the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable
    OK = True

    # print the header
    if is_menu_call:
        clib.clear_screen()
        clib.print_headers_with_environment(f'Configuration - Recreate {xlib.get_project_name()} config file')

    # get current region and zone names
    region_name = xconfiguration.get_current_region_name()
    zone_name = xconfiguration.get_current_zone_name()

    # get basic AWS data and contact e-mail address from NGScloud config file
    (user_id, access_key_id, secret_access_key) = xconfiguration.get_basic_aws_data()
    email = xconfiguration.get_contact_data()

    # confirm or change the AWS data and contact e-mail address
    print(xlib.get_separator())
    user_id = cinputs.input_user_id(user_id)
    access_key_id = cinputs.input_access_key_id(access_key_id)
    secret_access_key = cinputs.input_secret_access_key(secret_access_key)
    email = cinputs.input_email(email)

    # check the AWS access key identification and the AWS secret access key   
    print(xlib.get_separator())
    print('Checking the AWS access key identification and the AWS secret access key')
    OK = xec2.check_aws_credentials(access_key_id, secret_access_key)
    if OK:
        print('The credentials are OK.')
    else:
        print('ERROR: The credentials are wrong. Please review your access key identification and secret access key in the AWS web.')
        if not is_menu_call:
            raise xlib.ProgramException('EXIT')

    # confirm the creation of the NGScloud config file
    if OK:
        if is_menu_call:
            print(xlib.get_separator())
            OK = clib.confirm_action(f'The {xlib.get_project_name()} config file is going to be created. The previous files will be lost.')

    # create the NGScloud config file corresponding to the environment
    if OK:
        print(xlib.get_separator())
        print(f'The file {xconfiguration.get_ngscloud_config_file()} is being created ...')
        (OK, error_list) = xconfiguration.create_ngscloud_config_file(user_id, access_key_id, secret_access_key, email)
        if OK:
            print('The config file is created with default values.')
            print()
            print('You can modify the conection data and contact e-mail address in:')
            print('    "Cloud control" -> "Configuration" -> "Update connection data and contact e-mail"')
            print()
            print(f'The assigned region and zone are {xconfiguration.get_default_region_name()} and {xconfiguration.get_default_zone_name()}, respectively. You can modify them in:')
            print('    "Cloud control" -> "Configuration" -> "Update region and zone data"')
        else:
            for error in error_list:
                print(error)
            raise xlib.ProgramException('C001')

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_view_ngscloud_config_file():
    '''
    List the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable
    OK = True

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # view the file
    text = f'Configuration - View {xlib.get_project_name()} config file'
    OK = clib.view_file(ngscloud_config_file, text)

    # show continuation message 
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_instance_types():
    '''
    List the characteristics of the instance types.
    '''

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Configuration - List instance types')

    # get the instance type dictionary
    instance_type_dict = xconfiguration.get_instance_type_dict(xconfiguration.get_cluster_mode_native())

    # set data width
    use_width = 17
    instance_type_width = 13
    vcpu_width = 5
    memory_width = 9
    nitro_width = 5
    starcluster_width = 13
    generation_width = 11

    # set line
    line = '{0:' + str(use_width) + '}   {1:' + str(instance_type_width) + '}   {2:' + str(vcpu_width) + '}   {3:>' + str(memory_width) + '}   {4:' + str(nitro_width) + '}   {5:' + str(starcluster_width) + '}   {6:' + str(generation_width) + '}'

    # print header
    print(line.format('Use', 'Instance type', 'vCPUs', 'Memory', 'Nitro', 'StarCluster', 'Generation'))
    print(line.format('=' * use_width, '=' * instance_type_width, '=' * vcpu_width, '=' * memory_width, '=' * nitro_width, '=' * starcluster_width, '=' * generation_width))

    # print detail lines
    for key in sorted(instance_type_dict.keys()):
        use = instance_type_dict[key]['use']
        instance_type = instance_type_dict[key]['id']
        vcpu = instance_type_dict[key]['vcpu']
        memory = f'{instance_type_dict[key]["memory"]} GiB'
        nitro = instance_type_dict[key]['nitro']
        starcluster = instance_type_dict[key]['starcluster']
        generation = instance_type_dict[key]['generation']
        print(line.format(use, instance_type, vcpu, memory, nitro, starcluster, generation))

    # show warnings about characteristics and pricing
    print(xlib.get_separator())
    print('You can consult the characteristics of the EC2 intance types in:')
    print('    https://aws.amazon.com/ec2/instance-types/')
    print('and the EC2 pricing is detailed in:')
    print('    https://aws.amazon.com/ec2/pricing/')

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_update_connection_data():
    '''
    Update the user id, access key id,  secret access key and contact e-mail address
    in the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Configuration - Update connection data')

    # get basic AWS data and contact e-mail address from NGScloud config file
    (user_id, access_key_id, secret_access_key) = xconfiguration.get_basic_aws_data()
    email = xconfiguration.get_contact_data()

    # input the new AWS data and the contact e-mail address
    print(xlib.get_separator())
    user_id = cinputs.input_user_id(user_id)
    access_key_id = cinputs.input_access_key_id(access_key_id)
    secret_access_key = cinputs.input_secret_access_key(secret_access_key)
    email = cinputs.input_email(email)

    # check the AWS access key identification and the AWS secret access key   
    print(xlib.get_separator())
    print('Checking the AWS access key identification and the AWS secret access key')
    OK = xec2.check_aws_credentials(access_key_id, secret_access_key)
    if OK:
        print('The credentials are OK.')
    else:
        print('ERROR: The credentials are wrong. Please review your access key identification and secret access key in the AWS web.')

    # get the NGScloud config file
    if OK:
        ngscloud_config_file = xconfiguration.get_ngscloud_config_file()
    
    # confirm the connection data update in the NGScloud config file
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The file {ngscloud_config_file} is going to be update with the new connection data.')

    # save the option dictionary in the NGScloud config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {ngscloud_config_file} is being update with the new connection data ...')
        (OK, error_list) = xconfiguration.update_connection_data(user_id, access_key_id, secret_access_key)
        if OK:
            (OK, error_list) = xconfiguration.update_contact_data(email)
        if OK:
            print('The config file has been update.')
        else:
            for error in error_list:
                print(error)
            raise xlib.ProgramException('C001')

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_update_region_zone():
    '''
    Update the current region and zone names in the NGScloud config file
    corresponding to the envoronment.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Configuration - Update region and zone')

    # input new current region and zone
    print(xlib.get_separator())
    region_name = cinputs.input_region_name(xconfiguration.get_current_region_name(), help=True)
    zone_name = cinputs.input_zone_name(region_name, xconfiguration.get_current_zone_name(), help=True)
  
    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()
  
    # confirm the region and zone update in the NGScloud config file
    print(xlib.get_separator())
    OK = clib.confirm_action(f'The file {ngscloud_config_file} is going to be update with the new region and zone.')

    # save the option dictionary in the NGScloud config file
    if OK:
        print(xlib.get_separator())
        print(f'The file {ngscloud_config_file} is being update with the new region and zone ...')
        (OK, error_list) = xconfiguration.update_region_zone_data(region_name, zone_name)
        if OK:
            for error in error_list:
                print(error)
            print('The config file has been update.')
        else:
            for error in error_list:
                print(error)
            raise xlib.ProgramException('C001')

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_link_volumes():
    '''
    Link volumes in the NGScloud config file corresponding to the environment.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Configuration - Link volumes')

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # check there are created volumes in the zone
    if xec2.get_created_volume_name_list(zone_name) == []:
        print(f'*** WARNING: There is not any volume created in the zone {zone_name}.')
        OK = False

    # get the dataset structure
    if OK:
        print(xlib.get_separator())
        dataset_structure = cinputs.input_code(text='Dataset structure', code_list=xconfiguration.get_dataset_structure_list(), default_code=xconfiguration.get_dataset_structure_singlevolume())

    # get volume names
    if OK:

        # when there is not any volume linked
        if dataset_structure.lower() == xconfiguration.get_dataset_structure_none():
            ngscloud_volume = ''
            app_volume = ''
            database_volume = ''
            read_volume = ''
            reference_volume = ''
            result_volume = ''

        # when the dataset structure is single-volume
        elif dataset_structure.lower() == xconfiguration.get_dataset_structure_singlevolume():
            ngscloud_volume = cinputs.input_volume_name(text='NGScloud volume name', zone_name=zone_name, type='created', allowed_none=False, help=True)
            app_volume = ''
            database_volume = ''
            read_volume = ''
            reference_volume = ''
            result_volume = ''

        # when the dataset structure is multi-volume
        elif dataset_structure.lower() == xconfiguration.get_dataset_structure_multivolume():

            # input volume names
            ngscloud_volume = ''
            app_volume = cinputs.input_volume_name(text='Application volume name', zone_name=zone_name, type='created', allowed_none=False, help=True)
            database_volume = cinputs.input_volume_name(text='Database volume name', zone_name=zone_name, type='created', allowed_none=True, help=True)
            read_volume = cinputs.input_volume_name(text='Read volume name', zone_name=zone_name, type='created', allowed_none=True, help=True)
            reference_volume = cinputs.input_volume_name(text='Reference volume name', zone_name=zone_name, type='created', allowed_none=True, help=True)
            result_volume = cinputs.input_volume_name(text='Result volume name', zone_name=zone_name, type='created', allowed_none=True, help=True)

        # check volume names
        used_volume_name_list = []

        if ngscloud_volume == 'NONE': ngscloud_volume = ''
        if app_volume == 'NONE': app_volume = ''
        if database_volume == 'NONE': database_volume = ''
        if read_volume == 'NONE': read_volume = ''
        if reference_volume == 'NONE': reference_volume = ''
        if result_volume == 'NONE': result_volume = ''

        if app_volume != '':
            used_volume_name_list.append(app_volume)
        if database_volume != '':
            if database_volume not in used_volume_name_list:
                used_volume_name_list.append(database_volume)
            else:
                OK = False
        if OK and read_volume != '':
            if read_volume not in used_volume_name_list:
                used_volume_name_list.append(read_volume)
            else:
                OK = False
        if OK and reference_volume != '':
            if reference_volume not in used_volume_name_list:
                used_volume_name_list.append(reference_volume)
            else:
                OK = False
        if OK and result_volume != '':
            if result_volume not in used_volume_name_list:
                used_volume_name_list.append(result_volume)
            else:
                OK = False

            if not OK:
                print('*** ERROR: A volume can be linked only once.')

    # confirm the linkage of volumes
    if OK:
        OK = clib.confirm_action('The dataset structure are going to be modified.')

    # link volumes
    if OK:
        devstdout = xlib.DevStdOut(xconfiguration.link_volumes.__name__)
        (OK, error_list) = xconfiguration.link_volumes(dataset_structure.lower(), ngscloud_volume, app_volume, database_volume, read_volume, reference_volume, result_volume, devstdout, function=None)
        if OK:
            for error in error_list:
                print(error)
            print('The data structure is modified.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_keypairs():
    '''
    List the key pairs of a region.
    '''

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Security - List key pairs')

    # get the key pair dictionary and the keypair names list
    keypairs_dict = xec2.get_keypair_dict(xconfiguration.get_current_region_name())
    keypair_keys_list = sorted(keypairs_dict.keys())
 
    # list keypairs
    print(xlib.get_separator())
    if keypair_keys_list == []:
        print(f'WARNING: There is not any keypair created in the region {xconfiguration.get_current_region_name()}.')
    else:
        # set data width
        keypair_name_width = 25
        fingerprint_width = 59
        # set line
        line = '{0:' + str(keypair_name_width) + '}   {1:' + str(fingerprint_width) + '}'
        # print header
        print(line.format('Key Pair Name', 'Fingerprint'))
        print(line.format('=' * keypair_name_width, '=' * fingerprint_width))
        # print detail lines
        for keypair_key in keypair_keys_list:
            keypair_name = keypairs_dict[keypair_key]['keypair_name']
            fingerprint = keypairs_dict[keypair_key]['fingerprint']
            print(line.format(keypair_name, fingerprint))

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_create_keypairs():
    '''
    Create the key pairs of a region.
    '''

    # initialize the control variable
    OK = True

    # print the header and get the cluster name
    clib.clear_screen()
    clib.print_headers_with_environment('Security - Create key pairs')

    # get current region name
    region_name = xconfiguration.get_current_region_name()

    # confirm the creation of the key pairs
    print(xlib.get_separator())
    OK = clib.confirm_action(f'The key pairs of the region {region_name} are going to be created.')

    # create key pairs
    if OK:
        print(xlib.get_separator())
        print(f'The key pairs of the region {region_name} are been created ...')
        (OK, error_list) = xec2.create_keypairs(region_name)
        if OK:
            print('The key pairs and their corresponding local files have been created.')
        else:
            for error in error_list:
                print(error)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_clusters():
    '''
    List clusters.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cluster operation - List clusters')

    # list clusters
    devstdout = xlib.DevStdOut(xcluster.list_clusters.__name__)
    OK = xcluster.list_clusters(devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_create_cluster():
    '''
    Create a cluster with one instance type.
    '''

    # initialize the control variable
    OK = True

    # get current region name
    region_name = xconfiguration.get_current_region_name()

    # get the NGScloud config file
    ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

    # get the option dictionary corresponding to the NGScloud config file
    ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

    # get the dataset structure and NGScloud_volume
    dataset_structure = ngscloud_options_dict['dataset info']['dataset_structure'].lower()

    # get the volume type dictionary
    volume_type_dict = xec2.get_volume_type_dict()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cluster operation - Create cluster')

    # check if there is a cluster running
    running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
    if len(running_cluster_list) > 0:
        print(xlib.get_separator())
        print('WARNING: There is already a cluster running in this environment.')
        OK = False

    # check the dataset structure
    if OK:
        if dataset_structure == xconfiguration.get_dataset_structure_none():
            print(xlib.get_separator())
            print(f'The dataset structure is not {xconfiguration.get_dataset_structure_singlevolume()} nor {xconfiguration.get_dataset_structure_multivolume()}, then datasets will be created in the root volume and they will be lost when the cluster is terminated.')
            print()
            print('The dataset structure can be modified in:')
            print('    "Cloud control" -> "Configuration" -> "Links volumes"')
            OK = clib.confirm_action('')

    # show sites related to EBS volumes
    if OK:
        print(xlib.get_separator())
        print('You can consult the characteristics of the EC2 intance types in:')
        print('    https://aws.amazon.com/ec2/instance-types/')
        print('and the EC2 pricing is detailed in:')
        print('    https://aws.amazon.com/ec2/pricing/')

    # check if StarCluster is installed
    if OK:
        is_starcluster_installed = True
        command = f'{xlib.get_starcluster()} --version'
        devstdout = xlib.DevStdOut('starcluster_version', print_stdout=False)
        rc = xlib.run_command(command, devstdout)
        if rc != 0:
            is_starcluster_installed = False
        else:
            with open(devstdout.get_log_file(), 'r') as log_command:
                version_found = False
                for line in log_command:
                    if line.startswith('0.95.6'):
                        version_found = True
                if not version_found:
                    is_starcluster_installed = False

    # get the cluster mode list
    if OK:
        cluster_mode_list = xconfiguration.get_cluster_mode_list()
        if not is_starcluster_installed:
            cluster_mode_list.remove(xconfiguration.get_cluster_mode_starcluster())

    # get the cluster mode
    if OK:
        print(xlib.get_separator())
        cluster_mode = cinputs.input_code(text='Cluster mode', code_list=cluster_mode_list, default_code=xconfiguration.get_cluster_mode_native())

    # check the AMI identification is create in the current zone
    if OK:
        if cluster_mode == xconfiguration.get_cluster_mode_native() and xec2.get_ubuntu_ami_id(region_name) == xec2.get_unknown_ami_id():
            print(f'*** ERROR: The AMI {xec2.get_ubuntu_ami_name()} is not found in the region {region_name}.')
            OK = False
        elif cluster_mode == xconfiguration.get_cluster_mode_starcluster() and xec2.get_starcluster_ami_id(region_name) == xec2.get_unknown_ami_id():
            print(f'*** ERROR: The AMI {xec2.get_starcluster_ami_name()} is not found in the region {region_name}.')
            OK = False

    # get the instance type
    if OK:
        instance_type = cinputs.input_instance_type(cluster_mode, help=True)

    # get the root volume type, volume size and pruchasing option
    if OK:
        if cluster_mode == xconfiguration.get_cluster_mode_native():
            volume_type = cinputs.input_code(text='Root volume type', code_list=xec2.get_volume_type_id_list(only_possible_root_disk=True), default_code='gp2')
            volume_size = cinputs.input_int(text='Root volume size', default=8, minimum=8, maximum=volume_type_dict[volume_type]['maximum_size'])
            purchasing_option = cinputs.input_code(code_list=xec2.get_purchasing_option_list(), default_code=xec2.get_purchasing_option_ondemand(), text='Purchasing option')

    # get the maximum splot price and interruption behavior
    if OK:
        if cluster_mode == xconfiguration.get_cluster_mode_native():
            if purchasing_option == xec2.get_purchasing_option_spot():
                max_spot_price = cinputs.input_float(text='Maximum spot_price (in $; 0.0 is on-demand price)', default=0.0, minimum=0.0)
                interruption_behavior = cinputs.input_code(code_list=xec2.get_interruption_behavior_list(), default_code=xec2.get_interruption_behavior_terminate(), text='Interruption behavior')
            else:
                max_spot_price = 0.0
                interruption_behavior = xec2.get_interruption_behavior_terminate()

    # confirm the creation of the cluster
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The cluster is going to be created.')

    # create the cluster
    if OK:
        if cluster_mode == xconfiguration.get_cluster_mode_native():
            devstdout = xlib.DevStdOut(xinstance.create_instance.__name__)
            OK = xinstance.create_instance(instance_type, log=sys.stdout, root_volume_type=volume_type, root_volumen_size=volume_size, purchasing_option=purchasing_option, max_spot_price=max_spot_price, interruption_behavior=interruption_behavior, function=None, is_menu_call=True)
        elif cluster_mode == xconfiguration.get_cluster_mode_starcluster():
            template_name = xconfiguration.build_cluster_name(instance_type)
            cluster_name = xconfiguration.build_cluster_name(instance_type)
            devstdout = xlib.DevStdOut(xcluster.create_cluster.__name__)
            (OK, master_state_code, master_state_name) = xcluster.create_cluster(template_name, cluster_name, devstdout, function=None, is_menu_call=True)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_terminate_cluster(force):
    '''
    Terminate a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    if not force:
        clib.print_headers_with_environment('Cluster operation - Terminate cluster')
    else:
        clib.print_headers_with_environment('Cluster operation - Force termination of a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if not force:
        if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
            print('WARNING: There is not any running cluster.')
            OK = False
        else:
            cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        instance_type = cinputs.input_instance_type(cluster_mode=xconfiguration.get_cluster_mode_native(), help=True)
        cluster_name = xconfiguration.build_cluster_name(instance_type)

    # confirm the termination of the cluster
    if OK:
        print(xlib.get_separator())
        if not force:
            OK = clib.confirm_action('The cluster is going to be terminated.')
        else:
            OK = clib.confirm_action('The cluster is going to be forced to terminate.')

    # terminate the cluster
    if OK:
        if not force:
            devstdout = xlib.DevStdOut(xcluster.terminate_cluster.__name__)
            if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native():
                OK = xinstance.terminate_instance(cluster_name, devstdout, function=None, is_menu_call=True)
            elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
                OK = xcluster.terminate_cluster(cluster_name, force, devstdout, function=None, is_menu_call=True)
        else:
            devstdout = xlib.DevStdOut(xcluster.terminate_cluster.__name__)
            if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native() or xec2.get_cluster_mode(cluster_name) is None:
                OK = xinstance.terminate_instance(cluster_name, devstdout, function=None, is_menu_call=True)
            elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
                OK = xcluster.terminate_cluster(cluster_name, force, devstdout, function=None, is_menu_call=True)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_show_cluster_composition():
    '''
    Show cluster information of every node: OS, CPU number and memory.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cluster operation - Show cluster composition')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # check the cluster mode
    if xec2.get_cluster_mode(cluster_name) != xconfiguration.get_cluster_mode_starcluster():
        print(f'WARNING: This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.')
        OK = False

    # show the status of batch jobs
    if OK:
        devstdout = xlib.DevStdOut(xcluster.show_cluster_composition.__name__)
        xcluster.show_cluster_composition(cluster_name, devstdout, function=None)

    # show continuation message
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_show_status_batch_jobs():
    '''
    Show the status of batch jobs in the cluster.
    '''

    # initialize the control variable
    OK = True

    # initialize the list of identification of the batch jobs
    batch_job_id_list = []

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cluster operation - Show status batch jobs')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # check the cluster mode
    if xec2.get_cluster_mode(cluster_name) != xconfiguration.get_cluster_mode_starcluster():
        print(f'WARNING: This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the batch job dictionary
    if OK:
        (OK, error_list, batch_job_dict) = xcluster.get_batch_job_dict(ssh_client)

    # build the list of identifications of the batch jobs
    if OK:
        for job_id in batch_job_dict.keys():
            batch_job_id_list.append(job_id)
        if batch_job_id_list != []:
            batch_job_id_list.sort()
        else:
            print('WARNING: There is not any batch job.')
            OK = False

    # print the batch jobs
    if OK:
        print(xlib.get_separator())
        # set data width
        job_id_width = 6
        job_name_width = 10
        state_width = 15
        start_date_width = 10
        start_time_width = 10
        # set line
        line = '{0:' + str(job_id_width) + '} {1:' + str(job_name_width) + '} {2:' + str(state_width) + '} {3:' + str(start_date_width) + '} {4:' + str(start_time_width) + '}'
        # print header
        print(line.format('Job id', 'Job name', 'State', 'Start date', 'Start time'))
        print(line.format('=' * job_id_width, '=' * job_name_width, '=' * state_width, '=' * start_date_width, '=' * start_time_width))
        # print detail lines
        for job_id in batch_job_id_list:
            job_name = batch_job_dict[job_id]['job_name']
            state = batch_job_dict[job_id]['state']
            start_date = batch_job_dict[job_id]['start_date']
            start_time = batch_job_dict[job_id]['start_time']
            print(line.format(job_id, job_name, state, start_date, start_time))

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_kill_batch_job():
    '''
    Kill a batch job in the cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cluster operation - Kill batch job')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) == []:
        print('WARNING: There is not any running cluster.')
        OK = False
    else:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)

    # check the cluster mode
    if xec2.get_cluster_mode(cluster_name) != xconfiguration.get_cluster_mode_starcluster():
        print(f'WARNING: This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.')
        OK = False

    # create the SSH client connection
    if OK:
        (OK, error_list, ssh_client) = xssh.create_ssh_client_connection(cluster_name)
        for error in error_list:
            print(error)

    # get the batch job identificaction
    if OK:
        batch_job_id = cinputs.input_batch_job_id(ssh_client, help=True)
        if batch_job_id == '':
            print('WARNING: There is not any batch job.')
            OK = False

    # confirm the kill of the batch job
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action(f'The batch job {batch_job_id} is going to be killed.')

    # kill the batch job
    if OK:
        devstdout = xlib.DevStdOut(xcluster.kill_batch_job.__name__)
        xcluster.kill_batch_job(cluster_name, batch_job_id, devstdout, function=None)

    # close the SSH client connection
    if OK:
        xssh.close_ssh_client_connection(ssh_client)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_nodes():
    '''
    List nodes running.
    '''

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Node operation - List nodes')

    # get the node dictionary and node key list
    node_dict = xec2.get_node_dict()
    node_key_list = sorted(node_dict.keys())

    # list nodes
    print(xlib.get_separator())
    if node_key_list == []:
        print('WARNING: There is not any node running.')
    else:
        # set data width
        security_group_name_width = 28
        zone_name_width = 20
        node_name_width = 20
        node_id_width = 19
        state_width = 20
        # set line
        line = '{0:' + str(security_group_name_width) + '}   {1:' + str(zone_name_width) + '}   {2:' + str(node_name_width) + '}   {3:' + str(node_id_width) + '}   {4:' + str(state_width) + '}'
        # print header
        print(line.format('Security Group', 'Zone', 'Node Name', 'Node Id', 'State'))
        print(line.format('=' * security_group_name_width, '=' * zone_name_width, '=' * node_name_width, '=' * node_id_width, '=' * state_width))
        # print detail lines
        for node_key in node_key_list:
            security_group_name = node_dict[node_key]['security_group_name']
            zone_name = node_dict[node_key]['zone_name']
            node_name = node_dict[node_key]['node_name']
            node_id = node_dict[node_key]['node_id'] 
            state = node_dict[node_key]['state']
            print(line.format(security_group_name, zone_name, node_name, node_id, state))

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_add_node():
    '''
    Add a node in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Node operation - Add node in a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # check the cluster mode
    if OK:
        if xec2.get_cluster_mode(cluster_name) != xconfiguration.get_cluster_mode_starcluster():
            print(f'WARNING: This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.')
            OK = False

    # check the instance number
    if OK:
        if len(xec2.get_cluster_node_list(cluster_name)) >= xec2.get_max_node_number():
            print(f'WARNING: The maximum number ({xec2.get_max_node_number()}) of instances is already running.')
            OK = False

    # get the node name
    if OK:
        node_name = cinputs.input_node_name(cluster_name, new=True, is_master_valid=False, help=True)

    # confirm the addition of the node in the cluster
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The node is going to be added.')

    # add node in cluster
    if OK:
        devstdout = xlib.DevStdOut(xnode.add_node.__name__)
        xnode.add_node(cluster_name, node_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_node():
    '''
    Remove a node in a cluster.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Node operation - Remove node in a cluster')

    # get the cluster name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # check the cluster mode
    if OK:
        if xec2.get_cluster_mode(cluster_name) != xconfiguration.get_cluster_mode_starcluster():
            print(f'WARNING: This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.')
            OK = False

    # get the node name
    if OK:
        node_name = cinputs.input_node_name(cluster_name, new=False, is_master_valid=False, help=True)
        if node_name == []:
            print('WARNING: There is not any running node besides the master.')
            OK = False

    # confirm the removal of the node in the cluster
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The node is going to be removed.')

    # remove node
    if OK:
        devstdout = xlib.DevStdOut(xnode.remove_node.__name__)
        xnode.remove_node(cluster_name, node_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_list_volumes():
    '''
    List volumes created.
    '''

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - List volumes')

    # get the volume dictionary and volume key list
    volumes_dict = xec2.get_volume_dict()
    volume_keys_list = sorted(volumes_dict.keys())

    # list volume
    print(xlib.get_separator())
    if volume_keys_list == []:
        print('WARNING: There is not any volume created.')
    else:
        # set data width
        zone_name_width = 20
        volume_name_width = 20
        volume_id_width = 21
        size_width = 10
        state_width = 10
        attachments_number_width = 11
        # set line
        line = '{0:' + str(zone_name_width) + '}   {1:' + str(volume_name_width) + '}   {2:' + str(volume_id_width) + '}   {3:' + str(size_width) + '}   {4:' + str(state_width) + '}   {5:' + str(attachments_number_width) + '}'
        # print header
        print(line.format('Zone', 'Volume Name', 'Volume Id', 'Size (GiB)', 'State', 'Attachments'))
        print(line.format('=' * zone_name_width, '=' * volume_name_width, '=' * volume_id_width, '=' * size_width, '=' * state_width, '=' * attachments_number_width))
        # print detail lines
        for volume_key in volume_keys_list:
            zone_name = volumes_dict[volume_key]['zone_name']
            volume_name = volumes_dict[volume_key]['volume_name'] 
            volume_id = volumes_dict[volume_key]['volume_id']
            size = volumes_dict[volume_key]['size']
            state = volumes_dict[volume_key]['state']
            attachments_number = volumes_dict[volume_key]['attachments_number']
            print(line.format(zone_name, volume_name, volume_id, size, state, attachments_number))

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_create_volume():
    '''
    Create a volume in the current zone.
    '''

    # initialize the control variable
    OK = True

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # get the volume type dictionary
    volume_type_dict = xec2.get_volume_type_dict()

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - Create volume')

    # show sites related to EBS volumes
    print(xlib.get_separator())
    print('You can consult the characteristics of the EBS volumes in:')
    print('    https://aws.amazon.com/ebs/details/')
    print('and the EBS pricing is detailed in:')
    print('    https://aws.amazon.com/ebs/pricing/')

    # get the cluster name, node name, volume name, volume type and volume size
    print(xlib.get_separator())
    volume_name = cinputs.input_volume_name(text='volume name', zone_name=zone_name, type='created', allowed_none=False, help=False)
    volume_type = cinputs.input_code(text='Volume type', code_list=xec2.get_volume_type_id_list(only_possible_root_disk=False), default_code='standard')
    minimum_size = volume_type_dict[volume_type]['minimum_size']
    maximum_size = volume_type_dict[volume_type]['maximum_size']
    volume_size = cinputs.input_int(text='Volume size', default=minimum_size , minimum=minimum_size , maximum=maximum_size)
    terminate_indicator = True if cinputs.input_code(text='Is necessary to terminate the volume creator?', code_list=['y', 'n'], default_code=None).lower() == 'y' else False

    # confirm the creation of the volume
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The volume is going to be created.')

    # create the volume
    if OK:
        devstdout = xlib.DevStdOut(xvolume.create_volume.__name__)
        OK = xvolume.create_volume(volume_name, volume_type, volume_size, terminate_indicator, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_remove_volume():
    '''
    Remove a volume in the current zone.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - Remove volume')

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # get the volume name
    print(xlib.get_separator())
    volume_name = cinputs.input_volume_name(text='volume name', zone_name=zone_name, type='created', allowed_none=False, help=True)

    # confirm the removal of the volume
    print(xlib.get_separator())
    OK = clib.confirm_action('The volume is going to be removed.')

    # remove the volume
    if OK:
        devstdout = xlib.DevStdOut(xvolume.remove_volume.__name__)
        OK = xvolume.remove_volume(volume_name, devstdout, function=None)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_terminate_volume_creator():
    '''
    Terminate de volume creator of the current zone.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - Terminate volume creator')

    # confirm the termination of the volume creator
    print(xlib.get_separator())
    OK = clib.confirm_action('The volume creator is going to be terminated.')

    # terminate the volume creator
    if OK:
        devstdout = xlib.DevStdOut(xcluster.terminate_cluster.__name__)
        OK = xinstance.terminate_volume_creator(devstdout, function=None, is_menu_call=False)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_mount_volume():
    '''
    Mount a volume in a node.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - Mount volume in a node')

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # get the cluster name and node name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
        node_name = cinputs.input_node_name(cluster_name, new=False, is_master_valid=True, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # get the volume name, AWS device file and directory path
    if OK:
        volume_name = cinputs.input_volume_name(text='volume name', zone_name=zone_name, type='created', allowed_none=False, help=True)
        aws_device_file = cinputs.input_device_file(node_name, volume_name)
        mount_path = cinputs.input_mount_path(node_name)

    # confirm the mounting of the volume
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The volume is going to be mounted.')

    # mount the volume in the node
    if OK:
        devstdout = xlib.DevStdOut(xvolume.mount_volume.__name__)
        xvolume.mount_volume(cluster_name, node_name, volume_name, aws_device_file, mount_path, devstdout, function=None, is_menu_call=True)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_unmount_volume():
    '''
    Unmount a volume in a node.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Volume operation - Unmount volume in a node')

    # get current zone name
    zone_name = xconfiguration.get_current_zone_name()

    # get the cluster name and node name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
        node_name = cinputs.input_node_name(cluster_name, new=False, is_master_valid=True, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # get the volume name and mount path
    if OK:
        volume_name = cinputs.input_volume_name(text='volume name', zone_name=zone_name, type='created', allowed_none=False, help=True)
        mount_path = cinputs.input_mount_path(node_name)

    # confirm the unmounting of the volume
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The volume is going to be unmounted.')

    # unmount the volume in the node
    if OK:
        devstdout = xlib.DevStdOut(xvolume.unmount_volume.__name__)
        xvolume.unmount_volume(cluster_name, node_name, volume_name, mount_path, devstdout, function=None, is_menu_call=True)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

def form_open_terminal():
    '''
    Open a terminal windows of a cluster node.
    '''

    # initialize the control variable
    OK = True

    # print the header
    clib.clear_screen()
    clib.print_headers_with_environment('Cloud control - Open a terminal')

    # get the cluster name and node name
    print(xlib.get_separator())
    if xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False) != []:
        cluster_name = cinputs.input_cluster_name(volume_creator_included=False, help=True)
        node_name = cinputs.input_node_name(cluster_name, new=False, is_master_valid=True, help=True)
    else:
        print('WARNING: There is not any running cluster.')
        OK = False

    # confirm the terminal opening
    if OK:
        print(xlib.get_separator())
        OK = clib.confirm_action('The terminal is going to be opened using StarCluster.')

    # open de terminal
    if OK:
        xcluster.open_terminal(cluster_name, node_name)

    # show continuation message 
    print(xlib.get_separator())
    input('Press [Intro] to continue ...')

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the functions related to forms corresponding to Cloud Control menu items in console mode.')
     sys.exit(0)

#-------------------------------------------------------------------------------
