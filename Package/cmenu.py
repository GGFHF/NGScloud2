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
This file contains the functions related to menus in console mode.
'''
#-------------------------------------------------------------------------------

import sys

import cbioinfoapp
import ccloud
import cdataset
import clib
import clog
import ctoa
import xddradseqtools
import xlib
import xngshelper
import xraddesigner
import xtoa

#-------------------------------------------------------------------------------

def build_menu_main():
    '''
    Build the menu Main.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Main')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Cloud control')
        print()
        print( '    2. De novo RNA-seq')
        print( '    3. Reference-based RNA-seq')
        print( '    4. RAD-seq')
        print( '    5. Taxonomy-oriented annotation')
        print()
        print( '    6. Datasets')
        print( '    7. Logs')
        print()
        print(f'    X. Exit {xlib.get_project_name()}')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cloud_control()
        elif option == '2':
            build_menu_denovo_rnaseq()
        elif option == '3':
            build_menu_reference_based_rnaseq()
        elif option == '4':
            build_menu_radseq()
        elif option == '5':
            build_menu_toa()
        elif option == '6':
            build_menu_datasets()
        elif option == '7':
            build_menu_logs()
        elif option == 'X':
            sure = ''
            print( '')
            while sure not in ['Y', 'N']:
                sure = input(f'Are you sure to exit {xlib.get_project_name()}? (y or n): ').upper()
            if sure == 'Y':
                break

#-------------------------------------------------------------------------------

def build_menu_cloud_control():
    '''
    Build the menu Cloud control.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Cloud control')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Set environment')
        print()
        print( '    2. Configuration')
        print( '    3. Security')
        print()
        print( '    4. Cluster operation')
        print( '    5. Node operation')
        print( '    6. Volume operation')
        print()
        print( '    7. Bioinfo software installation')
        print()
        print( '    8. Open a terminal')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_set_environment()
        elif option == '2':
            build_menu_configuration()
        elif option == '3':
            build_menu_security()
        elif option == '4':
            build_menu_cluster_operation()
        elif option == '5':
            build_menu_node_operation()
        elif option == '6':
            build_menu_volume_operation()
        elif option == '7':
            build_menu_bioinfo_software_installation()
        elif option == '8':
            ccloud.form_open_terminal()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_configuration():
    '''
    Build the menu Configuration.
    '''


    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Configuration')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. Recreate {xlib.get_project_name()} config file')
        print(f'    2. View {xlib.get_project_name()} config file')
        print()
        print( '    3. List instance types')
        print()
        print( '    4. Update connection data and contact e-mail')
        print( '    5. Update region and zone')
        print()
        print( '    6. Link volumes')
        print()
        print( '    X. Return to menu Cloud control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_create_ngscloud_config_file(is_menu_call=True)
        elif option == '2':
            ccloud.form_view_ngscloud_config_file()
        elif option == '3':
            ccloud.form_list_instance_types()
        elif option == '4':
            ccloud.form_update_connection_data()
        elif option == '5':
            ccloud.form_update_region_zone()
        elif option == '6':
            ccloud.form_link_volumes()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_security():
    '''
    Build the menu Security.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Security')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. List key pairs')
        print( '    2. Create key pairs')
        print()
        print( '    3. List cluster security groups (coming soon!)')
        print( '    4. Force removal of a cluster security group (coming soon!)')
        print()
        print( '    X. Return to menu Cloud control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_list_keypairs()
        elif option == '2':
            ccloud.form_create_keypairs()
        elif option == '3':
            pass
        elif option == '3':
            pass
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cluster_operation():
    '''
    Build the menu Cluster operation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Cluster operation')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. List clusters')
        print()
        print( '    2. Create cluster')
        print( '    3. Terminate cluster')
        print()
        print( '    4. Force termination of a cluster')
        print()
        print( '    5. Show cluster composition')
        print()
        print( '    6. Show status of batch jobs')
        print( '    7. Kill batch job')
        print()
        print( '    X. Return to menu Cloud Control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_list_clusters()
        elif option == '2':
            ccloud.form_create_cluster()
        elif option == '3':
            ccloud.form_terminate_cluster(force=False)
        elif option == '4':
            ccloud.form_terminate_cluster(force=True)
        elif option == '5':
            ccloud.form_show_cluster_composition()
        elif option == '6':
            ccloud.form_show_status_batch_jobs()
        elif option == '7':
            ccloud.form_kill_batch_job()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_node_operation():
    '''
    Build the menu Node operation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Node operation')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. List nodes')
        print()
        print( '    2. Add node in a cluster')
        print( '    3. Remove node in a cluster')
        print()
        print( '    X. Return to menu Cloud Control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_list_nodes()
        elif option == '2':
            ccloud.form_add_node()
        elif option == '3':
            ccloud.form_remove_node()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_volume_operation():
    '''
    Build the menu Volume operation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Volume operation')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. List volumes')
        print()
        print( '    2. Create volume')
        print( '    3. Remove volume')
        print()
        print( '    4. Terminate volume creator')
        print()
        print( '    5. Mount volume in a node')
        print( '    6. Unmount volume in a node')
        print()
        print( '    X. Return to menu Cloud Control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ccloud.form_list_volumes()
        elif option == '2':
            ccloud.form_create_volume()
        elif option == '3':
            ccloud.form_remove_volume()
        elif option == '4':
            ccloud.form_terminate_volume_creator()
        elif option == '5':
            ccloud.form_mount_volume()
        elif option == '6':
            ccloud.form_unmount_volume()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_bioinfo_software_installation():
    '''
    Build the menu Bioinfo software installation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Bioinfo software installation')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    0. {xlib.get_miniconda3_name()} (Bioconda infrastructure)')
        print()
        print(f'    1. {xlib.get_bcftools_name()}')
        print(f'    2. {xlib.get_bedtools_name()}')
        print(f'    3. {xlib.get_blastplus_name()}')
        print(f'    4. {xlib.get_bowtie2_name()}')
        print(f'    5. {xlib.get_busco_name()}')
        print(f'    6. {xlib.get_cd_hit_name()}')
        print(f'    7. {xlib.get_cufflinks_name()}')
        print(f'    8. {xlib.get_cutadapt_name()}')
        print(f'    9. {xlib.get_ddradseqtools_name()}')
        print(f'   10. {xlib.get_detonate_name()}')
        print(f'   11. {xlib.get_diamond_name()}')
        print(f'   12. {xlib.get_entrez_direct_name()}')
        print(f'   13. {xlib.get_express_name()}')
        print(f'   14. {xlib.get_fastqc_name()}')
        print(f'   15. {xlib.get_gmap_gsnap_name()}')
        print(f'   16. {xlib.get_hisat2_name()}')
        print(f'   17. {xlib.get_htseq_name()}')
        print(f'   18. {xlib.get_ipyrad_name()}')
        print(f'   19. {xlib.get_kallisto_name()}')
        print(f'   20. {xlib.get_ngshelper_name()}')
        print(f'   21. {xlib.get_quast_name()}')
        print(f'   22. {xlib.get_raddesigner_name()}')
        print(f'   23. {xlib.get_rnaquast_name()}')
        print(f'   24. {xlib.get_rsem_name()}')
        print(f'   25. {xlib.get_samtools_name()}')
        print(f'   26. {xlib.get_soapdenovo2_name()}')
        print(f'   27. {xlib.get_soapdenovotrans_name()}')
        print(f'   28. {xlib.get_star_name()}')
        print(f'   29. {xlib.get_starcode_name()}')
        print(f'   30. {xlib.get_toa_name()}')
        print(f'   31. {xlib.get_tophat_name()}')
        print(f'   32. {xlib.get_transabyss_name()}')
        print(f'   33. {xlib.get_transdecoder_name()}')
        print(f'   34. {xlib.get_transrate_name()}')
        print(f'   35. {xlib.get_trimmomatic_name()}')
        print(f'   36. {xlib.get_trinity_name()}')
        print(f'   37. {xlib.get_vcftools_name()}')
        print(f'   38. {xlib.get_vcftools_perl_libraries_name()}')
        print(f'   39. {xlib.get_vsearch_name()}')
        print()
        print(f'   40. {xlib.get_r_name()} & analysis packages')
        print()
        print( '    X. Return to menu Cloud Control')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '0':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_miniconda3_code())
        elif option == '1':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_bcftools_code())
        elif option == '2':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_bedtools_code())
        elif option == '3':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_blastplus_code())
        elif option == '4':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_bowtie2_code())
        elif option == '5':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_busco_code())
        elif option == '6':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_cd_hit_code())
        elif option == '7':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_cufflinks_code())
        elif option == '8':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_cutadapt_code())
        elif option == '9':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_ddradseqtools_code())
        elif option == '10':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_detonate_code())
        elif option == '11':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_diamond_code())
        elif option == '12':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_entrez_direct_code())
        elif option == '13':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_express_code())
        elif option == '14':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_fastqc_code())
        elif option == '15':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_gmap_gsnap_code())
        elif option == '16':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_hisat2_code())
        elif option == '17':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_htseq_code())
        elif option == '18':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_ipyrad_code())
        elif option == '19':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_kallisto_code())
        elif option == '20':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_ngshelper_code())
        elif option == '21':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_quast_code())
        elif option == '22':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_raddesigner_code())
        elif option == '23':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_rnaquast_code())
        elif option == '24':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_rsem_code())
        elif option == '25':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_samtools_code())
        elif option == '26':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_soapdenovo2_code())
        elif option == '27':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_soapdenovotrans_code())
        elif option == '28':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_star_code())
        elif option == '29':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_starcode_code())
        elif option == '30':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_toa_code())
        elif option == '31':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_tophat_code())
        elif option == '32':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_transabyss_code())
        elif option == '33':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_transdecoder_code())
        elif option == '34':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_transrate_code())
        elif option == '35':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_trimmomatic_code())
        elif option == '36':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_trinity_code())
        elif option == '37':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_vcftools_code())
        elif option == '38':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_vcftools_perl_libraries_code())
        elif option == '39':
            cbioinfoapp.form_installation_bioinfo_app(xlib.get_vsearch_code())
        elif option == '40':
             cbioinfoapp.form_installation_bioinfo_app(xlib.get_r_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq():
    '''
    Build the menu De novo RNA-seq.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Read quality')
        print( '    2. Trimming')
        print( '    3. Digital normalization')
        print()
        print( '    4. Assembly')
        print()
        print( '    5. Read alignment')
        print()
        print( '    6. Transcriptome quality assessment')
        print( '    7. Transcriptome filtering')
        print()
        print( '    8. Quantitation')
        print()
        print( '    9. Variant calling')
        print()
        print( '    A. Annotation')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_denovo_rnaseq_read_quality()
        elif option == '2':
            build_menu_denovo_rnaseq_trimming()
        elif option == '3':
            build_menu_denovo_rnaseq_digital_normalization()
        elif option == '4':
            build_menu_denovo_rnaseq_assembly()
        elif option == '5':
            build_menu_denovo_rnaseq_read_alignment()
        elif option == '6':
            build_menu_denovo_rnaseq_transcriptome_quality_assessment()
        elif option == '7':
            build_menu_denovo_rnaseq_transcriptome_transcriptome_filtering()
        elif option == '8':
            build_menu_denovo_rnaseq_quantitation()
        elif option == '9':
            build_menu_denovo_rnaseq_variant_calling()
        elif option == 'A':
            build_menu_denovo_rnaseq_annotation()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_read_quality():
    '''
    Build the menu De novo RNA-seq - Read quality.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Read quality')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_fastqc_name()}')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_fastqc()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_trimming():
    '''
    Build the menu De novo RNA-seq - Trimming.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Trimming')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cutadapt_name()}')
        print(f'    2. {xlib.get_trimmomatic_name()}')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cutadapt()
        elif option == '2':
            build_menu_trimmomatic()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_digital_normalization():
    '''
    Build the menu De novo RNA-seq - Digital normalization.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Digital normalization')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_insilico_read_normalization_name()} ({xlib.get_trinity_name()} package)')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_insilico_read_normalization()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_assembly():
    '''
    Build the menu De novo RNA-seq - Assembly.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Assembly')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_soapdenovotrans_name()}')
        print(f'    2. {xlib.get_transabyss_name()}')
        print(f'    3. {xlib.get_trinity_name()}')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_soapdenovotrans()
        elif option == '2':
            build_menu_transabyss()
        elif option == '3':
            build_menu_trinity()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_read_alignment():
    '''
    Build the menu denovo RNA-seq - Read alignment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Read alignment')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_bowtie2_name()}')
        print(f'    2. {xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_bowtie2()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_transcriptome_quality_assessment():
    '''
    Build the menu De novo RNA-seq- Transcriptome quality assessment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Transcriptome quality assessment')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_busco_name()}')
        print(f'    2. {xlib.get_quast_name()}')
        print(f'    3. {xlib.get_rnaquast_name()}')
        print(f'    4. {xlib.get_rsem_eval_name()} ({xlib.get_detonate_name()} package)')
        print(f'    5. {xlib.get_transrate_name()}')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_busco()
        elif option == '2':
            build_menu_quast()
        if option == '3':
            build_menu_rnaquast()
        elif option == '4':
            build_menu_rsem_eval()
        elif option == '5':
            build_menu_transrate()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_transcriptome_transcriptome_filtering():
    '''
    Build the menu De novo RNA-seq - Transcriptome filtering.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Transcriptome filtering')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cd_hit_est_name()} ({xlib.get_cd_hit_name()} package)')
        print(f'    2. {xlib.get_transcript_filter_name()} ({xlib.get_ngshelper_name()} package)')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cd_hit_est()
        elif option == '2':
            build_menu_transcript_filter()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_quantitation():
    '''
    Build the menu De novo RNA-seq - Quantitation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Quantitation')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_express_name()}')
        print(f'    2. {xlib.get_kallisto_name()}')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_express()
        elif option == '':
            build_menu_kallisto()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_variant_calling():
    '''
    Build the menu De novo RNA-seq - Variant calling.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Variant calling')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_variant_calling_name()} ({xlib.get_ddradseqtools_name()} package)')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_variant_calling()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_denovo_rnaseq_annotation():
    '''
    Build the menu De novo RNA-seq - Annotation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('De novo RNA-seq - Annotation')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_transcriptome_blastx_name()} ({xlib.get_ngshelper_name()} package)')
        print()
        print( '    X. Return to menu De novo RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_transcriptome_blastx()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq():
    '''
    Build the menu Reference-based RNA-seq.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Read quality')
        print( '    2. Trimming')
        print()
        print( '    3. Read alignment')
        print()
        print( '    4. Assembly')
        print()
        print( '    5. Transcriptome alignment')
        print()
        print( '    6. Quantitation')
        print( '    7. Differential expression')
        print()
        print( '    8. Variant calling')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_reference_based_rnaseq_read_quality()
        elif option == '2':
            build_menu_reference_based_rnaseq_trimming()
        elif option == '3':
            build_menu_reference_based_rnaseq_read_alignment()
        elif option == '4':
            build_menu_reference_based_rnaseq_assembly()
        elif option == '5':
            build_menu_reference_based_rnaseq_transcriptome_alignment()
        elif option == '6':
            build_menu_reference_based_rnaseq_quantitation()
        elif option == '7':
            build_menu_reference_based_rnaseq_differential_expression()
        elif option == '8':
            build_menu_reference_based_rnaseq_variant_calling()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_read_quality():
    '''
    Build the menu Reference-based RNA-seq - Read quality.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Read quality')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_fastqc_name()}')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_fastqc()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_trimming():
    '''
    Build the menu Reference-based RNA-seq - Trimming.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Trimming')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cutadapt_name()}')
        print(f'    2. {xlib.get_trimmomatic_name()}')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cutadapt()
        elif option == '2':
            build_menu_trimmomatic()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_read_alignment():
    '''
    Build the menu Reference-based RNA-seq - Read alignment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Read alignment')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_bowtie2_name()}')
        print(f'    2. {xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)')
        print(f'    3. {xlib.get_hisat2_name()}')
        print(f'    4. {xlib.get_star_name()}')
        print(f'    5. {xlib.get_tophat_name()}')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_bowtie2()
        if option == '2':
            build_menu_gsnap()
        elif option == '3':
            build_menu_hisat2()
        elif option == '4':
            build_menu_star()
        elif option == '5':
            build_menu_tophat()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_assembly():
    '''
    Build the menu Reference-based RNA-seq - Assembly.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Assembly')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cufflinks_cuffmerge_name()} ({xlib.get_cufflinks_name()} package)')
        print(f'    2. {xlib.get_ggtrinity_name()} ({xlib.get_trinity_name()} package)')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cufflinks_cuffmerge()
        elif option == '2':
            build_menu_ggtrinity()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_transcriptome_alignment():
    '''
    Build the menu Reference-based RNA-seq - Transcriptome alignment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Transcriptome alignment')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_gmap_name()} ({xlib.get_gmap_gsnap_name()} package)')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_gmap()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_quantitation():
    '''
    Build the menu Reference-based RNA-seq - Quantitation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Quantitation')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cuffquant_name()} ({xlib.get_cufflinks_name()} package)')
        print(f'    2. {xlib.get_htseq_count_name()} ({xlib.get_htseq_name()} package)')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cuffquant()
        elif option == '2':
            build_menu_htseq_count()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_differential_expression():
    '''
    Build the menu Reference-based RNA-seq - Differential expression.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Differential expression')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cuffdiff_name()} ({xlib.get_cufflinks_name()} package)')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cuffdiff()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_based_rnaseq_variant_calling():
    '''
    Build the menu Reference-based RNA-seq - Variant calling.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference-based RNA-seq - Variant calling')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_variant_calling_name()} ({xlib.get_ddradseqtools_name()} package)')
        print()
        print( '    X. Return to menu Reference-based RNA-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_variant_calling()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq():
    '''
    Build the menu RAD-seq.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Maintenance of data files')
        print()
        print( '    2. Enzyme analysis')
        print( '    3. In silico simulations')
        print( '    4. RAD design')
        print()
        print( '    5. Read quality')
        print( '    6. Trimming')
        print()
        print( '    7. Pseudo assembly')
        print()
        print( '    8. Read alignment')
        print()
        print( '    9. Variant calling')
        print()
        print( '    A. Pipelines')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_radseq_data_files_maintenance()
        elif option == '2':
            build_menu_radseq_enzyme_analysis()
        elif option == '3':
            build_menu_radseq_in_silico_simulations()
        elif option == '4':
            build_menu_radseq_rad_design()
        elif option == '5':
            build_menu_radseq_read_quality()
        elif option == '6':
            build_menu_radseq_trimming()
        elif option == '7':
            build_menu_radseq_pseudo_assembly()
        elif option == '8':
            build_menu_radseq_read_alignment()
        elif option == '9':
            build_menu_radseq_variant_calling()
        elif option == 'A':
            build_menu_radseq_pipelines()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_data_files_maintenance():
    '''
    Build the menu RAD-seq - Maintenance of data files.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Maintenance of data files')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Restriction site file')
        print( '    2. End file')
        print( '    3. Individual file')
        print( '    4. VCF sample file')
        print( '    5. RADdesigner condition file')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_restriction_site_file()
        elif option == '2':
            build_menu_end_file()
        elif option == '3':
            build_menu_individual_file()
        elif option == '4':
            build_menu_vcf_sample_file()
        elif option == '5':
            build_raddesinger_condition_file()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_enzyme_analysis():
    '''
    Build the menu RAD-seq - Enzyme analysis.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Enzyme analysis')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_rsitesearch_name()} ({xlib.get_ddradseqtools_name()} package)')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_rsitesearch()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_in_silico_simulations():
    '''
    Build the menu RAD-seq - In silico simulation.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - In silico simulations')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_ddradseq_simulation_name()} ({xlib.get_ddradseqtools_name()} package)')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_ddradseq_simulation()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_rad_design():
    '''
    Build the menu RAD-seq - RAD design.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - RAD design')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_raddesigner_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_raddesigner()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_read_quality():
    '''
    Build the menu RAD-seq - Read quality.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Read quality')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_fastqc_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_fastqc()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_trimming():
    '''
    Build the menu RAD-seq - Trimming.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq- Trimming')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_cutadapt_name()}')
        print(f'    2. {xlib.get_trimmomatic_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_cutadapt()
        elif option == '2':
            build_menu_trimmomatic()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_read_alignment():
    '''
    Build the menu RAD-seq - Read alignment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Read alignment')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_bowtie2_name()}')
        print(f'    2. {xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)')
        print(f'    3. {xlib.get_hisat2_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_gsnap()
        elif option == '2':
            build_menu_gsnap()
        elif option == '3':
            build_menu_hisat2()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_pseudo_assembly():
    '''
    Build the menu RAD-seq - Pseudo assembly.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Pseudo assembly')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_soapdenovo2_name()}')
        print(f'    2. {xlib.get_starcode_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_soapdenovo2()
        elif option == '2':
            build_menu_starcode()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_variant_calling():
    '''
    Build the menu RAD-seq - Variant calling.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Variant calling')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_variant_calling_name()} ({xlib.get_ddradseqtools_name()} package)')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_variant_calling()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_radseq_pipelines():
    '''
    Build the menu RAD-seq - Pipelines.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RAD-seq - Pipelines')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_ipyrad_name()}')
        print()
        print( '    X. Return to menu RAD-seq')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_ipyrad()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_bowtie2():
    '''
    Build the menu Bowtie2.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_bowtie2_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Read alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_bowtie2_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_bowtie2_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_bowtie2_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_busco():
    '''
    Build the menu BUSCO.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_busco_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly quality assessment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_busco_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_busco_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_busco_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cd_hit_est():
    '''
    Build the menu CD-HIT-EST.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_cd_hit_est_name()} ({xlib.get_cd_hit_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run transcriptome filtering process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Transcriptome filtering')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_cd_hit_est_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_cd_hit_est_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_cd_hit_est_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cuffdiff():
    '''
    Build the menu Cuffdiff.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_cuffdiff_name()} ({xlib.get_cufflinks_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run differential expression process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Differential expression')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_cuffdiff_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_cuffdiff_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_cuffdiff_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cufflinks_cuffmerge():
    '''
    Build the menu Cufflinks-Cuffmerge.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_cufflinks_cuffmerge_name()} ({xlib.get_cufflinks_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_cufflinks_cuffmerge_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_cufflinks_cuffmerge_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_cufflinks_cuffmerge_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cuffquant():
    '''
    Build the menu Cuffquant.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_cuffquant_name()} ({xlib.get_cufflinks_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run quantitation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Quantitation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_cuffquant_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_cuffquant_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_cuffquant_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_cutadapt():
    '''
    Build the menu cutadapt.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_cutadapt_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run trimming process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Trimming')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_cutadapt_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_cutadapt_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_cutadapt_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_ddradseq_simulation():
    '''
    Build the menu ddRADseqTools.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_ddradseq_simulation_name()} ({xlib.get_ddradseqtools_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run simulation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart simulation process')
        print()
        print( '    X. Return to menu In silico simulations')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_ddradseq_simulation_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_ddradseq_simulation_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_ddradseq_simulation_code())
        elif option == '4':
             cbioinfoapp.form_restart_bioinfo_process(xlib.get_ddradseq_simulation_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_express():
    '''
    Build the menu eXpress.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_express_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run quantitation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Quantitation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_express_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_express_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_express_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_fastqc():
    '''
    Build the menu FastQC.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_fastqc_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read quality process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_fastqc_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_fastqc_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_fastqc_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_ggtrinity():
    '''
    Build the menu Genome-guided Trinity.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_ggtrinity_name()} ({xlib.get_trinity_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart assembly process')
        print()
        print( '    X. Return to menu Assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_ggtrinity_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_ggtrinity_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_ggtrinity_code())
        elif option == '4':
             cbioinfoapp.form_restart_bioinfo_process(xlib.get_ggtrinity_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_gmap():
    '''
    Build the menu GMAP.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_gmap_name()} ({xlib.get_gmap_gsnap_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run transcriptome alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Transcriptome alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_gmap_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_gmap_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_gmap_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_gsnap():
    '''
    Build the menu GSNAP.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_gsnap_name()} ({xlib.get_gmap_gsnap_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Read alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_gsnap_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_gsnap_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_gsnap_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_hisat2():
    '''
    Build the menu HISAT2.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_hisat2_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Read alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_hisat2_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_hisat2_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_hisat2_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_htseq_count():
    '''
    Build the menu htseq-count.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_htseq_count_name()} ({xlib.get_htseq_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run quantitation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Quantitation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_htseq_count_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_htseq_count_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_htseq_count_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_insilico_read_normalization():
    '''
    Build the menu insilico_read_normalization (Trinity package).
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_insilico_read_normalization_name()} ({xlib.get_trinity_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run digital normalization process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart digital normalization process')
        print()
        print( '    X. Return to menu Digital normalization')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_insilico_read_normalization_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_insilico_read_normalization_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_insilico_read_normalization_code())
        elif option == '4':
             cbioinfoapp.form_restart_bioinfo_process(xlib.get_insilico_read_normalization_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_ipyrad():
    '''
    Build the menu ipyrad.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_ipyrad_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run pipeline process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Pipelines')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_ipyrad_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_ipyrad_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_ipyrad_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_kallisto():
    '''
    Build the menu kallisto.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_kallisto_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run quantitation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Quantitation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_kallisto_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_kallisto_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_kallisto_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_quast():
    '''
    Build the menu QUAST.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_quast_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly quality assessment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_quast_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_quast_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_quast_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_raddesigner():
    '''
    Build the menu RADdesigner.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RADdesigner')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run design process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart process')
        print()
        print( '    X. Return to menu RAD design')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_raddesigner_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_raddesigner_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_raddesigner_code())
        elif option == '4':
            cbioinfoapp.form_restart_bioinfo_process(xlib.get_raddesigner_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_rnaquast():
    '''
    Build the menu rnaQUAST.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_rnaquast_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly quality assessment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_rnaquast_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_rnaquast_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_rnaquast_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_rsem_eval():
    '''
    Build the menu RSEM-EVAL (reference-free evaluation of DETONATE package).
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_rsem_eval_name()} ({xlib.get_detonate_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly quality assessment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_rsem_eval_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_rsem_eval_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_rsem_eval_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_rsitesearch():
    '''
    Build the menu rsitesearch.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_rsitesearch_name()} ({xlib.get_ddradseqtools_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run enzyme analysis process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Enzyme analysis of a RAD-seq or ddRADseq experiment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_rsitesearch_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_rsitesearch_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_rsitesearch_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_soapdenovo2():
    '''
    Build the menu SOAPdenovo2.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_soapdenovo2_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run pseudo assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart assembly process')
        print()
        print( '    X. Return to menu Pseudo assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_soapdenovo2_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_soapdenovo2_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_soapdenovo2_code())
        elif option == '4':
             cbioinfoapp.form_restart_bioinfo_process(xlib.get_soapdenovo2_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_soapdenovotrans():
    '''
    Build the menu SOAPdenovo-Trans.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_soapdenovotrans_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart assembly process')
        print()
        print( '    X. Return to menu Assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_soapdenovotrans_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_soapdenovotrans_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_soapdenovotrans_code())
        elif option == '4':
             cbioinfoapp.form_restart_bioinfo_process(xlib.get_soapdenovotrans_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_star():
    '''
    Build the menu STAR.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_star_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Read alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_star_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_star_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_star_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_starcode():
    '''
    Build the menu starcode.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_starcode_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run pseudo assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Pseudo assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_starcode_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_starcode_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_starcode_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_tophat():
    '''
    Build the menu TopHat.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_tophat_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run read alignment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Read alignment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_tophat_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_tophat_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_tophat_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_transabyss():
    '''
    Build the menu Trans-ABySS.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_transabyss_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_transabyss_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_transabyss_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_transabyss_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_transcript_filter():
    '''
    Build the menu transcript-filter.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_transcript_filter_name()} ({xlib.get_ngshelper_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run transcriptome filtering process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Transcriptome filtering')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_transcript_filter_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_transcript_filter_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_transcript_filter_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_transcriptome_blastx():
    '''
    Build the menu CD-HIT-EST.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_transcriptome_blastx_name()} ({xlib.get_ngshelper_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run annotation process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Annotation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_transcriptome_blastx_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_transcriptome_blastx_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_transcriptome_blastx_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_transrate():
    '''
    Build the menu Transrate.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_transrate_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly quality assessment process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Assembly quality assessment')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_transrate_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_transrate_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_transrate_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_trimmomatic():
    '''
    Build the menu Trimmomatic.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_trimmomatic_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run trimming process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print()
        print( '    X. Return to menu Trimming')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_trimmomatic_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_trimmomatic_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_trimmomatic_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_trinity():
    '''
    Build the menu Trinity.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_trinity_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run assembly process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart assembly process')
        print()
        print( '    X. Return to menu De novo assembly')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_trinity_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_trinity_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_trinity_code())
        elif option == '4':
            cbioinfoapp.form_restart_bioinfo_process(xlib.get_trinity_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_variant_calling():
    '''
    Build the menu Variant calling process.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_variant_calling_name()} ({xlib.get_ddradseqtools_name()} package)')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run variant calling process')
        print( '       (CAUTION: before running a process, the config file should be updated)')
        print( '    4. Restart process')
        print()
        print( '    X. Return to menu Variant calling')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_variant_calling_config_file()
        elif option == '2':
            cbioinfoapp.form_edit_bioinfo_config_file(xlib.get_variant_calling_code())
        elif option == '3':
            cbioinfoapp.form_run_bioinfo_process(xlib.get_variant_calling_code())
        elif option == '4':
            cbioinfoapp.form_restart_bioinfo_process(xlib.get_variant_calling_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_end_file():
    '''
    Build the menu End file.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('End file')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate data file')
        print( '    2. Edit data file')
        print()
        print( '    X. Return to menu Maintenance of data files')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xddradseqtools.get_end_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xddradseqtools.get_end_file())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_individual_file():
    '''
    Build the menu Individual file.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Individual file')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate data file')
        print( '    2. Edit data file')
        print()
        print( '    X. Return to menu Maintenance of data files')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xddradseqtools.get_individual_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xddradseqtools.get_individual_file())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_raddesinger_condition_file():
    '''
    Build the menu RADdesigner condition file.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('RADdesigner condition file')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate data file')
        print( '    2. Edit data file')
        print()
        print( '    X. Return to menu Maintenance of data files')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xraddesigner.get_condition_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xraddesigner.get_condition_file())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_restriction_site_file():
    '''
    Build the menu Restriction site file.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Restriction site file')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate data file')
        print( '    2. Edit data file')
        print()
        print( '    X. Return to menu Maintenance of data files for simulations')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xddradseqtools.get_restriction_site_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xddradseqtools.get_restriction_site_file())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_vcf_sample_file():
    '''
    Build the menu VCF sample file.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('VCF sample file')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate data file')
        print( '    2. Edit data file')
        print()
        print( '    X. Return to menu Maintenance of data files')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xngshelper.get_vcf_sample_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xngshelper.get_vcf_sample_file())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa():
    '''
    Build the menu RAD-seq.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Taxonomy-oriented Annotation')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_toa_name()} configuration')
        print()
        print( '    2. Genomic databases')
        print()
        print( '    3. Annotation pipelines')
        print( '    4. Statistics')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_toa_configuration()
        elif option == '2':
            build_menu_toa_databases()
        elif option == '3':
            build_menu_toa_pipelines()
        elif option == '4':
            build_menu_toa_stats()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_configuration():
    '''
    Build the menu TOA configuration.
    '''


    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Configuration')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. Recreate {xlib.get_toa_name()} config file')
        print(f'    2. View {xlib.get_toa_name()} config file')
        print()
        print(f'    3. Recreate {xlib.get_toa_name()} database')
        # -- print(f'    4. Rebuild {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Taxonomy-oriented Annotation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_create_toa_config_file()
        elif option == '2':
            ctoa.form_view_toa_config_file()
        elif option == '3':
            ctoa.form_manage_toa_database(xlib.get_toa_type_recreate())
        #elif option == '4':
        #    ctoa.form_manage_toa_database(xlib.get_toa_type_rebuild())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_databases():
    '''
    Build the menu TOA - Genomic databases.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Genomic databases')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_toa_data_basic_data_name()}')
        print()
        print(f'    2. {xlib.get_toa_data_gymno_01_name()}')
        print(f'    3. {xlib.get_toa_data_dicots_04_name()}')
        print(f'    4. {xlib.get_toa_data_monocots_04_name()}')
        print()
        print(f'    5. {xlib.get_toa_data_refseq_plant_name()}')
        # -- print(f'    X. {xlib.get_toa_data_taxonomy_name()}')
        print(f'    6. {xlib.get_toa_data_nt_name()}')
        # -- print(f'    X. {xlib.get_toa_data_viridiplantae_nucleotide_gi_name()}')
        print(f'    7. {xlib.get_toa_data_nr_name()}')
        # -- print(f'    X. {xlib.get_toa_data_viridiplantae_protein_gi_name()}')
        print(f'    8. {xlib.get_toa_data_gene_name()}')
        print()
        print(f'    9. {xlib.get_toa_data_interpro_name()}')
        print()
        print(f'    A. {xlib.get_toa_data_go_name()}')
        print()
        print( '    X. Return to menu Taxonomy-oriented Annotation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_toa_basic_data()
        elif option == '2':
            build_menu_toa_gymno_01()
        elif option == '3':
            build_menu_toa_dicots_04()
        elif option == '4':
            build_menu_toa_monocots_04()
        elif option == '5':
            build_menu_toa_refseq_plant()
        # -- elif option == 'X':
        # --     build_menu_toa_taxonomy()
        elif option == '6':
            build_menu_toa_nt()
        # -- elif option == 'X':
        # --     build_menu_toa_nucleotide_gi()
        elif option == '7':
            build_menu_toa_nr()
        # -- elif option == 'X':
        # --     build_menu_toa_protein_gi()
        elif option == '8':
            build_menu_toa_gene()
        elif option == '9':
            build_menu_toa_interpro()
        elif option == 'A':
            build_menu_toa_go()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_basic_data():
    '''
    Build the menu Basic data.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_basic_data_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate genomic dataset file')
        print( '    2. Edit genomic file')
        print()
        print( '    3. Recreate species file')
        print( '    4. Edit species file')
        print()
        print( '    5. Download other basic data')
        print()
        print(f'    6. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cbioinfoapp.form_recreate_data_file(xtoa.get_dataset_file())
        elif option == '2':
            cbioinfoapp.form_edit_data_file(xtoa.get_dataset_file())
        elif option == '3':
            cbioinfoapp.form_recreate_data_file(xtoa.get_species_file())
        elif option == '4':
            cbioinfoapp.form_edit_data_file(xtoa.get_species_file())
        elif option == '5':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_basic_data_code())
        elif option == '6':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_basic_data_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_gymno_01():
    '''
    Build the menu Gymno PLAZA 1.0.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_gymno_01_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build proteome')
        print()
        print( '    2. Download functional annotations from PLAZA server')
        print(f'    3. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_proteome(), xlib.get_toa_data_gymno_01_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_gymno_01_code())
        elif option == '3':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_gymno_01_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_dicots_04():
    '''
    Build the menu Dicots PLAZA 4.0.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_dicots_04_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build proteome')
        print()
        print( '    2. Download functional annotations from PLAZA server')
        print(f'    3. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_proteome(), xlib.get_toa_data_dicots_04_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_dicots_04_code())
        elif option == '3':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_dicots_04_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_monocots_04():
    '''
    Build the menu Monocots PLAZA 4.0.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_monocots_04_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build proteome')
        print()
        print( '    2. Download functional annotations from PLAZA server')
        print(f'    3. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_proteome(), xlib.get_toa_data_monocots_04_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_monocots_04_code())
        elif option == '3':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_monocots_04_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_refseq_plant():
    '''
    Build the menu NCBI RefSeq Plant.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_refseq_plant_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build proteome')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_proteome(), xlib.get_toa_data_refseq_plant_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_taxonomy():
    '''
    Build the menu NCBI Taxonomy.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_taxonomy_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Download taxonomy data from NCBI server')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_taxonomy_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_nt():
    '''
    Build the menu NCBI BLAST database NT.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_nt_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build database for BLAST+')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_blastplus_db(), xlib.get_toa_data_nt_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_nucleotide_gi():
    '''
    Build the menu NCBI Nucleotide GenInfo identifier list.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_viridiplantae_nucleotide_gi_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build identifier list using NCBI server')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_gilist(), xlib.get_toa_data_viridiplantae_nucleotide_gi_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_nr():
    '''
    Build the menu NCBI BLAST database NR.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_nr_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build database for BLAST+')
        print( '    2. Build database for DIAMOND')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_blastplus_db(), xlib.get_toa_data_nr_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_diamond_db(), xlib.get_toa_data_nr_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_protein_gi():
    '''
    Build the menu NCBI Protein GenInfo identifier lists.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_viridiplantae_protein_gi_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Build identifier list using NCBI server')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_build_gilist(), xlib.get_toa_data_viridiplantae_protein_gi_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_gene():
    '''
    Build the menu NCBI Gene.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_gene_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Download functional annotations from NCBI server')
        print(f'    2. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_gene_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_gene_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_interpro():
    '''
    Build the menu InterPro.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_interpro_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Download functional annotations from InterPro server')
        print(f'    2. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_interpro_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_interpro_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_go():
    '''
    Build the menu Gene Ontology.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(xlib.get_toa_data_go_name())

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Download functional annotations from Gene Ontology server')
        print(f'    2. Load data into {xlib.get_toa_name()} database')
        print()
        print( '    X. Return to menu Genomic databases')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_download_data(), xlib.get_toa_data_go_code())
        elif option == '2':
            ctoa.form_manage_genomic_database(xlib.get_toa_type_load_data(), xlib.get_toa_data_go_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_pipelines():
    '''
    Build the menu TOA - Pipelines.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} - Pipelines')

        # print the menu options
        print( 'Options:')
        print()
        print(f'    1. {xlib.get_toa_name()} {xlib.get_toa_process_pipeline_nucleotide_name()}')
        print(f'    2. {xlib.get_toa_name()} {xlib.get_toa_process_pipeline_aminoacid_name()}')
        print()
        print(f'    3. Annotation merger of {xlib.get_toa_name()} pipelines')
        print()
        print( '    X. Return to menu Taxonomy-oriented Annotation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_toa_nucleotide_pipeline()
        elif option == '2':
            build_menu_toa_aminoacid_pipeline()
        elif option == '3':
            build_menu_toa_annotation_merger()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_nucleotide_pipeline():
    '''
    Build the menu TOA nucleotide pipeline.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} {xlib.get_toa_process_pipeline_nucleotide_name()}')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run pipeline')
        print( '    4. Restart pipeline')
        print()
        print(f'    X. Return to menu {xlib.get_toa_name()} pipelines')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_recreate_pipeline_config_file(xlib.get_toa_process_pipeline_nucleotide_code())
        elif option == '2':
            ctoa.form_edit_pipeline_config_file(xlib.get_toa_process_pipeline_nucleotide_code())
        elif option == '3':
            ctoa.form_run_pipeline_process(xlib.get_toa_process_pipeline_nucleotide_code())
        elif option == '4':
            ctoa.form_restart_pipeline_process(xlib.get_toa_process_pipeline_nucleotide_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_aminoacid_pipeline():
    '''
    Build the menu TOA amino acid pipeline.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'{xlib.get_toa_name()} {xlib.get_toa_process_pipeline_aminoacid_name()}')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run pipeline')
        print( '    4. Restart pipeline')
        print()
        print(f'    X. Return to menu {xlib.get_toa_name()} pipelines')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_recreate_pipeline_config_file(xlib.get_toa_process_pipeline_aminoacid_code())
        elif option == '2':
            ctoa.form_edit_pipeline_config_file(xlib.get_toa_process_pipeline_aminoacid_code())
        elif option == '3':
            ctoa.form_run_pipeline_process(xlib.get_toa_process_pipeline_aminoacid_code())
        elif option == '4':
            ctoa.form_restart_pipeline_process(xlib.get_toa_process_pipeline_aminoacid_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_annotation_merger():
    '''
    Build the menu annotation merger of TOA pipelines.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment(f'Annotation merger of {xlib.get_toa_name()} pipelines')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run process')
        print()
        print(f'    X. Return to menu {xlib.get_toa_name()} pipelines')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_recreate_annotation_merger_config_file()
        elif option == '2':
            ctoa.form_edit_pipeline_config_file(xlib.get_toa_process_merge_annotations_code())
        elif option == '3':
            ctoa.form_run_pipeline_process(xlib.get_toa_process_merge_annotations_code())
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_toa_stats():
    '''
    Build the menu Statistics.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Alignment')
        print()
        print( '    2. Annotation datasets')
        print()
        print( '    3. Species')
        print( '    4. Family')
        print( '    5. Phylum')
        print()
        print( '    6. EC')
        print( '    7. Gene Ontology')
        print( '    8. InterPro')
        print( '    9. KEGG')
        print( '    A. MapMan')
        print( '    B. MetaCyc')
        print()
        print( '    X. Return to menu Taxonomy-oriented Annotation')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            build_menu_alignment_stats()
        elif option == '2':
            build_menu_annotation_dataset_stats()
        elif option == '3':
            build_menu_species_stats()
        elif option == '4':
            build_menu_family_stats()
        elif option == '5':
            build_menu_phylum_stats()
        elif option == '6':
            build_menu_ec_stats()
        elif option == '7':
            build_menu_go_stats()
        elif option == '8':
            build_menu_interpro_stats()
        elif option == '9':
            build_menu_kegg_stats()
        elif option == 'A':
            build_menu_mapman_stats()
        elif option == 'B':
            build_menu_metacyc_stats()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_alignment_stats():
    '''
    Build the menu Statistics - Alignment.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Alignment')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. # HITs per # HSPs data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_x_per_y_data(stats_code='hit_per_hsp')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_annotation_dataset_stats():
    '''
    Build the menu Statistics - Annotation datases.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Annotation datases')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_dataset_data_frecuency()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_species_stats():
    '''
    Build the menu Statistics - Species.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Species')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_phylogenic_data_frecuency(stats_code='species')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_family_stats():
    '''
    Build the menu Statistics - Family.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Family')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_phylogenic_data_frecuency(stats_code='family')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_phylum_stats():
    '''
    Build the menu Statistics - Phylum.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Phylum')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_phylogenic_data_frecuency(stats_code='phylum')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_ec_stats():
    '''
    Build the menu Statistics - EC.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - EC')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print( '    2. # sequences per # ids data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_ontologic_data_frecuency(stats_code='ec')
        elif option == '2':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_ec')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_go_stats():
    '''
    Build the menu Statistics - Gene Ontology.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Gene Ontology')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data per term')
        print( '    2. Frecuency distribution data per namespace')
        print( '    3. # sequences per # terms data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_go_data_frecuency()
        elif option == '2':
            ctoa.form_view_phylogenic_data_frecuency(stats_code='namespace')
        elif option == '3':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_go')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_interpro_stats():
    '''
    Build the menu Statistics - InterPro.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - InterPro')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print( '    2. # sequences per # ids data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_ontologic_data_frecuency(stats_code='interpro')
        elif option == '2':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_interpro')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_kegg_stats():
    '''
    Build the menu Statistics - KEGG.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - KEGG')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print( '    2. # sequences per # ids data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_ontologic_data_frecuency(stats_code='kegg')
        elif option == '2':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_kegg')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_mapman_stats():
    '''
    Build the menu Statistics - Mapman.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - Mapman')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print( '    2. # sequences per # ids data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_ontologic_data_frecuency(stats_code='mapman')
        elif option == '2':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_mapman')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_metacyc_stats():
    '''
    Build the menu Statistics - MetaCyc.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Statistics - MetaCyc')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Frecuency distribution data')
        print( '    2. # sequences per # ids data')
        print()
        print( '    X. Return to menu Statistics')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            ctoa.form_view_ontologic_data_frecuency(stats_code='metacyc')
        elif option == '2':
            ctoa.form_view_x_per_y_data(stats_code='seq_per_metacyc')
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_datasets():
    '''
    Build the menu Datasets.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Datasets')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. List dataset (coming soon!)')
        print()
        print( '    2. Reference dataset file transfer')
        print( '    3. Reference dataset file compression/decompression')
        print( '    4. Remove reference dataset')
        print()
        print( '    5. Database file transfer')
        print( '    6. Database file compression/decompression')
        print( '    7. Remove database')
        print()
        print( '    8. Read dataset file transfer')
        print( '    9. Read dataset file compression/decompression')
        print( '    A. Remove read dataset')
        print()
        print( '    B. Result dataset file transfer')
        print( '    C. Result dataset file compression/decompression')
        print( '    D. Remove result dataset')
        print()
        print( '    E. Remove experiment')
        print()
        print( '    X. Return to menu Main')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            pass
        elif option == '2':
            build_menu_reference_file_transfer()
        elif option == '3':
            build_menu_reference_file_compression_decompression()
        elif option == '4':
            cdataset.form_remove_reference_dataset()
        elif option == '5':
            build_menu_database_file_transfer()
        elif option == '6':
            build_menu_database_file_compression_decompression()
        elif option == '7':
            cdataset.form_remove_database_dataset()
        elif option == '8':
            build_menu_read_file_transfer()
        elif option == '9':
            build_menu_read_file_compression_decompression()
        elif option == 'A':
            cdataset.form_remove_read_dataset()
        elif option == 'B':
            build_menu_result_file_transfer()
        elif option == 'C':
            build_menu_result_file_compression_decompression()
        elif option == 'D':
            cdataset.form_remove_result_dataset()
        elif option == 'E':
            cdataset.form_remove_experiment()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_file_transfer():
    '''
    Build the menu Reference dataset file transfer.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference dataset file transfer')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Upload dataset to a cluster')
        print( '       (CAUTION: before running a upload process, the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_reference_transfer_config_file()
        elif option == '2':
            cdataset.form_edit_reference_transfer_config_file()
        elif option == '3':
            cdataset.form_upload_reference_dataset()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_reference_file_compression_decompression():
    '''
    Build the menu Reference dataset file compression/decompression.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Reference dataset file compression/decompression')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run compression/decompression process')
        print( '       (CAUTION: before running a compression/decompression process,')
        print( '                 the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_reference_gzip_config_file()
        elif option == '2':
            cdataset.form_edit_reference_gzip_config_file()
        elif option == '3':
            cdataset.form_run_reference_gzip_process()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_database_file_transfer():
    '''
    Build the menu Database file transfer.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Database file transfer')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Upload dataset to a cluster')
        print( '       (CAUTION: before running a upload process, the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_database_transfer_config_file()
        elif option == '2':
            cdataset.form_edit_database_transfer_config_file()
        elif option == '3':
            cdataset.form_upload_database_dataset()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_database_file_compression_decompression():
    '''
    Build the menu Database file compression/decompression.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Database file compression/decompression')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run compression/decompression process')
        print( '       (CAUTION: before running a compression/decompression process,')
        print( '                 the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_database_gzip_config_file()
        elif option == '2':
            cdataset.form_edit_database_gzip_config_file()
        elif option == '3':
            cdataset.form_run_database_gzip_process()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_read_file_transfer():
    '''
    Build the menu Read dataset file transfer.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Read dataset file transfer')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Upload dataset to a cluster')
        print( '       (CAUTION: before running a upload process, the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_read_transfer_config_file()
        elif option == '2':
            cdataset.form_edit_read_transfer_config_file()
        elif option == '3':
            cdataset.form_upload_read_dataset()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_read_file_compression_decompression():
    '''
    Build the menu Read dataset file compression/decompression.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Read dataset file compression/decompression')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run compression/decompression process')
        print( '       (CAUTION: before running a compression/decompression process,')
        print( '                 the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_read_gzip_config_file()
        elif option == '2':
            cdataset.form_edit_read_gzip_config_file()
        elif option == '3':
            cdataset.form_run_read_gzip_process()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_result_file_transfer():
    '''
    Build the menu Result dataset file transfer.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Result dataset file transfer')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Download dataset from a cluster')
        print( '       (CAUTION: before running a download process, the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_result_transfer_config_file()
        elif option == '2':
            cdataset.form_edit_result_transfer_config_file()
        elif option == '3':
            cdataset.form_download_result_dataset()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_result_file_compression_decompression():
    '''
    Build the menu Result dataset file compression/decompression.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Result dataset file compression/decompression')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. Recreate config file')
        print( '    2. Edit config file')
        print()
        print( '    3. Run compression/decompression process')
        print( '       (CAUTION: before running a compression/decompression process,')
        print( '                 the corresponding config file should be updated)')
        print()
        print( '    X. Return to menu Datasets')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            cdataset.form_recreate_result_gzip_config_file()
        elif option == '2':
            cdataset.form_edit_result_gzip_config_file()
        elif option == '3':
            cdataset.form_run_result_gzip_process()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

def build_menu_logs():
    '''
    Build the menu Logs.
    '''

    while True:

        # print headers
        clib.clear_screen()
        clib.print_headers_with_environment('Cluster logs')

        # print the menu options
        print( 'Options:')
        print()
        print( '    1. View the cluster start log')
        print()
        print( '    2. List submission logs in the local computer')
        print( '    3. View a submission log in the local computer')
        print()
        print( '    4. List result logs in the cluster')
        print( '    5. View a result log in the cluster')
        print()
        print( '    X. Return to menu Logs')
        print()

        # get the selected option
        option = input('Input the selected option: ').upper()

        # process the selected option
        if option == '1':
            clog.form_view_cluster_start_log()
        elif option == '2':
            clog.form_list_submission_logs()
        elif option == '3':
            clog.form_view_submission_log()
        elif option == '4':
            clog.form_list_cluster_experiment_processes()
        elif option == '5':
            clog.form_view_cluster_experiment_process_log()
        elif option == 'X':
            break

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print( 'This file contains the functions related to menus in console mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
