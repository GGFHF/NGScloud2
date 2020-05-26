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
This file contains the classes related to forms corresponding to Cloud Control
menu items in gui mode.
'''

#-------------------------------------------------------------------------------

import os
import PIL.Image
import PIL.ImageTk
import sys
import threading
import tkinter
import tkinter.filedialog
import tkinter.ttk

import gdialogs
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
import xrnaquast
import xsoapdenovo2
import xsoapdenovotrans
import xssh
import xstar
import xstarcode
import xread
import xreference
import xresult
import xtoa
import xtophat
import xtransabyss
import xtransrate
import xtrimmomatic
import xtrinity
import xvolume

#-------------------------------------------------------------------------------

class FormSetEnvironment(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormSetEnvironment" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Set environment'

        # create the wrappers to track changes in inputs
        self.wrapper_environment = tkinter.StringVar()
        self.wrapper_environment.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormSetEnvironment".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_instructions" and register it with the grid geometry manager
        self.label_instructions = tkinter.Label(self, text='Select an environment or type a new one:')
        self.label_instructions.grid(row=0, column=0, columnspan=2, padx=(15,5), pady=(75,5), sticky='w')

        # create "label_environment" and register it with the grid geometry manager
        self.label_environment = tkinter.Label(self, text='Environment')
        self.label_environment.grid(row=1, column=0, padx=(15,5), pady=(25,5), sticky='e')

        # create "combobox_environment" and register it with the grid geometry manager
        self.combobox_environment = tkinter.ttk.Combobox(self, width=40, height=4, state='normal', textvariable=self.wrapper_environment)
        self.combobox_environment.grid(row=1, column=1, padx=(5,5), pady=(25,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(55+xlib.get_os_size_fix()))
        self.label_fit.grid(row=2, column=2, padx=(0,0), pady=(25,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_exit" and register it with the grid geometry manager
        self.button_exit = tkinter.ttk.Button(self, text='Exit', command=self.exit)
        self.button_exit.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_environment.bind('<<ComboboxSelected>>', self.combobox_environment_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_environment()

    #---------------

    def populate_combobox_environment(self):
        '''
        Populate data in "combobox_environment".
        '''

        # clear the value selected in the combobox
        #self.wrapper_environment.set('')

        # check if there are any running clusters
        if xconfiguration.get_environments_list() == []:
            message = 'There is not any environment recorded. You have to type a new one.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_environment['values'] = xconfiguration.get_environments_list()

    #---------------

    def combobox_environment_selected_item(self, event=None):
        '''
        Process the event when an environment item has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormSetEnvironment" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_continue" has to be enabled or disabled
        if self.wrapper_environment.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Set the environment.
        '''

        # if the environment is new, check it is right
        if self.wrapper_environment.get() not in xconfiguration.get_environments_list():
            message = f'{self.wrapper_environment.get()} is not an environment recorded. Do you like to record it?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - Exit', message)
            if not OK:
                return

        # set the current environment
        xconfiguration.environment = self.wrapper_environment.get()

        # record the curren environment if it is not recorded
        if xconfiguration.environment not in xconfiguration.get_environments_list():
            (OK, error_list) = xconfiguration.add_environment(xconfiguration.environment)
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                raise xlib.ProgramException('C002')

        # check if it is necesary to create the NGScloud config file corresponding to the environment
        if not xconfiguration.is_ngscloud_config_file_created():

            # clear the label of the current process name
            self.main.label_process['text'] = ''

            # close the current form
            self.main.close_current_form()

            # create and register "form_create_config_files" in "container" of "Main with the grid geometry manager
            form_create_config_files = FormCreateConfigFiles(self.main.container, self.main)
            form_create_config_files.grid(row=0, column=0, sticky='nsew')

            # set "form_create_config_files" as current form and add it in the forms dictionary
            self.main.current_form = 'form_create_config_files'
            self.main.forms_dict[self.main.current_form] = form_create_config_files

            # raise "form_create_config_files" to front
            form_create_config_files.tkraise()

        # in case of the config file of the environment is created
        else:

            # set the environment variables corresponding to the NGScloud config file, the AWS access key identification, AWS secret access key and the current region name
            xconfiguration.set_environment_variables()

            # set the current region, zone and environment in the corresponding widgets of "frame_information"
            self.main.label_region_value['text'] = xconfiguration.get_current_region_name()
            self.main.label_zone_value['text'] = xconfiguration.get_current_zone_name()
            self.main.label_environment_value['text'] = xconfiguration.environment

            # clear the label of the current process name
            self.main.label_process['text'] = ''

            # build the full menu of the application
            self.main.build_full_menu()

            # close the current form
            self.main.close_current_form()

    #---------------

    def exit(self, event=None):
        '''
        Exit the application.
        '''

        # initialize the control variable
        OK = True

        # confirm the exit of NGScloud
        message = f'Are you sure to exit {xlib.get_project_name()}?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - Exit', message)

        # exit NGScloud
        if OK:
            # destroy "FormStart"
            self.destroy()
            # destroy "Main"
            self.main.destroy()

    #---------------

#-------------------------------------------------------------------------------

class FormCreateConfigFiles(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormCreateConfigFiles" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Create config files'

        # create the wrappers to track changes in inputs
        self.wrapper_user_id = tkinter.StringVar()
        self.wrapper_user_id.trace('w', self.check_inputs)
        self.wrapper_access_key_id = tkinter.StringVar()
        self.wrapper_access_key_id.trace('w', self.check_inputs)
        self.wrapper_secret_access_key = tkinter.StringVar()
        self.wrapper_secret_access_key.trace('w', self.check_inputs)
        self.wrapper_email = tkinter.StringVar()
        self.wrapper_email.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormCreateConfigFiles".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_user_id" and register it with the grid geometry manager
        self.label_user_id = tkinter.Label(self, text='User id')
        self.label_user_id.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "entry_user_id" and register it with the grid geometry manager
        self.entry_user_id = tkinter.Entry(self, textvariable=self.wrapper_user_id, width=15, validatecommand=self.check_inputs)
        self.entry_user_id.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_access_key_id" and register it with the grid geometry manager
        self.label_access_key_id = tkinter.Label(self, text='Access key id')
        self.label_access_key_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_access_key_id" and register it with the grid geometry manager
        self.entry_access_key_id = tkinter.Entry(self, textvariable=self.wrapper_access_key_id, width=25, validatecommand=self.check_inputs)
        self.entry_access_key_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_secret_access_key" and register it with the grid geometry manager
        self.label_secret_access_key = tkinter.Label(self, text='Secret access key')
        self.label_secret_access_key.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_secret_access_key" and register it with the grid geometry manager
        self.entry_secret_access_key = tkinter.Entry(self, textvariable=self.wrapper_secret_access_key, width=50, validatecommand=self.check_inputs)
        self.entry_secret_access_key.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email" and register it with the grid geometry manager
        self.label_email = tkinter.Label(self, text='Contact e-mail')
        self.label_email.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_email" and register it with the grid geometry manager
        self.entry_email = tkinter.Entry(self, textvariable=self.wrapper_email, width=50, validatecommand=self.check_inputs)
        self.entry_email.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email_warning" and register it with the grid geometry manager
        self.label_email_warning = tkinter.Label(self, text='')
        self.label_email_warning.grid(row=4, column=1, padx=(5,5), pady=(5,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(1+xlib.get_os_size_fix()))
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(25,5), sticky='e')

        # create "button_create" and register it with the grid geometry manager
        self.button_create = tkinter.ttk.Button(self, text='Create', command=self.create, state='disabled')
        self.button_create.grid(row=5, column=3, padx=(5,5), pady=(25,5), sticky='e')

        # create "button_exit" and register it with the grid geometry manager
        self.button_exit = tkinter.ttk.Button(self, text='Exit', command=self.exit)
        self.button_exit.grid(row=5, column=4, padx=(5,5), pady=(25,5), sticky='w')

    #---------------

    def create(self):
        '''
        Create the config files.
        '''

        # initialize the control variable
        OK = True

        # check the AWS access key identification and the AWS secret access key   
        OK = xec2.check_aws_credentials(self.wrapper_access_key_id.get(), self.wrapper_secret_access_key.get())
        if not OK:
            message = 'ERROR: The credentials are wrong. Please review your access key identification and secret access key in the AWS web.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # create the NGScloud config file corresponding to the environment
        if OK:
            (OK, error_list) = xconfiguration.create_ngscloud_config_file(self.wrapper_user_id.get(), self.wrapper_access_key_id.get(), self.wrapper_secret_access_key.get(), self.wrapper_email.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                raise xlib.ProgramException('C001')

        # create the key pairs directory
        if OK:
            if not os.path.exists(xlib.get_keypairs_dir()):
                os.makedirs(xlib.get_keypairs_dir())

        # create the config files
        if OK:

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

            # create the TopHat config file
            (OK, error_list) = xtophat.create_tophat_config_file()

            # create the Trans-ABySS config file
            (OK, error_list) = xtransabyss.create_transabyss_config_file()

            # create the starcode config file
            (OK, error_list) = xstarcode.create_starcode_config_file()

            # create the transcript-filter config file
            (OK, error_list) = xngshelper.create_transcriptome_blastx_config_file()

            # create the transcriptome-blastx config file
            (OK, error_list) = xngshelper.create_transcript_filter_config_file()

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

        # show the creation message
        if OK:
            message = 'The config files are created with default values.'
            tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - Start', message)

        # set the environment variables corresponding to the NGScloud config file
        if OK:
            xconfiguration.set_environment_variables()

        # set the current region, zone and environment in the corresponding widgets of "frame_information"
        if OK:
            self.main.label_region_value['text'] = xconfiguration.get_current_region_name()
            self.main.label_zone_value['text'] = xconfiguration.get_current_zone_name()
            self.main.label_environment_value['text'] = xconfiguration.environment

        # clear the label of the current process name
        if OK:
            self.main.label_process['text'] = ''

        # build the full menu of the application
        if OK:
            self.main.build_full_menu()

        # close the current form
        if OK:
            self.main.close_current_form()

    #---------------

    def exit(self, event=None):
        '''
        Exit the application.
        '''

        # initialize the control variable
        OK = True

        # confirm the exit of NGScloud
        message = f'Are you sure to exit {xlib.get_project_name()}?'
        OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - Exit', message)

        # exit NGScloud
        if OK:
            # destroy "FormStart"
            self.destroy()
            # destroy "Main"
            self.main.destroy()

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormCreateConfigFiles" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_email"
        if not self.check_entry_email():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_user_id.get() != '' and self.wrapper_access_key_id.get() != '' and self.wrapper_secret_access_key.get() != '' and self.wrapper_email.get() != '':
            self.button_create['state'] = 'enable'
        else:
            self.button_create['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_email(self):
        '''
        Check the content of "entry_email"
        '''

        # initialize the control variable
        OK = True

        # check if the e-mail address in "entry_email" is valid
        if not xlib.is_email_address_valid(self.wrapper_email.get()):
            self.label_email_warning['text'] = 'The e-mail address is not OK.'
            self.label_email_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_email_warning['text'] = ''
            self.label_email_warning['foreground'] = 'black'

        # return the control variable
        return OK

   #---------------

#-------------------------------------------------------------------------------

class FormRecreateNGScloudConfigFile(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormRecreateNGScloudConfigFile" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = f'Configuration - Recreate {xlib.get_project_name()} config file'

        # create the wrappers to track changes in inputs
        self.wrapper_user_id = tkinter.StringVar()
        self.wrapper_user_id.trace('w', self.check_inputs)
        self.wrapper_access_key_id = tkinter.StringVar()
        self.wrapper_access_key_id.trace('w', self.check_inputs)
        self.wrapper_secret_access_key = tkinter.StringVar()
        self.wrapper_secret_access_key.trace('w', self.check_inputs)
        self.wrapper_email = tkinter.StringVar()
        self.wrapper_email.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial value to inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRecreateNGScloudConfigFile".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_user_id" and register it with the grid geometry manager
        self.label_user_id = tkinter.Label(self, text='User id')
        self.label_user_id.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "entry_user_id" and register it with the grid geometry manager
        self.entry_user_id = tkinter.Entry(self, textvariable=self.wrapper_user_id, width=15, validatecommand=self.check_inputs)
        self.entry_user_id.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_access_key_id" and register it with the grid geometry manager
        self.label_access_key_id = tkinter.Label(self, text='Access key id')
        self.label_access_key_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_access_key_id" and register it with the grid geometry manager
        self.entry_access_key_id = tkinter.Entry(self, textvariable=self.wrapper_access_key_id, width=25, validatecommand=self.check_inputs)
        self.entry_access_key_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_secret_access_key" and register it with the grid geometry manager
        self.label_secret_access_key = tkinter.Label(self, text='Secret access key')
        self.label_secret_access_key.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_secret_access_key" and register it with the grid geometry manager
        self.entry_secret_access_key = tkinter.Entry(self, textvariable=self.wrapper_secret_access_key, width=50, validatecommand=self.check_inputs)
        self.entry_secret_access_key.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email" and register it with the grid geometry manager
        self.label_email = tkinter.Label(self, text='Contact e-mail')
        self.label_email.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_email" and register it with the grid geometry manager
        self.entry_email = tkinter.Entry(self, textvariable=self.wrapper_email, width=50, validatecommand=self.check_inputs)
        self.entry_email.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email_warning" and register it with the grid geometry manager
        self.label_email_warning = tkinter.Label(self, text='')
        self.label_email_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(29+xlib.get_os_size_fix()))
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # get basic AWS data and contact e-mail address from NGScloud config file
        (user_id, access_key_id, secret_access_key) = xconfiguration.get_basic_aws_data()
        email = xconfiguration.get_contact_data()

        # load initial data in inputs
        self.wrapper_user_id.set(user_id)
        self.wrapper_access_key_id.set(access_key_id)
        self.wrapper_secret_access_key.set(secret_access_key)
        self.wrapper_email.set(email)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRecreateNGScloudConfigFile" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_email"
        if not self.check_entry_email():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_user_id.get() != '' and self.wrapper_access_key_id.get() != '' and self.wrapper_secret_access_key.get() != '' and self.wrapper_email.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_email(self):
        '''
        Check the content of "entry_email"
        '''

        # initialize the control variable
        OK = True

        # check if the e-mail address in "entry_email" is valid
        if not xlib.is_email_address_valid(self.wrapper_email.get()):
            self.label_email_warning['text'] = 'The e-mail address is not OK.'
            self.label_email_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_email_warning['text'] = ''
            self.label_email_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the recreation of the NGScloud config file.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the recreation of the NGScloud config file
        if OK:
            message = f'The file {xconfiguration.get_ngscloud_config_file()} is going to be recreated. The previous file will be lost.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # recreate the NGScloud config file corresponding to the environment
        if OK:
            (OK, error_list) = xconfiguration.create_ngscloud_config_file(self.wrapper_user_id.get(), self.wrapper_access_key_id.get(), self.wrapper_secret_access_key.get(), self.wrapper_email.get())
            if OK:
                message = f'The file {xconfiguration.get_ngscloud_config_file()} is created with default values.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False

        # set the current region and zone in the corresponding widgets of "frame_information"
        if OK:
            self.main.label_region_value['text'] = xconfiguration.get_current_region_name()
            self.main.label_zone_value['text'] = xconfiguration.get_current_zone_name()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormRecreateNGScloudConfigFile".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

class FormUpdateConnectionData(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormUpdateConnectionData" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Configuration - Update connection data and contact e-mail'

        # create the wrappers to track changes in inputs
        self.wrapper_user_id = tkinter.StringVar()
        self.wrapper_user_id.trace('w', self.check_inputs)
        self.wrapper_access_key_id = tkinter.StringVar()
        self.wrapper_access_key_id.trace('w', self.check_inputs)
        self.wrapper_secret_access_key = tkinter.StringVar()
        self.wrapper_secret_access_key.trace('w', self.check_inputs)
        self.wrapper_email = tkinter.StringVar()
        self.wrapper_email.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial value to inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormUpdateConnectionData".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_user_id" and register it with the grid geometry manager
        self.label_user_id = tkinter.Label(self, text='User id')
        self.label_user_id.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "entry_user_id" and register it with the grid geometry manager
        self.entry_user_id = tkinter.Entry(self, textvariable=self.wrapper_user_id, width=15, validatecommand=self.check_inputs)
        self.entry_user_id.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_access_key_id" and register it with the grid geometry manager
        self.label_access_key_id = tkinter.Label(self, text='Access key id')
        self.label_access_key_id.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_access_key_id" and register it with the grid geometry manager
        self.entry_access_key_id = tkinter.Entry(self, textvariable=self.wrapper_access_key_id, width=25, validatecommand=self.check_inputs)
        self.entry_access_key_id.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_secret_access_key" and register it with the grid geometry manager
        self.label_secret_access_key = tkinter.Label(self, text='Secret access key')
        self.label_secret_access_key.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_secret_access_key" and register it with the grid geometry manager
        self.entry_secret_access_key = tkinter.Entry(self, textvariable=self.wrapper_secret_access_key, width=50, validatecommand=self.check_inputs)
        self.entry_secret_access_key.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email" and register it with the grid geometry manager
        self.label_email = tkinter.Label(self, text='Contact e-mail')
        self.label_email.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_email" and register it with the grid geometry manager
        self.entry_email = tkinter.Entry(self, textvariable=self.wrapper_email, width=50, validatecommand=self.check_inputs)
        self.entry_email.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_email_warning" and register it with the grid geometry manager
        self.label_email_warning = tkinter.Label(self, text='')
        self.label_email_warning.grid(row=3, column=2,columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(29+xlib.get_os_size_fix()))
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''
        # get basic AWS data and contact e-mail address from NGScloud config file
        (user_id, access_key_id, secret_access_key) = xconfiguration.get_basic_aws_data()
        email = xconfiguration.get_contact_data()

        # load initial data in inputs
        self.wrapper_user_id.set(user_id)
        self.wrapper_access_key_id.set(access_key_id)
        self.wrapper_secret_access_key.set(secret_access_key)
        self.wrapper_email.set(email)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormUpdateConnectionData" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_email"
        if not self.check_entry_email():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_user_id.get() != '' and self.wrapper_access_key_id.get() != '' and self.wrapper_secret_access_key.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_email(self):
        '''
        Check the content of "entry_email"
        '''

        # initialize the control variable
        OK = True

        # check if the e-mail address in "entry_email" is valid
        if not xlib.is_email_address_valid(self.wrapper_email.get()):
            self.label_email_warning['text'] = 'The e-mail address is not OK.'
            self.label_email_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_email_warning['text'] = ''
            self.label_email_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the update of the connection data in the NGScloud config file.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check the AWS access key identification and the AWS secret access key
        if OK:
            OK = xec2.check_aws_credentials(self.wrapper_access_key_id.get(), self.wrapper_secret_access_key.get())
            if not OK:
                message = 'ERROR: The credentials are wrong. Please review your access key identification and secret access key in the AWS web.'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the NGScloud config file
        if OK:
            ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

        # confirm the connection data update in the NGScloud config file
        if OK:
            message = f'The file {ngscloud_config_file} is going to be update with the new connection data.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # update the connection data in the NGScloud config file corresponding to the environment
        if OK:
            (OK, error_list) = xconfiguration.update_connection_data(self.wrapper_user_id.get(), self.wrapper_access_key_id.get(), self.wrapper_secret_access_key.get())
            if OK:
                (OK, error_list) = xconfiguration.update_contact_data(self.wrapper_email.get())
            if OK:
                message = f'The file {ngscloud_config_file} has been updated.'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                raise xlib.ProgramException('C001')

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormUpdateConnectionData".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

class FormUpdateRegionZone(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormUpdateRegionZone" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Configuration - Update region and zone'

        # set the region name and zone name
        self.region_name = None
        self.zone_name = None

        # create the wrappers to track changes in inputs
        self.wrapper_region_name = tkinter.StringVar()
        self.wrapper_region_name.trace('w', self.check_inputs)
        self.wrapper_zone_name = tkinter.StringVar()
        self.wrapper_zone_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormUpdateRegionZone".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_region_name" and register it with the grid geometry manager
        self.label_region_name = tkinter.Label(self, text='Region name')
        self.label_region_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_region_name" and register it with the grid geometry manager
        self.combobox_region_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_region_name)
        self.combobox_region_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_zone_name" and register it with the grid geometry manager
        self.label_zone_name = tkinter.Label(self, text='Zone name')
        self.label_zone_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_zone_name" and register it with the grid geometry manager
        self.combobox_zone_name = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_zone_name)
        self.combobox_zone_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(94+xlib.get_os_size_fix()))
        self.label_fit.grid(row=2, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_region_name.bind('<<ComboboxSelected>>', self.combobox_region_name_selected_item)
        self.combobox_zone_name.bind('<<ComboboxSelected>>', self.combobox_zone_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_region_name()

    #---------------

    def populate_combobox_region_name(self):
        '''
        Populate data in "combobox_region_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_region_name.set('')

        # load the region names list in the combobox
        self.combobox_region_name['values'] = xec2.get_available_region_list()

    #---------------

    def populate_combobox_zone_name(self):
        '''
        Populate data in "combobox_zone_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_zone_name.set('')

        # load the region name list in the combobox
        self.combobox_zone_name['values'] = xec2.get_available_zone_list(self.wrapper_region_name.get())

    #---------------

    def combobox_region_name_selected_item(self, parent, event=None):
        '''
        Process the event when a region name item has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # load data in "combobox_zone_name"
        self.populate_combobox_zone_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_zone_name_selected_item(self, parent, event=None):
        '''
        Process the event when a region name item has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormUpdateRegionZone" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_region_name.get() != '' and self.wrapper_zone_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the update of the region and zone in the NGScloud config file.
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {self.head}', message)

        # get the NGScloud config file
        if OK:
            ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

        # confirm the region and zone update in the NGScloud config file
        if OK:
            message = f'The file {ngscloud_config_file} is going to be update with the new the new region and zone.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # update the connection data in the NGScloud config file corresponding to the environment
        if OK:
            (OK, error_list) = xconfiguration.update_region_zone_data(self.wrapper_region_name.get(), self.wrapper_zone_name.get())
            if OK:
                if error_list == []:
                    message = f'The file {ngscloud_config_file} is updated.'
                else:
                    message = ''
                    for error in error_list:
                        message = f'{message}{error}\n'
                tkinter.messagebox.showinfo(f'{xlib.get_project_name()} - {self.head}', message)
            else:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {self.head}', message)

        # set the current region and zone in the corresponding widgets of "frame_information"
        if OK:
            self.main.label_region_value['text'] = xconfiguration.get_current_region_name()
            self.main.label_zone_value['text'] = xconfiguration.get_current_zone_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormUpdateRegionZone".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

class FormLinkVolumes(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormLinkVolumes" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Configuration - Link volumes'

        # get current zone name
        self.zone_name = xconfiguration.get_current_zone_name()

        # get created volume name list in the zone
        self.created_volume_name_list = xec2.get_created_volume_name_list(self.zone_name)

        # create the wrappers to track changes in inputs
        self.wrapper_dataset_structure = tkinter.StringVar()
        self.wrapper_dataset_structure.trace('w', self.check_inputs)
        self.wrapper_ngscloud_volume = tkinter.StringVar()
        self.wrapper_ngscloud_volume.trace('w', self.check_inputs)
        self.wrapper_app_volume = tkinter.StringVar()
        self.wrapper_app_volume.trace('w', self.check_inputs)
        self.wrapper_database_volume = tkinter.StringVar()
        self.wrapper_database_volume.trace('w', self.check_inputs)
        self.wrapper_read_volume = tkinter.StringVar()
        self.wrapper_read_volume.trace('w', self.check_inputs)
        self.wrapper_reference_volume = tkinter.StringVar()
        self.wrapper_reference_volume.trace('w', self.check_inputs)
        self.wrapper_result_volume = tkinter.StringVar()
        self.wrapper_result_volume.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # check there are created volumes in the zone
        if self.created_volume_name_list == []:
            message = f'There is not any volume created in the zone {self.zone_name}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check there are created volumes in the zone
        if self.created_volume_name_list == []:
            message = f'There is not any volume created in the zone {self.zone_name}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormLinkVolumes".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_dataset_structure" and register it with the grid geometry manager
        self.label_dataset_structure = tkinter.Label(self, text='Dataset structure')
        self.label_dataset_structure.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_dataset_structure" and register it with the grid geometry manager
        self.combobox_dataset_structure = tkinter.ttk.Combobox(self, width=20, height=4, state='readonly', textvariable=self.wrapper_dataset_structure)
        self.combobox_dataset_structure.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_ngscloud_volume" and register it with the grid geometry manager
        self.label_ngscloud_volume = tkinter.Label(self, text='NGScloud volume')
        self.label_ngscloud_volume.grid(row=1, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_ngscloud_volume" and register it with the grid geometry manager
        self.combobox_ngscloud_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_ngscloud_volume)
        self.combobox_ngscloud_volume.grid(row=1, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_app_volume" and register it with the grid geometry manager
        self.label_app_volume = tkinter.Label(self, text='Application volume')
        self.label_app_volume.grid(row=2, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_app_volume" and register it with the grid geometry manager
        self.combobox_app_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_app_volume)
        self.combobox_app_volume.grid(row=2, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_database_volume" and register it with the grid geometry manager
        self.label_database_volume = tkinter.Label(self, text='Database volume')
        self.label_database_volume.grid(row=3, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_database_volume" and register it with the grid geometry manager
        self.combobox_database_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_database_volume)
        self.combobox_database_volume.grid(row=3, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_read_volume" and register it with the grid geometry manager
        self.label_read_volume = tkinter.Label(self, text='Read volume')
        self.label_read_volume.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_read_volume" and register it with the grid geometry manager
        self.combobox_read_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_read_volume)
        self.combobox_read_volume.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_reference_volume" and register it with the grid geometry manager
        self.label_reference_volume = tkinter.Label(self, text='Reference volume')
        self.label_reference_volume.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_reference_volume" and register it with the grid geometry manager
        self.combobox_reference_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_reference_volume)
        self.combobox_reference_volume.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_result_volume" and register it with the grid geometry manager
        self.label_result_volume = tkinter.Label(self, text='Result volume')
        self.label_result_volume.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_result_volume" and register it with the grid geometry manager
        self.combobox_result_volume = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_result_volume)
        self.combobox_result_volume.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(40+xlib.get_os_size_fix()))
        self.label_fit.grid(row=7, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=7, column=3, padx=(5,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=7, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_dataset_structure.bind('<<ComboboxSelected>>', self.combobox_dataset_structure_selected_item)
        self.combobox_ngscloud_volume.bind('<<ComboboxSelected>>', self.combobox_ngscloud_volume_selected_item)
        self.combobox_app_volume.bind('<<ComboboxSelected>>', self.combobox_app_volume_selected_item)
        self.combobox_database_volume.bind('<<ComboboxSelected>>', self.combobox_database_volume_selected_item)
        self.combobox_read_volume.bind('<<ComboboxSelected>>', self.combobox_read_volume_selected_item)
        self.combobox_reference_volume.bind('<<ComboboxSelected>>', self.combobox_reference_volume_selected_item)
        self.combobox_result_volume.bind('<<ComboboxSelected>>', self.combobox_result_volume_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_dataset_structure()
        self.populate_combobox_ngscloud_volume()
        self.populate_combobox_app_volume()
        self.populate_combobox_database_volume()
        self.populate_combobox_read_volume()
        self.populate_combobox_reference_volume()
        self.populate_combobox_result_volume()

        # load initial data in inputs
        self.wrapper_dataset_structure.set(xconfiguration.get_dataset_structure_singlevolume())
        self.wrapper_ngscloud_volume.set('')
        self.wrapper_app_volume.set('')
        self.combobox_app_volume['state'] = 'disabled'
        self.wrapper_database_volume.set('')
        self.combobox_database_volume['state'] = 'disabled'
        self.wrapper_read_volume.set('')
        self.combobox_read_volume['state'] = 'disabled'
        self.wrapper_reference_volume.set('')
        self.combobox_reference_volume['state'] = 'disabled'
        self.wrapper_result_volume.set('')
        self.combobox_result_volume['state'] = 'disabled'

    #---------------

    def populate_combobox_dataset_structure(self):
        '''
        Populate data in "combobox_dataset_structure".
        '''

        # clear the value selected in the combobox
        self.wrapper_dataset_structure.set('')

        # load the mounting points list in the combobox
        self.combobox_dataset_structure['values'] = xconfiguration.get_dataset_structure_list()

    #---------------

    def populate_combobox_ngscloud_volume(self):
        '''
        Populate data in "combobox_ngscloud_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_ngscloud_volume.set('')

        # load the volume names list in the combobox
        self.combobox_ngscloud_volume['values'] = self.created_volume_name_list

    #---------------

    def populate_combobox_app_volume(self):
        '''
        Populate data in "combobox_app_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_app_volume.set('')

        # load the mounting points list in the combobox
        self.combobox_app_volume['values'] = self.created_volume_name_list

    #---------------

    def populate_combobox_database_volume(self):
        '''
        Populate data in "combobox_database_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_database_volume.set('')

        # load the mounting points list in the combobox
        self.combobox_database_volume['values'] = ['NONE'] + self.created_volume_name_list

    #---------------

    def populate_combobox_read_volume(self):
        '''
        Populate data in "combobox_read_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_read_volume.set('')

        # load the mounting points list in the combobox
        self.combobox_read_volume['values'] = ['NONE'] + self.created_volume_name_list

    #---------------

    def populate_combobox_reference_volume(self):
        '''
        Populate data in "combobox_reference_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_reference_volume.set('')

        # load the mounting points list in the combobox
        self.combobox_reference_volume['values'] = ['NONE'] + self.created_volume_name_list

    #---------------

    def populate_combobox_result_volume(self):
        '''
        Populate data in "combobox_result_volume".
        '''

        # clear the value selected in the combobox
        self.wrapper_result_volume.set('')

        # load the mounting points list in the combobox
        self.combobox_result_volume['values'] = ['NONE'] + self.created_volume_name_list

    #---------------

    def combobox_dataset_structure_selected_item(self, event=None):
        '''
        Process the event when a dataset structure item has been selected
        '''

        if self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_none():
            self.combobox_ngscloud_volume['state'] = 'readonly'
            self.wrapper_ngscloud_volume.set('')
            self.combobox_ngscloud_volume['state'] = 'disabled'
            self.combobox_app_volume['state'] = 'readonly'
            self.wrapper_app_volume.set('')
            self.combobox_app_volume['state'] = 'disabled'
            self.combobox_database_volume['state'] = 'readonly'
            self.wrapper_database_volume.set('')
            self.combobox_database_volume['state'] = 'disabled'
            self.combobox_read_volume['state'] = 'readonly'
            self.wrapper_read_volume.set('')
            self.combobox_read_volume['state'] = 'disabled'
            self.combobox_reference_volume['state'] = 'readonly'
            self.wrapper_reference_volume.set('')
            self.combobox_reference_volume['state'] = 'disabled'
            self.combobox_result_volume['state'] = 'readonly'
            self.wrapper_result_volume.set('')
            self.combobox_result_volume['state'] = 'disabled'
        elif self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_singlevolume():
            self.combobox_ngscloud_volume['state'] = 'readonly'
            self.wrapper_ngscloud_volume.set('')
            self.combobox_app_volume['state'] = 'readonly'
            self.wrapper_app_volume.set('')
            self.combobox_app_volume['state'] = 'disabled'
            self.combobox_database_volume['state'] = 'readonly'
            self.wrapper_database_volume.set('')
            self.combobox_database_volume['state'] = 'disabled'
            self.combobox_read_volume['state'] = 'readonly'
            self.wrapper_read_volume.set('')
            self.combobox_read_volume['state'] = 'disabled'
            self.combobox_reference_volume['state'] = 'readonly'
            self.wrapper_reference_volume.set('')
            self.combobox_reference_volume['state'] = 'disabled'
            self.combobox_result_volume['state'] = 'readonly'
            self.wrapper_result_volume.set('')
            self.combobox_result_volume['state'] = 'disabled'
        elif self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_multivolume():
            self.combobox_ngscloud_volume['state'] = 'readonly'
            self.wrapper_ngscloud_volume.set('')
            self.combobox_ngscloud_volume['state'] = 'disabled'
            self.combobox_app_volume['state'] = 'readonly'
            self.wrapper_app_volume.set('')
            self.combobox_database_volume['state'] = 'readonly'
            self.wrapper_database_volume.set('')
            self.combobox_read_volume['state'] = 'readonly'
            self.wrapper_read_volume.set('')
            self.combobox_reference_volume['state'] = 'readonly'
            self.wrapper_reference_volume.set('')
            self.combobox_result_volume['state'] = 'readonly'
            self.wrapper_result_volume.set('')

    #---------------

    def combobox_ngscloud_volume_selected_item(self, event=None):
        '''
        Process the event when a volume name item has been selected
        '''

        pass

    #---------------

    def combobox_app_volume_selected_item(self, event=None):
        '''
        Process the event when a application volume item has been selected
        '''

        pass

    #---------------

    def combobox_database_volume_selected_item(self, event=None):
        '''
        Process the event when a database volume item has been selected
        '''

        pass

    #---------------

    def combobox_read_volume_selected_item(self, event=None):
        '''
        Process the event when a read volume item has been selected
        '''

        pass

    #---------------

    def combobox_reference_volume_selected_item(self, event=None):
        '''
        Process the event when a reference volume item has been selected
        '''

        pass

    #---------------

    def combobox_result_volume_selected_item(self, event=None):
        '''
        Process the event when a result volume item has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormLinkVolumes" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_none() or self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_singlevolume() and self.wrapper_ngscloud_volume.get() != '' or self.wrapper_dataset_structure.get() == xconfiguration.get_dataset_structure_multivolume() and self.wrapper_app_volume.get() != '' and self.wrapper_database_volume.get() != '' and self.wrapper_read_volume.get() != '' and self.wrapper_reference_volume.get() != '' and self.wrapper_result_volume.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the link volumes to the configuration.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # prepare and check volume names
        if OK:
            used_volume_name_list = []

            ngscloud_volume = self.wrapper_ngscloud_volume.get()
            if ngscloud_volume == 'NONE': ngscloud_volume = ''
            app_volume = self.wrapper_app_volume.get()
            if app_volume == 'NONE': app_volume = ''
            database_volume = self.wrapper_database_volume.get()
            if database_volume == 'NONE': database_volume = ''
            read_volume = self.wrapper_read_volume.get()
            if read_volume == 'NONE': read_volume = ''
            reference_volume = self.wrapper_reference_volume.get()
            if reference_volume == 'NONE': reference_volume = ''
            result_volume = self.wrapper_result_volume.get()
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
                message = 'A volume can be linked only once.'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the region and zone update in the NGScloud config file
        if OK:
            message = 'The dataset structure are going to be modified.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # link volumes in the NGScloud config file corresponding to the environment
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xconfiguration.link_volumes.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xconfiguration.link_volumes, args=(self.wrapper_dataset_structure.get(), ngscloud_volume, app_volume, database_volume, read_volume, reference_volume, result_volume, dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormLinkVolumes".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

class FormCreateCluster(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormCreateCluster" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Cluster operation - Create cluster'

        # get current region name
        self.region_name = xconfiguration.get_current_region_name()

        # get the NGScloud config file
        ngscloud_config_file = xconfiguration.get_ngscloud_config_file()

        # get the option dictionary corresponding to the NGScloud config file
        ngscloud_options_dict = xlib.get_option_dict(ngscloud_config_file)

        # get the dataset structure and NGScloud_volume
        self.dataset_structure = ngscloud_options_dict['dataset info']['dataset_structure'].lower()

        # check the dataset structure
        if self.dataset_structure == xconfiguration.get_dataset_structure_none():
            message = f'The dataset structure is not {xconfiguration.get_dataset_structure_singlevolume()} nor {xconfiguration.get_dataset_structure_multivolume()}, then datasets will be created in the root volume and they will be lost when the cluster is terminated.\n\nThe dataset structure can be modified in:\n\n"Cloud control" -> "Configuration" -> "Links volumes"'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the volume type dictinary
        self.volume_type_dict = xec2.get_volume_type_dict()

        # create the wrappers to track changes in inputs
        self.wrapper_cluster_mode = tkinter.StringVar()
        self.wrapper_cluster_mode.trace('w', self.check_inputs)
        self.wrapper_instance_type = tkinter.StringVar()
        self.wrapper_instance_type.trace('w', self.check_inputs)
        self.wrapper_volume_type_text = tkinter.StringVar()
        self.wrapper_volume_type_text.trace('w', self.check_inputs)
        self.wrapper_volume_size = tkinter.StringVar()
        self.wrapper_volume_size.trace('w', self.check_inputs)
        self.wrapper_purchasing_option = tkinter.StringVar()
        self.wrapper_purchasing_option.trace('w', self.check_inputs)
        self.wrapper_max_spot_price = tkinter.StringVar()
        self.wrapper_max_spot_price.trace('w', self.check_inputs)
        self.wrapper_interruption_behavior = tkinter.StringVar()
        self.wrapper_interruption_behavior.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

        # show warnings about characteristics and pricing
        message = 'You can consult the characteristics of the EC2 intance types in:\n\n'
        message += 'https://aws.amazon.com/ec2/instance-types/\n\n'
        message += 'and the EC2 pricing is detailed in:\n\n'
        message += 'https://aws.amazon.com/ec2/pricing/'
        tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormCreateCluster".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_open_folder"
        image_select_instance_type = PIL.Image.open('./image_select_instance_type.png')
        imagetk_select_instance_type = PIL.ImageTk.PhotoImage(image_select_instance_type)  

        # create "label_cluster_mode" and register it with the grid geometry manager
        self.label_cluster_mode = tkinter.Label(self, text='Cluster mode')
        self.label_cluster_mode.grid(row=0, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_cluster_mode" and register it with the grid geometry manager
        self.combobox_cluster_mode = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_cluster_mode)
        self.combobox_cluster_mode.grid(row=0, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_instance_type" and register it with the grid geometry manager
        self.label_instance_type = tkinter.Label(self, text='Instance type')
        self.label_instance_type.grid(row=1, column=0, padx=(15,5), pady=(30,5), sticky='e')

        # create "entry_instance_type" and register it with the grid geometry manager
        self.entry_instance_type = tkinter.Entry(self, width=30, state='disabled', textvariable=self.wrapper_instance_type)
        self.entry_instance_type.grid(row=1, column=1, padx=(5,5), pady=(30,5), sticky='w')

        # create "button_select_instance_type" and register it with the grid geometry manager
        self.button_select_instance_type = tkinter.ttk.Button(self, image=imagetk_select_instance_type, command=self.select_instance_type)
        self.button_select_instance_type.image = imagetk_select_instance_type
        self.button_select_instance_type.grid(row=1, column=2, padx=(5,0), pady=(30,5), sticky='w')

        # create "label_instance_type_warning" and register it with the grid geometry manager
        self.label_instance_type_warning = tkinter.Label(self, text='')
        self.label_instance_type_warning.grid(row=2, column=1, columnspan=4, padx=(5,5), pady=(0,5), sticky='w')

        # create "label_volume_type_text" and register it with the grid geometry manager
        self.label_volume_type_text = tkinter.Label(self, text='Root volume type')
        self.label_volume_type_text.grid(row=3, column=0, padx=(15,5), pady=(5,5), sticky='e')

        # create "combobox_volume_type_text" and register it with the grid geometry manager
        self.combobox_volume_type_text = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_volume_type_text)
        self.combobox_volume_type_text.grid(row=3, column=1, padx=(5,5), pady=(5,5), sticky='w')

        # create "label_volume_size" and register it with the grid geometry manager
        self.label_volume_size = tkinter.Label(self, text='Root size (in GiB)')
        self.label_volume_size.grid(row=4, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_volume_size" and register it with the grid geometry manager
        self.entry_volume_size = tkinter.Entry(self, width=10, textvariable=self.wrapper_volume_size, validatecommand=self.check_inputs)
        self.entry_volume_size.grid(row=4, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_volume_size_warning" and register it with the grid geometry manager
        self.label_volume_size_warning = tkinter.Label(self, text='')
        self.label_volume_size_warning.grid(row=4, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_purchasing_option" and register it with the grid geometry manager
        self.label_purchasing_option = tkinter.Label(self, text='Purchasing option')
        self.label_purchasing_option.grid(row=5, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_purchasing_option" and register it with the grid geometry manager
        self.combobox_purchasing_option = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_purchasing_option)
        self.combobox_purchasing_option.grid(row=5, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_max_spot_price" and register it with the grid geometry manager
        self.label_max_spot_price = tkinter.Label(self, text='Maximum spot price (in $)')
        self.label_max_spot_price.grid(row=6, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "entry_max_spot_price" and register it with the grid geometry manager
        self.entry_max_spot_price = tkinter.Entry(self, width=10, textvariable=self.wrapper_max_spot_price, validatecommand=self.check_inputs)
        self.entry_max_spot_price.grid(row=6, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_max_spot_price_warning" and register it with the grid geometry manager
        self.label_max_spot_price_warning = tkinter.Label(self, text='')
        self.label_max_spot_price_warning.grid(row=6, column=2, columnspan=3, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_interruption_behavior" and register it with the grid geometry manager
        self.label_interruption_behavior = tkinter.Label(self, text='Interruption behavior')
        self.label_interruption_behavior.grid(row=7, column=0, padx=(15,5), pady=(35,5), sticky='e')

        # create "combobox_interruption_behavior" and register it with the grid geometry manager
        self.combobox_interruption_behavior = tkinter.ttk.Combobox(self, width=30, height=4, state='readonly', textvariable=self.wrapper_interruption_behavior)
        self.combobox_interruption_behavior.grid(row=7, column=1, padx=(5,5), pady=(35,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(50+xlib.get_os_size_fix()))
        self.label_fit.grid(row=8, column=2, padx=(0,0), pady=(35,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=8, column=3, padx=(5,5), pady=(35,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=8, column=4, padx=(5,5), pady=(35,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_mode.bind('<<ComboboxSelected>>', self.combobox_cluster_mode_selected_item)
        self.combobox_volume_type_text.bind('<<ComboboxSelected>>', self.combobox_volume_type_text_selected_item)
        self.combobox_purchasing_option.bind('<<ComboboxSelected>>', self.combobox_purchasing_option_selected_item)
        self.combobox_interruption_behavior.bind('<<ComboboxSelected>>', self.combobox_interruption_behavior_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # initialize the minimum and maximum volume size
        self.minimum_size = 0
        self.maximum_size = 0

        # populate data in comboboxes
        self.populate_combobox_cluster_mode()
        self.populate_combobox_volume_type_text()
        self.populate_combobox_purchasing_option()
        self.populate_combobox_interruption_behavior()

        # load initial data in inputs
        self.wrapper_cluster_mode.set('native')
        self.wrapper_volume_type_text.set('general purpose SSD')
        self.volume_type_id = xec2.get_volume_type_id(self.wrapper_volume_type_text.get())
        self.wrapper_volume_size.set('8')
        self.minimum_size = 8
        self.maximum_size = self.volume_type_dict[self.volume_type_id]['maximum_size']
        self.wrapper_purchasing_option.set(xec2.get_purchasing_option_ondemand())
        self.wrapper_max_spot_price.set('')
        self.entry_max_spot_price['state'] = 'disabled'
        self.wrapper_interruption_behavior.set('')
        self.combobox_interruption_behavior['state'] = 'disabled'

    #---------------

    def populate_combobox_cluster_mode(self):
        '''
        Populate data in "combobox_cluster_mode".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_mode.set('')

        # check if StarCluster is installed
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
        cluster_mode_list = xconfiguration.get_cluster_mode_list()
        if not is_starcluster_installed:
            cluster_mode_list.remove(xconfiguration.get_cluster_mode_starcluster())

        # load the cluster mode list in the combobox
        self.combobox_cluster_mode['values'] = cluster_mode_list

    #---------------

    def populate_combobox_volume_type_text(self):
        '''
        Populate data in "combobox_volume_type_text".
        '''

        # clear the value selected in the combobox
        self.wrapper_volume_type_text.set('')

        # load the volume types in the combobox
        self.combobox_volume_type_text['values'] = xec2.get_volume_type_text_list(only_possible_root_disk=True)

    #---------------

    def populate_combobox_purchasing_option(self):
        '''
        Populate data in "combobox_purchasing_option".
        '''

        # clear the value selected in the combobox
        self.wrapper_purchasing_option.set('')

        # load the volume types in the combobox
        self.combobox_purchasing_option['values'] = xec2.get_purchasing_option_list()

    #---------------

    def populate_combobox_interruption_behavior(self):
        '''
        Populate data in "combobox_interruption_behavior".
        '''

        # clear the value selected in the combobox
        self.wrapper_interruption_behavior.set('')

        # load the volume types in the combobox
        self.combobox_interruption_behavior['values'] = xec2.get_interruption_behavior_list()

    #---------------

    def select_instance_type(self):
        '''
        Select a instance type and assign it to "entry_instance_type".
        '''

        # get the instance type dictionary
        instance_type_dict = xconfiguration.get_instance_type_dict(self.wrapper_cluster_mode.get())

        # build the data list
        data_list = ['use', 'id', 'vcpu', 'memory', 'processor', 'speed', 'nitro', 'starcluster', 'generation']

        # build the data dictionary
        data_dict = {}
        data_dict['use'] = {'text': 'Use', 'width': 135, 'alignment': 'left'}
        data_dict['id'] = {'text': 'Instance Type', 'width': 120, 'alignment': 'left'}
        data_dict['vcpu'] = {'text': 'vCPUs', 'width': 50, 'alignment': 'right'}
        data_dict['memory'] = {'text': 'Memory (GiB)', 'width': 100, 'alignment': 'right'}
        data_dict['processor'] = {'text': 'Processor', 'width': 290, 'alignment': 'left'}
        data_dict['speed'] = {'text': 'Clock Speed', 'width': 100, 'alignment': 'right'}
        data_dict['nitro'] = {'text': 'Nitro System', 'width': 90, 'alignment': 'left'}
        data_dict['starcluster'] = {'text': 'StarCluster', 'width': 110, 'alignment': 'left'}
        data_dict['generation'] = {'text': 'Generation', 'width': 90, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        dialog_table = gdialogs.DialogTable(self, 'Instance type selection', 400, 1100, data_list, data_dict, instance_type_dict, sorted(instance_type_dict.keys()), action='select_instace_type')
        self.wait_window(dialog_table)

    #---------------

    def combobox_cluster_mode_selected_item(self, event=None):
        '''
        Process the event when a type instance has been selected.
        '''

        # clear the value selected in the combobox
        self.wrapper_instance_type.set('')

        # clear the value shown in the warning
        self.label_instance_type_warning['text'] = ''


        # set data of "entry_max_spot_price" and "combobox_interruption_behavior"
        if self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_native():
            self.combobox_volume_type_text['state'] = 'readonly'
            self.wrapper_volume_type_text.set('general purpose SSD')
            self.volume_type_id = xec2.get_volume_type_id(self.wrapper_volume_type_text.get())
            self.entry_volume_size['state'] = 'normal'
            self.wrapper_volume_size.set('8')
            self.minimum_size = 8
            self.maximum_size = self.volume_type_dict[self.volume_type_id]['maximum_size']
            self.combobox_purchasing_option['state'] = 'readonly'
            self.wrapper_purchasing_option.set(xec2.get_purchasing_option_ondemand())
            self.wrapper_max_spot_price.set('')
            self.entry_max_spot_price['state'] = 'disabled'
            self.wrapper_interruption_behavior.set('')
            self.combobox_interruption_behavior['state'] = 'disabled'
        elif self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_starcluster():
            self.combobox_volume_type_text['state'] = 'readonly'
            self.wrapper_volume_type_text.set('')
            self.combobox_volume_type_text['state'] = 'disable'
            self.entry_volume_size['state'] = 'normal'
            self.wrapper_volume_size.set('')
            self.entry_volume_size['state'] = 'disable'
            self.combobox_purchasing_option['state'] = 'readonly'
            self.wrapper_purchasing_option.set('')
            self.combobox_purchasing_option['state'] = 'disable'
            self.entry_max_spot_price['state'] = 'normal'
            self.wrapper_max_spot_price.set('')
            self.entry_max_spot_price['state'] = 'disabled'
            self.combobox_interruption_behavior['state'] = 'normal'
            self.wrapper_interruption_behavior.set('')
            self.combobox_interruption_behavior['state'] = 'disabled'

    #---------------

    def combobox_volume_type_text_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_volume_type_text" has been selected.
        '''

        # get the volume type identification
        self.volume_type_id = xec2.get_volume_type_id(self.wrapper_volume_type_text.get())

        # get the minimum and maximun size of the volume type
        self.minimum_size = 8
        self.maximum_size = self.volume_type_dict[self.volume_type_id]['maximum_size']

        # set the minimum value in volume size "entry_volume_size"
        self.wrapper_volume_size.set(self.minimum_size)

    #---------------

    def combobox_purchasing_option_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_purchasing_option" has been selected.
        '''

        # set data of "entry_max_spot_price" and "combobox_interruption_behavior"
        if self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_ondemand():
            self.entry_max_spot_price['state'] = 'normal'
            self.combobox_interruption_behavior['state'] = 'readonly'
            self.wrapper_max_spot_price.set('')
            self.wrapper_interruption_behavior.set('')
            self.entry_max_spot_price['state'] = 'disabled'
            self.combobox_interruption_behavior['state'] = 'disable'
        elif self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_spot():
            self.entry_max_spot_price['state'] = 'normal'
            self.combobox_interruption_behavior['state'] = 'readonly'
            self.wrapper_max_spot_price.set('')
            self.wrapper_interruption_behavior.set(xec2.get_interruption_behavior_terminate())

    #---------------

    def combobox_interruption_behavior_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_interruption_behavior" has been selected.
        '''

        pass

    #---------------

    def write_label_instance_type_warning(self, event=None):
        '''
        Print data of the instance type that has been selected.
        '''

        # get the data dictoonary corresponding to the instance type
        instance_type_data_dict = xconfiguration.get_instance_type_data_dict(self.wrapper_instance_type.get())

        # set the instance type description
        description = f'Use: {instance_type_data_dict["use"]} - vCPU: {instance_type_data_dict["vcpu"]} - Memory: {instance_type_data_dict["memory"]} GiB - Generation: {instance_type_data_dict["generation"]}'

        # show the description in the instance type warning
        self.label_instance_type_warning['text'] = description

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormCreateCluster" and do the actions linked to its value.
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_volume_size"
        if not self.check_entry_volume_size():
            OK = False

        # check the content of "entry_max_spot_price"
        if not self.check_entry_max_spot_price():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_instance_type.get() != '' and (self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_starcluster() or self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_native() and self.wrapper_volume_type_text.get() != '' and self.wrapper_volume_size.get() != '' and (self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_ondemand() or self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_spot() and self.wrapper_interruption_behavior.get() != '')):
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_volume_size(self):
        '''
        Check the content of "entry_volume_size".
        '''

        # initialize the control variable
        OK = True

        # check that "entry_volume_size" is an integer value between the minimum and maximum value size
        if self.wrapper_volume_type_text.get() != '':
            self.label_volume_size_warning['text'] = f'Integer number between {self.minimum_size} and {self.maximum_size}'
            self.label_volume_size_warning['foreground'] = 'black'
            if self.wrapper_volume_size.get() != '':
                try:
                    volume_size = int(self.wrapper_volume_size.get())
                except:
                    self.label_volume_size_warning['foreground'] = 'red'
                    OK = False
                else:
                    if volume_size < self.minimum_size or volume_size > self.maximum_size:
                        self.label_volume_size_warning['foreground'] = 'red'
                        OK = False
        else:
            self.label_volume_size_warning['text'] = ''
            self.label_volume_size_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def check_entry_max_spot_price(self):
        '''
        Check the content of "entry_max_spot_price".
        '''

        # initialize the control variable
        OK = True

        # check that "entry_max_spot_price" a float number greater than 0.0
        if self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_ondemand():
            self.label_max_spot_price_warning['text'] = ''
            self.label_max_spot_price_warning['foreground'] = 'black'
        elif self.wrapper_purchasing_option.get() == xec2.get_purchasing_option_spot():
            self.label_max_spot_price_warning['text'] = 'Float number greater than 0.0 (empty if on-demand price).'
            self.label_max_spot_price_warning['foreground'] = 'black'
            if self.wrapper_max_spot_price.get() != '':
                try:
                    max_spot_price = float(self.wrapper_max_spot_price.get())
                except:
                    self.label_max_spot_price_warning['foreground'] = 'red'
                    OK = False
                else:
                    if max_spot_price <= 0.0:
                        self.label_max_spot_price_warning['foreground'] = 'red'
                        OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the cluster creation process.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # check the AMI identification is create in the current zone
        if OK:
            if self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_native() and xec2.get_ubuntu_ami_id(self.region_name) == xec2.get_unknown_ami_id():
                message = f'*** ERROR: The AMI {xec2.get_ubuntu_ami_name()} is not found in the region {self.region_name}.'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False
            elif self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_starcluster() and xec2.get_starcluster_ami_id(self.region_name) == xec2.get_unknown_ami_id():
                message = f'*** ERROR: The AMI {xec2.get_starcluster_ami_name()} is not found in the region {self.region_name}.'
                tkinter.messagebox.showerror(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False

        # confirm the creation of the cluster
        if OK:
            if self.dataset_structure == xconfiguration.get_dataset_structure_none():
                message = f'The cluster with instace type {self.wrapper_instance_type.get()} is going to be created.\n\nThe dataset structure is not {xconfiguration.get_dataset_structure_singlevolume()} nor {xconfiguration.get_dataset_structure_multivolume()}, then datasets will be created in the root volume and they will be lost when the cluster is terminated.\n\nAre you sure to continue?\n\n'.format(self.wrapper_instance_type.get())
            else:

                message = f'The cluster with instace type {self.wrapper_instance_type.get()} is going to be created.\n\nAre you sure to continue?\n\n'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # create the cluster
        if OK:
            if self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_native():
                dialog_log = gdialogs.DialogLog(self, self.head, xinstance.create_instance.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                max_spot_price = 0.0 if self.wrapper_volume_size.get() == '' else float(self.wrapper_volume_size.get())
                interruption_behavior = xec2.get_interruption_behavior_terminate() if self.wrapper_interruption_behavior.get() == '' else self.wrapper_interruption_behavior.get()
                threading.Thread(target=xinstance.create_instance, args=(self.wrapper_instance_type.get(), dialog_log, self.volume_type_id, int(self.wrapper_volume_size.get()), self.wrapper_purchasing_option.get(), max_spot_price, interruption_behavior, lambda: dialog_log.enable_button_close(), True)).start()
            elif self.wrapper_cluster_mode.get() == xconfiguration.get_cluster_mode_starcluster():
                template_name = xconfiguration.build_cluster_name(self.wrapper_instance_type.get())
                cluster_name = xconfiguration.build_cluster_name(self.wrapper_instance_type.get())
                dialog_log = gdialogs.DialogLog(self, self.head, xcluster.create_cluster.__name__)
                threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                threading.Thread(target=xcluster.create_cluster, args=(template_name, cluster_name, dialog_log, lambda: dialog_log.enable_button_close(), True)).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormCreateCluster".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormTerminateCluster(tkinter.Frame):

    #---------------

    def __init__(self, parent, main, force):
        '''
        Execute actions correspending to the creation of a "FormTerminateCluster" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main
        self.force = force

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = f'Cluster operation - {"Terminate cluster" if not self.force else "Force termination of a cluster"}'

        # create the wrappers to track changes in the inputs
        if not self.force:
            self.wrapper_cluster_name = tkinter.StringVar()
            self.wrapper_cluster_name.trace('w', self.check_inputs)
        else:
            self.wrapper_instance_type = tkinter.StringVar()
            self.wrapper_instance_type.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormTerminateCluster".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "image_open_folder"
        image_select_instance_type = PIL.Image.open('./image_select_instance_type.png')
        imagetk_select_instance_type = PIL.ImageTk.PhotoImage(image_select_instance_type)  

        # create "label_cluster_name" and register it with the grid geometry manager
        if not self.force:
            self.label_cluster_name = tkinter.Label(self, text='Cluster name')
            self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        if not self.force:
            self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
            self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_instance_type" and register it with the grid geometry manager
        if self.force:
            self.label_instance_type = tkinter.Label(self, text='Instance type')
            self.label_instance_type.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "entry_instance_type" and register it with the grid geometry manager
        if self.force:
            self.entry_instance_type = tkinter.Entry(self, width=30, state='disabled', textvariable=self.wrapper_instance_type)
            self.entry_instance_type.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "button_select_instance_type" and register it with the grid geometry manager
        if self.force:
            self.button_select_instance_type = tkinter.ttk.Button(self, image=imagetk_select_instance_type, command=self.select_instance_type)
            self.button_select_instance_type.image = imagetk_select_instance_type
            self.button_select_instance_type.grid(row=0, column=2, padx=(5,0), pady=(75,5), sticky='w')

        # create "label_instance_type_warning" and register it with the grid geometry manager
        if self.force:
            self.label_instance_type_warning = tkinter.Label(self, text='')
            self.label_instance_type_warning.grid(row=1, column=1, columnspan=4, padx=(5,5), pady=(5,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        if not self.force:
            self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        else:
            self.label_fit = tkinter.Label(self, text=' '*(75+xlib.get_os_size_fix()))
        self.label_fit.grid(row=1, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        if not self.force:
            self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        if not self.force:
            self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def select_instance_type(self):
        '''
        Select a instance type and assign it to "entry_instance_type".
        '''

        # get the instance type dictionary
        instance_type_dict = xconfiguration.get_instance_type_dict(xconfiguration.get_cluster_mode_native())

        # build the data list
        data_list = ['use', 'id', 'vcpu', 'memory', 'processor', 'speed', 'nitro', 'starcluster', 'generation']

        # build the data dictionary
        data_dict = {}
        data_dict['use'] = {'text': 'Use', 'width': 135, 'alignment': 'left'}
        data_dict['id'] = {'text': 'Instance Type', 'width': 120, 'alignment': 'left'}
        data_dict['vcpu'] = {'text': 'vCPUs', 'width': 50, 'alignment': 'right'}
        data_dict['memory'] = {'text': 'Memory (GiB)', 'width': 100, 'alignment': 'right'}
        data_dict['processor'] = {'text': 'Processor', 'width': 290, 'alignment': 'left'}
        data_dict['speed'] = {'text': 'Clock Speed', 'width': 100, 'alignment': 'right'}
        data_dict['nitro'] = {'text': 'Nitro System', 'width': 90, 'alignment': 'left'}
        data_dict['starcluster'] = {'text': 'StarCluster', 'width': 110, 'alignment': 'left'}
        data_dict['generation'] = {'text': 'Generation', 'width': 90, 'alignment': 'left'}

        # create the dialog Table to show the nodes running
        dialog_table = gdialogs.DialogTable(self, 'Instance type selection', 400, 1100, data_list, data_dict, instance_type_dict, sorted(instance_type_dict.keys()), action='select_instace_type')
        self.wait_window(dialog_table)

    #---------------

    def write_label_instance_type_warning(self, event=None):
        '''
        Process the event when a instance type item has been selected
        '''

        # get the data dictoonary corresponding to the instance type
        instance_type_data_dict = xconfiguration.get_instance_type_data_dict(self.wrapper_instance_type.get())

        # set the instance type description
        description = f'Use: {instance_type_data_dict["use"]} - vCPU: {instance_type_data_dict["vcpu"]} - Memory: {instance_type_data_dict["memory"]} GiB - Generation: {instance_type_data_dict["generation"]}'

        # show the instance type description in the instance type warning
        self.label_instance_type_warning['text'] = description

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormTerminateCluster" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if not self.force:
            if self.wrapper_cluster_name.get() != '':
                self.button_execute['state'] = 'enable'
            else:
                self.button_execute['state'] = 'disabled'
        else:
            if self.wrapper_instance_type.get() != '':
                self.button_execute['state'] = 'enable'
            else:
                self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the cluster termination process.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the termination of the cluster
        if OK:
            if not self.force:
                message = f'The cluster {self.wrapper_cluster_name.get()} is going to be terminated.\n\nAre you sure to continue?'
            else:
                message = f'The cluster  {xconfiguration.build_cluster_name(self.wrapper_instance_type.get())} is going to be forced to terminate..\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # terminate the cluster and initialize inputs
        if OK:
            if not self.force:
                if xec2.get_cluster_mode(self.wrapper_cluster_name.get()) == xconfiguration.get_cluster_mode_native():
                    dialog_log = gdialogs.DialogLog(self, self.head, xinstance.terminate_instance.__name__)
                    threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                    threading.Thread(target=xinstance.terminate_instance, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()
                elif xec2.get_cluster_mode(self.wrapper_cluster_name.get()) == xconfiguration.get_cluster_mode_starcluster():
                    dialog_log = gdialogs.DialogLog(self, self.head, xcluster.terminate_cluster.__name__)
                    threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                    threading.Thread(target=xcluster.terminate_cluster, args=(self.wrapper_cluster_name.get(), self.force, dialog_log, lambda: dialog_log.enable_button_close())).start()
            else:
                cluster_name = xconfiguration.build_cluster_name(self.wrapper_instance_type.get())
                if xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_native() or xec2.get_cluster_mode(cluster_name) is None:
                    dialog_log = gdialogs.DialogLog(self, self.head, xinstance.terminate_instance.__name__)
                    threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                    threading.Thread(target=xinstance.terminate_instance, args=(cluster_name, dialog_log, lambda: dialog_log.enable_button_close())).start()
                elif xec2.get_cluster_mode(cluster_name) == xconfiguration.get_cluster_mode_starcluster():
                    dialog_log = gdialogs.DialogLog(self, self.head, xcluster.terminate_cluster.__name__)
                    threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
                    threading.Thread(target=xcluster.terminate_cluster, args=(cluster_name, self.force, dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormTerminateCluster".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormShowClusterComposition(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormShowClusterComposition" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Cluster operation - Show cluster composition'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormShowClusterComposition".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=1, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=1, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=1, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running in mode StarCluster, set its cluster name by default
        if len(running_cluster_list) == 1 and xec2.get_cluster_mode(running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormShowClusterComposition" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        OK = self.check_combobox_cluster_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the show of cluster information of every node: OS, CPU number and memory.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # show the status of batch jobs in the cluster
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xcluster.show_cluster_composition.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xcluster.show_cluster_composition, args=(self.wrapper_cluster_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormShowClusterComposition".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormShowStatusBatchJobs(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormShowStatusBatchJobs" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Cluster operation - Show status of batch jobs'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # initialize the cluster name previously selected
        self.cluster_name_ant = None

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormShowStatusBatchJobs".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=1, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=1, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=1, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running in mode StarCluster, set its cluster name by default
        if len(running_cluster_list) == 1 and xec2.get_cluster_mode(running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def check_inputs(self, *args):
        '''
        check the content of each input of "FormShowStatusBatchJobs" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        OK = self.check_combobox_cluster_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the show of the status of batch jobs in the cluster.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # get the batch job dictionary
        if OK:
            (OK, error_list, batch_job_dict) = xcluster.get_batch_job_dict(self.ssh_client)

        # check if there are any batch jobs
        if OK:
            if batch_job_dict == {}:
                message = 'There is not any batch job.'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False

        # build the data list
        if OK:
            data_list = ['job_id', 'job_name', 'state', 'start_date', 'start_time']

        # build the data dictionary
        if OK:
            data_dict = {}
            data_dict['job_id'] = {'text': 'Job id', 'width': 50, 'alignment': 'right'}
            data_dict['job_name'] = {'text': 'Job name', 'width': 100, 'alignment': 'left'}
            data_dict['state'] = {'text': 'State', 'width': 150, 'alignment': 'left'}
            data_dict['start_date'] = {'text': 'Start date', 'width': 100, 'alignment': 'right'}
            data_dict['start_time'] = {'text': 'Start time', 'width': 100, 'alignment': 'right'}

        # create the dialog Table to show the nodes running
        if OK:
            dialog_table = gdialogs.DialogTable(self, self.head, 400, 900, data_list, data_dict, batch_job_dict, sorted(batch_job_dict.keys()))
            self.wait_window(dialog_table)

        # close the form
        self.close()

    #---------------

    def close(self):
        '''
        Close "FormShowStatusBatchJobs".
        '''

        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormKillBatchJob(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormKillBatchJob" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Cluster operation - Kill batch job'

        # initialize the SSH client connection and previous cluster name
        self.ssh_client = None
        self.cluster_name_ant = None

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_job = tkinter.StringVar()
        self.wrapper_job.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormKillBatchJob".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_job" and register it with the grid geometry manager
        self.label_job = tkinter.Label(self, text='Job')
        self.label_job.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_job" and register it with the grid geometry manager
        self.combobox_job = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_job)
        self.combobox_job.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=2, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_job.bind('<<ComboboxSelected>>', self.combobox_job_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_job['values'] = []
        self.wrapper_job.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox_cluster_name
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running in mode StarCluster, set its cluster name by default
        if len(running_cluster_list) == 1 and xec2.get_cluster_mode(running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_job(self):
        '''
        Populate data in "combobox_job".
        '''

        # clear the value selected in the combobox
        self.wrapper_job.set('')

        # get the batch job dictionary
        (OK, error_list, batch_job_dict) = xcluster.get_batch_job_dict(self.ssh_client)
        if batch_job_dict == {}:
            message = 'There is not any batch job.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # build the list of batch jobs
        batch_job_id_list = []
        for job_id in batch_job_dict.keys():
            job_name = batch_job_dict[job_id]['job_name']
            state_name = batch_job_dict[job_id]['state_name']
            start_date = batch_job_dict[job_id]['start_date']
            start_time = batch_job_dict[job_id]['start_time']
            batch_job_id_list.append('{0:>3} ({1}; started at: {2} {3}; state: {4})'.format(job_id, job_name, start_date, start_time, state_name))
        batch_job_id_list.sort()

        # load the names of cluster nodes created in the combobox
        self.combobox_job['values'] = batch_job_id_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # check if the cluster name selected is different to the previous cluster name
        if self.wrapper_cluster_name.get() != self.cluster_name_ant:

            # close SSH client connection
            if self.cluster_name_ant is not None:
                xssh.close_ssh_client_connection(self.ssh_client)

            # create the SSH client connection
            (OK, error_list, self.ssh_client) = xssh.create_ssh_client_connection(self.wrapper_cluster_name.get())
            if not OK:
                message = ''
                for error in error_list:
                    message = f'{message}{error}\n'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                self.close()

            # save current cluster name as previous cluster name
            self.cluster_name_ant = self.wrapper_cluster_name.get()

        # load data in "combobox_job"
        self.populate_combobox_job()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_job_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_job" has been selected
        '''

        try:
            self.job_id = int(self.wrapper_job.get()[0:3])
        except:
            self.job_id = 0

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormKillBatchJob" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        OK = self.check_combobox_cluster_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_job.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the removal of a node in a cluster.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the removal of the volume in the node
        if OK:
            message = f'The batch job {self.wrapper_job.get()} is going to be killed in the cluster {self.wrapper_cluster_name.get()}.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # remove the volume and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xcluster.kill_batch_job.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xcluster.kill_batch_job, args=(self.wrapper_cluster_name.get(), self.job_id, dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormKillBatchJob".
        '''

        # close SSH client connection
        if self.cluster_name_ant is not None:
            xssh.close_ssh_client_connection(self.ssh_client)

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormAddNode(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormAddNode" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Node operation - Add node in a cluster'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_node_name = tkinter.StringVar()
        self.wrapper_node_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # initialize the cluster node list
        self.cluster_node_list = []

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormAddNode".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_node_name" and register it with the grid geometry manager
        self.label_node_name = tkinter.Label(self, text='Node name')
        self.label_node_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_node_name" and register it with the grid geometry manager
        self.entry_node_name = tkinter.Entry(self, textvariable=self.wrapper_node_name, width=40, validatecommand=self.check_inputs)
        self.entry_node_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_node_name_warning" and register it with the grid geometry manager
        self.label_node_name_warning = tkinter.Label(self, text='')
        self.label_node_name_warning.grid(row=1, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=2, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=2, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=2, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.wrapper_node_name.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox_node_name
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running in mode StarCluster, set its cluster name by default
        if len(running_cluster_list) == 1 and xec2.get_cluster_mode(running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''

        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # get the cluster node list
        self.cluster_node_list = xec2.get_cluster_node_list(self.wrapper_cluster_name.get())

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormAddNode" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        OK = self.check_combobox_cluster_name()

        # check the content of "entry_node_name"
        OK = self.check_entry_node_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_node_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # check if the maximum number of instances is already running
        if len(self.cluster_node_list) >= xec2.get_max_node_number():
            message = f'The maximum number ({xec2.get_max_node_number()}) of instances is already running.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def check_entry_node_name(self):
        '''
        Check the content of "entry_node_name"
        '''

        # initialize the control variable
        OK = True

        # check if the node name is valid
        if self.wrapper_node_name.get() == 'master':
            self.label_node_name_warning['text'] = 'The name master cannot be used.'
            self.label_node_name_warning['foreground'] = 'red'
            OK = False
        elif self.wrapper_node_name.get() != '' and not self.wrapper_node_name.get().isalnum():
            self.label_node_name_warning['text'] = 'It is not an alphanumeric string'
            self.label_node_name_warning['foreground'] = 'red'
            OK = False
        elif self.wrapper_node_name.get() != '' and not self.wrapper_node_name.get()[0].isalpha():
            self.label_node_name_warning['text'] = 'The first character is not alphabetic'
            self.label_node_name_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_node_name_warning['text'] = ''
            self.label_node_name_warning['foreground'] = 'black'

        # check if the a node with the node name is instances is already running
        if OK:
            if self.wrapper_node_name.get() in self.cluster_node_list:
                message = f'The {self.wrapper_node_name.get()} is already running.'
                tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
                OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the detetion moval of a node in a cluster.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the addition of the volume in the node
        if OK:
            message = f'The node {self.wrapper_node_name.get()} is going to be added in the cluster {self.wrapper_cluster_name.get()}.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # add the node
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xnode.add_node.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xnode.add_node, args=(self.wrapper_cluster_name.get(), self.wrapper_node_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormAddNode".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormRemoveNode(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormRemoveNode" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Node operation - Remove node in a cluster'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_node_name = tkinter.StringVar()
        self.wrapper_node_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRemoveNode".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_node_name" and register it with the grid geometry manager
        self.label_node_name = tkinter.Label(self, text='Node name')
        self.label_node_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_node_name" and register it with the grid geometry manager
        self.combobox_node_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_node_name)
        self.combobox_node_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_node_name.bind('<<ComboboxSelected>>', self.combobox_node_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_node_name['values'] = []
        self.wrapper_node_name.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox_cluster_name
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        self.running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if self.running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = self.running_cluster_list

        # if there is only one cluster running in mode StarCluster, set its cluster name by default
        if len(self.running_cluster_list) == 1 and xec2.get_cluster_mode(self.running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.wrapper_cluster_name.set(self.running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_node_name(self):
        '''
        Populate data in "combobox_node_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_node_name.set('')

        # get cluster node list
        cluster_node_list = xec2.get_cluster_node_list(self.wrapper_cluster_name.get())

        # remove master in cluster node list
        cluster_node_list.remove('master')

        # check if there are some nodes besides master
        if cluster_node_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running node besides the master.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of cluster nodes created in the combobox
        self.combobox_node_name['values'] = cluster_node_list

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # load data in "combobox_node_name"
        if xec2.get_cluster_mode(self.running_cluster_list[0]) == xconfiguration.get_cluster_mode_starcluster():
            self.populate_combobox_node_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_node_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRemoveNode" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "combobox_cluster_name"
        OK = self.check_combobox_cluster_name()

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_node_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_combobox_cluster_name(self):
        '''
        Check the content of "combobox_cluster_name"
        '''

        # initialize the control variable
        OK = True

        # check the cluster mode
        if self.wrapper_cluster_name.get() != '' and xec2.get_cluster_mode(self.wrapper_cluster_name.get()) != xconfiguration.get_cluster_mode_starcluster():
            message = f'This option is only available for clusters started in mode {xconfiguration.get_cluster_mode_starcluster()}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the removal of a node in a cluster.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the removal of the volume in the node
        if OK:
            message = f'The node {self.wrapper_node_name.get()} is going to be removed in the cluster {self.wrapper_cluster_name.get()}.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # remove the volume and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xnode.remove_node.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xnode.remove_node, args=(self.wrapper_cluster_name.get(), self.wrapper_node_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormRemoveNode".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormCreateVolume(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormCreateVolume" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Volume operation - Create volume'

        # get the volume type dictinary
        self.volume_type_dict = xec2.get_volume_type_dict()

        # create the wrappers to track changes in inputs
        self.wrapper_volume_name = tkinter.StringVar()
        self.wrapper_volume_name.trace('w', self.check_inputs)
        self.wrapper_volume_type_text = tkinter.StringVar()
        self.wrapper_volume_type_text.trace('w', self.check_inputs)
        self.wrapper_volume_size = tkinter.StringVar()
        self.wrapper_volume_size.trace('w', self.check_inputs)
        self.wrapper_terminate_indicator = tkinter.IntVar()
        self.wrapper_terminate_indicator.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

        # show warnings about characteristics and pricing
        message = 'You can consult the characteristics of the EBS volumes in:\n\n'
        message += 'https://aws.amazon.com/ebs/details/\n\n'
        message += 'and the EBS pricing is detailed in:\n\n'
        message += '    https://aws.amazon.com/ebs/pricing/'
        tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormCreateVolume".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_volume_name" and register it with the grid geometry manager
        self.label_volume_name = tkinter.Label(self, text='Volume name')
        self.label_volume_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "entry_volume_name" and register it with the grid geometry manager
        self.entry_volume_name = tkinter.Entry(self, textvariable=self.wrapper_volume_name, width=40, validatecommand=self.check_inputs)
        self.entry_volume_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_volume_type_text" and register it with the grid geometry manager
        self.label_volume_type_text = tkinter.Label(self, text='Volume type')
        self.label_volume_type_text.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_volume_type_text" and register it with the grid geometry manager
        self.combobox_volume_type_text = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_volume_type_text)
        self.combobox_volume_type_text.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_volume_size" and register it with the grid geometry manager
        self.label_volume_size = tkinter.Label(self, text='Volume size (in GiB)')
        self.label_volume_size.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_volume_size" and register it with the grid geometry manager
        self.entry_volume_size = tkinter.Entry(self, textvariable=self.wrapper_volume_size, width=10, validatecommand=self.check_inputs)
        self.entry_volume_size.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_volume_size_warning" and register it with the grid geometry manager
        self.label_volume_size_warning = tkinter.Label(self, text='')
        self.label_volume_size_warning.grid(row=2, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "checkbutton_terminate_indicator" and register it with the grid geometry manager
        self.checkbutton_terminate_indicator = tkinter.ttk.Checkbutton(self, text='Terminate volume creator?', variable=self.wrapper_terminate_indicator)
        self.checkbutton_terminate_indicator.grid(row=3, column=1, columnspan=2, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(43+xlib.get_os_size_fix()))
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(0,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_volume_type_text.bind('<<ComboboxSelected>>', self.combobox_volume_type_text_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # initialize the minimum and maximum volume size
        self.minimum_size = 0
        self.maximum_size = 0

        # populate data in comboboxes
        self.populate_combobox_volume_type_text()

        # load initial data in inputs
        self.wrapper_volume_name.set('')
        self.wrapper_volume_type_text.set('standard HDD')
        self.volume_type_id = xec2.get_volume_type_id(self.wrapper_volume_type_text.get())
        self.minimum_size = self.volume_type_dict[self.volume_type_id]['minimum_size']
        self.maximum_size = self.volume_type_dict[self.volume_type_id]['maximum_size']
        self.wrapper_volume_size.set(self.minimum_size)
        self.wrapper_terminate_indicator.set(1)

    #---------------

    def populate_combobox_volume_type_text(self):
        '''
        Populate data in "combobox_volume_type_text".
        '''

        # clear the value selected in the combobox
        self.wrapper_volume_type_text.set('')

        # load the volume types in the combobox
        self.combobox_volume_type_text['values'] = xec2.get_volume_type_text_list(only_possible_root_disk=False)

    #---------------

    def combobox_volume_type_text_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_volume_type_text" has been selected
        '''

        # get the volume type identification
        self.volume_type_id = xec2.get_volume_type_id(self.wrapper_volume_type_text.get())

        # get the minimum and maximun size of the volume type
        self.minimum_size = self.volume_type_dict[self.volume_type_id]['minimum_size']
        self.maximum_size = self.volume_type_dict[self.volume_type_id]['maximum_size']

        # set the minimum value in volume size "entry_volume_size"
        self.wrapper_volume_size.set(self.minimum_size)

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormCreateVolume" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_volume_size"
        if not self.check_entry_volume_size():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_volume_name.get() != '' and self.wrapper_volume_type_text.get() != '' and self.wrapper_volume_size.get() != '' and self.wrapper_terminate_indicator.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_volume_size(self):
        '''
        Check the content of "entry_volume_size"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_volume_size" an integer value and it is  greater than or equal to the minimum value size
        self.label_volume_size_warning['text'] = f'It has to be an integer number between {self.minimum_size} and {self.maximum_size}'
        self.label_volume_size_warning['foreground'] = 'black'
        try:
            volume_size = int(self.wrapper_volume_size.get())
        except:
            self.label_volume_size_warning['foreground'] = 'red'
            OK = False
        else:
            if volume_size < self.minimum_size or volume_size > self.maximum_size:
                self.label_volume_size_warning['foreground'] = 'red'
                OK = False

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the creation of the volume in the currrent zone.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the creation of the cluster
        if OK:
            message = f'The volume {self.wrapper_volume_name.get()} is going to be created.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # create the cluster and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xvolume.create_volume.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xvolume.create_volume, args=(self.wrapper_volume_name.get(), self.volume_type_id, int(self.wrapper_volume_size.get()), self.wrapper_terminate_indicator.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormCreateVolume".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

   #---------------

#-------------------------------------------------------------------------------

class FormRemoveVolume(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormRemoveVolume" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Volume operation - Remove volume'

        # create the wrappers to track changes in the inputs
        self.wrapper_volume_name = tkinter.StringVar()
        self.wrapper_volume_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormRemoveVolume".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_volume_name" and register it with the grid geometry manager
        self.label_volume_name = tkinter.Label(self, text='Volume name')
        self.label_volume_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_volume_name" and register it with the grid geometry manager
        self.combobox_volume_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_volume_name)
        self.combobox_volume_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=1, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=1, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=1, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_volume_name.bind('<<ComboboxSelected>>', self.combobox_volume_name_selected_item)

    #---------------

    def combobox_volume_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # populate data in comboboxes
        self.populate_combobox_volume_name()

    #---------------

    def populate_combobox_volume_name(self):
        '''
        Populate data in "combobox_volume_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_volume_name.set('')

        # get the current zone name
        zone_name = xconfiguration.get_current_zone_name()

        # get the volume name list
        volume_names_list = xec2.get_created_volume_name_list(zone_name)

        # check if there are any volumes created
        if volume_names_list == []:
            self.combobox_volume_name['values'] = []
            message = 'There is not any volume created.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_volume_name['values'] = volume_names_list

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormRemoveVolume" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_volume_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the read file transfer process.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the removal of the volume
        if OK:
            message = f'The volume {self.wrapper_volume_name.get()} is going to be removed.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # remove the volume and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xvolume.remove_volume.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xvolume.remove_volume, args=(self.wrapper_volume_name.get(), dialog_log, lambda: dialog_log.enable_button_close())).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormRemoveVolume".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormMountVolume(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormMountVolume" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Volume operation - Mount volume in a node'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_node_name = tkinter.StringVar()
        self.wrapper_node_name.trace('w', self.check_inputs)
        self.wrapper_volume_name = tkinter.StringVar()
        self.wrapper_volume_name.trace('w', self.check_inputs)
        self.wrapper_aws_device_file = tkinter.StringVar()
        self.wrapper_aws_device_file.trace('w', self.check_inputs)
        self.wrapper_mount_path = tkinter.StringVar()
        self.wrapper_mount_path.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormMountVolume".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_node_name" and register it with the grid geometry manager
        self.label_node_name = tkinter.Label(self, text='Node name')
        self.label_node_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_node_name" and register it with the grid geometry manager
        self.combobox_node_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_node_name)
        self.combobox_node_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_volume_name" and register it with the grid geometry manager
        self.label_volume_name = tkinter.Label(self, text='Volume name')
        self.label_volume_name.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_volume_name" and register it with the grid geometry manager
        self.combobox_volume_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_volume_name)
        self.combobox_volume_name.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_aws_device_file" and register it with the grid geometry manager
        self.label_aws_device_file = tkinter.Label(self, text='Device file')
        self.label_aws_device_file.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_aws_device_file" and register it with the grid geometry manager
        self.entry_aws_device_file = tkinter.Entry(self, textvariable=self.wrapper_aws_device_file, width=40, validatecommand=self.check_inputs)
        self.entry_aws_device_file.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_aws_device_file_warning" and register it with the grid geometry manager
        self.label_aws_device_file_warning = tkinter.Label(self, text='')
        self.label_aws_device_file_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_mount_path" and register it with the grid geometry manager
        self.label_mount_path = tkinter.Label(self, text='Mount path')
        self.label_mount_path.grid(row=4, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_mount_path" and register it with the grid geometry manager
        self.entry_mount_path = tkinter.Entry(self, textvariable=self.wrapper_mount_path, width=40, validatecommand=self.check_inputs)
        self.entry_mount_path.grid(row=4, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_mount_path_warning" and register it with the grid geometry manager
        self.label_mount_path_warning = tkinter.Label(self, text='')
        self.label_mount_path_warning.grid(row=4, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(50+xlib.get_os_size_fix()))
        self.label_fit.grid(row=5, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=5, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=5, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_node_name.bind('<<ComboboxSelected>>', self.combobox_node_name_selected_item)
        self.combobox_volume_name.bind('<<ComboboxSelected>>', self.combobox_volume_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_node_name['values'] = []
        self.wrapper_node_name.set('')
        self.wrapper_aws_device_file.set('/dev/sdm')
        self.wrapper_mount_path.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_volume_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_node_name(self):
        '''
        Populate data in "combobox_node_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_node_name.set('')

        # load the names of cluster nodes created in the combobox
        self.combobox_node_name['values'] = xec2.get_cluster_node_list(self.wrapper_cluster_name.get())

    #---------------

    def populate_combobox_volume_name(self):
        '''
        Populate data in "combobox_volume_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_volume_name.set('')

        # get current zone name
        zone_name = xconfiguration.get_current_zone_name()

        # check if there are any linked volumes
        if xec2.get_created_volume_name_list(zone_name) == []:
            message = f'There is not any volume created in the zone {zone_name}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the volume names list in the combobox
        self.combobox_volume_name['values'] = xec2.get_created_volume_name_list(zone_name)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # load data in "combobox_node_name"
        self.populate_combobox_node_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_node_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def combobox_volume_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormMountVolume" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "entry_aws_device_file"
        if not self.check_entry_aws_device_file():
            OK = False

        # check the content of "self.entry_mount_path"
        if not self.check_entry_mount_path():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_node_name.get() != '' and self.wrapper_volume_name.get() != '' and self.wrapper_aws_device_file.get() != '' and self.entry_mount_path.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_aws_device_file(self):
        '''
        Check the content of "entry_aws_device_file"
        '''

        # initialize the control variable
        OK = True

        # initialize the device file pattern
        device_file_pattern = '/dev/sd[m-p]'

        # check that "entry_aws_device_file" value is valid
        if not xlib.is_device_file(self.wrapper_aws_device_file.get(), device_file_pattern):
            self.label_aws_device_file_warning['text'] = f'It has to have a pattern {device_file_pattern}.'
            self.label_aws_device_file_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_aws_device_file_warning['text'] = ''
            self.label_aws_device_file_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def check_entry_mount_path(self):
        '''
        Check the content of "entry_mount_path"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_mount_path" value is an absolute path
        if not xlib.is_absolute_path(self.entry_mount_path.get(), 'linux'):
            self.label_mount_path_warning['text'] = 'It has to be a Linux absolute path.'
            self.label_mount_path_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_mount_path_warning['text'] = ''
            self.label_mount_path_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the mounting of a volume in a node.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the mounting of the volume in the node
        if OK:
            message = f'The volume {self.wrapper_volume_name.get()} is going to be mounted in the node {self.wrapper_node_name.get()}.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # terminate the cluster and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xvolume.mount_volume.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xvolume.mount_volume, args=(self.wrapper_cluster_name.get(), self.wrapper_node_name.get(), self.wrapper_volume_name.get(), self.wrapper_aws_device_file.get(), self.entry_mount_path.get(), dialog_log, lambda: dialog_log.enable_button_close(), True)).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormMountVolume".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormUnmountVolume(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormUnmountVolume" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Volume operation - Unmount volume in a node'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_node_name = tkinter.StringVar()
        self.wrapper_node_name.trace('w', self.check_inputs)
        self.wrapper_volume_name = tkinter.StringVar()
        self.wrapper_volume_name.trace('w', self.check_inputs)
        self.wrapper_mount_path = tkinter.StringVar()
        self.wrapper_mount_path.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormUnmountVolume".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_node_name" and register it with the grid geometry manager
        self.label_node_name = tkinter.Label(self, text='Node name')
        self.label_node_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_node_name" and register it with the grid geometry manager
        self.combobox_node_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_node_name)
        self.combobox_node_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_volume_name" and register it with the grid geometry manager
        self.label_volume_name = tkinter.Label(self, text='Volume name')
        self.label_volume_name.grid(row=2, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_volume_name" and register it with the grid geometry manager
        self.combobox_volume_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_volume_name)
        self.combobox_volume_name.grid(row=2, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_mount_path" and register it with the grid geometry manager
        self.label_mount_path = tkinter.Label(self, text='Mount path')
        self.label_mount_path.grid(row=3, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "entry_mount_path" and register it with the grid geometry manager
        self.entry_mount_path = tkinter.Entry(self, textvariable=self.wrapper_mount_path, width=40, validatecommand=self.check_inputs)
        self.entry_mount_path.grid(row=3, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_mount_path_warning" and register it with the grid geometry manager
        self.label_mount_path_warning = tkinter.Label(self, text='')
        self.label_mount_path_warning.grid(row=3, column=2, columnspan=3, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=4, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=4, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=4, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_node_name.bind('<<ComboboxSelected>>', self.combobox_node_name_selected_item)
        self.combobox_volume_name.bind('<<ComboboxSelected>>', self.combobox_volume_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_node_name['values'] = []
        self.wrapper_node_name.set('')
        self.wrapper_mount_path.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()
        self.populate_combobox_volume_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_node_name(self):
        '''
        Populate data in "combobox_node_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_node_name.set('')

        # load the names of cluster nodes created in the combobox
        self.combobox_node_name['values'] = xec2.get_cluster_node_list(self.wrapper_cluster_name.get())

    #---------------

    def populate_combobox_volume_name(self):
        '''
        Populate data in "combobox_volume_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_volume_name.set('')

        # get current zone name
        zone_name = xconfiguration.get_current_zone_name()

        # check if there are any volumes linked
        if xec2.get_created_volume_name_list(zone_name) == []:
            message = f'There is not any volume created in the zone {zone_name}.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the volume names list in the combobox
        self.combobox_volume_name['values'] = xec2.get_created_volume_name_list(zone_name)

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # load data in "combobox_node_name"
        self.populate_combobox_node_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_node_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def combobox_volume_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormUnmountVolume" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check the content of "self.entry_mount_path"
        if not self.check_entry_mount_path():
            OK = False

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_node_name.get() != '' and self.wrapper_volume_name.get() != '' and self.wrapper_mount_path.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def check_entry_mount_path(self):
        '''
        Check the content of "entry_mount_path"
        '''

        # initialize the control variable
        OK = True

        # check that "entry_mount_path" value is an absolute path
        if not xlib.is_absolute_path(self.entry_mount_path.get(), 'linux'):
            self.label_mount_path_warning['text'] = 'It has to be a Linux absolute path.'
            self.label_mount_path_warning['foreground'] = 'red'
            OK = False
        else:
            self.label_mount_path_warning['text'] = ''
            self.label_mount_path_warning['foreground'] = 'black'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the unmounting of a volume in a node.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the unmounting of the volume in the node
        if OK:
            message = f'The volume {self.wrapper_volume_name.get()} is going to be unmounted in the node {self.wrapper_node_name.get()}.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # terminate the cluster and initialize inputs
        if OK:
            dialog_log = gdialogs.DialogLog(self, self.head, xvolume.unmount_volume.__name__)
            threading.Thread(target=self.wait_window, args=(dialog_log,)).start()
            threading.Thread(target=xvolume.unmount_volume, args=(self.wrapper_cluster_name.get(), self.wrapper_node_name.get(), self.wrapper_volume_name.get(), self.entry_mount_path.get(), dialog_log, lambda: dialog_log.enable_button_close(), True)).start()

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormUnmountVolume".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

class FormOpenTerminal(tkinter.Frame):

    #---------------

    def __init__(self, parent, main):
        '''
        Execute actions correspending to the creation of a "FormOpenTerminal" instance.
        '''

        # save initial parameters in instance variables
        self.parent = parent
        self.main = main

        # call the init method of the parent class
        tkinter.Frame.__init__(self, self.parent)

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # assign the text of the "head"
        self.head = 'Cluster operation - Open a terminal'

        # create the wrappers to track changes in the inputs
        self.wrapper_cluster_name = tkinter.StringVar()
        self.wrapper_cluster_name.trace('w', self.check_inputs)
        self.wrapper_node_name = tkinter.StringVar()
        self.wrapper_node_name.trace('w', self.check_inputs)

        # build the graphical user interface
        self.build_gui()

        # load initial data in inputs
        self.initialize_inputs()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def build_gui(self):
        '''
        Build the graphical user interface of "FormOpenTerminal".
        '''

        # assign the text to the label of the current process name
        self.main.label_process['text'] = self.head

        # create "label_cluster_name" and register it with the grid geometry manager
        self.label_cluster_name = tkinter.Label(self, text='Cluster name')
        self.label_cluster_name.grid(row=0, column=0, padx=(15,5), pady=(75,5), sticky='e')

        # create "combobox_cluster_name" and register it with the grid geometry manager
        self.combobox_cluster_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_cluster_name)
        self.combobox_cluster_name.grid(row=0, column=1, padx=(5,5), pady=(75,5), sticky='w')

        # create "label_node_name" and register it with the grid geometry manager
        self.label_node_name = tkinter.Label(self, text='Node name')
        self.label_node_name.grid(row=1, column=0, padx=(15,5), pady=(45,5), sticky='e')

        # create "combobox_node_name" and register it with the grid geometry manager
        self.combobox_node_name = tkinter.ttk.Combobox(self, width=40, height=4, state='readonly', textvariable=self.wrapper_node_name)
        self.combobox_node_name.grid(row=1, column=1, padx=(5,5), pady=(45,5), sticky='w')

        # create "label_fit" and register it with the grid geometry manager
        self.label_fit = tkinter.Label(self, text=' '*(53+xlib.get_os_size_fix()))
        self.label_fit.grid(row=3, column=2, padx=(0,0), pady=(45,5), sticky='e')

        # create "button_execute" and register it with the grid geometry manager
        self.button_execute = tkinter.ttk.Button(self, text='Execute', command=self.execute, state='disabled')
        self.button_execute.grid(row=3, column=3, padx=(5,5), pady=(45,5), sticky='e')

        # create "button_close" and register it with the grid geometry manager
        self.button_close = tkinter.ttk.Button(self, text='Close', command=self.close)
        self.button_close.grid(row=3, column=4, padx=(5,5), pady=(45,5), sticky='w')

        # link a handler to events
        self.combobox_cluster_name.bind('<<ComboboxSelected>>', self.combobox_cluster_name_selected_item)
        self.combobox_node_name.bind('<<ComboboxSelected>>', self.combobox_node_name_selected_item)

    #---------------

    def initialize_inputs(self):
        '''
        Load initial data in inputs.
        '''

        # load initial data in inputs
        self.combobox_node_name['values'] = []
        self.wrapper_node_name.set('')

        # populate data in comboboxes
        self.populate_combobox_cluster_name()

    #---------------

    def populate_combobox_cluster_name(self):
        '''
        Populate data in "combobox_cluster_name".
        '''

        # clear the value selected in the combobox_node_name
        self.wrapper_cluster_name.set('')

        # check if there are some running clusters
        running_cluster_list = xec2.get_running_cluster_list(only_environment_cluster=True, volume_creator_included=False)
        if running_cluster_list == []:
            self.combobox_cluster_name['values'] = []
            message = 'There is not any running cluster.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)
            return

        # load the names of clusters which are running in the combobox
        self.combobox_cluster_name['values'] = running_cluster_list

        # if there is only one cluster running, set its cluster name by default
        if len(running_cluster_list) == 1:
            self.wrapper_cluster_name.set(running_cluster_list[0])
            self.combobox_cluster_name['state'] = 'disabled'
            self.combobox_cluster_name_selected_item()

    #---------------

    def populate_combobox_node_name(self):
        '''
        Populate data in "combobox_node_name".
        '''

        # clear the value selected in the combobox
        self.wrapper_node_name.set('')

        # load the names of cluster nodes created in the combobox
        self.combobox_node_name['values'] = xec2.get_cluster_node_list(self.wrapper_cluster_name.get())

    #---------------

    def combobox_cluster_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        # set cursor to show busy status
        self.main.config(cursor='watch')
        self.main.update()

        # load data in "combobox_node_name"
        self.populate_combobox_node_name()

        # set cursor to show normal status
        self.main.config(cursor='')
        self.main.update()

    #---------------

    def combobox_node_name_selected_item(self, event=None):
        '''
        Process the event when an item of "combobox_cluster_name" has been selected
        '''

        pass

    #---------------

    def check_inputs(self, *args):
        '''
        Check the content of each input of "FormOpenTerminal" and do the actions linked to its value
        '''

        # initialize the control variable
        OK = True

        # check if "button_execute" has to be enabled or disabled
        if self.wrapper_cluster_name.get() != '' and self.wrapper_node_name.get() != '':
            self.button_execute['state'] = 'enable'
        else:
            self.button_execute['state'] = 'disabled'

        # return the control variable
        return OK

    #---------------

    def execute(self):
        '''
        Execute the opening of window of a node cluster.
        '''

        # check inputs
        OK = self.check_inputs()
        if not OK:
            message = 'Some input values are not OK.'
            tkinter.messagebox.showwarning(f'{xlib.get_project_name()} - {self.head}', message)

        # confirm the opening of a windows terminal of a node cluster
        if OK:
            message = f'A terminal window of the node {self.wrapper_node_name.get()} in the cluster {self.wrapper_cluster_name.get()} is going to be opened.\n\nAre you sure to continue?'
            OK = tkinter.messagebox.askyesno(f'{xlib.get_project_name()} - {self.head}', message)

        # open the terminal windows
        if OK:
            xcluster.open_terminal(self.wrapper_cluster_name.get(), self.wrapper_node_name.get())

        # close the form
        if OK:
            self.close()

    #---------------

    def close(self):
        '''
        Close "FormOpenTerminal".
        '''

        # clear the label of the current process name
        self.main.label_process['text'] = ''

        # close the current form
        self.main.close_current_form()

    #---------------

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    print('This file contains the classes related to forms corresponding to Cloud Control menu items in gui mode.')
    sys.exit(0)

#-------------------------------------------------------------------------------
