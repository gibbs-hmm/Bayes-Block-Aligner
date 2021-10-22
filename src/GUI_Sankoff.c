/*


Bayes Aligner: A bayesian sequence alignment program
(Zhu, J.;Liu, J.S.; Lawrence, C.E. Bayesian adaptive sequence alignment algorithms. Bioinformatics, 14:25-39, 1998)
Please acknowledge the program authors on any publication of scientific 
results based in part on use of the program and cite the above article in 
which the program was described.

Copyright (C) 2006   Health Research Inc. 
                                                                    
HEALTH RESEARCH INCORPORATED (HRI),                   
ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                  
                                                                   
Email:  Bayes-Block@wadsworth.org
                                       
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



*/

/**************************************************************************/
/*                                                                        */
/* $Id: GUI_Sankoff.c,v 1.14 1997/07/28 15:49:12 junzhu Exp junzhu $      */
/*                                                                        */
/* Author: Jun Zhu 1996/7/15                                              */
/* Description: Create GUI for Sankoff algorithm (Bayes aligner)          */
/**************************************************************************/
#include <Xm/Xm.h>
#include <Xm/MainW.h>
#include <Xm/Text.h>
#include <Xm/RowColumn.h>
#include <Xm/ToggleBG.h>
#include <Xm/ToggleB.h>       /* radio button */
#include <Xm/FileSB.h>
#include <X11/cursorfont.h>
#include "GUI_common.h"
#include "common.h"                             
#include "routines.h" 

static void DisplayBusyCursor ( Widget w );
static void RemoveBusyCursor ( Widget w );

Widget helpshell;

/* file scope variables */
    Widget        iw_Sankoff = (Widget)0;
    Widget        iws_rowcol;         /* rowcol to contain buttons        */
    Widget        iws_label;          /* label for function               */
    Widget        buts_stop;          /* button for stop                  */
    Widget        buts_accept;        /* button for accept                */
    Widget        buts_help;          /* button for help                  */
    Widget        buts_plot_marginal; /* button for plot_marginal         */ 
    Widget        buts_plot_profile;  /* button for plot_profile          */
    Widget        buts_plot_cell;     /* button for plot_cell             */
    Widget        buts_browse;        /* button for browsing              */
    Widget        iws_queryfile;      /* for query file                   */
    Widget        iws_datafile;       /* for data file                    */
    Widget        iws_max_blocks;     /* upper limit on number of blocks  */
    Widget        iws_blocks;         /* for number of blocks             */
    Widget        iws_outfile;        /* for output file                  */
    Widget        iws_outfile_name;   /* for output file                  */
    Widget        iws_output_sequence;/* for output sequence              */
    Widget        iws_best_alignment; /* flag may turn off best alingment */
    Widget        iws_backsampling_number;
    Widget        iws_backsampling_sequence;
    Widget        iws_backsampling_cutoff;
    Widget        iws_backsampling_cutoff_value;
    Widget        iws_dna_matrix;     /* PAM index range for DNA          */
    Widget        iws_dna_stepsize;   /* PAM index stepsize for DNA       */
    Widget        iws_debug;          /* flag for debug                   */
    Widget        iws_alpha;          /* Protein or DNA                   */
    Widget        iws_special_site;   /* flag for exist special sites     */
    Widget        iws_sites;          /* enumerating special sites        */
    Widget        iws_site_weight;    /* extra weight of special sites    */
                                      /* weight for normal site is 1      */
    Widget        iws_advanced;       /* advanced options button          */
    Widget        iws_radio_protein;  /* protein radio button             */
    Widget        iws_radio_dna;      /* dna radio button                 */
    Widget        iws_exactorback;    /* exact or back                    */
    Widget        iws_radio_exact;    /* exact posterior probability RB   */
    Widget        iws_radio_back;     /* backsampling radio button        */
    Widget        iws_speedormemory;  /* speed or memory                  */
    Widget        iws_radio_speed;    /* speed radiobutton                */
    Widget        iws_radio_min_memory;/* min_memory radio button         */
    Widget        iws_proportionalorfixed; /* fixed or proportional       */
    Widget        iws_radio_proportional; /* proportional radio button    */
    Widget        iws_radio_fixed;     /* fixed radio button              */   
    Widget        iws_radio_profile;     /* profile radio button          */
    Widget        radio_profile_full;
    Widget        radio_profile_sparse;
    Widget        radio_profile_no;

    int           Sankoff_matrix_num;
    int           run_done_flag=0;

    static char*  s_queryfile=" ";    /* input file name                  */
    static char*  s_datafile=" ";     /* input file name                  */
/* option strings */

   int   opradio_protein_dna;
    char  opstr_relation_matrix[MAX_STRING_LENGTH];
    char  opstr_number_of_sankoff_blocks[MAX_NUMBER_LENGTH];
    int   opcheck_output_file;
    char  opstr_output_file[MAX_FILE_NAME_LENGTH];
    int   opcheck_advanced_options;
    int   opradio_back_exact;
    char  opstr_sample_size[MAX_NUMBER_LENGTH];
    int   opradio_backprop_backfixed;
    int   opcheck_output_sample;
    int   opcheck_find_all;
    char  opstr_cutoff[MAX_NUMBER_LENGTH];
    int   opstr_base_blocksize[MAX_NUMBER_LENGTH];
    int   opcheck_special_site;
    char  opstr_special_site[MAX_STRING_LENGTH];
    char  opstr_special_weight[MAX_NUMBER_LENGTH];
    int   opcheck_best_align;
    int   opcheck_output_sequence;
    int   opradio_minmemory_maxspeed;
    int   opcheck_debug;
     
/* routine define here */
void GUI_Sankoff_accept(Widget, void*, void*);
void GUI_Sankoff_option_menu(Widget, void*, void*);
void GUI_Sankoff_stop(Widget, void*, void*);
void GUI_Sankoff_help(Widget, void*, void*);
void GUI_Sankoff_plot_marginal(Widget, void*, void*);
void GUI_Sankoff_plot_profile (Widget, void*, void*);
void GUI_Sankoff_plot_cell (Widget, void*, void*);
void GUI_Sankoff_toggle_true(Widget, void*, void*);
void GUI_Sankoff_toggle_false(Widget, void*, void*);
void GUI_Sankoff_toggle_radio_true(Widget, void*, void*);
void GUI_Sankoff_toggle_radio_false(Widget, void*, void*);
void nested_radio(Widget parent, void* data, void* call_data);
void GUI_Sankoff_off_hide(Widget parent, void* data, void* call_data);

void SelectFileCallback ( Widget    w, 
                          Widget cbWidget, 
                          XtPointer callData );
void CancelCallback ( Widget    w, 
                      XtPointer clientData, 
                      XtPointer callData );

void OKCallback ( Widget    w, 
                  Widget cbWidget, 
                  XtPointer callData );



static void ValueChangedCallback ( Widget    w, 
                                   XtPointer clientData, 
                                   XtPointer callData );



/********************   GUI_Sankoff   ******************************/

void GUI_Sankoff(Widget iw_temp, void* data, void* call_data) 
{
    Widget rowcol_tmp,rowcol_tmp1;   /* tmp row_col widget */
    Widget option_menu;
    char cval[50];

    Arg           args[10];         /* arg list */
    register      int  n;           /* arg count */
    /* set up default values for options */

    /*    strncpy(opstr_query_file,"",MAX_FILE_NAME_LENGTH-1);  
    strncpy(opstr_data_file,"",MAX_FILE_NAME_LENGTH-1);
    */
    
    

    if(iw_Sankoff!=(Widget)0) { /* already created */
    
        XtManageChild(iw_Sankoff);
	return;
    }
    
    /* otherwise create a menu widget first */        

     iw_Sankoff = wid_dialog(iw_temp,iw_Sankoff,"Bayesian Adaptive Alignment",
			    30,30);
    
   if(iws_rowcol!=(Widget)0) XtUnmanageChild(iws_rowcol); 
    iws_rowcol   = wid_rowcol(iw_Sankoff ,'v',-1,-1);
    iws_label=wid_labelg(iws_rowcol,0,
			 "Bayesian Adaptive Alignment for Protein and DNA"
			 ,0,0);

    /* widget for input file */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1); 
    sprintf(cval,"%s",s_queryfile);
    iws_queryfile=wid_textboxb(rowcol_tmp,0,"Query file:",cval,50);       
    buts_browse=wid_pushg(rowcol_tmp,0," Query  Browse ",
                         SelectFileCallback,iws_queryfile,0,0); 
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1); 
    sprintf(cval,"%s",s_datafile);
    iws_datafile=wid_textboxb(rowcol_tmp,0,"Data file:",cval,40); 
    buts_browse=wid_pushg(rowcol_tmp,0," Data Browse  ",
                         SelectFileCallback,iws_datafile,0,0);

    /* maximum number of blocks and relation matrix */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
  
    /*  iws_protein=wid_toggleg(rowcol_tmp,0,
			"Protein:",TRUE,
			GUI_Sankoff_toggle_true,NULL,-1,-1); */


    /***********************************************************/


  iws_alpha = XtVaCreateManagedWidget("A Simple Toggle Button",
			      xmRowColumnWidgetClass,
			      rowcol_tmp,
			      XmNentryAlignment, XmALIGNMENT_CENTER,
			      XmNorientation, XmHORIZONTAL,
			      XmNpacking, XmPACK_TIGHT,
			      XmNradioBehavior, True,
			      XmNnumColumns, 1,
			      NULL);



 
   /*
    * Create the children of the XmRowColumn widget.
    */ 

   
        iws_radio_protein = XtVaCreateManagedWidget ( "PROTEIN", 
                                           xmToggleButtonWidgetClass,
                                           iws_alpha,
                                           NULL);
                                          

        XtAddCallback ( iws_radio_protein, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
      
        iws_radio_dna = 
                      XtVaCreateManagedWidget ( "DNA", 
                                              xmToggleButtonWidgetClass,
                                              iws_alpha,

                              NULL);
                                          

        XtAddCallback ( iws_radio_dna, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
    
     XmToggleButtonSetState(iws_radio_protein,TRUE,FALSE);
	/***************************************************************/

    Sankoff_matrix_num=5; 
    rowcol_tmp1=wid_rowcol(rowcol_tmp,'h',-1,-1);
    option_menu=XmVaCreateSimpleOptionMenu(rowcol_tmp1,"option_menu",
	   XmStringCreateLocalized("Relation matrix:"),'R',Sankoff_matrix_num,
					        GUI_Sankoff_option_menu,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum30"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum35"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum40"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum45"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum50"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum62"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum80"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("Blosum100"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("All BLOSUM"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM40"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM50"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM60"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM70"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM80"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM90"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM100"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM110"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM120"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM130"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM140"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM150"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM160"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM170"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM180"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM190"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM200"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM210"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM220"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM230"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM240"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM250"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM260"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM270"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM280"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM290"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("PAM300"),' ',NULL,NULL,
	   XmVaPUSHBUTTON,XmStringCreateLocalized("All PAM"),' ',NULL,NULL,
	   NULL);
    XtManageChild(option_menu);

    /* score matrix option for DNA */
    sprintf(cval,"%d-%d",40,100);
    iws_dna_matrix=wid_textboxb(rowcol_tmp,0,"Relation matrix:PAM",
				    cval,8);
    sprintf(cval,"%d",10);
    iws_dna_stepsize=wid_textboxb(XtParent(iws_dna_matrix),0,
			      "incremental stepsize:",cval,3);
    XtUnmanageChild(XtParent(iws_dna_matrix));

    XtRemoveAllCallbacks(iws_radio_protein, XmNvalueChangedCallback);
    XtAddCallback(iws_radio_protein,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_true,option_menu);
    XtAddCallback(iws_radio_protein,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_false,iws_dna_matrix);


    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
    sprintf(cval,"%d",DEFAULT_NUM_BLOCKS);
    iws_max_blocks=wid_textboxb(rowcol_tmp,0,"Maximum number of blocks:",
				    cval,3);
  
    /* output file */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
    iws_outfile=wid_toggleg(rowcol_tmp,0,
			"Output file?",FALSE,
			GUI_Sankoff_toggle_true,NULL,-1,-1);
     sprintf(cval,"%s"," ");
     iws_outfile_name=wid_textboxb(rowcol_tmp,0,
				  "Output file name:",cval,30);  
     XtUnmanageChild(XtParent(iws_outfile_name));
     XtRemoveAllCallbacks(iws_outfile, XmNvalueChangedCallback);
     XtAddCallback(iws_outfile,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_outfile_name);


     /* advanced options */
   rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
    iws_advanced=wid_toggleg(rowcol_tmp,0, "Advanced Options",FALSE,
			  GUI_Sankoff_toggle_true,NULL,-1,-1);



    /* exact posterior or backsampling */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
     /****************************************************************/
  iws_exactorback = XtVaCreateManagedWidget("A Simple Toggle Button",
			      xmRowColumnWidgetClass,
			      rowcol_tmp,
			      XmNentryAlignment, XmALIGNMENT_CENTER,
			      XmNorientation, XmHORIZONTAL,
			      XmNpacking, XmPACK_TIGHT,
			      XmNradioBehavior, True,
			      XmNnumColumns, 1,
			      NULL);



 
   /*
    * Create the children of the XmRowColumn widget.
    */ 

   
        iws_radio_exact = XtVaCreateManagedWidget ( "Exact Posterior Probability)", 
                                           xmToggleButtonWidgetClass,
                                           iws_exactorback,
                                           NULL);
                                          

        XtAddCallback ( iws_radio_exact, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
      
        iws_radio_back = 
                      XtVaCreateManagedWidget ( "Back Sampling", 
                                              xmToggleButtonWidgetClass,
                                              iws_exactorback,

                              NULL);
                                          

        XtAddCallback ( iws_radio_back, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
    
     XmToggleButtonSetState(iws_radio_back,TRUE,FALSE);
	/***************************************************************/

    /* back sampling */
     sprintf(cval,"%s","1000");
    
     iws_backsampling_number=wid_textboxb(rowcol_tmp,0,
				  "Sample size:",cval,5);  

       /* backsampling proportional  option */
   /***********************************************************/
 rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);

  iws_proportionalorfixed = XtVaCreateManagedWidget("A Simple Toggle Button",
			      xmRowColumnWidgetClass,
			      rowcol_tmp,
			      XmNentryAlignment, XmALIGNMENT_CENTER,
			      XmNorientation, XmHORIZONTAL,
			      XmNpacking, XmPACK_TIGHT,
			      XmNradioBehavior, True,
			      XmNnumColumns, 1,
			      NULL);



 
   /*
    * Create the children of the XmRowColumn widget.
    */ 

   
        iws_radio_proportional = XtVaCreateManagedWidget ( "Back Sample proportional", 
                                           xmToggleButtonWidgetClass,
                                           iws_proportionalorfixed,
                                           NULL);
                                          

        XtAddCallback (iws_radio_proportional , XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
      
        iws_radio_fixed = 
                      XtVaCreateManagedWidget ( "Back Sample Fixed", 
                                              xmToggleButtonWidgetClass,
                                              iws_proportionalorfixed,

                              NULL);
                                          

        XtAddCallback ( iws_radio_fixed, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
    
     XmToggleButtonSetState(iws_radio_proportional,TRUE,FALSE);
	/***************************************************************/

     rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
     iws_backsampling_sequence=wid_toggleg(rowcol_tmp,0,
			"Output sampling sequence?",FALSE,
			GUI_Sankoff_toggle_true,NULL,-1,-1);
   
     /* trace back all the alignment above a cutoff */
     rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
     iws_backsampling_cutoff=wid_toggleg(rowcol_tmp,0,
			"Find all alignment above a cutoff?",FALSE,
			GUI_Sankoff_toggle_true,NULL,-1,-1);
     sprintf(cval,"%s"," ");
     iws_backsampling_cutoff_value=wid_textboxb(rowcol_tmp,0,
				  "Cutoff value:",cval,5);  
      XtUnmanageChild(XtParent(iws_backsampling_cutoff_value)); 

     XtRemoveAllCallbacks(iws_backsampling_cutoff, XmNvalueChangedCallback);
     XtAddCallback(iws_backsampling_cutoff,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_backsampling_cutoff_value);
  

     XtRemoveAllCallbacks(iws_radio_back, XmNvalueChangedCallback);

     XtAddCallback(iws_radio_back,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_true,iws_proportionalorfixed); 
      XtAddCallback(iws_radio_back,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_true,iws_backsampling_number);
     XtAddCallback(iws_radio_back,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_true,iws_backsampling_sequence);
     XtAddCallback(iws_radio_back,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_radio_true,iws_backsampling_cutoff);

     rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1); 
     sprintf(cval,"%d",BASE_BLOCKSIZE);
     iws_blocks=wid_textboxb(rowcol_tmp,0,"Base block size:",
				    cval,3);

 iws_radio_profile = XtVaCreateManagedWidget("A Simple Toggle Button",
			      xmRowColumnWidgetClass,
			      rowcol_tmp,
			      XmNentryAlignment, XmALIGNMENT_CENTER,
			      XmNorientation, XmHORIZONTAL,
			      XmNpacking, XmPACK_TIGHT,
			      XmNradioBehavior, True,
			      XmNnumColumns, 1,
			      NULL);

      radio_profile_full = XtCreateManagedWidget( "Full Profile",
						   xmToggleButtonWidgetClass,
						   iws_radio_profile,NULL,0);

   radio_profile_sparse = XtCreateManagedWidget( "Sparse Profile",
						   xmToggleButtonWidgetClass,
						   iws_radio_profile,NULL,0);
 XmToggleButtonSetState(radio_profile_sparse,TRUE,FALSE);

    radio_profile_no = XtCreateManagedWidget( "No Profile",
						   xmToggleButtonWidgetClass,
						   iws_radio_profile,NULL,0);
  

   XtUnmanageChild(XtParent(radio_profile_sparse));
     /*********************************************************************************/

    /* special sites */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
    iws_special_site=wid_toggleg(rowcol_tmp,0,
			"Special sites:",FALSE,
			GUI_Sankoff_toggle_true,NULL,-1,-1);
    sprintf(cval,"%s"," ");
    iws_sites=wid_textboxb(rowcol_tmp,0," ",cval,30);  
    XtUnmanageChild(XtParent(iws_sites)); 
    iws_site_weight=wid_textboxb(rowcol_tmp,0,"weight:",cval,5);  
    XtUnmanageChild(XtParent(iws_site_weight));
    XtRemoveAllCallbacks(iws_special_site, XmNvalueChangedCallback);
    XtAddCallback(iws_special_site,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_sites);
    XtAddCallback(iws_special_site,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_site_weight);
 
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);
    iws_best_alignment=wid_toggleg(rowcol_tmp,0,
				    "Best alignment?",FALSE,
				    GUI_Sankoff_toggle_true,NULL,-1,-1);
    iws_output_sequence=wid_toggleg(rowcol_tmp,0,
				    "Output matching sequence?",FALSE,
				    GUI_Sankoff_toggle_true,NULL,-1,-1);


    /*************** min memory or speed **************************/

    /***********************************************************/
 rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1);

  iws_speedormemory = XtVaCreateManagedWidget("A Simple Toggle Button",
			      xmRowColumnWidgetClass,
			      rowcol_tmp,
			      XmNentryAlignment, XmALIGNMENT_CENTER,
			      XmNorientation, XmHORIZONTAL,
			      XmNpacking, XmPACK_TIGHT,
			      XmNradioBehavior, True,
			      XmNnumColumns, 1,
			      NULL);



 
   /*
    * Create the children of the XmRowColumn widget.
    */ 

   
        iws_radio_min_memory = XtVaCreateManagedWidget ( "Minimum Memory", 
                                           xmToggleButtonWidgetClass,
                                           iws_speedormemory,
                                           NULL);
                                          

        XtAddCallback (iws_radio_min_memory , XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
      
        iws_radio_speed = 
                      XtVaCreateManagedWidget ( "Maximum Speed", 
                                              xmToggleButtonWidgetClass,
                                              iws_speedormemory,

                              NULL);
                                          

        XtAddCallback ( iws_radio_speed, XmNvalueChangedCallback, 
                        ValueChangedCallback, NULL );
    
     XmToggleButtonSetState(iws_radio_min_memory,TRUE,FALSE);
	/***************************************************************/

    /* debug information */
    iws_debug=wid_toggleg(rowcol_tmp,0, "Debug",FALSE,
			  GUI_Sankoff_toggle_true,NULL,-1,-1);

     /* unmanage for the advance options */

 XtRemoveAllCallbacks(iws_advanced, XmNvalueChangedCallback);
 
  XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_exactorback);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_output_sequence);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_best_alignment);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_special_site);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_radio_min_memory); 
XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,radio_profile_sparse); 

 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_debug); 
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_toggle_true,iws_blocks); 

 /* hide back sampling options */


 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_off_hide,iws_proportionalorfixed); 
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_off_hide,iws_backsampling_number);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_off_hide,iws_backsampling_sequence);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_off_hide,iws_backsampling_cutoff);
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,
		   GUI_Sankoff_off_hide,iws_blocks); 

XtCallCallbacks(iws_advanced,XmNvalueChangedCallback,NULL); /* hide it all */
 /* show back sampling options */
 XtAddCallback(iws_advanced,XmNvalueChangedCallback,nested_radio,NULL);
 

    /* button for stop and accept */
    rowcol_tmp=wid_rowcol(iws_rowcol,'h',-1,-1); 
   
    buts_accept=wid_pushg(rowcol_tmp,0,"  Accept  ",
                         GUI_Sankoff_accept,NULL,0,0);
    buts_plot_profile=wid_pushg(rowcol_tmp,0," Histo-plot profile  ",
                         GUI_Sankoff_plot_profile,NULL,0,0);
    buts_plot_cell=wid_pushg(rowcol_tmp,0,"  Cell-plot profile  ",
                         GUI_Sankoff_plot_cell,NULL,0,0);
    buts_plot_marginal=wid_pushg(rowcol_tmp,0,"  Plot marginal  ",
                         GUI_Sankoff_plot_marginal,NULL,0,0);  
    buts_help=wid_pushg(rowcol_tmp,0,"  Help  ",
                         GUI_Sankoff_help,NULL,0,0);
    buts_stop=wid_pushg(rowcol_tmp,0,"   Exit  ",
                         GUI_Sankoff_stop,NULL,0,0);
    XtManageChild(iw_Sankoff);
    /*  XtAppMainLoop(app_context); */
}

/************ accept button callback *********************************/

void GUI_Sankoff_accept(Widget iw_temp, void* data, void* call_data)
{
  /*  Model M;
  Sank_S SK;*/
  int i,j;
  int backsampling, backsampling_blocks,backsampling_number;
  int Sankoff_max_blocks;
  int backsampling_cutoff;
  int  exact_posterior;
  float backsampling_cutoff_value;
  int Sankoff_blocksize;
  int val;
  int protein_flag;           /* flag for protein or DNA sequences */
  int matrix_start,matrix_end;/* start and end matrix for DNA      */
  int matrix_increment_step;  /* DNA matrix incremental step       */
  char *cval;
  char input[20];
  char *outputfile;
  FILE *query_fptr,*data_fptr;
  int int_number;
  double double_number;
  /* char command_string[MAX_STRING_LENGTH]; */
  char *opstr_query_file;
  char *opstr_data_file;
  char *opstr_base_blocksize;
  char *opstr_number_of_sankoff_blocks;
  char option_string[MAX_STRING_LENGTH];
  void build_command(char *command_string,char *option_string);
  FILE *fp;

  NInput_N=0;       /* no special sites */
  Input_alpha=0;    /* no special site weight */
  /* initialize command string nulls */ 
  strncpy(command_string,"",1);
  sargc=1;
  sargv[0]="bayes_aligner";
  /* sargv[1]="-PBayes_align";
  sargc++;*/
  strncpy(option_string,"",MAX_STRING_LENGTH-1);
  

  /* input data file */
  opstr_query_file=XmTextGetString(iws_queryfile);
   GUI_strip_space(opstr_query_file); 
  if((query_fptr=fopen(opstr_query_file,"r"))==NULL)
    {
      wid_error(iw_Sankoff,0,"Check query file name!!!",0,0);
      return;
    }
  fclose(query_fptr);
  /*  printf("just before the build command\n");
  printf("opstr_query_file= %s\n",opstr_query_file); */
  sprintf(option_string,"-QUERY=%s",opstr_query_file);
  build_command(command_string,option_string);
  

  opstr_data_file=XmTextGetString(iws_datafile);
  GUI_strip_space(opstr_data_file);
  if((data_fptr=fopen(opstr_data_file,"r"))==NULL)
  {
      wid_error(iw_Sankoff,0,"Check data file name!!!",0,0);
      return;
  }
  fclose(data_fptr);
  sprintf(option_string,"-DATA=%s",opstr_data_file);
  build_command(command_string,option_string);
  /* check base Block size */

  opstr_base_blocksize=XmTextGetString(iws_blocks);
  GUI_strip_space(opstr_base_blocksize); 
  if(sscanf(opstr_base_blocksize,"%d", &int_number)!= 1)
  { 
      wid_error(iw_Sankoff,0,"Check input of base block size",
		0,0);
      return;        
  }
  sprintf(option_string,"-BASE_BLOCKSIZE=%s",opstr_base_blocksize);
  build_command(command_string,option_string); 
  /* check maximum number of blocks */

  opstr_number_of_sankoff_blocks=XmTextGetString(iws_max_blocks);
  GUI_strip_space(opstr_number_of_sankoff_blocks);
  if(sscanf(opstr_number_of_sankoff_blocks,"%d", &int_number)!= 1)
  { 
      wid_error(iw_Sankoff,0,"Check input of maximum number of blocks",
		0,0);
      return;        
  }    
  if (int_number > 200 ||   int_number < 1)
  { 
      wid_error(iw_Sankoff,0,"Check input of maximum number of blocks",
		0,0);
      return;        
   }
 sprintf(option_string,"-MAX_NUM_BLOCKS=%s",opstr_number_of_sankoff_blocks);
  build_command(command_string,option_string);

  /* check protein or DNA sequence flag */
 if(XmToggleButtonGetState(iws_radio_protein)) {
     sprintf(cval,"%d",Sankoff_matrix_num);
     sprintf(option_string,"-MATRIX_NUM=%s",cval);
     build_command(command_string,option_string);
      protein_flag=1;
  }
  else{ /*  DNA sequence */
      protein_flag=0;
      build_command(command_string," -NON_PROTEIN");
      /* check input for PAM indeces */
      cval=XmTextGetString(iws_dna_matrix);      
      GUI_strip_space(cval);
      if(strlen(cval)-1==0){
	  wid_error(iw_Sankoff,0, 
		    "Input matrices as start-end (e.g. 40-100)",0,0);
	  return;        
      }
      /* scan input */
      j=0;
      for(i=0;i<=strlen(cval);i++){
	  if(cval[i]=='-'){
	      /* add deliminator */
	      input[j]=' ';
	  }
	  else{
	      input[j]=cval[i];
	  }
	  j++;
      }
      sscanf(input,"%d %d", &matrix_start,&matrix_end);
      if(matrix_start<=0 || matrix_start>matrix_end) {
	  wid_error(iw_Sankoff,0, 
		    "Input matrices as start-end (e.g. 40-100)\nstart >0, and start<=end",0,0);
	  return;        
      }

      cval=XmTextGetString(iws_dna_stepsize);      
      GUI_strip_space(cval);
      if(strlen(cval)==0) {
	  wid_error(iw_Sankoff,0, 
		    "check incremental step!!",0,0);
	  return;        
      }
      sscanf(cval,"%d", &matrix_increment_step);
      if(matrix_increment_step<=0) {
	  wid_error(iw_Sankoff,0, 
		    "Incremental step should be greater than zero",0,0);
	  return;        
      }
      sprintf(cval,"%d-%d,%d",matrix_start,matrix_end,matrix_increment_step);  
      sprintf(option_string,"-DNA_MATRIX=%s",cval);
      build_command(command_string,option_string);
             
  }
  

  /* check special sites */
  if(XmToggleButtonGadgetGetState(iws_special_site)){
      cval=XmTextGetString(iws_sites);
      GUI_strip_space(cval);
      if(strlen(cval)-1==0){
	  wid_error(iw_Sankoff,0, "Input special sites",0,0);
	  return;        
      }
      /* scan sites */
      NInput_N=0;
      j=0;
      for(i=0;i<=strlen(cval);i++){
	  if(cval[i]==',' || cval[i]==NULL){
	      /* add terminator the end of option string */
	      input[j]=NULL;
	      sscanf(input,"%d", &Input_N[NInput_N]);
	      NInput_N++;
	      j=0;
	  }
	  else{
	      input[j]=cval[i];
	      j++;
	  }
      }
      printf("NInput_N=%d\n",NInput_N);
       build_command(command_string," -IMPORT_SITE=");
       build_command(command_string,input);
       sprintf(option_string,"-IMPORT_SITE=%s",input);
       build_command(command_string,option_string);
      
      cval=XmTextGetString(iws_site_weight);
      GUI_strip_space(cval);
      if(sscanf(cval,"%f", &Input_alpha)!=1){
	wid_error(iw_Sankoff,0, "Check special site weight",0,0);
	  return;        
      }       
     sprintf(option_string,"-SITE_WEIGHT=%s",cval);
     build_command(command_string,option_string);
       
      
  }   
  /* check back sampling */
 if(XmToggleButtonGetState(iws_radio_back)) 
    {
      backsampling=TRUE;
      exact_posterior=FALSE;    	  

      /* check number of sequence for back sampling */
      cval=XmTextGetString(iws_backsampling_number);
      GUI_strip_space(cval); 
      if(sscanf(cval,"%d", &backsampling_number)!= 1)
	{ 
	  wid_error(iw_Sankoff,0,
		    "Check number of sequence for back sampling!!",0,0);
	  return;        
	} 
      sprintf(option_string,"-BACK_S_SIZE=%s",cval);
      build_command(command_string,option_string);
       
             
      /* check trace back all alignment above a cutoff */
      if(XmToggleButtonGadgetGetState(iws_backsampling_cutoff))
	{
	  backsampling_cutoff=TRUE;
          build_command(command_string,"-BACK_CUTOFF");
	  cval=XmTextGetString(iws_backsampling_cutoff_value);
	  GUI_strip_space(cval); 
	  if(sscanf(cval,"%f", &backsampling_cutoff_value)!= 1)
	    { 
	      wid_error(iw_Sankoff,0,
			"Check cutoff value for back sampling!!",0,0);
	      return;        
	    } 
          sprintf(option_string,"-BACK_C_VALUE=%s",cval);
          build_command(command_string,option_string);
        
	}       
      else
	backsampling_cutoff=FALSE;
    }
  else
    {
      backsampling=FALSE;
      exact_posterior=TRUE;
      if(XmToggleButtonGetState(iws_radio_min_memory)) {
           wid_error(iw_Sankoff,0,
	   "Can not use the MINUMIM MEMORY option with \n EXACT POSTERIOR PROBABILITY option!!",0,0);
	   return;
      }
  
      build_command(command_string,"-EXACT_POST" );  
 
    }

  /* output file */
 if(XmToggleButtonGadgetGetState(iws_outfile))
    {
      outputfile=XmTextGetString(iws_outfile_name);
      GUI_strip_space(outputfile); 
      if (strlen(outputfile) < 1) {
          wid_error(iw_Sankoff,0,
	   "Error, Unspecified output file!!",0,0);
	   return;
      }

      sprintf(option_string,"-OUTFILE=%s",outputfile);
      build_command(command_string,option_string);
    
    }

  /* flags for best alignment and output alignment sequences */
  if(XmToggleButtonGadgetGetState(iws_best_alignment)){
     if(XmToggleButtonGetState(iws_radio_min_memory)) {
           wid_error(iw_Sankoff,0,
	   "Can not use the MINUMIM MEMORY option with \n BEST ALIGNMENT option!!",0,0);
	   return;
      }
     build_command(command_string, "-BEST_ALIGN" );     
  }
 
  if(XmToggleButtonGadgetGetState(iws_output_sequence)){   
      build_command(command_string,"-OUT_SEQ" );            
  }
  if(!XmToggleButtonGetState(iws_radio_min_memory)) {
     build_command(command_string,"-MAX_SPEED" );
  }
 if(XmToggleButtonGetState(radio_profile_full)){   
      build_command(command_string,"-FULL_PROFILE" );          
  }
 if(XmToggleButtonGetState(radio_profile_no)){   
      build_command(command_string,"-NO_PROFILE" );          
  }
  if(XmToggleButtonGadgetGetState(iws_debug)){   
      build_command(command_string,"-DEBUG" );          
  }
  /*   XtUnmanageChild(iw_Sankoff);
    XtUnmanageChild(iw_temp); */
  XFlush(idispl);

  /* Sankoff(M,SK); */
   printf("command string = %s\n",command_string); 
  for (i=0;i<sargc;i++)
    printf("%s ",sargv[i]);
  printf("\n");
 DisplayBusyCursor (toplevel ); 
 
 if (fp=fopen("Bayes_aligner.cmd","w")) { 
  for (i=0;i<sargc;i++)
    fprintf(fp,"%s ",sargv[i]);
  fprintf(fp,"\n");
  /*  fprintf(fp,"%s\n",command_string); */
    fclose(fp);
 }
 else {
    printf("unable to open Bayes_aligner.cmd\n");
 }
 Gibbs_NON_GUI(sargc,sargv);   
 RemoveBusyCursor ( toplevel);
 run_done_flag=1; /* I have done at least one run so I can plot */ 
  
}


/*************** option menu call back for relation matrix ************************/
void GUI_Sankoff_option_menu(Widget parent, void* data, void* call_data)
{

  Sankoff_matrix_num=(int) data;

}


/*************** toggle button click *********************************/
void GUI_Sankoff_toggle_true(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;

  if(tmp!=NULL){
      if(XmToggleButtonGadgetGetState(parent)){
	  /* pop up the following widget */
	  XtManageChild(XtParent(tmp));
           XtUnmanageChild(XtParent(XtParent(tmp)));
           XtManageChild(XtParent(XtParent(tmp)));
            
      }
      else{
	  /* pop down the following widget */
	 XtUnmanageChild(XtParent(tmp)); 
          XtUnmanageChild(XtParent(XtParent(tmp)));
          XtManageChild(XtParent(XtParent(tmp)));
    
      }
  }
}



/*************** hide button click *********************************/
void GUI_Sankoff_off_hide(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;

  if(tmp!=NULL){
      if(!XmToggleButtonGadgetGetState(parent)){
	  /* pop up the following widget */
	  XtUnmanageChild(XtParent(tmp));
      } else {
       XtCallCallbacks(iws_radio_proportional,XmNvalueChangedCallback,call_data);
      }
 
  }
}

void GUI_Sankoff_toggle_false(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;

  if(tmp!=NULL){
      if(!XmToggleButtonGadgetGetState(parent)){
	  /* pop up the following widget */
	  XtManageChild(XtParent(tmp));
      }
      else{
	  /* pop down the following widget */
	  XtUnmanageChild(XtParent(tmp));
      }
  }
}


/************ stop button callback *********************************/

void GUI_Sankoff_stop(Widget iw_temp, void* data, void* call_data)
{
      XtUnmanageChild(iw_Sankoff); 
      exit(0); 
}



/************ plot_profile button callback *********************************/

void GUI_Sankoff_plot_profile(Widget iw_temp, void* data, void* call_data)
{ 
      int return_code,maxi,maxj;
      FILE *fptr;
      char line[100],command[200],filename[100], xlabel[100],ylabel[100];
      if(!run_done_flag) {
           wid_error(iw_Sankoff,0,
	   "No current plot_profile file\n You must first run an alignment\n(press Accept to an alignment)",0,0);
      }
      else {
      DisplayBusyCursor (toplevel ); 
      if ((fptr=fopen("bayes.info","r"))==NULL) {
         fprintf(stderr,"unable to open bayes.info for read\n");
         return;
      }
      fgets(line,99,fptr);
      sscanf(line,"%d",&maxi);
      fgets(line,99,fptr);
      sscanf(line,"%d",&maxj);
      fgets(line,99,fptr); /* read past the alignment filename */
      fgets(line,99,fptr); /* get profile file name */
      line[strlen(line)-1]='\0';
      strcpy(filename,line);
      fgets(line,99,fptr);
      line[strlen(line)-1]='\0'; 
      strcpy(xlabel,line);
      fgets(line,99,fptr);
      line[strlen(line)-1]='\0'; 
      strcpy(ylabel,line);
      
      fclose(fptr);

      sprintf(command,"histodisp -FILENAME=%s -MAXI=%d -MAXJ=%d -LABELX=%s -LABELY=%s -LABELZ='Probability of alignment'",
          filename,maxi,maxj,xlabel,ylabel);
      if (return_code=system(command)) {
          printf("problem with histodisp command line=%s , return code =%d\n",command, return_code); 
       wid_error(iw_Sankoff,0,
	   "Unable to plot\n Please check to see if histodisp is on your path",0,0);
      }
      RemoveBusyCursor ( toplevel); 
      }
}


/************ plot_profile button callback *********************************/

void GUI_Sankoff_plot_cell(Widget iw_temp, void* data, void* call_data)
{ 
      int return_code,maxi,maxj;
      FILE *fptr;
      char line[100],command[200],filename[100], xlabel[100],ylabel[100];
      if(!run_done_flag) {
           wid_error(iw_Sankoff,0,
	   "No current plot_profile file\n You must first run an alignment\n(press Accept to an alignment)",0,0);
      }
      else {
      DisplayBusyCursor (toplevel ); 
      if ((fptr=fopen("bayes.info","r"))==NULL) {
         fprintf(stderr,"unable to open bayes.info for read\n");
         return;
      }
      fgets(line,99,fptr);
      sscanf(line,"%d",&maxi);
      fgets(line,99,fptr);
      sscanf(line,"%d",&maxj);
      fgets(line,99,fptr); /* read past the alignment filename */
      fgets(line,99,fptr); /* get profile file name */
      line[strlen(line)-1]='\0';
      strcpy(filename,line);
      fgets(line,99,fptr);
      line[strlen(line)-1]='\0'; 
      strcpy(xlabel,line);
      fgets(line,99,fptr);
      line[strlen(line)-1]='\0'; 
      strcpy(ylabel,line);
      
      fclose(fptr);

      sprintf(command,"celldisp -FILENAME=%s -MAXI=%d -MAXJ=%d -LABELX=%s -LABELY=%s",
          filename,maxi,maxj,xlabel,ylabel);
      if (return_code=system(command)) {
          printf("problem with histodisp command line=%s , return code =%d\n",command, return_code); 
       wid_error(iw_Sankoff,0,
	   "Unable to plot\n Please check to see if histodisp is on your path",0,0);
      }
      RemoveBusyCursor ( toplevel); 
      }
}

/************ plot_marginal button callback *********************************/

void GUI_Sankoff_plot_marginal(Widget iw_temp, void* data, void* call_data)
{ 
      int return_code;
      if(!run_done_flag) {
           wid_error(iw_Sankoff,0,
	   "No current plot_marginal file\n You must first run an alignment\n(press Accept to an alignment)",0,0);
      }
      else {
      DisplayBusyCursor (toplevel ); 
      if (return_code=system("sciplot marginal.txt")) {
          printf("problem with sciplot, return code =%d\n",return_code); 
       wid_error(iw_Sankoff,0,
	   "Unable to plot\n Please check to see if sciplot is on your path",0,0);
      }
      RemoveBusyCursor ( toplevel); 
      }
}

/**********************help button callback **********************/
void GUI_Sankoff_help(Widget iw_temp, void* data, void* call_data)
{
     Widget helpwindow, text;
     Arg        args[10];    /*  arg list        */ 
     register int    n;        /*  arg count        */ 
 
     n=0;
     XtSetArg (args[n], XmNwidth, 400); n++;
     XtSetArg (args[n], XmNheight, 500); n++;
     XtSetArg (args[n], XtNlabel, "Help"); n++;
     
    helpshell = XtCreatePopupShell ( "help", topLevelShellWidgetClass,
                                  toplevel, args, n );
     XtPopup ( helpshell, XtGrabNone );
 
     n=0;
     XtSetArg (args[n], XmNwidth, 200); n++;
     XtSetArg (args[n], XmNheight, 100); n++;
     XtSetArg (args[n], XmNeditMode, XmMULTI_LINE_EDIT); n++; 
     XtSetArg(args[n], XmNeditable,False);n++;     
     text = XmCreateScrolledText ( helpshell, "text", args, n );
        XtManageChild (text ); 
    
XmTextSetString(text,"\
INTRODUCTION:\n\
   This is the distribution package for the Bayesian Aligner.\n\
   The detail theory of Bayes aligner can be found in\n\
\n\
   (1) Bayesian Adaptive Alignment and Inference\n\
       Jun Zhu, Jun S. Liu and Charles E. Lawrence\n\
       ISMB-97 pp.358-368.\n\
   (2) Bayesian adaptive sequence alignment algorithm\n\
       Jun Zhu, Jun S. Liu and Charles E. Lawrence\n\
       Bioinformatics (in press) 2/98.\n\
\n\
INSTALLATION:\n\
   The package is a tar.gz file.\n\
   The directory structure of the tar file is:\n\
\n\
   Bayesaligner/     Top directory\n\
      examples/     data and command files for examples\n\
      source/       c source code and makefiles\n\
      unix/\n\
         bin/       solaris binary executable\n\
         tools/     binary display tools for unix\n\
      win95/\n\
         bin/       win95 binary executable\n\
         tools/     binary display tools for win95\n\
\n\
   To install:\n\
      1) copy the distribution file(Bayesaligner.tar.gz) \n\
        to your favorite installation directory.\n\
      2) cd to your favorite installation directory.\n\
      3) gunzip Bayesaligner.tar.gz\n\
      4) tar -xvf  Bayesaligner\n\
     \n\
USING THE BAYESALIGNER:\n\
   There are two ways to run the program, \n\
      1) graphic user interface (GUI)(on unix only) or\n\
      2) plain command line\n\
\n\
    By default, it starts in GUI mode if no command line options are given.\n\
    The GUI mode provides  a user friendly interface for specifying options.\n\
    (An equivalent NON-GUI command line will be generated and placed in a \n\
     file name \"Bayes_aligner.cmd\".)\n\
\n\
    Make sure that the Bayesaligner bin/ and tools/ directories are on your path. \n\
\n\
    The query and data files must contain sequences in FASTA format. \n\
 \n\
    There are multiple relation matrices (PAM and BLOSUM) from which you can choose.\n\
\n\
COMMAND LINE OPTIONS:\n\
    -QUERY=: filename of sequence to align.\n\
    -DATA=: The filename of the 2nd sequence to align.\n\
    -NON_PROTEIN: switch between protein/DNA (Protein is the default).\n\
    -MAX_NUM_BLOCKS=: The number of Sankoff blocks to align\n\
    -DNA_MATRIX,STEP: DNA relation matrix number: [1-500] PAM,step (ie use a range).\n\
    -MATRIX_NUM: protein relation matrix number:[0-7] BLOSUM30-100, [8] ALL BLOSUM, \n\
                  [9-35] PAM40-300, step 10,[36] ALL PAM.\n\
\n\
     Advanced options: (These options should rarely be used, please send us email \n\
                         (steve@wadsworth.org) if you plan to use them)\n\
         -DEBUG=: gives you detail information.\n\
         -OUTFILE=: The output filename (Default is just the screen).\n\
         -BEST_ALIGN: a straight forward Sankoff method.\n\
         -BASE_BLOCKSIZE=: the base size of a block used in calculation of the\n\
                           number of Sankoff blocks.(Advanced).\n\
         -EXACT_BLOCK: Don't backsample proportional (Default is false).\n\
         -EXACT_POST: Calculate the exact posterior probability (Default is false).\n\
         -OUT_SEQ:output matching sequence if best alignments are calculated.\n\
         -NO_PROFILE: Don't create a full(AKA big) .profile file .\n\
         -FULL_PROFILE:  Create a full(AKA big) .profile file (Default is a sparse file) .\n\
         -BACK_S_SIZE=: The number of backsamples.\n\
         -BACK_S_OUTSEQ: Output all alignments\n\
         -BACK_CUTOFF: find all alignments above a cutoff value.\n\
         -BACK_C_VALUE=: Cutoff value\n\
         -BACK_S_BLOCKS=: find the [N] most likely alignments (Default is 20)\n\
\n\
ANOTATED EXAMPLE: \n\
     Bayesaligner  -QUERY=50A.human.aldolas.pM.fta -DATA=50C.mouse.j05517.same.fta  \n\
                   -OUTFILE=demo.out -NON_PROTEIN -DNA_MATRIX=1-1,1 -MAX_NUM_BLOCKS=20\n\
                   -BACK_S_SIZE=1000  \n\
\n\
     A brief description of the above command line:\n\
          -QUERY=50A.human.aldolas.pM.fta -DATA=50C.mouse.j05517.same.fta :\n\
           The two fasta files containing the sequences\n\
\n\
           -OUTFILE=demo.out:  Send the most of the screen output to demo.out\n\
\n\
           -NON_PROTEIN:  This is a DNA sequence\n\
\n\
           -DNA_MATRIX=1-1,1:\n\
                 use just the PAM1 maxtrix. Another example 10-30,10,\n\
                 would use PAM10 through PAM30 in increments of 10 \n\
                 (so it would use 3 matricies PAM10,PAM20,PAM30).\n\
           -BACK_S_SIZE=1000:\n\
                   backsample 1000 times (a 1000 will give you enough backsamples \n\
                   to \"smooth out\" the results of this stochatisic process)\n\
\n\
UNDERSTANDING THE RESULTS:\n\
     The Bayes_aligner is a \"statistical\" procedure that returns the entire\n\
     alignment space (not just the maximum). The products of this program are \n\
     the posterior probabilities of alignment variables and their visualizations.\n\
     P(k|R)---posterior probability of k blocks given this pair of sequences.\n\
     P(matrix|R)---posterior probability that this relation matrix correctly\n\
             relates these sequences.\n\
\n\
     The posterior probability of each pair of bases aligning can be obtained using\n\
     back sampling or by exact calculation.  Backsampling is the default method\n\
     because it uses less memory and less compute time . Backsampling will generate\n\
     three files ****.profile.*, ***.marginal.* and ****.alignment.*. The posterior \n\
     probability of each base aligning (***.profile) can be viewed using gnuplot(proview)\n\
     or the GUI Histo-plot button. This is a 3D plot, the \"peaks\" represent the locations\n\
     in each sequence where there is a strong alignment. The marginal posterior probability\n\
     (***.marginal.*) of a particular base in the query sequence aligning with any base \n\
     in the data sequence can be view by gnuplot(margview) or the GUI Plot marginal button.\n\
     The marginal plot is useful in determining what are the conserved regions of the\n\
     query sequence in the data sequence(Phylogenetic footprinting). The alignment\n\
     (****.alignment.*) displays the most likely aligned blocks between the query and\n\
     data sequences.    \n\
\n\
TOOLS:\n\
     The GUI tools associated with the buttons are avaiable only on X-window/unix systems.\n\
     The gnuplot tools(margview and proview) are platform independent. \n\
          How to use margview:\n\
              1) type    gnuplot margview\n\
          How to use proview:\n\
              1) copy your (*****.profile.*) file to profile.txt\n\
              2) type    gnuplot proview\n\
\n\
WARNINGS:\n\
     Short sequence(<20) or identical sequences tend to cause underflow and/or overflow\n\
     problems. Please avoid these type of trivial sequence alignments.\n\
\n\
MAKING NEW BINARIES:\n\
     How to make new binaries (this should not be necessary):\n\
         1)To compile this on your unix machine , just type \"make\" in the source directory.\n\
         2)To compile this on your win95, just type makeit. \n\
     The GUI version of this program requires X windows and Motif(or Lesstif) libraries.\n\
\n\
\n\
\n"); 
}





static void ValueChangedCallback ( Widget    w, 
                                   XtPointer clientData, 
                                   XtPointer callData ) 
{
    XmToggleButtonCallbackStruct *cbs =
                             ( XmToggleButtonCallbackStruct * )  callData;

   
    /*   printf ( "%s: %s\n", XtName ( w ) , cbs->set ? "set" : "unset" ); */
}



/*************** toggle button click *********************************/
void GUI_Sankoff_toggle_radio_true(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;

  if(tmp!=NULL){
      if(XmToggleButtonGetState(parent)){
	  /* pop up the following widget */
	  XtManageChild(XtParent(tmp));
           XtUnmanageChild(XtParent(XtParent(tmp)));
           XtManageChild(XtParent(XtParent(tmp)));
      }
      else{
	  /* pop down the following widget */
	  XtUnmanageChild(XtParent(tmp));
          XtUnmanageChild(XtParent(XtParent(tmp)));
          XtManageChild(XtParent(XtParent(tmp)));
      }
  }
}



/*************** nested radio button *********************************/
void nested_radio(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;
if(XmToggleButtonGetState(parent))
  XtCallCallbacks(iws_radio_back,XmNvalueChangedCallback,call_data);
}

void GUI_Sankoff_toggle_radio_false(Widget parent, void* data, void* call_data)
{
  Widget tmp=(Widget) data;

  if(tmp!=NULL){
      if(!XmToggleButtonGetState(parent)){
	  /* pop up the following widget */
	  XtManageChild(XtParent(tmp));
          XtUnmanageChild(XtParent(XtParent(tmp)));
           XtManageChild(XtParent(XtParent(tmp)));
      }
      else{
	  /* pop down the following widget */
	  XtUnmanageChild(XtParent(tmp));
          XtUnmanageChild(XtParent(XtParent(tmp)));
          XtManageChild(XtParent(XtParent(tmp)));
      }
  }
}



                                
void SelectFileCallback ( Widget    w,
                          Widget cbWidget,
                          XtPointer callData)
{
    static Widget dialog = NULL;
    Arg           args[3];          /* arg list */
    int           n;                /* arg count */
    XmString     xmstr;        /* String label */
    char         *str,*title;

    if ( dialog )
    {
       XtRemoveAllCallbacks(dialog,XmNokCallback); 
       XtRemoveAllCallbacks(dialog,XmNokCallback);
    }
    else  
      {
         dialog = XmCreateFileSelectionDialog ( w, "openFileDialog", 
                                                NULL,(Cardinal) 0 );
	 /* XtFree(str); */
      }
   
    XtAddCallback ( dialog, XmNokCallback,
                    OKCallback,cbWidget );
    XtAddCallback ( dialog, XmNcancelCallback,
                     CancelCallback, NULL );

    XtManageChild ( dialog );
}
    
void OKCallback ( Widget w, Widget cbWidget, XtPointer callData )
{
    XmFileSelectionBoxCallbackStruct *cbs = 
                         ( XmFileSelectionBoxCallbackStruct * ) callData;
    char *fileName;
    Widget tmp; 

   /*
    * Remove the widget from the screen.
    */

    XtUnmanageChild ( w );

   /*
    * Retrieve the character string from the compound string format.
    */

    XmStringGetLtoR ( cbs->value, XmFONTLIST_DEFAULT_TAG, &fileName );

   /*
    * For this demo, just echo the selected string to standard out.
    */

    /*   printf ( "Selected file = %s\n", fileName ); */
   
   
    /* iws_queryfile=wid_textboxb(0,iws_queryfile,"Query file:",fileName,40);*/
  
     tmp = wid_textboxb(0,cbWidget,"",fileName,40);
}      
                
void CancelCallback ( Widget    w,
                      XtPointer clientData,
                      XtPointer callData )
{
    XtUnmanageChild ( w );
}
        
 void build_command(char *command_string,char *option_string)
{
  char *com;
  /*  printf("string before cat %s\n",command_string);
    printf("option_string = %s\n",option_string);*/
  
  if ((strlen(command_string) + strlen(option_string)) < MAX_STRING_LENGTH) {
    strncat(command_string,option_string,MAX_STRING_LENGTH -1);
    /* printf("string after cat %s\n",command_string); */
   GUI_strip_space(option_string);
   if (!(sargv[sargc]=(char *)malloc((strlen(option_string)+1)*sizeof(char)))){
       printf("malloc error in build_command\n");
       exit(1);
   }
   strncpy(sargv[sargc],option_string,strlen(option_string)+1);
   sargc++;
  }
 else
   {
        printf(stderr, "MAXIMUM string length exceeded in build_command, aborting!\n"); 
    exit(1);
   } 
}



            
void DisplayBusyCursor ( Widget w ) 
{
    static cursor = NULL;

    if ( !cursor ) 
        cursor = XCreateFontCursor ( XtDisplay ( w ), XC_watch );  

    XDefineCursor ( XtDisplay ( w ), XtWindow ( w ), cursor );

    XFlush ( XtDisplay ( w ) );

    XtAppAddWorkProc ( XtWidgetToApplicationContext ( w ), 
                       RemoveBusyCursor, ( XtPointer ) w );
}
    
void RemoveBusyCursor ( Widget w) 
{

    XUndefineCursor ( XtDisplay ( w ), XtWindow ( w ) );

}
        
      
