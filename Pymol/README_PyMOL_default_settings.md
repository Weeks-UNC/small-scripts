# PYMOL DEFAULT SETTINGS AND ALIASES GUIDE

### Custom default settings and aliases (custom functions) for PyMOL are stored in a configuration file titled pymolrc. The pymolrc file will initalize everytime you start PyMOL to incorrporate your custom setting automatically.  

### This tutorial will cover:

1. How to edit your PyMOL pymolrc file with default settings and/or custom aliases.
2. What to include in your PyMOL pymolrc file
3. How to write an basic PyMOL alias (custom function).
4. Example pymolrc files

# 1. Edit your pymolrc file

## QUICK SET-UP: Edit your pymolrc file in PyMOL

1. Open a window of PyMOL and choose, **File > Edit pymolrc**
2. Add your desired commands and aliases and save.  

    **If you plan to use the Weeks lab's default [pymolrc.pml](pymolrc.pml) file:**
    - Copy the contents of the [pymolrc.pml](pymolrc.pml) file from the Weeks-UNC/small-scripts/Pymol github directory.
    - Paste the contents of the [pymolrc.pml](pymolrc.pml) file into the open notepad file and save.

## ADVANCED SET-UP: Edit your pymolrc file within the command line

### WINDOWS USERS:
1. Open command prompt window and paste
    > notepad "%HOMEDRIVE%%HOMEPATH%\pymolrc.pml"
2. This will open your pymolrc.pml file using notepad allowing you to make edits.
3. Add your desired commands and aliases and save.  

    **If you plan to use the Weeks lab's default pymolrc.pml file:**
    - Copy the contents of the pymolrc.pml file from the Weeks-UNC/small-scripts/Pymol github directory.  
    - Paste the contents of the pymolrc.pml file into your open notepad file and save.

### MAC/LINUX USERS:
1. Open your terminal and type  
    > nano ~/.pymolrc

    or (depending on your preference)

    > vim ~/.pymolrc
2. This will open your pymolrc.pml file in your perfered text editor allowing you to make edits.
3. Add your desired commands and aliases and save. 

    **If you plan to use the Weeks lab's default pymolrc.pml file:**  
    - Copy the contents of the pymolrc.pml file from the Weeks-UNC/small-scripts/Pymol github directory.  
    - Paste the contents of the pymolrc.pml file into your open nano session and save.

### MULTIPLE PYMOLRC FILES:

PyMOL will even load multiple pymolrc files, however only from the same directory.  
You can have multiple scripts: e.g. pymolrc-settings.pml and pymolrc-misc.pml in your home directory.

- pymolrc-settings.pml can e.g. be used to define 'permanent' custom Settings that you rarely change
- pymolrc-misc.pml can e.g. be used to define more transient custom Settings, such as Working Directory or Fetch Path

You can query which pymolrc files have been loaded:

> PyMOL> print invocation.options.deferred

**NOTE:** You do NOT need to use the PYMOL API format in pymolrc.pml files.

# 2. What to include in your pymolrc file

### Your pymolrc file should contain settings that you always apply to your PyMOL sessions. 
### Examples:

- ### Apply basic settings to PyMOL 

      # Change backround color to white for easier viewing.
      bg_color white

      # Display sequence viewer.
      set seq_view, 1

      # Toggles whether PyMOL loops movies.
      set movie_loop, 1

      # Save fetched PDB files here.
      set fetch_path, /your/fetch/path

- ### Apply consistent styles
        
      # Adjust the thickness of atomic bonds viewed as sticks.
      set stick_radius, 0.3

      # Adjust the radius of spheres.
      set sphere_scale, 1

      # Settings related to surface features.
      set surface_quality, 1
      set solvent_radius, 1.6
      set transparency, 0.5
      set surface_color, grey80

      # Settings related to labels.
      set label_size, 60
      set label_outline_color, 1
      set label_color, 0
      set label_position, [0, 0, 10]
        
- ### Apply raytrace settings

      # Enable multi-thread processing for faster raytrace performance.
      set max_threads, 16
    
      # Increase memory allowance for faster raytracing performance.
      set hash_max, 2048

      # Sets specular reflection intensity to lowwer value for aesthetics. 
      set specular, 0.1

      # Removes distracing shadows from raytraced images.
      set ray_shadows, off

      # Improves antiliasing/edge smoothing for raytraced images.
      set antialias, 2

- ### Apply aliases (custom function)

    Adding an alias (custom function) to your pymolrc file will ensure that it can be used by default anytime you use PyMOL.  
    Creating an alias is covered below in section 3.


Read more at: https://pymolwiki.org/index.php/Pymolrc

# 3. Create your own PyMOL alias

**alias** allows you to bind a commonly used command (or series of commands) to a single PyMOL keyword of your choice.  
In essence, aliases allow you to add customized functions to PyMOL.

## USAGE
 
> PyMOL> alias *name*, *command-sequence*

alias: name of the PyMOL command for creating a new alias.  
*name*: name of your new custom function.  
*command-sequence*: list of your desired PyMOL> command line commands. 

## EXAMPLES

### Aliases can be used to simplify a commonly used command or series of commands into a single function name.
    
- remove_protein: alias which removes all amino acid residues.
    > PyMOL> alias remove_protein, remove (byres polymer & name CA)

- remove_rna: alias which removes all RNA residues.
    > Pymol> alias remove_rna, remove (byres polymer & name O2')


### Aliases are speficically useful for consistantly applying style settings for publication.

- rna: alias which implements a cartoon style for RNAs. (removes protein)
    > PyMOL> alias rna, set cartoon_ring_mode, 3; set cartoon_ring_finder, 1; remove resn hoh; remove inorganic and not resn STP; cartoon oval; set cartoon_oval_length, 0.75; set cartoon_oval_width, 0.25; color greyish, (polymer); util.cbao (organic); remove (byres polymer & name CA)

- rna_surface: alias which displays the surface of RNAs.
    > PyMOL> alias rna_surface, show surface, byres polymer & name O2'


|Style: PyMOL Default|Style: rna Alias|Style: rna + rna_surface Aliases|
|:-:|:-:|:-:|
|![Style: PyMOL Default](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/3E5C_RNA_Default_Style.png?raw=true)|![Style: rna Alias](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/3E5C_RNA_rna_Style.png?raw=true)|![Style: rna + rna_surface Alias](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/3E5C_RNA_rna_s_Style.png?raw=true)

- rna_protein: alias which implements a cartoon style for RNA and protein.
    > Pymol> alias rna_protein, set cartoon_ring_mode, 3; set cartoon_ring_finder, 1; remove resn hoh; remove inorganic and not resn STP; cartoon oval; set cartoon_oval_length, 0.75; set cartoon_oval_width, 0.25; color greyish, (polymer); util.cbao (organic); color lightpink, (byres polymer & name CA); cartoon automatic, (byres polymer & name CA)

- rna_ribbon: alias which changes the RNA cartoon style to display a ribbon only.
    > Pymol> alias rna_ribbon, set cartoon_ring_finder, 0


|Style: PyMOL Default|Style: rna_protein Alias|Style: rna_protein + rna_ribbon Aliases|
|:-:|:-:|:-:|
|![Style: PyMOL Default](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/5GX2_RNA-Protein_default_Style.png?raw=true)|![Style: rna_protein Alias](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/5GX2_RNA_rna_protein_Style.png?raw=true)|![Style: rrna_protein + rna_ribbon Aliases](https://github.com/Weeks-UNC/small-scripts/blob/master/Pymol/Images/5GX2_RNA_rna+rna_ribbon_Style.png?raw=true)


## NOTES

Multiple commands should be seperated by a semi-colon.

After defining an alias, you can implement your alias simply by entering the alias name into the PyMOL> command line.
	
Read more at: https://pymolwiki.org/index.php/Alias

# 4. Example pymolrc files

### From https://pymolwiki.org/index.php/Pymolrc

    # simple test: change background color of PyMOL window
    bg blue

    # this will run the script in the specified location
    run /path/to/home/pymol/load_sep.py

    # your favorite settings
    set movie_loop, 0
    set two_sided_lighting, 1

    set label_size, 60
    set label_outline_color, 1
    set label_color, 0
    set label_position, [0, 0, 10]

    # for images:
    #   antialias =1 smooths jagged edges, 0 turns it off
    set antialias, 1

    #   stick_radius -adjust thickness of atomic bonds
    set stick_radius, 0.3

    # save fetched PDB files here
    set fetch_path, /your/fetch/path

    # Personal short-cut to color_obj function
    import color_obj
    cmd.extend("co",color_obj.color_obj)

### From https://betainverse.wordpress.com/2017/04/14/pymol-default-settings-pymolrc/

    # default path that PyMOL uses to load files from before it tries to download them from the PDB.
    set fetch_path, /home/edmonds/PDBs/
    
    # use white background
    bg_color white
    
    # set ray tracing settings
    set ray_opaque_background, off
    set ray_shadows,off
    
    # run custom pymol scripts
    run /home/edmonds/scripts/pymolscripts/data2bfactor.py
    run /home/edmonds/scripts/pymolscripts/spectrumany.py
    
### From http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/

    # use white background
    set bg_rgb, 1 1 1

    # set orthoscopic, rather than perspective as the default
    set orthoscopic, 1

    # set size of graphics window
    #viewport 800, 600
    viewport 1100, 900

    set hash_max, 170

    # turn off auto zoom
    set auto_zoom, 0

    set cgo_line_radius, 0.05
    set line_width, 2
    set cgo_line_width, 2
    set ribbon_radius, .2
    set ribbon_width, 2
    set ribbon_sampling, 1
    set cartoon_smooth_loop, 0

    set stick_radius, .2
    set fog, 1.

    # set spec_refl to 2. for bright reflections, .5 for dimmer
    set spec_refl, 1.5
    # set spec_power to 200 for tight reflections, 40 for broader
    set spec_power, 100

    # default to antialiased rendering (on or 1)
    set antialias, on
    # to set background to transparent for ray tracing, set the following off, or 0
    set ray_opaque_background, on

    # The following may not be necessary now as one can set the preferred
    # color space see the "Display -> Color Space" menu item in the main GUI
    #optimized rgb values for cmyk output:
    set_color _dblue= [0.05 , 0.19 , 0.57]
    set_color _blue=  [0.02 , 0.50 , 0.72]
    set_color gold, [1.0, 0.8, 0.0]
    set_color gold2, [0.90, 0.6, 0.0]
    set_color lightblue, [0.6, 0.8, 1.0]
    set_color lightgreen, [0.,1.0, 0.3]
    set_color darkgreen, [0.,0.7,0.0]


    # run some of my scripts (see http://adelie.biochem.queensu.ca/~rlc/work/pymol/ for these)
    run c:\Program Files\DeLano Scientific\PyMOL\load_best.py
    run c:\Program Files\DeLano Scientific\PyMOL\load_models.py
    run c:\Program Files\DeLano Scientific\PyMOL\load_sep.py
    run c:\Program Files\DeLano Scientific\PyMOL\load_list.py

### From Seth Veenbaas.
### Settings designed for RNA.

    # Created by: Seth Veenbaas
    # Weeks Lab UNC Chapel Hill
    # Last date modified: 10/06/22

    # Installation instructions at https://github.com/Weeks-UNC/small-scripts/tree/master/Pymol


    #### Startup settings ####
    # set size of PYMOL window on boot
    viewport 1500, 1000

    # enable multi-thread processing
    set max_threads, 16

    # increase raytracing memory allowance
    set hash_max, 2048

    # change backround color
    bg_color white

    # display sequence viewer
    set seq_view, 1

    # sets specular reflection intensity
    set specular, 0.1

    # controls appearence of shadows for ray-traced images
    set ray_shadows, off

    # controls antiliasing/edge smoothing for ray-traced images
    set antialias, 2

    # settings related to surface features
    set surface_quality, 1
    set solvent_radius, 1.6
    set transparency, 0.5
    set surface_color, grey80

    # settings related to rendering meshes
    set mesh_quality, 2
    set mesh_type, 0
    set mesh_width, 0.5
    set mesh_radius, 0.015

    # orthoscopic turns on and off the perspective handling
    set orthoscopic, off

    # define "greyish" color
    set_color greyish, [0.625, 0.7, 0.7]
    set_color novel, [0.0, 0.4314322351168238, 0.1118361643280874]
    set_color known, [0.908605075491407, 0.3955005147576708, 0.0]

    #### Custom Aliases #####

    # RNA cartoon command
    # Sets RNA cartoon style, removes ions and water, colors RNA lightteal and ligand brightorange.
    alias rna, set cartoon_ring_mode, 3; set cartoon_ring_finder, 1; remove resn hoh; remove inorganic and not resn STP; cartoon oval; set cartoon_oval_length, 0.75; set cartoon_oval_width, 0.25; color greyish, (polymer); util.cbao (organic); remove (byres polymer & name CA)

    # RNA cartoon command + extract ligand as an object
    alias rna_ligand_obj, rna; extract Ligand, organic; util.cbao (organic)

    # RNA cartoon command + colors proeteins lightpink instead of removing.
    alias rna_protein, set cartoon_ring_mode, 3; set cartoon_ring_finder, 1; remove resn hoh; remove inorganic and not resn STP; cartoon oval; set cartoon_oval_length, 0.75; set cartoon_oval_width, 0.25; color greyish, (polymer); util.cbao (organic); color lightpink, (byres polymer & name CA); cartoon automatic, (byres polymer & name CA)

    # RNA & Protein command + extracts objects for proteins and ligands
    alias rna_protein_obj, rna_protein; extract Protein, (byres polymer & name CA); extract Ligand, organic; orient

    # Ribosome
    alias ribosome, set surface_quality, 0; rna_protein_obj

    # RNA ribbon command
    alias rna_ribbon, set cartoon_ring_finder, 0

    # show surface of RNA only
    alias rna_surface, show surface, byres polymer & name O2'

    # show surface of RNA and Protein
    alias rna_protein_surface, show surface, byres polymer & name O2'; show surface, byres polymer & name CA; set surface_color, lightpink, byres polymer & name CA

    # ligand colored by element
    alias ligand, util.cbao (organic)

    # shows details of ligand binding site.
    alias binding_site, ligand; select Binding_Site, byresidue polymer within 4.5 of organic; set cartoon_ring_finder, 0; show sticks, Binding_Site; util.cnc Binding_Site; contacts Ligand, Binding_Site; disable contacts_all & contacts_aa & contacts_dd; hide labels, contacts; color purple, contacts_polar; color purple, contacts_polar_ok; color tv_yellow, contacts_all; zoom (Ligand)

    # Select ligands
    alias select_ligand, select ligand, (organic) and not (byres polymer & name O2')

    # remove organic ligands
    alias remove_ligand, remove organic and not (byres polymer & name O2')

    # select all amino acids
    alias select_protein, select protein, (byres polymer & name CA)

    # remove all amino acids
    alias remove_protein, remove (byres polymer & name CA)

    # select all RNAs
    alias select_rna, select RNA, (byres polymer & name O2')

    # remove all rnas
    alias remove_rna, remove (byres polymer & name O2')

    # sorts pocket predictions based on if they contact RNA-only, protein-only, or RNA-Protein.
    alias sort_pockets, extract Protein_Pockets, byresidue resn STP within 4 of (byres polymer & name CA); extract RNA-Protein_Pockets, byresidue resn STP and Protein_Pockets and within 4 of (byres polymer & name O2'); extract RNA_Pockets, byresidue resn STP and not Protein_Pockets; color purple, Protein_Pockets; color skyblue, RNA-Protein_Pockets color tv_orange, RNA_Pockets

    # set a-sphere radius based on data in b factor column. (for use with *real_sphere.pdb files)
    alias fpocket_radius, cmd.alter('resn STP', 'vdw = b - 1.65')



        
