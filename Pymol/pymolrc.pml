# Created by: Seth Veenbaas
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

# Show surface of RNA only
alias rna_surface, show surface, byres polymer & name O2'

# Show surface of RNA and Protein
alias rna_protein_surf, show surface, byres polymer & name O2'; show surface, byres polymer & name CA; set surface_color, lightpink, byres polymer & name CA

# ligand colored by element
alias ligand, util.cbao (organic)

# shows details of ligand binding site.
# requires installing a pymol plug-in
### In Pymol navigate to: Plugin > Plugin Manager > Install New Plugin 
### paste into URL field: https://raw.githubusercontent.com/dkoes/show_contacts/master/show_contacts.py > click Fetch
alias binding_site, ligand; select Binding_Site, byresidue polymer within 4.5 of (organic); set cartoon_ring_finder, 0; show sticks, Binding_Site; util.cnc Binding_Site; contacts (organic), Binding_Site; disable contacts_all & contacts_aa & contacts_dd; hide labels, contacts; color deeppurple, contacts_polar; color violet, contacts_polar_ok; color tv_yellow, contacts_all

# Select ligands
alias select_ligand, select ligand, (organic) and not (byres polymer & name O2')

# Remove ligands
alias remove_ligand, remove organic and not (byres polymer & name O2')

# Select all amino acids
alias select_protein, select protein, (byres polymer & name CA)

# Remove all amino acids
alias remove_protein, remove (byres polymer & name CA)

# Selects all RNAs
alias select_rna, select RNA, (byres polymer & name O2')

# Remove all rnas
alias remove_rna, remove (byres polymer & name O2')

# Sorts pocket predictions based on if they contact RNA-only, protein-only, or RNA-Protein.
alias sort_pockets, extract Protein_Pockets, byresidue resn STP within 4 of (byres polymer & name CA); extract RNA-Protein_Pockets, byresidue resn STP and Protein_Pockets and within 4 of (byres polymer & name O2'); extract RNA_Pockets, byresidue resn STP and not Protein_Pockets; color purple, Protein_Pockets; color skyblue, RNA-Protein_Pockets color tv_orange, RNA_Pockets

# Sets a-sphere radius based on data in b factor column. (for use with *real_sphere.pdb files)
alias fpocket_radius, cmd.alter('resn STP', 'vdw = b - 1.65')
