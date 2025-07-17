from pymol import cmd

HIDE_CYCLINS = True


#----reinitialize the session
cmd.reinitialize()

#-----LOADING STRUCTURES-------

# fetch cdk7 thz-1 structure (chain J)
cmd.fetch("6XD3", "CDK7")
# fetch cdk8 structure in complex with compound 22 (chain A)
cmd.fetch("5CEI", "CDK8")
# load cdk19 Alphafold structure (no available PDB structures)
cmd.load("AF-Q9BWU1-F1-model_v4.pdb", "CDK19")
# fetch cdk9 KB-0742 structure (chain A)
cmd.fetch("8k5r", "CDK9")
# fetch cdk13 thz-1 structure (chain A)
cmd.fetch("7NXJ", "CDK13")
# fetch cdk12 hq461 structure (chain B)
cmd.fetch("8BUG", "CDK12")

#------COLOR-CODING STRUCTURES--------
def color_by_hex(selection, hex_code, color_name=None):
    # Remove '#' if present
    hex_code = hex_code.lstrip("#")
    # Convert to RGB tuple (0.0 – 1.0 scale)
    rgb = tuple(int(hex_code[i:i+2], 16)/255.0 for i in (0, 2, 4))
    # Optional name for the new color
    if not color_name:
        color_name = f"custom_{hex_code}"
    # Register the color in PyMOL
    cmd.set_color(color_name, rgb)
    # Apply the color to the selection
    cmd.color(color_name, selection)


color_dict = {"CDK7": "#00ffff",
"CDK8": "#1f77b4", 
"CDK19": "#6495ed",
"CDK9": "#9400d3",
"CDK12": "#2ca02c",
"CDK13": "#556b2f"}


#-----CATALYTIC SITE SELECTION AND ALIGNMENT-----
#select the ligands
cmd.select("CDK7_ligand", "CDK7 and resn V0G and resi 401") #THZ1
cmd.select("CDK8_ligand", "CDK8 and resn 50R and resi 508") #Compound 22
cmd.select("CDK9_ligand", "CDK9 and resn VQE and resi 401") #KB-0742
cmd.select("CDK12_ligand", "CDK12 and resn RPW and resi 1101") #HQ461
cmd.select("CDK13_ligand", "CDK13 and chain A and resn 5I1 and resi 2000") #THZ1


#selecting the CDK
cmd.select("CDK7_full", "CDK7 and chain J and not CDK7_ligand")
cmd.select("CDK8_full", "CDK8 and chain A and not CDK8_ligand")
cmd.select("CDK19_full", "CDK19")
cmd.select("CDK9_full", "CDK9 and chain A and not CDK9_ligand")
cmd.select("CDK12_full", "CDK12 and chain B and not CDK12_ligand")
cmd.select("CDK13_full", "CDK13 and chain A and not CDK13_ligand")

#color them
for i in color_dict.keys():
    color_by_hex(f'{i}_full', color_dict[i])

#selecting residues ranging from beginning to end of KLIFS catalytic site
cmd.select("cdk7_catalytic", "CDK7 and chain J and resi 16-159")
cmd.select("cdk8_catalytic", "CDK8 and chain A and resi 25-177")
cmd.select("cdk19_catalytic", "CDK19 and resi 25-177")
cmd.select("cdk9_catalytic", "CDK9 and chain A and resi 23-171")
cmd.select("cdk12_catalytic", "CDK12 and chain B and resi 731-881")
cmd.select("cdk13_catalytic", "CDK13 and chain A and resi 709-859")


#select the associated cyclins...
#NOTE: not present for current CDK19 AlphaFold structure (could include with multimer)
cmd.select("CDK7_cyclin", "CDK7 and chain I")
cmd.select("CDK8_cyclin", "CDK8 and chain B")
cmd.select("CDK9_cyclin", "CDK9 and chain B")
cmd.select("CDK12_cyclin", "CDK12 and chain C")
cmd.select("CDK13_cyclin", "CDK13 and chain B")


#-----HIDE THINGS THAT AREN'T IN THE CATALYTIC DOMAIN------
set1 = ["cdk9_catalytic", "cdk7_catalytic", "cdk8_catalytic", "cdk13_catalytic", "cdk12_catalytic", "cdk19_catalytic"]

for sel in set1:
    n = sel.split('_')[0].upper()

    if n=='CDK19':
        cmd.select(f"everything_{n}", f"{n} and not {n}_full")
    else:
        #hide other stuff and keep ligands + cyclin
        cmd.select(f"everything_{n}", f"{n} and not {n}_full and not {n}_ligand and not {n}_cyclin")

        #hide other organics
        cmd.hide("everything", f"{n} and organic and not {n}_ligand")

        if HIDE_CYCLINS==True:
            cmd.hide("everything", f"{n}_cyclin")


    cmd.hide("everything", f"everything_{n}")
    

#also hide waters
cmd.hide("everything", "resn HOH")


#-------SET TRANSPARENCY-----
#setting things outside of catalytic core to higher transparency

for sel in set1:
    n = sel.split('_')[0].upper()

    # Create the inverse selection (everything else)
    cmd.select(f"non_catalytic_{n}", f"{n}_full and not {sel}")
    # Set transparency on the non-catalytic selection
    cmd.set("cartoon_transparency", 0.8, f"non_catalytic_{n}")  # 0 = opaque, 1 = fully transparent


#----ALIGNMENT-------
#align these residues
reference = "cdk9_catalytic"
targets = ["cdk7_catalytic", "cdk8_catalytic", "cdk13_catalytic", "cdk12_catalytic", "cdk19_catalytic"]

#and write the alignment stats to txt file
with open("alignment_stats.txt", "w") as f:
    f.write("Target\tReference\tRMSD\tAlignedAtoms\tAlignmentCycles\tRMSD_pre_refinement\tNo_atoms_pre_refinement\tRaw_score\tNo_residues_aligned\n")
    for target in targets:
        rmsd, n_aligned, n_cycles, rmsd_pre_align, n_aligned_pre, raw_score, n_residues = cmd.align(target, reference)
        f.write(f"{target}\t{reference}\t{rmsd}\t{n_aligned}\t{n_cycles}\t{rmsd_pre_align}\t{n_aligned_pre}\t{raw_score}\t{n_residues}\n")

#----and set view and other plotting things----
cmd.set_view((
    -0.326536447,  0.859916389,  0.392320514,
     0.126545385,  0.451109022, -0.883449554,
    -0.936673164, -0.238833144, -0.256121993,
     0.000075411, -0.000149459, -205.771041870,
    25.209163666, 41.874084473,  8.689048767,
 -75678.609375000, 76090.328125000, -20.000000000
))

cmd.set("ambient", 0.5)
cmd.set("ray_trace_mode", 1)