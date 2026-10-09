"""Molecule sets for the four enhanced-stereo specs.

Each entry is a dict:
  title   - short caption shown under the picture
  smiles  - (CX)SMILES; RDKit lays it out and writes a V3000 molblock
  wedges  - [(begin, end, cfg)] extra wedge bonds by 0-based atom index
            (cfg 1 = up, 3 = down); the narrow end is 'begin'
  either  - [k] k-th non-aromatic double bond drawn crossed (CFG=2)
  coll    - [(name, [items])] extra V3000 collections; items are 0-based
            atom indices for STE*, or for STEB* either ('db', k) = k-th
            non-aromatic double bond (C=O included) or ('bond', i) = 1-based bond index
  rotate  - rotate the drawing by this many radians
  sdf     - (file, record) instead of smiles: a molblock from a file
  tests   - (file, anchor) instead of smiles: the first C string literal
            after 'anchor' in a unit test
  edit    - [(old, new)] text replacements applied to an sdf/tests molblock
  relayout - draw with a fresh 2D layout (source drawing overlaps)
"""

UP, DOWN = 1, 3

ENH_DATA = "../InChI-specs/enhanced-stereo/datasets"
FIXTURES = "../InChI/INCHI-1-TEST/tests/test_unit/fixtures"
ATROP_SDF = FIXTURES + "/atropisomers_test_file_1_v2.sdf"
ATROP_TESTS = "../InChI/INCHI-1-TEST/tests/test_unit/test_atropisomers.cpp"
ENH_TESTS = "../InChI/INCHI-1-TEST/tests/test_unit/test_enhancedStereo.cpp"

# Threonine-like skeleton used for the group combinations: centres 1 and 3
THR = "C[C@@H](O)[C@H](N)C(=O)O"
BUTANOL = "C[C@H](O)CC"

ENHANCED = [
    dict(title="(S)-butan-2-ol, ABS group", smiles=BUTANOL + " |a:1|"),
    dict(title="butan-2-ol, OR group: one enantiomer, unknown which",
         smiles=BUTANOL + " |o1:1|"),
    dict(title="butan-2-ol, AND group: racemate", smiles=BUTANOL + " |&1:1|"),
    dict(title="ungrouped wedge = absolute", smiles=BUTANOL),
    dict(title="two centres, one OR group (relative configuration)",
         smiles=THR + " |o1:1,3|"),
    dict(title="two centres, one AND group (racemic diastereomer)",
         smiles=THR + " |&1:1,3|"),
    dict(title="two independent OR groups", smiles=THR + " |o1:1,o2:3|"),
    dict(title="two independent AND groups", smiles=THR + " |&1:1,&2:3|"),
    dict(title="ABS + OR", smiles=THR + " |a:1,o1:3|"),
    dict(title="ABS + AND", smiles=THR + " |a:1,&1:3|"),
    dict(title="OR + AND", smiles=THR + " |o1:1,&1:3|"),
    dict(title="three classes: ABS, OR, AND",
         smiles="C[C@@H](O)[C@H](Cl)[C@@H](Br)C |a:1,o1:3,&1:5|"),
    dict(title="edge: meso-tartaric acid in an AND group",
         smiles="O[C@@H](C(=O)O)[C@H](O)C(=O)O |&1:1,5|"),
    dict(title="trans-1,2-dimethylcyclohexane, OR group",
         smiles="C[C@@H]1CCCC[C@H]1C |o1:1,6|"),
    dict(title="edge: ungrouped wedge beside an OR group",
         smiles="C[C@@H](O)CC[C@H](C)Cl |o1:1|"),
    dict(title="two components, one OR group each",
         smiles="C[C@H](O)CC.C[C@H](O)CCC |o1:1,o2:6|"),
    dict(title="two components, both absolute: uniform /s",
         smiles="C[C@H](O)CC.C[C@@H](O)CCC |a:1,6|"),
    dict(title="edge: AND group spanning two components (dropped)",
         smiles="C[C@H](O)CC.C[C@H](O)CCC |&1:1,6|"),
    dict(title="edge: centre chiral only by isotope (ethanol-1-d), AND",
         smiles="C[C@H]([2H])O |&1:1|"),
    dict(title="edge: ABS label on a CH2 (not a stereocentre)",
         smiles=BUTANOL + " |&1:1|", coll=[("STEABS", [3])]),
    dict(title="ibuprofen, racemate",
         smiles="CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O |&1:10|"),
    dict(title="thalidomide, racemate",
         smiles="O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1 |&1:4|"),
    dict(title="ephedrine, relative configuration (OR)",
         smiles="CN[C@@H](C)[C@H](O)c1ccccc1 |o1:2,4|"),
    dict(title="BIOVIA white paper, Fig. 22", sdf=(ENH_DATA + "/biovia_white_paper_fig22.mol", 0)),
    dict(title="cyclosporin, AND group (OpenChemLib)",
         sdf=(ENH_DATA + "/openchemlib_cyclosporin_rac.mol", 0)),
    dict(title="two AND groups (RDKit)", sdf=(ENH_DATA + "/rdkit_two_and_groups.mol", 0)),
    dict(title="OR and AND groups (Indigo)", sdf=(ENH_DATA + "/indigo_rel_and_rac.mol", 0)),
    dict(title="OR enantiomer (Ketcher)", sdf=(ENH_DATA + "/ketcher_or_enantiomer.mol", 0)),
    dict(title="ChEMBL record, OR only (CDPKit)",
         sdf=(ENH_DATA + "/cdpkit_chembl_sterel_only.mol", 0)),
    dict(title="EOS49757, ABS only", sdf=(ENH_DATA + "/EOS49757_steabs_only.mol", 0)),
]

# 1-bromo-1-fluoro-3-chlorobuta-1,2-diene (deck slide 18):
# 0 F, 1 C, 2 Br, 3 C (central), 4 C, 5 CH3, 6 Cl
DECK = "FC(Br)=C=C(C)Cl"
DECK_RA = [(1, 2, UP), (4, 6, DOWN)]

ALLENES = [
    dict(title="1,3-dichloroallene, two wedges", smiles="ClC=C=CCl",
         wedges=[(1, 0, UP), (3, 4, DOWN)]),
    dict(title="1,3-dichloroallene, mirror image", smiles="ClC=C=CCl",
         wedges=[(1, 0, DOWN), (3, 4, UP)]),
    dict(title="1,3-dichloroallene, single wedge", smiles="ClC=C=CCl",
         wedges=[(1, 0, UP)]),
    dict(title="edge: flat drawing, axis undefined", smiles="ClC=C=CCl"),
    dict(title="deck slide 18, two wedges", smiles=DECK, wedges=DECK_RA),
    dict(title="deck slide 18, mirror image", smiles=DECK,
         wedges=[(1, 2, DOWN), (4, 6, UP)]),
    dict(title="deck slide 18, single wedge", smiles=DECK, wedges=[(4, 6, UP)]),
    dict(title="edge: two wedges at one end contradict", smiles=DECK,
         wedges=[(1, 0, UP), (1, 2, UP)]),
    dict(title="deck slide 18, rotated drawing", smiles=DECK, wedges=DECK_RA, rotate=1.0),
    dict(title="ABS group on the central atom", smiles=DECK, wedges=DECK_RA,
         coll=[("STEABS", [3])]),
    dict(title="OR group on the central atom", smiles=DECK, wedges=DECK_RA,
         coll=[("STEREL1", [3])]),
    dict(title="AND group on the central atom", smiles=DECK, wedges=DECK_RA,
         coll=[("STERAC1", [3])]),
    dict(title="AND group, single wedge", smiles=DECK, wedges=[(4, 6, UP)],
         coll=[("STERAC1", [3])]),
    dict(title="edge: AND group names a terminal atom: ignored, but /s is lost (bug)", smiles=DECK,
         wedges=DECK_RA, coll=[("STERAC1", [1])]),
    dict(title="edge: symmetric end, no axis", smiles="CC(C)=C=CCl",
         wedges=[(1, 0, UP), (4, 5, DOWN)]),
    dict(title="penta-2,3-diene", smiles="CC=C=CC", wedges=[(1, 0, UP), (3, 4, DOWN)]),
    dict(title="glutinic acid (penta-2,3-dienedioic acid)", smiles="OC(=O)C=C=CC(=O)O",
         wedges=[(3, 1, UP), (5, 6, DOWN)]),
    dict(title="1,3-di-tert-butylallene", smiles="CC(C)(C)C=C=CC(C)(C)C",
         wedges=[(4, 1, UP), (6, 7, DOWN)]),
    dict(title="allene AND + stereocentre ABS", smiles="C[C@H](O)C=C=CCl |a:1|",
         wedges=[(5, 6, UP)], coll=[("STERAC1", [4])]),
    dict(title="two allene axes, two OR groups", smiles="ClC=C=CCC=C=CCl",
         wedges=[(1, 0, UP), (3, 4, DOWN), (5, 4, UP), (7, 8, DOWN)],
         coll=[("STEREL1", [2]), ("STEREL2", [6])]),
    dict(title="allene next to an E double bond", smiles="Cl/C=C/C=C=CC",
         wedges=[(3, 2, UP), (5, 6, DOWN)]),
    dict(title="cyclonona-1,2-diene", smiles="C1CCCC=C=CCC1",
         wedges=[(4, 3, UP), (6, 7, DOWN)]),
    dict(title="edge: butatriene (even cumulene, E/Z not axial)", smiles="ClC=C=C=CCl"),
    dict(title="edge: pentatetraene (odd cumulene, length 4)", smiles="ClC=C=C=C=CCl",
         wedges=[(1, 0, UP), (5, 6, DOWN)]),
]

# BrClC=CClF (deck slides 25-28) in both geometries; the C=C is double bond 0
BROMO_Z = "Cl/C(Br)=C(/Cl)F"
BROMO_E = "Cl/C(Br)=C(\\Cl)F"
# Two independent C=C: F(Br)C=C(Cl)-CH2-C(Cl)=C(I)F
DIENE = "F/C(Br)=C(/Cl)C/C(Cl)=C(\\I)F"
DIENE_FLIP = "F/C(Br)=C(\\Cl)C/C(Cl)=C(\\I)F"
# CH3-CH(OH)-C(Br)=C(Cl)F (deck slide 30): centre 1, double bond 0
BOND_CENTRE = "C[C@@H](O)/C(Br)=C(/Cl)F"

EZ = [
    dict(title="BrClC=CClF, no collection", smiles=BROMO_Z),
    dict(title="BrClC=CClF, other isomer", smiles=BROMO_E),
    dict(title="STEBABS: absolute", smiles=BROMO_Z, coll=[("STEBABS", [("db", 0)])]),
    dict(title="STEBREL1: one isomer, unknown which", smiles=BROMO_Z,
         coll=[("STEBREL1", [("db", 0)])]),
    dict(title="STEBRAC1: E/Z mixture", smiles=BROMO_Z, coll=[("STEBRAC1", [("db", 0)])]),
    dict(title="OR group, drawn as the other isomer (normalised)", smiles=BROMO_E,
         coll=[("STEBREL1", [("db", 0)])]),
    dict(title="edge: crossed double bond (CFG=2, RDKit draws it plain) is unknown, not OR", smiles=BROMO_Z,
         either=[0], coll=[("STEBREL1", [("db", 0)])]),
    dict(title="diene, no collection", smiles=DIENE),
    dict(title="diene, both bonds in one AND group", smiles=DIENE,
         coll=[("STEBRAC1", [("db", 0), ("db", 1)])]),
    dict(title="diene, one bond flipped, one AND group", smiles=DIENE_FLIP,
         coll=[("STEBRAC1", [("db", 0), ("db", 1)])]),
    dict(title="diene, two AND groups", smiles=DIENE,
         coll=[("STEBRAC2", [("db", 0)]), ("STEBRAC1", [("db", 1)])]),
    dict(title="diene, one OR bond, the other absolute", smiles=DIENE,
         coll=[("STEBREL1", [("db", 1)])]),
    dict(title="diene, ABS + AND", smiles=DIENE,
         coll=[("STEBABS", [("db", 0)]), ("STEBRAC1", [("db", 1)])]),
    dict(title="diene, AND + OR", smiles=DIENE,
         coll=[("STEBRAC1", [("db", 0)]), ("STEBREL1", [("db", 1)])]),
    dict(title="OR bond + absolute stereocentre", smiles=BOND_CENTRE + " |a:1|",
         coll=[("STEBREL2", [("db", 0)])]),
    dict(title="absolute bond + racemic stereocentre", smiles=BOND_CENTRE + " |&1:1|",
         coll=[("STEBABS", [("db", 0)])]),
    dict(title="edge: single bond in a bond collection (ignored)", smiles=BROMO_Z,
         coll=[("STEBREL1", [("bond", 2)])]),
    dict(title="edge: bond in two collections (dropped)", smiles=DIENE,
         coll=[("STEBREL1", [("db", 0)]), ("STEBRAC1", [("db", 0), ("db", 1)])]),
    dict(title="edge: bond index out of range (dropped)", smiles=DIENE,
         coll=[("STEBRAC1", [("db", 0), ("bond", 99)])]),
    dict(title="two equal components, OR on one, AND on the other",
         smiles=BROMO_Z + "." + BROMO_Z,
         coll=[("STEBREL1", [("db", 0)]), ("STEBRAC1", [("db", 1)])]),
    dict(title="stilbene, OR group", smiles="c1ccc(cc1)/C=C/c1ccccc1",
         coll=[("STEBREL1", [("db", 0)])]),
    dict(title="fumaric acid, AND group", smiles="OC(=O)/C=C/C(=O)O",
         coll=[("STEBRAC1", [("db", 1)])]),
    dict(title="citral (geranial), AND group", smiles="CC(C)=CCC/C(C)=C/C=O",
         coll=[("STEBRAC1", [("db", 1)])]),
    dict(title="acetophenone oxime (C=N), OR group", smiles="C/C(=N/O)c1ccccc1",
         coll=[("STEBREL1", [("db", 0)])]),
    dict(title="all-trans retinal, AND + OR on two of four bonds",
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=O",
         coll=[("STEBRAC1", [("db", 2)]), ("STEBREL1", [("db", 4)])]),
]

ATROPISOMERS = [
    dict(title="(M)-6,6'-dinitrodiphenic acid", sdf=(ATROP_SDF, 0), relayout=True),
    dict(title="(P)-6,6'-dinitrodiphenic acid", sdf=(ATROP_SDF, 1), relayout=True),
    dict(title="(M)-2,2'-dibromo-6,6'-dichlorobiphenyl", sdf=(ATROP_SDF, 6)),
    dict(title="same, other drawing", sdf=(ATROP_SDF, 7)),
    dict(title="(P)-2,2'-dibromo-6,6'-dichlorobiphenyl", sdf=(ATROP_SDF, 35)),
    dict(title="3D input (projection), 44° torsion", tests=(ATROP_TESTS, "k_biaryl3d_p44 =")),
    dict(title="3D input (projection), -44° torsion", tests=(ATROP_TESTS, "k_biaryl3d_m44 =")),
    dict(title="2D, two wedges", tests=(ATROP_TESTS, "k_dummy1_molblock =")),
    dict(title="2D, single wedge", tests=(ATROP_TESTS, "k_dummy1_molblock ="),
         edit=[("  8  7  1  1", "  8  7  1  0")]),
    dict(title="edge: 2D flat, axis undefined", tests=(ATROP_TESTS, "k_dummy1_molblock ="),
         edit=[("  2  3  1  1", "  2  3  1  0"), ("  8  7  1  1", "  8  7  1  0")]),
    dict(title="2D V3000, ABS collection on the axis",
         tests=(ENH_TESTS, "test_EnhancedStereochemistry_2_atropisomer")),
    dict(title="2D V3000, AND collection on the axis",
         tests=(ENH_TESTS, "test_EnhancedStereochemistry_2_atropisomer"),
         edit=[("MDLV30/STEABS", "MDLV30/STERAC1")]),
    dict(title="edge: both wedges on one side (degenerate)", sdf=(ATROP_SDF, 14)),
    dict(title="(M)-BINOL", sdf=(ATROP_SDF, 16)),
    dict(title="(P)-BINOL", sdf=(ATROP_SDF, 18)),
    dict(title="lignan, biaryl in 8-ring (M); labelled atropisomer, no axis since 2ea128d", sdf=(ATROP_SDF, 21)),
    dict(title="same lignan (P)", sdf=(ATROP_SDF, 25)),
    dict(title="edge: bare bridged biaryl, 7-ring (fast ring flip)",
         tests=(ATROP_TESTS, "k_dummy12_molblock =")),
    dict(title="edge: acyclic C14H26 axis", sdf=(ATROP_SDF, 42)),
    dict(title="aryl N,N-dimethylamide, aryl-carbonyl axis", sdf=(ATROP_SDF, 43)),
    dict(title="edge: C-N axis (out of scope, v1 carbon only)", sdf=(ATROP_SDF, 26)),
    dict(title="metolachlor-like herbicide, C-N axis", sdf=(ATROP_SDF, 44)),
    dict(title="edge: N-N axis (M)", sdf=(ATROP_SDF, 28)),
    dict(title="edge: N-N axis (P)", sdf=(ATROP_SDF, 29)),
    dict(title="edge: dihydrobiindole, no axis", sdf=(ATROP_SDF, 30)),
    dict(title="edge: 3,3'-bipyridine, too little hindrance", sdf=(ATROP_SDF, 33)),
    dict(title="edge: overcrowded alkene (Z,M)", sdf=(ATROP_SDF, 45)),
    dict(title="edge: overcrowded alkene (E,P)", sdf=(ATROP_SDF, 47)),
]

TABLES = [
    ("enhanced", "Enhanced stereochemistry (ABS / OR / AND)", ENHANCED),
    ("allenes", "Allenes", ALLENES),
    ("ez", "Double-bond E/Z", EZ),
    ("atropisomers", "Atropisomers", ATROPISOMERS),
]
