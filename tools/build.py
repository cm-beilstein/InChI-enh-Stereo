"""Build the static page: data/<table>.sdf, svg/<table>-<id>.svg, index.html.

Usage (from the repo root):
  INCHI_OLD=/path/to/dev/inchi-1 INCHI_NEW=/path/to/atropisomers/inchi-1 \\
      python tools/build.py

  INCHI_OLD  inchi-1 built from IUPAC-InChI/InChI 'dev' (standard options)
  INCHI_NEW  inchi-1 built from cm-beilstein's 'atropisomers' branch,
             run with -EnhancedStereochemistry
"""

import hashlib
import html
import math
import os
import re
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from molecules import TABLES, UP, DOWN

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.dirname(ROOT)
OLD_OPTS = []
NEW_OPTS = ["-EnhancedStereochemistry"]
SVG_W, SVG_H = 320, 240
CANON_FONT_PX = 9
CANON_GAP_PX = 8          # number to bare carbon
CANON_GAP_SYMBOL_PX = 14  # number to an atom symbol (O, Cl, ...)
ATOM_FONT_MAX_PX = 18     # atom symbols on small molecules

# V3000 collection name -> label drawn next to the atom or bond
#   STEABS -> abs, STEREL1 -> or1, STERAC2 -> and2, STEBREL1 -> or1
COLL_RE = re.compile(r"MDLV30/STEB?(ABS|REL|RAC)(\d*) (ATOMS|BONDS)=\((\d+) ([\d ]+)\)")
COLL_LABEL = {"ABS": "abs", "REL": "or", "RAC": "and"}


# --------------------------------------------------------------- molblocks

def read_sdf_record(path, index):
    """Molblock 'index' (0-based) of an SDF, without its data fields."""
    records = open(os.path.join(CODE, path[3:])).read().split("$$$$\n")
    return records[index].split("M  END")[0] + "M  END\n"


def read_test_literal(path, anchor):
    """First run of C string literals after 'anchor' in a unit test."""
    src = open(os.path.join(CODE, path[3:])).read()
    tail = src[src.index(anchor):]
    tail = tail[:tail.index(";")]
    parts = re.findall(r'"((?:[^"\\]|\\.)*)"', tail)
    return "".join(parts).encode().decode("unicode_escape")


def double_bonds(mol):
    """Non-aromatic double bonds (C=O too) in bond order: the 'db' references."""
    return [b.GetIdx() for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE]


def from_smiles(entry):
    """Lay out a (CX)SMILES and write V3000 with the entry's extras."""
    mol = Chem.MolFromSmiles(entry["smiles"])
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    # Rotate the drawing: same molecule, different coordinates
    angle = entry.get("rotate", 0.0)
    if angle:
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x * math.cos(angle) - p.y * math.sin(angle),
                                     p.x * math.sin(angle) + p.y * math.cos(angle), 0))

    Chem.WedgeMolBonds(mol, mol.GetConformer())
    block = Chem.MolToV3KMolBlock(mol)
    lines = block.split("\n")
    dbs = double_bonds(mol)

    # Extra wedges (allene ends) and crossed double bonds, by bond line
    bond_cfg = {}
    for begin, end, cfg in entry.get("wedges", []):
        bond_cfg[mol.GetBondBetweenAtoms(begin, end).GetIdx() + 1] = (begin + 1, end + 1, cfg)
    for k in entry.get("either", []):
        b = mol.GetBondWithIdx(dbs[k])
        bond_cfg[dbs[k] + 1] = (b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1, 2)

    in_bonds = False
    for i, line in enumerate(lines):
        if line.startswith("M  V30 BEGIN BOND"):
            in_bonds = True
            continue
        if line.startswith("M  V30 END BOND"):
            break
        if not in_bonds:
            continue
        f = line.split()
        idx = int(f[2])
        if idx in bond_cfg:
            a1, a2, cfg = bond_cfg[idx]
            lines[i] = f"M  V30 {idx} {f[3]} {a1} {a2} CFG={cfg}"

    # Extra collections; atoms 1-based, bonds by double-bond rank or index
    extra = []
    for name, items in entry.get("coll", []):
        if name.startswith("STEB"):
            ids = [dbs[v] + 1 if kind == "db" else v for kind, v in items]
            extra.append(f"M  V30 MDLV30/{name} BONDS=({len(ids)} {' '.join(map(str, ids))})")
        else:
            ids = [a + 1 for a in items]
            extra.append(f"M  V30 MDLV30/{name} ATOMS=({len(ids)} {' '.join(map(str, ids))})")

    if extra:
        if "M  V30 END COLLECTION" in lines:
            at = lines.index("M  V30 END COLLECTION")
            lines[at:at] = extra
        else:
            at = lines.index("M  V30 END CTAB")
            lines[at:at] = ["M  V30 BEGIN COLLECTION"] + extra + ["M  V30 END COLLECTION"]

    lines[0] = entry["title"][:70]
    return "\n".join(lines)


def molblock(entry):
    if "smiles" in entry:
        return from_smiles(entry)
    if "sdf" in entry:
        block = read_sdf_record(*entry["sdf"])
    else:
        block = read_test_literal(*entry["tests"])
    for old, new in entry.get("edit", []):
        assert old in block, (entry["title"], old)
        block = block.replace(old, new)
    return block


# ------------------------------------------------------------------ InChI

def run_cli(exe, block, opts):
    """Output and log lines of one inchi-1 run."""
    with tempfile.TemporaryDirectory() as tmp:
        mol, out, log = (os.path.join(tmp, n) for n in ("in.mol", "out.txt", "log.txt"))
        open(mol, "w").write(block if block.endswith("\n") else block + "\n")
        subprocess.run([exe, mol, out, log, os.devnull, "-NoLabels"] + opts,
                       capture_output=True, timeout=60)
        return open(out).read().splitlines(), open(log).read().splitlines()


def run_inchi(exe, block, opts):
    """InChI string and canonical numbers {0-based input atom: number};
    the log's error line instead of the InChI when there is none.

    AuxInfo /N: lists input atoms in canonical order, per component:
      AuxInfo=1/1/N:4,1,2,3;5,6  ->  atom 4 is 1, atom 1 is 2, ... atom 6 is 2
    """
    out, log = run_cli(exe, block, opts)
    inchi = next((l for l in out if l.startswith("InChI=")), None)
    if inchi is None:
        error = next((l for l in log if l.startswith("Error")), "no output")
        return "(" + error + ")", {}

    canon = {}
    aux = next((l for l in out if l.startswith("AuxInfo=")), "")
    numbers = re.search(r"/N:([^/]*)", aux)
    for component in (numbers.group(1).split(";") if numbers else []):
        for rank, atom in enumerate(component.split(","), 1):
            canon[int(atom) - 1] = rank
    return inchi, canon


# -------------------------------------------------------------------- SVG

def collection_notes(block):
    """{('atom'|'bond', 0-based idx): label} from the V3000 collections."""
    notes = {}
    for kind, num, what, _, ids in COLL_RE.findall(block):
        label = COLL_LABEL[kind] + num
        target = "atom" if what == "ATOMS" else "bond"
        for i in ids.split():
            notes[(target, int(i) - 1)] = label
    return notes


def draw(block, canon, relayout=False):
    mol = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False, strictParsing=False)
    if mol is None:
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=False)
        mol.UpdatePropertyCache(strict=False)

    # 3D records (any z, whatever the header says): draw the x/y projection,
    # which shows the twist; RDKit does not wedge an axis from 3D
    conf = mol.GetConformer()
    if any(abs(conf.GetAtomPosition(i).z) > 1e-4 for i in range(mol.GetNumAtoms())):
        pass
    elif relayout:
        # Source drawing overlaps: new 2D layout, RDKit re-wedges the
        # stereo it perceived (atropisomer axes included)
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)
        Chem.WedgeMolBonds(mol, mol.GetConformer())
    else:
        Chem.ReapplyMolBlockWedging(mol, allBondTypes=True)

    # Explicit H after all heavy atoms: hide them, atom numbers stay valid
    hs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    heavy = mol.GetNumAtoms() - len(hs)
    if hs and min(hs) >= heavy:
        mol = Chem.RemoveHs(mol, sanitize=False)

    # Labels come from the file's collections, not RDKit's stereo groups,
    # which miss allene axes and bond groups
    mol = Chem.RWMol(mol)
    mol.SetStereoGroups([])
    for (target, i), label in collection_notes(block).items():
        if target == "atom" and i < mol.GetNumAtoms():
            mol.GetAtomWithIdx(i).SetProp("atomNote", label)
        if target == "bond" and i < mol.GetNumBonds():
            mol.GetBondWithIdx(i).SetProp("bondNote", label)

    d = rdMolDraw2D.MolDraw2DSVG(SVG_W, SVG_H)
    o = d.drawOptions()
    o.addStereoAnnotation = False
    o.annotationFontScale = 0.8
    o.clearBackground = False
    o.padding = 0.1  # room for the canonical numbers at the edge
    o.maxFontSize = ATOM_FONT_MAX_PX
    mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True, addChiralHs=False,
                                           wedgeBonds=False)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    svg = svg[svg.index("<svg"):]
    return svg.replace("</svg>", canon_labels(d, mol, canon) + "</svg>")


def has_symbol(atom, pos):
    """Does RDKit print a symbol for this atom? Heteroatoms, lone atoms
    and linear carbons (allene centre) get one."""
    nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
    if atom.GetAtomicNum() != 6 or not nbrs:
        return True
    if len(nbrs) != 2:
        return False
    a, b, c = pos[nbrs[0]], pos[atom.GetIdx()], pos[nbrs[1]]
    u, v = (a[0] - b[0], a[1] - b[1]), (c[0] - b[0], c[1] - b[1])
    cos = (u[0] * v[0] + u[1] * v[1]) / (math.hypot(*u) * math.hypot(*v) or 1)
    return cos < -0.98


def canon_labels(d, mol, canon):
    """Small grey canonical numbers beside each atom, on the side away
    from its bonds:

            OH 5
            |
        1 - C 4 ...      offset 8 px (14 px next to a symbol), less on small drawings
    """
    pos = {a.GetIdx(): tuple(d.GetDrawCoords(a.GetIdx())) for a in mol.GetAtoms()}
    lengths = [math.dist(pos[b.GetBeginAtomIdx()], pos[b.GetEndAtomIdx()])
               for b in mol.GetBonds()]
    bond = sorted(lengths)[len(lengths) // 2] if lengths else 30.0

    out = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        if i not in canon:
            continue

        # Direction: away from the mean of the bond directions; a linear
        # or lone atom goes perpendicular / down-right
        x, y = pos[i]
        dx = dy = 0.0
        for n in atom.GetNeighbors():
            nx, ny = pos[n.GetIdx()]
            length = math.hypot(nx - x, ny - y) or 1
            dx -= (nx - x) / length
            dy -= (ny - y) / length
        if math.hypot(dx, dy) < 0.3:
            if atom.GetDegree():
                nx, ny = pos[atom.GetNeighbors()[0].GetIdx()]
                dx, dy = -(ny - y), nx - x
            else:
                dx, dy = 1.0, 1.0
        norm = math.hypot(dx, dy)
        r = min(bond * (0.45 if has_symbol(atom, pos) else 0.15),
                CANON_GAP_SYMBOL_PX if has_symbol(atom, pos) else CANON_GAP_PX)
        lx = min(max(x + dx / norm * r, CANON_FONT_PX), SVG_W - CANON_FONT_PX)
        ly = min(max(y + dy / norm * r, CANON_FONT_PX), SVG_H - CANON_FONT_PX)
        out.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="{CANON_FONT_PX}" '
                   'fill="#666" font-family="sans-serif" text-anchor="middle" '
                   f'dominant-baseline="central">{canon[i]}</text>')
    return "".join(out)


# ------------------------------------------------------------------- HTML

def mark_diff(old, new):
    """Bold the part of 'new' after the common prefix with 'old' (after
    the version prefix), so the change is visible at a glance."""
    if old == new:
        return html.escape(new)
    a, b = old[9:], new[9:]
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    return html.escape(new[:9 + n]) + "<b>" + html.escape(new[9 + n:]) + "</b>"


PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>InChI enhanced stereochemistry: old vs. new</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2rem auto; max-width: 1200px;
       padding: 0 1rem; color: #222; }}
nav a {{ margin-right: 1rem; }}
table {{ border-collapse: collapse; width: 100%; margin-bottom: 3rem; }}
th, td {{ border: 1px solid #ccc; padding: .4rem .6rem; vertical-align: middle; }}
th {{ background: #f3f3f3; text-align: left; }}
td.id {{ text-align: right; width: 3rem; }}
td.pic {{ width: {w}px; text-align: center; font-size: .85rem; color: #555; }}
td.pic img {{ display: block; margin: 0 auto .3rem; }}
.inchi {{ font-family: ui-monospace, monospace; font-size: .85rem; word-break: break-all; }}
.inchi div {{ margin: .2rem 0; }}
.tag {{ display: inline-block; width: 2.6rem; color: #888; }}
.same {{ color: #888; }}
b {{ color: #b00; }}
</style>
</head>
<body>
<h1>InChI enhanced stereochemistry: old vs. new</h1>
<p>Four extensions specified in
<a href="https://github.com/cm-beilstein/InChI-specs">InChI-specs</a> (private).
<b>old</b>: InChI <code>dev</code> ({old_ver}), default options.
<b>new</b>: branch <code>atropisomers</code> ({new_ver}, contains all four features),
<code>-EnhancedStereochemistry</code>. The changed tail of the new InChI is in red.
Small grey numbers in the pictures: canonical atom numbers of the new InChI
(per component). Collection labels <code>abs</code>, <code>or<i>n</i></code>,
<code>and<i>n</i></code>. Generated by <code>tools/build.py</code>.</p>
<nav>{nav}</nav>
{tables}
</body>
</html>
"""


def build():
    old_exe, new_exe = os.environ["INCHI_OLD"], os.environ["INCHI_NEW"]
    tables, nav = [], []
    ident = 0

    for key, title, entries in TABLES:
        rows, sdf = [], []
        for entry in entries:
            ident += 1
            block = molblock(entry)
            old, _ = run_inchi(old_exe, block, OLD_OPTS)
            new, canon = run_inchi(new_exe, block, NEW_OPTS)

            svg_name = f"svg/{key}-{ident}.svg"
            svg = draw(block, canon, entry.get("relayout", False))
            open(os.path.join(ROOT, svg_name), "w").write(svg)

            # Content hash in the URL: browsers reload a changed picture
            # despite GitHub Pages' 10-minute cache
            svg_url = svg_name + "?v=" + hashlib.sha1(svg.encode()).hexdigest()[:8]
            sdf.append(block.rstrip("\n") + f"\n>  <ID>\n{ident}\n\n>  <OLD>\n{old}\n\n"
                       f">  <NEW>\n{new}\n\n$$$$")

            same = old[9:] == new[9:]
            rows.append(
                f'<tr><td class="id">{ident}</td>'
                f'<td class="pic"><img src="{svg_url}" width="{SVG_W}" height="{SVG_H}" '
                f'alt="molecule {ident}">{html.escape(entry["title"])}</td>'
                f'<td class="inchi"><div><span class="tag">old</span>{html.escape(old)}</div>'
                f'<div{" class=same" if same else ""}><span class="tag">new</span>'
                f'{mark_diff(old, new)}</div></td></tr>')
            print(f"{ident:3} {key:13} {new}", file=sys.stderr)

        open(os.path.join(ROOT, f"data/{key}.sdf"), "w").write("\n".join(sdf) + "\n")
        nav.append(f'<a href="#{key}">{html.escape(title)}</a>')
        tables.append(f'<h2 id="{key}">{html.escape(title)} ({len(entries)})</h2>\n'
                      '<table><tr><th>ID</th><th>Molecule</th><th>InChI</th></tr>\n'
                      + "\n".join(rows) + "\n</table>")

    page = PAGE.format(w=SVG_W + 20, nav="".join(nav), tables="\n".join(tables),
                       old_ver=os.environ.get("INCHI_OLD_REV", "dev"),
                       new_ver=os.environ.get("INCHI_NEW_REV", "atropisomers"))
    open(os.path.join(ROOT, "index.html"), "w").write(page)


if __name__ == "__main__":
    build()
