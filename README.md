# InChI enhanced stereochemistry: old vs. new

Static page comparing InChIs from the current InChI (`dev`) with the
prototype for four stereo extensions:

1. Enhanced stereochemistry (ABS / OR / AND groups)
2. Allenes
3. Double-bond E/Z
4. Atropisomers

plus a table of meso compounds, where a group changes nothing
([#286](https://github.com/IUPAC-InChI/InChI/issues/286)), with its
known limitation.

Page: https://cm-beilstein.github.io/InChI-enh-Stereo/

Specs: https://github.com/cm-beilstein/InChI-specs (private)

## Explanations for chemists

What the labels mean and how they reach the InChI, one page per feature
([source](https://github.com/cm-beilstein/InChI-enh-Stereo/tree/main/explain)):

- [Enhanced stereochemistry](https://cm-beilstein.github.io/InChI-enh-Stereo/explain/enhanced.html)
  — incl. [meso compounds](https://cm-beilstein.github.io/InChI-enh-Stereo/explain/enhanced.html#meso)
- [Allenes](https://cm-beilstein.github.io/InChI-enh-Stereo/explain/allenes.html)
- [Double-bond E/Z](https://cm-beilstein.github.io/InChI-enh-Stereo/explain/ez.html)
- [Atropisomers](https://cm-beilstein.github.io/InChI-enh-Stereo/explain/atropisomers.html)

## Contents

- `index.html` — the page
- `svg/` — molecule pictures (RDKit)
- `data/<table>.sdf` — input molblocks with old and new InChI
- `tools/molecules.py` — molecule sets
- `tools/build.py` — generator
- `explain/<table>.html` — walkthrough per spec for chemists; copied from
  InChI-Specifications `<spec>/explain.html` (built there by `tools/explain.py`)
- `workflow/<table>.html` — how the new InChI is built, per spec, for chemists;
  rendered from `workflow/src/*.json` with [archify](https://github.com/tt-a1i/archify) (MIT)

## Rebuild

Needs RDKit and two `inchi-1` binaries:

```
INCHI_OLD=/path/dev/inchi-1 INCHI_NEW=/path/atropisomers/inchi-1 \
    python tools/build.py
```

`INCHI_OLD` runs with default options, `INCHI_NEW` with
`-EnhancedStereochemistry`. `tools/molecules.py` reads some inputs
from sibling checkouts `../InChI` (branch `atropisomers`) and
`../InChI-specs`.

Rebuild a diagram (archify checked out next to this repo):

```
node ../archify/archify/bin/archify.mjs finalize workflow \
    workflow/src/ez.json workflow/ez.html --quality showcase
```
