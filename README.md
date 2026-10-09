# InChI enhanced stereochemistry: old vs. new

Static page comparing InChIs from the current InChI (`dev`) with the
prototype for four stereo extensions:

1. Enhanced stereochemistry (ABS / OR / AND groups)
2. Allenes
3. Double-bond E/Z
4. Atropisomers

Page: https://cm-beilstein.github.io/InChI-enh-Stereo/

Specs: https://github.com/cm-beilstein/InChI-specs (private)

## Contents

- `index.html` — the page
- `svg/` — molecule pictures (RDKit)
- `data/<table>.sdf` — input molblocks with old and new InChI
- `tools/molecules.py` — molecule sets
- `tools/build.py` — generator

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
