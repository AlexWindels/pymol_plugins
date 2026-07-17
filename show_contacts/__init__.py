"""
Show Contacts - PyMOL plugin: draw interface contacts (polar/H-bond,
hydrophobic C-C, salt-bridge, pi-stacking) between the two chains (A and H)
of an object, via a dialog (or a command).

Assumes every object contains exactly two chains, A (target) and H (binder).
Contacts are always computed BETWEEN chain A and chain H of the same object,
so they are interface-only (never intra-chain / whole-protein).

Optionally highlights the residues involved (sticks + residue labels) and
labels each contact with its distance.

Install: Plugin > Plugin Manager > Install New Plugin > File... and pick this
         __init__.py (PyMOL installs the whole folder). Or zip the folder.
Menu:    Plugin > "Show Contacts..."   -> opens the options dialog
Command: show_contacts [obj=all][, polar=1][, hydrophobic=1][, salt=1][, pi=1]
                       [, polar_cutoff=3.6][, hydrophobic_cutoff=4.0]
                       [, salt_cutoff=4.0][, pi_cutoff=5.5]
                       [, show_res=1][, show_dist=1][, chain1=A][, chain2=H]
         e.g.  show_contacts                 # all objects
               show_contacts obj=2734_0_RP   # one object

Per object <O>, results are grouped under <O>_contacts:
  <O>_contacts_polar        yellow  - polar/H-bond contacts (PyMOL distance mode 2)
  <O>_contacts_hydrophobic  orange  - nonpolar carbon-carbon pairs within cutoff
  <O>_contacts_saltbridge   blue    - charged N (Lys/Arg/His) <-> O (Asp/Glu)
  <O>_contacts_pistack      purple  - aromatic ring-centroid pairs within cutoff
and the involved residues are selected as <O>_contacts_res (shown as sticks).
"""

from pymol import cmd

# Nonpolar carbons: carbons NOT bonded to N/O (drops backbone C=O and CA, and
# polar side-chain carbons), leaving aliphatic/aromatic hydrophobic carbons.
_HYDROPHOBIC = "(elem C and not (neighbor (elem N+O)))"

# Charged groups for salt bridges.
_CATION = ("((resn LYS and name NZ) or (resn ARG and name NH1+NH2+NE) "
           "or (resn HIS and name ND1+NE2))")
_ANION = "((resn ASP and name OD1+OD2) or (resn GLU and name OE1+OE2))"
_AROMATIC = "(resn PHE+TYR+TRP+HIS)"

# Aromatic ring atoms per residue (for pi-stacking centroids).
_RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # 6-membered ring
}

_SUFFIXES = ("_polar", "_hydrophobic", "_saltbridge", "_saltbridge_rev",
             "_pistack", "_pi_a", "_pi_b")

_ALL = "(all objects)"


def _ring_centroids(sel):
    """Return [(chain, resi, (x,y,z)), ...] ring centroids for aromatics in sel."""
    model = cmd.get_model("(%s) and %s" % (sel, _AROMATIC))
    byres = {}
    for at in model.atom:
        byres.setdefault((at.segi, at.chain, at.resi, at.resn), {})[at.name] = at.coord
    out = []
    for (segi, chain, resi, resn), atoms in byres.items():
        coords = [atoms[n] for n in _RINGS.get(resn, []) if n in atoms]
        if len(coords) < 3:
            continue
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        out.append((chain, resi, (cx, cy, cz)))
    return out


def _style(obj, color, gap=None, width=2, labels=False):
    cmd.set("dash_color", color, obj)
    cmd.set("dash_width", width, obj)
    if gap is not None:
        cmd.set("dash_gap", gap, obj)
    (cmd.show if labels else cmd.hide)("labels", obj)  # distance value labels


def _within(a, cutoff, b):
    return "byres ((%s) within %g of (%s))" % (a, cutoff, b)


def _clear(name):
    rsel = name + "_res"
    if rsel in cmd.get_names("selections"):
        cmd.label("(%s) and name CA" % rsel, "")
        cmd.hide("sticks", rsel)
        cmd.delete(rsel)
    for sfx in _SUFFIXES:
        cmd.delete(name + sfx)
    cmd.delete(name)


def _contacts_between(sel1, sel2, name, polar, hydrophobic, salt, pi,
                      pc, hc, sc, pic, show_res, show_dist):
    """Draw requested contact types between sel1/sel2 (disjoint chains)."""
    _clear(name)
    made = []
    res_exprs = []  # selections of residues involved in the contacts

    if polar:
        d1, d2 = ("(%s) and (donor or acceptor)" % s for s in (sel1, sel2))
        pn = name + "_polar"
        cmd.distance(pn, sel1, sel2, pc, mode=2)  # PyMOL H-bond detection
        _style(pn, "yellow", width=3, labels=show_dist)
        made.append(pn)
        res_exprs += [_within(d1, pc, d2), _within(d2, pc, d1)]

    if hydrophobic:
        h1, h2 = ("(%s) and %s" % (s, _HYDROPHOBIC) for s in (sel1, sel2))
        hn = name + "_hydrophobic"
        cmd.distance(hn, h1, h2, hc, mode=0)
        _style(hn, "orange", gap=0.5, labels=show_dist)
        made.append(hn)
        res_exprs += [_within(h1, hc, h2), _within(h2, hc, h1)]

    if salt:
        c1, a1 = ("(%s) and %s" % (sel1, _CATION), "(%s) and %s" % (sel1, _ANION))
        c2, a2 = ("(%s) and %s" % (sel2, _CATION), "(%s) and %s" % (sel2, _ANION))
        for obj, x, y in [(name + "_saltbridge", c1, a2),
                          (name + "_saltbridge_rev", a1, c2)]:
            if cmd.count_atoms(x) and cmd.count_atoms(y):
                cmd.distance(obj, x, y, sc, mode=0)
                _style(obj, "marine", gap=0.3, width=3, labels=show_dist)
                made.append(obj)
        res_exprs += [_within(c1, sc, a2), _within(a2, sc, c1),
                      _within(a1, sc, c2), _within(c2, sc, a1)]

    if pi:
        cen1, cen2 = _ring_centroids(sel1), _ring_centroids(sel2)
        if cen1 and cen2:
            a, b = name + "_pi_a", name + "_pi_b"
            for chain, resi, pos in cen1:
                cmd.pseudoatom(a, pos=list(pos), name="CEN", chain=chain, resi=resi)
            for chain, resi, pos in cen2:
                cmd.pseudoatom(b, pos=list(pos), name="CEN", chain=chain, resi=resi)
            pn = name + "_pistack"
            cmd.distance(pn, a, b, pic, mode=0)
            _style(pn, "purple", width=3, labels=show_dist)
            cmd.hide("everything", a)
            cmd.hide("everything", b)
            made += [pn, a, b]
        ar1 = "(%s) and %s" % (sel1, _AROMATIC)
        ar2 = "(%s) and %s" % (sel2, _AROMATIC)
        res_exprs += [_within(ar1, pic, ar2), _within(ar2, pic, ar1)]

    if made:
        cmd.group(name, " ".join(made))

    # highlight the residues involved (sticks) + label them by residue
    if show_res and res_exprs:
        rsel = name + "_res"
        cmd.select(rsel, " or ".join("(%s)" % e for e in res_exprs))
        if cmd.count_atoms(rsel):
            cmd.show("sticks", rsel)
            cmd.label("(%s) and name CA" % rsel, "resn+resi")
            cmd.disable(rsel)  # hide the pink selection dots; keep the sticks
        else:
            cmd.delete(rsel)

    return made


def show_contacts(obj="all", polar=1, hydrophobic=1, salt=1, pi=1,
                  polar_cutoff=3.6, hydrophobic_cutoff=4.0,
                  salt_cutoff=4.0, pi_cutoff=5.5,
                  show_res=1, show_dist=1, chain1="A", chain2="H", quiet=0):
    """Interface contacts between chain `chain1` and `chain2` within object(s)."""
    polar, hydrophobic, salt, pi = int(polar), int(hydrophobic), int(salt), int(pi)
    show_res, show_dist, quiet = int(show_res), int(show_dist), int(quiet)
    pc, hc = float(polar_cutoff), float(hydrophobic_cutoff)
    sc, pic = float(salt_cutoff), float(pi_cutoff)

    objs = _molecule_objects() if obj in ("all", _ALL, "") else [obj]
    if not objs:
        print("show_contacts: no molecule objects loaded.")
        return

    cmd.set("label_distance_digits", 2)  # 2-decimal distances

    for o in objs:
        sel1 = "(%s) and chain %s" % (o, chain1)
        sel2 = "(%s) and chain %s" % (o, chain2)
        if cmd.count_atoms(sel1) == 0 or cmd.count_atoms(sel2) == 0:
            print("show_contacts: %s missing chain %s/%s (%d/%d atoms) - skipped."
                  % (o, chain1, chain2, cmd.count_atoms(sel1), cmd.count_atoms(sel2)))
            continue
        made = _contacts_between(sel1, sel2, o + "_contacts",
                                 polar, hydrophobic, salt, pi,
                                 pc, hc, sc, pic, show_res, show_dist)
        if not quiet:
            shown = ", ".join(m for m in made if not m.endswith(("_pi_a", "_pi_b")))
            print("show_contacts: %s (chain %s <-> %s) -> %s"
                  % (o, chain1, chain2, shown or "nothing"))


cmd.extend("show_contacts", show_contacts)


# --------------------------------------------------------------------------- #
# GUI
# --------------------------------------------------------------------------- #
dialog = None  # module ref so the dialog isn't garbage-collected


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Show Contacts...", run_plugin_gui)


def _molecule_objects():
    """Only real structures, so generated contact objects don't show up."""
    return [o for o in cmd.get_object_list()
            if cmd.get_type(o) == "object:molecule"]


def _refresh_objects():
    combo = dialog._obj
    keep = combo.currentText()
    combo.clear()
    combo.addItem(_ALL)
    for o in _molecule_objects():
        combo.addItem(o)
    i = combo.findText(keep)
    combo.setCurrentIndex(i if i >= 0 else 0)


def run_plugin_gui():
    from pymol.Qt import QtWidgets
    global dialog
    if dialog is None:
        dialog = _build_dialog(QtWidgets)
    _refresh_objects()
    dialog.show()
    dialog.raise_()


def _contact_row(QtWidgets, label, default_cut, cut_min=1.0, cut_max=8.0):
    chk = QtWidgets.QCheckBox(label)
    chk.setChecked(True)
    spin = QtWidgets.QDoubleSpinBox()
    spin.setRange(cut_min, cut_max)
    spin.setSingleStep(0.1)
    spin.setValue(default_cut)
    row = QtWidgets.QHBoxLayout()
    row.addWidget(chk)
    row.addWidget(QtWidgets.QLabel("cutoff (A):"))
    row.addWidget(spin)
    return row, chk, spin


def _build_dialog(QtWidgets):
    d = QtWidgets.QDialog()
    d.setWindowTitle("Show Contacts")
    lay = QtWidgets.QVBoxLayout(d)

    lay.addWidget(QtWidgets.QLabel(
        "Interface contacts between chain A and chain H of the chosen object\n"
        "(or all objects). Each object is assumed to hold chains A and H."))

    d._obj = QtWidgets.QComboBox()
    form = QtWidgets.QFormLayout()
    form.addRow("Object:", d._obj)
    lay.addLayout(form)

    prow, polar, polar_cut = _contact_row(QtWidgets, "Polar / H-bonds (yellow)", 3.6)
    hrow, hydro, hydro_cut = _contact_row(QtWidgets, "Hydrophobic C-C (orange)", 4.0)
    srow, salt, salt_cut = _contact_row(QtWidgets, "Salt bridges (blue)", 4.0)
    pirow, pi, pi_cut = _contact_row(QtWidgets, "Pi-stacking (purple)", 5.5)
    for r in (prow, hrow, srow, pirow):
        lay.addLayout(r)

    show_res = QtWidgets.QCheckBox("Show contact residues (sticks + residue labels)")
    show_res.setChecked(True)
    show_dist = QtWidgets.QCheckBox("Label distances on contacts")
    show_dist.setChecked(True)
    lay.addWidget(show_res)
    lay.addWidget(show_dist)

    btns = QtWidgets.QHBoxLayout()
    show_btn = QtWidgets.QPushButton("Show")
    clear_btn = QtWidgets.QPushButton("Clear")
    close_btn = QtWidgets.QPushButton("Close")
    for b in (show_btn, clear_btn, close_btn):
        btns.addWidget(b)
    lay.addLayout(btns)

    def _targets():
        sel = d._obj.currentText()
        return _molecule_objects() if sel in ("", _ALL) else [sel]

    def on_show():
        show_contacts(
            obj=d._obj.currentText() or "all",
            polar=1 if polar.isChecked() else 0,
            hydrophobic=1 if hydro.isChecked() else 0,
            salt=1 if salt.isChecked() else 0,
            pi=1 if pi.isChecked() else 0,
            polar_cutoff=polar_cut.value(),
            hydrophobic_cutoff=hydro_cut.value(),
            salt_cutoff=salt_cut.value(),
            pi_cutoff=pi_cut.value(),
            show_res=1 if show_res.isChecked() else 0,
            show_dist=1 if show_dist.isChecked() else 0,
        )

    def on_clear():
        for o in _targets():
            _clear(o + "_contacts")

    show_btn.clicked.connect(on_show)
    clear_btn.clicked.connect(on_clear)
    close_btn.clicked.connect(d.close)
    return d
