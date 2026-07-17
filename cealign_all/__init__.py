"""
CEalign All - PyMOL plugin: CE-align every object in the session onto a
reference, via a dialog (or a command).

Install: Plugin > Plugin Manager > Install New Plugin > File... and pick this
         __init__.py (PyMOL installs the whole folder). Or zip this folder and
         install the .zip.
Menu:    Plugin > "CEalign All..."   -> opens the options dialog
Command: cealign_all [ref=<obj>][, sel=<selection>][, quiet=0/1]
         e.g.  cealign_all sel=chain A
"""

from pymol import cmd


# --------------------------------------------------------------------------- #
# Core: also usable from the PyMOL command line as `cealign_all`
# --------------------------------------------------------------------------- #
def cealign_all(ref="", sel="", quiet=0):
    """Superpose all objects onto `ref` (default: first object) via cealign.

    ref  : reference object name; empty = first object in the session.
    sel  : optional selection to align on, applied to both, e.g. "chain A".
    quiet: 1 to suppress per-object output.
    """
    quiet = int(quiet)
    objects = cmd.get_object_list()
    if not objects:
        print("cealign_all: no objects loaded.")
        return

    ref = ref or objects[0]
    if ref not in objects:
        print("cealign_all: reference '%s' not in session (%s)."
              % (ref, ", ".join(objects)))
        return

    n = 0
    for obj in objects:
        if obj == ref:
            continue
        target = "(%s) and (%s)" % (ref, sel) if sel else ref
        mobile = "(%s) and (%s)" % (obj, sel) if sel else obj
        try:
            r = cmd.cealign(target, mobile)
            n += 1
            if not quiet:
                print("cealign_all: %-32s RMSD=%.3f over %d res"
                      % (obj, r["RMSD"], r["alignment_length"]))
        except Exception as e:
            print("cealign_all: failed on %s: %s" % (obj, e))

    if not quiet:
        print("cealign_all: aligned %d object(s) onto '%s'." % (n, ref))
    return n


cmd.extend("cealign_all", cealign_all)


# --------------------------------------------------------------------------- #
# GUI
# --------------------------------------------------------------------------- #
dialog = None  # module-level ref so the dialog isn't garbage-collected


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("CEalign All...", run_plugin_gui)


def run_plugin_gui():
    from pymol.Qt import QtWidgets

    global dialog
    if dialog is None:
        dialog = _build_dialog(QtWidgets)

    # refresh the reference dropdown with the currently loaded objects
    combo = dialog._ref_combo
    keep = combo.currentText()
    combo.clear()
    combo.addItem("(first object)")
    for o in cmd.get_object_list():
        combo.addItem(o)
    i = combo.findText(keep)
    if i >= 0:
        combo.setCurrentIndex(i)

    dialog.show()
    dialog.raise_()


def _build_dialog(QtWidgets):
    d = QtWidgets.QDialog()
    d.setWindowTitle("CEalign All")
    lay = QtWidgets.QVBoxLayout(d)

    lay.addWidget(QtWidgets.QLabel("Reference (align everything onto this):"))
    ref_combo = QtWidgets.QComboBox()
    lay.addWidget(ref_combo)

    chain_a = QtWidgets.QCheckBox("Align using chain A only (shared target)")
    chain_a.setChecked(True)
    lay.addWidget(chain_a)

    lay.addWidget(QtWidgets.QLabel("Custom selection (optional; overrides the checkbox):"))
    sel_edit = QtWidgets.QLineEdit()
    sel_edit.setPlaceholderText("blank = whole object; e.g. chain A and resi 1-120")
    lay.addWidget(sel_edit)

    quiet = QtWidgets.QCheckBox("Quiet (suppress per-object output)")
    lay.addWidget(quiet)

    btns = QtWidgets.QHBoxLayout()
    run_btn = QtWidgets.QPushButton("Align")
    close_btn = QtWidgets.QPushButton("Close")
    btns.addWidget(run_btn)
    btns.addWidget(close_btn)
    lay.addLayout(btns)

    def on_run():
        ref = ref_combo.currentText()
        ref = "" if ref == "(first object)" else ref
        sel = sel_edit.text().strip()
        if not sel and chain_a.isChecked():
            sel = "chain A"
        cealign_all(ref=ref, sel=sel, quiet=1 if quiet.isChecked() else 0)

    run_btn.clicked.connect(on_run)
    close_btn.clicked.connect(d.close)

    d._ref_combo = ref_combo  # stash for refresh in run_plugin_gui
    return d
