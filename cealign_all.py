"""
cealign_all - CE-align every object in the PyMOL session onto a reference
(default: the first loaded object).

Install: PyMOL > Plugin > Plugin Manager > Install New Plugin > (this file).
Or just: run "C:/.../cealign_all.py"

Command-line usage in PyMOL:
    cealign_all                     # align all objects onto the first one
    cealign_all ref=2vyr_1          # align all onto object 'ref'
    cealign_all sel=chain A         # superpose using only chain A (e.g. the shared target)
    cealign_all ref=2vyr_1, sel=chain A
"""

from pymol import cmd


def cealign_all(ref="", sel="", quiet=0):
    """Superpose all objects onto `ref` (default: first object) via cealign.

    ref  : reference object name; empty = first object in the session.
    sel  : optional selection to align on (applied to both), e.g. "chain A".
    quiet: 1 to suppress per-object output.
    """
    quiet = int(quiet)
    objects = cmd.get_object_list()
    if not objects:
        print("cealign_all: no objects loaded.")
        return

    ref = ref or objects[0]
    if ref not in objects:
        print("cealign_all: reference '%s' not in session (%s)." % (ref, ", ".join(objects)))
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
                print("cealign_all: %-32s RMSD=%.3f over %d residues"
                      % (obj, r["RMSD"], r["alignment_length"]))
        except Exception as e:
            print("cealign_all: failed on %s: %s" % (obj, e))

    if not quiet:
        print("cealign_all: aligned %d object(s) onto '%s'." % (n, ref))


cmd.extend("cealign_all", cealign_all)


def __init_plugin__(app=None):
    """Add a Plugin menu entry (PyMOL 2.x Qt)."""
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("CEalign all onto first", lambda: cealign_all())
    except Exception:
        pass
