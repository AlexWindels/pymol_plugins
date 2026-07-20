## show_contacts

What it does: a PyMOL plugin that visualizes the interaction interface between the two chains (A = target, H = binder) of a structure. For a selected object — or all loaded objects at once — it draws the inter-chain contacts as color-coded dashes: polar/H-bonds (yellow), hydrophobic C–C (orange), salt bridges (blue), and π-stacking (purple), each with an adjustable distance cutoff. Contacts are computed strictly between chains A and H, so only interface interactions are shown. It optionally highlights the participating residues as labelled sticks and annotates each contact with its distance. 

Usage: Use it from the Plugin → Show Contacts… dialog or the show_contacts command.

## cealign_all

What it does: CE-aligns every object in the session onto a reference (default = the first loaded object).

Usage:
- cealign_all — align all objects onto the first one
- cealign_all ref=<object> — align all onto a chosen reference
- cealign_all sel=chain A — superpose using only chain A
