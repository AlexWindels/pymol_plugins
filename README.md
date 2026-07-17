## show_contacts

Shows **polar** contacts.

Good hydrogen bonds (as determined by PyMOL) are shown in yellow. 
Electrostatic clashes (donor-donor or acceptor-acceptor) are shown in red. 
Close (<4.0 A) but not ideal contacts are shown in purple. 

Based on [show contacts](https://pymolwiki.org/index.php/Show_contacts)

## cealign_all

What it does: CE-aligns every object in the session onto a reference (default = the first loaded object).

Usage:
- cealign_all — align all objects onto the first one
- cealign_all ref=<object> — align all onto a chosen reference
- cealign_all sel=chain A — superpose using only chain A
