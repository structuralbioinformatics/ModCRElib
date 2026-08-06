#!/usr/bin/env python3
"""
mutate_bases.py - standalone reimplementation of the classic X3DNA
`mutate_bases` utility (the compiled binary no longer works against
modern X3DNA/DSSR installs).

USAGE
-----
Generate a mutation-entry template listing every nucleotide residue
in a structure:

    python3 mutate_bases.py -e input.pdb [template.dat]

Apply mutations described in a control file to a structure:

    python3 mutate_bases.py input.pdb [mutations.dat] -o output.pdb

MUTATION FILE FORMAT
---------------------
One entry per line. Lines starting with '#' (after stripping leading
whitespace) are comments and are ignored, as are blank lines.

    chain=A snum=2 name=C mutation=G
    chain=B snum=15 m=A

Fields (key=value, whitespace separated):
    chain     PDB chain ID                                    (required)
    snum      PDB residue/sequence number  (alias: seqnum)     (required)
    icode     PDB insertion code                               (optional)
    name      original base name - sanity-checked if given     (optional)
    mutation  target base: A, C, G, T or U   (alias: m)        (required)

ALGORITHM
---------
All five standard ribo/deoxyribo bases (A, C, G, T, U) share atom
names N1, C2, N3, C4, C5, C6 for their six-membered ring (a purine's
ring is fused to the imidazole ring N7/C8/N9, but it still carries
this same six-atom set). That means a single, base-type-agnostic
procedure works for every mutation, including purine<->pyrimidine
swaps:

  1. Read the six shared ring atoms from the residue being mutated.
  2. Kabsch-superpose an idealized copy of the TARGET base's ring
     atoms onto those same six atoms.
  3. Apply that rotation + translation to every atom of the idealized
     target base.
  4. Remove the original base atoms (keep sugar/phosphate untouched),
     insert the transformed target-base atoms, rename the residue.

NOTE ON GEOMETRY
----------------
IDEAL_BASES below holds approximate idealized planar ring geometry
(generic aromatic bond lengths/angles), good enough for most
downstream uses (visualization, distance-based analysis, further
refinement/energy minimization). For crystallographic-grade ideal
coordinates, swap IDEAL_BASES for the PDB Chemical Component
Dictionary's ideal coordinates for DA/DC/DG/DT/DU (or A/C/G/U for
RNA) - freely downloadable from RCSB - keeping the same atom-name
keys used here.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Idealized base geometry (approximate planar ring coordinates, angstroms).
# Origin/orientation are arbitrary here - only *relative* geometry matters,
# since everything gets rigidly transformed onto the real residue's frame.
# Ring atoms (N1,C2,N3,C4,C5,C6) MUST be present for every base, since the
# superposition is always done on that shared set.
# --------------------------------------------------------------------------

def _hexagon(bond=1.39):
    """Regular six-membered ring in the xy-plane, centered at origin."""
    names = ["N1", "C2", "N3", "C4", "C5", "C6"]
    r = bond / (2 * math.sin(math.pi / 6))
    pts = {}
    for i, name in enumerate(names):
        theta = math.radians(60 * i)
        pts[name] = np.array([r * math.cos(theta), r * math.sin(theta), 0.0])
    return pts


def _add(base, name, anchor, bond):
    """Place an exocyclic substituent `bond` A outward from `anchor`,
    along the radial direction from the ring center (origin) through
    `anchor`. Fine for a planar, roughly-regular idealized ring."""
    a = base[anchor]
    d = a / np.linalg.norm(a)
    base[name] = a + d * bond


def _pyrimidine_ring():
    return _hexagon()


def _purine_ring():
    """Six-membered ring + fused five-membered imidazole ring
    (N7, C8, N9) sharing the C4-C5 edge, built in a local frame
    (x = along the C4->C5 edge, y = outward, away from the ring
    center) so bond lengths come out chemically sane by construction:
    C4-N9 = C5-N7 ~= 1.37 A, N9-C8 = N7-C8 ~= 1.35 A, N9...N7 = 2.23 A.
    """
    base = _hexagon()
    c4, c5 = base["C4"], base["C5"]
    edge_dir = (c5 - c4)
    edge_dir = edge_dir / np.linalg.norm(edge_dir)
    center = np.zeros(3)
    outward = ((c4 + c5) / 2) - center
    outward = outward / np.linalg.norm(outward)
    perp = np.array([-edge_dir[1], edge_dir[0], 0.0])
    if np.dot(perp, outward) < 0:
        perp = -perp
    mid = (c4 + c5) / 2

    def local(x, y):
        return mid + edge_dir * x + perp * y

    base["N9"] = local(-1.115, 1.304)
    base["N7"] = local(1.115, 1.304)
    base["C8"] = local(0.0, 2.065)
    return base


def _build_ideal_bases():
    bases = {}

    # --- Cytosine ---
    c = _pyrimidine_ring()
    _add(c, "O2", "C2", 1.24)
    _add(c, "N4", "C4", 1.33)
    bases["C"] = c

    # --- Uracil ---
    u = _pyrimidine_ring()
    _add(u, "O2", "C2", 1.22)
    _add(u, "O4", "C4", 1.23)
    bases["U"] = u

    # --- Thymine (uracil + C5 methyl) ---
    t = _pyrimidine_ring()
    _add(t, "O2", "C2", 1.22)
    _add(t, "O4", "C4", 1.23)
    _add(t, "C7", "C5", 1.50)
    bases["T"] = t

    # --- Adenine ---
    a = _purine_ring()
    _add(a, "N6", "C6", 1.34)
    bases["A"] = a

    # --- Guanine ---
    g = _purine_ring()
    _add(g, "O6", "C6", 1.24)
    _add(g, "N2", "C2", 1.34)
    bases["G"] = g

    return bases


IDEAL_BASES = _build_ideal_bases()
RING_ATOMS = ["N1", "C2", "N3", "C4", "C5", "C6"]
PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T", "U"}
BASE_ATOM_NAMES = {
    b: set(coords.keys()) for b, coords in IDEAL_BASES.items()
}


def normalize_base_name(name):
    """Map a residue name like 'DA', 'DG', 'RC', 'A5', 'DT3' etc. down to
    a single-letter base code."""
    n = name.strip().upper()
    for b in "ACGTU":
        if n == b or n.replace("D", "").replace("R", "") == b:
            return b
    # last resort: last alphabetic char that's a valid base letter
    for ch in reversed(n):
        if ch in "ACGTU":
            return ch
    raise ValueError(f"Could not interpret base name: {name!r}")


# --------------------------------------------------------------------------
# Minimal PDB I/O (ATOM/HETATM only - sufficient for nucleic acid mutation)
# --------------------------------------------------------------------------

class Atom:
    __slots__ = ("record", "serial", "name", "altloc", "resname", "chain",
                 "resseq", "icode", "xyz", "occ", "bfac", "element", "raw")

    def __init__(self, line):
        self.record = line[0:6].strip()
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altloc = line[16]
        self.resname = line[17:20].strip()
        self.chain = line[21]
        self.resseq = int(line[22:26])
        self.icode = line[26]
        self.xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        self.occ = line[54:60]
        self.bfac = line[60:66]
        self.element = line[76:78] if len(line) >= 78 else ""
        self.raw = line

    def formatted(self, serial):
        return (
            f"ATOM  {serial:5d} {self.name:<4s}{self.altloc}{self.resname:>3s} "
            f"{self.chain}{self.resseq:4d}{self.icode}   "
            f"{self.xyz[0]:8.3f}{self.xyz[1]:8.3f}{self.xyz[2]:8.3f}"
            f"{self._occ():6s}{self._bfac():6s}"
            f"          {self.element:>2s}\n"
        )

    def _occ(self):
        try:
            return f"{float(self.occ):6.2f}"
        except ValueError:
            return "  1.00"

    def _bfac(self):
        try:
            return f"{float(self.bfac):6.2f}"
        except ValueError:
            return "  0.00"


def read_pdb(path):
    atoms = []
    others = []  # non ATOM/HETATM lines, kept for passthrough on write
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(Atom(line))
            elif not line.startswith("END"):
                others.append(line)
    return atoms, others


def residue_key(atom):
    return (atom.chain, atom.resseq, atom.icode)


def group_residues(atoms):
    """chain,resseq,icode -> list[Atom] in file order."""
    residues = {}
    order = []
    for a in atoms:
        k = residue_key(a)
        if k not in residues:
            residues[k] = []
            order.append(k)
        residues[k].append(a)
    return residues, order


# --------------------------------------------------------------------------
# Kabsch superposition
# --------------------------------------------------------------------------

def kabsch(mobile, target):
    """Return (R, t) such that R @ mobile[i] + t ~= target[i] (least squares)."""
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    cm = mobile.mean(axis=0)
    ct = target.mean(axis=0)
    mc = mobile - cm
    tc = target - ct
    H = mc.T @ tc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ct - R @ cm
    return R, t


# --------------------------------------------------------------------------
# Mutation-file parsing
# --------------------------------------------------------------------------

class MutationEntry:
    def __init__(self, chain, snum, icode=" ", name=None, mutation=None, lineno=0, raw=""):
        self.chain = chain
        self.snum = snum
        self.icode = icode
        self.name = name
        self.mutation = mutation
        self.lineno = lineno
        self.raw = raw


def parse_mutation_file(path):
    entries = []
    n_total = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            n_total += 1
            fields = {}
            for tok in stripped.split():
                if "=" not in tok:
                    continue
                k, v = tok.split("=", 1)
                fields[k.strip().lower()] = v.strip()

            chain = fields.get("chain")
            snum = fields.get("snum") or fields.get("seqnum")
            mutation = fields.get("mutation") or fields.get("m")

            if chain is None or snum is None or mutation is None:
                print(f"{lineno}\tmiss required fields: chain, seqnum or mutation")
                continue
            try:
                snum = int(snum)
            except ValueError:
                print(f"{lineno}\tinvalid snum: {snum!r}")
                continue

            icode = fields.get("icode", " ") or " "
            name = fields.get("name")

            try:
                mutation = normalize_base_name(mutation)
            except ValueError:
                print(f"unrecognized mutation: {lineno} {stripped!r}")
                continue

            entries.append(MutationEntry(chain, snum, icode, name, mutation, lineno, stripped))

    n_valid = len(entries)
    print(f"{n_valid} out of {n_total} mutation entries are valid.")
    if n_total and n_valid == 0:
        print(f"{path} contains no mutation entry")
    return entries


# --------------------------------------------------------------------------
# Core mutation logic
# --------------------------------------------------------------------------

def mutate_residue(residue_atoms, target_base):
    ring_coords = {}
    for a in residue_atoms:
        if a.name in RING_ATOMS:
            ring_coords[a.name] = a.xyz

    missing = [n for n in RING_ATOMS if n not in ring_coords]
    if missing:
        raise ValueError(f"residue is missing ring atoms {missing}; cannot compute frame")

    mobile = np.array([IDEAL_BASES[target_base][n] for n in RING_ATOMS])
    target = np.array([ring_coords[n] for n in RING_ATOMS])
    R, t = kabsch(mobile, target)

    ideal = IDEAL_BASES[target_base]

    all_possible_base_atoms = set()
    for coords in IDEAL_BASES.values():
        all_possible_base_atoms |= set(coords.keys())
    kept = [a for a in residue_atoms if a.name not in all_possible_base_atoms]

    template = kept[0] if kept else residue_atoms[0]

    # Preserve the file's existing naming convention (DA/DC/DG/DT for DNA,
    # RA/RC/RG/RU or bare A/C/G/U for RNA) and apply it to EVERY atom in
    # the residue, not just the newly built base atoms -- otherwise the
    # residue ends up with a mixed resname (old code on kept backbone
    # atoms, new code only on the swapped-in base atoms).
    orig_resname = template.resname.strip().upper()
    if orig_resname.startswith("D") and len(orig_resname) > 1:
        target_resname = "D" + target_base
    elif orig_resname.startswith("R") and len(orig_resname) > 1:
        target_resname = "R" + target_base
    else:
        target_resname = target_base

    new_atoms = []
    for a in kept:
        a.resname = target_resname
        new_atoms.append(a)

    for name, xyz in ideal.items():
        new_xyz = R @ xyz + t
        a = Atom.__new__(Atom)
        a.record = "ATOM"
        a.serial = 0  # renumbered on write
        a.name = name
        a.altloc = " "
        a.resname = target_resname
        a.chain = template.chain
        a.resseq = template.resseq
        a.icode = template.icode
        a.xyz = new_xyz
        a.occ = "  1.00"
        a.bfac = "  0.00"
        a.element = name[0]
        new_atoms.append(a)

    return new_atoms


def apply_mutations(atoms, entries):
    residues, order = group_residues(atoms)
    applied = []
    for e in entries:
        key = (e.chain, e.snum, e.icode if e.icode != " " else " ")
        # try exact match, then relax icode
        match_key = None
        for k in order:
            if k[0] == e.chain and k[1] == e.snum and (e.icode.strip() == "" or k[2] == e.icode):
                match_key = k
                break
        if match_key is None:
            print(f"Mutation entry {e.lineno} {e.raw!r} has no PDB residue match")
            continue

        residue_atoms = residues[match_key]
        orig_resname = residue_atoms[0].resname
        try:
            orig_base = normalize_base_name(orig_resname)
        except ValueError:
            orig_base = None

        if e.name and orig_base and normalize_base_name(e.name) != orig_base:
            print(f"Mutation entry {e.lineno} {e.raw!r}: name={e.name} does not match "
                  f"residue {orig_resname} at {match_key}")
            continue

        if orig_base == e.mutation:
            print(f"{e.lineno}: mutation is the same as original")
            continue

        try:
            new_residue_atoms = mutate_residue(residue_atoms, e.mutation)
        except ValueError as ex:
            print(f"skip {match_key}: {ex}")
            continue

        residues[match_key] = new_residue_atoms
        applied.append((match_key, orig_resname, e.mutation))
        print(f"mutated: chain={match_key[0]} snum={match_key[1]} "
              f"{orig_resname} -> {e.mutation}")

    out_atoms = []
    for k in order:
        out_atoms.extend(residues[k])
    return out_atoms, applied


def write_pdb(path, atoms, others, remarks):
    with open(path, "w") as fh:
        fh.write("REMARK    Mutated PDB File from 3DNA-compatible 'mutate_bases.py'\n")
        for r in remarks:
            fh.write(r)
        serial = 1
        for a in atoms:
            fh.write(a.formatted(serial))
            serial += 1
        fh.write("END\n")


# --------------------------------------------------------------------------
# Template ("-e") generation
# --------------------------------------------------------------------------

def write_template(pdb_path, out_path):
    atoms, _ = read_pdb(pdb_path)
    residues, order = group_residues(atoms)
    with open(out_path, "w") as fh:
        fh.write("# Empty or comment (starting with #s) lines are ignored\n\n")
        fh.write("#           chain=A snum=2 name=C mutation=G\n")
        fh.write("#   to mutate base C (on chain A and with residue number 2) to G\n")
        fh.write("# add m=BASE_NAME (up-to three letters) to an entry for mutation\n")
        for k in order:
            chain, snum, icode = k
            resname = residues[k][0].resname
            fh.write(f"chain={chain}   snum={snum:<4d}  name={resname:<3s} icode={icode}\n")
    print(f"Number of mutations: {len(order)}")
    print(f"Template written to {out_path}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def run_apply(pdb_path, mut_path, out_path):
    """Parse `mut_path`, apply its mutations to `pdb_path`, write `out_path`.
    Shared by both the legacy (-l) and the standalone CLI entry points."""
    if not pdb_path.exists():
        sys.exit(f"PDB file <{pdb_path}> does not exist!")
    if not mut_path.exists():
        sys.exit(f"Mutation file <{mut_path}> does not exist!")

    entries = parse_mutation_file(mut_path)
    if not entries:
        return

    atoms, others = read_pdb(pdb_path)
    new_atoms, applied = apply_mutations(atoms, entries)

    remarks = [f"REMARK    Mutation#{i+1} chain={k[0]} snum={k[1]} {orig} to [{new}]\n"
               for i, (k, orig, new) in enumerate(applied)]
    write_pdb(out_path, new_atoms, others, remarks)
    print(f"\n[SUCCESS] wrote mutated structure -> {out_path}")


def main():
    argv = sys.argv[1:]

    # ---------------------------------------------------------------
    # Legacy X3DNA-compatible invocation, as used by existing callers
    # (e.g. model_dna.py):
    #     mutate_bases -l mutation_file pdb_file output_file
    # ---------------------------------------------------------------
    if argv and argv[0] == "-l":
        if len(argv) != 4:
            sys.exit("usage: mutate_bases.py -l mutation_file pdb_file output_file")
        _, mut_arg, pdb_arg, out_arg = argv
        run_apply(Path(pdb_arg), Path(mut_arg), Path(out_arg))
        return

    # ---------------------------------------------------------------
    # Standalone CLI: mutate_bases.py [-e] pdb_file [mutation_file] [-o out]
    # ---------------------------------------------------------------
    ap = argparse.ArgumentParser(
        description="Reimplementation of X3DNA's mutate_bases",
        add_help=True,
    )
    ap.add_argument("-e", action="store_true",
                     help="generate a mutation-entry template instead of applying mutations")
    ap.add_argument("pdb_file")
    ap.add_argument("mutation_or_template_file", nargs="?", default=None)
    ap.add_argument("-o", "--output", default=None,
                     help="output PDB path (mutation mode only; default: <input>_mutated.pdb)")
    args = ap.parse_args(argv)

    pdb_path = Path(args.pdb_file)
    if not pdb_path.exists():
        sys.exit(f"PDB file <{pdb_path}> does not exist!")

    if args.e:
        out = Path(args.mutation_or_template_file or "mutations.dat")
        write_template(pdb_path, out)
        return

    mut_path = Path(args.mutation_or_template_file or "mutations.dat")
    out_path = Path(args.output) if args.output else pdb_path.with_name(pdb_path.stem + "_mutated.pdb")
    run_apply(pdb_path, mut_path, out_path)


if __name__ == "__main__":
    main()
