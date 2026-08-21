#!/usr/bin/env python
"""QA: randomized ground-truth testing of pileup_parse.

pileup_parse carries 8 hand-written self-tests. Hand-written cases only cover the
failure modes their author already thought of - and bug 4 WAS a case its author had not
thought of. So instead of adding more hand cases, this CONSTRUCTS pileup strings from
tokens whose correct answer is known by construction, then checks the parser against it.
Ground truth comes from the generator, not from a second parser that could repeat the
same mistake.

Token grammar exercised (samtools mpileup):
  . ,            ref match, fwd / rev                       -> depth +1
  A C G T a c g t   mismatch, fwd / rev                     -> depth +1, alt if matching
  N n            ambiguous                                  -> depth +1
  *              deleted base in this read                  -> depth +0
  ^X             read start, X = ANY ascii mapping-quality  -> not a base
  $              read end                                   -> not a base
  +N<seq> -N<seq>  indel payload, seq may contain . , $ ^ ACGT and digits
Adversarial choices: mapping-quality char drawn from the FULL printable range so it
lands on '.', 'T', 'a', '^', '$', '+', '-' and digits; indel payloads use the real nucleotide alphabet,
with multi-digit lengths to exercise the length parser.
"""
import sys, random, string
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else "/mnt/data/a3a")
import pileup_parse as pp

rng = random.Random(20260821)
BASES_F = "ACGT"
MAPQ = string.printable[:95]           # includes . , ^ $ + - digits ACGT acgt
PAYLOAD = "ACGTNacgtn*"   # REAL mpileup indel payloads are nucleotides only.
# A payload containing digits or +/- would make the token genuinely ambiguous
# (".-1" + "8" is indistinguishable from a 18-base deletion). samtools never emits
# that, so testing it would manufacture a false alarm rather than find a real bug.


def make_case(alt):
    s, f, r, depth = [], 0, 0, 0
    for _ in range(rng.randint(0, 14)):
        if rng.random() < 0.25:                       # read start, adversarial mapq
            s.append("^" + rng.choice(MAPQ))
        kind = rng.random()
        if kind < 0.40:
            c = rng.choice(".,"); s.append(c); depth += 1
        elif kind < 0.80:
            b = rng.choice(BASES_F)
            up = rng.random() < 0.5
            c = b if up else b.lower()
            s.append(c); depth += 1
            if b == alt:
                f += up; r += (not up)
        elif kind < 0.88:
            s.append(rng.choice("Nn")); depth += 1
        else:
            s.append("*")                              # deleted base: NOT depth
        if rng.random() < 0.20:                        # indel payload after the base
            n = rng.randint(1, 15)
            s.append(rng.choice("+-") + str(n) +
                     "".join(rng.choice(PAYLOAD) for _ in range(n)))
        if rng.random() < 0.20:
            s.append("$")
    return "".join(s), f, r, depth


bad = 0
N = 20000
for i in range(N):
    alt = rng.choice(BASES_F)
    s, ef, er, ed = make_case(alt)
    gf, gr = pp.count_alt(s, alt)
    gd = pp.count_depth(s)
    if (gf, gr, gd) != (ef, er, ed):
        bad += 1
        if bad <= 5:
            print(f"  MISMATCH alt={alt} str={s!r}")
            print(f"    expected fwd={ef} rev={er} depth={ed}")
            print(f"    got      fwd={gf} rev={gr} depth={gd}")
print(f"randomized ground-truth cases: {N:,}   mismatches: {bad}")

# targeted adversarial singles that the random generator may under-sample
targeted = [
    ("^T.,.$",            "T", (0, 0), 3),   # mapq char is 'T'
    ("^a,,",              "A", (0, 0), 2),   # mapq char is 'a'; 2 bases follow
    ("^^AT",              "T", (1, 0), 2),   # mapq char is '^'
    ("^$AT",              "T", (1, 0), 2),   # mapq char is '$'
    ("^+AT",              "T", (1, 0), 2),   # mapq char is '+'
    ("^5AT",              "T", (1, 0), 2),   # mapq char is a digit
    (".+12ACGTACGTACGT,", "A", (0, 0), 2),   # two-digit indel length
    (".+2NN,",            "A", (0, 0), 2),   # indel payload is not a base
    (".-3TTT,t",          "T", (0, 1), 3),   # deletion payload is not a base
    ("*,.T",              "T", (1, 0), 3),   # '*' is a read but NOT a usable base
    ("",                  "T", (0, 0), 0),   # empty pileup (zero coverage)
    ("^",                 "T", (0, 0), 0),   # truncated read-start at end
]
tbad = 0
for s, alt, exp_ar, exp_d in targeted:
    got = pp.count_alt(s, alt); gd = pp.count_depth(s)
    ok = (got == exp_ar and gd == exp_d)
    if not ok:
        tbad += 1
        print(f"  TARGETED FAIL {s!r} alt={alt}: got {got}/{gd}, expected {exp_ar}/{exp_d}")
print(f"targeted adversarial cases: {len(targeted)}   failures: {tbad}")
print("VERDICT:", "PASS" if bad == 0 and tbad == 0 else "*** FAIL ***")
