"""Correct samtools mpileup base-string parsing.

THREE bugs this replaces (found by QA before s6 ever ran):
 1. counting only uppercase "T" loses reverse-strand supporting reads -- and loses
    them STRAND-ASYMMETRICALLY. That is the mechanism behind the 21.8:1 skew that
    contaminated the earlier Selict stream in this project.
 2. same defect mirrored on minus-strand sites (count("a") misses "A").
 3. raw .count() is corrupted by indel tokens (+2AG, -1T) and read-start markers
    (^X, where X is an arbitrary ASCII char that may itself be "T" or "a").

Returns (n_fwd, n_rev) for the requested alt base so strand balance stays
auditable instead of assumed.
"""
import re

_INDEL = re.compile(r"[+-](\d+)")


def clean(bases):
    """Strip read-start/end markers and indel payloads from a pileup string."""
    out = []
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == "^":          # read start: skip the mapping-quality char too
            i += 2
            continue
        if c == "$":          # read end
            i += 1
            continue
        if c in "+-":         # indel: +<len><seq> / -<len><seq>
            m = _INDEL.match(bases, i)
            if m:
                i = m.end() + int(m.group(1))
                continue
        out.append(c)
        i += 1
    return "".join(out)


def count_alt(bases, alt_base):
    """(forward, reverse) read counts supporting alt_base. alt_base upper-case."""
    b = clean(bases)
    return b.count(alt_base.upper()), b.count(alt_base.lower())


def count_depth(bases):
    """Reads usable at the site after cleaning (matches + mismatches)."""
    b = clean(bases)
    return sum(1 for c in b if c in ".,ACGTNacgtn")


if __name__ == "__main__":
    # self-tests: each case targets one of the three bugs
    assert count_alt(".,.T,t", "T") == (1, 1), "must count BOTH strands"
    assert count_alt("^T.,.$", "T") == (0, 0), "^X mapping-qual char is not a base"
    assert count_alt(".+2AG.,", "A") == (0, 0), "indel payload is not a base"
    assert count_alt(".-1T.,t", "T") == (0, 1), "deletion payload is not a base"
    assert count_alt("TTtt", "T") == (2, 2)
    assert count_depth(".,.T,t") == 6
    assert count_depth("^T.,.$") == 3
    print("pileup_parse self-tests PASS")
