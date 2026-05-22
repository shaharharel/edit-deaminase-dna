"""Hybrid editor / enzyme feature schema for DeaminaFormer.

Combines:
  - Categorical editor_id (50 named editors + UNK)
  - Decomposed slots that let the model extrapolate to novel editors:
      deaminase_family, deaminase_mutations (multi-hot),
      cas_variant, pam_class, linker, ugi_count, fusions,
      localization, window_center, window_width, edit_type

Final embedding dim: 128 (after concat + MLP)
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import torch
import torch.nn as nn

# --- Vocabularies ---

DEAMINASE_FAMILIES = [
    "rAPOBEC1", "hA3A", "hA3B", "hA3G", "hAID", "PmCDA1",
    "TadA_WT", "TadA_evolved", "hADAR2", "evoCDA1", "evoAPOBEC1",
    "evoFERNY", "Anc689", "UNK",
]

# Activity-altering mutations (curated from literature)
DEAMINASE_MUTATIONS = [
    # rAPOBEC1
    "W90A", "W90Y", "R33A", "K34A", "R126E", "R132E",
    # hA3A
    "Y130F", "N57G",
    # TadA7.10 (ABE7.10 mutations from Gaudelli 2017)
    "H36L", "W23R", "P48A", "R51L", "L84F", "A106V", "D108N", "H123Y",
    "S146C", "D147Y", "R152P", "E155V", "I156F", "K157N",
    # TadA8e directed-evolution (Richter 2020)
    "A109S", "T111R", "D119N", "H122N", "Y147D", "F149Y", "T166I", "D167N",
    # Other TadA mutations / engineered variants
    "V82G", "V82S", "Q154R", "Y147R", "Y123H", "A158K",
    "F148A", "V106W", "E59A", "N108Q", "del153",
    # ADAR
    "E488Q", "T375G",
    # placeholder
    "OTHER",
]

CAS_VARIANTS = [
    "dCas9", "nCas9_D10A", "nCas9_H840A", "nCas9_NG", "SpRY",
    "xCas9", "enAsCas12a", "LbCas12a", "dCas13b", "none",
]

PAM_CLASSES = ["NGG", "NG", "NRN", "NYN", "TTTV", "none"]

LINKER_TYPES = ["XTEN16", "XTEN32", "GS3", "GS6", "intein_split", "none", "OTHER"]

EXTRA_FUSIONS = ["UGI", "UNG", "XRCC1", "Mu_Gam", "MS2", "PUF", "SunTag", "repressor"]

LOCALIZATIONS = ["NLS", "NES", "both", "none"]

EDIT_TYPES = ["C_to_T", "A_to_G", "C_to_G", "RNA_A_to_I", "RNA_C_to_U"]


@dataclass
class EditorConfig:
    """One row in the BE catalog. Describes a single named editor."""
    editor_id: str                       # e.g. "BE4max", "ABE8e"; "UNK" for novel
    deaminase_family: str                 # from DEAMINASE_FAMILIES
    deaminase_mutations: list[str] = field(default_factory=list)
    cas_variant: str = "nCas9_D10A"
    pam_class: str = "NGG"
    linker_type: str = "XTEN16"
    linker_aa_count: int = 16
    ugi_count: int = 2
    extra_fusions: list[str] = field(default_factory=list)
    localization: str = "NLS"
    nls_copies: int = 2
    window_center: float = 6.0   # protospacer position
    window_width: float = 4.0    # editing window width
    edit_type: str = "C_to_T"


# --- The catalog (~35 named editors) ---
# Filled from Data Agent D's output
CATALOG: dict[str, EditorConfig] = {
    "BE1": EditorConfig("BE1", "rAPOBEC1", cas_variant="dCas9", ugi_count=0,
                        linker_type="XTEN16", localization="none", nls_copies=0,
                        window_center=6.0, window_width=4.0),
    "BE2": EditorConfig("BE2", "rAPOBEC1", cas_variant="dCas9", ugi_count=1,
                        localization="none", nls_copies=0),
    # BE3: Komor 2016 — single C-terminal NLS (Addgene 73021) — QA C6 fix
    "BE3": EditorConfig("BE3", "rAPOBEC1", cas_variant="nCas9_D10A", ugi_count=1,
                        linker_type="XTEN16", localization="NLS", nls_copies=1,
                        window_center=6.0, window_width=4.0),
    "BE4": EditorConfig("BE4", "rAPOBEC1", cas_variant="nCas9_D10A", ugi_count=2,
                        localization="NLS", nls_copies=1),
    "BE4max": EditorConfig("BE4max", "rAPOBEC1", cas_variant="nCas9_D10A", ugi_count=2,
                           localization="NLS", nls_copies=2),
    "AncBE4max": EditorConfig("AncBE4max", "Anc689", cas_variant="nCas9_D10A", ugi_count=2,
                              localization="NLS", nls_copies=2),
    "YE1_BE3": EditorConfig("YE1_BE3", "rAPOBEC1", ["W90Y", "R126E"],
                            cas_variant="nCas9_D10A", ugi_count=1,
                            window_center=5.5, window_width=2.0),
    "YE1_BE4": EditorConfig("YE1_BE4", "rAPOBEC1", ["W90Y", "R126E"],
                            cas_variant="nCas9_D10A", ugi_count=2,
                            window_center=5.5, window_width=2.0,
                            localization="NLS", nls_copies=2),
    "YEE_BE3": EditorConfig("YEE_BE3", "rAPOBEC1", ["W90Y", "R126E", "R132E"],
                            cas_variant="nCas9_D10A", ugi_count=1,
                            window_center=5.5, window_width=1.5),
    "evoAPOBEC1_BE4max": EditorConfig("evoAPOBEC1_BE4max", "evoAPOBEC1",
                                       cas_variant="nCas9_D10A", ugi_count=2,
                                       localization="NLS", nls_copies=2),
    "evoCDA1_BE4max": EditorConfig("evoCDA1_BE4max", "evoCDA1",
                                    cas_variant="nCas9_D10A", ugi_count=2,
                                    window_center=7.5, window_width=11.0,  # wide!
                                    localization="NLS", nls_copies=2),
    "evoFERNY_BE4max": EditorConfig("evoFERNY_BE4max", "evoFERNY",
                                     cas_variant="nCas9_D10A", ugi_count=2,
                                     localization="NLS", nls_copies=2),
    # hA3A_BE3: Wang 2018 — active window C4–C10, not C2–C13 — QA C6 fix
    "hA3A_BE3": EditorConfig("hA3A_BE3", "hA3A", cas_variant="nCas9_D10A", ugi_count=1,
                             window_center=7.0, window_width=6.0),
    "hA3A_eBE_Y130F": EditorConfig("hA3A_eBE_Y130F", "hA3A", ["Y130F"],
                                    cas_variant="nCas9_D10A", ugi_count=1),
    "eA3A_BE3": EditorConfig("eA3A_BE3", "hA3A", ["N57G"],
                              cas_variant="nCas9_D10A", ugi_count=1,
                              window_center=6.0, window_width=2.0),
    "Target_AID": EditorConfig("Target_AID", "PmCDA1", cas_variant="nCas9_D10A", ugi_count=1,
                                window_center=3.0, window_width=2.0),
    "SECURE_BE3_R33A": EditorConfig("SECURE_BE3_R33A", "rAPOBEC1", ["R33A"],
                                     cas_variant="nCas9_D10A", ugi_count=1),
    "SECURE_BE3_R33A_K34A": EditorConfig("SECURE_BE3_R33A_K34A", "rAPOBEC1",
                                          ["R33A", "K34A"],
                                          cas_variant="nCas9_D10A", ugi_count=1),
    "ABE7_10": EditorConfig("ABE7_10", "TadA_evolved",
                             ["W23R", "H36L", "P48A", "R51L", "L84F", "A106V",
                              "D108N", "H123Y", "S146C", "D147Y", "R152P",
                              "E155V", "I156F", "K157N"],
                             cas_variant="nCas9_D10A", ugi_count=0,
                             edit_type="A_to_G", window_center=5.5, window_width=3.0),
    "ABEmax": EditorConfig("ABEmax", "TadA_evolved",
                            ["W23R", "H36L", "P48A", "R51L", "L84F", "A106V",
                             "D108N", "H123Y", "S146C", "D147Y", "R152P",
                             "E155V", "I156F", "K157N"],
                            cas_variant="nCas9_D10A", ugi_count=0,
                            edit_type="A_to_G", window_center=5.5, window_width=3.0,
                            localization="NLS", nls_copies=2),
    # ABE8e: Richter 2020 — ABE7.10 TadA* mutations + 8 directed-evolution mutations
    # Citing Richter 2020 Nat Biotechnol 38:892. CRITICAL — verify with Levanon's biologist
    "ABE8e": EditorConfig("ABE8e", "TadA_evolved",
                           ["A109S", "T111R", "D119N", "H122N",
                            "Y147D", "F149Y", "T166I", "D167N"],
                           cas_variant="nCas9_D10A", ugi_count=0,
                           edit_type="A_to_G", window_center=6.0, window_width=6.0,
                           localization="NLS", nls_copies=2),
    "ABE8e_F148A": EditorConfig("ABE8e_F148A", "TadA_evolved",
                                 ["A109S", "T111R", "D119N", "H122N",
                                  "Y147D", "F148A", "F149Y", "T166I", "D167N"],
                                 cas_variant="nCas9_D10A", ugi_count=0,
                                 edit_type="A_to_G", window_center=6.0, window_width=6.0,
                                 localization="NLS", nls_copies=2),
    # ABE9: Chen 2022 — TadA8e narrowed-window with V82S in addition
    "ABE9": EditorConfig("ABE9", "TadA_evolved",
                          ["A109S", "T111R", "D119N", "H122N",
                           "Y147D", "F149Y", "T166I", "D167N", "V82S", "F148A"],
                          cas_variant="nCas9_D10A", ugi_count=0,
                          edit_type="A_to_G", window_center=5.5, window_width=1.5,
                          localization="NLS", nls_copies=2),
    "ABE8e_V106W": EditorConfig("ABE8e_V106W", "TadA_evolved",
                                 ["A109S", "T111R", "D119N", "H122N",
                                  "Y147D", "F149Y", "T166I", "D167N", "V106W"],
                                 cas_variant="nCas9_D10A", ugi_count=0,
                                 edit_type="A_to_G", localization="NLS", nls_copies=2),
    "CGBE1": EditorConfig("CGBE1", "rAPOBEC1", ["R33A"], cas_variant="nCas9_D10A",
                           ugi_count=0, extra_fusions=["XRCC1"], edit_type="C_to_G",
                           window_center=6.0, window_width=3.0),
    "GBE_AID": EditorConfig("GBE_AID", "hAID", cas_variant="nCas9_D10A", ugi_count=0,
                             extra_fusions=["UNG"], edit_type="C_to_G",
                             window_center=6.0, window_width=3.0),
    "REPAIR_v2": EditorConfig("REPAIR_v2", "hADAR2", ["E488Q", "T375G"],
                               cas_variant="dCas13b", pam_class="none", ugi_count=0,
                               edit_type="RNA_A_to_I", linker_type="none"),
    "RESCUE": EditorConfig("RESCUE", "hADAR2",
                            cas_variant="dCas13b", pam_class="none", ugi_count=0,
                            edit_type="RNA_C_to_U", linker_type="none"),
    "BEACON1": EditorConfig("BEACON1", "hA3A", cas_variant="enAsCas12a",
                             pam_class="TTTV", ugi_count=1),
    "BE_NG": EditorConfig("BE_NG", "rAPOBEC1", cas_variant="nCas9_NG",
                           pam_class="NG", ugi_count=2,
                           localization="NLS", nls_copies=2),
    "SpRY_CBE": EditorConfig("SpRY_CBE", "rAPOBEC1", cas_variant="SpRY",
                              pam_class="NRN", ugi_count=2,
                              localization="NLS", nls_copies=2),
    "SpRY_ABE8e": EditorConfig("SpRY_ABE8e", "TadA_evolved",
                                ["V82G", "Q154R", "Y147R", "Y123H", "A158K"],
                                cas_variant="SpRY", pam_class="NRN", ugi_count=0,
                                edit_type="A_to_G", localization="NLS", nls_copies=2),
    # placeholder for novel editor designs
    "UNK": EditorConfig("UNK", "UNK"),
}


# Lock the editor_id → index mapping. CRITICAL (QA C4): this must stay stable
# across catalog edits. To add a new editor, append to CATALOG but DO NOT
# re-insert in the middle, otherwise downstream embedding indices shift and
# saved checkpoints become silently wrong.
EDITOR_ID_TO_IDX = {k: i for i, k in enumerate(CATALOG)}
NUM_EDITORS = len(CATALOG)


# --- Embedding module ---

class EditorEmbedding(nn.Module):
    """Encodes an EditorConfig into a fixed-dim embedding.

    Hybrid: editor_id (16d) + decomposed slots (~80d) → concat → MLP → 128d.
    """
    def __init__(self, d_out: int = 128):
        super().__init__()
        # Vocabularies — drop padding_idx (was zeroing BE1's gradient) — QA M22 fix
        self.editor_id_emb = nn.Embedding(NUM_EDITORS, 16)
        self.deam_family_emb = nn.Embedding(len(DEAMINASE_FAMILIES), 16)
        self.deam_mutations_emb = nn.Linear(len(DEAMINASE_MUTATIONS), 16, bias=False)  # multi-hot
        self.cas_emb = nn.Embedding(len(CAS_VARIANTS), 8)
        self.pam_emb = nn.Embedding(len(PAM_CLASSES), 4)
        self.linker_type_emb = nn.Embedding(len(LINKER_TYPES), 4)
        self.fusions_emb = nn.Linear(len(EXTRA_FUSIONS), 8, bias=False)  # multi-hot
        self.localization_emb = nn.Embedding(len(LOCALIZATIONS), 4)
        self.edit_type_emb = nn.Embedding(len(EDIT_TYPES), 4)
        # Continuous slots: linker_aa, ugi_count, nls_copies, window_center, window_width = 5
        raw_dim = 16 + 16 + 16 + 8 + 4 + 4 + 8 + 4 + 4 + 5
        self.proj = nn.Sequential(
            nn.Linear(raw_dim, 256), nn.GELU(), nn.Dropout(0.1),
            nn.Linear(256, d_out),
        )

    def forward(self, cfgs: list[EditorConfig]) -> torch.Tensor:
        """cfgs: list of EditorConfig of length B. Returns (B, d_out)."""
        device = next(self.parameters()).device
        B = len(cfgs)
        # Build indices/multi-hots on the fly
        # Use locked EDITOR_ID_TO_IDX (QA C4 fix) — UNK is always at the end
        editor_id_idx = torch.tensor(
            [EDITOR_ID_TO_IDX.get(c.editor_id, EDITOR_ID_TO_IDX["UNK"]) for c in cfgs],
            device=device, dtype=torch.long)
        df_idx = torch.tensor([DEAMINASE_FAMILIES.index(c.deaminase_family) for c in cfgs],
                              device=device, dtype=torch.long)
        # multi-hot mutations
        mh_mut = torch.zeros(B, len(DEAMINASE_MUTATIONS), device=device)
        for i, c in enumerate(cfgs):
            for m in c.deaminase_mutations:
                if m in DEAMINASE_MUTATIONS:
                    mh_mut[i, DEAMINASE_MUTATIONS.index(m)] = 1.0
        cas_idx = torch.tensor([CAS_VARIANTS.index(c.cas_variant) for c in cfgs],
                               device=device, dtype=torch.long)
        pam_idx = torch.tensor([PAM_CLASSES.index(c.pam_class) for c in cfgs],
                               device=device, dtype=torch.long)
        lt_idx = torch.tensor([LINKER_TYPES.index(c.linker_type) for c in cfgs],
                              device=device, dtype=torch.long)
        # multi-hot fusions
        mh_fus = torch.zeros(B, len(EXTRA_FUSIONS), device=device)
        for i, c in enumerate(cfgs):
            for f in c.extra_fusions:
                if f in EXTRA_FUSIONS:
                    mh_fus[i, EXTRA_FUSIONS.index(f)] = 1.0
        loc_idx = torch.tensor([LOCALIZATIONS.index(c.localization) for c in cfgs],
                               device=device, dtype=torch.long)
        et_idx = torch.tensor([EDIT_TYPES.index(c.edit_type) for c in cfgs],
                              device=device, dtype=torch.long)
        cont = torch.tensor(
            [[c.linker_aa_count / 32.0, c.ugi_count / 2.0, c.nls_copies / 2.0,
              c.window_center / 10.0, c.window_width / 10.0] for c in cfgs],
            device=device, dtype=torch.float32)
        raw = torch.cat([
            self.editor_id_emb(editor_id_idx),
            self.deam_family_emb(df_idx),
            self.deam_mutations_emb(mh_mut),
            self.cas_emb(cas_idx),
            self.pam_emb(pam_idx),
            self.linker_type_emb(lt_idx),
            self.fusions_emb(mh_fus),
            self.localization_emb(loc_idx),
            self.edit_type_emb(et_idx),
            cont,
        ], dim=-1)
        return self.proj(raw)
