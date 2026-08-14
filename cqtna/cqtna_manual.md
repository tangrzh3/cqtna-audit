# CQTNA — the three diagnostics no tool can run for you

`cqtna.py` automates five of the eight diagnostics. These three it cannot, and
the report names them in every run rather than omitting them, because a
nomination audited on five of eight is audited on five of eight.

Each entry states what is needed, what the check is, and what a failure looks
like — the last of these matters most, since a check whose failure mode is
undefined is not a check.

---

## (ii) Colocalisation, with multiple-signal modelling and matched LD

**Needs** regional summary statistics for both exposure and outcome, and an LD
reference matched to the outcome cohort. A candidate list does not contain these.

**Do** run colocalisation rather than treating SMR/HEIDI as sufficient; report
how many variants entered each HEIDI test; and prefer a method that admits more
than one causal variant per region, with in-sample LD where the outcome provider
publishes it.

**Failure looks like** a region assigned to distinct causal variants by
colocalisation while HEIDI passes it, which is the signature of a neighbouring
gene's eQTL being tagged by a large local signal it does not share. Also treat
as failure any credible-set count that keeps rising as you raise the allowed
number of signals: proxy LD splits one signal into many, and a count that never
saturates is a diagnostic of the LD reference, not of the region.

**Do not** report a colocalisation posterior for a region in which the outcome
GWAS resolves no credible set at all. That posterior is computed where there is
nothing to colocalise with.

---

## (vii) Cell-level matching on lineage composition

**Needs** the single-cell data and the split itself.

**Do** match on lineage composition as well as sequencing depth, **at the cell
level**, whenever cells are divided by a gene score.

**Failure looks like** an apparent phenotype that disappears once lineage is
matched. A cluster-level control can pass while the cell-level control fails,
because contaminating cells are distributed within clusters rather than forming
their own; passing the cluster-level version is therefore not evidence.

**Note** that any score containing lineage-associated genes will split lineage
as well as state, so removing lineage genes from the score is part of the check
and not an alternative to it.

---

## (viii) Endpoint definitions, before comparing across releases

**Needs** the phenotype definitions themselves, code by code. Release notes do
not reliably summarise them.

**Do** verify that the endpoint is the same one before treating two releases of
a GWAS resource as a power series. Where the definition differs, the comparison
is still worth making, but it is a transfer test rather than a power point, and
the direction of every difference should be written down before the result is
seen.

**Failure looks like** phenotype drift scored as instability of the candidate
list — the very quantity the comparison exists to measure, moved in the direction
that makes a list look less reproducible than it is. It also runs the other way:
a definition that silently widens can make a list look more stable.

**Worked example of getting this wrong.** A release described its successor
endpoint as "including Hilmo", the hospital discharge register, which reads as a
newly widened case definition. The code-level definitions showed the cases were
identically defined and that what had changed was the rule screening cancers out
of the controls. The wrong conclusion was drawn from that one line before the
definitions were checked.

---

## What the tool reports about itself

Every `cqtna.py` run ends with this list and with the sentence that no gene can
reach `target-supported` from the tool alone. That is deliberate. The three
checks above are the ones that would license the word, and they are the three a
candidate list cannot answer.
