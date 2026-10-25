# Scientific and Computational Decisions

## Project: EBV Antigen Region Conservation Analysis

This document records the main methodological decisions, quality-control criteria, reference definitions, and limitations of the computational analysis.

The purpose is to make the workflow transparent and to support the interpretation and reproduction of the reported results.

## 1. Selection of Antigen Regions

Three antigen regions were investigated because of their relevance to published studies of EBV-specific immune responses associated with multiple sclerosis.

| Protein | Reference Accession | Positions | Reference Peptide |
|---|---|---|---|
| EBNA1 | YP_401677.1 | 386–405 | SQSSSSGSPPRRPPPGRRPF |
| EBNA3A | YP_401669.1 | 325–333 | FLRGRAYGI |
| BZLF1 | YP_401673.1 | 54–64 | EPLPQGQLTAY |

**Decision:** Use fixed reference accessions and one-based reference-relative positions throughout the analysis.

**Rationale:** Explicit reference definitions reduce ambiguity when comparing protein sequences and reporting amino acid substitutions.

## 2. EBNA3A Reference Discrepancy

The peptide `FLRGRAYGL` has been reported in EBV-specific immunological studies, including research involving individuals with MS.

However, the selected reference protein `YP_401669.1` contains `FLRGRAYGI` at positions 325–333.

Positions 193–201 in this reference correspond to a different peptide, `FVHLQATLG`.

**Decision:** Analyze positions 325–333 relative to `YP_401669.1` and explicitly document the difference between the reference peptide and the published peptide.

**Rationale:** Reference-dependent coordinates must not be silently treated as interchangeable.

**Remaining limitation:** The discrepancy in published and reference-relative numbering has not been fully resolved. Direct equivalence between the literature coordinates and the selected reference coordinates should not be assumed without further sequence-level verification.

### Reference Sequence Verification

The selected EBNA3A reference sequence (`YP_401669.1`, 944 amino acids) was directly examined using Biopython.

The following observations were confirmed:

- Positions 193–201: `FVHLQATLG`
- Positions 325–333: `FLRGRAYGI`
- The published peptide `FLRGRAYGL` differs from the reference peptide at amino acid position 9 (I to L).

## 3. Initial Protein Retrieval

Candidate EBV protein records were retrieved from NCBI Protein using Biopython Entrez.

Initial retrieval counts:

| Protein | Records |
|---|---:|
| EBNA1 | 696 |
| EBNA3A | 148 |
| BZLF1 | 591 |
| Total | 1,435 |

**Decision:** Preserve the downloaded FASTA files and accession metadata.

**Rationale:** Retaining the original candidate collection allows downstream inclusion decisions and results to be traced back to specific database records.

**Important consideration:** Search results were not automatically treated as confirmed records of the target protein.

## 4. Protein Identity Screening

Protein descriptions were initially checked using target-specific name patterns.

Records without a recognized target name were flagged for further identity review.

Additional screening considered sequence length and nonstandard amino acid symbols.

**Decision:** Exclude unresolved identity records from the primary analysis rather than treating all NCBI search hits as valid target proteins.

**Rationale:** Database searches may retrieve proteins other than the intended target, potentially introducing unrelated sequences into multiple sequence alignment and conservation calculations.

## 5. BZLF1 Identity Validation

The initial BZLF1 screening classified:

- 418 records as included.
- 170 records as requiring identity review.
- 2 records as too short.
- 1 record as containing nonstandard residues.

Further examination showed that some records had ambiguous descriptions despite strong sequence similarity to the BZLF1 reference.

Pairwise protein alignment was used to prioritize records with at least 80% reference coverage and 80% amino acid identity under the selected screening procedure.

Of 170 reviewed records, 99 met these exploratory similarity thresholds.

GenBank gene annotations were then examined for these 99 candidates.

- 97 records contained an explicit BZLF1 gene annotation.
- 2 records remained unresolved by this annotation check.

**Decision:** Incorporate the 97 gene-annotated records into the primary BZLF1 sequence collection while retaining the unresolved records outside the primary analysis.

**Rationale:** Combining sequence similarity with available gene annotation provides stronger evidence of protein identity than an ambiguous protein description alone.

**Limitation:** Gene annotations were not independently experimentally verified. The original name-matched records were also not all subjected to the same additional validation procedure.

## 6. Sequence Quality and Alignment Eligibility

The screened multiple sequence alignment workflow used sequences containing the 20 standard amino acid symbols.

Sequences with nonstandard symbols were flagged for review.

Identical sequence strings were deduplicated before multiple sequence alignment.

**Decision:** Preserve accession-level inclusion records while aligning unique protein sequences.

**Rationale:** This reduces redundant alignment input without losing the link between aligned sequences and their original accessions.

**Limitation:** Excluding records with ambiguous amino acid symbols, even outside the target region, may remove otherwise informative records.

## 7. Multiple Sequence Alignment

PyFAMSA was used for multiple sequence alignment.

The reference protein was explicitly included in each alignment.

**Decision:** Map target antigen regions according to the positions of nongap reference residues in the alignment.

**Rationale:** Reference-relative mapping provides a consistent coordinate system when comparing sequences with different lengths or alignment gaps.

**Limitation:** Gaps in aligned target regions may reflect incomplete sequence coverage, biological insertions or deletions, or alignment uncertainty. The current mapping-based exclusion procedure does not fully distinguish these possibilities.

## 8. Conservation Metrics

Two complementary levels of analysis were used:

**Accession-record level:** Each eligible database record contributes to the exact-match calculation.

**Unique-protein-sequence level:** Identical protein sequences are counted once.

A peptide was considered an exact match only when its entire extracted target region matched the reference peptide.

Position-wise amino acid frequencies and Shannon entropy were used to describe variability.

**Decision:** Report both accession-level and unique-sequence exact-match frequencies.

**Rationale:** Repeated identical sequences can affect descriptive frequency estimates.

**Limitation:** Neither accession counts nor deduplicated sequence counts represent independent patient samples or unbiased estimates of viral population frequencies.

## 9. Final Computed Results

| Protein | Comparable Records | Exact Matches | Exact-Match Frequency |
|---|---:|---:|---:|
| EBNA1 | 59 | 59 | 100.00% |
| EBNA3A | 58 | 7 | 12.07% |
| BZLF1 | 443 | 421 | 95.03% |

After deduplication of identical protein sequences:

| Protein | Unique Eligible Sequences | Exact-Match Frequency |
|---|---:|---:|
| EBNA1 | 35 | 100.00% |
| EBNA3A | 39 | 12.82% |
| BZLF1 | 63 | 82.54% |

These values describe the specific sequences retained by the computational workflow.

## 10. Biological Interpretation

The project investigates conservation and variability of selected EBV antigen regions.

The association between EBV and MS provides scientific context for the selection of these regions.

**Decision:** Interpret the findings as descriptive viral sequence analysis rather than evidence that a specific peptide variant causes or protects against MS.

**Rationale:** The dataset does not contain patient-level MS status, immune-response measurements, or functional antigen-recognition experiments.

## 11. Reproducibility

The principal analysis pipeline consists of:

1. `align_screened.py`
2. `map_regions.py`
3. `variant_summary.py`
4. `conservation.py`
5. `visuals.py`

These scripts can be executed sequentially using `run_pipeline.py`.

The pipeline was successfully rerun using the existing local input files, reproducing the reported conservation summaries and regenerating the figures.

**Environment used for the verified local run:**

- Python 3.15.0rc3
- Biopython 1.88
- PyFAMSA 0.7.0
- Matplotlib 3.11.2

**Remaining work:** Verify clean-environment installation and execution, finalize accession provenance documentation, and review the handling of incomplete and gap-containing target regions before final release.