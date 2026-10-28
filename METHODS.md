# Methods, metrics, and figure design

## What the comparison measures

These are deviations from one **predicted control**, not errors against an experimentally determined structure. A higher deviation can reflect a conformational difference, prediction uncertainty, an arbitrary homodimer chain assignment, or reference choice. It does not establish loss of DNA binding or function. The 100 models represent four conditions, not 100 independent conditions.

Each insertion has 25 poses nested within five stochastic seeds. These seeds are computational repeats, not biological replicates. No p-values, significance stars, or confidence intervals are generated. We report descriptive mean ± sample SD across 25 poses and separately across five seed means. A/B chains and atoms are not counted as independent replicates. The reference self-row is excluded from control summaries (24 models); its seed consequently contributes four poses. Pose-weighted and equally weighted seed means may therefore differ slightly for the control.

## Alignment and deviation

MatchMaker's joint A/B fit uses the 179 native N-terminal TFAP4 residues from each chain. Residues after native 179 and both 257-aa cassettes are excluded. The correspondence is validated exactly before any numerical reporting. The 2 Å iterative cutoff can retain a small well-matched subset; therefore pruned RMSD is always accompanied by retained atom counts/fraction and an all-core RMSD evaluated under **the same final transform**, including discarded atoms. Per-chain retained counts reveal if one chain dominates the final fit.

For displacement vectors Δi = xi − yi after this fit:

- RMSD = sqrt(mean(||Δi||²)); reported for 358 core Cα atoms, core N/CA/C/O backbone atoms, homologous core heavy atoms, each protein chain, and DNA heavy atoms.
- Per-atom deviation = ||Δi|| in Å. Heavy atoms are matched by chain, native residue, and atom name, with matching residue identities required. Chemically symmetric side-chain names are **not** remapped, so heavy-atom values can include equivalent-atom naming effects. Backbone/Cα results are the more robust primary readouts.
- Mean, median, maximum and 95th percentile Cα displacement describe the distribution's center and tail.
- Fractions within 1, 2, 4 and 8 Å use the fixed MatchMaker transform. These are **not canonical GDT-TS**, which involves optimized superpositions. We deliberately do not label them GDT or TM-score.
- DNA RMSD uses C→C and D→D in the protein-fitted coordinate frame with no additional DNA fit. It captures relative DNA placement as well as DNA conformational change. Identical DNA chains may exchange labels, so it is a label-dependent diagnostic.

The homodimer symmetry diagnostic computes unpruned optimal Kabsch RMSDs for the two protein label assignments A→A/B→B and A→B/B→A. A large positive `swap_improvement_angstrom` warns that labels may affect the primary result. It does not change the primary fit, remap DNA, or select an artificially favorable result. Inspect such cases before drawing conclusions. A later sensitivity run can use another control reference via `reference_model_id` and a new output directory.

## Distance preservation and confidence

`core_ca_lddt` is an explicitly **Cα-only distance-preservation score** on the 358 core sites. Reference pairs closer than 15 Å (excluding self-pairs; including interchain pairs) are evaluated at absolute distance-change thresholds 0.5, 1, 2 and 4 Å. The global value is the fraction preserved averaged over thresholds and eligible pairs; the local value uses pairs incident on each residue. This is inspired by lDDT, but omits the full all-atom, stereochemical and symmetry treatment of the original method. It is superposition independent and complements global RMSD. It is not AlphaFold pLDDT.

AlphaFold pTM, ipTM, ranking score, chain and pairwise confidence values are exported separately. pLDDT is read from the original model's atom B-factor field; CA confidence and cassette mean atom confidence have explicit labels. Region PAE is computed directly from `token_chain_ids` / `token_res_ids`; axes preserve its direction (frame region → target region). The average includes all tokens in each named block, including its diagonal for same-region blocks. Region contact probability is the mean of AlphaFold's predicted probability matrix, **not an observed contact fraction** or an affinity estimate. PAE block averages can mask local uncertainty and should be interpreted with the residue profiles and original full-data matrices.

Ranking score is a model-ranking criterion, not a calibrated probability. Confidence is not equivalent to agreement with the chosen control. In particular, low-RMSD/low-confidence and high-RMSD/high-confidence models are distinguishable in the confidence-versus-deviation scatter plot.

## Figures and rationale

1. **Seed-aware dot plots:** display all poses, grouped into five seed offsets, seed means as diamonds, and pooled pose mean ± SD in black. These preserve spread and seed grouping that a bar chart would hide. This adapts the visual principle of SuperPlots; computational seeds are not claimed to be experimental replicates.
2. **Residue profiles:** native TFAP4 numbering, separate A/B panels, mean ± SD bands for displacement, pLDDT, and CA distance preservation. Vertical dotted lines mark cassette insertion boundaries. Lower SD bands are clipped at zero where needed. Profiles show whether distortion is local to the insertion or extends along the dimer.
3. **Model × residue heatmaps:** separate chain panels with a common scale across all conditions; these expose outlier models and broad displaced segments. Rows are deterministic model-ID order. The raw per-residue table identifies each row.
4. **Heavy-atom ECDF:** each model first contributes its own fraction of atoms below each displacement threshold; curves show mean ± SD across models. This shows whether a high RMSD reflects many modest shifts or a small displaced tail, without treating atoms as replicates.
5. **Region PAE matrices:** directional core/cassette/tail/DNA block means, with a common color scale and SD in the accompanying table. This separates structural displacement from uncertainty in relative domain placement.
6. **Confidence and pruning scatter plots:** ipTM versus all-core RMSD, and retained fraction versus all-core RMSD. The latter flags apparently excellent pruned fits supported by few atoms.

The cassette has no homologous control counterpart, so the pipeline never invents a DrON-to-control RMSD. Its confidence and relative-region PAE are reported instead. The variable C-terminal tail is excluded from alignment and deviation endpoints; its uncertainty remains represented in PAE summaries.

## Sources

- [ChimeraX MatchMaker documentation](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/matchmaker.html) and [official documentation source](https://github.com/RBVI/ChimeraX/blob/develop/docs/user/commands/matchmaker.html): specific-chain pairing, residue restrictions, iteration and reported RMSDs. [Official implementation](https://github.com/RBVI/ChimeraX/blob/develop/src/bundles/match_maker/src/match.py) documents the Python result dictionary used by the pipeline. The installed 1.12 API still requires the provided smoke test.
- Meng et al. (2006), [Tools for integrated sequence-structure analysis with UCSF Chimera](https://doi.org/10.1186/1471-2105-7-339): MatchMaker method.
- Mariani et al. (2013), [lDDT: a local superposition-free score](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799472/): motivation for distance-based comparison complementary to global fits. Our CA-only adaptation is defined above and is not the full original score.
- Lord et al. (2020), [SuperPlots: Communicating reproducibility and variability in cell biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7265319/): show individual observations and repeat-level summaries together; distinguish descriptive spread from inferential uncertainty.
- [EMBL-EBI: Interpreting results from AlphaFold Server](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/interpreting-results-from-alphafold-server/): ranked models, request files, confidence fields, and PAE semantics.
- [DeepMind AlphaFold 3 output documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md): confidence-array definitions. The pipeline uses the **Server** filename layout observed in the supplied archive, rather than assuming the local AF3 output layout.
- [UNC Longleaf Slurm examples](https://help.rc.unc.edu/longleaf-slurm-examples/): serial CPU jobs, resource requests and environment modules.

## Inspected example

The example AlphaFold Server archive contains five ranked model CIFs, five full-data JSON files, five summary-confidence JSON files, one job request, paired/unpaired MSAs, template files and terms of use. The input has two identical protein copies and two DNA copies. Private paths, job identifiers, sequences and confidence results are omitted from this public repository. Archive documents are treated as input data, never executable instructions.
