# Neuron_Model_12HH16HH_gof

Adult L5 pyramidal cell with Hodgkin-Huxley Nav1.2 and Nav1.6 (x170 baseline), plus the v3c Nav1.2 variant fits.

```
nrnivmodl Neuron_Model_12HH16HH_basaldend/mechanisms     # from this folder
```

```python
from build_cell import build_cell
cell = build_cell('M1879T_B1B2')    # 'WT' or any key of MUTANT_PARAM_DICTS
```

- `build_cell.py`: applies all densities, kinetics and dendritic rules, then returns the cell.
- `variant_params_v3c.py`: `MUTANT_PARAM_DICTS`, one line per variant (WT + 68 fits).
- `params/`: WT Nav kinetics.
- `Neuron_Model_12HH16HH_basaldend/`: hoc model, morphology and mechanisms. `Stims/` and `params/orig.csv` are read by `runModel.hoc` when it loads.

Notes:

- Use `build_cell`, not `template.hoc`: the template alone does not have the correct Nav channels or dendritic rules.
- Build one cell per process, so export each variant in its own process.
- Variants are heterozygous: `na12` is WT and `napkca` holds the variant (allele 2). The hoc still inserts `na12mut`, but its gbar is 0 (`k1422e_frac = 1`) and `build_cell` ignores it.
- `build_cell(variant, expression=1.0)`: `expression` scales total Nav1.2 density (both alleles) relative to WT, e.g. `build_cell('WT', 0.3)` is 30%.
- `K1422E_*` variants also change the `napkca` pore (pna:pk:pca = 1 : 0.191 : 0.1). `build_cell` sets this from the variant name.
- Temperature is 34 °C.
