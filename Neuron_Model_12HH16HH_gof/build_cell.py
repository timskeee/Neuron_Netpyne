"""Build the adult L5 pyramidal cell (12HH16HH, x170 baseline), WT or a v3c Nav1.2 heterozygote.

    nrnivmodl Neuron_Model_12HH16HH_basaldend/mechanisms     # once, from this folder
    from build_cell import build_cell
    cell = build_cell('M1879T_B1B2')                         # any key of MUTANT_PARAM_DICTS

The cell is a singleton (h.cell): build one per process.
"""
import json
import math
import os

from neuron import h

from variant_params_v3c import MUTANT_PARAM_DICTS

HERE = os.path.dirname(os.path.abspath(__file__))
_built = False

# Whole-type density multipliers
NAV12 = 1.21 * 1.3 * 0.75
NAV16 = 1.43 * 1.2 * 1.5
# Regional multipliers on the hoc base densities
DEND_NAV12 = 1 * 1.3 * 0.5 * 3.3 * 0.9
SOMA_NAV12 = 2.56 * 1.3 * 0.9
SOMA_NAV16 = 0.8
AIS_NAV12 = 5.76 * 1.5 * 0.6
AIS_NAV16 = 0.648 * 3 * 1.5 * 0.75
DEND_NAV16_GBAR = 5.05e-06 * 3.3                 # Nav1.6 gbar within 20 um of the soma, 0 beyond


def _set_globals():
    h.dend_na12 = 2.48e-03 * DEND_NAV12 * DEND_NAV12
    h.dend_na16 = 0
    h.dend_k = 0.0043685576 * (1 * 5 * 0.005)
    h.soma_na12 = 3.24e-02 * SOMA_NAV12 * SOMA_NAV12
    h.soma_na16 = 7.88e-02 * SOMA_NAV16 * SOMA_NAV16
    h.soma_K = 0.21330453 * 1.25 * (0.022 * 10 * 10)
    h.ais_na12 = AIS_NAV12 / NAV12
    h.ais_na16 = AIS_NAV16 / NAV16
    h.ais_ca = 0.0010125926 * 0.0455 * 21.5
    h.ais_KCa = 0.0009423347 * 0.25 * 0.0709 * (0.25 * 2)
    h.node_na = 0.9934221
    h.axon_KP = 0.048675
    h.axon_KT = 0.004465
    h.axon_K = 0.89699364 * 2.1 * (1 * 0.4)
    h.axon_LVA = 0.00034828275
    h.axon_HVA = 1.05e-05 * 0.0455
    h.axon_KCA = 0.4008224 * 0.25 * 0.0709
    h.gpas_all = 1.34e-05
    h.cm_all = 1.6171424


def _load(name):
    p = json.load(open(os.path.join(HERE, 'params', name + '.txt')))
    p.pop('_frame')
    return p


def _write(mech, params):
    for sec in h.allsec():
        if h.ismembrane(mech, sec=sec):
            for seg in sec:
                m = getattr(seg, mech)
                for k, v in params.items():
                    if hasattr(m, k):            # Ena exists only on na12mut
                        setattr(m, k, v)


def _scale_gbar(mech, fac):
    for sec in h.allsec():
        if h.ismembrane(mech, sec=sec):
            for seg in sec:
                m = getattr(seg, mech)
                m.gbar *= fac


def _refine_dendrites(d_lambda=0.1, f=100.0):
    for sec in h.allsec():
        if 'apic' in sec.name() or 'dend' in sec.name():
            lam = 1e5 * math.sqrt(sec.diam / (4.0 * math.pi * f * sec.Ra * sec.cm))
            sec.nseg = max(int(sec.L / (d_lambda * lam) + 0.9) | 1, 1)


def _na_taper(old, d, m=24.0, a=150.0, b=300.0):
    if d <= a:
        return old * m
    if d >= b:
        return old
    return old * (1.0 + (m - 1.0) * ((b - d) / (b - a)))


def _im_band(old, d, m=300, a=150.0, b=600.0, s=40.0):
    w = (1.0 / (1.0 + math.exp(-(d - a) / s))) * (1.0 / (1.0 + math.exp((d - b) / s)))
    return old * (1.0 + (m - 1.0) * w)


def _apical_rules():
    apic = [sec for sec in h.allsec() if 'apic' in sec.name()]
    base = {sec.name(): (sec(0.5).SK_E2.gSK_E2bar, sec(0.5).Im.gImbar) for sec in apic}
    _refine_dendrites()
    h.distance(sec=h.cell.soma[0])
    for sec in apic:
        sk0, im0 = base[sec.name()]
        for seg in sec:
            d = h.distance(seg.x, sec=sec)
            for mech in ('na12', 'na12mut'):
                m = getattr(seg, mech)
                m.gbar = _na_taper(m.gbar, d)
            seg.Im.gImbar = _im_band(im0, d) * 0.7
            seg.SK_E2.gSK_E2bar = sk0 * 0.30 * (0.35 if d > 800.0 else 1.0)


def build_cell(variant='WT', expression=1.0):
    """Returns the cell. WT allele on na12; variant allele on na12mut.

    expression: total Nav1.2 density (both alleles) as a fraction of the WT cell, e.g. 0.3 = 30%.
    """
    global _built
    if _built:
        raise RuntimeError('build_cell can run only once per process')
    _built = True
    allele2 = dict(MUTANT_PARAM_DICTS[variant])
    gs = allele2.pop('gbar_scaling')
    wt = _load('na12annaTFHH')

    cwd = os.getcwd()
    os.chdir(os.path.join(HERE, 'Neuron_Model_12HH16HH_basaldend'))
    h.load_file(1, 'runModel.hoc')
    _set_globals()
    h.working()
    os.chdir(cwd)

    for sec in h.allsec():
        if h.ismembrane('SKv3_1', sec=sec):
            sec.mtaumul_SKv3_1 = 6 * 0.85
    h.working()

    # Nav kinetics and densities
    _write('na12', wt)
    _write('na12mut', {**wt, **allele2})
    _scale_gbar('na12', NAV12 * expression)
    _scale_gbar('na12mut', NAV12 * gs * expression)
    _write('na16', _load('na16HH_TF'))
    _write('na16mut', _load('na16mutHH_TF'))
    _scale_gbar('na16', NAV16)
    _scale_gbar('na16mut', NAV16)

    # Nav1.6 in dendrites only near the soma
    h.distance(0, 0.5, sec=h.cell.soma[0])
    for sec in h.allsec():
        if 'dend' in sec.name() or 'apic' in sec.name():
            for seg in sec:
                gbar = DEND_NAV16_GBAR if h.distance(seg) <= 20 else 0
                for mech in ('na16', 'na16mut'):
                    if hasattr(seg, mech):
                        getattr(seg, mech).gbar = gbar

    _apical_rules()

    # Allele 2 pore: na12mut is Na-selective by default (pna:pk:pca = 1:0:0). Only K1422E changes it.
    if variant.split('_')[0] == 'K1422E':
        for sec in h.allsec():
            for seg in sec:
                if hasattr(seg, 'na12mut'):
                    seg.na12mut.pna, seg.na12mut.pk, seg.na12mut.pca = (1.0, 0.1910, 0.1)
    return h.cell
