"""Guideline-informed research routes. No model, labels, files or old A1/B verdicts."""
import math
import re
import unicodedata


FEATURES = ['recorded_pneumonia', 'acute_respiratory_deterioration', 'infectious_symptoms',
    'immunocompromising_condition', 'immunosuppressant_exposure', 'structural_lung_disorder',
    'compatible_pulmonary_imaging', 'diffuse_groundglass', 'nodular_or_cavitary_imaging']
METHODS = ['culture', 'targeted_pcr', 'panel_pcr', 'antigen', 'microscopy_histology', 'antibody', 'unknown']
COMPARTMENTS = ['lower_respiratory', 'upper_respiratory', 'blood', 'urine', 'csf', 'tissue', 'other', 'unknown']
CLINICAL_ALIASES = {
    'ecoli': 'escherichiacoli', 'pjp': 'pneumocystisjirovecii',
    'hsv1': 'humanalphaherpesvirus1', 'herpessimplexvirustype1': 'humanalphaherpesvirus1',
    'herpessimplexvirus1': 'humanalphaherpesvirus1',
    'cmv': 'humanbetaherpesvirus5', 'humancytomegalovirus': 'humanbetaherpesvirus5',
    'vzv': 'humanalphaherpesvirus3', 'varicellazostervirus': 'humanalphaherpesvirus3',
    'sarscov2': 'sarscov2', 'severeacuterespiratorysyndromerelatedcoronavirus': 'sarscov2',
    'mrsa': 'staphylococcusaureus', 'candidaglabrata': 'nakaseomycesglabratus',
}


def key(value):
    return re.sub('[^a-z0-9]', '', unicodedata.normalize('NFKC', str(value)).casefold())


def clinical_key(value):
    k = key(value)
    return CLINICAL_ALIASES.get(k, k)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def present(ledger, name, primary=False):
    f = next(f for f in ledger['features'] if f['feature'] == name)
    return f['status'] == 'present' and bool(f['evidence']) and (not primary or any(
        ledger['fact_roles'][e['fact_id']] == 'primary_window' for e in f['evidence']))


def event_gate(ledger, mode):
    # Same phenotype cannot establish organism attribution. This is event only.
    recorded = present(ledger, 'recorded_pneumonia')
    signs = present(ledger, 'acute_respiratory_deterioration') or present(ledger, 'infectious_symptoms')
    image = present(ledger, 'compatible_pulmonary_imaging', primary=True)
    return image and signs if mode == 'strict' else recorded or (image and signs)


def profile(name, category):
    k = clinical_key(name)
    if k.startswith(('candida', 'nakaseomyces')):
        return 'respiratory_candida_caution'
    if k.startswith(('humanalphaherpesvirus', 'humanbetaherpesvirus', 'humangammaherpesvirus')):
        return 'latent_herpes_caution'
    if 'fung' in category.casefold() or k == 'pneumocystisjirovecii':
        return 'fungal_host_route'
    return 'unclassified_no_special_prior'


def quant_branch(b, reads, percentile):
    return (b['compartment'] == 'lower_respiratory' and b['temporal_role'] == 'primary_window'
        and b['qc_status'] in ['evaluable', 'partial_evaluable']
        and finite(b['reads']) and b['reads'] >= reads
        and b['code_priority'] in [1, 2] and finite(b['percentile_all_positive'])
        and b['percentile_all_positive'] >= percentile)


def matched_observations(name, ledger):
    # Group-only data never identify an exact species/type. Explicitly group-named
    # candidates may match that same group, but cannot type one of its members.
    is_group = any(s in name.casefold() for s in ['group', 'complex', 'spp.', 'rhinovirus/enterovirus'])
    return [o for o in ledger['observations'] if clinical_key(o['organism_as_written']) == clinical_key(name)
        and (o['organism_scope'] == 'species' or is_group and o['organism_scope'] == 'group')]


def decide(candidate, ledger, spec):
    observations = matched_observations(candidate['organism'], ledger)
    positives = [o for o in observations if o['detection'] == 'positive'
        and ledger['fact_roles'][o['fact_id']] == 'primary_window']
    direct = [o for o in positives if o['compartment'] == 'lower_respiratory'
        and o['method'] in ['culture', 'targeted_pcr', 'panel_pcr', 'microscopy_histology']]
    invasive = [o for o in positives if o['invasive_or_cytopathic_evidence']
        and o['method'] == 'microscopy_histology' and o['compartment'] in ['lower_respiratory', 'tissue']]
    blood_bridge = any(o['compartment'] == 'blood' and o['method'] in ['culture', 'microscopy_histology'] for o in positives)
    organism_profile = profile(candidate['organism'], candidate['category'])
    event = event_gate(ledger, spec['event_mode'])
    q = [b['seq_id'] for b in candidate['branches'] if quant_branch(b, spec['min_reads'], spec['min_percentile'])]
    moderate = [b['seq_id'] for b in candidate['branches'] if quant_branch(b, 100, .8)]
    host = present(ledger, 'immunocompromising_condition') or present(ledger, 'immunosuppressant_exposure')
    phenotype = present(ledger, 'diffuse_groundglass', True) or present(ledger, 'nodular_or_cavitary_imaging', True)
    host_route = bool(spec['host_route'] and organism_profile == 'fungal_host_route' and host and phenotype and moderate)
    routes = []
    if spec['direct_route'] and direct: routes.append('D_direct_respiratory_test')
    if spec['quant_route'] and q: routes.append('Q_same_branch_quantitative')
    if host_route: routes.append('H_host_phenotype_fungal')
    # Prespecified caution is a research gate; not a claim that absent tissue = no
    # infection. It is ablated on/off in the entire frozen parameter grid.
    caution_block = False
    if spec['cautions'] and organism_profile == 'respiratory_candida_caution':
        caution_block = not (invasive or blood_bridge and direct)
    if spec['cautions'] and organism_profile == 'latent_herpes_caution':
        caution_block = not invasive
    add = bool(event and routes and not caution_block)
    reasons = ([] if event else ['pulmonary_event_not_established_for_this_route'])
    if not routes: reasons.append('no_enabled_evidence_route_met')
    if caution_block: reasons.append('low_specificity_detection_requires_more_attribution')
    return {'add': add, 'routes': routes, 'event_gate': event, 'organism_profile': organism_profile,
        'caution_block': caution_block, 'direct_fact_ids': sorted({o['fact_id'] for o in direct}),
        'invasive_fact_ids': sorted({o['fact_id'] for o in invasive}), 'quant_seq_ids': q,
        'host_route': host_route, 'blocking_or_missing': reasons,
        'interpretation': 'research rescue/defer, not proven infection or negative diagnosis'}


def specifications():
    common = {'min_reads': 1000, 'min_percentile': .8, 'event_mode': 'recorded',
        'host_route': False, 'cautions': True, 'direct_route': True, 'quant_route': False}
    rows = [{'id': 'C0_baseline', **common, 'direct_route': False},
            {'id': 'C1_direct', **common},
            {'id': 'C2_direct_quant', **common, 'quant_route': True},
            {'id': 'C3_primary', **common, 'quant_route': True, 'host_route': True}]
    for reads in [10, 100, 1000]:
        for percentile in [.8, .95]:
            for mode in ['recorded', 'strict']:
                for host in [False, True]:
                    for caution in [False, True]:
                        rows.append({'id': f'G_r{reads}_p{int(percentile*100)}_{mode}_h{int(host)}_c{int(caution)}',
                            **common, 'min_reads': reads, 'min_percentile': percentile,
                            'event_mode': mode, 'host_route': host, 'cautions': caution, 'quant_route': True})
    return rows
