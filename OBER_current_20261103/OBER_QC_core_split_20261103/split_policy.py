"""Conservative separation of an agreed ineligible core from auxiliary QC."""
import copy

CORE_FIELDS = ('scope', 'human_any', 'qualified', 'basis', 'nonhuman_only')
DETAIL_FIELDS = ('noncausal', 'secondary_other')


def core_profile(record):
    if record['status'] == 'unresolved_qc':
        return ['unresolved_qc']
    return ['resolved'] + [record.get(k) for k in CORE_FIELDS]


def split_consensus(left, right, baseline):
    """Never qualify a new positive; retain individual outputs and open details."""
    result = copy.deepcopy(baseline)
    can_resolve = (
        baseline['status'] == 'unresolved_qc'
        and left['status'] == right['status'] == 'resolved'
        and left.get('qualified') is False
        and right.get('qualified') is False
        and not left.get('internal_issues')
        and not right.get('internal_issues')
        and all(left.get(k) == right.get(k) for k in CORE_FIELDS)
        and all(k in left and k in right for k in CORE_FIELDS)
    )
    differences = [k for k in DETAIL_FIELDS if left.get(k) != right.get(k)]
    unknowns = [k for k in DETAIL_FIELDS if left.get(k) is None or right.get(k) is None]
    if can_resolve:
        result = copy.deepcopy(left)
        result['status'] = 'qc_core_consensus_negative'
        for field in differences:
            result[field] = None
            result['source_summaries'][field] = '兩次細節分類不同，保留細節QC；未選定其中一次。'
        result['signature'] = [result.get(k) for k in CORE_FIELDS + DETAIL_FIELDS]
        for field in ('human_ids', 'qualifying_ids', 'identity_ids', 'design_ids'):
            result[field] = sorted(set(left.get(field, []) + right.get(field, [])))
        result['open_nondecision_fields'] = sorted(set(differences + unknowns))
    result['core_recovered_from_detail_disagreement'] = can_resolve
    result['detail_QC_pending'] = bool(differences or unknowns)
    result['detail_QC_fields'] = sorted(set(differences + unknowns))
    result['detail_QC_pair'] = [
        {k: copy.deepcopy(x.get(k)) for k in DETAIL_FIELDS} for x in (left, right)
    ]
    result['original_full_consensus_status'] = baseline['status']
    result['full_signature_agreement'] = left.get('signature') == right.get('signature')
    return result
