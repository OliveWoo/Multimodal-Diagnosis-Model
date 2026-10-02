from tools import mngs_big_agent as agent


def test_alias_normalization_and_display_names_are_canonical():
    assert agent.normalize_organism_name("CMV") == agent.normalize_organism_name(
        "Cytomegalovirus"
    )
    assert agent.normalize_organism_name("HHV-7") == agent.normalize_organism_name(
        "Human betaherpesvirus 7"
    )
    assert agent.normalize_organism_name("HSV-1") == agent.normalize_organism_name(
        "Human alphaherpesvirus 1"
    )
    assert agent.canonical_display_name("Cytomegalovirus") == "CMV"
    assert agent.canonical_display_name("Human betaherpesvirus 7") == "HHV-7"
    assert agent.canonical_display_name("Human alphaherpesvirus 1") == "HSV-1"
    assert agent.is_oral_upper_airway_flora_name("Veillonella atypica") is True
    assert agent.is_oral_upper_airway_flora_name("Gemella morbillorum") is False


def test_postprocess_keeps_gemella_m2_not_protected_level3():
    ranked_mngs = {
        "records": [
            {
                "specimen_code": "241114_002A",
                "candidates": [
                    {
                        "organism_name": "Human alphaherpesvirus 1",
                        "reads": 5375,
                        "source_category": "3.Virus",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.9287,
                            "possibility_level": "high",
                        },
                    },
                    {
                        "organism_name": "Gemella morbillorum",
                        "reads": 19,
                        "source_category": "1.Bac",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.501,
                            "possibility_level": "high",
                        },
                    },
                ],
            }
        ]
    }
    payload = {
        "best_available_summary": {
            "selection_mode": "Confirmed",
            "picked_count": 2,
            "picked_pathogens": [
                {
                    "organism_name": "Human alphaherpesvirus 1",
                    "classification": "Viral",
                    "picked_role": "Primary",
                    "basis_level": "Level 1",
                    "mngs_signal_tier": "M1_strong",
                    "rank_priority": "1",
                    "reads": 5375,
                },
                {
                    "organism_name": "Gemella morbillorum",
                    "classification": "Bacterial",
                    "picked_role": "Protected_Level3",
                    "basis_level": "Level 3",
                    "mngs_signal_tier": "M3_protected",
                    "rank_priority": "1",
                    "reads": 9,
                },
            ],
        },
        "final_infection_likelihood": "Possible",
        "pathogen_candidates": [
            {
                "organism_name": "Human alphaherpesvirus 1",
                "classification": "Viral",
                "integrated_causative_level": "Level 1",
                "mngs_signal_tier": "M1_strong",
                "rank_priority": "1",
                "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                "reads": 5375,
                "reads_tier": "R3_high",
                "reads_percentile": 0.9287,
                "specimen_alignment": "Aligned",
                "is_protected_pathogen": True,
                "protected_retention_condition": True,
                "module_support_summary": {"mNGS": "M1_strong"},
                "key_evidence": {"guardrail_rule": "", "level_rule": "R-S5-L1"},
            },
            {
                "organism_name": "Gemella morbillorum",
                "classification": "Bacterial",
                "integrated_causative_level": "Level 3",
                "mngs_signal_tier": "M3_protected",
                "rank_priority": "1",
                "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                "reads": 9,
                "reads_tier": "R0_trace",
                "reads_percentile": 0.501,
                "specimen_alignment": "Aligned",
                "is_protected_pathogen": True,
                "protected_retention_condition": True,
                "module_support_summary": {"mNGS": "M3_protected"},
                "key_evidence": {"guardrail_rule": "", "level_rule": "R-S5-L3"},
            },
        ],
        "excluded_candidates": [
            {
                "organism_name": "Human betaherpesvirus 7",
                "classification": "Viral",
                "observed_level": "Level 5",
            }
        ],
    }

    fixed = agent.postprocess_mngs_max_output(payload, ranked_mngs=ranked_mngs)
    gemella = fixed["pathogen_candidates"][1]
    gemella_pick = fixed["best_available_summary"]["picked_pathogens"][1]

    assert fixed["pathogen_candidates"][0]["organism_name"] == "HSV-1"
    assert fixed["excluded_candidates"][0]["organism_name"] == "HHV-7"
    assert gemella["reads"] == 19
    assert gemella["reads_tier"] == "R1_low"
    assert gemella["is_protected_pathogen"] is False
    assert gemella["protected_retention_condition"] is False
    assert gemella["mngs_signal_tier"] == "M2_moderate"
    assert gemella["integrated_causative_level"] == "Level 2"
    assert gemella_pick["basis_level"] == "Level 2"
    assert gemella_pick["mngs_signal_tier"] == "M2_moderate"
    assert gemella_pick["picked_role"] == "Secondary"


def test_postprocess_removes_invalid_gemella_oral_guardrail_and_adds_pick():
    ranked_mngs = {
        "records": [
            {
                "specimen_code": "241114_002A",
                "candidates": [
                    {
                        "organism_name": "Human alphaherpesvirus 1",
                        "reads": 5375,
                        "source_category": "3.Virus",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.9287,
                            "possibility_level": "high",
                        },
                    },
                    {
                        "organism_name": "Enterococcus faecium",
                        "reads": 83,
                        "source_category": "1.Bac",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.6891,
                            "possibility_level": "high",
                        },
                    },
                    {
                        "organism_name": "Veillonella atypica",
                        "reads": 56,
                        "source_category": "1.Bac",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.6515,
                            "possibility_level": "high",
                        },
                    },
                    {
                        "organism_name": "Gemella morbillorum",
                        "reads": 19,
                        "source_category": "1.Bac",
                        "ranking": {
                            "rank_priority": 1,
                            "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                            "reads_percentile": 0.501,
                            "possibility_level": "high",
                        },
                    },
                ],
            }
        ]
    }
    payload = {
        "best_available_summary": {
            "selection_mode": "Confirmed",
            "picked_count": 2,
            "picked_pathogens": [
                {
                    "organism_name": "HSV-1",
                    "classification": "Viral",
                    "basis_level": "Level 2",
                    "mngs_signal_tier": "M1_strong",
                    "rank_priority": "1",
                    "reads": 5375,
                },
                {
                    "organism_name": "Enterococcus faecium",
                    "classification": "Bacterial",
                    "basis_level": "Level 2",
                    "mngs_signal_tier": "M2_moderate",
                    "rank_priority": "1",
                    "reads": 83,
                },
            ],
        },
        "pathogen_candidates": [
            {
                "organism_name": "HSV-1",
                "classification": "Viral",
                "integrated_causative_level": "Level 2",
                "mngs_signal_tier": "M1_strong",
                "rank_priority": "1",
                "reads": 5375,
                "reads_tier": "R3_high",
                "specimen_alignment": "Aligned",
                "is_likely_colonizer_or_background": False,
                "module_support_summary": {"mNGS": "M1_strong"},
                "key_evidence": {"guardrail_rule": "", "level_rule": "R-S5-L2"},
            },
            {
                "organism_name": "Enterococcus faecium",
                "classification": "Bacterial",
                "integrated_causative_level": "Level 2",
                "mngs_signal_tier": "M2_moderate",
                "rank_priority": "1",
                "reads": 83,
                "reads_tier": "R1_low",
                "specimen_alignment": "Aligned",
                "is_likely_colonizer_or_background": False,
                "module_support_summary": {"mNGS": "M2_moderate"},
                "key_evidence": {"guardrail_rule": "", "level_rule": "R-S5-L2"},
            },
            {
                "organism_name": "Veillonella atypica",
                "classification": "Bacterial",
                "integrated_causative_level": "Level 4",
                "mngs_signal_tier": "M4_weak",
                "rank_priority": "1",
                "reads": 56,
                "reads_tier": "R1_low",
                "specimen_alignment": "Aligned",
                "is_likely_colonizer_or_background": True,
                "module_support_summary": {"mNGS": "M4_weak"},
                "key_evidence": {
                    "guardrail_rule": "R-S4-ORAL-UPPER_AIRWAY_HARD",
                    "level_rule": "R-S5-L4",
                },
            },
            {
                "organism_name": "Gemella morbillorum",
                "classification": "Bacterial",
                "integrated_causative_level": "Level 4",
                "mngs_signal_tier": "M4_weak",
                "rank_priority": "1",
                "reads": 19,
                "reads_tier": "R1_low",
                "specimen_alignment": "Aligned",
                "is_likely_colonizer_or_background": True,
                "module_support_summary": {"mNGS": "M4_weak"},
                "key_evidence": {
                    "guardrail_rule": "R-S4-ORAL-UPPER_AIRWAY_HARD",
                    "level_rule": "R-S5-L4",
                },
            },
        ],
        "excluded_candidates": [
            {
                "organism_name": "Veillonella atypica",
                "classification": "Bacterial",
                "observed_level": "Level 4",
                "rank_priority": "1",
                "reads": 56,
                "exclusion_reason_code": "colonizer_or_background",
            },
            {
                "organism_name": "Gemella morbillorum",
                "classification": "Bacterial",
                "observed_level": "Level 4",
                "rank_priority": "1",
                "reads": 19,
                "exclusion_reason_code": "colonizer_or_background",
            },
        ],
    }

    fixed = agent.postprocess_mngs_max_output(payload, ranked_mngs=ranked_mngs)
    candidates = {
        item["organism_name"]: item for item in fixed["pathogen_candidates"]
    }
    picked_names = [
        item["organism_name"]
        for item in fixed["best_available_summary"]["picked_pathogens"]
    ]
    excluded_names = [item["organism_name"] for item in fixed["excluded_candidates"]]

    assert candidates["Gemella morbillorum"]["mngs_signal_tier"] == "M2_moderate"
    assert candidates["Gemella morbillorum"]["integrated_causative_level"] == "Level 2"
    assert candidates["Gemella morbillorum"]["key_evidence"]["guardrail_rule"] == ""
    assert candidates["Veillonella atypica"]["mngs_signal_tier"] == "M4_weak"
    assert "Gemella morbillorum" in picked_names
    assert "Gemella morbillorum" not in excluded_names
    assert fixed["best_available_summary"]["picked_count"] == 3
