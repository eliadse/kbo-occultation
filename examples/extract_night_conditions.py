"""
Condense each night's MAGIC CC ``.rep`` reports (weather, pyrometer, AMC,
starguider, pointing; ~300 MB per telescope) into small CSVs under
``data/observations/conditions/<night>/``, restricted to the night's
observing span +- PAD_S. The SNR study notebook reads only these caches.

Re-run after adding a night (or its .rep files) to observation_log.yaml.
"""

import os

from kbo_occultation.aux_reports import read_cc_reports, save_conditions
from kbo_occultation.observation_log import load_observation_log

PAD_S = 1800.0


def main():
    log = load_observation_log()
    for night in log.nights.values():
        if not night.cc_reports or night.conditions_dir is None:
            print(f"{night.name}: no CC reports listed, skipping")
            continue
        spans = [obs.time_span() for obs in night.observations.values() if os.path.exists(obs.path)]
        t_start = min(s[0] for s in spans) - PAD_S
        t_stop = max(s[1] for s in spans) + PAD_S
        for telescope, rep_path in night.cc_reports.items():
            if not os.path.exists(rep_path):
                print(f"{night.name} {telescope}: {rep_path} not found, skipping")
                continue
            reports = read_cc_reports(rep_path, t_start, t_stop)
            save_conditions(reports, night.conditions_dir, telescope)
            counts = ", ".join(f"{k}={len(v)}" for k, v in reports.items())
            print(f"{night.name} {telescope}: {counts}")


if __name__ == "__main__":
    main()
