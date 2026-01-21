#!/usr/bin/env python3
import yaml
import os
from copy import deepcopy

VALID_QMD_FILES = {
    'Home': 'index.qmd',
    'QC Filtering': 'qc_filter.qmd',
}

def configure_quarto(navbar_entries: list):
    return {
        'project': {
            'type': 'website',
            'output-dir': 'report',
            'render': '*.qmd',
        },
        'format': {
            'html': {
                'theme': {
                    'light': ['flatly', 'styles.scss'],
                },
                'toc': True,
                'toc-location': 'left',
                'toc-depth': 2,
                'code-link': True,
                'code-fold': True,
                'code-line-numbers': True,
            },
        },
        'website': {
            'title': 'scRNAvigator Report',
            'navbar': {
                'background': 'primary',
                'left': navbar_entries
            }
        },
    }

def construct_navbar(qmd_files: list = []):
    qmd_basenames = [
        os.path.basename(f)
        for f in qmd_files
    ]

    navbar_list = []
    i = 0
    for title, f in VALID_QMD_FILES.items():
        if f not in qmd_basenames:
            continue
        navbar_list.append({
            'text': f'{i}. {title}',
            'href': f,
        })
        i += 1

    return navbar_list

if __name__ == '__main__':
    all_files = os.listdir()
    qmd_files = [
        f for f in all_files
        if f.endswith('.qmd')
    ]

    navbar_entries = construct_navbar(qmd_files=qmd_files)

    qmd_config = configure_quarto(navbar_entries=navbar_entries)

    config_file = '_quarto.qmd'

    with open(config_file, 'w') as f:
        yaml.safe_dump(qmd_config, f)