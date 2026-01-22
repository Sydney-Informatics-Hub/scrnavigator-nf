#!/usr/bin/env python3
import yaml
import os

VALID_QMD_FILES_NESTED = {
    'Quality Control': [
        ('Filtering', 'qc_filter.qmd'),
        ('Clustering', 'qc_cluster.qmd'),
        ('Doublet Detection', 'qc_doublets.qmd'),
    ]
}

def configure_quarto(navbar_entries: list):
    return {
        'project': {
            'type': 'website',
            'output-dir': 'report',
            'render': ['*.qmd'],
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

    # Start with the index page
    navbar_list = [
        {'text': 'Home', 'href': 'index.qmd'}
    ]

    # Iterate through valid nested QMD files and add them if they exist
    i = 1
    for menu in ['Quality Control', 'Integration', 'Annotation', 'Analysis']:
        valid_qmd_files = VALID_QMD_FILES_NESTED.get(menu)
        if valid_qmd_files is None:
            continue
        navbar_nested_dict = {
            'text': f'{i}. {menu}',
            'menu': []
        }
        j = 1
        for title, f in valid_qmd_files:
            if f not in qmd_basenames:
                continue
            navbar_nested_dict['menu'].append({
                'text': f'{j}. {title}',
                'href': f,
            })
            j += 1
        navbar_list.append(navbar_nested_dict)

    return navbar_list

if __name__ == '__main__':
    all_files = os.listdir()
    qmd_files = [
        f for f in all_files
        if f.endswith('.qmd')
    ]

    navbar_entries = construct_navbar(qmd_files=qmd_files)

    qmd_config = configure_quarto(navbar_entries=navbar_entries)

    config_file = '_quarto.yml'

    with open(config_file, 'w') as f:
        yaml.safe_dump(qmd_config, f)