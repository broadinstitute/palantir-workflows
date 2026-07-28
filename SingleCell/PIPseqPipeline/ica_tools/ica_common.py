"""Shared helpers for the ICA API scripts in this directory (export_pipeline_to_ica.py, start_analysis.py)."""

import os

API_URL = 'https://ica.illumina.com/ica/rest/api'

PROJECT_NAMES_AND_IDS = {
    'BCL Shared Development': '4ab33fe6-c169-4dc9-928d-6ce7a8062d34',
    'MDL Single Cell Dev': 'de460091-a137-4848-8e23-c59b851d7425',
}


def load_api_key():
    api_key_path = os.path.expanduser('~/.icav2/api_key.txt')
    return open(api_key_path).read().strip()


def prompt_choice(prompt, options):
    """Prints a numbered list of options and asks the user to pick one.

    Returns the chosen 0-based index, or None if the user aborted.
    """
    print(prompt)
    for i, option in enumerate(options, start=1):
        print(f'  {i}. {option}')
    choice = input('Enter the number of your choice (or anything else to abort): ')
    if not choice.isdigit() or int(choice) < 1 or int(choice) > len(options):
        return None
    return int(choice) - 1
