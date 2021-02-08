#!/usr/bin/env python3

import os
from difflib import SequenceMatcher


path_to_input_folder = input('path to folder with input files: ')
path_to_output_folder = input('path to folder with output files: ')
cap = input('cap sequence: ')

def search_similar_str(text, substr):
    """Поиск наиболее похожей подстроки в тексте.
    Принимает на вход текст и искомую подстроку.
    Возвращает коэффициент схожести и список похожих подстрок"""
    if substr in text:
        return 1, [substr]

    leader = (0, set())
    for window in (text[i:i + len(substr)] for i in
                   range(len(text) - len(substr) + 1)):
        seq = SequenceMatcher(a=window, b=substr)
        score = seq.ratio()
        if score > leader[0]:
            leader = (score, {window})
        elif score == leader[0]:
            leader[1].add(window)
        if leader[0] == 1:
            return leader[0], list(leader[1])

    return leader[0], list(leader[1])

# list of file names in directory
for _, _, f in os.walk(path_to_input_folder):
    input_folder = f
    break
for file in input_folder:
    with open(f'{path_to_input_folder}/{file}') as infile:
        content = infile.readlines()
    for i in range(1, len(content) + 1, 2):
        found = list(search_similar_str(content[i], cap))
        if found[1][0] != 0:
            content[i] = content[i].replace(found[1][0], 'N' * len(found[1][0]), 1)
    with open(f'{path_to_output_folder}/masked_{file}', mode='w') as outfile:
        outfile.writelines(content)












