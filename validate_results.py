import os
import tqdm

for file in tqdm.tqdm([x for x in os.listdir('outputs') if x.endswith('.txt')]):
    with open(f'outputs/{file}', 'r') as f:
        mine = f.readlines()
    with open(f'outputs/{file.replace(".txt", "")}', 'r') as f:
        prof = f.readlines()
    equal = True
    for line_m, line_prof in zip(mine, prof):
        line_m = line_m.replace(' \n', '\n')
        equal = equal and (line_m==line_prof)
    print(f'{file} -> {equal}')