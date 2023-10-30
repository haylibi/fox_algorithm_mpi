import os

for file in [x for x in os.listdir('outputs') if x.endswith('.txt')]:
    with open(f'outputs/{file}', 'r') as f:
        mine = f.readlines()
    prof_file = file.split("-")[0]
    prof_file = prof_file.replace('_input', '')
    with open(f'outputs/{prof_file}', 'r') as f:
        prof = f.readlines()
    equal = True
    for line_m, line_prof in zip(mine, prof):
        line_m = line_m.replace(' \n', '\n')
        equal = equal and (line_m==line_prof)
    if 'input5' in file: continue
    print(f'{file:30s} -> {equal}')