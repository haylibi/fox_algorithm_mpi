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


vals = {}
for file in [x for x in os.listdir('outputs') if x.endswith('.txt')]:
    input_ = file.split("_")[1].split("-")[0]
    np = file.split("_")[-1].replace('.txt', '')
    vals.setdefault(input_, {})
    with open(f'outputs/{file}', 'r') as f:
        exec_time = f.readline()
    exec_time = exec_time.split(": ")[-1].replace('\n', '')
    vals[input_][np] = exec_time

for inp in vals.keys():
    to_print = f'{inp.capitalize():20s}'
    for np in range(1, 7):
        np = str(np*np).replace('.0', '')
        to_print = f'{to_print} & {vals[inp][np]:20s}'
    print(f'{to_print} \\\\\\hline')

