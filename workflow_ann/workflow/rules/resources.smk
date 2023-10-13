'''
Utilities for resources allocation.
'''

import collections
import re
import yaml


def allocate_resources():
    '''
    Get resources requirements for all rules from the resources config file
    specified in config.yaml. Overrides these values with values defined directly
    in the config file.
    '''
    parsed_resources = collections.defaultdict(dict)
    cfg_resources = yaml.safe_load(open(config['resources_info']))
    threads_default = config['presets']['threads'][cfg_resources['default']['threads']]
    mem_mb_default = config['presets']['mem_mb'][cfg_resources['default']['mem_mb']]
    runtime_default = config['presets']['runtime'][cfg_resources['default']['runtime']]
    for rule, cfg in cfg_resources.items():
        if rule == 'default':
            continue
        try:
            parsed_resources[rule]['threads'] = config['presets']['threads'][cfg['threads']]
        except KeyError:
            print(f'Could not find threads resources preset <{cfg["threads"]}> for rule <{rule}>. Using default ({threads_default})')
            parsed_resources[rule]['threads'] = threads_default
        try:
            parsed_resources[rule]['mem_mb'] = config['presets']['mem_mb'][cfg['mem_mb']]
        except KeyError:
            print(f'Could not find mem_mb resources preset <{cfg["mem_mb"]}> for rule <{rule}>. Using default ({mem_mb_default})')
            parsed_resources[rule]['mem_mb'] = mem_mb_default
        try:
            parsed_resources[rule]['runtime'] = config['presets']['runtime'][cfg['runtime']]
        except KeyError:
            print(f'Could not find runtime resources preset <{cfg["runtime"]}> for rule <{rule}>. Using default ({runtime_default})')
            parsed_resources[rule]['runtime'] = runtime_default
    if config['resources']:
        for rule, specs in config['resources'].items():
            try:
                for spec, value in specs.items():
                    parsed_resources[rule][spec] = value
            except KeyError:
                print(f'Invalid resource <{spec}> or rule name <{rule}> in resources section of config file')
    config['resources'] = parsed_resources
    config['resources']['default'] = cfg_resources['default']


def get_threads(rule):
    '''
    Get number of threads for a rule from the config dictionary.
    '''
    try:
        threads = config['resources'][rule]['threads']
    except KeyError:
        threads = config['resources']['default']['threads']
    return threads


def get_mem(rule, attempt):
    '''
    Get memory requirement for a rule from the config dictionary.
    Memory is increased 1.5x per attempt.
    '''
    try:
        mem_mb = config['resources'][rule]['mem_mb']
    except KeyError:
        mem_mb = config['presets']['mem_mb'][config['resources']['default']['mem_mb']]
    if isinstance(mem_mb, (int, float)):
        mem_mb = int(mem_mb * (1 + (attempt - 1) / 2))
    elif mem_mb.isdigit():
        mem_mb = int(int(mem_mb) * (1 + (attempt - 1) / 2))
    elif mem_mb[-1] in ('G', 'M', 'K'):  # Careful, this cannot be used in resources (not an int)
        tmp = float(mem_mb[:-1]) * (1 + (attempt - 1) / 2)
        mem_mb = f'{tmp}{mem_mb[-1]}'
    return mem_mb


def get_runtime(rule, attempt):
    '''
    Get runtime requirement for a rule from the config dictionary.
    Runtime is increased 1.5x per attempt.
    '''
    try:
        runtime = config['resources'][rule]['runtime']
    except KeyError:
        runtime = config['presets']['runtime'][config['resources']['default']['runtime']]
    if isinstance(runtime, int):
        time = runtime
    else:
        try:
            d, h, m, s = (int(f) for f in re.split(':|-', runtime))
            time = ((((d * 24) + h) * 60) + m) * 60 + s
        except ValueError:
            print(f'Invalid runtime format for rule <{rule}>: <{runtime}>')
            time = 3600
    # time = int(time * (1 + (attempt - 1) / 2))  # Remove this at the moment since curnagl has such a low max time limit
    return time


# Parse all resources specification sources and populate the config dictionary with resources for each rule
allocate_resources()
