import argparse
from ruamel.yaml import YAML
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("--changed-workflows", help="Name of workflow(s) which have changed", nargs='*')
parser.add_argument("--watt-config", help="WATT config yaml")
parser.add_argument("--womtool-jar")

#copied from watt
def resolve_relative_path(rel_path: str) -> str:
    """
    Given test path, check if it should be interpreted as a relative path inside a repo or absolute path on system.
    Return value is an absolute path on the host system pointing to the file at the given path.
    This should allow users running tests on local machines to resolve the correct paths, even if the repo has a
    root given by a non-root dir on the local system.
    """
    in_repo = False
    in_root = False
    working_dir = os.path.curdir
    while not in_repo and not in_root:
        if os.path.exists(os.path.join(working_dir, os.path.join('.git'))):
            in_repo = True
        else:
            new_working_dir = os.path.abspath(os.path.join(working_dir, '..'))
            in_root = working_dir == new_working_dir
            working_dir = new_working_dir
    if in_repo:
        return os.path.join(working_dir, rel_path.removeprefix('/'))
    else:
        return rel_path

def get_wdl_dependencies(womtool_run: subprocess.CompletedProcess):
    womtool_stdout = womtool_run.stdout
    prior_bytes = b'List of Workflow dependencies is:\n'
    start_byte_index = womtool_stdout.find(prior_bytes) + len(prior_bytes)
    dependencies = womtool_stdout[start_byte_index:].decode().split('\n')
    return [d for d in dependencies if d.endswith(".wdl")] 
    

    

if __name__ == '__main__':
    args = parser.parse_args()

    try:
        yaml = YAML(typ='safe')
        with open(args.watt_config) as file:
            config = yaml.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Cannot find configuration file at path: {args.config}")
    
    tests_to_run = set()
    
    if args.changed_workflows: 
        changed_wdls = {resolve_relative_path(path) for path in args.changed_workflows}

        for test, test_info in config.items():
            wdl_path = resolve_relative_path(test_info['path'])
            if wdl_path in changed_wdls:
                tests_to_run.add(test)
            else:
                womtool_run = subprocess.run(['java', "-jar", args.womtool_jar, "validate", "-l", wdl_path], capture_output=True)
                dependencies = get_wdl_dependencies(womtool_run)
                for dependency in dependencies:
                    if dependency in changed_wdls:
                        tests_to_run.add(test)
                        break
        
        print(" ".join(tests_to_run))
