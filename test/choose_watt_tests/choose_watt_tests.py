import argparse
from ruamel.yaml import YAML
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("--changed-workflows", help="Name of workflow(s) which have changed", nargs='*')
parser.add_argument("--changed-jsons", help="Name of json(s) which have changed", nargs='*')
parser.add_argument("--watt-config", help="WATT config yaml")
parser.add_argument("--womtool-jar")

#copied from watt
def resolve_path(path: str) -> str:
    """
    Given test path, return an absolute path corresponding to the intended path.
    If the given path is a relative path (does not begin with / on unix), it is taken to be relative to the current working directory.
    If the given path is an absolute path (begins with / on unix), then it is taken to either be relative to a git repo if the current
    working directory is inside a git repo, or absolute on the host system if the current working directory is not inside a git repo or
    the corresponding file in the git repo does not exist.
    This should allow users running tests on local machines to resolve the correct paths, even if the repo has a
    root given by a non-root dir on the local system.
    """
    # If the path is not absolute, we assume it is relative to current working directory
    if not os.path.isabs(path):
        return os.path.abspath(path)

    # If the path is absolute, we need to check if it is inside a git repo
    in_repo = False
    in_root = False
    working_dir = os.getcwd()
    while not in_repo and not in_root:
        if os.path.exists(os.path.join(working_dir, os.path.join('.git'))):
            in_repo = True
        else:
            new_working_dir = os.path.abspath(os.path.join(working_dir, '..'))
            in_root = working_dir == new_working_dir
            working_dir = new_working_dir
    if in_repo:
        in_repo_path = os.path.join(working_dir, path.removeprefix('/'))
        if os.path.exists(in_repo_path):
            return in_repo_path
        else:
            return path
    else:
        return path

def get_wdl_dependencies(womtool_run: subprocess.CompletedProcess):
    return [line.decode() for line in womtool_run.stdout.splitlines() if line.endswith(b'.wdl')]

if __name__ == '__main__':
    args = parser.parse_args()

    try:
        yaml = YAML(typ='safe')
        with open(args.watt_config) as file:
            config = yaml.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Cannot find configuration file at path: {args.config}")
    
    tests_to_run = set()
    
    if args.changed_workflows or args.changed_jsons: 
        changed_wdls = {resolve_path(path) for path in args.changed_workflows} if args.changed_workflows else set()
        changed_jsons = {resolve_path(path) for path in args.changed_jsons} if args.changed_jsons else set()

        for wf, wf_tests_info in config.items():
            wdl_path = resolve_path(wf_tests_info['path'])
            if wdl_path in changed_wdls:
                tests_to_run.add(wf)
            else:
                for test, test_info in wf_tests_info['tests'].items():
                    if resolve_path(test_info['test_inputs']) in changed_jsons or resolve_path(test_info['expected_outputs']) in changed_jsons:
                        tests_to_run.add(wf)
            if wf not in tests_to_run:
                womtool_run = subprocess.run(['java', "-jar", args.womtool_jar, "validate", "-l", wdl_path], capture_output=True)
                dependencies = get_wdl_dependencies(womtool_run)
                for dependency in dependencies:
                    if dependency in changed_wdls:
                        tests_to_run.add(wf)
                        #don't need to check other dependencies if wf already added to tests_to_run
                        break
        
    print(" ".join(tests_to_run))
