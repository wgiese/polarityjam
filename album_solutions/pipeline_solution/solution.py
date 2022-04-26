from album.runner.api import setup, get_args


env_file = """name: Vascuec-Solution
channels:
  - conda-forge
  - defaults
  - fastai
  - pytorch
dependencies:
  - python=3.8.13
  - git=2.34
  - typing-extensions=4.1.1
  - six=1.16
  - pyparsing=3.0.8
  - pillow=9.1.0
  - numpy=1.21.6
  - llvmlite=0.38
  - colorama=0.4.4
  - tqdm=4.64.0
  - tifffile=2022.4.8
  - scipy=1.8.0
  - PyWavelets=1.3.0
  - pytz=2022.1
  - python-dateutil=2.8.2
  - packaging=21.3
  - numba=0.55.1
  - networkx=2.8
  - natsort=8.1.0
  - kiwisolver=1.4.2
  - imageio=2.16.2
  - fonttools=4.32.0
  - fastremap=1.12.2
  - cycler=0.11.0
  - scikit-image=0.19.2
  - pandas=1.4.2
  - matplotlib=3.5.1
  - pytorch=1.11
  - opencv-python-headless=4.5
  - cellpose=2.0
  - cmocean=2.0
  - PyYAML=6.0
  - pip
  - pip:
    - pyunpack==0.2.2
    - GitPython==3.1
"""


def install():
    from album.runner.api import setup, get_args , get_package_path
    from git import Repo
    from pathlib import Path
    import subprocess
    import sys
    import os
    print("Installing vascu_ec and dependencies")
    vascu_repo = Path(get_package_path()).joinpath('vascu_repo')
    os.mkdir(vascu_repo)
    Repo.clone_from("https://github.com/wgiese/vascu-ec-app",vascu_repo)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", Path(vascu_repo).joinpath('vascu_ec')])


def run():
    from album.runner.api import setup, get_args , get_package_path
    from pathlib import Path
    import os
    
    #collecting argument as strings for better readability and to build the vascu_ec call string  
    run_string = str(get_args().run_command)
    parameter_string = str(get_args().parameter_file)
    input_string = str(get_args().input_file)
    output_string = str(get_args().output_path)
    filename_prefix_string = str(get_args().filename_prefix)
    key_string = str(get_args().key_csv)
    vascu_string = "vascu_ec"
    
    
    #Determine the mode to run Vascu_ec
    
    if run_string == "test":
        vascu_call = "vascu_ec_test"
        os.system(vascu_call)
    else:
    
        #Normal run
        if run_string == "run":
            #run with prefix
            if filename_prefix_string != "None":
                vascu_call = '%s %s %s %s %s --filename_prefix %s' %(vascu_string,run_string,parameter_string,input_string,output_string,filename_prefix_string)
                os.system(vascu_call)
                
            #run without prefix
            else:
                vascu_call = '%s %s %s %s %s' %(vascu_string,run_string,parameter_string,input_string,output_string)
                os.system(vascu_call)
                
        #run with keys csv
        if run_string == "run-key":
            if key_string != "None":
                vascu_call = '%s %s %s %s %s %s' %(vascu_string,run_string,parameter_string,input_string,key_string,output_string)
                os.system(vascu_call)
            else:
                print("You tried to run vascu_ec in the key mode without providing a key.csv file!")
        #Run on multiple files
        if run_string == "run-stack":
            vascu_call = '%s %s %s %s %s' %(vascu_string,run_string,parameter_string,input_string,output_string)
            os.system(vascu_call)


def prepare_test():
    return {}


def test():
    from album.runner.api import setup, get_args
    import os
    vascu_call = "vascu_ec_test"
    os.system("vascu_ec_test")


def close():
    from album.runner.api import setup, get_args
    print("closing")


setup(
    group="album",
    name="Vascuec-Solution",
    version="0.1.0",
    title="A Solution to run the Vascuec App",
    description="A Solution to run the Vascuec App.",
    authors=["Lucas Rieckert", "Jan Philipp Albrecht"],
    cite=[{
        "text": "Your first citation text",
        "doi": "your first citation doi"
    }],
    tags=["vascuec", "test"],
    license="UNLICENSE",
    documentation="",
    covers=[{
        "description": "Dummy cover image.",
        "source": "cover.png"
    }],
    album_api_version="0.3.1",
    args=[{
        "name": "run_command",
        "type": "string",
        "default": "run",
        "description": "How do you want to run vascuec? run for single tiff file, run-stack for a directory containing multiple tiff files, run-key for a csv file containing a list of  directorys containing tiff files or test for the test-suit",
		"required": True
		},
        {
        "name": "input_file",
        "type": "file",
        "description": "Path to the input tiff file (or directory of tiff files for run-stack).",
		"required": True
        },
        {
        "name": "parameter_file",
        "type": "file",
        "description": "Path to the parameter file.",
		"required": True
        },
		{
		"name": "output_path",
		"type": "string",
		"description": "Path to the directory where the output will get saved.",
		"required": True
		},
        {
		"name": "filename_prefix",
		"type": "string",
		"description": "Optional file prefix to rename the input files.",
		},
        {
		"name": "key_csv",
		"type": "file",
		"description": "CSV file containing the keys for the input.",
		}

    ],
    install=install,
    run=run,
    close=close,
    pre_test=prepare_test,
    test=test,
    dependencies={'environment_file': env_file}
)