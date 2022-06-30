# Pynteny: synteny hmm searches made easy

```
    ____              __                  
   / __ \__  ______  / /____  ____  __  __
  / /_/ / / / / __ \/ __/ _ \/ __ \/ / / /
 / ____/ /_/ / / / / /_/  __/ / / / /_/ / 
/_/    \__, /_/ /_/\__/\___/_/ /_/\__, /  
      /____/                     /____/   
```


## Installation

Download (or fork) repo to local directory:

```bash
git clone https://github.com/Robaina/Pynteny.git
```

cd to downloaded repo and install conda enviroment:

```bash
cd Pynteny
conda env create -f environment.yml    
```

Install pynteny in conda enviroment:

```bash
conda activate pynteny
python setup.py install
```

Check that installation worked fine:

```bash
pynteny --help
```